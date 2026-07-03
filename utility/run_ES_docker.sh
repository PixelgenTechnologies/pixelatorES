#!/bin/bash

# Exit on error
set -e

CONTAINER_IMAGE="ghcr.io/pixelgentechnologies/pixelatores:main" # Default container image
MOUNT_DIRS=()
NAME=""
CSV_FILE=""
RUN_TERMINAL="false"    # Explicitly initialize to false
QUARTO_PARAMS=()        # Initialize an array to store multiple -P parameters

# Usage function (updated)
usage() {
  echo ""
  echo "Usage: $0 -d <pixelator-dir> -n <output-name> -s <samplesheet-csv> [-i <container-image>] [-P <param>] [-D]"
  echo ""
  echo "  -d <pixelator-dir>    : The Pixelator output directory to mount in the container (required)."
  echo "  -n <output-name>      : Name or path of the output HTML file (required)."
  echo "  -s <samplesheet-csv>  : The ES samplesheet CSV file (required)."
  echo "  -i <ghcr-image>       : The Docker image to use (optional)."
  echo "  -P <parameter>        : Additional parameters to pass to Quarto (can be used multiple times)."
  echo "  -D <debug>            : Run a bash terminal inside the container instead of experiment-summary."
  exit 1
}

# Parse command-line options
while getopts "d:n:s:i:P:D" opt; do
  case $opt in
    P)
      # Store -P and the argument as separate elements for clean bash expansion
      QUARTO_PARAMS+=("-P" "$OPTARG") 
      ;;
    d)
      MOUNT_DIRS+=("$OPTARG")
      ;;
    n)
      NAME="$OPTARG"
      ;;
    s)
      CSV_FILE="$OPTARG"
      ;;
    i)
      CONTAINER_IMAGE="$OPTARG"
      ;;
    D)
      RUN_TERMINAL="true"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
      usage
      ;;
  esac
done

# Check for required arguments only if not running in terminal mode
if [ "$RUN_TERMINAL" == "false" ]; then
  if [ ${#MOUNT_DIRS[@]} -eq 0 ] || [ -z "$NAME" ] || [ -z "$CSV_FILE" ]; then
    RED='\033[0;31m'
    NC='\033[0m'
    printf "${RED}Error${NC}: When not using -D, all required arguments (-d, -n, -s) must be provided.\n" >&2
    usage
  fi
fi

# Robustly handle output paths and names (resolves relative path issues)
if [ -n "$NAME" ]; then
  OUTPUT_DIR=$(realpath "$(dirname "$NAME")")
  OUTPUT_BASE=$(basename "$NAME")
else
  OUTPUT_DIR=$(pwd)
  OUTPUT_BASE="summary"
fi

# Build base docker args
docker_args=(
  --rm
  -it
  -v "$OUTPUT_DIR":/workspace/output
)

# Only attempt to mount CSV if it was actually provided
if [ -n "$CSV_FILE" ]; then
  docker_args+=(--mount type=bind,source="$(realpath "$CSV_FILE")",target=/workspace/samplesheet.csv,readonly)
fi

# Add each data directory as a subfolder in /workspace/data
for dir in "${MOUNT_DIRS[@]}"; do
  base=$(basename "$(realpath "$dir")")
  docker_args+=("-v" "$(realpath "$dir"):/workspace/data/$base")
done

echo "Pulling image: $CONTAINER_IMAGE..."
docker pull "$CONTAINER_IMAGE"

# Run the container
if [ "$RUN_TERMINAL" == "true" ]; then
  echo "Running bash terminal in container..."
  docker run \
    "${docker_args[@]}" \
    "$CONTAINER_IMAGE" \
    bash
else
  echo "Running experiment-summary in container..."
  docker run \
    "${docker_args[@]}" \
    "$CONTAINER_IMAGE" \
    bash -c "experiment-summary \
      -P sample_sheet=/workspace/samplesheet.csv -P data_folder=/workspace/data \
      \"${QUARTO_PARAMS[@]}\" && \
      mv /workspace/inst/quarto/pixelatorES.html /workspace/output/${OUTPUT_BASE}.html && \
      chown $(id -u):$(id -g) /workspace/output/${OUTPUT_BASE}.html"
fi

exit 0