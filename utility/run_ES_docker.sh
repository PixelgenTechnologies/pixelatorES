#!/bin/bash

# Exit on error
set -e

CONTAINER_IMAGE="ghcr.io/pixelgentechnologies/pixelatorES:main" # Default container image
MOUNT_DIRS=()
OUTPUT_PATH=""
CSV_FILE=""
SKIP_CLEANUP="false"  # Default: do cleanup
MISSING_FLAGS="false"   # Flag to track missing required arguments
QUARTO_PARAMS=() # Initialize an array to store multiple -P parameters


# Usage function (updated)
usage() {
  echo ""
  echo "Usage: $0 -i <container> -d <pixelator-dir> -n <output-name> -s <samplesheet-csv> [-c <container-name>] [-skip]"
  echo ""
  echo "  -d <pixelator-dir>	  : The Pixelator output directory to mount in the container (required)."
  echo "  -n <output-name>  	  :  Path of the output HTML file without extension (required)."
  echo "  -s <samplesheet-csv> 	: The ES samplesheet CSV file (required)."
  echo "  -i <ghcr-image>	      : The Docker image to use."
  echo "  -P <parameter>	      : Additional parameters to pass to Quarto (can be used multiple times)."
  echo "  -D <debug>            : Run a bash terminal inside the container instead of experiment-summary."
  exit 1
}

# Parse command-line options
while getopts "d:n:s:i:P:D" opt; do
  case $opt in
    P)
      QUARTO_PARAMS+=("-P $OPTARG") # Append each -P parameter to the array
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
      MISSING_FLAGS="true" # Set the flag
      usage
      ;;
    :)
      echo "Option -$OPTARG requires an argument." >&2
       MISSING_FLAGS="true" # Set the flag
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
    MISSING_FLAGS="true" # Set flag BEFORE calling usage
    usage
    exit 1 #Must exit after usage, otherwise trap will run
  fi
fi


# Build base docker args
docker_args=(
  --rm
  -it
  -v "$(dirname "$NAME")":/workspace/output
  --mount type=bind,source="$(realpath "$CSV_FILE")",target=/workspace/samplesheet.csv,readonly
)

# Add each data directory as a subfolder in /workspace/data
for dir in "${MOUNT_DIRS[@]}"; do
  base=$(basename "$(realpath "$dir")")
  docker_args+=("-v" "$(realpath "$dir"):/workspace/data/$base")
done

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
      ${QUARTO_PARAMS[@]} \
      --output "$NAME.html" && \
      mv /workspace/inst/quarto/$NAME.html /workspace/output/ && \
      chown $(id -u):$(id -g) output/$NAME.html"
fi

exit 0
