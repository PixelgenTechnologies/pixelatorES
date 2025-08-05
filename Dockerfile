FROM ghcr.io/pixelgentechnologies/pixelatorr:main
ARG QUARTO_VERSION="1.5.54"
ARG GITHUB_PAT=
# Set the environment variable for the GitHub Personal Access Token
# This is needed to access private repositories during the build process
ENV GITHUB_PAT=${GITHUB_PAT}

RUN apt-get update && apt-get install -y curl && rm -rf /var/lib/apt/lists/* &&\
    export QUARTO_VERSION=${QUARTO_VERSION} &&\
    mkdir -p /opt/quarto/${QUARTO_VERSION} &&\
    curl -o quarto.tar.gz -L "https://github.com/quarto-dev/quarto-cli/releases/download/v${QUARTO_VERSION}/quarto-${QUARTO_VERSION}-linux-amd64.tar.gz" &&\
    tar -zxvf quarto.tar.gz -C "/opt/quarto/${QUARTO_VERSION}" --strip-components=1 &&\
    rm quarto.tar.gz
ENV PATH=/opt/quarto/$QUARTO_VERSION/bin:$PATH

WORKDIR /workspace

# Install dependencies
COPY DESCRIPTION /workspace/
RUN R -e "pak::pak()"

# Move the quarto files into the container
COPY inst /workspace/inst

# Install the package
COPY R /workspace/R
COPY NAMESPACE /workspace/NAMESPACE
RUN R -e "devtools::install()"

# Setup runable entrypoint
COPY experiment-summary /usr/local/bin/experiment-summary
