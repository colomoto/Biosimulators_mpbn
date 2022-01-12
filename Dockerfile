# Base OS
FROM python:3.9-slim-buster

ARG VERSION="0.1"
ARG SIMULATOR_VERSION="1.6"

# metadata
LABEL \
    org.opencontainers.image.title="mpbn" \
    org.opencontainers.image.version="${SIMULATOR_VERSION}" \
    org.opencontainers.image.description="Brief Python implementation of Most Permissive Boolean Networks" \
    org.opencontainers.image.url="https://mpbn.readthedocs.io/" \
    org.opencontainers.image.documentation="https://mpbn.readthedocs.io/" \
    org.opencontainers.image.authors="BioSimulators Team <info@biosimulators.org>" \
    org.opencontainers.image.vendor="BioSimulators Team" \
    org.opencontainers.image.licenses="SPDX:GPL-3.0-only" \
    \
    base_image="python:3.9-slim-buster" \
    version="${VERSION}" \
    software="mpbn" \
    software.version="${SIMULATOR_VERSION}" \
    about.summary="Brief Python implementation of Most Permissive Boolean Networks" \
    about.home="https://mpbn.readthedocs.io/" \
    about.documentation="https://mpbn.readthedocs.io/" \
    about.tags="BioSimulators,mathematical model,logical model,systems biology,computational biology,SBML,SED-ML,COMBINE,OMEX" \
    maintainer="Loïc Paulevé <loic.pauleve@labri.fr>"

# Install GINsim
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends \
        openjdk-11-jre-headless \
    && pip install ginsim \
    && python -m ginsim_setup \
    && apt-get autoremove -y \
    && apt clean -y \
    && rm -rf /var/lib/apt/lists/*

# Install mpbn
RUN pip install mpbn==${SIMULATOR_VERSION}

# fonts for matplotlib
RUN apt-get update -y \
    && apt-get install -y --no-install-recommends libfreetype6 \
    && rm -rf /var/lib/apt/lists/*

# Copy code for command-line interface into image and install it
COPY . /root/Biosimulators_mpbn
RUN pip install /root/Biosimulators_mpbn \
    && rm -rf /root/Biosimulators_mpbn

ENV VERBOSE=0 \
    MPLBACKEND=PDF

# Declare the environment variables that the simulation tool supports (e.g., ALGORITHM_SUBSTITUTION_POLICY) and their default values

# Entrypoint
ENTRYPOINT ["biosimulators-{my-simulator}"]
CMD []
