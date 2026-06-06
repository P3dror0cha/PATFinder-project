FROM mambaorg/micromamba:1.5-bullseye-slim

ARG ENV_FILE

COPY ${ENV_FILE} /tmp/env.yml

RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"