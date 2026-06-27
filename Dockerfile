FROM mambaorg/micromamba:1.5-bullseye-slim

ARG ENV_FILE
COPY ${ENV_FILE} /tmp/environment.yml

ARG MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH="$MAMBA_ROOT_PREFIX/bin:$PATH"

RUN micromamba install -c defaults -c bioconda -c conda-forge -y -n base git gcc_linux-64 gxx_linux-64 && \
    micromamba install -y -n base -f /tmp/environment.yml && \
    micromamba clean --all --yes

WORKDIR /home/mambauser/app

RUN git clone https://github.com/Rinoahu/POEM_py3k.git . \
    && git checkout 09aa0973ccc9f44db916040e5115de7660e26173

