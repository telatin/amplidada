# ── Stage 1: Install pixi and resolve environment ────────────────────────────
FROM ubuntu:24.04 AS builder

ARG PIXI_VERSION=v0.41.4
ARG TARGETARCH

# Install curl + ca-certs (needed to fetch pixi and conda packages)
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Install pixi
RUN curl -fsSL https://pixi.sh/install.sh | bash \
    && mv /root/.pixi/bin/pixi /usr/local/bin/pixi

WORKDIR /app

# Copy the pixi manifest (and lock file if present)
COPY pixi.* .

# Install the environment (--frozen uses lock file if present, otherwise solves)
RUN pixi install --manifest-path pixi.toml

# ── Stage 2: Lean runtime image ───────────────────────────────────────────────
FROM ubuntu:24.04

RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 \
        libgfortran5 \
        ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy pixi binary
COPY --from=builder /usr/local/bin/pixi /usr/local/bin/pixi

# Copy the resolved pixi project (manifest + installed .pixi/ environment)
WORKDIR /app
COPY --from=builder /app /app

# Expose the pixi-managed R on PATH
ENV PATH="/app/.pixi/envs/default/bin:${PATH}"
ENV R_HOME="/app/.pixi/envs/default/lib/R"

# Smoke-test: confirm R and the three Bioconductor packages load
RUN Rscript -e "\
    pkgs <- c('dada2', 'phyloseq', 'microbiome'); \
    for (p in pkgs) { \
        library(p, character.only = TRUE); \
        cat(p, as.character(packageVersion(p)), 'OK\n') \
    }"

CMD ["R"]
