# Epicurus Neo full pipeline image: brings every external tool + the epicurus-neo package.
# Build:  docker build -t epicurus-neo .
# Run:    docker run --rm -v "$PWD":/work epicurus-neo run-pipeline --config /work/patient.yaml --output-dir /work/out
#
# Reference bundles are large and licensing varies, so they are NOT baked in; mount
# them at run time and pass --bundle-dir / references.bundle_dir.
FROM mambaorg/micromamba:1.5-jammy

COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml
RUN micromamba install -y -n base -f /tmp/environment.yml && micromamba clean --all --yes

# Install the package itself into the base env.
COPY --chown=$MAMBA_USER:$MAMBA_USER . /opt/epicurus-neo
ARG MAMBA_DOCKERFILE_ACTIVATE=1
RUN pip install --no-deps -e /opt/epicurus-neo

WORKDIR /work
ENTRYPOINT ["micromamba", "run", "-n", "base", "epicurus-neo"]
CMD ["--help"]
