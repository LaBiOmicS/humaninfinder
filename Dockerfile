# HumaninFinder container image.
#
# Installs the *released* package from PyPI, so the image holds exactly what `pip install
# humanin-finder==<version>` gives users, with nothing from the working tree. Published by
# .github/workflows/container.yml as ghcr.io/labiomics/humaninfinder:<version>.
#
#   docker build --build-arg VERSION=1.0.0b4 -t humaninfinder:1.0.0b4 .
#   docker run --rm -v "$PWD":/data humaninfinder:1.0.0b4 predict -i genome.fasta -o results
#
# Files are read from and written to /data (the working directory); mount your directory there.

FROM python:3.12-slim

ARG VERSION
RUN test -n "$VERSION" || (echo "build with --build-arg VERSION=<released version>" >&2; exit 1)

LABEL org.opencontainers.image.title="HumaninFinder" \
      org.opencontainers.image.description="Homology-first location and annotation of Humanin and SHLP1-6 in mitogenomes and NUMTs" \
      org.opencontainers.image.version="${VERSION}" \
      org.opencontainers.image.source="https://github.com/LaBiOmicS/humaninfinder" \
      org.opencontainers.image.url="https://labiomics.github.io/humaninfinder/" \
      org.opencontainers.image.licenses="MIT"

ENV PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1

# Wheels only: pyhmmer ships manylinux wheels for amd64 and arm64, so nothing is compiled here.
RUN pip install --only-binary=:all: "humanin-finder==${VERSION}" \
 && humanin-finder setup

# Run as an unprivileged user. With `--user "$(id -u):$(id -g)"` the outputs belong to you.
RUN useradd --create-home --uid 1000 humanin
USER humanin
WORKDIR /data

ENTRYPOINT ["humanin-finder"]
CMD ["--help"]
