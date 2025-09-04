FROM python:3.13-bookworm AS rbuild
WORKDIR /build
RUN apt-get update && apt-get install -y libdeflate0 r-base r-cran-littler
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install()"
RUN R -e "BiocManager::install(c('Biostrings', 'Rsamtools', 'GenomicAlignments'))"
RUN R -e "install.packages(c('data.table', 'dplyr', 'mltools', 'randomForest', 'xgboost', 'knitr', 'remotes', 'MTPS', 'R.utils'))"
RUN R -e "remotes::install_github('virtool/iimi', ref='1bfac4d429a72e01a6899a288cb6d6a511fb303c')"

FROM python:3.13-bookworm AS deps
WORKDIR /app
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/bowtie2/2.5.4/bowtie* /usr/local/bin/
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/pigz/2.8/pigz /usr/local/bin/
COPY --from=ghcr.io/virtool/tools:1.1.0 /tools/samtools/1.22.1/bin/samtools /usr/local/bin/

FROM python:3.13-bookworm AS uv
WORKDIR /app
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:${PATH}" \
    UV_CACHE_DIR='/tmp/uv_cache'
COPY uv.lock pyproject.toml README.md ./
RUN uv sync

FROM deps AS base
WORKDIR /app
RUN apt-get update && apt-get install -y libdeflate0 r-base r-cran-littler
ENV PATH="/app/.venv/bin:/root/.local/bin:${PATH}"
COPY --from=uv /app/.venv /app/.venv
COPY --from=rbuild /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY run.r utils.py workflow.py VERSION* ./

FROM deps AS test
WORKDIR /app
RUN apt-get update && apt-get install -y libdeflate0 r-base r-cran-littler
RUN curl -LsSf https://astral.sh/uv/install.sh | sh
ENV PATH="/root/.local/bin:$PATH" \
    UV_CACHE_DIR='/tmp/uv_cache'
COPY uv.lock pyproject.toml ./
COPY README.md ./
RUN uv sync
COPY --from=rbuild /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY conftest.py ./
COPY tests ./tests
COPY example ./example
COPY run.r utils.py workflow.py VERSION* ./
