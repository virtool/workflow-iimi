FROM debian:buster as prep
WORKDIR /build
RUN apt-get update && apt-get install -y cmake gcc g++ make unzip wget zlib1g-dev
RUN wget https://zlib.net/pigz/pigz-2.8.tar.gz
RUN tar -xvf pigz-2.8.tar.gz
WORKDIR /build/pigz-2.8
RUN make
WORKDIR /build
RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
RUN unzip fastqc_v0.11.9.zip
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.5.2/bowtie2-2.5.2-linux-x86_64.zip
RUN unzip bowtie2-2.5.2-linux-x86_64.zip
RUN mkdir bowtie2
RUN cp bowtie2-2.5.2-linux-x86_64/bowtie* bowtie2
RUN wget https://github.com/relipmoc/skewer/archive/0.2.2.tar.gz
RUN tar -xf 0.2.2.tar.gz
WORKDIR /build/skewer-0.2.2
RUN make
RUN mv skewer /build
WORKDIR /build
RUN apt-get install -y bzip2 libcurses-ocaml-dev libbz2-dev liblzma-dev
RUN wget https://github.com/samtools/samtools/releases/download/1.19/samtools-1.19.tar.bz2
RUN tar -xjf samtools-1.19.tar.bz2
WORKDIR /build/samtools-1.19
RUN ./configure --prefix /build/samtools
RUN make
RUN make install

FROM python:3.12.3-bookworm as rbuild
WORKDIR /build
RUN apt-get update && apt-get install -y libdeflate0 r-base r-cran-littler
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install()"
RUN R -e "BiocManager::install(c('Biostrings', 'Rsamtools', 'GenomicAlignments'))"
RUN R -e "install.packages(c('data.table', 'dplyr', 'mltools', 'randomForest', 'xgboost', 'knitr'))"
COPY ./iimi ./iimi
RUN R -e "install.packages('iimi', repos = NULL, type = 'source')"

FROM python:3.12.3-bookworm as base
WORKDIR /app
RUN apt-get update && apt-get install -y libdeflate0 r-base r-cran-littler
COPY --from=prep /build/bowtie2/* /usr/local/bin/
COPY --from=prep /build/pigz-2.8/pigz /usr/local/bin/pigz
COPY --from=prep /build/samtools/bin/samtools /usr/local/bin/samtools
COPY --from=prep /build/skewer /usr/local/bin/
RUN curl -sSL https://install.python-poetry.org | python3 - --version 1.8.3
ENV PATH="/root/.local/bin:$PATH"
COPY poetry.lock pyproject.toml ./
RUN poetry install --without dev --no-root
COPY --from=rbuild /usr/local/lib/R/site-library /usr/local/lib/R/site-library
COPY run.r utils.py workflow.py VERSION* ./
RUN poetry install
ENTRYPOINT ["poetry", "run"]

FROM base as test
RUN poetry install
COPY conftest.py ./
COPY tests ./tests
COPY example ./example
ENTRYPOINT ["poetry", "run"]
