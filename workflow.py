import asyncio
from operator import itemgetter
from pathlib import Path

from pyfixtures import fixture
from virtool_workflow import step, RunSubprocess
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.indexes import WFIndex
from virtool_workflow.data.ml import WFMLModelRelease
from virtool_workflow.data.samples import WFSample

from utils import (
    write_iimi_nucleotide_info,
    untar,
    load_and_format_prediction_results,
    write_all_otu_fasta,
)


@fixture
async def all_otu_index_path(work_path: Path):
    path = work_path / "all_otu_index"
    await asyncio.to_thread(path.mkdir)

    return path / "index"


@fixture
async def nucleotide_info_path(work_path: Path) -> Path:
    return work_path / "nucleotide_info.csv"


@fixture
async def output_path(work_path: Path) -> Path:
    path = work_path / "output"
    await asyncio.to_thread(path.mkdir)

    return path


@step(name="Build all-OTU index")
async def build_all_otu_index(
    all_otu_index_path: Path,
    index: WFIndex,
    logger,
    nucleotide_info_path: Path,
    proc: int,
    run_subprocess: RunSubprocess,
):
    """Map reads against all OTU sequences in the configured index."""

    all_otu_fasta_path = all_otu_index_path.parent / "all.fa"

    await asyncio.to_thread(
        write_iimi_nucleotide_info,
        index.path / "reference.json.gz",
        nucleotide_info_path,
        logger,
    )

    await asyncio.to_thread(
        write_all_otu_fasta,
        index.path / "reference.json.gz",
        all_otu_fasta_path,
    )

    await run_subprocess(
        [
            "bowtie2-build",
            "--threads",
            proc,
            "-f",
            all_otu_fasta_path,
            all_otu_index_path,
        ]
    )


@step(name="Map all OTUs")
async def map_all_otus(
    all_otu_index_path: Path,
    proc: int,
    run_subprocess: RunSubprocess,
    sample: WFSample,
    work_path: Path,
):
    sam_path = work_path / "mapped.sam"

    await run_subprocess(
        [
            "bowtie2",
            "--threads",
            proc,
            "-a",
            "-x",
            all_otu_index_path,
            "-U",
            *sample.read_paths,
            "-S",
            sam_path,
        ],
    )

    await run_subprocess(
        [
            "samtools",
            "view",
            "-bS",
            sam_path,
            "-o",
            work_path / "mapped.bam",
        ]
    )


@step
async def predict(
    analysis: WFAnalysis,
    bam_path: Path,
    index: WFIndex,
    logger,
    ml: WFMLModelRelease,
    nucleotide_info_path: Path,
    run_subprocess: RunSubprocess,
    output_path: Path,
    work_path: Path,
):
    """Run prediction using Iimi."""

    async def handler(line: bytes):
        logger.info("stdout", line=line.decode())

    assert ml.file_path.exists()
    assert ml.file_path.name == "model.tar.gz"

    await asyncio.to_thread(untar, ml.file_path, work_path)

    logger.info(
        "Running R script", model_contents=list((work_path / "model").iterdir())
    )

    await run_subprocess(
        [
            "Rscript",
            "./run.r",
            bam_path,
            work_path / "model" / "mappability_profile.rds",
            work_path / "model" / "trained_rf.rds",
            nucleotide_info_path,
            output_path,
            "--verbose",
        ],
        stdout_handler=handler,
    )

    result = await asyncio.to_thread(
        load_and_format_prediction_results,
        index.path / "reference.json.gz",
        output_path,
    )

    await analysis.upload_result(
        {"hits": sorted([h.dict() for h in result], key=itemgetter("name"))}
    )
