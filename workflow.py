import asyncio
from operator import itemgetter
from pathlib import Path

from pyfixtures import fixture
from virtool_workflow import RunSubprocess, step
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.ml import WFMLModelRelease
from virtool_workflow.data.samples import WFSample

from utils import (
    load_and_format_prediction_results,
    write_all_otu_fasta,
    write_iimi_nucleotide_info,
)


@fixture
async def all_otu_index_path(work_path: Path):
    path = work_path / "all_otu_index"
    await asyncio.to_thread(path.mkdir)

    return path / "index"


@fixture
async def output_path(work_path: Path) -> Path:
    path = work_path / "output"
    await asyncio.to_thread(path.mkdir)

    return path


@step(name="Build all-OTU index")
async def build_all_otu_index(
    all_otu_index_path: Path,
    ml: WFMLModelRelease,
    logger,
    proc: int,
    run_subprocess: RunSubprocess,
):
    """Map reads against all OTU sequences in the configured index."""

    logger.info("creating bowtie2 mapping index")

    await run_subprocess(
        [
            "bowtie2-build",
            "--threads",
            proc,
            "-f",
            ml.path / "reps.fa",
            all_otu_index_path,
        ],
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
        ],
    )


@step
async def predict(
    analysis: WFAnalysis,
    logger,
    ml: WFMLModelRelease,
    run_subprocess: RunSubprocess,
    output_path: Path,
    work_path: Path,
):
    """Run prediction using Iimi."""

    async def handler(line: bytes):
        logger.info("stdout", line=line.decode())

    logger.info("running rscript")

    print( [
            "Rscript",
            "./run.r",
            work_path / "mapped.bam",
            ml.path / "unreliable_regions.csv",
            ml.path / "trained_xgb.rds",
            ml.path / "sequence_info.csv",
            output_path,
            "--verbose",
        ])

    await run_subprocess(
        [
            "Rscript",
            "./run.r",
            work_path / "mapped.bam",
            ml.path / "unreliable_regions.csv",
            ml.path / "trained_xgb.rds",
            ml.path / "sequence_info.csv",
            output_path,
            "--verbose",
        ],
        stdout_handler=handler,
    )

    result = await asyncio.to_thread(
        load_and_format_prediction_results,
        ml.path / "reference.json.gz",
        output_path,
    )

    await analysis.upload_result(
        {"hits": sorted([h.dict() for h in result], key=itemgetter("name"))},
    )
