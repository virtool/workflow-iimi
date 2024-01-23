import json
import os
import shutil
from pathlib import Path

import pytest
from structlog import get_logger
from syrupy import SnapshotAssertion
from virtool_core.models.analysis import AnalysisSample
from virtool_core.models.enums import LibraryType
from virtool_core.models.index import IndexNested
from virtool_core.models.reference import ReferenceNested
from virtool_core.models.subtraction import SubtractionNested
from virtool_core.utils import decompress_file
from virtool_workflow import RunSubprocess
from virtool_workflow.data.analyses import WFAnalysis
from virtool_workflow.data.indexes import WFIndex
from virtool_workflow.data.ml import WFMLModelRelease
from virtool_workflow.data.samples import WFSample
from virtool_workflow.pytest_plugin.data import Data
from virtool_workflow.utils import untar

from workflow import build_all_otu_index, map_all_otus, predict


@pytest.fixture
def analysis(data: Data, mocker) -> WFAnalysis:
    data.analysis.ready = False

    analysis = mocker.Mock(spec=WFAnalysis)

    analysis.id = data.analysis.id
    analysis.index = IndexNested.parse_obj(data.index)
    analysis.ml = data.ml
    analysis.reference = ReferenceNested.parse_obj(data.index.reference)
    analysis.sample = AnalysisSample.parse_obj(data.sample)
    analysis.subtractions = [SubtractionNested.parse_obj(data.subtraction)]
    analysis.workflow = data.analysis.workflow

    return analysis


@pytest.fixture
def index(data: Data, example_path: Path, work_path: Path) -> WFIndex:
    index_path = work_path / "indexes" / data.index.id

    shutil.copytree(example_path / "index", index_path)

    decompress_file(index_path / "otus.json.gz", index_path / "otus.json")

    index_path / "otus.json.gz"

    with open(example_path / "index" / "sequence_otu_map.json") as f:
        sequence_otu_map = json.load(f)

    return WFIndex(
        id=data.index.id,
        path=index_path,
        manifest={},
        reference=data.index.reference,
        sequence_lengths={},
        sequence_otu_map=sequence_otu_map,
    )


@pytest.fixture
def ml(data: Data, example_path: Path, work_path: Path) -> WFMLModelRelease:
    path = work_path / "ml" / str(data.ml.id) / str(data.ml.model.id)
    path.mkdir(parents=True)

    obj = WFMLModelRelease(id=data.ml.id, name=data.ml.name, path=path)

    shutil.copyfile(example_path / "model.tar.gz", obj.file_path)
    untar(obj.file_path, path)

    for filepath in (path / "model").iterdir():
        os.rename(filepath, path / filepath.name)

    logger = get_logger("ml")
    logger.info("Extracted ML model", contents=list(path.iterdir()))

    return obj


@pytest.fixture
def sample(data: Data, example_path: Path, work_path: Path) -> WFSample:
    sample_path = work_path / "samples" / data.sample.id
    sample_path.mkdir(parents=True)

    shutil.copyfile(
        example_path / "sample" / "large.fq.gz",
        sample_path / "reads_1.fq.gz",
    )

    return WFSample(
        id=data.sample.id,
        library_type=LibraryType.normal,
        name=data.sample.name,
        paired=False,
        quality=data.sample.quality,
        read_paths=(sample_path / "reads_1.fq.gz",),
    )


@pytest.mark.asyncio
async def test_workflow(
    analysis: WFAnalysis,
    example_path: Path,
    index: WFIndex,
    logger,
    ml: WFMLModelRelease,
    proc: int,
    run_subprocess: RunSubprocess,
    sample: WFSample,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    """Make sure the workflow runs without error."""

    all_otu_index_path = work_path / "all_otu_index"
    all_otu_index_path.mkdir()

    nucleotide_info_path = work_path / "nucleotide_info.csv"

    # Test Step: Build all-OTU index
    await build_all_otu_index(
        all_otu_index_path / "all",
        index,
        logger,
        nucleotide_info_path,
        proc,
        run_subprocess,
    )

    assert set([path.name for path in all_otu_index_path.iterdir()]) == {
        "all.1.bt2",
        "all.2.bt2",
        "all.3.bt2",
        "all.4.bt2",
        "all.fa",
        "all.rev.1.bt2",
        "all.rev.2.bt2",
    }

    # Test Step: Map all OTUs
    await map_all_otus(
        work_path / "all_otu_index" / "all",
        proc,
        run_subprocess,
        sample,
        work_path,
    )

    output_path = work_path / "output"
    output_path.mkdir()

    await predict(
        analysis,
        logger,
        ml,
        nucleotide_info_path,
        run_subprocess,
        output_path,
        work_path,
    )

    assert analysis.upload_result.call_args[0] == snapshot


@pytest.mark.asyncio
async def test_predict(
    analysis: WFAnalysis,
    example_path: Path,
    index: WFIndex,
    logger,
    ml: WFMLModelRelease,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    bam_path = work_path / "mapped.bam"
    shutil.copyfile(example_path / "mapped.bam", bam_path)

    nucleotide_info_path = work_path / "nucleotide_info.csv"
    shutil.copyfile(example_path / "nucleotide_info.csv", nucleotide_info_path)

    output_path = work_path / "output"
    output_path.mkdir()

    await predict(
        analysis,
        logger,
        ml,
        nucleotide_info_path,
        run_subprocess,
        output_path,
        work_path,
    )

    assert analysis.upload_result.call_args[0] == snapshot
