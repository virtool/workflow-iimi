"""Tests for the Iimi workflow."""

import json
import shutil
from pathlib import Path

import pytest
from pytest_mock import MockerFixture
from structlog import get_logger
from syrupy import SnapshotAssertion
from virtool.analyses.models import AnalysisSample
from virtool.indexes.models import IndexNested
from virtool.models.enums import LibraryType
from virtool.references.models import ReferenceNested
from virtool.subtractions.models import SubtractionNested
from virtool.utils import decompress_file, decompress_tgz
from virtool.workflow import RunSubprocess
from virtool.workflow.data.analyses import WFAnalysis
from virtool.workflow.data.indexes import WFIndex
from virtool.workflow.data.ml import WFMLModelRelease
from virtool.workflow.data.samples import WFSample
from virtool.workflow.pytest_plugin import WorkflowData

from workflow import build_all_otu_index, map_all_otus, predict


@pytest.fixture
def analysis(workflow_data: WorkflowData, mocker: MockerFixture) -> WFAnalysis:
    """A fake analysis for testing."""
    workflow_data.analysis.ready = False

    analysis = mocker.Mock(spec=WFAnalysis)

    analysis.id = workflow_data.analysis.id
    analysis.index = IndexNested.parse_obj(workflow_data.index)
    analysis.ml = workflow_data.ml
    analysis.reference = ReferenceNested.parse_obj(workflow_data.index.reference)
    analysis.sample = AnalysisSample.parse_obj(workflow_data.sample)
    analysis.subtractions = [SubtractionNested.parse_obj(workflow_data.subtraction)]
    analysis.workflow = workflow_data.analysis.workflow

    return analysis


@pytest.fixture
def index(workflow_data: WorkflowData, example_path: Path, work_path: Path) -> WFIndex:
    """A fake index for testing."""
    index_path = work_path / "indexes" / workflow_data.index.id

    shutil.copytree(example_path / "index", index_path)

    decompress_file(index_path / "otus.json.gz", index_path / "otus.json")

    with (example_path / "index" / "sequence_otu_map.json").open() as f:
        sequence_otu_map = json.load(f)

    return WFIndex(
        id=workflow_data.index.id,
        path=index_path,
        manifest={},
        reference=workflow_data.index.reference,
        sequence_lengths={},
        sequence_otu_map=sequence_otu_map,
    )


@pytest.fixture
def ml(
    example_path: Path,
    work_path: Path,
    workflow_data: WorkflowData,
) -> WFMLModelRelease:
    path = work_path / "ml" / str(workflow_data.ml.id) / str(workflow_data.ml.model.id)
    path.mkdir(parents=True)

    obj = WFMLModelRelease(
        id=workflow_data.ml.id,
        name=workflow_data.ml.name,
        path=path,
    )

    shutil.copyfile(example_path / "model.tar.gz", obj.file_path)
    decompress_tgz(obj.file_path, path)

    for filepath in (path / "model").iterdir():
        filepath.rename(path / filepath.name)

    return obj


@pytest.fixture
def sample(
    example_path: Path,
    work_path: Path,
    workflow_data: WorkflowData,
) -> WFSample:
    sample_path = work_path / "samples" / workflow_data.sample.id
    sample_path.mkdir(parents=True)

    shutil.copyfile(
        example_path / "sample" / "large.fq.gz",
        sample_path / "reads_1.fq.gz",
    )

    return WFSample(
        id=workflow_data.sample.id,
        library_type=LibraryType.normal,
        name=workflow_data.sample.name,
        paired=False,
        quality=workflow_data.sample.quality,
        read_paths=(sample_path / "reads_1.fq.gz",),
    )


@pytest.mark.asyncio
async def test_workflow(
    analysis: WFAnalysis,
    index: WFIndex,
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
        get_logger("test_workflow"),
        nucleotide_info_path,
        proc,
        run_subprocess,
    )

    assert {path.name for path in all_otu_index_path.iterdir()} == {
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
        ml,
        nucleotide_info_path,
        run_subprocess,
        output_path,
        work_path,
    )

    assert analysis.upload_result.call_args[0] == snapshot


async def test_predict(
    analysis: WFAnalysis,
    example_path: Path,
    ml: WFMLModelRelease,
    run_subprocess: RunSubprocess,
    snapshot: SnapshotAssertion,
    work_path: Path,
):
    """Test that the prediction step works as expected."""
    bam_path = work_path / "mapped.bam"
    shutil.copyfile(example_path / "mapped.bam", bam_path)

    nucleotide_info_path = work_path / "nucleotide_info.csv"
    shutil.copyfile(example_path / "nucleotide_info.csv", nucleotide_info_path)

    output_path = work_path / "output"
    output_path.mkdir()

    await predict(
        analysis,
        get_logger("test_predict"),
        ml,
        nucleotide_info_path,
        run_subprocess,
        output_path,
        work_path,
    )

    assert analysis.upload_result.call_args[0] == snapshot
