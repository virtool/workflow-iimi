import shutil
from pathlib import Path

from syrupy import SnapshotAssertion

from utils import write_iimi_nucleotide_info, load_and_format_prediction_results


def test_write_iimi_nucleotide_info(
    example_path: Path, logger, snapshot: SnapshotAssertion, work_path: Path
):
    shutil.copyfile(
        example_path / "index" / "reference.json.gz", work_path / "reference.json.gz"
    )

    output_path = work_path / "nucleotide_info.csv"

    write_iimi_nucleotide_info(work_path / "reference.json.gz", output_path, logger)

    assert output_path.read_text() == snapshot(name="file")


def test_load_prediction_results(snapshot: SnapshotAssertion, work_path: Path):
    output_path = work_path / "output"
    output_path.mkdir()

    shutil.copyfile("example/index/reference.json.gz", work_path / "reference.json.gz")

    for filename in (
        "coverage.csv",
        "untrustworthy.csv",
        "prediction_sequence.csv",
        "prediction_virus.csv",
    ):
        shutil.copyfile(f"example/{filename}", output_path / filename)

    assert [
        obj.dict()
        for obj in load_and_format_prediction_results(
            work_path / "reference.json.gz", output_path
        )
    ] == snapshot
