"""Test for workflow utilities."""

import shutil
from pathlib import Path

from syrupy import SnapshotAssertion

from utils import load_and_format_prediction_results


def test_load_prediction_results(snapshot: SnapshotAssertion, work_path: Path):
    """Test that the load_prediction_results function works as expected."""
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
            work_path / "reference.json.gz",
            output_path,
        )
    ] == snapshot
