from pathlib import Path

import pytest


@pytest.fixture
def proc(request: pytest.FixtureRequest) -> int:
    """An arbitrary `proc` config value for testing."""
    return int(request.config.getoption("--proc"))


@pytest.fixture
def work_path(tmp_path: Path) -> Path:
    path = tmp_path / "work"
    path.mkdir()

    return path


# Add an option to configure a process count.
# This is used to test the parallelism of the workflow.
def pytest_addoption(parser: pytest.Parser) -> None:
    """Add command line options for pytest."""
    parser.addoption(
        "--proc",
        action="store",
        default=1,
        help="Number of processes allow for workflow tools",
    )
