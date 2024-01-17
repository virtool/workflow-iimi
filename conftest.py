from pathlib import Path

import pytest
import arrow
from structlog import get_logger

pytest_plugins = [
    "virtool_workflow.pytest_plugin",
]


@pytest.fixture
def example_path():
    return Path(__file__).parent / "example"


@pytest.fixture
def logger():
    return get_logger("workflow")


@pytest.fixture
def proc(request):
    return request.config.getoption("--proc")


@pytest.fixture
def static_datetime():
    return arrow.get(2020, 1, 1, 1, 1, 1).naive


@pytest.fixture
def work_path(tmpdir) -> Path:
    path = Path(tmpdir) / "work"
    path.mkdir()

    return path


# Add an option to configure a process count.
# This is used to test the parallelism of the workflow.
def pytest_addoption(parser):
    parser.addoption(
        "--proc",
        action="store",
        default=1,
        help="Number of processes allow for workflow tools",
    )
