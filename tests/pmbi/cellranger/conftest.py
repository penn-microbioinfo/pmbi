import pytest
from pathlib import Path

def pytest_addoption(parser):
    parser.addoption(
        "--cellranger_atac_path",
        action="store",
        required=True,
        help="Path to cellranger-atac package directory"
    )

@pytest.fixture
def runopts(request):
    pkg = Path(request.config.getoption("--cellranger_atac_path"))
    return {
        "binary": pkg.joinpath("bin", "cellranger-atac"),
        "id": "test",
        "reference": pkg.joinpath("external/arc_testrun_files/reference"),
        "fastqs": pkg.joinpath("external/cellranger_atac_tiny_fastq/1.0.0"),
        "sample": "testfastq_ATAC",
    }
