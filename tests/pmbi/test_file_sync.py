import pytest
from pathlib import Path

from pmbi.sync import PathSync


@pytest.fixture
def syncdirs(tmp_path):
    """
    Pytest fixture to set up temporary source and destination directories
    with dummy files for testing PathSync.
    """
    from_dir = tmp_path / "source_dir"
    to_dir = tmp_path / "dest_dir"

    from_dir.mkdir()
    to_dir.mkdir()

    # Create dummy files in from_dir
    (from_dir / "file1.txt").write_text("content of file1")
    (from_dir / "file2.csv").write_text("a,b,c\n1,2,3")
    (from_dir / "data").mkdir()
    (from_dir / "data" / "nested_file.json").write_text('{"key": "value"}')
    (from_dir / "data" / "sub_data").mkdir()
    (from_dir / "data" / "sub_data" / "deep_file.log").write_text("log message")

    return {
        "from_dir": from_dir,
        "to_dir": to_dir,
        "files_created": [
            "file1.txt",
            "file2.csv",
            "data/nested_file.json",
            "data/sub_data/deep_file.log",
        ],
    }

def test_sync_dir(syncdirs):
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], ["data"])
    ps.sync()
    expected_relpaths = ["data/nested_file.json", "data/sub_data/deep_file.log"]
    for relpath in expected_relpaths:
        fullpath=syncdirs["to_dir"] / relpath
        assert fullpath.is_file()
        assert fullpath.read_text() == (syncdirs["from_dir"] / fullpath).read_text()

def test_init_valid_dirs(syncdirs):
    """Test that PathSync initializes correctly with valid directories."""
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], ["file1.txt"])
    assert ps.src_dir == syncdirs["from_dir"]
    assert ps.dest_dir == syncdirs["to_dir"]
    assert ps.files_to_sync == [Path("file1.txt")]

def test_init_non_existent_from_dir(tmp_path):
    """Test that __init__ raises NotADirectoryError if from_dir does not exist."""
    non_existent_dir = tmp_path / "non_existent"
    with pytest.raises(NotADirectoryError, match="Source directory does not exist"):
        PathSync(non_existent_dir, tmp_path / "dest", ["some_file.txt"])

def test_sync_empty_list(syncdirs):
    """Test that sync() does nothing and doesn't raise errors for an empty files_to_sync list."""
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], [])
    ps.sync()
    # Assert that dest_dir is still empty (except for its own creation)
    assert not list(syncdirs["to_dir"].iterdir())

def test_sync_basic_file(syncdirs):
    """Test syncing a single file at the root of src_dir."""
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], ["file1.txt"])
    ps.sync()

    expected_dest_path = syncdirs["to_dir"] / "file1.txt"
    assert expected_dest_path.is_file()
    assert expected_dest_path.read_text() == (syncdirs["from_dir"] / "file1.txt").read_text()

def test_sync_nested_file(syncdirs):
    """Test syncing a file located in a subdirectory of src_dir."""
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], ["data/nested_file.json"])
    ps.sync()

    expected_dest_path = syncdirs["to_dir"] / "data" / "nested_file.json"
    assert expected_dest_path.is_file()
    assert expected_dest_path.parent.is_dir()  # Check if data/ was created
    assert expected_dest_path.read_text() == (
        syncdirs["from_dir"] / "data" / "nested_file.json"
    ).read_text()

def test_sync_multi_files(syncdirs):
    """Test syncing multiple files, including root and nested ones."""
    files_to_sync = ["file1.txt", "data/sub_data/deep_file.log", "file2.csv"]
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], files_to_sync)
    ps.sync()

    for rel_path in files_to_sync:
        expected_dest_path = syncdirs["to_dir"] / rel_path
        assert expected_dest_path.is_file()
        assert expected_dest_path.read_text() == (syncdirs["from_dir"] / rel_path).read_text()

def test_sync_file_not_found(syncdirs):
    """Test that sync() raises FileNotFoundError if a file in files_to_sync does not exist."""
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], ["non_existent.txt"])
    with pytest.raises(FileNotFoundError, match=f"Path not found in src_dir: non_existent.txt"):
        ps.sync()

def test_sync_overwrites_existing_file(syncdirs):
    """Test that sync() correctly overwrites a file if it already exists in dest_dir."""
    file_to_sync = "file1.txt"
    ps = PathSync(syncdirs["from_dir"], syncdirs["to_dir"], [file_to_sync])
    
    # First sync
    ps.sync()
    original_content = (syncdirs["from_dir"] / file_to_sync).read_text()
    (syncdirs["to_dir"] / file_to_sync).write_text("old content") # Manually change content in dest

    # Modify source file content
    new_content = "updated content of file1"
    (syncdirs["from_dir"] / file_to_sync).write_text(new_content)

    # Second sync should overwrite
    ps.sync()

    assert (syncdirs["to_dir"] / file_to_sync).read_text() == new_content

