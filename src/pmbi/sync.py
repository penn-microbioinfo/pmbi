from pmbi.logging import streamLogger
import shutil
from pathlib import Path
from typing import List, Union

logger = streamLogger(__name__)

class PathSync:
    """
    Synchronizes a specific list of files from one directory to another.
    """

    def __init__(
        self,
        src_dir: Union[str, Path],
        dest_dir: Union[str, Path],
        files_to_sync: List[Union[str,Path]],
    ):
        """
        Prepares to sync a list of files from a source to a destination.

        Args:
            src_dir: The source directory.
            dest_dir: The destination directory.
            files_to_sync: A list of *relative* paths of files to copy.

        Raises:
            NotADirectoryError: if src_dir does not exist
        """
        self.src_dir = Path(src_dir)
        self.dest_dir = Path(dest_dir)
        self.files_to_sync = [Path(p) for p in files_to_sync]

        if not self.src_dir.is_dir():
            mes = f"Source directory does not exist or is not a directory: {self.src_dir}"
            logger.critical(mes)
            raise NotADirectoryError(mes)

    def sync(self):
        """
        Executes the file synchronization.

        Copies the list of files from the source directory to the destination
        directory, preserving the relative path structure.

        Raises:
            FileNotFoundError: If a file in the sync list does not exist in the source.
        """
        logger.info(
            f"Syncing {len(self.files_to_sync)} files "
            f"from '{self.src_dir}' to '{self.dest_dir}'."
        )
        if not self.files_to_sync:
            logger.warning("files_to_sync is empty; nothing to do.")
            return

        for relpath in self.files_to_sync:
            src = self.src_dir / relpath
            dest = self.dest_dir / relpath

            if not src.exists():
                raise FileNotFoundError(
                    f"Path not found in src_dir: {relpath}"
                )

            # Create dest directory structure that doesn't exist
            dest.parent.mkdir(parents=True, exist_ok=True)

            # Sync it
            logger.debug(f"Copying '{src}' to '{dest}'")
            if src.is_file():
                shutil.copy(src, dest)
            elif src.is_dir():
                shutil.copytree(src, dest, dirs_exist_ok=True)
            else:
                raise OSError("Path is not file or dir: {src}")

