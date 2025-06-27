"""Utilities for handling subprocess output streaming"""

import subprocess
from typing import Iterator, Optional
from io import StringIO

def Popen_stream(
    cmd: list[str],
    **kwargs
) -> Iterator[str]:
    """
    Yeild subprocess.Popen stdout/stderr 
    
    Args:
        cmd: Command list to execute
        **kwargs: Additional kwargs passed to subprocess.Popen
        
    Yields:
        Output lines as they become available
    """
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        **kwargs
    ) as proc:
        
        while True:
            line = proc.stdout.readline()
            if not line and proc.poll() is not None:
                break
                
            if line:
                yield line.rstrip('\n')
