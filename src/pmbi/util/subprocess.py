"""Utilities for handling subprocess output streaming"""

import subprocess
from typing import Iterator, Optional
from io import StringIO

def stream_output(
    cmd: list[str],
    buffer: Optional[StringIO] = None,
    **kwargs
) -> Iterator[str]:
    """
    Stream subprocess output to both a buffer and yield lines.
    
    Args:
        cmd: Command list to execute
        buffer: Optional StringIO buffer to capture output
        **kwargs: Additional kwargs passed to subprocess.Popen
        
    Yields:
        Output lines as they become available
    """
    if buffer is None:
        buffer = StringIO()
        
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
                buffer.write(line)
                yield line.rstrip('\n')
