import subprocess
from io import StringIO
from pmbi.logging import streamLogger

def run_and_log(cmd, logger = streamLogger(__name__)):
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,  # Affirms that proc.stdout will not be None
        stderr=subprocess.STDOUT,
        text=True,
    ) as proc:
        while True:
            line = proc.stdout.readline()
            if not line and proc.poll() is not None:
                break

            if line:
                logger.info(line.strip())

        if proc.returncode != 0:
            mes = f"Process failed with returncode: {proc.returncode}"
            logger.critical(mes)
            raise subprocess.SubprocessError(mes)

def run_and_capture(cmd, logger = streamLogger(__name__)) -> StringIO:
    with subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,  # Affirms that proc.stdout will not be None
        stderr=None,
        text=True,
    ) as proc:
        buf = StringIO()
        while True:
            line = proc.stdout.readline()
            if not line and proc.poll() is not None:
                break

            if line:
                buf.write(f"{line}")

        if proc.returncode != 0:
            mes = f"Process failed with returncode: {proc.returncode}"
            logger.critical(mes)
            raise subprocess.SubprocessError(mes)

        buf.seek(0)
        return buf
