import subprocess
from pmbi.logging import streamLogger

def run_and_log(cmd, logger):
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
