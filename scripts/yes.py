from joblib import Parallel, delayed
import subprocess

def _yes():
    x=0
    while True:
        if x == 1e12:
            x = 0
        x = x+1
    # p = subprocess.call(["yes"], stdout = subprocess.PIPE, stderr=subprocess.PIPE)
    # out,err = p.communicate()

n_jobs = 354

out = Parallel(n_jobs=n_jobs, timeout=1200)(
    delayed(_yes)() for f in list(range(0,n_jobs))
)

