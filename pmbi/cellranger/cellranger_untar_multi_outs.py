import tarfile
import logging
import argparse
import pathlib
import re
from cellranger_command import super_sample_p
import sys
from pathos.multiprocessing import ProcessPool
import pathos.helpers
class Pooled(object):
    def __init__ (self, func = None, nproc = pathos.helpers.cpu_count()):
        pass
def pool_it(func):
    def wrapper(*args, **kwargs):
        with ProcessPool(nodes = 2) as p:
            p.map(func, *args, **kwargs)
    return wrapper

@pool_it
def untar_cellranger_outs(fname):
    tf = tarfile.open(fname, 'r')
    ss = re.search(super_sample_p, fname).group(0)

    headdirs = {pathlib.PurePath(m.name).parts[0] for m in tf.getmembers()}
    if len(headdirs) == 1:
        ### Had some early archives that were made with the super_sample as the headdir, versus `outs`
        if headdirs.pop() == ss:
            extract_to = pathlib.Path('.')
        else:
            extract_to = pathlib.Path(ss)
            extract_to.mkdir()
        logging.warning(f"Extracting {tf}")    
        tf.extractall(path=extract_to)
        logging.warning(f"Done extracting {tf}")    
        pathlib.Path(fname).unlink()
    else:
        raise ValueError("More than one head directory in tar: {fname}")

    tf.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("archives", nargs='+', help="TAR archives to be untarred.")
    args = parser.parse_args()

    untar_cellranger_outs(args.archives)
