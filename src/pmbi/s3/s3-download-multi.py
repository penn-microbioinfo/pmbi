import argparse
import re
import os
import logging
import boto3
from pathos.multiprocessing import ProcessPool
from pathos.helpers import cpu_count

'''
def list_files(bucket, prefix):
    olist = boto3.client("s3").list_objects(Bucket = bucket, Prefix = prefix)
    if olist["IsTruncated"]:
        raise OSError("Response is truncated - you should use a paginator here instead.")

    if "Contents" in olist:
        return olist["Contents"]
    else:
        raise OSError(f"`{prefix}` does not match any objects in bucket `{bucket}`")
'''
def list_files(bucket, prefix):
    objnames = []
    pag = boto3.client("s3").get_paginator("list_objects")
    page_itr = pag.paginate(Bucket = bucket, Prefix = prefix)

    for page in page_itr:
        for i in page["Contents"]:
            yield i
            #objnames.append(i["Key"])
    #return objnames

def _get_file(key, bucket, chunksize = 1.024e6, do_not_overwrite = False):
    getobj = boto3.client("s3").get_object(Bucket = bucket, Key = key)
    base = os.path.basename(key)
    if os.path.isfile(base) and do_not_overwrite:
        pass
    else:
        with open(base, 'wb') as local:
            for chunk in getobj["Body"].iter_chunks():
                local.write(chunk)
    
def get_files(keylist, bucket, nproc, chunksize = 1.024e6, do_not_overwrite=False):
    with ProcessPool(nodes = nproc) as p:
        f_args = (
                [bucket]*len(keylist),
                [chunksize]*len(keylist),
                [do_not_overwrite]*len(keylist)
                )
        p.map(_get_file, keylist, *f_args)
        #p.map(_get_file, keylist, buckets, do_not_overwrite=do_not_overwrite)

if __name__ == "__main__":
    # Python > 3.8
    #logging.basicConfig(encoding='utf-8', level=logging.INFO)
    # Python <= 3.8
    logging.basicConfig(level=logging.INFO)
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bucket", action = "store", help = "")
    parser.add_argument("--do_not_overwrite", action = "store_true", help = "Don't overwrite files that exist. Default: overwrite")
    parser.add_argument("-k", "--prefix", action = "store", help = "Prefix for objects to download (i.e., the bucket directory).")
    parser.add_argument("-p", "--nproc", action = "store", default = cpu_count(), type = int, help = "")
    parser.add_argument("--chunksize",  action = "store", default = 1.024e6, type = float, help = "")
    parser.add_argument("--pattern", default = ".", required = False)
    args = parser.parse_args()
    
    p = re.compile(args.pattern)
    oblist = [x["Key"] for x in list_files(args.bucket, args.prefix) if re.search(p, x["Key"]) is not None]
    print(oblist)
    get_files(oblist, args.bucket, nproc=args.nproc, chunksize=args.chunksize, do_not_overwrite = args.do_not_overwrite)



