import boto3
import time
import uuid
import io
import botocore
import argparse
import subprocess
import os
import logging
import re
import sys
import hashlib
import base64
from pmbi.s3.lib import split_s3_uri
from pathos.multiprocessing import ProcessPool
from pathos.helpers import cpu_count

PARTS_FNAME_PATTERN="^x[a-z]{2}$"

def pooled(f, ):
    f()
    
class Part(object):
    def __init__(self, partnumber, part_fname, part_digest = None, part_hexdigest = None):
        self.number = partnumber
        self.fname = part_fname
        self.digest = part_digest
        self.digest = part_hexdigest
    def _etag(self):
        if self.digest is not None:
            return md5hash(io.BytesIO(self.digest)).hexdigest()
            #self.etag =  f"{md5hash( io.BytesIO(b''.join(self.part_digests)) ).hexdigest()}-{len(self.part_digests)}"
        else:
            raise ValueError("No part digest to make etag out of.")

class S3MultiPartUpload(object):
    def __init__(self, bucket, key, abort_on_failure_to_complete = False, nproc = 1):
        self.client = None 
        self.bucket = bucket
        self.key = key
        self.upload_id = None
        #self.part_fnames = None
        self.parts = None
        self.nparts = None
        #self.part_digests = None
        self.etag = None
        self.abort_on_failure_to_complete = abort_on_failure_to_complete
        self.nproc = nproc

        self.upload_id = self._create_upload()

    def _create_upload(self):
        return boto3.client("s3").create_multipart_upload(Bucket = self.bucket, Key = self.key)["UploadId"]
    
    def _upload_part(self, partnumber, partpath, parthash):
        start = time.time()
        with open(partpath, 'rb') as fopen:
            try:
                boto3.client("s3").upload_part(Bucket = self.bucket, Key = self.key, UploadId = self.upload_id, PartNumber = partnumber, Body = fopen, ContentMD5 = parthash)
            except botocore.exceptions.ClientError:
                self._abort()
                raise
        logging.info(f"Uploaded part #{partnumber}: {partpath} - {time.time()-start} seconds")

    def _part_fnames(self):
        return [p.fname for p in self.parts]
    
    def _part_numbers(self):
        return [p.number for p in self.parts]

    def _part_digests(self, as_bytes = True):
        if as_bytes:
            return [p.digest for p in self.parts]
        else:
            return [p.hexdigest for p in self.parts]

    def upload_parts(self, part_fnames):
        self.nparts = len(part_fnames)
        self.parts = [Part(pn, pf) for pn,pf in zip(range(1,self.nparts+1), part_fnames)]
        if self.nparts < self.nproc:
            logging.info(f"Only using as many processes as there are parts: {self.nparts}")
            self.nproc = self.nparts
        logging.info(f"Part filenames (IN THIS ORDER): {self._part_fnames()}")
        self.part_digests = S3MultiPartUpload._md5_digest_parts(self.parts, self.nproc)

        logging.info(f"Started uploading {self.nparts} parts.")
        start = time.time()
        print(self._part_numbers(), self._part_fnames())
        with ProcessPool(nodes=self.nproc) as p:
            p.map(self._upload_part, self._part_numbers(), self._part_fnames(), S3MultiPartUpload._encode_b64([p.digest for p in self.parts]))
   
        logging.info(f"Total upload time ({self.nparts} parts): {time.time()-start} seconds")

        self.delete_part_files()

    def delete_part_files(self):
        for pf in self._part_fnames():
            os.remove(pf)
        logging.info(f"Removed {len(self.parts)} part files: {self._part_fnames()})")

    def _etag():
        pass

    def _generate_final_etag (self):
        if self.parts is not None and all([p.digest is not None for p in self.parts]):
            self.etag =  f"{md5hash( io.BytesIO(b''.join(self._part_digests())) ).hexdigest()}-{len(self._part_digests())}"
        else:
            raise ValueError("Cannot generate ETag without parts.")
    
    def _matches_etag(self, other_etag):
        if self.etag == other_etag:
            return True
        else:
            return False

    def _encode_b64(part_digests):
        return [base64.b64encode(x).decode("utf-8") for x in part_digests]
    
    def _md5_digest_parts(parts, local_nproc = 1):
        def both_digests(x):
            x_hash = md5hash(x)
            x_digest = x_hash.digest()
            x_hexdigest = x_hash.hexdigest()
            return (x_digest, x_hexdigest)

        part_fnames = [p.fname for p in parts]

        logging.info(f"Generating md5 digests for {len(part_fnames)} file parts.")
        streams = [open(x, 'rb') for x in [p.fname for p in parts]]
        with ProcessPool(local_nproc) as p:
            digests = p.map(lambda x: both_digests(x), streams)
        _close_streams = [x.close() for x in streams]
        for i,p in enumerate(parts):
            p.digest = digests[i][0]
            p.hexdigest = digests[i][1]
        return [d[0] for d in digests]
        
        #return [md5hash(open(x, 'rb')).digest() for x in part_fnames]

    def list_parts(self):
        try:
            r = boto3.client("s3").list_parts(Bucket = self.bucket, Key = self.key, UploadId = self.upload_id)
        except botocore.exceptions.ClientError:
            self._abort()
            raise
        return r["Parts"]

    def _uploaded_parts_in_correct_order(self):
        # Sort remote parts by part number
        remote_parts = {p["PartNumber"]: p["ETag"].replace('"', '') for p in self.list_parts()}
        remote_parts = {i: remote_parts[i] for i in sorted(remote_parts.keys())}
        
        # Sort local parts by filename
        local_parts = {p.fname: (p.number, p.hexdigest) for p in self.parts} 
        local_parts = {i: local_parts[i] for i in sorted(local_parts.keys())} 
      
        # To make sure that the remote parts are numbered correctly, compare the hexdigests
        # from each list (sorted differently)
        local_fnames = list(local_parts.keys())
        remote_etags = list(remote_parts.values())
        local_etags = [x[1] for x in local_parts.values()]
        for idx,local in enumerate(local_etags):
            logging.info(f"{local_fnames[idx]}\t{local}\t{remote_etags[idx]}")
            print(local_fnames[idx])
            print(local, remote_etags[idx])
            if local != remote_etags[idx]:
                pass
                #return False 
        return False
        #return True

    def complete_upload(self):
        
        if not self._uploaded_parts_in_correct_order():
            raise ValueError("Local and remote parts mismatched.")

        self._generate_final_etag()

        #check list_parts against stored parts()

        multiparts = {"Parts": [{"PartNumber": p["PartNumber"], "ETag": p["ETag"]} for p in self.list_parts()]}

        try:
            r = boto3.client("s3").complete_multipart_upload(Bucket = self.bucket, Key = self.key, UploadId = self.upload_id, MultipartUpload=multiparts) 
        except:
            self._abort()
        remote_etag = r["ETag"].replace('"', '') 
        if not self._matches_etag(remote_etag):
            logging.critical(f"Local and remote ETags do not match: {self.etag}")
            if self.abort_on_failure_to_complete:
                self._abort()
            else:
                logging.critical(f"Unable to complete MultipartUpload. It will have to be done manually:\n{'-'*25}\naws s3api complete-multipart-upload --bucket {self.bucket} --key {self.key} --upload-id {self.upload_id} --multipart-upload {multiparts}")

    def _abort(self):
        if self.upload_id is None:
            logging.critical("Unable to abort prior to creation.")
            raise ValueError
        try:
            boto3.client("s3").abort_multipart_upload(Bucket = self.bucket, Key = self.key, UploadId = self.upload_id)
            logger.critical("Upload aborted successfully.")
        except botocore.exceptions.ClientError:
            logging.critical(f"Unable to abort MultipartUpload. It will have to be done manually:\n{'-'*25}\naws s3api abort-multipart-upload --bucket {self.bucket} --key {self.key} --upload-id {self.upload_id}")
            raise

    def upload_md5(self, md5_key, md5_bytes):
        assert os.path.dirname(md5_key) == os.path.dirname(self.key), "Dirnames for largefile key and md5 do not match"
        assert f"{os.path.basename(self.key)}.md5" == os.path.basename(md5_key), "MD5 base key does not match the large file name + .md5"
        boto3.client("s3").put_object(Bucket = self.bucket, Key = md5_key, Body = md5_bytes)

def get_upload_id(response):
    return response["UploadId"]

def make_file_parts(largefile, part_size = 4950000000, parts_outdir = None): 
    cmd = ["split", "-b", str(part_size), largefile] 
    p = subprocess.run(cmd, capture_output = True)
    try:
        p.check_returncode()
    except subprocess.CalledProcessError:
        logging.critical(f"Call to `split` failed with: {p.stderr}")
        sys.exit(1)

    pat = re.compile(PARTS_FNAME_PATTERN)
    part_fnames = [x for x in os.listdir(".") if re.match(pat, x)]
    
    # This sort is very important, becasue python reads the file names out
    # of sort order
    part_fnames.sort()
    
    check_file_part_sizes(largefile, part_fnames)

    return part_fnames

def check_file_part_sizes(largefile, part_fnames):
    file_part_sizes = [os.stat(x).st_size for x in part_fnames]
    file_parts_total_size = sum(file_part_sizes)
    largefile_size = os.stat(largefile).st_size
    if file_parts_total_size != largefile_size:
        logging.critical(f"Sum of file part sizes is not equal to file size. {file_parts_total_size} != {largefile_size}")
        sys.exit(1)

def md5hash(f, block_size = 2**20):
    md5 = hashlib.new("md5")
    while True:
        block = f.read(block_size)
        if not block:
            break
        md5.update(block)
    return md5


def check_final_etag(local_etag, remote_etag):
    if local_etag != remote_etag:
        logging.critical(f"Local ETag does not equal remote ETag: {local_etag} != {remote_etag}")
        sys.exit(1)

def main(uri: str,
         largefile: str,
         partsize: int = 4950000000,
         nproc: int = cpu_count()
         ):

    bucket, key = split_s3_uri(uri)

    logging.basicConfig(level=logging.INFO)

    logging.info(f"Uploading with up to {nproc} processes (default = detected {cpu_count()} CPUs))")

    logging.info(f"Generating MD5 for {largefile} in background.")
    logging.info(f"Splitting {largefile} ({float(os.stat(largefile).st_size)/float(1e9)} Gb) into parts.")
   
    largefile_path = os.path.realpath(largefile)
    unique_dirname = str(uuid.uuid4())
    logging.info(f"Creating temporary directory to populate with file parts: {unique_dirname}")
    os.mkdir(unique_dirname)
    os.chdir(unique_dirname)

    md5p = subprocess.Popen(["md5sum", largefile_path], stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    part_fnames = make_file_parts(largefile_path, part_size = partsize)
    md5_bytes,err = md5p.communicate()
    md5_key = f"{key}.md5"
    md5_basename = f"{os.path.basename(md5_key)}"
    if md5p.returncode != 0:
        raise OSError(err)
    else:
        with open(md5_basename, 'wb') as md5:
            md5.write(md5_bytes)

    # Begin interactions with S3
    mpu = S3MultiPartUpload(bucket, key, nproc = nproc)
    mpu.upload_parts(part_fnames)
    mpu.upload_md5(md5_key, md5_bytes)
    mpu.complete_upload()

    logging.shutdown()


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-u", "--uri", action = "store", required = True, type = str, help = "S3 URI to upload large file as.")
    parser.add_argument("-f", "--largefile", action = "store", type = str, help = "Path to large file to upload.")
    parser.add_argument("-s", "--partsize", action = "store", required = True, type = int, help = "", default = 4950000000)
    parser.add_argument("-p", "--nproc", action = "store", default = cpu_count(), type = int, help = "")
    args = parser.parse_args()

    main(**vars(args))
 
    # On AWS, a .nfsXXXX file is created that doesn't get removed until
    # after script closes, so can't do anything about it here - yet...
    # Just call script from a different directory and then you can delete 
    # all the uuid-named directories
    #os.chdir("../")
    #os.rmdir(unique_dirname)

