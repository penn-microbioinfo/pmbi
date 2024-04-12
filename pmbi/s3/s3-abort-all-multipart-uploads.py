import boto3
import json
import sys

bucket = sys.argv[1]
client = boto3.client("s3")
active_uploads = client.list_multipart_uploads(Bucket = bucket)
for upload in active_uploads["Uploads"]:
   client.abort_multipart_upload(Bucket = bucket, Key = upload["Key"], UploadId = upload["UploadId"]) 

