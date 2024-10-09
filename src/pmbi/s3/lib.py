import boto3
import re

# %%
def split_s3_uri(uri: str) -> tuple[str, str]:
    uri_pattern = re.compile("^s3[:][/]{2}([^/]+)[/](.*$)")
    s = re.match(uri_pattern, uri)
    if s is not None:
        Bucket = s.group(1)
        Prefix = s.group(2)
        return Bucket, Prefix
    else:
        raise ValueError(f"Invalid URI provided: {uri}")

# %%
def object_key_list(uri: str) -> list[str]:
    Bucket,Prefix = split_s3_uri(uri)
    s3 = boto3.resource("s3")
    bucket = s3.Bucket(Bucket)
    return [o.key for o in bucket.objects.filter(Prefix = Prefix)]
