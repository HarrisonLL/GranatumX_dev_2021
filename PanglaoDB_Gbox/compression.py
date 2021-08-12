import json
import gzip
from io import StringIO, BytesIO
import base64



def compress_chunk(param1):
    param1["matrix"] = param1["matrix"].tolist()
    tmp = json.dumps(param1)
    bio = BytesIO()
    bio.write(tmp.encode('utf-8'))
    bio.seek(0)
    stream = BytesIO()
    compressor = gzip.GzipFile(fileobj=stream, mode = 'w')
    while True:
        chunk = bio.read(8192)
        if not chunk:
            compressor.close()
            break
        compressor.write(chunk)
    encoded = base64.b64encode(stream.getvalue())
    return(encoded.decode('utf-8'))
