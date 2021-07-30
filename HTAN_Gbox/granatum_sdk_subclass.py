import granatum_sdk
import base64
from io import StringIO, BytesIO
import gzip

class granatum_extended(granatum_sdk):

    """
    chunk data model:
    {   "origin data size": [#, #]
        "current chunk": ["last_gbox_name", "col/row"],
        "suggested chunk": {
                "deepimpute": ["col", ##],
                "log-transform":["row", ##],
                ..
        },
        "chunk1": base64string,
        "chunk2": base64string,
        ...
    }

    """

    def __init__(self, gbox_name):
        self.gbox_name = gbox_name


    def check_chunk(self, json_dict):
        return "current chunk" in json_dict


    def decode_json(jsonstr):
        bio = BytesIO()
        stream = BytesIO(base64.b64decode(jsonstr.encode('utf-8')))
        decompressor = gzip.GzipFile(fileobj=stream, mode='r')

        while True:
            chunk = decompressor.read(8192)
            if not chunk:
                decompressor.close()
                bio.seek(0)
                output = bio.read().decode("utf-8")
                break
            bio.write(chunk)
        return json.loads(output)


    def adjust_transform(self, json_dict):
        # user rerun the gbox
        if json_dict["current chunk"][0] == self.gbox_name:
            return
        # four situations
        if json_dict["current chunk"][-1] == json_dict["suggested chunk"].get(self.gbox_name)[0] == "col":
            pass
        elif json_dict["current chunk"][-1] == json_dict["suggested chunk"].get(self.gbox_name)[0] == "row":
            pass
        elif json_dict["current chunk"][-1] == "col" and json_dict["suggested chunk"].get(self.gbox_name)[0] == "row":
            for i in range(len(json_dict)-3):
                # decompress
            pass
        elif json_dict["current chunk"][-1] == "row" and json_dict["suggested chunk"].get(self.gbox_name)[0] == "col":
            pass


    def iterate_chunks(self, json_dict):
        for k,v in json_dict.items():
            tmp = base64.b64decode(v)
            bio = BytesIO()
            stream = BytesIO(base64.b64decode(b.encode('utf-8')))
            decompressor = gzip.GzipFile(fileobj=stream, mode='r')
            while True:  # until EOF
                chunk = decompressor.read(8192)
                if not chunk:
                    decompressor.close()
                    bio.seek(0)
                    output = bio.read().decode("utf-8")
                    break
                bio.write(chunk)
            # yield an assay type
            yield output
