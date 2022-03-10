import numpy as np
import pandas as pd
import os
import shutil
from tqdm import tqdm
import string
from scipy import sparse
from scipy.sparse import rand
import random
import gc
import json
from io import StringIO, BytesIO
import gzip
import base64


def compress_assay(exported_assay):
    tmp = json.dumps(exported_assay)
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


def main():
     # destroy if exits, and then create ./datasets directory
    dirpath = './datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    cell_sizes = [100, 500,1000, 2000, 4000,6000,8000]
    gene_sizes = [int(1e4), int(2e4), int(3e4), int(4e4), int(5e4),int(6e4)]
    percent_nz = [0.1,0.5,1,2,5,10,20,30]
    print("Generating datasets and saving them..", flush=True)
    count = 0
    for csize in tqdm(cell_sizes):
        for gsize in tqdm(gene_sizes,leave=False):
            for percent in tqdm(percent_nz, leave=False):
                if count < 50:
                    output = {}
                    percent /= 100
                
                    matrix = rand(gsize, csize, density=percent, dtype=np.uint8, format="csr", random_state=42)
                    sample_IDs = [''.join(random.choices(string.ascii_uppercase + string.digits,k = 10)) for _ in range(csize)]
                    gene_IDs = [''.join(random.choices(string.ascii_uppercase + string.digits,k = 20)) for _ in range(gsize)]
                    assay = {
                        "matrix":matrix.todense().tolist(),
                        "sampleIds":sample_IDs,
                        "geneIds":gene_IDs
                        }
            
                    name = "_".join((str(gsize), str(csize), str(percent*100)))
                    output[name] = compress_assay(assay)
                    with open(os.path.join("./datasets",name),"wt") as f:
                        json.dump(output, f)
                    del matrix, sample_IDs, gene_IDs, assay, output
                    gc.collect()
                count += 1

                #name = "_".join((str(gsize), str(csize), str(percent*100))) + ".npz"
                #sparse.save_npz(os.path.join('./datasets',name), matrix)
                #del (matrix)
                #gc.collect()

if __name__ == "__main__":
    main()
© 2022 GitHub, Inc.
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About
