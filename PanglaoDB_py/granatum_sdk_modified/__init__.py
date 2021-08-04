import io
import json
import pickle
import random
import string
import sys
from base64 import b64encode
from glob import glob
from os import environ, path
import zipfile
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import scanpy as sc

from .gene_id_helpers import *

from scipy.sparse import coo_matrix


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def random_string(n=5):
    return "".join([random.choice(list(string.ascii_uppercase + string.digits)) for _ in range(n)])


class Granatum:

    def __init__(self, swd="/data"):
        self.swd = environ.get("GRANATUM_SWD", swd)
        self.uploaded_files_dir = path.join(self.swd, "uploaded_files")
        self.exports_dir = path.join(self.swd, "exports")
        self.exports_anno_file = path.join(self.swd, "exports_anno.json")
        self.imports_dir = path.join(self.swd, "imports")
        self.args_file = path.join(self.swd, "args.json")
        self.results_file = path.join(self.swd, "results.json")
        self.dynamic_exports = []
        self.results = []

        self.debug_dir = path.join(self.swd, "debug")

        with open(path.join(self.swd, "args.json")) as f:
            try:
                self.args = json.load(f)
            except FileNotFoundError:
                self.args = {}

    # -- uploaded_files -------------------------------------------------

    def get_uploaded_file_path(self, inject_into):
        try:
            return glob(path.join(self.uploaded_files_dir, inject_into, "*"))[0]
        except IndexError:
            return None

    # -- imports -------------------------------------------------

    def get_import(self, inject_into):
        import_file = path.join(self.imports_dir, inject_into)
        with open(import_file, "r") as f:
            return json.load(f)

    ''' return a list of imported files '''
    def get_import_from_zip(self, inject_into):
        import_file = path.join(self.imports_dir, inject_into)
        new_file_path = path.join((self.imports_dir), inject_into.replace("zip", "")

        with zipfile.ZipFile(import_file,"r") as zip_ref:
            zip_ref.extractall(new_file_path)

        ret = []
        for file in os.listdir(new_file_path):
            ret.append(self.get_import(file))
        return ret
    # -- args -------------------------------------------------

    def get_arg(self, inject_into, default=None):
        return self.args.get(inject_into, default)

    # -- exports  -------------------------------------------------

    def export_statically(self, data, extract_from, raw=False):
        if raw:
            with open(path.join(self.exports_dir, extract_from), "w") as f:
                f.write(data)
        else:
            with open(path.join(self.exports_dir, extract_from), "w") as f:
                json.dump(data, f)

    def export(self, data, extract_from, kind=None, meta=None, dynamic=True, raw=False):
        if dynamic:
            self.dynamic_exports.append({"extractFrom": extract_from, "kind": kind, "meta": meta})

        self.export_statically(data, extract_from, raw)

    def export_current_figure(self, extract_from, format_='pdf', zoom=2, width=750, height=650, dpi=100, **kwargs):
        fig = plt.gcf()
        fig.set_figheight(height / dpi)
        fig.set_figwidth(width / dpi)
        bs = io.BytesIO()
        fig.savefig(bs, dpi=zoom * dpi, format=format_)

        self.export(bs, extract_from, **kwargs, raw=True)

    # -- results  -------------------------------------------------

    def add_current_figure_to_results(
        self, description=None, zoom=2, width=750, height=650, dpi=100, savefig_kwargs={}
    ):
        save_filepath = path.join("/tmp", random_string() + ".png")

        fig = plt.gcf()
        fig.set_figheight(height / dpi)
        fig.set_figwidth(width / dpi)
        fig.savefig(save_filepath, dpi=zoom * dpi, **savefig_kwargs)

        with open(save_filepath, "rb") as f:
            image_b64 = b64encode(f.read()).decode("utf-8")

        self.results.append(
            {
                "type": "png",
                "width": width,
                "height": height,
                "description": description,
                "data": image_b64,
            }
        )

    def add_pandas_df(self, df, description=None):
        self.results.append(
            {
                "type": "table",
                "data": {
                    "caption": description,
                    "columns": json.loads(json.dumps(df.columns.values.tolist())),
                    "data": json.loads(df.to_json(orient='records')),
                }
            }
        )

    def add_result(self, data, data_type="json", **kwargs):
        self.results.append({"type": data_type, "data": data, **kwargs})

    def add_markdown(self, md, **kwargs):
        self.results.append({"type": "markdown", "data": md, **kwargs})

    # -- commit  -------------------------------------------------

    def commit(self):
        with open(self.exports_anno_file, "w") as f:
            json.dump(self.dynamic_exports, f)

        with open(self.results_file, "w") as f:
            json.dump(self.results, f)

    # -- helpers -------------------------------------------------

    def ann_data_from_assay(self, assay):

        sparse_matrix = coo_matrix(assay.get("matrix")).tocsc()
        adata = sc.AnnData(sparse_matrix.transpose())

        adata.var_names = assay.get("geneIds")
        adata.obs_names = assay.get("sampleIds")
        return adata

    def assay_from_ann_data(self, adata):
        return {
            "matrix": adata.X.T.toarray().tolist(),
            "sampleIds": adata.obs_names.tolist(),
            "geneIds": adata.var_names.tolist(),
        }

    def pandas_from_assay(self, assay, samples_as_rows=False):

        sparse_matrix = coo_matrix(assay.get("matrix"))
        df = pd.DataFrame.sparse.from_spmatrix(sparse_matrix, index=assay.get("geneIds"), columns=assay.get("sampleIds"))
        
        if samples_as_rows:
            df = df.transpose()
        return df

    def assay_from_pandas(self, df, samples_as_rows=False):
        if samples_as_rows:
            df = df.transpose()
        return {
            "matrix": df.where((pd.notnull(df)), 0).values.tolist(),
            "geneIds": df.index.values.tolist(),
            "sampleIds": df.columns.values.tolist(),
        }

    # -- error  -------------------------------------------------

    def error(self, message):
        raise Exception(message)

    def _pickle(self, data, filename="debug.pickle"):
        with open(path.join(self.debug_dir, filename), "wb") as f:
            pickle.dump(data, f)
