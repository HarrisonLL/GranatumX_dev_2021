import re
import pandas as pd
import pkg_resources
from collections import Counter


def get_all_genes(species, column='Symbol'):
    biomart_tb = pd.read_csv(pkg_resources.resource_stream(__name__, f"data/biomart_{species}.csv"))
    return biomart_tb[column].dropna().unique().tolist()


def guess_gene_id_type(ids_):
    if type(ids_) is list:
        return Counter(map(guess_gene_id_type, ids_)).most_common(1)[0][0]
    elif type(ids_) is int:
        return 'entrez'
    elif type(ids_) is str:
        if re.match(r'^[0-9]+$', ids_):
            return 'entrez'
        elif re.match(r'^ENS[A-Z]+[0-9]{11}', ids_):
            return 'ensembl'
        else:
            return 'symbol'
    else:
        return 'unknown'


biomart_col_dict = {
    'ensembl': 'Ensembl Gene ID',
    'ensembl_trans': 'Ensembl Transcript ID',
    'entrez': 'NCBI gene ID',
    'symbol': 'Symbol',
}


def convert_gene_ids(ids_, from_type, to_type, species, return_new_meta=False):
    from_col = biomart_col_dict[from_type]
    to_col = biomart_col_dict[to_type]

    biomart_tb = pd.read_csv(pkg_resources.resource_stream(__name__, f"data/biomart_{species}.csv"))

    if from_col != to_col:
        lookup_tb = biomart_tb.groupby(from_col).first()
        info_tb = lookup_tb.loc[ids_, :]

        new_ids_ = info_tb[to_col].values.tolist()
        new_meta = info_tb.reset_index().set_index(to_col)
    else:
        lookup_tb = biomart_tb.groupby(from_col).first()
        new_meta = lookup_tb.loc[ids_, :]

        new_ids_ = ids_

    if return_new_meta:
        return (new_ids_, new_meta)
    else:
        return new_ids_
