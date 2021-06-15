import granatum_sdk
import requests
import os
from tqdm import tqdm
import shutil

def download_file(url, output_path):
    url = url.replace('/fetch', '')  # Work around https://github.com/DataBiosphere/azul/issues/2908
    
    response = requests.get(url, stream=True)
    response.raise_for_status()
    
    total = int(response.headers.get('content-length', 0))
    print(f'Downloading to: {output_path}', flush=True)
    
    with open(output_path, 'wb') as f:
        with tqdm(total=total, unit='B', unit_scale=True, unit_divisor=1024) as bar:
            for chunk in response.iter_content(chunk_size=1024):
                size = f.write(chunk)
                bar.update(size)


def iterate_matrices_tree(tree, keys=()):
    if isinstance(tree, dict):
        for k, v in tree.items():
            yield from iterate_matrices_tree(v, keys=(*keys, k))
    elif isinstance(tree, list):
        for file in tree:
            yield keys, file
    else:
        assert False


def download_data(ProjectID):
    print(os.getcwd())
    # destroy if exits, and then create ./tmp_datasets directory
    dirpath = './tmp_datasets'
    if os.path.exists(dirpath) and os.path.isdir(dirpath):
        shutil.rmtree(dirpath)
    os.mkdir(dirpath)

    project_uuid = ProjectID
    catalog = 'dcp5'
    endpoint_url = f'https://service.azul.data.humancellatlas.org/index/projects/{project_uuid}'
    save_location = dirpath
    response = requests.get(endpoint_url, params={'catalog': catalog})
    response.raise_for_status()

    response_json = response.json()
    project = response_json['projects'][0]

    file_urls = set()
    for key in ('matrices', 'contributorMatrices'):
        tree = project[key]
        for path, file_info in iterate_matrices_tree(tree):
            url = file_info['url']
            if url not in file_urls:
                dest_path = os.path.join(save_location, file_info['name'])
                download_file(url, dest_path)
                file_urls.add(url)



def main():
    gn = granatum_sdk.Granatum()

    #assay = gn.get_import('assay')
    #assay_df = gn.pandas_from_assay(assay)

    #checkbox_value = gn.get_arg('someCheckbox')
    #number_value = gn.get_arg('someNumber')
    #seed_value = gn.get_arg('someSeed')

    # markdown_str = f"""\
    # * checkbox_value = **{checkbox_value}**
    # * number_value = **{number_value}**
    # * seed_value = **{seed_value}**
    # * Shape of the assay = **{assay_df.shape}**"""

    #gn.add_result(markdown_str, data_type='markdown')

    #gn.commit()
    ProjectID = gn.get_arg('PID', '4a95101c-9ffc-4f30-a809-f04518a23803')
    download_data(ProjectID)
    print(os.listdir('./tmp_datasets'))






if __name__ == "__main__":
    main()
