import requests
import os
from tqdm import tqdm

# 4a95101c-9ffc-4f30-a809-f04518a23803
# 3089d311-f9ed-44dd-bb10-397059bad4dc
project_uuid = '3089d311-f9ed-44dd-bb10-397059bad4dc'
catalog = 'dcp6'
endpoint_url = f'https://service.azul.data.humancellatlas.org/index/projects/{project_uuid}'
save_location = './'

def download_file(url, output_path):
    url = url.replace('/fetch', '')  # Work around https://github.com/DataBiosphere/azul/issues/2908
    
    response = requests.get(url, params={'catalog': catalog})
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


response = requests.get(endpoint_url)
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

            # only download DCP-processed files with selected organ type and species
            name = file_info['name']
            if  (name.endswith('.loom')) and 'blood' in name.lower() and 'human' in name.lower():
                download_file(url, dest_path)

            file_urls.add(url)
print('Downloads Complete.')
