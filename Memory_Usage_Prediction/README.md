Run the data_download.py to get the loom file.

Then install deepimpute:
$ pip install deepimpute

Then Run the mup.py:
$ python mup.py --input_file_path ./t-cell-activation-human-blood-10XV2.loom --gene_size 3000 \
The script will output sugggested cell size and a regression png file.

