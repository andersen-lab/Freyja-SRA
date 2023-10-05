import subprocess
import pandas as pd

new_samples = pd.read_csv('data/new_samples.csv', index_col=0)

cmd = [
    'nextflow run main.nf',
    '--input data/new_samples.csv',
    '--num_samples 5',
    '-profile docker',
    '-entry fetch_sra'
]

sample_to_run = len(new_samples)
samples_finished = 0
while samples_finished < sample_to_run:
    subprocess.run(cmd, shell=True)
    subprocess.run('rm -rf work', shell=True)