#!/usr/bin/env python3
import argparse
import subprocess
import pandas as pd
import os

parser = argparse.ArgumentParser(description='Aggregate covariants outputs')
parser.add_argument('--baseDir', type=str, help='Base directory')

args = parser.parse_args()
agg_df = pd.DataFrame(columns=['Covariants','Sample', 'Count', 'Max_count', 'Freq', 'Coverage_start', 'Coverage_end'])

for covar_path in os.listdir('${baseDir}/outputs/covariants'):
    df = pd.read_csv('${baseDir}/outputs/covariants/' + covar_path,sep='\t')
    sample_name = covar_path.split('/')[-1].split('.')[0]
    df['Sample'] = sample_name
    agg_df = pd.concat((agg_df,df), axis=0)

agg_df = agg_df.set_index(['Covariants','Sample'])
agg_df.to_csv('aggregate_covariants.tsv',sep='\t')
subprocess.run(["gzip", "aggregate_covariants.tsv"])