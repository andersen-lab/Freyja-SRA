import subprocess
import pandas as pd
import time
new_samples = pd.read_csv('data/new_samples.csv', index_col=0)

# Get the first 5 samples from new_samples, if there are less than 5 samples, get all of them

while new_samples.shape[0] > 0:
    
    subprocess.run(['bash', 'scripts/run_pipeline.sh'], shell=True)
    time.sleep(3600) # sleep for 1 hour
    new_samples = new_samples.iloc[5:]
    new_samples.to_csv('data/new_samples.csv', index=True)
