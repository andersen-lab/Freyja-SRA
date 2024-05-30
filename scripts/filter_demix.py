import os

rerun = os.listdir('demix_updated')
existing = os.listdir('outputs/demix')

# If a file in rerun is in existing, remove it from existing

for file in rerun:
    if file in existing:
        os.remove(f'outputs/demix/{file}')

# Copy all files in rerun to outputs/demix
for file in rerun:
    os.system(f'cp demix_updated/{file} outputs/demix/{file}')
