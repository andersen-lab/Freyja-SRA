from elasticsearch import Elasticsearch
import eland as ed
import pandas as pd

USER = 'elastic'
PASS = 'g*tfsy4uJIXUtyBrFhHs'
es = Elasticsearch(['https://localhost:9200'], basic_auth=(USER, PASS), verify_certs=False)

df = pd.read_csv('output/demix/aggregate.tsv', sep='\t')
df.rename(columns={'Unnamed: 0': 'accession'}, inplace=True)
df['accession'] = df['accession'].str.split('.').str[0]

df.index = df['accession']
df.drop(columns=['accession'], inplace=True)

print(df.head())

INDEX = 'demix'

df = ed.pandas_to_eland(
    pd_df=df,
    es_client=es,
    # Where the data will live in Elasticsearch
    es_dest_index="demix",
    es_if_exists="append",
    # Wait for data to be indexed before returning
    es_refresh=True,
)


