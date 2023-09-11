import eland as ed
from elasticsearch import Elasticsearch

USER = 'elastic'
PASS = 'changeme'
CERTIFICATE = '/etc/ssl/certs/ca-certificates.crt'

es = Elasticsearch(['https://localhost:9200'], basic_auth=(USER, PASS), verify_certs=False, ca_certs=CERTIFICATE)
df = ed.DataFrame(es_client=es, es_index_pattern="demix")
df = df.to_pandas()


print(df.loc['SRR25726357'])
