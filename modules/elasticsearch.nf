/*
 * Check which samples are present in the database (if any)
 */
process CHECK_SAMPLES_IN_ES {
    container = 'dylanpilz/elastic-eland-amd64:latest'

    input:
    val host
    val user
    val password

    output:
    path 'acc_in_db.csv'

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import eland as ed
    from elasticsearch import Elasticsearch

    es = Elasticsearch(['${host}'], basic_auth=('${user}', '${password}'), verify_certs=False)
    try:
        df = ed.DataFrame(es, es_index_pattern='demix').to_pandas()
    except:
        df = pd.DataFrame()
    pd.Series(df.index).to_csv('acc_in_db.csv', index=False, header=False)
    """
}

process PUSH_TO_ES {
    container = 'dylanpilz/elastic-eland-amd64:latest'

    input:
    val host
    val user
    val password
    path variants
    path demix
    path covariants

    script:
    """
    #!/usr/bin/env python3

    import pandas as pd
    import eland as ed
    from elasticsearch import Elasticsearch

    es = Elasticsearch(['${host}'], basic_auth=('${user}', '${password}'), verify_certs=False)

    for df, index in zip([pd.read_csv('${variants}', sep='\t'), pd.read_csv('${demix}', sep='\t'), pd.read_csv('${covariants}', sep='\t')], ['variants', 'demix', 'covariants']):
        if index == 'demix':
            df.rename(columns={'Unnamed: 0': 'accession'}, inplace=True)
            df['accession'] = df['accession'].str.split('.').str[0]
            df.index = df['accession']
            df.drop(columns=['accession'], inplace=True)

        print(df.head())

        df = ed.pandas_to_eland(
            pd_df=df,
            es_client=es,
            # Where the data will live in Elasticsearch
            es_dest_index=index,
            es_if_exists="append",
            # Wait for data to be indexed before returning
            es_refresh=True,
        )
    """
}