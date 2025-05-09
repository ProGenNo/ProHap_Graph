import argparse
import pandas as pd
import time
from neo4j import GraphDatabase

parser = argparse.ArgumentParser(
    description='Updates the overview statistics in the graph database.')

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")

parser.add_argument("-uri", dest="neo4j_uri", required=True,
                    help="the URI address of the neo4j instance")

parser.add_argument("-usr", dest="neo4j_usr", required=True,
                    help="the username for the neo4j instance")

parser.add_argument("-pwd", dest="neo4j_pwd", required=True,
                    help="the password for the neo4j instance")

args = parser.parse_args()

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)
all_gene_ids = gene_id_df['GeneID'].drop_duplicates().tolist()
gene_id_df.set_index('TranscriptID', inplace=True)

print ('Connecting to DB')
DB_connection = GraphDatabase.driver(uri=args.neo4j_uri, auth=(args.neo4j_usr, args.neo4j_pwd))
session = DB_connection.session()

print ('Updating overview statistics.')

for i,gID in enumerate(all_gene_ids):
    query_str = 'MATCH (g:Gene {id: \'' + gID + '\'}) '
    query_str += "CALL apoc.path.subgraphAll(g, {relationshipFilter:'<TRANSCRIPT_OF|<EXON_PART_OF|<VARIANT_MAPS_TO|<INCLUDES_ALT_ALLELE|INCLUDES_EXON>|HAPLO_FORM_OF|<ENCODED_BY_TRANSCRIPT|<ENCODED_BY_HAPLOTYPE|<MAPS_TO'}) YIELD nodes RETURN nodes, [ v in nodes | labels(v) ] as node_types;"
    query_response = session.run(query_str).data()
    time.sleep(0.3)

    if (len(query_response) > 0):
        matched_variants = len([ elem for elem in query_response[0]['nodes'] if ('overlapping_peptide' in elem) and (elem['overlapping_peptide']) ])

        # Get all the protein sequences, chop off UTRs
        all_proteoforms = [ elem['sequence'][elem['start_aa']:].split('*')[0] for i,elem in enumerate(query_response[0]['nodes']) if (query_response[0]['node_types'][i][0] == 'Proteoform') ]
        total_proteoforms = len(list(dict.fromkeys(all_proteoforms)))

        total_peptides = len([ elem for elem in query_response[0]['node_types'] if elem[0] == 'Peptide'])
        variant_peptides = len([ elem for elem in query_response[0]['nodes'] if ('pep_class_1' in elem) and (('variant' in elem['pep_class_1']) or ('frameshift' in elem['pep_class_1'])) ])
        

        query_str = 'MATCH (g:Gene {id: \'' + gID + '\'}) '
        query_str += "SET g.matched_vars = " + str(matched_variants) 
        query_str += " SET g.total_proteoforms = " + str(total_proteoforms) 
        query_str += " SET g.total_peptides = " + str(total_peptides) 
        query_str += " SET g.variant_peptides = " + str(variant_peptides) 
        session.run(query_str)
        print('Processed:', i, '/', len(all_gene_ids), end='\r')
        time.sleep(0.3)
    

session.close()
print ('Done')
