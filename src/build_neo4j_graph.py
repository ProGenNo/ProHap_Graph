import argparse
import pandas as pd
from neo4j import GraphDatabase
from neo4j_commands import Neo4jCommands
from common import read_fasta
import gffutils

parser = argparse.ArgumentParser(
    description='Builds a neo4j graph from provided input files.')

parser.add_argument("-hap_fa", dest="haplo_fasta", required=True,
                    help="haplotypes fasta file")

parser.add_argument("-ref_fa", dest="reference_fasta", required=True,
                    help="reference fasta file")

parser.add_argument("-cdna_fa", dest="cdna_fasta", required=True,
                    help="cDNA fasta file (Ensembl)")

parser.add_argument("-hap_db", dest="haplo_tsv", required=True,
                    help="haplotypes TSV file")

parser.add_argument("-annot_db", dest="annotation_db", required=True,
                    help="annotations database file")

parser.add_argument("-tr_id", dest="transcript_ids", required=True,
                    help="csv file mapping protein IDs to transcript IDs")

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")

parser.add_argument("-uri", dest="neo4j_uri", required=True,
                    help="the URI address of the neo4j instance")

parser.add_argument("-usr", dest="neo4j_usr", required=True,
                    help="the username for the neo4j instance")

parser.add_argument("-pwd", dest="neo4j_pwd", required=True,
                    help="the password for the neo4j instance")

args = parser.parse_args()

print ('Reading', args.annotation_db)
# Load the annotations database
annotations_db = gffutils.FeatureDB(args.annotation_db)

print ("Reading", args.transcript_ids)
tr_id_df = pd.read_csv(args.transcript_ids, header=0)
all_transcript_ids = tr_id_df['TranscriptID'].apply(lambda x: x.split('.',1)[0]).drop_duplicates().tolist()
tr_id_df.set_index('ProteinID', inplace=True)

cDNA_sequences = read_fasta(args.cdna_fasta)

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)
all_gene_ids = gene_id_df['GeneID'].drop_duplicates().tolist()
gene_id_df.set_index('TranscriptID', inplace=True)

unique_haplotype_count = 0

print ('Connecting to DB')
DB_connection = GraphDatabase.driver(uri=args.neo4j_uri, auth=(args.neo4j_usr, args.neo4j_pwd))
session = DB_connection.session()

print ('Adding genes and transcripts')

exon_match_commands = {}
gene_match_commands = {}
transcript_match_commands = {}

session.run('CREATE TEXT INDEX gID_index for (n:Gene) on (n.id)')
session.run('CREATE TEXT INDEX trID_index for (n:Transcript) on (n.id)')
session.run('CREATE TEXT INDEX exID_index for (n:Exon) on (n.id)')

# create all the genes, exons and transcripts, edges between them
for i,geneID in enumerate(all_gene_ids[:]):
    gene_feature = annotations_db[geneID]

    if (geneID not in gene_match_commands):
        query_str = 'CREATE ' 
        query_str += Neo4jCommands.create_gene_command(gene_feature)
        gene_match_commands[geneID] = "MATCH (" + geneID + ':Gene {id: \'' + geneID + '\'})'

        matched_transcripts = []
        matched_exons = []

        for transcript_feature in annotations_db.children(gene_feature, featuretype='transcript'):
            transcriptID = transcript_feature.id

            # transcript not used to create haplotype database -> ignore
            if (transcriptID not in all_transcript_ids):
                continue

            # transcript already created -> don't create again
            if (transcriptID not in transcript_match_commands):
                transcript_cDNA = cDNA_sequences[transcriptID]['sequence']

                query_str += ', ' +  Neo4jCommands.create_transcript_command(transcript_feature, annotations_db, transcript_cDNA)
                query_str += ', (' + transcriptID + ')-[:TRANSCRIPT_OF]->(' + geneID + ')'
                transcript_match_commands[transcriptID] = "MATCH (" + transcriptID + ':Transcript {id: \'' + transcriptID + '\'})'
                matched_transcripts.append(transcriptID)                                # make sure not to match this node twice

            # transcript already created or matched in this query
            elif transcriptID not in matched_transcripts:
                query_str = transcript_match_commands[transcriptID] + ' ' + query_str
                matched_transcripts.append(transcriptID)                                # make sure not to match this node twice

            for exon_feature in annotations_db.children(transcript_feature, featuretype='exon'):
                exonID = exon_feature.attributes['exon_id'][0]   
                exon_number = exon_feature.attributes['exon_number'][0]             

                # exon already created -> don't create again, connect to the gene only when it's created
                if (exonID not in exon_match_commands):
                    query_str += ', ' + Neo4jCommands.create_exon_command(exon_feature)
                    exon_match_commands[exonID] = "MATCH (" + exonID + ':Exon {id: \'' + exonID + '\'})'
                    matched_exons.append(exonID)        # make sure not to match this node twice

                # exon already created or matched in this query
                elif exonID not in matched_exons:
                    query_str = exon_match_commands[exonID] + ' ' + query_str
                    matched_exons.append(exonID)        # make sure not to match this node twice

                # connect the exon to this transcript
                query_str += ', (' + transcriptID + ')' + '-[:INCLUDES_EXON {exon_number: ' + exon_number + '}]->(' + exonID + ')'

        #session.run(query_str)
        print('Processed:', i, '/', len(all_gene_ids), end='\r')

# add the haplotypes and protein sequences
haplotype_match_commands = {}
haplotype_transcript_edges = {}
proteoform_match_commands = {}
variant_match_commands = {}
skipped_haplotypes = 0

print ()
print ("Reading", args.haplo_tsv)
haplo_df = pd.read_csv(args.haplo_tsv, header=0, sep='\t')
haplo_df.set_index('HaplotypeID', inplace=True)

print ('Adding proteins and haplotypes')
session.run("CREATE TEXT INDEX varID_index for (n:Variant) on (n.id)")
session.run("CREATE TEXT INDEX haploID_index for (n:Haplotype) on (n.id)")
session.run('CREATE TEXT INDEX protID_index for (n:Proteoform) on (n.id)')

haplo_proteins = read_fasta(args.haplo_fasta)
for prot_idx,fasta_elem in enumerate(haplo_proteins.values()):
    matching_haplotypes = [ elem[0] for elem in fasta_elem['matching_proteins'] ]

    for hap_idx,haploID in enumerate(matching_haplotypes):
        #try:
        haplotype_row = haplo_df.loc[haploID]
        #except:
        #    continue

        # skip haplotypes that don't have the corresponding transcript
        if (haplotype_row['TranscriptID'] not in transcript_match_commands):
            skipped_haplotypes += 1
            continue

        geneID = gene_id_df.loc[haplotype_row['TranscriptID']]['GeneID']    # Gene ID
        chromosome = haplotype_row['chromosome']                            # chromosome
        DNA_changes = haplotype_row['DNA_changes']                          # String: DNA changes
        reading_frame = int(fasta_elem['reading_frames'][hap_idx][0])             # Reading frame
        protein_changes = haplotype_row['all_protein_changes']              # String: all protein changes
        foo = haplotype_row['frequency']                                    # FoO: frequency of occurrence

        # In case of 3-frame translation pick the protein changes with the correct RF
        if (haplotype_row['reading_frame'] == -1):
            protein_changes = ';'.join([ ch.split('|')[reading_frame] for ch in protein_changes.split(';') ])

        # init the query string - match the corresponding gene and transcript
        query_str = gene_match_commands[geneID] + ' ' + transcript_match_commands[haplotype_row['TranscriptID']] + ' CREATE '
        matched_variants = []   # IDs of variants that have already been matched in this query

        # create the haplotype node
        haplo_hash = 'hap_' + hex(hash(geneID + DNA_changes))[2:]
        is_new_haplotype = haplo_hash not in haplotype_match_commands
        if is_new_haplotype:
            unique_haplotype_count += 1
            query_str += Neo4jCommands.create_haplotype_command(unique_haplotype_count, haplo_hash, geneID)
            haplotype_match_commands[haplo_hash] = Neo4jCommands.match_haplotype_command(unique_haplotype_count, haplo_hash, geneID)
            haplotype_transcript_edges[haplo_hash] = [haplotype_row['TranscriptID']]

            # connect haplotype node to transcript
            query_str += ', (' + haplo_hash + ')-[:HAPLO_FORM_OF {frequency: ' + str(foo) + '}]->(' + haplotype_row['TranscriptID'] + ')'
        else:
            query_str = haplotype_match_commands[haplo_hash] + ' ' + query_str

            # this haplotype might need to be connected to another transcript
            if haplotype_row['TranscriptID'] not in haplotype_transcript_edges[haplo_hash]:
                if not query_str.endswith('CREATE '):
                    query_str += ', '
                query_str += '(' + haplo_hash + ')-[:HAPLO_FORM_OF {frequency: ' + str(foo) + '}]->(' + haplotype_row['TranscriptID'] + ')'
                haplotype_transcript_edges[haplo_hash].append(haplotype_row['TranscriptID'])

        # create the proteoform node
        proteoform_id = haploID + '_rf' + str(reading_frame)

        if proteoform_id not in proteoform_match_commands:
            if not query_str.endswith('CREATE '):
                query_str += ', '
            query_str += Neo4jCommands.create_proteoform_haplo_command(proteoform_id, fasta_elem['sequence'], reading_frame, haplotype_row, protein_changes)
            proteoform_match_commands[proteoform_id] = 'MATCH (' + proteoform_id + ':Proteoform {id: \'' + proteoform_id + '\'})'
        
            query_str += ', (' + proteoform_id + ')-[:ENCODED_BY_HAPLOTYPE]->(' + haplo_hash + ')'
            query_str += ', (' + proteoform_id + ')-[:ENCODED_BY_TRANSCRIPT]->(' + haplotype_row['TranscriptID'] + ')'

        # create the variant nodes
        for i,change in enumerate(DNA_changes.split(';')):
            change_id = str(chromosome) + ':' + change
            change_id_safe = ('chr' + str(chromosome) + ':' + change).replace(':', '_').replace('>', 'x')   # for compatibility with the query language

            if change_id_safe not in variant_match_commands:
                query_str += ', ' + Neo4jCommands.create_variant_command(change_id, change_id_safe)
                query_str += ', (' + change_id_safe + ')-[:VARIANT_MAPS_TO]->(' + geneID + ')'
                variant_match_commands[change_id_safe] = 'MATCH (' + change_id_safe + ':Variant {id: \'' + change_id + '\'})'
                matched_variants.append(change_id_safe)
            
            elif change_id_safe not in matched_variants:
                query_str = variant_match_commands[change_id_safe] + ' ' + query_str

            # connect variants to haplotype only if they haven't been connnected already
            if is_new_haplotype:
                query_str += ', (' + haplo_hash + ')-[:INCLUDES_ALT_ALLELE {var_order: ' + str(i) + '}]->(' + change_id_safe + ')'

        session.run(query_str)
    print('Processed:', prot_idx, '/', len(haplo_proteins), end='\r')

print()
print('Adding canonical proteins')

ref_proteins = read_fasta(args.reference_fasta)
for i,fasta_elem in enumerate(ref_proteins.values()):
    protID = fasta_elem['accession'].replace('.', '_')
    transcriptID = fasta_elem['description'].split('transcript:',1)[1].split('.',1)[0]

    if (protID not in proteoform_match_commands) and (transcriptID in transcript_match_commands):
        proteoform_match_commands[protID] = 'MATCH (' + protID + ':Proteoform {id: \'' + protID + '\'})'

        query_str = transcript_match_commands[transcriptID] + ' CREATE '

        query_str += Neo4jCommands.create_proteoform_reference_command(protID, fasta_elem)
        query_str += ', (' + protID + ')-[:ENCODED_BY_TRANSCRIPT]->(' + transcriptID + ')'

        session.run(query_str)

    print('Processed:', i, '/', len(ref_proteins), end='\r')

print ()
print (skipped_haplotypes, 'haplotypes skipped')

session.close()
print ('Done')
