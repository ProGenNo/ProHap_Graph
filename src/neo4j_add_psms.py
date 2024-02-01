import argparse
import pandas as pd
import requests
from neo4j import GraphDatabase
from neo4j_commands import Neo4jCommands

parser = argparse.ArgumentParser(
    description='Adds peptides and spectra from provided input files to the graph.')

parser.add_argument("-psm", dest="psm_file", required=True,
                    help="PSM CSV file")

parser.add_argument("-qval_col", dest="qval_column", required=True,
                    help="q-value column name in the PSM table")

parser.add_argument("-qval_thr", dest="qval_threshold", required=False, type=int,
                    help="maximum acceptable q-value, default: 0.01", default=0.01)

parser.add_argument("-mf", dest="metadata_file", required=True,
                    help="sample metadata CSV file")
                
parser.add_argument("-tr_id", dest="transcript_ids", required=True,
                    help="csv file mapping protein IDs to transcript IDs")

parser.add_argument("-hap_csv", dest="haplo_tsv", required=True,
                    help="haplotypes TSV file")

parser.add_argument("-frag_col", dest="frag_column", required=True,
                    help="fragmentation technique column in metadata")      

parser.add_argument("-prot_col", dest="proteases_column", required=True,
                    help="proteases column in metadata")                    

parser.add_argument("-instr_col", dest="instrument_column", required=True,
                    help="spectrometer name column in metadata")
                    
parser.add_argument("-ID_col", dest="ID_column", required=True,
                    help="rawfile ID column in metadata")                  

parser.add_argument("-tissue_col", dest="tissue_column", required=True,
                    help="tissue name column in metadata")            

parser.add_argument("-age_col", dest="age_column", required=False,
                    help="donor age column in metadata")    

parser.add_argument("-sex_col", dest="sex_column", required=False,
                    help="donor sex column in metadata")           

parser.add_argument("-sample_id", dest="sample_id_column", required=True,
                    help="sample or experiment ID column in metadata")      

parser.add_argument('-pride_acc', dest="pride_accession", required=False,
                    help="accession number of the Pride dataset")

parser.add_argument("-uri", dest="neo4j_uri", required=True,
                    help="the URI address of the neo4j instance")

parser.add_argument("-usr", dest="neo4j_usr", required=True,
                    help="the username for the neo4j instance")

parser.add_argument("-pwd", dest="neo4j_pwd", required=True,
                    help="the password for the neo4j instance")

args = parser.parse_args()

print ("Reading", args.psm_file)
psm_df = pd.read_csv(args.psm_file, header=0)
psm_df = psm_df[psm_df[args.qval_column] <= args.qval_threshold]

print ("Reading", args.metadata_file)
meta_df = pd.read_csv(args.metadata_file, header=0)
meta_df.set_index(args.ID_column, inplace=True)

print ("Reading", args.transcript_ids)
tr_id_df = pd.read_csv(args.transcript_ids, header=0)
tr_id_df['TranscriptID'] = tr_id_df['TranscriptID'].apply(lambda x: x.split('.',1)[0])
tr_id_df.set_index('ProteinID', inplace=True)

print ("Reading", args.haplo_tsv)
haplo_df = pd.read_csv(args.haplo_tsv, header=0, sep='\t')
haplo_df.set_index('HaplotypeID', inplace=True)

print ('Connecting to DB')
DB_connection = GraphDatabase.driver(uri=args.neo4j_uri, auth=(args.neo4j_usr, args.neo4j_pwd))
session = DB_connection.session()

peptide_match_commands = {}
spectrum_match_commands = {}
sample_match_commands = {}

psms_log_file = open('tmp/psms_added.log', 'w')
psm_processed_count = 0
total_psm_count = len(psm_df)


print ('Adding', total_psm_count, 'PSMs')

for index,psm_row in psm_df.iterrows():
    rawfile_ID = psm_row['rawfile_ID']
    peptide_seq = psm_row['Sequence']
    peptide_hash = 'pep_' + hex(hash(peptide_seq))[2:]

    query_str = ""

    if (psm_row['psm_type1'] == 'decoy') or (psm_row['psm_type1'] == 'contaminant'):
        psm_processed_count += 1
        continue

    # create or match peptide
    if peptide_hash not in peptide_match_commands:
        query_str = "CREATE " + Neo4jCommands.create_peptide_command(peptide_hash, peptide_seq, psm_row['psm_type1'], psm_row['psm_type2'])
        peptide_match_commands[peptide_hash] = 'MATCH (' + peptide_hash + ':Peptide {id:\'' + peptide_hash + '\'})'

        matching_proteins = psm_row['matching_proteins'].split(';')
        matching_RFs = psm_row['reading_frames'].split(';')
        protein_positions = psm_row['positions_in_proteins'].split(';')

        canonical_transcripts = [ tr_id_df.loc[prot_name]['TranscriptID'] for prot_name in matching_proteins if prot_name.startswith('ENSP') ]

        for i,prot_name in enumerate(matching_proteins):
            
            if prot_name.startswith('haplo'):
                # Check if there is a canonical match to a protein product of this transcript. If so - don't create an edge to the haplotypes for efficiency
                trID = haplo_df.loc[prot_name]['TranscriptID']
                if (trID in canonical_transcripts):
                    continue

                protID = prot_name + '_rf' + matching_RFs[i]
            else:
                protID = prot_name.replace('.', '_')

            query_str = 'MATCH (' + protID + ':Proteoform {id: \'' + protID + '\'}) ' + query_str

            query_str += ', (' + peptide_hash + ')-[:MAPS_TO {position: ' + protein_positions[i] + '}]->(' + protID + ')'
    else:
        query_str = peptide_match_commands[peptide_hash] + ' '

    # create or match spectrum
    spectrum_hash = 'spec_' + hex(hash(rawfile_ID + psm_row['SpectrumTitle']))[2:]

    if spectrum_hash not in spectrum_match_commands:
        # get metadata
        metadata_row = meta_df.loc[rawfile_ID]

        frag_tech = metadata_row[args.frag_column]
        proteases = metadata_row[args.proteases_column]
        instrument = metadata_row[args.instrument_column]
        indiv_age = metadata_row[args.age_column]
        if (args.sex_column):
            indiv_sex = metadata_row[args.sex_column]
        else:
            indiv_sex = '-'
        tissue_name = metadata_row[args.tissue_column]
        sample_ID = metadata_row[args.sample_id_column]

        if 'CREATE' in query_str:
            query_str += ', '
        else: 
            query_str += 'CREATE '

        query_str += Neo4jCommands.create_spectrum_command(spectrum_hash, psm_row['SpectrumTitle'], frag_tech, proteases, instrument, rawfile_ID, psm_row['measured_mz'], -1, psm_row['measured_rt'])
        spectrum_match_commands[spectrum_hash] = 'MATCH (' + spectrum_hash + ':Spectrum {id:\'' + spectrum_hash + '\'})'

        # create or match sample
        if sample_ID not in sample_match_commands:
            query_str += ', ' + Neo4jCommands.create_sample_command(sample_ID, args.pride_accession, tissue_name, indiv_age, indiv_sex)
            sample_match_commands[sample_ID] = 'MATCH (' + sample_ID + ':Sample {id:\'' + sample_ID + '\'})'
        else:
            query_str = sample_match_commands[sample_ID] + ' ' + query_str

        query_str += ', (' + spectrum_hash + ')-[:MEASURED_FROM]->(' + sample_ID + ')'
    else:
        continue
        #query_str = spectrum_match_commands[spectrum_hash] + ' ' + query_str
        #print ('[Warning]: Spectrum', spectrum_hash, 'matched to multiple peptides!')

    query_str += ', ' + Neo4jCommands.connect_peptide_spectrum_command(psm_row, peptide_hash, spectrum_hash)
    session.run(query_str)
    psms_log_file.write(psm_row['PSMId'] + '\n')
    psm_processed_count += 1
    print('Processed:', psm_processed_count, '/', total_psm_count, end='\r')

psms_log_file.close()

# Create index on peptide sequence for faster search and close session
#session.run('CREATE TEXT INDEX pep_seq_index FOR (n:Peptide) ON (n.sequence)')
session.close()
print ('')
print ('Done')
