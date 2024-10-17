import argparse
import pandas as pd
import time
import re
from neo4j import GraphDatabase
from neo4j_commands import Neo4jCommands

parser = argparse.ArgumentParser(
    description='Adds peptides and spectra from provided input files to the graph.')

parser.add_argument("-psm", dest="psm_file", required=True,
                    help="PSM CSV file")

parser.add_argument("-rawfile_id", dest="rawfile_id", required=False,
                    help="Rawfile name in case results are not concatenated (skip otherwise)", default=None)

parser.add_argument("-qval_col", dest="qval_column", required=False,
                    help="q-value column name in the PSM table", default="q-value")

parser.add_argument("-qval_thr", dest="qval_threshold", required=False, type=float,
                    help="maximum acceptable q-value, default: 0.01", default=0.01)

parser.add_argument("-mf", dest="metadata_file", required=True,
                    help="sample metadata CSV file")

parser.add_argument("-added_peps", dest="added_peps", required=True,
                    help="file containing peptide sequences that are already added")

parser.add_argument("-g_id", dest="gene_ids", required=True,
                    help="csv file mapping transcript IDs to gene IDs")
                
parser.add_argument("-tr_id", dest="transcript_ids", required=True,
                    help="csv file mapping protein IDs to transcript IDs")

parser.add_argument("-hap_tsv", dest="haplo_tsv", required=True,
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

parser.add_argument("-pheno_col", dest="phenotype_column", required=False,
                    help="donor phenotype / disease column in metadata")        

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
psm_df['USI'] = psm_df['USI'].apply(lambda x: '\'' + x + '\'')

print ("Reading", args.metadata_file)
meta_df = pd.read_csv(args.metadata_file, header=0, sep='\t')
meta_df.set_index(args.ID_column, inplace=True)

print ("Reading", args.gene_ids)
gene_id_df = pd.read_csv(args.gene_ids, header=0)
all_transcript_ids = gene_id_df['TranscriptID'].tolist()

print ("Reading", args.transcript_ids)
tr_id_df = pd.read_csv(args.transcript_ids, header=0)
tr_id_df['TranscriptID'] = tr_id_df['TranscriptID'].apply(lambda x: x.split('.',1)[0])
#tr_id_df = tr_id_df[tr_id_df["TranscriptID"].isin(all_transcript_ids)]
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

peps_file = open(args.added_peps, 'r')
added_peptides = [ l for l in peps_file.readlines() if (len(l) > 1)]
peps_file.close()
for pep in added_peptides:
    peptide_match_commands[pep[:-1]] = 'MATCH (pep_' + pep[:-1] + ':Peptide {id:\'pep_' + pep[:-1] + '\'})'

peps_file = open(args.added_peps, 'a')
psm_processed_count = 0
total_psm_count = len(psm_df)

#PSM_features_to_add = ['rt_Abs_error', 'rt_apex_dist', 'spectra_cos_similarity', 'spectra_angular_similarity', 'posterior_error_prob', 'q-value']
PSM_features_to_add = ['posterior_error_probability', 'q-value', 'USI']

print ('Adding', total_psm_count, 'PSMs')

if (len(added_peptides) == 0):
    try:
        session.run('CREATE TEXT INDEX pep_seq_index FOR (n:Peptide) ON (n.sequence)')
    except:
        print('Could not create index on peptide sequence - most likely already exists')

    try:
        session.run('CREATE TEXT INDEX pep_seq_index FOR (n:Peptide) ON (n.id)')
    except:
        print('Could not create index on peptide ID - most likely already exists')

for index,psm_row in psm_df.iterrows():
    rawfile_ID = args.rawfile_id if args.rawfile_id is not None else psm_row['rawfile_ID']
    peptide_seq = psm_row['sequence']

    query_str = ""
    set_query_str = ""  # separate query to set variant properties

    if (psm_row['psm_type1'] == 'decoy') or (psm_row['psm_type1'] == 'contaminant'):
        psm_processed_count += 1
        continue
    
    matched_protein = False

    # create or match peptide
    if peptide_seq not in peptide_match_commands:
        query_str = "CREATE " + Neo4jCommands.create_peptide_command('pep_'+peptide_seq, peptide_seq, psm_row['psm_type1'], psm_row['psm_type2'], psm_row['expected_maximum_frequency'])

        matching_proteins = psm_row['matching_proteins'].split(';')
        matching_RFs = psm_row['reading_frames'].split(';')
        protein_positions = psm_row['positions_in_proteins'].split(';')

        canonical_transcripts = [ tr_id_df.loc[prot_name]['TranscriptID'] for prot_name in matching_proteins if prot_name.startswith('ENSP') ]

        for i,prot_name in enumerate(matching_proteins):
            
            if prot_name.startswith('haplo'):
                # Check if there is a canonical match to a protein product of this transcript. If so - don't create an edge to the haplotypes for efficiency
                trID = haplo_df.loc[prot_name]['TranscriptID']
                if ((trID not in all_transcript_ids) or (trID in canonical_transcripts)):
                    continue

                protID = prot_name + '_rf' + matching_RFs[i]
            else:
                trID = tr_id_df.loc[prot_name]['TranscriptID']
                if (trID not in all_transcript_ids):
                    continue
                protID = prot_name.replace('.', '_')

            query_str = 'MATCH (' + protID + ':Proteoform {id: \'' + protID + '\'}) ' + query_str

            query_str += ', (pep_' + peptide_seq + ')-[:MAPS_TO {position: ' + protein_positions[i] + '}]->(' + protID + ')'
            matched_protein = True

        if (matched_protein):     

            # in case of variant peptides - flag the variant as matched
            if ('variant' in psm_row['psm_type1']):
                for varID in re.split(r"[,;|]", psm_row['covered_alleles_dna']):
                    if ('>' in varID):
                        varID_safe = varID.replace(':', '_').replace('>', 'x')
                        set_query_str = 'MATCH (var_' + varID_safe + ':Variant {id:\'' + varID + '\'}) SET var_' + varID_safe + '.overlapping_peptide = TRUE'

            peptide_match_commands[peptide_seq] = 'MATCH (pep_' + peptide_seq + ':Peptide {id:\'pep_' + peptide_seq + '\'})'       
            peps_file.write(peptide_seq + '\n')

    else:
        query_str = peptide_match_commands[peptide_seq] + ' '
        matched_protein = True

    if not (matched_protein):
        psm_processed_count += 1
        continue

    # create or match spectrum
    spectrum_hash = 'spec_' + hex(hash(rawfile_ID + psm_row['SpectrumTitle']))[2:]

    if spectrum_hash not in spectrum_match_commands:
        # get metadata
        metadata_row = meta_df.loc[rawfile_ID]

        frag_tech = metadata_row[args.frag_column] if args.frag_column != "_" else "-"
        proteases = metadata_row[args.proteases_column] if args.proteases_column != "_" else "-"
        instrument = metadata_row[args.instrument_column] if args.instrument_column != "_" else "-"
        indiv_age = metadata_row[args.age_column] if args.age_column != "_" else "-"
        indiv_sex = metadata_row[args.sex_column] if args.sex_column != "_" else "-"
        tissue_name = metadata_row[args.tissue_column] if args.tissue_column != "_" else "-"
        phenotype = metadata_row[args.phenotype_column] if args.phenotype_column != "_" else "-"
        sample_ID = metadata_row[args.sample_id_column]

        if 'CREATE' in query_str:
            query_str += ', '
        else: 
            query_str += 'CREATE '

        #query_str += Neo4jCommands.create_spectrum_command(spectrum_hash, psm_row['SpectrumTitle'], frag_tech, proteases, instrument, rawfile_ID, psm_row['measured_mz'], -1, psm_row['measured_rt'])
        query_str += Neo4jCommands.create_spectrum_command(spectrum_hash, psm_row['SpectrumTitle'], frag_tech, proteases, instrument, rawfile_ID, -1, -1, -1)
        spectrum_match_commands[spectrum_hash] = 'MATCH (' + spectrum_hash + ':Spectrum {id:\'' + spectrum_hash + '\'})'

        # create or match sample
        if sample_ID not in sample_match_commands:
            query_str += ', ' + Neo4jCommands.create_sample_command(sample_ID, args.pride_accession, tissue_name, indiv_age, indiv_sex, phenotype)
            sample_match_commands[sample_ID] = 'MATCH (' + sample_ID + ':Sample {id:\'' + sample_ID + '\'})'
        else:
            query_str = sample_match_commands[sample_ID] + ' ' + query_str

        query_str += ', (' + spectrum_hash + ')-[:MEASURED_FROM]->(' + sample_ID + ')'
    else:
        psm_processed_count += 1
        # query_str = spectrum_match_commands[spectrum_hash] + ' ' + query_str
        print('[Warning]: Spectrum', spectrum_hash, 'matched to multiple peptides!')
        print()
        continue

    query_str += ', ' + Neo4jCommands.connect_peptide_spectrum_command(psm_row, 'pep_' + peptide_seq, spectrum_hash, psm_row['ID'], PSM_features_to_add)
    try:
        session.run(query_str)
        time.sleep(1)
    except:
        time.sleep(10)
        session.run(query_str)
        time.sleep(5)

    if (len(set_query_str) > 0):
        try:
            session.run(set_query_str)
            time.sleep(1)
        except:
            time.sleep(10)
            session.run(set_query_str)
            time.sleep(2)

    psm_processed_count += 1
    print('Processed:', psm_processed_count, '/', total_psm_count, end='\r')

peps_file.close()

# Create index on peptide sequence for faster search and close session
#session.run('CREATE TEXT INDEX pep_seq_index FOR (n:Peptide) ON (n.sequence)')
session.close()
print ('')
print ('Done')
