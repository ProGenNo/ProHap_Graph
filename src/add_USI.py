import argparse
import pandas as pd
import xmltodict

parser = argparse.ArgumentParser(description='Adds the Universal Spectrum Identifier column')

parser.add_argument("-psm", dest="psm_file", required=True,
                    help="PSM TSV file")

parser.add_argument("-qval_col", dest="qval_column", required=False,
                    help="q-value column name in the PSM table", default="q-value")

parser.add_argument("-qval_thr", dest="qval_threshold", required=False, type=float,
                    help="maximum acceptable q-value, default: 0.01", default=0.01)

parser.add_argument('-pride_acc', dest="pride_accession", required=False,
                    help="accession number of the Pride dataset")

parser.add_argument("-o", dest="output_file", required=True,
                    help="putput TSV file")

args = parser.parse_args()

mods = {}

with open('data/unimod.xml','r') as file:
    for mod in xmltodict.parse(file.read())['umod:unimod']['umod:modifications']['umod:mod']:
        mods[mod['@record_id']] = mod

def replace_mods(pep_str):
    split_str = pep_str.split('[')
    result_str = split_str[0]

    for substr in split_str[1:]:
        mod_id = substr.split(']',1)[0]
        mod_mass_shift = float(mods[mod_id]['umod:delta']['@mono_mass'])
        result_str += '[' + ('+' if not mod_mass_shift.startswith('-') else '') + "{:.4f}".format(mod_mass_shift) + ']' + substr.split(']',1)[1]

    return result_str

print ("Reading", args.psm_file)
psm_df = pd.read_csv(args.psm_file, header=0, sep='\t')

# Apply FDR threshold and remove contaminants and decoys
psm_df = psm_df[(psm_df[args.qval_column] <= args.qval_threshold) & (~psm_df['psm_type1'].isin(['decoy', 'contaminant']))]  #  & (psm_df['matching_proteins'].str.count(';') < 500)

psm_df['USI'] = psm_df.apply(lambda row: 'mzspec:' + args.pride_accession + ':' + row['SpectrumFilename'] + ':' + row['SpectrumTitle'].split('scan=',1)[1] + '/' + ('2' if row['charge_2'] else ('3' if row['charge_3'] else '4')), axis=1)

psm_df.to_csv(args.output_file, sep='\t', index=False)