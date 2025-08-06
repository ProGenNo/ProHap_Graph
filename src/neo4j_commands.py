class Neo4jCommands:
    def create_gene_command(gene_feature):
        command_str = '(' + gene_feature.id + ':Gene { '
        command_str += 'id: \'' + gene_feature.id + '\', '
        try:
            command_str += 'name: \'' + gene_feature.attributes['gene_name'][0] + '\', '
        except:
            command_str += 'name: \'-\', '
        command_str += 'chrom: \'' + str(gene_feature.chrom) + '\', '
        command_str += 'version: \'' + gene_feature.attributes['gene_version'][0] + '\', '
        command_str += 'biotype: \'' + gene_feature.attributes['gene_biotype'][0] + '\', '
        command_str += 'strand: \'' + str(gene_feature.strand) + '\', '
        command_str += 'bp_from: ' + str(gene_feature.start) + ', '
        command_str += 'bp_to: ' + str(gene_feature.end) + ' })'
        return command_str

    def create_transcript_command(transcript_feature, annotations_db, cDNA):
        start_codons = [ str(sc.start) for sc in annotations_db.children(transcript_feature, featuretype='start_codon', order_by='start') ]    # there should be only one, but just in case...
        stop_codons = [ str(sc.start) for sc in annotations_db.children(transcript_feature, featuretype='stop_codon', order_by='start') ] 
        command_str = '(' + transcript_feature.id + ':Transcript { '
        command_str += 'id: \'' + transcript_feature.id + '\', '
        command_str += 'cDNA_sequence: \'' + cDNA + '\', '
        command_str += 'start: [' + ', '.join(start_codons) + '], '
        command_str += 'stop: [' + ', '.join(stop_codons) + '], '
        command_str += 'Ensembl_canonical: ' + str(int(('tag' in transcript_feature.attributes) and ('Ensembl_canonical' in transcript_feature.attributes['tag']))) + ', '
        command_str += 'MANE_select: ' + str(int(('tag' in transcript_feature.attributes) and ('MANE_Select' in transcript_feature.attributes['tag']))) + ', '
        command_str += 'version: \'' + transcript_feature.attributes['transcript_version'][0] + '\', '
        command_str += 'biotype: \'' + transcript_feature.attributes['transcript_biotype'][0] + '\' })'
        return command_str

    def create_exon_command(exon_feature):
        command_str = '(' + exon_feature.attributes['exon_id'][0] + ':Exon { '
        command_str += 'id: \'' + exon_feature.attributes['exon_id'][0] + '\', '
        command_str += 'bp_from: ' + str(exon_feature.start) + ', '
        command_str += 'bp_to: ' + str(exon_feature.end) + ' })'
        return command_str

    def create_haplotype_command(count, haplo_hash, geneID):
        command_str = '(' + haplo_hash + ':Haplotype { '
        command_str += 'id: \'' + geneID + '_hap_' + hex(count)[2:] + '\' })'
        return command_str

    def match_haplotype_command(count, haplo_hash, geneID):
        command_str = 'MATCH (' + haplo_hash + ':Haplotype { '
        command_str += 'id: \'' + geneID + '_hap_' + hex(count)[2:] + '\' })'
        return command_str

    def create_variant_command(changeID, changeIDsafe):
        command_str = '(' + changeIDsafe + ':Variant {'
        command_str += 'id: \'' + changeID + '\', '
        command_str += 'location: ' + changeID.split(':',2)[1] + ', '
        command_str += 'overlapping_peptide: FALSE, '
        command_str += 'ref: \'' + changeID.split(':', 2)[2].split('>',1)[0] + '\', '
        command_str += 'alt: \'' + changeID.split('>',1)[1] + '\' })'
        return command_str

    def create_proteoform_haplo_command(proteoformID, sequence, reading_frame, haplo_row, protein_changes):
        command_str = '(' + proteoformID + ':Proteoform {'
        command_str += 'id: \'' + proteoformID + '\', '
        command_str += 'length: ' + str(len(sequence)) + ', '
        command_str += 'sequence: \'' + sequence + '\', '
        command_str += 'reading_frame: ' + str(reading_frame) + ', '
        command_str += 'protein_changes: \'' + protein_changes + '\', '
        command_str += 'cDNA_changes: \'' + haplo_row['cDNA_changes'] + '\', '
        command_str += 'start_aa: ' + str(haplo_row['protein_prefix_length']) + ', '
        command_str += 'start_lost: ' + str(int(haplo_row['start_lost'])) + ', '
        command_str += 'splice_sites_affected: \'' + str(haplo_row['splice_sites_affected']) + '\' })'
        return command_str

    def create_proteoform_reference_command(proteoformID, fasta_element):
        command_str = '(' + proteoformID + ':Proteoform {'
        command_str += 'id: \'' + proteoformID + '\', '
        command_str += 'length: ' + str(len(fasta_element['sequence'])) + ', '
        command_str += 'sequence: \'' + fasta_element['sequence'] + '\', '
        command_str += 'reading_frame: -1, '
        command_str += 'protein_changes: \'-\', '
        command_str += 'cDNA_changes: \'-\', '
        command_str += 'start_aa: 0, '
        command_str += 'start_lost: 0, '
        command_str += 'splice_sites_affected: \'-\' })'
        return command_str

    def create_peptide_command(pep_id, sequence, changes, class_1, class_2, freq, contam_matches):
        command_str = '(' + pep_id + ':Peptide {'
        command_str += 'id: \'' + pep_id + '\', '
        command_str += 'length: ' + str(len(sequence)) + ', '
        command_str += 'sequence: \'' + sequence + '\', '
        command_str += 'pep_changes: \'' + changes + '\', '
        command_str += 'expected_max_freq: \'' + freq + '\', '
        command_str += 'contaminant_matches: \'' + contam_matches + '\', '
        command_str += 'pep_class_1: \'' + class_1 + '\', '
        command_str += 'pep_class_2: \'' + class_2 + '\' })'
        return command_str

    def create_spectrum_command(spec_id, spec_title, frag_tech, proteases, instrument, fraction_ID, precursor_mz, precursor_intens, retention_time):
        command_str = '(' + spec_id + ':Spectrum {'
        command_str += 'id: \'' + spec_id + '\', '
        command_str += 'title: \'' + spec_title + '\', '
        command_str += 'frag_technique: \'' + frag_tech + '\', '
        command_str += 'proteases: \'' + proteases + '\', '
        command_str += 'spectrometer: \'' + instrument + '\', '
        command_str += 'fraction_id: \'' + fraction_ID + '\', '
        command_str += 'retention_time: ' + str(retention_time) + ', '
        command_str += 'precursor_mz: ' + str(precursor_mz) + ', '
        command_str += 'precursor_intensity: ' + str(precursor_intens) + ' })'
        return command_str

    def create_sample_command(sample_id, pride_id, tissue, age, sex, phenotype):
        command_str = '(' + sample_id + ':Sample {'
        command_str += 'id: \'' + sample_id + '\', '
        command_str += 'pride_project_accession: \'' + pride_id + '\', '
        command_str += 'tissue_name: \'' + tissue + '\', '
        command_str += 'phenotype: \'' + phenotype + '\', '
        command_str += 'individual_age: \'' + str(age) + '\', '
        command_str += 'individual_sex: \'' + sex + '\' })'
        return command_str

    def connect_peptide_spectrum_command(psm_row, pep_id, spec_id, psm_id, features):
        command_str = "(" + pep_id + ')-[:MATCHED_TO {'
        command_str += 'id: \'' + psm_id + '\', '
        command_str += ', '.join([ (feature.replace('-', '_') + ': ' + str(psm_row[feature])) for feature in features ]) + '}'
        command_str += ']->(' + spec_id + ')'
        return command_str

