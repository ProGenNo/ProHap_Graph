# Specify the configuration file below:
configfile: "config.yaml"

Ensembl_FTP_URL = "ftp.ensembl.org/pub/release-" + str(config['ensembl_release']) + "/"
annotationFilename = "Homo_sapiens.GRCh38." + str(config['ensembl_release']) + ".chr_patch_hapl_scaff"

sample_names_file = open(config['psm_file_list'], 'r')
lines = sample_names_file.readlines()
sample_names_file.close()
SAMPLES = [l[:-1] for l in lines[:]]

rule all:
    input:
        base="flag_baseDB_ready",
        psms=expand("tmp/{sample}_added", sample=SAMPLES),
        overview_stats="tmp/overview_stats"

rule download_reference_proteome:
    output:
        temp("data/fasta/Homo_sapiens.GRCh38.pep.all.fa")
    shell:
        "mkdir -p data/fasta ; "
        "wget " + Ensembl_FTP_URL + "fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O {output}.gz && gunzip {output}.gz; "

rule reference_filter_format:
    input:
        "data/fasta/Homo_sapiens.GRCh38.pep.all.fa"
    output:
        "data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa"
    conda: "condaenv.yaml"
    shell:
        "python3 src/fasta_format_headers.py -i {input} -o {output} -t _ensref -use_ENST 1 "

rule download_cdnas_fasta:
    output:
        out1=temp("data/fasta/Homo_sapiens.GRCh38.ncrna.fa"),
        out2=temp("data/fasta/Homo_sapiens.GRCh38.cdna.all.fa")
    shell:
        "mkdir -p data/fasta ; "
        "wget " + Ensembl_FTP_URL + "fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz -O {output.out1}.gz && gunzip {output.out1}.gz; "
        "wget " + Ensembl_FTP_URL + "fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz -O {output.out2}.gz && gunzip {output.out2}.gz; "

rule merge_cdnas_fasta:
	input:
		in1="data/fasta/Homo_sapiens.GRCh38.ncrna.fa",
		in2="data/fasta/Homo_sapiens.GRCh38.cdna.all.fa"
	output:
		temp("data/fasta/total_cdnas_" + str(config['ensembl_release']) + ".fa")
	shell:
		"cat {input.in1} > {output}; cat {input.in2} >> {output}"

rule download_gtf:
    output:
        temp("data/gtf/" + annotationFilename + ".gtf")
    shell:
        "mkdir -p data/gtf ; "
        "wget " + Ensembl_FTP_URL + "gtf/homo_sapiens/" + annotationFilename + ".gtf.gz -O {output}.gz && gunzip {output}.gz; "

rule parse_gtf_whole:
    input:
        "data/gtf/" + annotationFilename + ".gtf"
    output:
        temp("data/gtf/" + annotationFilename + ".db")
    conda: "condaenv.yaml"
    shell:
        "python src/parse_gtf.py -i {input} -o {output}"

rule build_database_graph:
    input:
        annot="data/gtf/" + annotationFilename + ".db",
        haplo_fasta=config['haplo_fasta'],
        ref_fasta="data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa",
        cdna_fasta="data/fasta/total_cdnas_" + str(config['ensembl_release']) + ".fa",
        haplo_table=config['haplo_table'],
        tr_ids='protein_transcript_ids_' + str(config['ensembl_release']) + '.csv',
        gene_ids='gene_transcript_ids_' + str(config['ensembl_release']) + '_full.csv',
    output:
        "flag_baseDB_ready"
    params:
        uri=config['neo4j_uri'],
        usr=config['neo4j_username'],
        pwd=config['neo4j_pwd']
    conda: "condaenv.yaml"
    shell:
        "python src/build_neo4j_graph.py -hap_fa {input.haplo_fasta} -ref_fa {input.ref_fasta} -cdna_fa {input.cdna_fasta} -hap_db {input.haplo_table} -annot_db {input.annot} -tr_id {input.tr_ids} -g_id {input.gene_ids} -uri {params.uri} -usr {params.usr} -pwd {params.pwd} && touch {output}"

rule fix_psm_report:
    input:
        psm=config['psm_file_location'] + "{sample}.tsv",
        haplo_table=config['haplo_table'],
        haplo_fasta=config['haplo_fasta'],
        ref_fasta="data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa",
        tr_ids='protein_transcript_ids_' + str(config['ensembl_release']) + '.csv',
        gene_ids='gene_transcript_ids_' + str(config['ensembl_release']) + '_full.csv',
    output:
        temp("tmp/{sample}_fixed.tsv")
    shell:
        "python src/psm_fix_annotation_newDB.py -i {input.psm} -o {output} -hap_tsv {input.haplo_table} -ref_fa {input.ref_fasta} -f {input.full_fasta} -tr_id {input.tr_ids} -g_id {input.gene_ids} -t {params.max_cores} -log /dev/null"

rule neo4j_add_psms:
    input:
        psm=config['psm_file_location'] + "{sample}.tsv",
        mf=config['meta_file'],
        haplo_table=config['haplo_table'],
        tr_ids='protein_transcript_ids_' + str(config['ensembl_release']) + '.csv',
        gene_ids='gene_transcript_ids_' + str(config['ensembl_release']) + '_full.csv',
        base_graph="flag_baseDB_ready"
    output:
        "tmp/{sample}_added"
    params:
        added_peps='tmp/already_added_peptides.txt',
        qval_thr=config['qval_thr'],
        rawfile_col=config['rawfile_col'],
        sample_col=config['sample_col'],
        frag_col=config['frag_col'],
        proteases_col=config['proteases_col'],
        instrument_col=config['instrument_col'],
        tissue_col=config['tissue_col'],
        age_col=config['age_col'],
        sex_col=config['sex_col'],
        phenotype_col=config['phenotype_col'],
        pride_acc=config['pride_acc'],
        uri=config['neo4j_uri'],
        usr=config['neo4j_username'],
        pwd=config['neo4j_pwd']
    conda: "condaenv.yaml"
    threads: 200
    shell:
        "mkdir -p tmp ; touch {params.added_peps} ; python src/neo4j_add_psms.py -psm {input.psm} -hap_tsv {input.haplo_table} -mf {input.mf} -tr_id {input.tr_ids} -g_id {input.gene_ids} -qval_thr {params.qval_thr} "
        "-sample_id {params.sample_col} -ID_col {params.rawfile_col} -frag_col {params.frag_col} -prot_col {params.proteases_col} -instr_col {params.instrument_col} -tissue_col {params.tissue_col} "
        "-age_col {params.age_col} -sex_col {params.sex_col} -pheno_col {params.phenotype_col} -pride_acc {params.pride_acc} -added_peps {params.added_peps} "
        "-uri {params.uri} -usr {params.usr} -pwd {params.pwd} && touch {output}"
    
rule get_overview_stats:
    input:
        psms=expand("tmp/{sample}_added", sample=SAMPLES),
        gene_ids='gene_transcript_ids_' + str(config['ensembl_release']) + '_full.csv',
        base_graph="flag_baseDB_ready"
    output:
        flag="tmp/overview_stats"
    params:
        uri=config['neo4j_uri'],
        usr=config['neo4j_username'],
        pwd=config['neo4j_pwd']
    conda: "condaenv.yaml"
    shell:
        "python src/refresh_overview_stats.py -g_id {input.gene_ids} -uri {params.uri} -usr {params.usr} -pwd {params.pwd} && touch {output.flag}"

