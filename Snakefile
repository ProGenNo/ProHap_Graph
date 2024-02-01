# Specify the configuration file below:
configfile: "config_example.yaml"

Ensembl_FTP_URL = "ftp.ensembl.org/pub/release-" + str(config['ensembl_release']) + "/"
annotationFilename = "Homo_sapiens.GRCh38." + str(config['ensembl_release']) + ".chr_patch_hapl_scaff"

rule all:
    input:
        "flag_baseDB_ready"

rule download_reference_proteome:
    output:
        temp("data/fasta/Homo_sapiens.GRCh38.pep.all.fa")
    shell:
        "mkdir -p data/fasta ; "
        "wget " + Ensembl_FTP_URL + "fasta/homo_sapiens/pep/Homo_sapiens.GRCh38.pep.all.fa.gz -O {output}.gz && gunzip {output}.gz; "

rule reference_fix_headers:
    input:
        "data/fasta/Homo_sapiens.GRCh38.pep.all.fa"
    output:
        temp("data/fasta/ensembl_reference_proteinDB_" + str(config['ensembl_release']) + "_tagged.fa")
    conda: "condaenv.yaml"
    shell:
        "python src/fix_headers.py -i {input} -o {output} -t _ensref "

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
        tr_ids='protein_transcript_ids_' + config['ensembl_release'] + '.csv',
        gene_ids='gene_transcript_ids_' + config['ensembl_release'] + '.csv',
    output:
        "flag_baseDB_ready"
    params:
        uri=config['neo4j_uri'],
        usr=config['neo4j_username'],
        pwd=config['neo4j_pwd']
    conda: "condaenv.yaml"
    shell:
        "python src/build_neo4j_graph.py -hap_fa {input.haplo_fasta} -ref_fa {input.ref_fasta} -cdna_fa {input_cdna_fasta} -hap_db {input.haplo_db} -annot_db {input.annot} -tr_id {input.tr_ids} -g_id {input.gene_ids} -uri {params.uri} -usr {params.usr} -p {params.pwd} && touch {output}"

