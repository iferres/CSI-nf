#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.in = "*.fasta"
params.outdir = "results"
params.blastdb = "/mnt/cive/ncbi/nt"
//params.blastn_name = "nt"

process CRISPRCASFINDER {
    label 'crisprcasfinder'

    input:
    tuple val(genome_id), path("${genome_id}.fasta")

    output:
    tuple val(genome_id), path("Result*/GFF/*.gff")

    script:
    """
    crispr-entrypoint.sh ${genome_id}.fasta
    """
}

process EXTRACT_SPACERS {
    publishDir "SPACERS", mode: 'copy'

    input:
    tuple val(genome_id), 
        val(contig_id), 
        path("${genome_id}_${contig_id}.gff")

    output:
    tuple val(genome_id), 
        val(contig_id), 
        path("${genome_id}_${contig_id}_sp.gff"), 
        emit: gff,
        optional: true
    tuple val(genome_id),
        val(contig_id),
        path("fasta_spacers/*.fasta"), 
        emit: fasta,
        optional: true

    script:
    """
    #!/usr/bin/env Rscript

    gff <- try(read.table("${genome_id}_${contig_id}.gff", stringsAsFactors=FALSE))

    if (class(gff)!="try-error"){
        spacers <- gff[gff\$V3=="CRISPRspacer",,drop=FALSE]
    
        spl <- strsplit(spacers\$V9, ";")
    
        spacers\$V10 <- sapply(spl, function(x){
            gp <- grep("sequence=", x, value = TRUE)
            sq <- sub("sequence=", "", gp, fixed = TRUE)
            sq
        })
    
        spacers\$V11 <- sapply(spl, function(x) {
            gp <- grep("ID=", x, value=TRUE)
            id <- sub("ID=", "", gp, fixed = TRUE)
            id
        })
 
        spacers\$V12 <- paste("${genome_id}", "${contig_id}", spacers\$V11, sep = "_")

        write.table(spacers, 
            file = "${genome_id}_${contig_id}_sp.gff", 
            quote = FALSE, 
            sep = "\\t", 
            row.names = FALSE, 
            col.names = FALSE)

        fastas <- paste(paste0(">", spacers\$V12), spacers\$V10, sep = "\\n")
        
        dir.create("fasta_spacers")
        for (i in 1:length(fastas)){
            fi <- paste0("fasta_spacers/", spacers\$V11[i], ".fasta")
            cat(fastas[i], file = fi)
        }
        
    }
    """
}

 db_name = file(params.blastdb).name
 db_path = file(params.blastdb).parent

process BLASTN {

    input:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        path("${genome_id}_${contig_id}_${spacer_id}.fasta")

    output:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}.tab")

    script:
    """
    blastn -task "megablast" \
        -evalue 1e-5 \
        -num_threads 1 \
        -db ${db_path}/${db_name} \
        -query ${genome_id}_${contig_id}_${spacer_id}.fasta \
        -out ${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}.tab \
        -outfmt 6
    """
}


process PARSE_BLASTN {

    input:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}.tab")

    output:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_GIDS.txt")

    script:
    """
    awk '{print \$2} ' ${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}.tab > \
        ${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_GIDS.txt
    """
}

process GET_TAXID {
    input:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_GIDS.txt")
    
    output:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        stdout, optional: true
        //path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_TAXIDS.txt")

    script:
    """
    blastdbcmd \
        -db ${db_path}/${db_name} \
        -outfmt "%T" \
        -entry_batch ${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_GIDS.txt # > \
        # ${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_TAXIDS.txt
    """
}

process LCA {
    label 'taxonkit'

    input:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        val(stdin)
        //path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_TAXIDS.txt")

    output:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        stdout
        //path("${genome_id}_${contig_id}_${spacer_id}_vs_${db_name}_LCA_RANK.txt")

    when:
    stdin != ""
    
    script:
    """
    echo "${stdin}" | tr '\\n' '\\t' | \
        taxonkit --data-dir ${db_path} lca | \
        awk '{print \$NF}' | \
        taxonkit --data-dir ${db_path} lineage | \
        taxonkit --data-dir ${db_path} reformat | \
        cut -f 3
    """
}

process WRITE_TSV {
    input:
    tuple val(genome_id),
        val(contig_id),
        val(spacer_id),
        val(lca_ranks)

    output:
    path("${genome_id}_${contig_id}_${spacer_id}.tsv")

    script:
    """
    echo -ne "${genome_id}\\t${contig_id}\\t${spacer_id}\\t${lca_ranks}"  \
        > ${genome_id}_${contig_id}_${spacer_id}.tsv
    """
}


process COLLECT_RESULTS {
    publishDir "LCA_RANKS", mode: 'copy'

    input:
    path("*.tsv")
    
    output:
    path("ranks.tsv")

    script:
    """
    echo -ne "genome\\tcontig\\tspacer\\trank\\n" > ranks
    cat *.tsv >> ranks
    mv ranks ranks.tsv
    """
}

Channel
  .fromPath(params.in, checkIfExists: true)
  .ifEmpty { exit 1, "Non fastq files found: ${params.in}" }
  .map{ file -> tuple(file.baseName, file)}
  .set{ in_ch }


workflow {
    
    CRISPRCASFINDER( in_ch )
    CRISPRCASFINDER.out
        .set{ crout_ch }

    crout_ch
        .transpose()
        .map{it -> tuple(it[0], it[1].baseName, it[1])}
        .set{ crout_contig_ch }

    EXTRACT_SPACERS( crout_contig_ch )

    EXTRACT_SPACERS.out.fasta
        .transpose()
        .map{it -> tuple(it[0], it[1], it[2].baseName, it[2])}
        .set{ spacer_ch }
        
        //spacer_ch.view()

    BLASTN( spacer_ch )

    PARSE_BLASTN( BLASTN.out )

    GET_TAXID( PARSE_BLASTN.out )

    LCA( GET_TAXID.out )

    LCA.out
        .set{ result_ranks }

    WRITE_TSV( result_ranks )

    WRITE_TSV.out
        .collect()
        .set{ collected_ch }

    //collected_ch.view()

    COLLECT_RESULTS( collected_ch )
}