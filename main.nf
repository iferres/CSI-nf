#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.in = "*.fasta"
params.outdir = "results"
params.blastdb = "path/to/nt"
//params.blastn_name = "nt"

process CRISPRCASFINDER {
    label 'crisprcasfinder'

    input:
    tuple val(genome_id), path("${genome_id}.fasta")

    output:
    tuple val(genome_id), path("Result*/GFF/*.gff")

    shell:
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
        path("${genome_id}_${contig_id}_sp.fasta"), 
        emit: fasta,
        optional: true

    shell:
    """
    #!/usr/bin/env Rscript

    gff <- try(read.table("${genome_id}_${contig_id}.gff"))

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


        fasta <- paste0(">", paste(spacers\$V12, spacers\$V10, sep = "\\n"), collapse = "\\n")
        cat(fasta,
            file = "${genome_id}_${contig_id}_sp.fasta")
    }
    """
}

 db_name = file(params.blastdb).name
 db_path = file(params.blastdb).parent

process BLASTN {

    input:
    tuple val(genome_id),
        val(contig_id),
        path("${genome_id}_${contig_id}_sp.fasta")

    output:
    tuple val(genome_id),
        val(contig_id),
        path("${genome_id}_${contig_id}_sp_vs_${db_name}.tab")

    shell:
    """
    blastn -task "megablast" \
        -evalue 1e-5 \
        -num_threads 1 \
        -db ${db_path}/${db_name} \
        -query ${genome_id}_${contig_id}_sp.fasta \
        -out ${genome_id}_${contig_id}_sp_vs_${db_name}.tab \
        -outfmt 6
    """
}


process PARSE_BLASTN {

    input:
    tuple val(genome_id),
        val(contig_id),
        path("${genome_id}_${contig_id}_sp_vs_${db_name}.tab")

    output:
    //tuple val(genome_id),
    //    val(contig_id),
    //    path("${genome_id}_${contig_id}_sp_vs_${db_name}_Filt.tab")
    tuple val(genome_id),
        val(contig_id),
        path("${genome_id}_${contig_id}_sp_vs_${db_name}_GIDS.txt")

    shell:
    """
    awk '{print \$2} ' ${genome_id}_${contig_id}_sp_vs_${db_name}.tab > \
        ${genome_id}_${contig_id}_sp_vs_${db_name}_GIDS.txt
    """
}

process LCA {
    label 'taxonkit'

    input:
    tuple val(genome_id),
        val(contig_id),
        path("${genome_id}_${contig_id}_sp_vs_${db_name}_GIDS.txt")

    shell:
    """
    cat ${genome_id}_${contig_id}_sp_vs_${db_name}_GIDS.txt | \
        taxonkit --data-dir ${db_path} 
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

    BLASTN( EXTRACT_SPACERS.out.fasta )

    PARSE_BLASTN( BLASTN.out )

    LCA( PARSE_BLASTN.out )
}