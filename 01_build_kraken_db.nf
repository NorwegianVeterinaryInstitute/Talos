/*
 * This is a nextflow workflow to build a kraken database using the following databases
 * Archaea
 * Bacteria
 * plasmids
 * viral sequences
 * human
 * fungi
 * protozoa
 */

/* 
 * pipeline input parameters 
 */

log.info """\
         METAGENOMICS - N F   P I P E L I N E    
         ===================================
         
         temporary - directory          : ${workDir}
         location of kraken2 directory  : ${params.kraken2.path}
         kraken2 database               : ${params.kraken2_dir}
         """
         .stripIndent()

/* 
 * process names  
 */



// list of taxa to download to build the database
//taxons = Channel.from('archaea', 'plasmid','bacteria', 'viral','human','fungi', 'protozoa')
taxons = Channel.value('viral')

// downloading the taxonomy
process download_taxonomy {
    conda 'conda_yml/kraken2_env.yml'
    publishDir "${params.kraken2.path}", mode: "${params.savemode}"

    executor='local'
    
    output:
    file ("${params.kraken2_dir}/taxonomy") into taxonomy_ch

    """          
    kraken2-build --download-taxonomy --threads 1 --db ${params.kraken2_dir}
    
    """
}

// downloading the taxa for the database
process download_taxa {
    conda 'conda_yml/kraken2_env.yml'
    //publishDir "${params.kraken2.path}", mode: "${params.savemode}"
    tag "$taxa"

    executor='local'


    input:
    val taxa from taxons


    output:
    file ("${params.kraken2_dir}/library/$taxa") into downloads_ch

    """
    kraken2-build --download-library $taxa --threads 2 --no-masking --db ${params.kraken2_dir}
    
    """
}

// masking low-complexity sequencing, overwriting the previous unmasked dataset.
process masking_taxa {
    conda 'conda_yml/kraken2_env.yml'
    publishDir "${params.kraken2.path}", mode: "${params.savemode}"
    tag "$taxa"


    executor='slurm'
    label 'small'

    input:
    val taxa from taxons
    file ("${params.kraken2_dir}/library/$taxa") from downloads_ch


    output:
    file ("${params.kraken2_dir}/library/$taxa") into database_ch
    

    """
    cd ${params.kraken2_dir}/library/$taxa
    ls -lath 
    
    #dustmasker commands
    dustmasker -in library.fna -outfmt fasta | sed -e '/^>/!s/[a-z]/x/g' > library.fna.tmp
    mv library.fna.tmp library.fna
    touch library.fna.masked
    cd -
    
    """
}




process build_Kraken2_db {
    conda 'conda_yml/kraken2_env.yml'
    
    executor='slurm'
    label 'large'
    label 'longtime'

    input:
    val taxa from taxons
    file ("${params.kraken2_dir}/taxonomy") from taxonomy_ch
    file ("${params.kraken2_dir}/library/$taxa") from database_ch.collect()

    """
    echo This is working

    ls ${params.kraken2.path}/${params.kraken2_dir}
    ls ${params.kraken2.path}/${params.kraken2_dir}/library
    kraken2-build build --threads=8 --db ${params.kraken2.path}/${params.kraken2_dir}
    #kraken2-build build --threads=24 --db ${params.kraken2.path}/${params.kraken2_dir}

    """
}

