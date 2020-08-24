/*
 * This is a nextflow workflow to do a first quality check of metagenomic datasets.
 * The steps involve running fastqc, multiqc, and an analysis of sequencing depth using nonpareil.
 */

/* 
 * pipeline input parameters 
 */

log.info """\
         TALOS - a metagenomics shotgun pipeline    
         ===================================
         This current workflow generates a multiqc report
         and it calculates the coverage of the metagenomic samples.

         input - reads                  : ${params.reads}
         files in read set              : ${params.setsize}
         output - directory             : ${params.outdir}
         temporary - directory          : ${workDir} 
         """
         .stripIndent()

// Needed to run on the SAGA CLUSTER !!!! NO LONGER NEEDED
// preCmd = """
// if [ -f /cluster/bin/jobsetup ];
// then set +u; source /cluster/bin/jobsetup; set -u; fi
// """

// Creating the channels needed for the first analysis step
Channel 
    .fromFilePairs( params.reads, size:params.setsize, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch } 
 
 // first process is to run fastqc on the raw datasets
 // using the fastqc conda environment

process fastqc {
    conda 'conda_yml/fastqc_env.yml'

    publishDir "${params.outdir}/01_fastqc", mode: "${params.savemode}"
    
    tag "FASTQC on $sample_id"
    
    executor='slurm'
    label 'small'

    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_raw_ch


    script:
    """
    
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    
    """  
}  
 
// running multiqc on the fastqc files from the channel: fastqc_raw_ch

process multiqc {
    conda 'conda_yml/multiqc_env.yml'
    publishDir "${params.outdir}/02_multiqc", mode: "${params.savemode}"
    tag "multiqc on fastqc"
    
     executor='slurm'
    label 'small'

    input:
    file('*') from fastqc_raw_ch.collect()
    
    output:
    file('raw_data.multiqc_report.html')  
     
    script:
    """
    
    multiqc . 
    mv multiqc_report.html raw_data.multiqc_report.html
    
    """
} 

