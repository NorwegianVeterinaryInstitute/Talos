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

// Needed to run on the Abel cluster
preCmd = """
if [ -f /cluster/bin/jobsetup ];
then set +u; source /cluster/bin/jobsetup; set -u; fi
"""

// Creating the channels needed for the first analysis step
Channel 
    .fromFilePairs( params.reads, size:params.setsize, checkIfExists: true )
    .into { read_pairs_ch; read_pairs2_ch } 
 
 // first process is to run fastqc on the raw datasets
 // using the fastqc conda environment

process fastqc {
    conda 'configuration_files/fastqc_env.yml'

    publishDir "${params.outdir}/01_fastqc", mode: "copy"
    
    tag "FASTQC on $sample_id"

    input:
    set sample_id, file(reads) from read_pairs_ch

    output:
    file("fastqc_${sample_id}_logs") into fastqc_raw_ch


    script:
    """
    ${preCmd}
    mkdir fastqc_${sample_id}_logs
    fastqc -o fastqc_${sample_id}_logs -f fastq -q ${reads}
    """  
}  
 
// running multiqc on the fastqc files from the channel: fastqc_raw_ch

process multiqc {
    conda 'configuration_files/multiqc_env.yml'
    publishDir "${params.outdir}/02_multiqc", mode: "${params.savemode}"
       
    input:
    file('*') from fastqc_raw_ch.collect()
    
    output:
    file('raw_data.multiqc_report.html')  
     
    script:
    """
    ${preCmd}
    multiqc . 
    mv multiqc_report.html raw_data.multiqc_report.html
    """
} 

/*
 * Calculate the sequence coverage of the metagenomes
 */
process run_coverage {
    conda 'configuration_files/nonpareil_env.yml'
    publishDir "${params.outdir}/03_nonpareil_data", mode: "${params.savemode}"
    tag { sample_id }

    input:
    set sample_id, file(reads) from read_pairs2_ch

    output:
    file("${sample_id}*.npo") into r_plotting_ch
    file "${sample_id}*.npa"
    file "${sample_id}*.npc"
    file "${sample_id}*.npl"
    file "${sample_id}*.npo"
    
    

    """
    ${preCmd}
    gunzip -f *.gz
    nonpareil -s *.R1.* -T kmer -f fastq -b ${sample_id}_R1 \
     -X ${params.query} -n ${params.subsample} -t $task.cpus
     sleep 10s
    """
}

/*
 * Create coverage calculations plots and combine into single html document
 */

 process plot_coverage {
    conda 'configuration_files/nonpareil_env.yml'
    publishDir "${params.outdir}/04_coverage_plots_raw_data", mode: "${params.savemode}"
    tag { "All_samples" }

    input:
    file('*') from r_plotting_ch.collect()

    output:
    file "*.png"
    file "single_plots"   // folder with single file results

    """
    ${preCmd}
    mkdir single_plots
    Rscript $baseDir/Rscripts/process_npo_files.r
    """
}
