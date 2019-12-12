/*
 *This is a nextflow workflow to do a  first quality check of metagenomic datasets.
 */

/* 
 * pipeline input parameters 
 */

log.info """\
         METAGENOMICS - N F   P I P E L I N E    
         ===================================
         
         input - reads                  : ${params.reads}
         files in read set              : ${params.setsize}
         output - directory             : ${params.outdir}
         temporary - directory          : ${workDir}
         Trimmomatic adapters           : ${params.adapters} 
         Trimmomatic adapters directory : ${params.adapter_dir}
         phix - directory               : ${params.phix_dir} 
         host - directory               : ${params.host_dir}
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
    .set { read_pairs_ch } 
 
/* running trimmomatic to remove adapters sequences
 * $task.cpus to specify cpus to use
 */ 

process run_trim {
    conda 'configuration_files/trimmomatic_env.yml'
    publishDir "${params.outdir}/05_fastq_trimmed", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from read_pairs_ch

    output:
    set pair_id, file("${pair_id}*.trimmed.fq.gz") into reads_trimmed_ch
    file "${pair_id}_trimmed.log"

    """
    ${preCmd}
    trimmomatic PE -threads $task.cpus -trimlog ${pair_id}_trimmed.log ${pair_id}*.gz \
    -baseout ${pair_id}_trimmed.fq.gz ILLUMINACLIP:${params.adapter_dir}/${params.adapters}:${params.illuminaClipOptions} \
    SLIDINGWINDOW:${params.slidingwindow} \
    LEADING:${params.leading} TRAILING:${params.trailing} \
    MINLEN:${params.minlen} &> ${pair_id}_run.log
    mv ${pair_id}_trimmed_1P.fq.gz ${pair_id}_R1.trimmed.fq.gz
    mv ${pair_id}_trimmed_2P.fq.gz ${pair_id}_R2.trimmed.fq.gz
    cat ${pair_id}_trimmed_1U.fq.gz ${pair_id}_trimmed_2U.fq.gz > ${pair_id}_S_concat_stripped_trimmed.fq.gz
    """
}

/*
 * remove low-complexity reads from datasets with bbduk
 */
process run_low_complex {
    conda 'configuration_files/bbmap_env.yml'
    publishDir "${params.outdir}/06_bbduk_highC", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from reads_trimmed_ch

    output:
    set pair_id, file("${pair_id}*.trimmed.highC.fq.gz") into reads_highC_ch
    file "${pair_id}_bbduk_output.log"

    """
    ${preCmd}
    bbduk.sh threads=$task.cpus entropy=0.7 entropywindow=50 entropyk=5 \
    in1=${pair_id}_R1.trimmed.fq.gz \
    in2=${pair_id}_R2.trimmed.fq.gz \
    outm=${pair_id}.lowC.reads.fq.gz \
    out1=${pair_id}_R1.trimmed.highC.fq.gz \
    out2=${pair_id}_R2.trimmed.highC.fq.gz \
    stats=stats.txt &> ${pair_id}_bbduk_output.log
    """
}

/*
 * remove reads matching to phiX with bbduk
 */
process remove_phiX {
    conda 'configuration_files/bbmap_env.yml'
    publishDir "${params.outdir}/07_bbduk_phix", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from reads_highC_ch

    output:
    set pair_id, file("${pair_id}*.trimmed.highC.phix.fq.gz") into reads_phix_ch
    file "${pair_id}_bbduk_output.log"

    """
    ${preCmd}
    bbduk.sh threads=$task.cpus ref=${params.phix_dir}/${params.phix_file} k=31 hdist=1 \
    in1=${pair_id}_R1.trimmed.highC.fq.gz \
    in2=${pair_id}_R2.trimmed.highC.fq.gz\
    outm=${pair_id}.phix.reads.fq.gz \
    out1=${pair_id}.R1.trimmed.highC.phix.fq.gz \
    out2=${pair_id}.R2.trimmed.highC.phix.fq.gz \
    stats=stats.txt &> ${pair_id}_bbduk_output.log
    """
}

/*
 * remove reads matching to human genome with bbmap
 */

process remove_host {
    conda 'configuration_files/bbmap_env.yml'
    publishDir "${params.outdir}/08_bbmap_host", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from reads_phix_ch

    output:
    set pair_id, file("${pair_id}*.clean.fq.gz") into reads_host_ch
    file "${pair_id}.*.human.fq.gz"
    file "${pair_id}_bbmap_output.log"

    """
    ${preCmd}
    bbmap.sh -Xmx15g threads=6 \
    minid=0.95 maxindel=3 bwr=0.16 bw=12 \
    quickmatch fast minhits=2 \
    path=${params.host_dir} \
    in=${pair_id}.R1.trimmed.highC.phix.fq.gz\
    in2=${pair_id}.R2.trimmed.highC.phix.fq.gz \
    outu=${pair_id}.R1.clean.fq.gz \
    outu2=${pair_id}.R2.clean.fq.gz \
    outm=${pair_id}.R1.human.fq.gz \
    outm2=${pair_id}.R2.human.fq.gz \
    statsfile=${pair_id}.human_result.txt &> ${pair_id}_bbmap_output.log
    """
}

/*
 * Calculate the sequence coverage of the metagenomes
 */
process run_coverage {
    conda 'configuration_files/nonpareil_env.yml'
    publishDir "${params.outdir}/09_nonpareil", mode: "${params.savemode}"
    tag { pair_id }

    input:
    set pair_id, file(reads) from reads_host_ch

    output:
    file("${pair_id}*.npo") into r_plotting_ch
    file "${pair_id}*.npa"
    file "${pair_id}*.npc"
    file "${pair_id}*.npl"
    file "${pair_id}*.npo"
    
    

    """
    ${preCmd}
    gunzip -f *.fq.gz
    nonpareil -s *.R1.clean.fq -T kmer -f fastq -b ${pair_id}_R1 \
     -X ${params.query} -n ${params.subsample} -t $task.cpus
     sleep 10s
    """
}

/*
 * Create coverage calculations plots and combine into single html document
 */

 process plot_coverage {
    conda 'configuration_files/nonpareil_env.yml'
    publishDir "${params.outdir}/10_coverage_plots_clean_data", mode: "${params.savemode}"
    tag { "all samples" }

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