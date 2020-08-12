/*
 * This is a nextflow workflow to set up the conda environments for the the Talos workflow
 * The environments will be created in the work directory in a seperate folder called: conda.
 * each environment is tested by calling the version number of the main tool installed in the 
 * environment.
 */

/* 
 * pipeline input parameters 
 */

log.info """\
         METAGENOMICS - N F   P I P E L I N E    
         ===================================
         
         temporary - directory          : ${workDir}
         """
         .stripIndent()

/* 
 * process names  
 */

process fastqc {
    conda 'conda_yml/fastqc_env.yml'
    
    label 'small'
    
    """
    fastqc -v
    """
}

process multiqc {
    conda 'conda_yml/multiqc_env.yml'
    
    label 'small'
    
    """
    multiqc --version
    """
}

process run_coverage {
    conda 'conda_yml/nonpareil_env.yml'
    
    label 'small'
    
    """
    nonpareil -V
    """
}

process run_trim {
    conda 'conda_yml/trimmomatic_env.yml'
    
    label 'small'
    
    """
    trimmomatic -version
    """
}

process run_low_complex {
    conda 'conda_yml/bbmap_env.yml'
    
    label 'small'
    
    """
    bbmap.sh -version
    bbduk.sh -version
    """
}

process Average_gsize {
    conda 'conda_yml/microbecensus_env.yml'
    
    label 'small'
    
    """
    run_microbe_census.py --version
    """
}

process hulk_calculation {
    conda 'conda_yml/hulk_env.yml'
    
    label 'small'
    
    """
    hulk version
    """

}

process Kraken_classification {
    conda 'conda_yml/kraken2_env.yml'
    
    label 'small'
    
    """
    kraken2 --version
    """
}