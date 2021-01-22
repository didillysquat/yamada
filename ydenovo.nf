#!/usr/bin/env nextflow
// pipeline for doing the denovo assembly of yamadas 6 transciptomes generated from:
// 1.	N. agnita (free-living diatom, two bio-duplicates NA1 and 2)
// 2.	D. capensis (dinotom containing N. agnita, two bio-duplicates DC1 and 2);
// 3.	D. kwazulunatalensis (dinotom containing another diatom, two bio-duplicates DK1 and 2).

params.raw_reads_dir = "/home/humebc/projects/20201217_yamada/raw_seq_files"
params.output_dir = "${workflow.launchDir}/outputs"
bin_dir = "${workflow.launchDir}/bin"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"

params.subsample = true
params.subsample_depth = 10000

// Publish dirs
if (params.subsample){
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
}else{
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
}


// Threads
params.trimmomatic_threads = 50
params.rcorrector_threads = 50
params.trinity_threads = 40

/* 
If subsample, create a channel that will pass into the subsampling process and then pass the subsampled
files into the trimming and fastqc
if not subsample, then create channel directly from the raw sequencing files
*/
if (params.subsample){
    Channel.fromFilePairs("${params.raw_reads_dir}/*_{1,2}.fastq.gz").set{ch_subsample}
}else{
    Channel.fromFilePairs("${params.raw_reads_dir}/*_{1,2}.fastq.gz").map{it[1]}.into{ch_fastqc_pre_trim; ch_trimmomatic_input}
}

process subsample{
    tag "${pair_id}"
    conda "envs/seqtk.yaml"
    publishDir "${params.output_dir}/${params.subsample_depth}_subsampled_reads"

    input:
    set pair_id, file(reads) from ch_subsample

    output:
    tuple val("${pair_id}_sub_${params.subsample_depth}"), file("${pair_id}_sub_${params.subsample_depth}*.fastq.gz") into ch_fastqc_pre_trim,ch_trimmomatic_input

    script:
	read_out_one = reads[0].getName().replaceAll("${pair_id}", "${pair_id}_sub_${params.subsample_depth}")
    read_out_two = reads[1].getName().replaceAll("${pair_id}", "${pair_id}_sub_${params.subsample_depth}")
    
    """
    seqtk sample -s100 ${reads[0]} ${params.subsample_depth} | gzip > ${read_out_one}
    seqtk sample -s100 ${reads[1]} ${params.subsample_depth} | gzip > ${read_out_two}
    """
}


process fastqc_pre_trim{
    tag "${fastq_file}"
    conda "envs/fastqc.yaml"
    publishDir fastqc_pre_trim_publish_dir

    input:
    file fastq_file from ch_fastqc_pre_trim.map{it[1]}.flatten()

    output:
    file "*_fastqc.html" into ch_fastqc_pre_trim_output

    script:
    """
    fastqc -o . $fastq_file
    """
}



/* 
    Trim reads with trimmomatic
*/
process trimmomatic{
    cache 'lenient'
	tag "${pair_id}"
    conda "envs/trimmomatic.yaml"
	
	input:
	tuple val(pair_id), file(fastqs) from ch_trimmomatic_input
	
	output:
	// Output that will be used for the error_correction
	// This is a list of tuples that are the 1P and 2P output files only
	tuple val(pair_id), file("${pair_id}*{1,2}P.fq.gz") into ch_rcorrect,ch_fastqc_post_trim

	script:
	outbase = fastqs[0].getName().replaceAll('_1.fastq.gz', '.trimmed.fq.gz')
	"""
	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${fastqs[0]} \\
		-baseout $outbase \\
		ILLUMINACLIP:${tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
	"""
}

process fastqc_post_trim{
    tag "${fastq_file}"
    conda "envs/fastqc.yaml"
    publishDir fastqc_post_trim_publish_dir

    input:
    file fastq_file from ch_fastqc_post_trim.map{it[1]}.flatten()

    output:
    file "*_fastqc.html" into ch_fastqc_post_trim_output

    script:
    """
    fastqc -o . $fastq_file
    """
}


// Now error correction
process rcorrector{
    cache 'lenient'
    tag "${trimmed_read_one}"
    conda "envs/rcorrector.yaml"
    cpus params.rcorrector_threads
    
    
    input:
    tuple val(pair_id), file(fastqs) from ch_rcorrect

    output:
    tuple val(pair_id), file("${pair_id}*{1,2}P.cor.fq.gz") into ch_trinity

    script:
    """
    run_rcorrector.pl -1 ${fastqs[0]} -2 ${fastqs[1]} -od . -t ${task.cpus}
    """
}


// Now do the trinity assembly
process trinity{
    cache 'lenient'
    tag "$pair_id"
    conda "envs/trinity.yaml"
    cpus params.trinity_threads
    storeDir "nf_trinity_assembly" 

    input:
    tuple val(pair_id), file(fastqs) from ch_trinity

    output:
    file "${pair_id}.Trinity.fasta" into ch_remove_short_iso_forms_trinity_input
    
    script:
    // NB that the output directory for each trinity assembly must have 'trinity' in it.
    """
    Trinity --left ${fastqs[0]} --right ${fastqs[1]} --seqType fq --max_memory 150G --CPU ${task.cpus} \\
    --min_contig_length 250 --output trinity --full_cleanup  --SS_lib_type RF
    mv trinity.Trinity.fasta ${pair_id}.Trinity.fasta
    """
}

//

