#!/usr/bin/env nextflow
// pipeline for doing the denovo assembly of yamadas 6 transciptomes generated from:
// 1.	N. agnita (free-living diatom, two bio-duplicates NA1 and 2)
// 2.	D. capensis (dinotom containing N. agnita, two bio-duplicates DC1 and 2);
// 3.	D. kwazulunatalensis (dinotom containing another diatom, two bio-duplicates DK1 and 2).

params.raw_reads_dir = "/home/humebc/projects/20201217_yamada/raw_seq_files"
output_dir = "${workflow.launchDir}/outputs"
bin_dir = "${workflow.launchDir}/bin"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"

params.subsample = true
params.subsample_depth = 10000

// Threads
params.trimmomatic_threads = 50
params.rcorrector_threads = 50
params.trinity_threads = 50

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
    storeDir "${params.subsample_depth}_subsampled_reads"

    input:
    set pair_id, file(reads) from ch_subsample

    output:
    tuple val(pair_id), file("*${pair_id}*_sub_${params.subsample_depth}*.fastq.gz") into ch_fastqc_pre_trim,ch_trimmomatic_input

    script:
	read_out_one = reads[0].getName().replaceAll(".fastq.gz", "_sub_${params.subsample_depth}.fastq.gz")
    read_out_two = reads[1].getName().replaceAll(".fastq.gz", "_sub_${params.subsample_depth}.fastq.gz")
    
    """
    seqtk sample -s100 ${reads[0]} ${params.subsample_depth} | gzip > ${read_out_one}
    seqtk sample -s100 ${reads[1]} ${params.subsample_depth} | gzip > ${read_out_two}
    """
}

ch_fastqc_pre_trim.view()

// process fastqc_pre_trim{
//     tag "${fastq_file}"
//     conda "envs/nf_general.yaml"
//     publishDir "nf_fastqc_pre_trim"

//     input:
//     file fastq_file from ch_fastqc_pre_trim.flatten()

//     output:
//     file "*_fastqc.html" into ch_fastqc_pre_trim_output

//     script:
//     """
//     fastqc -o . $fastq_file
//     """
// }

// /* 
//     Trim reads with trimmomatic
// */
// process trimmomatic{
//     cache 'lenient'
// 	tag "${fastq_file_one}"
//     conda "envs/nf_general.yaml"
// 	storeDir "nf_trimmed"

// 	input:
// 	tuple file(fastq_file_one), file(fastq_file_two) from ch_trimmomatic_input
	
// 	output:
// 	// Output that will be used for the error_correction
// 	// This is a list of tuples that are the 1P and 2P output files only
// 	tuple file("${fastq_file_one.getName().replaceAll('_1.fastq.gz', '')}*1P.fq.gz"), file("${fastq_file_one.getName().replaceAll('_1.fastq.gz', '')}*2P.fq.gz") into ch_rcorrect_input

// 	script:
// 	outbase = fastq_file_one.getName().replaceAll('_1.fastq.gz', '.trimmed.fq.gz')
// 	"""
// 	trimmomatic PE -threads ${params.trimmomatic_threads} -basein ${fastq_file_one} \\
// 		-baseout $outbase \\
// 		ILLUMINACLIP:${tru_seq_pe_fasta_path}:2:30:10:2:keepBothReads \\
// 		LEADING:3 TRAILING:3 MINLEN:36 HEADCROP:11
// 	"""
// }



// // Now error correction
// process rcorrector{
//     cache 'lenient'
//     tag "${trimmed_read_one}"
//     conda "envs/nf_general.yaml"
//     cpus params.rcorrector_threads
    
//     storeDir "nf_error_corrected"
    
//     input:
//     tuple file(trimmed_read_one), file(trimmed_read_two) from ch_rcorrect_input

//     output:
//     tuple file("${trimmed_read_one.getName().replaceAll('.trimmed_1P.fq.gz', '')}*1P.cor.fq.gz"), file("${trimmed_read_one.getName().replaceAll('.trimmed_1P.fq.gz', '')}*2P.cor.fq.gz") into ch_trinity_input

//     script:
//     """
//     run_rcorrector.pl -1 $trimmed_read_one -2 $trimmed_read_two -od . -t ${task.cpus}
//     """
// }




// // Now do the trinity assembly
// // NB I have changed the code in here subtly without re running it in order to implement the
// // storeDir directive. I have added a final mv line to the script to rename the ambiguous trinity fasta
// // file so that the name is specific to the sample SRR base name.
// // I have also changed the output dir from being a specific folder for each of the samples (due to the 
// // ambiguity in the previous naming system) to all being held in the nf_trinity_assembly dir.
// // I will make the changes by hand now to the outputs that are already in this directory.
// process trinity{
//     cache 'lenient'
//     tag "${corrected_read_one}"
//     conda "envs/nf_general.yaml"
//     cpus params.trinity_threads
//     storeDir "nf_trinity_assembly"

//     input:
//     tuple file(corrected_read_one), file(corrected_read_two) from ch_trinity_input

//     output:
//     file "${corrected_read_one.getName().replaceAll('.trimmed_1P.cor.fq.gz', '')}*Trinity.fasta" into ch_remove_short_iso_forms_trinity_input
    
//     script:
//     // NB that the output directory for each trinity assembly must have 'trinity' in it.
//     """
//     Trinity --left $corrected_read_one --right $corrected_read_two --seqType fq --max_memory 150G --CPU ${task.cpus} \\
//     --min_contig_length 250 --output trinity --full_cleanup
//     mv trinity.Trinity.fasta ${corrected_read_one.getName().replaceAll('.trimmed_1P.cor.fq.gz','')}.trinity.Trinity.fasta
//     """
// }

