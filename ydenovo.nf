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

// Evalue for taxonomic analyses
// NB there is currently a bug in mmseqs where by they are using a FLOAT to store
// the evalue number. As such the eval will break if you go much above 1e-35.
params.tax_evalue = '1e-35'
// Publish dirs
// TODO consider making the differentiatl publishDir via the saveAs option of publishDir
if (params.subsample){
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_correct_publish_dir = [params.output_dir, "fastqc_post_correct_sub_${params.subsample_depth}"].join(File.separator)
    mmseqs_read_query_db_store_dir = [params.output_dir, "mmseqs_read_query_dbs_sub_${params.subsample_depth}"].join(File.separator)
    mmseqs_taxonomy_publish_dir = [params.output_dir, "mmseqs_taxonomy_sub_${params.subsample_depth}_${params.tax_evalue}"].join(File.separator)
    trinity_assembly_publish_dir = [params.output_dir, "trinity_assembly_sub_${params.subsample_depth}"].join(File.separator)
    trinity_stats_publish_dir = [params.output_dir, "trinity_basic_stats_sub_${params.subsample_depth}"].join(File.separator)
    busco_stats_publish_dir = [params.output_dir, "busco_stats_sub_${params.subsample_depth}"].join(File.separator)
    assembly_bowtie2_indexed_db_storeDir = [params.output_dir, "assembly_bowtie2_indexed_dbs_sub_${params.subsample_depth}"].join(File.separator)
    read_mapping_stats_publish_dir = [params.output_dir, "read_mapping_stats_sub_${params.subsample_depth}"].join(File.separator)
}else{
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
    fastqc_post_correct_publish_dir = [params.output_dir, "fastqc_post_correct"].join(File.separator)
    mmseqs_read_query_db_store_dir = [params.output_dir, "mmseqs_read_query_dbs"].join(File.separator)
    mmseqs_taxonomy_publish_dir = [params.output_dir, "mmseqs_taxonomy_${params.tax_evalue}"].join(File.separator)
    trinity_assembly_publish_dir = [params.output_dir, "trinity_assembly"].join(File.separator)
    trinity_stats_publish_dir = [params.output_dir, "trinity_basic_stats"].join(File.separator)
    busco_stats_publish_dir = [params.output_dir, "busco_stats"].join(File.separator)
    assembly_bowtie2_indexed_db_storeDir = [params.output_dir, "assembly_bowtie2_indexed_dbs"].join(File.separator)
    read_mapping_stats_publish_dir = [params.output_dir, "read_mapping_stats"].join(File.separator)
}

// DB paths
params.mmseqs_nt_path = '/home/humebc/nt/nt.fnaDB'
params.path_to_bowtie2_silva_db = '/home/humebc/silva/SILVA_SSUParc_LSUParc_tax'

// Threads
params.trimmomatic_threads = 50
params.rcorrector_threads = 50
params.trinity_threads = 50
params.mmseqs_threads = 50
params.busco_threads = 50
params.bowtie2_threads = 50
params.rsem_threads = 50

/* 
If subsample, create a channel that will pass into the subsampling process and then pass the subsampled
files into the trimming and fastqc
if not subsample, then create channel directly from the raw sequencing files
*/
if (params.subsample){
    Channel.fromFilePairs("${params.raw_reads_dir}/*_{1,2}.fastq.gz").set{ch_subsample}
    
    process subsample{
        tag "${pair_id}"
        conda "envs/seqtk.yaml"
        publishDir "${params.output_dir}/${params.subsample_depth}_subsampled_reads", mode: 'copy'
        storeDir "${params.output_dir}/${params.subsample_depth}_subsampled_reads"

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
}else{
    Channel.fromFilePairs("${params.raw_reads_dir}/*_{1,2}.fastq.gz").into{ch_fastqc_pre_trim; ch_trimmomatic_input}
}

process fastqc_pre_trim{
    tag "${fastq_file}"
    conda "envs/fastqc.yaml"
    publishDir fastqc_pre_trim_publish_dir, mode: 'copy'
    storeDir fastqc_pre_trim_publish_dir

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_pre_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${pair_id}.fastqc.html" into ch_fastqc_pre_trim_output

    script:
    """
    fastqc -o . $fastq_file
    mv *.html ${pair_id}.fastqc.html
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
    publishDir fastqc_post_trim_publish_dir, mode: 'copy'
    storeDir fastqc_post_trim_publish_dir

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_post_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${pair_id}.fastqc.html" into ch_fastqc_post_trim_output

    script:
    """
    fastqc -o . $fastq_file
    mv *.html ${pair_id}.fastqc.html
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
    tuple val(pair_id), file("${pair_id}*{1,2}P.cor.fq.gz") into ch_bowtie2_silva_mapping

    script:
    """
    run_rcorrector.pl -1 ${fastqs[0]} -2 ${fastqs[1]} -od . -t ${task.cpus}
    """
}

// // Discard any read pairs where one or both reads map to the SILVA rRNA database
// // per https://informatics.fas.harvard.edu/best-practices-for-de-novo-transcriptome-assembly-with-trinity.html
// // Out put the metrics of the mapping
// // We will conservatively only keep those samples that are unaligned in either read to the SILVA db
process bowtie2_silva_mapping{
    tag "${fastq_file}"
    conda "envs/bowtie2.yaml"
    cpus params.rcorrector_threads

    input:
    tuple val(pair_id), file(fastqs) from ch_bowtie2_silva_mapping

    output:
    tuple val(pair_id), file("${pair_id}_paired_unaligned_{R1,R2}.fq.gz") into ch_trinity,ch_fastqc_post_correct,ch_mmseq_create_query_dbs,ch_bowtie2_mapping_stats_reads,ch_rsem_reads

    script:
    metric_file = "${pair_id}.bowtie_metrics"
    """
    bowtie2 --quiet --very-sensitive-local --phred33  -x ${params.path_to_bowtie2_silva_db} -1 ${fastqs[0]} -2 ${fastqs[1]} --threads ${task.cpus} \\
    --al-conc-gz ${pair_id}_paired_aligned.fq.gz \\
    --un-conc-gz ${pair_id}_paired_unaligned.fq.gz  --al-gz ${pair_id}_unpaired_aligned.fq.gz \\
    --un-gz ${pair_id}_unpaired_unaligned.fq.gz
    mv ${pair_id}_paired_unaligned.fq.1.gz ${pair_id}_paired_unaligned_R1.fq.gz
    mv ${pair_id}_paired_unaligned.fq.2.gz ${pair_id}_paired_unaligned_R2.fq.gz
    """
}

process fastqc_post_correct{
    tag "${fastq_file}"
    conda "envs/fastqc.yaml"
    publishDir fastqc_post_correct_publish_dir, mode: 'copy'
    storeDir fastqc_post_correct_publish_dir

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_post_correct.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${pair_id}.fastqc.html" into ch_fastqc_post_correct_output

    script:
    """
    fastqc -o . $fastq_file
    mv *.html ${pair_id}.fastqc.html
    """
}

/* At this point we want to implement some taxonomic checks to see that we are working
   with the orgnisms that we expect
   An mmseqs taxonomy database has been create as per section "Taxonomy assignment" of:
   https://mmseqs.com/latest/userguide.pdf
   We will run an mmseqs query against this database for each of the seq files
   We will break this up into the following processes
   1 - mmseqs creatdb --> queryDB
   2 - mmseqs taxonomy --> taxonomyResult
        We will limit matches to very confident matches with evalues below 1e-50 
        per http://www.metagenomics.wiki/tools/blast/evalue#:~:text=Blast%20hits%20with%20an%20E,good%20hit%20for%20homology%20matches.
   3 - mmseqs createtsv --> taxonomyResult.tsv
   3 - mmseqs taxonomyreport --> taxonomyResult_report (--report-mode 1 can be used)
        We will produce the standard taxonomyResult_report and the html krona report.    
   We will want to process on a per readfile basis
   see if we can explicitly code in the files
*/
process mmseq_create_query_dbs{
    cache 'lenient'
    tag "$fastq_file"
    container 'soedinglab/mmseqs2:latest'
    containerOptions '-u $(id -u):$(id -g)'
    storeDir mmseqs_read_query_db_store_dir

    input:
    tuple val(pair_id), file(fastq_file) from ch_mmseq_create_query_dbs.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    tuple val(pair_id), file("${pair_id}.queryDB{,.dbtype,_h,_h.dbtype,_h.index,.index,.lookup,.source}") into ch_mmseq_taxonomy

    script:
    out_name = "${pair_id}.queryDB"
    """
    mmseqs createdb $fastq_file $out_name --dbtype 2
    """
}

// NB there is currently a bug in mmseqs where by they are using a FLOAT to store
// the evalue number. As such the eval will break if you go much above 1e-35.
// process mmseq_taxonomy{
//     cache 'lenient'
//     tag "${db_files[0]}"
//     container 'soedinglab/mmseqs2:latest'
//     containerOptions '-u $(id -u):$(id -g)'
//     cpus params.mmseqs_threads
//     publishDir mmseqs_taxonomy_publish_dir, mode: 'copy'
//     storeDir mmseqs_taxonomy_publish_dir
//     stageInMode 'copy'

//     input:
//     tuple val(pair_id), file(db_files) from ch_mmseq_taxonomy

//     output:
//     tuple file("${pair_id}.taxonomyResult.tsv"),\
//     file("${pair_id}.taxonomyResult_report"),\
//     file("${pair_id}.report.html") into ch_mmseq_taxonomy_out

//     script:
//     output_tsv = "${pair_id}.taxonomyResult.tsv"
//     output_report = "${pair_id}.taxonomyResult_report"
//     output_html = "${pair_id}.report.html"
//     """
//     mmseqs taxonomy  ${pair_id}.queryDB ${params.mmseqs_nt_path} taxonomyResult tmp --threads ${task.cpus} -e ${params.tax_evalue}
//     mmseqs createtsv ${pair_id}.queryDB taxonomyResult $output_tsv
//     mmseqs taxonomyreport ${params.mmseqs_nt_path} taxonomyResult $output_report
//     mmseqs taxonomyreport ${params.mmseqs_nt_path} taxonomyResult $output_html --report-mode 1
//     """
// }

// //TODO do two kmer sizes for trinity. Check what the default is but look at
// // possibly using 25 and 32 based on
// // https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0232005
// // Also, spadesrna seems like a smart choice for a second assembler.
// // Now do the trinity assembly. We want to export the assembly and the gene to transcript map
process trinity{
    cache 'lenient'
    tag "$pair_id"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    cpus params.trinity_threads
    storeDir trinity_assembly_publish_dir

    input:
    tuple val(pair_id), file(fastqs) from ch_trinity

    output:
    tuple val(pair_id), file("${pair_id}.Trinity.fasta") into ch_trinity_stats,ch_busco,ch_assembly_bowtie_db,ch_rsem_assemblies
    tuple val(pair_id), file("${pair_id}.gene_trans_map") into ch_trin_gene_trans_map

    script:
    // NB that the output directory for each trinity assembly must have 'trinity' in it.
    """
    Trinity --left ${fastqs[0]} --right ${fastqs[1]} --seqType fq --max_memory 150G --CPU ${task.cpus} \\
    --min_contig_length 250 --output trinity --full_cleanup  --SS_lib_type RF
    mv trinity.Trinity.fasta ${pair_id}.Trinity.fasta
    mv trinity.Trinity.fasta.gene_trans_map ${pair_id}.gene_trans_map
    """
}

// // TODO first do some basic stats using the Trinity stats package.
process trinity_stats{
    cache 'lenient'
    tag "$pair_id"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir trinity_stats_publish_dir

    input:
    tuple val(pair_id), file(trin_assembly_fasta) from ch_busco

    output:
    file("${pair_id}_trinity_assembly_stats.txt") into ch_trinity_stats_out

    script:
    """
    /usr/local/bin/trinityrnaseq/util/TrinityStats.pl $trin_assembly_fasta >& ${pair_id}_trinity_assembly_stats.txt
    """

}

// // busco
// process busco{
//     cache 'lenient'
//     tag "$pair_id"
//     container '11cfd93a8245'
//     containerOptions '-u $(id -u):$(id -g) -v $(pwd):/busco_wd'
//     publishDir busco_stats_publish_dir, mode: 'copy', saveAs: {filename -> filename.replaceAll('.txt', "_${pair_id}_.txt")}
//     stageInMode 'copy'
//     errorStrategy 'ignore'
//     cpus params.busco_threads

//     input:
//     tuple val(pair_id), file(trin_assembly_fasta) from ch_trinity_stats

//     output:
//     file("busco_results/*.txt") into ch_busco_out

//     script:
//     """
//     busco -i $trin_assembly_fasta -o busco_results -m transcriptome --auto-lineage-euk -c ${task.cpus} -f
//     mv busco_results/*.txt .
//     """
// }


// quantify reads mapping back to the assembly
// first, make indexed bowtie2 db of the assembly
process assembly_bowtie_db{
    cache 'lenient'
    tag "$pair_id"
    cpus params.bowtie2_threads
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    storeDir assembly_bowtie2_indexed_db_storeDir

    input:
    tuple val(pair_id), file(trin_assembly_fasta) from ch_assembly_bowtie_db

    output:
    tuple val(pair_id), file("${pair_id}_bt_db*") into ch_bowtie2_mapping_stats_indexed_assembly

    script:
    """
    bowtie2-build $trin_assembly_fasta ${pair_id}_bt_db --threads ${task.cpus}
    """

}

// ch_bowtie2_mapping_stats_indexed_assembly.join(ch_bowtie2_mapping_stats_reads).view()

// second, map the reads back to the bowtie2 indexed assembly db
process bowtie2_mapping_stats_reads{
    cache 'lenient'
    tag "$pair_id"
    publishDir read_mapping_stats_publish_dir, mode: 'copy'
    cpus params.bowtie2_threads
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'

    input:
    tuple val(pair_id), file(assembly_db), file(reads) from ch_bowtie2_mapping_stats_indexed_assembly.join(ch_bowtie2_mapping_stats_reads)

    output:
    file("${pair_id}_align_stats.txt")

    script:
    """
    bowtie2 -p ${task.cpus} -q --no-unal -k 20 -x ${pair_id}_bt_db -1 ${reads[0]} \\
    -2 ${reads[1]} 2>${pair_id}_align_stats.txt
    """
}

// Finally let's get the Ex90N50 and Ex90 Gene Count
process rsem{
    cache 'lenient'
    tag "$pair_id"
    cpus params.rsem_threads
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'

    input:
    tuple val(pair_id), file(assembly), file(reads) from ch_rsem_assemblies.join(ch_rsem_reads)

    output:
    tuple file("${pair_id}.RSEM.isoforms.results"), file("${pair_id}.RSEM.genes.results")

    script:
    """
    /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts $assembly \\
    --left ${reads[0]} --right ${reads[1]} \\
    --seqType fq --est_method RSEM --output_dir out_dir --aln_method bowtie2 --SS_lib_type RF \\
    --thread_count ${task.cpus} --trinity_mode --prep_reference
    mv out_dir/RSEM.isoforms.results ${pair_id}.RSEM.isoforms.results
    mv out_dir/RSEM.genes.results ${pair_id}.RSEM.genes.results
    """
}
// 
// TODO as an initial set of quality metrics
// take a look at the below suggestion from 
// https://academic.oup.com/gigascience/article/8/5/giz039/5488105
/*
On the basis of our observations, we suggest initially using reference-free metrics as provided 
by the TransRate [42] software. In general, TransRate’s optimal assembly score seems to be a 
good measure of the quality of an assembly. Assemblies that needed fewer contigs for a comprehensive 
description of the whole transcriptome also achieved in most cases good TransRate scores (Table S6). 
However, this score can be calculated only for paired-end RNA-Seq data at the moment.

If biological/reference-based metrics should be included, the 95%-assembled isoforms statistics calculated 
by rnaQUAST [39], as well as the scores calculated by BUSCO [43,42] and the number of fully reconstructed 
protein-coding transcripts, are good metrics for the evaluation of the best assembly results.


*/

