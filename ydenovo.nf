#!/usr/bin/env nextflow

// TODO consider working with two sets of vals that represent the assembly id and the sample id.

bin_dir = "${workflow.launchDir}/bin"
launch_dir = "${workflow.launchDir}"
tru_seq_pe_fasta_path = "${workflow.launchDir}/TruSeq3-PE.fa"


// Publish dirs
// TODO consider making the differentiatl publishDir via the saveAs option of publishDir
// TODO put all of the subsampled outputs into a separate folder to prevent having to clean up
// the main output folder prior to delivery
if (params.subsample){
    params.output_dir = "${workflow.launchDir}/outputs_sub_sampled/"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim_sub_${params.subsample_depth}"].join(File.separator)
    fastqc_post_correct_publish_dir = [params.output_dir, "fastqc_post_correct_sub_${params.subsample_depth}"].join(File.separator)
    mmseqs_read_query_db_store_dir = [params.output_dir, "mmseqs_read_query_dbs_sub_${params.subsample_depth}"].join(File.separator)
    mmseqs_taxonomy_publish_dir = [params.output_dir, "mmseqs_taxonomy_sub_${params.subsample_depth}_${params.tax_evalue}"].join(File.separator)
    cleaned_reads_publish_dir = [params.output_dir, "cleaned_reads_sub_${params.subsample_depth}"].join(File.separator)
    trinity_assembly_publish_dir = [params.output_dir, "trinity_assembly_sub_${params.subsample_depth}"].join(File.separator)
    trinity_stats_publish_dir = [params.output_dir, "trinity_basic_stats_sub_${params.subsample_depth}"].join(File.separator)
    busco_stats_publish_dir = [params.output_dir, "busco_stats_sub_${params.subsample_depth}"].join(File.separator)
    assembly_bowtie2_indexed_db_storeDir = [params.output_dir, "assembly_bowtie2_indexed_dbs_sub_${params.subsample_depth}"].join(File.separator)
    read_mapping_stats_publish_dir = [params.output_dir, "read_mapping_stats_sub_${params.subsample_depth}"].join(File.separator)
    expression_quantification_rsem_publish_dir = [params.output_dir, "expression_quantification_rsem_sub_${params.subsample_depth}"].join(File.separator)
}else{
    params.output_dir = "${workflow.launchDir}/outputs"
    fastqc_pre_trim_publish_dir = [params.output_dir, "fastqc_pre_trim"].join(File.separator)
    fastqc_post_trim_publish_dir = [params.output_dir, "fastqc_post_trim"].join(File.separator)
    fastqc_post_correct_publish_dir = [params.output_dir, "fastqc_post_correct"].join(File.separator)
    mmseqs_read_query_db_store_dir = [params.output_dir, "mmseqs_read_query_dbs"].join(File.separator)
    mmseqs_taxonomy_publish_dir = [params.output_dir, "mmseqs_taxonomy_${params.tax_evalue}"].join(File.separator)
    cleaned_reads_publish_dir = [params.output_dir, "cleaned_reads"].join(File.separator)
    trinity_assembly_publish_dir = [params.output_dir, "trinity_assembly"].join(File.separator)
    trinity_stats_publish_dir = [params.output_dir, "trinity_basic_stats"].join(File.separator)
    busco_stats_publish_dir = [params.output_dir, "busco_stats"].join(File.separator)
    assembly_bowtie2_indexed_db_storeDir = [params.output_dir, "assembly_bowtie2_indexed_dbs"].join(File.separator)
    read_mapping_stats_publish_dir = [params.output_dir, "read_mapping_stats"].join(File.separator)
    expression_quantification_rsem_publish_dir = [params.output_dir, "expression_quantification_rsem"].join(File.separator)
}

// DB paths
// For the docker mounting we need the directory that contains the SILVA database
File nt_file = new File(params.mmseqs_nt_path);
nt_db_dir = nt_file.getParent();
File b2_file = new File(params.path_to_bowtie2_silva_db);
bowtie2_silva_db_dir = b2_file.getParent();

// Threads
params.trimmomatic_threads = 50
params.rcorrector_threads = 50
params.trinity_threads = 50
params.mmseqs_threads = 20
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
        container 'biocontainers/seqtk:v1.3-1-deb_cv1'
        containerOptions '-u $(id -u):$(id -g)'
        publishDir "${params.output_dir}/${params.subsample_depth}_subsampled_reads", mode: 'copy', overwrite: true
        
        input:
        tuple val(pair_id), file(reads) from ch_subsample

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
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir fastqc_pre_trim_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_pre_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${pair_id}.pre_trim.fastqc.html" into ch_fastqc_pre_trim_output

    script:
    """
    fastqc -o . $fastq_file
    mv *.html ${pair_id}.pre_trim.fastqc.html
    """
}


/* 
    Trim reads with trimmomatic
*/
process trimmomatic{
    cache 'lenient'
	tag "${pair_id}"
    container 'davelabhub/trimmomatic:0.39--1'
    containerOptions '-u $(id -u):$(id -g)'
	
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
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir fastqc_post_trim_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_post_trim.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${pair_id}.post_trim.fastqc.html" into ch_fastqc_post_trim_output

    script:
    """
    fastqc -o . $fastq_file
    mv *.html ${pair_id}.post_trim.fastqc.html
    """
}


// Now error correction
process rcorrector{
    cache 'lenient'
    tag "${trimmed_read_one}"
    container 'tabotaab/rcorrector:latest'
    containerOptions '-u $(id -u):$(id -g)'
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
// ch_bowtie2_silva_mapping.view()
process bowtie2_silva_mapping{
    tag "${pair_id}"
    container 'biocontainers/bowtie2:v2.4.1_cv1'
    containerOptions "-u \$(id -u):\$(id -g) -v ${bowtie2_silva_db_dir}:${bowtie2_silva_db_dir}"
    cpus params.bowtie2_threads
    publishDir cleaned_reads_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), file(fastqs) from ch_bowtie2_silva_mapping

    output:
    tuple val(pair_id), file("${pair_id}_paired_unaligned_{R1,R2}.fq.gz") into ch_trinity,ch_fastqc_post_correct,ch_sub_for_tax,ch_bowtie2_mapping_stats_reads,ch_rsem_reads

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
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir fastqc_post_correct_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), file(fastq_file) from ch_fastqc_post_correct.flatMap{[["${it[0]}_1", it[1][0]], ["${it[0]}_2", it[1][1]]]}

    output:
    file "${pair_id}.post_correct.fastqc.html" into ch_fastqc_post_correct_output

    script:
    """
    fastqc -o . $fastq_file
    mv *.html ${pair_id}.post_correct.fastqc.html
    """
}

// If not already subsampled or if the subsample is greater than 10000
// subsample before running the taxonomy
if (!params.subsample || params.subsample_depth > 10000){
    process sub_for_tax{
        tag "${pair_id}"
        container 'biocontainers/seqtk:v1.3-1-deb_cv1'
        containerOptions '-u $(id -u):$(id -g)'
        
        input:
        tuple val(pair_id), file(reads) from ch_sub_for_tax

        output:
        tuple val("${pair_id}"), file("${pair_id}_sub_10000*.fq.gz") into ch_mmseq_create_query_dbs

        script:
        read_out_one = reads[0].getName().replaceAll("${pair_id}", "${pair_id}_sub_10000")
        read_out_two = reads[1].getName().replaceAll("${pair_id}", "${pair_id}_sub_10000")
        
        """
        seqtk sample -s100 ${reads[0]} 10000 | gzip > ${read_out_one}
        seqtk sample -s100 ${reads[1]} 10000 | gzip > ${read_out_two}
        """
    }
}else{
    // No script. This merely redirects the channel into a new channel with no processing
    process no_sub_for_tax{
        cache 'lenient'
        tag "$fastq_file"

        input:
        tuple val(pair_id), file(fastq_file) from ch_sub_for_tax

        output:
        tuple val(pair_id), file(fastq_file) into ch_mmseq_create_query_dbs

        script:
        """
        echo redirecting
        """
    }
}
// // /* At this point we want to implement some taxonomic checks to see that we are working
// //    with the orgnisms that we expect
// //    An mmseqs taxonomy database has been create as per section "Taxonomy assignment" of:
// //    https://mmseqs.com/latest/userguide.pdf
// //    We will run an mmseqs query against this database for each of the seq files
// //    We will break this up into the following processes
// //    1 - mmseqs creatdb --> queryDB
// //    2 - mmseqs taxonomy --> taxonomyResult
// //         We will limit matches to very confident matches with evalues below 1e-50 
// //         per http://www.metagenomics.wiki/tools/blast/evalue#:~:text=Blast%20hits%20with%20an%20E,good%20hit%20for%20homology%20matches.
// //    3 - mmseqs createtsv --> taxonomyResult.tsv
// //    3 - mmseqs taxonomyreport --> taxonomyResult_report (--report-mode 1 can be used)
// //         We will produce the standard taxonomyResult_report and the html krona report.    
// //    We will want to process on a per readfile basis
// //    see if we can explicitly code in the files
// // */
process mmseq_create_query_dbs{
    cache 'lenient'
    tag "$fastq_file"
    container 'soedinglab/mmseqs2:latest'
    containerOptions '-u $(id -u):$(id -g)'

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

// // NB there is currently a bug in mmseqs where by they are using a FLOAT to store
// // the evalue number. As such the eval will break if you go much above 1e-35.

// // of this on run time.
// // https://github.com/soedinglab/MMseqs2
process mmseq_taxonomy{
    cache 'lenient'
    tag "${db_files[0]}"
    container 'soedinglab/mmseqs2:latest'
    containerOptions "-u \$(id -u):\$(id -g) -v ${nt_db_dir}:${nt_db_dir}"
    cpus params.mmseqs_threads
    publishDir mmseqs_taxonomy_publish_dir, mode: 'copy', overwrite: true
    maxForks 1

    input:
    tuple val(pair_id), file(db_files) from ch_mmseq_taxonomy

    output:
    tuple file("${pair_id}.taxonomyResult.tsv"),\
    file("${pair_id}.taxonomyResult_report"),\
    file("${pair_id}.report.html") into ch_mmseq_taxonomy_out

    script:
    output_tsv = "${pair_id}.taxonomyResult.tsv"
    output_report = "${pair_id}.taxonomyResult_report"
    output_html = "${pair_id}.report.html"
    """
    mmseqs taxonomy  ${pair_id}.queryDB ${params.mmseqs_nt_path} taxonomyResult tmp --threads ${task.cpus} -e ${params.tax_evalue} -s 7.0
    mmseqs createtsv ${pair_id}.queryDB taxonomyResult $output_tsv
    mmseqs taxonomyreport ${params.mmseqs_nt_path} taxonomyResult $output_report
    mmseqs taxonomyreport ${params.mmseqs_nt_path} taxonomyResult $output_html --report-mode 1
    rm -r tmp
    """
}

// // FUTURE do two kmer sizes for trinity. Check what the default is but look at
// // possibly using 25 and 32 based on
// // https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0232005
// // Also, spadesrna seems like a smart choice for a second assembler.
// // Now do the trinity assembly. We want to export the assembly and the gene to transcript map


// // channel manipulation below does:
// // 1 - concatenate the sample_ids to only their base e.g. NG-25417_DK_1_lib406032_6953_2 --> NG-25417_DK
// // 2 - There are then 2 identical keys for each sample type which we group with groupTuple
// // 3 - Finally, we join the two lists each containing the two sequencing files from the original key files pairs
// // https://github.com/trinityrnaseq/trinityrnaseq/wiki
process trinity{
    cache 'lenient'
    tag "$pair_id"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    cpus params.trinity_threads
    publishDir trinity_assembly_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), file(fastqs) from ch_trinity.map{ key, files -> tuple( key[0..10], files ) }.groupTuple().map{key, files -> tuple(key, files[0] + files[1])}

    output:
    tuple val(pair_id), file("${pair_id}.Trinity.fasta") into ch_trinity_stats,ch_busco,ch_assembly_bowtie_db,ch_rsem_assemblies,ch_ExN50_assembly
    tuple val(pair_id), file("${pair_id}.gene_trans_map") into ch_abund_to_matrix_gene_map

    script:
    // NB that the output directory for each trinity assembly must have 'trinity' in it.
    """
    Trinity --left ${fastqs[0]},${fastqs[2]} --right ${fastqs[1]},${fastqs[3]} --seqType fq --max_memory 150G --CPU ${task.cpus} \\
    --min_contig_length 250 --output trinity --full_cleanup  --SS_lib_type RF
    mv trinity.Trinity.fasta ${pair_id}.Trinity.fasta
    mv trinity.Trinity.fasta.gene_trans_map ${pair_id}.gene_trans_map
    """
}

// // first do some basic stats using the Trinity stats package.
process trinity_stats{
    cache 'lenient'
    tag "$pair_id"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir trinity_stats_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_id), file(trin_assembly_fasta) from ch_busco

    output:
    file("${pair_id}_trinity_assembly_stats.txt") into ch_trinity_stats_out

    script:
    """
    /usr/local/bin/trinityrnaseq/util/TrinityStats.pl $trin_assembly_fasta >& ${pair_id}_trinity_assembly_stats.txt
    """
}

// busco
// https://busco.ezlab.org/busco_userguide.html
process busco{
    cache 'lenient'
    tag "$pair_id"
    container 'ezlabgva/busco:v5.0.0_cv1'
    containerOptions '-u $(id -u):$(id -g) -v $(pwd):/busco_wd'
    publishDir busco_stats_publish_dir, mode: 'copy', saveAs: {filename -> filename.replaceAll('.txt', "_${pair_id}_.txt")}, overwrite: true
    errorStrategy 'ignore'
    cpus params.busco_threads

    input:
    tuple val(pair_id), file(trin_assembly_fasta) from ch_trinity_stats

    output:
    tuple file("*eukaryota_odb10*.txt"), file("*alveolata_odb10*.txt"), file("*stramenopiles_odb10*.txt") into ch_busco_out

    script:
    """
    busco -i $trin_assembly_fasta -o busco_results -m transcriptome -l eukaryota_odb10 -c ${task.cpus} -f
    mv busco_results/*eukaryota_odb10*.txt .
    busco -i $trin_assembly_fasta -o busco_results -m transcriptome -l alveolata_odb10 -c ${task.cpus} -f
    mv busco_results/*alveolata_odb10*.txt .
    busco -i $trin_assembly_fasta -o busco_results -m transcriptome -l stramenopiles_odb10 -c ${task.cpus} -f
    mv busco_results/*stramenopiles_odb10*.txt .
    rm -r busco_downloads
    """
}


// quantify reads mapping back to the assembly
// first, make indexed bowtie2 db of the assembly
process assembly_bowtie_db{
    cache 'lenient'
    tag "$pair_id"
    cpus params.bowtie2_threads
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'

    input:
    tuple val(pair_id), file(trin_assembly_fasta) from ch_assembly_bowtie_db

    output:
    tuple val(pair_id), file("${pair_id}_bt_db*") into ch_bowtie2_mapping_stats_indexed_assembly

    script:
    """
    bowtie2-build $trin_assembly_fasta ${pair_id}_bt_db --threads ${task.cpus}
    """

}

// // Here we need to reassociate each of the assemblies back to its original reads
// // second, map the reads back to the bowtie2 indexed assembly db
process bowtie2_mapping_stats_reads{
    cache 'lenient'
    tag "${pair_ids[0]}"
    publishDir read_mapping_stats_publish_dir, mode: 'copy', overwrite: true
    cpus params.bowtie2_threads
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'

    input:
    tuple val(pair_ids), file(reads), file(assembly_db) from ch_bowtie2_mapping_stats_indexed_assembly.mix(ch_bowtie2_mapping_stats_reads).collect().flatMap{ items_list -> println "THIS is starting the pair_read_to_assembly"
        println "THIS is the items_list ${items_list}"
        def return_list = []
        // Map the val of the assembly objects to the val of the seq file objects
        def assembly_val_to_file_val_map = [:]
        for (i = 0; i <items_list.size(); i+=2) {
            println "items_list[i] is ${items_list[i]}"
            for (j = 0; j <items_list.size(); j+=2) {
                println "items_list[j] is ${items_list[j]}"
                try{
                    items_list[j].startsWith(items_list[i])
                }catch(Exception ex){
                    println "ERRORHERE i: ${i}, j: ${j}, items_list[i]:${items_list[i]}, items_list[j]:${items_list[j]}"
                }
                if(items_list[j].startsWith(items_list[i]) && items_list[i] != items_list[j]){
                    // Then items_list[i] is substring of items_list[j]
                    // Then items_list[i] is assembly related to reads items_list[j]
                    if(assembly_val_to_file_val_map.containsKey(items_list[i])){
                        // Then we already have one of the reads mapped to items_list[i]
                        // So add this second one to the list
                        current_read_val_list = assembly_val_to_file_val_map.get(items_list[i])
                        new_read_val_list = current_read_val_list + items_list[j]
                        assembly_val_to_file_val_map[items_list[i]] = new_read_val_list
                    }else{
                        // Then we don't have either of the reads mapped to the assembly val yet
                        assembly_val_to_file_val_map[items_list[i]] = [items_list[j]]
                    }
                }
            }
        }
        
        // Map each val object to its list of files
        def val_to_files_list_map = [:]
        for (i = 0; i <items_list.size(); i+=2) {
            
            val_to_files_list_map.put(items_list[i], items_list[i+1])
        }

        // Finally populate the return_list according to the mappings we have made
        // We want the strucutre to be tuples (2 items) with a tuple for every read set.
        // First item is the val of the read set
        // Second item is a list containing the R1, R2, followed by assembly files in that order
        assembly_val_to_file_val_map.each { assembly_key, read_val_list ->
            read_val_list.each{read_name ->
                def tup = tuple(read_name, val_to_files_list_map[read_name] + val_to_files_list_map[assembly_key])
                return_list << tuple([read_name, assembly_key], val_to_files_list_map[read_name], val_to_files_list_map[assembly_key])
            }
        }
        return return_list
        }

    output:
    file("${pair_ids[0]}_align_stats.txt")

    script:
    """
    bowtie2 -p ${task.cpus} -q --no-unal -k 20 -x ${pair_ids[1]}_bt_db -1 ${reads[0]} \\
    -2 ${reads[1]} 2>${pair_ids[0]}_align_stats.txt
    """
}

// // Finally let's get the Ex90N50 and Ex90 Gene Count
process rsem{
    cache 'lenient'
    tag "$pair_id"
    cpus params.rsem_threads
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir expression_quantification_rsem_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(pair_ids), file(reads), file(assembly) from ch_rsem_assemblies.mix(ch_rsem_reads).collect().flatMap{ items_list -> println "THIS is starting the pair_read_to_assembly"
        println "THIS is the items_list ${items_list}"
        def return_list = []
        // Map the val of the assembly objects to the val of the seq file objects
        def assembly_val_to_file_val_map = [:]
        for (i = 0; i <items_list.size(); i+=2) {
            println "items_list[i] is ${items_list[i]}"
            for (j = 0; j <items_list.size(); j+=2) {
                println "items_list[j] is ${items_list[j]}"
                try{
                    items_list[j].startsWith(items_list[i])
                }catch(Exception ex){
                    println "ERRORHERE i: ${i}, j: ${j}, items_list[i]:${items_list[i]}, items_list[j]:${items_list[j]}"
                }
                if(items_list[j].startsWith(items_list[i]) && items_list[i] != items_list[j]){
                    // Then items_list[i] is substring of items_list[j]
                    // Then items_list[i] is assembly related to reads items_list[j]
                    if(assembly_val_to_file_val_map.containsKey(items_list[i])){
                        // Then we already have one of the reads mapped to items_list[i]
                        // So add this second one to the list
                        current_read_val_list = assembly_val_to_file_val_map.get(items_list[i])
                        new_read_val_list = current_read_val_list + items_list[j]
                        assembly_val_to_file_val_map[items_list[i]] = new_read_val_list
                    }else{
                        // Then we don't have either of the reads mapped to the assembly val yet
                        assembly_val_to_file_val_map[items_list[i]] = [items_list[j]]
                    }
                }
            }
        }
        
        // Map each val object to its list of files
        def val_to_files_list_map = [:]
        for (i = 0; i <items_list.size(); i+=2) {
            
            val_to_files_list_map.put(items_list[i], items_list[i+1])
        }

        // Finally populate the return_list according to the mappings we have made
        // We want the strucutre to be tuples (2 items) with a tuple for every read set.
        // First item is the val of the read set
        // Second item is a list containing the R1, R2, followed by assembly files in that order
        assembly_val_to_file_val_map.each { assembly_key, read_val_list ->
            read_val_list.each{read_name ->
                def tup = tuple(read_name, val_to_files_list_map[read_name] + val_to_files_list_map[assembly_key])
                return_list << tuple([read_name, assembly_key], val_to_files_list_map[read_name], val_to_files_list_map[assembly_key])
            }
        }
        return return_list}

    output:
    tuple val("${pair_ids[1]}"), file("${pair_ids[0]}.RSEM.isoforms.results") into ch_abund_to_matrix_results
    tuple val("${pair_ids[1]}"), file("${pair_ids[0]}.RSEM.genes.results") into ch_abund_to_matrix_gen_results_out

    script:
    """
    /usr/local/bin/trinityrnaseq/util/align_and_estimate_abundance.pl --transcripts ${assembly} \\
    --left ${reads[0]} --right ${reads[1]} \\
    --seqType fq --est_method RSEM --output_dir out_dir --aln_method bowtie2 --SS_lib_type RF \\
    --thread_count ${task.cpus} --trinity_mode --prep_reference
    mv out_dir/RSEM.isoforms.results ${pair_ids[0]}.RSEM.isoforms.results
    mv out_dir/RSEM.genes.results ${pair_ids[0]}.RSEM.genes.results
    """
}

// TODO consider working with two sets of vals that represent the assembly id and the sample id.
process abund_to_matrix{
    cache 'lenient'
    tag "$assembly_id"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir expression_quantification_rsem_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(assembly_id), file(abund_results), file(gene_map) from ch_abund_to_matrix_results.groupTuple().join(ch_abund_to_matrix_gene_map)

    output:
    tuple val(assembly_id), file("${assembly_id}.RSEM.isoform.counts.matrix"), file("${assembly_id}.RSEM.isoform.TPM.not_cross_norm"), file("${assembly_id}.RSEM.isoform.TMM.EXPR.matrix") into ch_abund_to_matrix_out_out
    tuple val(assembly_id), file("${assembly_id}.RSEM.isoform.TMM.EXPR.matrix") into ch_ExN50_matrix

    script:
    """
    /usr/local/bin/trinityrnaseq/util/abundance_estimates_to_matrix.pl --est_method RSEM --gene_trans_map ${gene_map} --out_prefix RSEM ${abund_results[0]} ${abund_results[1]}
    mv RSEM.isoform.counts.matrix ${assembly_id}.RSEM.isoform.counts.matrix
    mv RSEM.isoform.TPM.not_cross_norm ${assembly_id}.RSEM.isoform.TPM.not_cross_norm
    mv RSEM.isoform.TMM.EXPR.matrix ${assembly_id}.RSEM.isoform.TMM.EXPR.matrix
    """
}

// Now compute the ExN50 stats per
// https://github.com/trinityrnaseq/trinityrnaseq/wiki/Transcriptome-Contig-Nx-and-ExN50-stats#contig-ex90n50-statistic-and-ex90-gene-count
process ExN50{
    cache 'lenient'
    tag "$assembly_id"
    container 'trinityrnaseq/trinityrnaseq:latest'
    containerOptions '-u $(id -u):$(id -g)'
    publishDir expression_quantification_rsem_publish_dir, mode: 'copy', overwrite: true

    input:
    tuple val(assembly_id), file(abund_matrix), file(assembly_fasta) from ch_ExN50_matrix.join(ch_ExN50_assembly)

    output:
    tuple file("${assembly_id}.RSEM.isoform.TMM.EXPR.matrix.E-inputs"), file("${assembly_id}.ExN50_plot.pdf"), file("${assembly_id}.ExN50.stats"), file("${assembly_id}.Ex90N50.txt") into ch_ExN50_out
    
    script:
    """
    /usr/local/bin/trinityrnaseq/util/misc/contig_ExN50_statistic.pl $abund_matrix $assembly_fasta | tee ExN50.stats
    /usr/local/bin/trinityrnaseq/util/misc/plot_ExN50_statistic.Rscript  ExN50.stats
    cat ${assembly_id}.RSEM.isoform.TMM.EXPR.matrix.E-inputs |  egrep -v ^\\# | awk '\$1 <= 90' | wc -l > Ex90N50.txt
    mv ExN50_plot.pdf ${assembly_id}.ExN50_plot.pdf
    mv ExN50.stats ${assembly_id}.ExN50.stats
    mv Ex90N50.txt ${assembly_id}.Ex90N50.txt
    """
}
// //
