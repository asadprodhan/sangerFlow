#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.samples_tsv = './samples.tsv'

params.db = ""
db_name = file(params.db).name
db_path = file(params.db).parent

params.evalue='0.1'
params.topHits = '5'
params.cpus='18'

params.publish_dir = './'


log.info """\

    s a n g e r F L O W
    A Sanger sequencing-based bioinformatics pipeline for pest and pathogen diagnosis 

    Asad Prodhan
    =======================================

    """
    .stripIndent()


//Processing R1
//rename the reads with their sample IDs
process renameR1 {

    tag "${read1.simpleName}"

    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sampleId), path(read1)

    output:
    tuple val(sampleId), path ("*.fasta") 

    script:
    """
    bioawk -c fastx '{ print ">"\$name;print (\$seq) }' ${read1} > ${sampleId}.R1.fasta
    """
}

//trimming low quality reads
process trimR1 {

    tag "${x.baseName}"
    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sampleId), path(x)	 				

    output:
    tuple val(sampleId), path ("*.trim.R1.fasta"), emit: R1_read 

    script:
    """
    sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D|N]/-/g' $x > ${x.simpleName}.trim.R1.fasta
    
    """
}

//Processing R2
//rename the reads with their sample IDs
process renameR2 {

    tag "${read2.simpleName}"

    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sampleId), path(read2)

    output:
    tuple val(sampleId), path ("*.fasta") 

    script:
    """
    bioawk -c fastx '{ print ">"\$name;print (\$seq) }' ${read2} > ${sampleId}.R2.fasta

    """
}

//R2 reverse complement
process revComp {

    tag "${x.baseName}"
    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sampleId), path(x) 	 				

    output:
    tuple val(sampleId), path ("*") 

    script:
    """
    bioawk -c fastx '{ print ">"\$name;print revcomp(\$seq) }' $x > ${x.simpleName}.revComp.R2.fasta    
    """
}
//trimming low quality reads
process trimR2 {

    tag "${x.baseName}"
    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sampleId), path(x)	 				

    output:
    tuple val(sampleId), path ("*.trim.R2.fasta"), emit: R2_read 

    script:
    """
    sed '/^[^>]/s/[R|Y|W|S|M|K|H|B|V|D|N]/-/g' $x > ${x.simpleName}.revComp.trim.R2.fasta 
        
    """
}

//concat the reads
process concatination {

    tag "${processed_read1.baseName}"
    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sample), path(processed_read1), path(processed_read2)				

    output:
    tuple val(sample), path("*.Merged.fasta")

    script:
    """
    cat ${processed_read1} ${processed_read2} > ${processed_read1.simpleName}.Merged.fasta 
                
    """
}


// Read 1 and Read2 alignment
process alignment {

    tag "${x.baseName}"
    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sample), path(x) 
    
    output:
    tuple val(sample), path ("*")

    script:
    """
    clustalo --in=$x --out=${x.simpleName}.clustal --force --outfmt=clu --wrap=80 --dealign                
    """
}

// Extracting the most representative sequence from the alignment file
process contig {

    tag "${x.baseName}"
    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sample), path(x) 
    
    output:
    tuple val(sample), path ("${x.simpleName}.fasta")

    script:
    """
    cons -sequence $x -outseq ${x.simpleName}.fasta              
    """
}

// Renaming the contig fasta header with actual sample name

process renameContigFastaHeader {

    tag "${x.simpleName}"

    publishDir (
        path: "${params.publish_dir}/temp",
        mode: 'copy',
        overwrite: 'true',
    )

    input:
    tuple val(sample), path(x) 

    output:
    tuple val(sample), path ("*.fasta"), emit: contig_seq

    script:
    """
    awk '{gsub(/EMBOSS_001/,"${x.simpleName}"); print}' $x > ${x.simpleName}.withHeader.fasta
    
    """
}

//Blastn
process blastn {
    tag "${query.baseName}"
	publishDir "${params.publish_dir}/results", mode:'copy'

    input:
	tuple val(query_id), path(query)
    path db

	output:
    path ("${query.baseName}_blast_sort.tsv"), emit: blastHits
    path ("${query.baseName}_blast.xml")
    path ("${query.baseName}_blast.html")

    script:
    """
    blastn \
        -query ${query} -db "${db}/${db_name}/nt" \
        -outfmt 11 -out ${query.baseName}_blast.asn \
        -evalue ${params.evalue} \
        -num_threads ${params.cpus}

	blast_formatter \
        -archive ${query.baseName}_blast.asn \
        -outfmt 5 -out ${query.baseName}_blast.xml
    
    blast_formatter \
        -archive ${query.baseName}_blast.asn \
        -html -out ${query.baseName}_blast.html

	blast_formatter \
        -archive ${query.baseName}_blast.asn \
        -outfmt "6 qaccver saccver pident length evalue bitscore mismatch gapopen qstart qend sstart send stitle" -out ${query.baseName}_blast_unsort.tsv

    sort -n -r -k 6 ${query.baseName}_blast_unsort.tsv > ${query.baseName}_blast_sort.tsv  

    """
}

//Blastn_topHits
process blastnTopHits {
    tag "${query.baseName}"
	publishDir "${params.publish_dir}/temp", mode:'copy'

    input:
	path(query)

	output:
    path ("${query.baseName}_blast_sort_topHits.tsv"), emit: topHits

    script:
    """
    head -n ${params.topHits} ${query} > ${query.baseName}_blast_sort_topHits.tsv 

    """
}

//concatinating sorted blastn files into a single file and adding headers
process blastHits_concat {

	publishDir "${params.publish_dir}/results", mode:'copy'

    input:
    file combinedFiles
    
    output:
    path "concatenatedHits_withHeaders.tsv" 
	
    script:
    
    """
    cat ${combinedFiles} >> concatenatedHits.tsv
    awk 'BEGIN{print "qaccver\tsaccver\tpident\tlength\tevalue\tbitscore\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tstitle"}1' concatenatedHits.tsv > concatenatedHits_withHeaders.tsv
	
    """
}

// //workflow definition
workflow {

    Channel.fromPath( params.samples_tsv )
        | splitCsv( header: true, sep: '\t' )
        | map { row -> tuple( row.sampleId, file(row.read1) ) }
        | renameR1
        | trimR1
        
    Channel.fromPath( params.samples_tsv )
        | splitCsv( header: true, sep: '\t' )
        | map { row -> tuple( row.sampleId, file(row.read2) ) }
        | renameR2
        | revComp
        | trimR2
    
    trimR1.out 
        | combine( trimR2.out, by: 0 ) 
        | map { sample, processed_read1, read2 ->
            processed_read2 = read2 instanceof Path ? [read2] : read2

            tuple( sample, processed_read1, processed_read2 )
        } 
        | transpose( by: 2 ) 
        | concatination 
        | alignment
        | contig
        | renameContigFastaHeader
	
    blastn (renameContigFastaHeader.out.contig_seq, db_path)
    blastnTopHits(blastn.out.blastHits)
    
    blastHits_concat(blastnTopHits.out.topHits.collect())
}
