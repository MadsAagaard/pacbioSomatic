#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
date2=new Date().format( 'yyMMdd HH:mm:ss' )
user="$USER"
runID="${date}.${user}"


//////////////////////////// SWITCHES ///////////////////////////////// 


// OUTPUT locations:

//lrsStorage="/lnx01_data3/storage/pacBioLRS/analyzedData/${params.assembly}/SOMATIC/"



/* --------------------------------- NAMING --------------------------------- 

single sample (preprocess, deepvariant, sawfish, hiphase, QC):

${meta.npn}.${meta.sampletype}.${params.genome_version}.{toolname}.{filetype}

metadata: [id, npnNormal, npnTumor, pcgr, type, npn]




T-N analysis (all somatic analysis requiring both T and N sample):
${meta.prefixTN}.{toolname}.{filetype}

Metadata: [id, npn, npnNormal, npnTumor, pcgr]
metamap can (should!) be used for joining


*/




///////////////////////////////////////////////////////////////////
////////////// --- Preprocessing --- /////////////////////////
///////////////////////////////////////////////////////////////////


process create_fofn {
    publishDir "${meta.id}/documents/", mode: 'copy',pattern: '*.fofn'
    input:
    tuple val(meta), path(data) //ubam

    output:
    tuple val(meta), path("${meta.npn}.${meta.sampletype}.fofn")
    script:
    """
    `realpath ${data} > ${meta.npn}.${meta.sampletype}.fofn`
    """
} 

process inputFiles_symlinks_ubam{
    errorStrategy 'ignore'
    publishDir "${meta.id}/inputSymlinks/", mode: 'link', pattern: '*.{bam,pbi}'

    input:

    tuple val(meta), path(data)   
    //data (default): 0:ubam, 1:ubam pbi
    output:
    tuple val(meta), path(data)
    script:
    """
    """
}


process pbmm2_align_mergedData {
    label "veryHigh"
    tag "$meta.npn"
    conda "${params.pbmm2}"


    publishDir "${meta.id}/documents/", mode: 'copy',pattern: '*.fofn'


    input:
    tuple val(meta), path(fofn)
    
    output:
    tuple val(meta), path("${meta.prefix}.pbmm2.bam"), path("${meta.prefix}.pbmm2*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.npn} \
    ${params.genome_mmi} \
    ${fofn} \
    ${meta.prefix}.pbmm2.bam
    """
}

/*
    samtools view \
    -T ${params.genome_fasta} \
    -C \
    -o ${meta.prefix}.pbmm2.cram ${meta.prefix}.pbmm2.bam

    samtools index ${meta.prefix}.pbmm2.cram
    */

///////////////////////////////////////////////////////////////////
////////////// --- Per sample processes --- /////////////////////////
///////////////////////////////////////////////////////////////////
process deepvariant{
    tag "$meta.id"
    label "veryHigh"

    input:
    tuple val(meta), val(data) // data from joined ch: [bamN:bamNormal,baiN:baiNormal,bamT:bamTumor,baiT:baiTumor] access as e.g. "data.bamN" 

    output:
    tuple val(meta), path("${meta.prefixNormal}.deepVariant.vcf.gz"), path("${meta.prefixNormal}.deepVariant.vcf.gz.tbi"), emit: dv_vcf
    //tuple val(meta), path("${meta.prefixTumor}.deepVariant.FAKE.vcf.gz"), path("${meta.prefixTumor}.deepVariant.FAKE.vcf.gz.tbi"), emit: dvTumorFAKE_vcf

    tuple val(meta), path("${meta.prefixNormal}.deepVariant.g.vcf.gz"), emit: dv_gvcf    
    //path("${meta.npn}.deepvariant.vcf_stats_report.txt")
    """
    singularity run -B ${params.s_bind} ${params.simgpath}/${params.deepvariant_image} /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref=${params.genome_fasta} \
    --reads=${data.bamNormal}  \
    --output_vcf=${meta.prefixNormal}.deepVariant.vcf.gz \
    --output_gvcf=${meta.prefixNormal}.deepVariant.g.vcf.gz \
    --num_shards=${task.cpus}

    """    
}


process sawFish2{
    tag "$meta.npnNormal"
    label "high"
    conda "${params.sawfish2}"

    //publishDir "${params.lrsStorageBase}/sawfish/", mode: 'copy', pattern:"*.sawfishSV.vcf.*"

    publishDir "${meta.id}/toolsOutputDNA/sawFish/sawfish_supporting_data/", mode: 'copy', pattern: "*.sawfishSV.*"

    input:
    tuple val(meta), val(data)
    
    output:
    //tuple val(meta), path("*.sawfishSV.*"),emit: all_sawfish_output

    tuple val(meta), path("${meta.prefixNormal}.sawfishSV.vcf.gz"), path("${meta.prefixNormal}.sawfishSV.vcf.gz.tbi"), emit:sv_vcf

    tuple val(meta), path("${meta.prefixNormal}.sawfishSV.supporting_reads.json.gz"), emit: sv_supporting_reads

    script:
    """
    sawfish discover \
    --threads ${task.cpus} \
    --ref ${params.genome_fasta} \
    --bam ${data.bamNormal} \
    --cnv-excluded-regions ${params.cnv_exclude_sawfish} \
    --output-dir ${meta.npnNormal}.normal.sawfishDiscover 

    sawfish joint-call \
    --threads ${task.cpus} \
    --report-supporting-reads \
    --sample ${meta.npnNormal}.normal.sawfishDiscover \
    --output-dir ${meta.npnNormal}.normal.sawfishSV 
    
    mv ${meta.npnNormal}.normal.sawfishSV/genotyped.sv.vcf.gz ${meta.prefixNormal}.sawfishSV.vcf.gz

    mv ${meta.npnNormal}.normal.sawfishSV/genotyped.sv.vcf.gz.tbi ${meta.prefixNormal}.sawfishSV.vcf.gz.tbi

   mv ${meta.npnNormal}.normal.sawfishSV/supporting_reads.json.gz ${meta.prefixNormal}.sawfishSV.supporting_reads.json.gz

    mv ${meta.npnNormal}.normal.sawfishSV/samples/*/gc_bias_corrected_depth.bw ${meta.prefixNormal}.sawfishSV.gc_bias_corrected_depth.bw

    mv ${meta.npnNormal}.normal.sawfishSV/samples/*/depth.bw ${meta.prefixNormal}.sawfishSV.depth.bw

    mv ${meta.npnNormal}.normal.sawfishSV/samples/*/copynum.bedgraph ${meta.prefixNormal}.sawfishSV.copynum.bedgraph

    mv ${meta.npnNormal}.normal.sawfishSV/samples/*/copynum.summary.json ${meta.prefixNormal}.sawfishSV.copynum.summary.json

    """
}

process svdb_SawFish {
    tag "$meta.id"
    label "low"
    conda "${params.svdb}"

    publishDir "${params.lrsStorageBase}/sawfish/", mode: 'copy',pattern: "*.sawfishSV.hiphase.svdb.vcf*"

    publishDir {"${meta.id}/toolsOutputDNA/sawFish/"}, mode: 'copy', pattern: "*.sawfishSV.hiphase.svdb.*"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: '*.sawfishSV.hiphase.svdb.vcf*'

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.sawfishSV.hiphase.svdb.*")
    tuple val(meta), path("${meta.prefixNormal}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz"),path("${meta.prefixNormal}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz.tbi"), emit: sawfishAF10
    script:
    """
    svdb --query \
    --query_vcf ${data.sawfish_vcf} \
    --sqdb ${params.sawfish_sqdb} > ${meta.prefixNormal}.sawfishSV.hiphase.svdb.vcf
    
    bgzip ${meta.prefixNormal}.sawfishSV.hiphase.svdb.vcf
    
    bcftools index -t ${meta.prefixNormal}.sawfishSV.hiphase.svdb.vcf.gz

    bcftools view -e 'INFO/FRQ>0.1' ${meta.prefixNormal}.sawfishSV.hiphase.svdb.vcf.gz -Oz -o ${meta.prefixNormal}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz

    bcftools index -t ${meta.prefixNormal}.sawfishSV.hiphase.svdb.AF_below10pct.vcf.gz

    """
}

process hiPhase {
    tag "$meta.id"
    label "intermediate"
    conda "${params.hiphase}"

    publishDir "${meta.id}/alignments/", mode: 'copy', pattern: "*.hiphase.ba*"
    publishDir "${meta.id}/toolsOutputDNA/deepVariant/", mode: 'copy', pattern: "*.hiphase.deepvariant.*"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: '*.hiphase.deepvariant.*'
    //publishDir "${meta.id}/toolsOutputDNA/sawfish_supporting_data/", mode: 'copy', pattern: "*.hiphase.sawfishSV.*"

    input:
    tuple val(meta), val(data), path(vcf), path(sv)
    
    output:
    tuple val(meta), path("${meta.prefixNormal}.hiphase.bam"), path("${meta.prefixNormal}.hiphase.bam.bai"), emit: hiphase_bam_normal                                                

    tuple val(meta), path("${meta.prefixTumor}.hiphase.bam"), path("${meta.prefixTumor}.hiphase.bam.bai"), emit: hiphase_bam_tumor       

    tuple val(meta), path("${meta.prefixNormal}.hiphase.deepvariant.vcf.gz"), path("${meta.prefixNormal}.hiphase.deepvariant.vcf.gz.tbi"), emit: hiphase_dv_vcf

    tuple val(meta), path("${meta.prefixNormal}.hiphase.deepvariant.WES_ROI.vcf.gz"), path("${meta.prefixNormal}.hiphase.deepvariant.WES_ROI.vcf.gz.tbi"), emit: hiphase_dv_roi_vcf

    tuple val(meta), path("${meta.prefixNormal}.hiphase.sawfishSV.vcf.gz"), path("${meta.prefixNormal}.hiphase.sawfishSV.vcf.gz.tbi"), emit: hiphase_sv_vcf
  

    script:
    """
    hiphase \
    --vcf ${vcf[0]} \
    --output-vcf ${meta.prefixNormal}.hiphase.deepvariant.vcf.gz \
    --vcf ${sv[0]} \
    --output-vcf ${meta.prefixNormal}.hiphase.sawfishSV.vcf.gz \
    --bam ${data.bamNormal} \
    --output-bam ${meta.prefixNormal}.hiphase.bam \
    --bam ${data.bamTumor} \
    --output-bam ${meta.prefixTumor}.hiphase.bam \
    --reference ${params.genome_fasta} \
    --threads ${task.cpus} \
    --ignore-read-groups \
    --io-threads ${task.cpus}

    bcftools index -t -f ${meta.prefixNormal}.hiphase.deepvariant.vcf.gz 

    ${params.gatk_exec} SelectVariants \
    -R ${params.genome_fasta} \
    -V  ${meta.prefixNormal}.hiphase.deepvariant.vcf.gz \
    -L ${params.ROI} \
    -O  ${meta.prefixNormal}.hiphase.deepvariant.WES_ROI.vcf.gz

    """
}

process pbCPGtools{
    tag "$meta.prefix"
    label "medium"
    conda "${params.pbCPGtools}"

    publishDir "${meta.id}/toolsOutputDNA/methylation/pbCPGtools/", mode: 'copy', pattern: "*.methylation.{hap1,hap2,combined}.*"
    //publishDir "${params.lrsStorageBase}/methylation/2025/${meta.id}/", mode: 'copy', pattern:"*.bed.*"


    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("${meta.prefix}.hiphase.methylation*")
    
    script:
    """
    aligned_bam_to_cpg_scores \
    --bam ${data.bam} \
    --output-prefix ${meta.prefix}.hiphase.methylation
    """
}

process methBat{
    tag "$meta.prefix"
    label "medium"
    conda "${params.methbat}"

    publishDir "${meta.id}/toolsOutputDNA/methylation/methBatOLD/", mode: 'copy'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.prefix}.met.*")
    tuple val(meta), path("${meta.prefix}.met.CelltypeEstimate.json"), emit: for_yaml_summary
    script:
    """
    methbat segment \
    --input-prefix ${meta.prefix}.hiphase.methylation \
    --output-prefix ${meta.prefix}.methBatSegments

    methbat profile \
    --input-prefix ${meta.prefix}.hiphase.methylation \
    --input-regions ${params.methylationBackground} \
    --output-region-profile ${meta.prefix}.met.profile

    methbat deconvolve \
    --input-prefix ${meta.prefix}.hiphase.methylation \
    --atlas-regions ${params.methbatAtlas} \
    --output-estimate ${meta.prefix}.met.CelltypeEstimate.json  
    
    
    """

}

//// WORKINPROGRESS


process methBatNEW_pileup{
    tag "$meta.id"
    label "intermediateCPU"
    conda "${params.methbat_v1}"

    publishDir {"${meta.id}/toolsOutputDNA/methylation/5mC_pileup/"},   mode: 'copy',   pattern: "*.5mC.bed.*"
    publishDir {"${meta.id}/toolsOutputDNA/methylation/5mC_bedgraphs/"},   mode: 'copy',   pattern: "*.5mC.bedgraph.*"

    publishDir "${params.lrsStorageBase}/methylationNEW/5mC_pileup/",   mode: 'copy',   pattern: "*.5mC.bed.*"

  //  publishDir {"${params.outBase(meta)}/specialAnalysis/methylation/5hmC/"},  mode: 'copy',   pattern: "*.5hmC.bed.*"
   // publishDir {"${params.outBase(meta)}/specialAnalysis/methylation/6mA/"},   mode: 'copy',   pattern: "*.6mA.bed.*"

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.met.*"), path("*.5mC.bedgraph.*")
    tuple val(meta), path("*.5mC.bed.gz"),  path("*.5mC.bed.gz.tbi"),   emit: met5mC
    tuple val(meta), path("*.5hmC.bed.gz"), path("*.5hmC.bed.gz.tbi"),  emit: met5hmC
    tuple val(meta), path("*.6mA.bed.gz"),  path("*.6mA.bed.gz.tbi"),   emit: met6mA
   
    script:
    """
    methbat pileup \
    --threads ${task.cpus} \
    --input-bam ${data.bam} \
    --output-prefix ${meta.prefix}.met.pileup

    zgrep "Total" ${meta.prefix}.met.pileup.5mC.bed.gz | \
    cut -f 1-3,7 | \
    bgzip > ${meta.prefix}.5mC.bedgraph.gz

    tabix -p bed ${meta.prefix}.5mC.bedgraph.gz

    """
}

process methBatNEW_profile_single {
    tag "$meta.prefix"
    label "low"
    conda "${params.methbat_v1}"

    publishDir "${meta.id}/toolsOutputDNA/methylation/5mC_profile/",   mode: 'copy',   pattern: "*.5mC.cpgIslands.profile.tsv"
    publishDir "${meta.id}/toolsOutputDNA/methylation/5mC_CGI_profiles/", mode: 'copy', pattern:"*.profile.tsv"
    publishDir "${params.lrsStorageBase}/methylationNEW/5mC_CGI_profiles/", mode: 'copy', pattern:"*.profile.tsv"
    
    input:
    tuple val(meta), path(data), path(tbi)
    
    output:
    tuple val(meta), path("*.5mC.cpgIslands.profile.tsv")

    script:
    """
    methbat profile \
    --input-regions ${params.methylationCpG_regions} \
    --input-pileup ${data} \
    --output-region-profile ${meta.prefix}.met.5mC.cpgIslands.profile.tsv
    
    methbat deconvolve \
    --input-pileup ${data} \
    --atlas-regions ${params.methbatAtlas} \
    --output-estimate ${meta.prefix}.met.CelltypeEstimate.json  
        
    
    
    """
}
















//// WORKINPROGRESS

///////////////////////////////////////////////////////////////////
/////////////////// --- QC processes --- //////////////////////////
///////////////////////////////////////////////////////////////////


process mosdepthROI {
    tag "$meta.prefix"
    label "low"
    conda "${params.mosdepth}"

    publishDir "${meta.id}/QC/mosdepth/", mode: 'copy'

    input: 
    tuple val(meta), val(data)  

    output:
    tuple val(meta), path("${meta.prefix}_roi.*"),emit: mosdepth_roi
    tuple val(meta), path("*.region.dist.txt"), emit:multiqc
    script:
    def callable=params.genome=="hg38" ? "--by ${params.CALLABLE_ROI}" : "--by 1000"
    """
    mosdepth \
    -t ${task.cpus} \
    $callable \
    ${meta.prefix}_roi \
    ${data.bam}
    """
}

process cramino {
    tag "$meta.prefix"
    label "low"
    conda "${params.cramino}"

    publishDir {"${meta.id}/QC/cramino/"}, mode: 'copy'

    input: 
    tuple val(meta), val(data)  // meta: [npn,datatype,sampletype,id], data: [cram,crai]

    output:
    tuple val(meta), path("${meta.prefix}.craminoQC.txt"), emit: for_yaml_summary

    script:
    """
    cramino \
    -t ${task.cpus} \
    --karyotype \
    --phased \
    ${data.bam} > ${meta.prefix}.craminoQC.txt
    """
}




// TO DO MULTIQC
process multiQC {
    
    label "low"
    publishDir "${meta.id}/QC/", mode: 'copy'


    conda "${params.multiqc}"

    input:
    //path(inputfiles)
    tuple val(meta),  path(data)  
    //   path("_fastqc.*").collect().ifEmpty([])
    // path("${meta.prefixTN}.samtools.sample.stats.txt").collect().ifEmpty([])
    // path("bamQC/*").collect().ifEmpty([]) 
    //path("${meta.prefixTN}.picardWGSmetrics.txt").collect().ifEmpty([]) 

    output:
    path ("*MultiQC*.html")

    script:
    """
    multiqc \
    -c ${params.multiqc_config} \
    -f -q ${launchDir}/${meta.id}/ \
    -n ${meta.prefixTN}.MultiQC.DNA.html
    """
}





///////////////////////////////////////////////////////////////////
/////////////////// --- Paired T-N analyses --- ///////////////////
///////////////////////////////////////////////////////////////////




process deepSomatic {
    label "high"
    tag "$meta.id"
    publishDir "${meta.id}/toolsOutputDNA/deepSomatic/", mode: 'copy'


    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.prefixTN}.deepSomatic.vcf.gz"), path("${meta.prefixTN}.deepSomatic.vcf.gz.tbi"), emit: vcf_idx

    tuple val(meta), path("${meta.prefixTN}.deepSomatic.vcf.gz"), emit: vcf

    tuple val(meta), path("${meta.prefixTN}.deepSomatic.PASS.vcf.gz"),path("${meta.prefixTN}.deepSomatic.PASS.vcf.gz.tbi"), emit: pcgr_vcf
    
    script:
    """
    singularity run -B ${params.s_bind} ${params.simgpath}/${params.deepsomatic_image} run_deepsomatic \
    --model_type=PACBIO \
    --ref=${params.genome_fasta} \
    --reads_normal=${data.bamNormal} \
    --reads_tumor=${data.bamTumor} \
    --sample_name_normal=${meta.npnNormal} \
    --sample_name_tumor=${meta.npnTumor} \
    --output_vcf=${meta.prefixTN}.deepSomatic.vcf.gz \
    --output_gvcf=${meta.prefixTN}.deepSomatic.g.vcf.gz \
    --num_shards=${task.cpus} \
    --regions ${params.ROI} \
    --logging_dir .

    bcftools view -f PASS ${meta.prefixTN}.deepSomatic.vcf.gz -Oz -o ${meta.prefixTN}.deepSomatic.PASS.vcf.gz
    bcftools index -t ${meta.prefixTN}.deepSomatic.PASS.vcf.gz
    """
}

process deepSomatic_edits {
    label "low"
    publishDir "${meta.id}/toolsOutputDNA/deepSomatic/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: '*.normalAdded.*'
    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("${meta.prefixTN}.deepSomatic.PASS.normalAdded.vcf.gz"), emit: vcf
    tuple val(meta), path("${meta.prefixTN}.deepSomatic.PASS.normalAdded.vcf.gz"), path("${meta.prefixTN}.deepSomatic.PASS.normalAdded.vcf.gz.tbi"), emit:vcf_idx
    
    script:
    """
    bcftools view -h -f PASS "${data}" > header.tmp

    awk 'BEGIN{OFS="\\t"}
        /^#CHROM/ {print \$0 "\\t${meta.npnNormal}"; next}
        /^[^#]/ {next}
        {print}
    ' header.tmp > header.modified.tmp


    bcftools view -H -f PASS "${data}" | \
    awk 'BEGIN{OFS="\\t"}
        {
            split(\$9, fmt_fields, ":")
            fake_normal=""
            for(i=1; i<=length(fmt_fields); i++) {
                if(fmt_fields[i]=="GT") fake_normal = fake_normal "0/0"
                else if(fmt_fields[i]=="GQ") fake_normal = fake_normal "."
                else if(fmt_fields[i]=="DP") fake_normal = fake_normal "0"
                else if(fmt_fields[i]=="AD") fake_normal = fake_normal "0,0"
                else if(fmt_fields[i]=="VAF") fake_normal = fake_normal "."
                else if(fmt_fields[i]=="PL") fake_normal = fake_normal "."
                else fake_normal = fake_normal "."
                if(i < length(fmt_fields)) fake_normal = fake_normal ":"
            }
            print \$0 "\\t" fake_normal
        }
    ' > variants.modified.tmp
    cat header.modified.tmp variants.modified.tmp | bgzip -c > ${meta.prefixTN}.deepSomatic.PASS.normalAdded.vcf.gz
    tabix -p vcf ${meta.prefixTN}.deepSomatic.PASS.normalAdded.vcf.gz
    """
}


/* 

IF EOF BGZF ERROR FROM HIPHASE:
samtools reheader -P <(samtools view -H file.bam) file.bam > fixed.bam

*/

process severus {
    label "high"
    tag "$meta.id"
    conda "${params.severus}"

    publishDir "${meta.id}/toolsOutputDNA/severus_somaticSV/", mode: 'copy'

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.prefixTN}.severus/"), emit: severusDir
    tuple val(meta),path("${meta.prefixTN}.severusSomaticSV.vcf"),emit:vcf

    script:
    """
    severus \
    --target-bam ${data.bamTumor} \
    --target-sample ${meta.npnTumor} \
    --control-bam ${data.bamNormal} \
    --control-sample ${meta.npnNormal} \
    --vntr-bed ${params.vntr_severus} \
    --phasing-vcf ${data.dv_vcf} \
    --threads ${task.cpus} \
    --use-supplementary-tag \
    --out-dir ${meta.prefixTN}.severus

    mv ${meta.prefixTN}.severus/somatic_SVs/severus_somatic.vcf ${meta.prefixTN}.severusSomaticSV.vcf    
    """

}

process severus_edits {
    label "low"
    tag "$meta.id"
    
    publishDir "${meta.id}/toolsOutputDNA/severus_somaticSV/", mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: '*.normalAdded.*'
    
    input:
    tuple val(meta), path(data) // severus output vcf

    output:
    tuple val(meta), path("${meta.prefixTN}.severusSomaticSV.normalAdded.vcf.gz"),emit:vcf

    tuple val(meta), path("${meta.prefixTN}.severusSomaticSV.normalAdded.vcf.gz"),path("${meta.prefixTN}.severusSomaticSV.normalAdded.vcf.gz.tbi"),emit:vcf_idx
    //tuple val(meta), path("${meta.prefixTN}.severusSomaticSV.BEDPE.bed"), emit:bedpe

    script:
    """
    bcftools view -h "${data}" \
    | awk 'BEGIN{OFS="\\t"} /^#CHROM/{\$0=\$0 "\\t${meta.npnNormal}"} {print}' > header_with_normal.txt
    
     bcftools view -H "${data}" \
      | awk 'BEGIN{OFS="\\t"}
        {
          split(\$9,fmt,":");
          fake="";
          for(i=1;i<=length(fmt);i++){
            if(fmt[i]=="GT")  fake=fake "0/0";
            else if(fmt[i]=="AD") fake=fake "0,0";
            else if(fmt[i]=="DP") fake=fake "0";
            else fake=fake ".";
            if(i<length(fmt)) fake=fake ":";
          }
          print \$0 "\\t" fake;
        }' > body_with_normal.txt

    cat header_with_normal.txt body_with_normal.txt > ${meta.prefixTN}.severusSomaticSV.normalAdded.vcf

    cat header_with_normal.txt body_with_normal.txt | bgzip -c > ${meta.prefixTN}.severusSomaticSV.normalAdded.vcf.gz
    tabix -p vcf ${meta.prefixTN}.severusSomaticSV.normalAdded.vcf.gz
    """
    
}
//    bash ${params.hrdetect_bedpe_script}/ ${meta.prefixTN}.severusSomaticSV.normalAdded.vcf.gz > ${meta.prefixTN}.severusSomaticSV.BEDPE.bed
process cobalt {
    label "high"
    tag "$meta.id"
    conda "${params.hmftools}"
    publishDir "${meta.id}/toolsOutputDNA/hmftools/", mode: 'copy'
  
    input: 
    tuple val(meta), val(data)

    output: 
    tuple val(meta), path("cobalt/"), emit: cobaltDir
    
    script:
    """
    cobalt "-Xmx16G" \
    -reference ${meta.npnNormal} \
    -reference_bam ${data.bamNormal} \
    -tumor ${meta.npnTumor} \
    -tumor_bam ${data.bamTumor} \
    -ref_genome ${params.genome_fasta} \
    -output_dir cobalt \
    -threads ${task.cpus} \
    -gc_profile ${params.hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp

    """
}
//    mv cobalt/${meta.npnTumor}.cobalt.gc.median.tsv cobalt/${meta.prefixTN}.cobalt.gc.median.tsv

process amber {
    label "medium"
    tag "$meta.id"
    conda "${params.hmftools}"

    publishDir "${meta.id}/toolsOutputDNA/hmftools/", mode: 'copy'

    input: 
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("amber/"), emit: amberDir
    tuple val(meta), path("amber/*.qc"), path("amber/*.baf.tsv.gz"),emit: for_yaml_summary   
    script:
    """
    amber "-Xmx16G" \
    -reference ${meta.npnNormal} \
    -reference_bam ${data.bamNormal} \
    -tumor ${meta.npnTumor} \
    -tumor_bam ${data.bamTumor} \
    -ref_genome ${params.genome_fasta} \
    -output_dir amber \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -loci ${params.hmftools_data_dir_v534}/copy_number/GermlineHetPon.38.vcf.gz
    """
}
//    #mv amber/${meta.npnTumor}.amber.qc  amber/${meta.prefixTN}.amber.qc

process purple {
    label "high"
    tag "$meta.id"
    conda "${params.purple}"

    publishDir "${meta.id}/toolsOutputDNA/hmftools/", mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.purity.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.segment.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.{html,png}"
    //publishDir "${outputDir}/${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.circos.png"
    //publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.driver.catalog.somatic.tsv"


    input:
    tuple val(meta), path(amber), path(cobalt), path(somatic_sv), path(somatic_smallvar)

    output:
    tuple val(meta), path("purple/"), emit: purpleDir
    tuple val(meta), path("purple/${meta.prefixTN}.purple.cnv.somatic.tsv"), emit: purple_pass_for_hrd
    path("${meta.prefixTN}.*")
    tuple val(meta), path("${meta.prefixTN}.purple.cnv.somatic.forPCGR.txt"), emit:cna_for_pcgr
    tuple val(meta), path("purple/${meta.prefixTN}.purple.purity.tsv"), path("purple/${meta.prefixTN}.purple.driver.catalog.somatic.tsv"), path("purple/${meta.prefixTN}.purple.cnv.somatic.tsv"), emit: for_yaml_summary
    script:
    """
    purple "-Xmx16G" \
    -reference ${meta.npnNormal} \
    -tumor ${meta.npnTumor} \
    -ref_genome ${params.genome_fasta} \
    -output_dir purple \
    -threads ${task.cpus} \
    -ref_genome_version 38 \
    -somatic_sv_vcf ${somatic_sv} \
    -somatic_vcf ${somatic_smallvar} \
    -gc_profile ${params.hmftools_data_dir_v534}/copy_number/GC_profile.1000bp.38.cnp \
    -ensembl_data_dir ${params.hmftools_data_dir_v534}/common/ensembl_data/ \
    -amber ${amber} \
    -cobalt ${cobalt} \
    -driver_gene_panel ${params.hmftools_data_dir_v534}/common/DriverGenePanel.38.tsv \
    -somatic_hotspots ${params.hmftools_data_dir_v534}/variants/KnownHotspots.somatic.38.vcf.gz \
    -circos /lnx01_data3/shared/programmer/miniconda3/envs/circos_purple/bin/circos


    mv purple/${meta.npnTumor}.purple.cnv.somatic.tsv purple/${meta.prefixTN}.purple.cnv.somatic.tsv
    mv purple/${meta.npnTumor}.purple.qc purple/${meta.prefixTN}.purple.qc
    mv purple/${meta.npnTumor}.purple.purity.tsv purple/${meta.prefixTN}.purple.purity.tsv
    mv purple/plot/${meta.npnTumor}.circos.png purple/${meta.prefixTN}.purple.PASS.circos.png

    mv purple/${meta.npnTumor}.purple.sv.vcf.gz purple/${meta.prefixTN}.purple.sv.vcf.gz
    mv purple/${meta.npnTumor}.purple.sv.vcf.gz.tbi purple/${meta.prefixTN}.purple.sv.vcf.gz.tbi

    mv purple/${meta.npnTumor}.purple.somatic.vcf.gz purple/${meta.prefixTN}.purple.somatic.vcf.gz
    mv purple/${meta.npnTumor}.purple.somatic.vcf.gz.tbi purple/${meta.prefixTN}.purple.somatic.vcf.gz.tbi
    
    mv purple/${meta.npnTumor}.purple.driver.catalog.somatic.tsv purple/${meta.prefixTN}.purple.driver.catalog.somatic.tsv

    mv purple/${meta.npnTumor}.purple.segment.tsv purple/${meta.prefixTN}.purple.segment.tsv

    cp purple/${meta.prefixTN}.* .

    cut -f 1,2,3,15,16 ${meta.prefixTN}.purple.cnv.somatic.tsv > ${meta.prefixTN}.purple.cnv.somatic.forPCGR.txt
    sed -i "1 s/majorAlleleCopyNumber/nMajor/;s/minorAlleleCopyNumber/nMinor/;s/chromosome/Chromosome/;s/start/Start/;s/end/End/" ${meta.prefixTN}.purple.cnv.somatic.forPCGR.txt
    """

}

process owl_msi {
    label "low"
    tag "$meta.id"

    publishDir "${meta.id}/toolsOutputDNA/owl_MSI/", mode: 'copy'

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.txt"), emit: owl_msi
    tuple val(meta), path("*.owl-scores.txt"), emit: for_yaml_summary
    script:
    """
    ${params.owl} profile \
    --bam ${data.bam} \
    --regions ${params.owl_markers} \
    --sample ${meta.npn} > ${meta.prefix}.owlMSI.txt

    ${params.owl} score \
    --file ${meta.prefix}.owlMSI.txt \
    --prefix ${meta.prefix}_MSI
    """
}

/*
process chord_hrd {
    label "low"
    tag "$meta.id"
    conda "${params.chord}"

    publishDir "${meta.id}/toolsOutputDNA/CHORD_HRD/", mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy', pattern: "*.chord_hrd.summary.txt"

    input:
    tuple val(meta), path(data)
    // data.deepSomaticVCF  : deepSomatic PASS VCF (from deepSomatic_edits)
    // data.severusVCF      : severus somatic SV VCF (from severus_edits)

    output:
    tuple val(meta), path("${meta.prefixTN}.chord_hrd.txt"),         emit: chord_full
    tuple val(meta), path("${meta.prefixTN}.chord_hrd.summary.txt"), emit: for_yaml_summary

    script:
    def (deepSomaticVCF, severusVCF) = data

    """
    Rscript ${params.chord_Rscript} \
        $deepSomaticVCF \
        $severusVCF \
        ${meta.npnTumor} \
        ${meta.prefixTN}
    """
}
*/
process hrd_scores {
    label "low"
    tag "$meta.id"
    conda "${params.sigrap02}"

    publishDir "${meta.id}/toolsOutputDNA/hrd_scores_sigrap/", mode: 'copy'

    input:
    tuple val(meta), path(data)
    // data.deepSomaticVCF  : deepSomatic PASS VCF (from deepSomatic_edits)
    // data.severusVCF      : severus somatic SV VCF (from severus_edits)

    output:
    tuple val(meta), path("${meta.prefixTN}.HRdetect.sigrap.json.gz"),        emit: hrdetect_json
    tuple val(meta), path("${meta.prefixTN}.CHORD.sigrap.json.gz"),           emit: chord_json

    tuple val(meta), path("${meta.prefixTN}.mutationalPattern/"),             emit: mutpat_dir
    tuple val(meta), path("${meta.prefixTN}_snv2020.json.gz"),                emit: snv2020_json
    tuple val(meta), path("${meta.prefixTN}_snv2015.json.gz"),                emit: snv2015_json
    tuple val(meta), path("${meta.prefixTN}_indel.json.gz"),                  emit: indel_json
    tuple val(meta), path("${meta.prefixTN}_dbs.json.gz"),                    emit: dbs_json



    script:
    def (deepSomaticVCF, severusVCF,purpleCNV) = data

    """

    sigrap.R hrdetect \
    --sample ${meta.prefixTN} \
    --snv $deepSomaticVCF \
    --sv $severusVCF \
    --cnv $purpleCNV \
    --out ${meta.prefixTN}.HRdetect.sigrap.json.gz

    sigrap.R chord \
    --sample ${meta.prefixTN} \
    --snv $deepSomaticVCF \
    --sv $severusVCF \
    --out ${meta.prefixTN}.CHORD.sigrap.json.gz

    sigrap.R mutpat \
    --sample ${meta.prefixTN} \
    --snv $deepSomaticVCF \
    --out ${meta.prefixTN}.mutationalPattern

    mv ${meta.prefixTN}.mutationalPattern/sigs/snv2020.json.gz ${meta.prefixTN}_snv2020.json.gz
    mv ${meta.prefixTN}.mutationalPattern/sigs/snv2015.json.gz ${meta.prefixTN}_snv2015.json.gz
    mv ${meta.prefixTN}.mutationalPattern/sigs/indel.json.gz ${meta.prefixTN}_indel.json.gz
    mv ${meta.prefixTN}.mutationalPattern/sigs/dbs.json.gz ${meta.prefixTN}_dbs.json.gz
    """
}

process scarhrd_purple {
    label "low"
    tag "$meta.id"
    conda "${params.scarhrd}"

    publishDir "${meta.id}/toolsOutputDNA/scarHRD/", mode: 'copy'

    input:
    tuple val(meta), path(purpleCNV)
    // purpleCNV: purple *.purple.cnv.somatic.tsv (from purple.out.purple_pass_for_hrd)

    output:
    tuple val(meta), path("${meta.prefixTN}.purple.scarHRD.txt"),         emit: scarhrd_full
    tuple val(meta), path("${meta.prefixTN}.purple.scarHRD.summary.txt"), emit: for_yaml_summary

    script:
    """
    Rscript ${params.scarhrd_Rscript_v2} \
        $purpleCNV \
        ${meta.npnTumor} \
        ${meta.prefixTN}.purple
    """
}

process scarhrd_wakhan {
    label "low"
    tag "$meta.id"
    conda "${params.scarhrd}"
    publishDir "${meta.id}/toolsOutputDNA/scarHRD/", mode: 'copy'

    input:
    tuple val(meta), path(wakhanVCF)

    output:
    tuple val(meta), path("${meta.prefixTN}.wakhan.scarHRD.txt"),         emit: scarhrd_full
    tuple val(meta), path("${meta.prefixTN}.wakhan.scarHRD.summary.txt"), emit: for_yaml_summary

    script:
    def wakhanPloidy = meta.wakhanPloidy ?: 'NA'
    def mincnq = params.scarhrd_wakhan_minCNQ ? "--min-cnq ${params.scarhrd_wakhan_minCNQ}" : ""
    """
    python3 ${params.wakhan_scarhrd_py} \
        --vcf ${wakhanVCF} \
        --sample ${meta.npnTumor} \
        --out ${meta.prefixTN}.wakhan.scarHRD_input.tsv \
        ${mincnq}

    Rscript ${params.scarhrd_Rscript_v2} \
        --input  ${meta.prefixTN}.wakhan.scarHRD_input.tsv \
        --sample ${meta.npnTumor} \
        --out    ${meta.prefixTN}.wakhan \
        --source table \
        --ploidy ${wakhanPloidy}
    """
}
/*
process hrdetect_hrd {
    label "low"
    tag "$meta.id"
    conda "${params.hrdetect}"

    publishDir "${meta.id}/toolsOutputDNA/HRD/",        mode: 'copy'
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/",    mode: 'copy', pattern: "*.hrdetect.summary.txt"

    input:
    tuple val(meta), val(data)
    // data.deepSomaticVCF  : deepSomatic PASS VCF
    // data.severusVCF      : severus somatic SV VCF (normalAdded)
    // data.purpleCNV       : purple CNV somatic TSV

    output:
    tuple val(meta), path("${meta.prefixTN}.hrdetect.txt"),         emit: hrdetect_full
    tuple val(meta), path("${meta.prefixTN}.hrdetect.summary.txt"), emit: hrdetect_summary

    script:

    """
    # Step 1: convert severus VCF to BEDPE
    bash ${params.hrdetect_bedpe_script} ${data.severusVCF} ${meta.prefixTN}.severus.bedpe

    # Step 2: run HRDetect
    Rscript ${params.hrdetect_Rrscript} \
        ${data.deepSomaticVCF} \
        ${meta.prefixTN}.severus.bedpe \
        ${data.purpleCNV} \
        ${meta.npnTumor} \
        ${meta.prefixTN}
    """
}

*/
process pcgr_v212_deepSomatic {
    tag "$meta.id"
    label 'medium'
    conda "${params.pcgr212}"

    publishDir "${meta.id}/toolsOutputDNA/PCGR212/", mode: 'copy', pattern: "*.pcgr.*"
    publishDir "${meta.id}/TUMORBOARDFILES/",mode: 'copy', pattern:"*.html"

    input:
    tuple val(meta),  path(data)        //meta, data: [vcf,idx,cna],

    output:
    path("*.pcgr.*")
    tuple val(meta), path("${meta.id}_pcgr212.pcgr.*.xlsx"), emit: for_yaml_summary
    script:

    //def rnaexp=!params.skipRNA       ? "--input_rna_expression ${data[2]}" : ""
    """

    pcgr \
    --input_vcf ${data[0]} \
    --refdata_dir ${params.pcgr_data_dir3} \
    --output_dir . \
    --vep_dir ${params.pcgr_VEP} \
    --genome_assembly ${params.pcgr_assembly} \
    --sample_id ${meta.id}_pcgr212 \
    --min_mutations_signatures 100 \
    --all_reference_signatures \
    --estimate_tmb \
    --tmb_display coding_non_silent \
    --estimate_msi \
    --exclude_dbsnp_nonsomatic \
    --assay WGS \
    --input_cna ${data[2]} \
    --pcgrr_conda ${params.pcgrr212} \
    --estimate_signatures \
    --tumor_site ${meta.pcgr}
    """
}

process collect_clinical_summary {
    label "low"
    tag "$meta.id"
    conda "${params.somaticSummaryEnv}"  // needs pyyaml, pandas, python-calamine

    publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.clinical_summaryDNA.html"
    publishDir "${meta.id}/summaryFiles/", mode: 'copy', pattern: "*.clinical_summaryDNA.*"
   
    publishDir "${params.lrsStorageBase}/clinicalSummaries/dna/json/", mode: 'copy', pattern: "*.clinical_summaryDNA.json"
    publishDir "${params.lrsStorageBase}/clinicalSummaries/dna/yaml/", mode: 'copy', pattern: "*.clinical_summaryDNA.yaml"

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.prefixTN}.clinical_summaryDNA.yaml"), emit: yaml
    tuple val(meta), path("${meta.prefixTN}.clinical_summaryDNA.json"), emit: json
    tuple val(meta), path("${meta.prefixTN}.clinical_summaryDNA.html"), emit: html
    script:
    """
    python3 ${params.clinical_summary_py} \
        --pcgr-code          ${meta.pcgr} \
        --case-id            ${meta.id} \
        --npn-tumor          ${meta.npnTumor} \
        --npn-normal         ${meta.npnNormal} \
        --gender             ${meta.gender} \
        --genome-version     ${params.genome_version} \
        --amber              ${data.amberQC} \
        --cramino-tumor      ${data.cramino_t} \
        --cramino-normal     ${data.cramino_n} \
        --owl-tumor          ${data.owl_t} \
        --owl-normal         ${data.owl_n} \
        --purple-drivers     ${data.purple_driver} \
        --purple-purity      ${data.purple_purity} \
        --cnv-plot-png       ${data.cnv_plot} \
        --pcgr-xlsx          ${data.pcgr} \
        --met-tumor          ${data.methbat_t} \
        --met-normal         ${data.methbat_n} \
        --wakhan             ${data.wakhan} \
        --hrdetect-json      ${data.hrdetectJson} \
        --chord-json         ${data.chordJson} \
        --scarhrd            ${data.scarhrd} \
        --mutpattern-snv2020 ${data.snv2020Json} \
        --mutpattern-snv2015 ${data.snv2015Json} \
        --mutpattern-indel   ${data.indelJson} \
        --mutpattern-dbs     ${data.dbsJson} \
        --cosmic-groups      ${params.cosmic_groups} \
        --html-template      ${params.clinical_summary_html} \
        --output             ${meta.prefixTN}.clinical_summaryDNA
    """
}
//  --chord              ${data.chord} \
process purple_genome_view {
    label "low"
    tag "$meta.id"
    conda "${params.somaticSummaryEnv}"  // needs pyyaml, pandas, python-calamine
   
    publishDir "${meta.id}/TUMORBOARDFILES/DNA/", mode: 'copy'

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.prefixTN}.purple_genome_view.html"), emit: html
    tuple val(meta), path("${meta.prefixTN}.purple_genome_view.png"),  emit: png

    script:
    """
    python3 ${params.somaticScripts}/purple_genome_view_v1.py \
        --cnv   ${data.purple_cnv} \
        --amber ${data.amberBAF} \
        --out   ${meta.prefixTN}.purple_genome_view.html \
        --png   ${meta.prefixTN}.purple_genome_view.png \
        --title "${meta.id} · ${params.genome_version}" \
        --dpi 180
    """
}

process wakhan {
    label "high"
    tag "$meta.id"
    conda "${params.wakhan}"  // needs pyyaml, pandas, python-calamine
    publishDir "${meta.id}/toolsOutputDNA/", mode: 'copy'

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("wakhan/"), emit: wakhanDir
    tuple val(meta), path("wakhan/${meta.prefixTN}.wakhan.solutions_ranks.tsv"), emit: wakhanTSV
    tuple val(meta), path("wakhan/solution_1/vcf_output/*_cna_integers.vcf.gz"), emit: vcf

    script:
    """
    wakhan all \
    --target-bam ${data.bamTumor} \
    --reference ${params.genome_fasta} \
    --threads ${task.cpus} \
    --normal-phased-vcf ${data.dv_vcf} \
    --breakpoints ${data.severusVCF} \
    --pdf-enable \
    --ploidy-range ${params.wakhan_ploidy_range} \
    --purity-range ${params.wakhan_purity_range} \
    --genome-name ${meta.prefixTN}.wakhan \
    --out-dir-plots wakhan

    mv wakhan/solutions_ranks.tsv wakhan/${meta.prefixTN}.wakhan.solutions_ranks.tsv
    """
}










/*
process owl_msi {
    label "low"
    tag "$meta.id"
    conda "${params.owl}"

    publishDir "${meta.id}/toolsOutputDNA/MSI/", mode: 'copy'

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("*.txt"), emit: owl_msi

    script:
    """
    ${params.owl} profile \
    --bam ${data.bamTumor} \
    --regions ${params.owl_markers} \
    --sample ${meta.npnT} > ${meta.prefixTumor}.owlMSI.txt

    ${params.owl} score \
    --file ${meta.prefixTumor}.owlMSI.txt \
    --prefix ${meta.prefixTumor}_MSI

    ${params.owl} profile \
    --bam ${data.bamNormal} \
    --regions ${params.owl_markers} \
    --sample ${meta.npnNormal} > ${meta.prefixNormal}.owlMSI.txt

    ${params.owl} score \
    --file  ${meta.prefixNormal}.owlMSI.txt \
    --prefix ${meta.prefixNormal}_MSI
    """
}


*/






























































/*
process pbmm2_align {
    errorStrategy 'ignore'
    tag "$meta.id"
   // publishDir "${meta.id}/alignments/", mode: 'copy',pattern: '*.{bam,bai}'
    //publishDir "${params.outdir}/${runfolder_basename}/fastq_symlinks/", mode: 'link', pattern:'*.{fastq,fq}.gz'
    cpus 24
    maxForks 8
    conda "${params.pbmm2}"

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta), path("${meta.prefix}.pbmm2.bam"), path("${meta.prefix}.pbmm2*bai"),  emit: bam
    
    script:
    """
    pbmm2 align \
    --preset HIFI \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.npn} \
    ${params.genome_mmi} \
    ${data[0]} \
    ${meta.prefix}.pbmm2.bam
    """
}




## Working with val(data) setup

process deepvariant{
    errorStrategy 'ignore'
    tag "$meta.id"
    publishDir "${meta.id}/deepVariant/", mode: 'copy'
    //publishDir "${params.lrsStorageBase}/variants/2025/deepVariant/gvcf", mode: 'copy', pattern: "*.deepVariant.g.vcf.*"    
    cpus 32
    maxForks 5

    input:
    tuple val(meta), val(data) // data from joined ch: [bamN:bamNormal,baiN:baiNormal,bamT:bamTumor,baiT:baiTumor] access as e.g. "data.bamN" 

    output:
    tuple val(meta), path("${meta.npnNormal}.${params.genome_version}.deepVariant.vcf.gz"),path("${meta.npnNormal}.${params.genome_version}.deepVariant.vcf.gz.tbi"), emit: dv_vcf
    tuple val(meta), path("${meta.npnNormal}.${params.genome_version}.deepVariant.g.vcf.gz"), emit: dv_gvcf    
    //path("${meta.npnNormal}.deepvariant.vcf_stats_report.txt")
    """
    singularity run -B ${params.s_bind} ${params.simgpath}/deepvariant180.sif /opt/deepvariant/bin/run_deepvariant \
    --model_type=PACBIO \
    --ref=${params.genome_fasta} \
    --reads=${data.bamN} \
    --output_vcf=${meta.npnNormal}.${params.genome_version}.deepVariant.vcf.gz \
    --output_gvcf=${meta.npnNormal}.${params.genome_version}.deepVariant.g.vcf.gz \
    --num_shards=${task.cpus}
    """    
}


process sawFish2{
    tag "$meta.id"
    errorStrategy 'ignore'
    cpus 16
    maxForks 5
    //publishDir "${params.lrsStorageBase}/structuralVariants/2025/sawfish/", mode: 'copy', pattern:"*.sawfishSV.vcf.*"

    publishDir "${meta.id}/CNV_and_SV/sawfish/", mode: 'copy', pattern: "*.sawfishSV.*"


    conda '/lnx01_data3/shared/programmer/miniconda3/envs/sawfish2/'

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.sawfishSV.*")

    tuple val(meta), path("${meta.prefixNormal}.sawfishSV.vcf.gz"), path("${meta.prefixNormal}.sawfishSV.vcf.gz.tbi"), emit:sv_vcf

    script:
    """
    sawfish discover \
    --threads ${task.cpus} \
    --ref ${params.genome_fasta} \
    --bam ${data.bamN} \
    --cnv-excluded-regions ${params.cnv_exclude_sawfish} \
    --output-dir ${meta.id}.sawfishDiscover 

    sawfish joint-call \
    --threads ${task.cpus} \
    --sample ${meta.id}.sawfishDiscover \
    --output-dir ${meta.id}.germline.sawfishSV 
    
    mv ${meta.id}.germline.sawfishSV/genotyped.sv.vcf.gz ${meta.prefixNormal}.sawfishSV.vcf.gz

    mv ${meta.id}.germline.sawfishSV/genotyped.sv.vcf.gz.tbi ${meta.prefixNormal}.sawfishSV.vcf.gz.tbi
    """
}





*/
