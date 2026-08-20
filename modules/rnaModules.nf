#!/usr/bin/env nextflow
nextflow.enable.dsl = 2


date=new Date().format( 'yyMMdd' )
date2=new Date().format( 'yyMMdd HH:mm:ss' )
user="$USER"
runID="${date}.${user}"


////////////////////////////////////////////
/////// ------- PREPROCESS + ALN ------- ///
////////////////////////////////////////////

process inputFiles_symlinks_ubamRNA {
    label "low"
    publishDir {"${meta.id}/documents/inputSymlinks/"}, mode: 'symlink', pattern: '*.{bam,pbi}'
    
    input:
    tuple val(meta), path(data)   

    output:
    tuple val(meta), path(data)
 
    script:
    """
    """

}

process merge_ubams {
    label "high"
    tag "$meta.id"
    conda "${params.isoseq_pbtk}"

    publishDir {"${meta.id}/inputSymlinks/merged_ubam/"}, mode: 'copy', pattern: '*.bam'

    input:
    tuple val(meta), path(bams)

    output:
    tuple val(meta), path("${meta.id}.${meta.npnRNA}.merged.bam")

    script:
    if (bams instanceof List && bams.size() > 1)
        """
        samtools merge -f -@ ${task.cpus} ${meta.id}.${meta.npnRNA}.merged.bam ${bams.join(' ')}
        """
    else
        """
        cp ${bams} ${meta.id}.${meta.npnRNA}.merged.bam
        """
}


process isoseq_refine_cluster {
    label "high"
    tag "$meta.id"
    conda "${params.isoseq_pbtk}"
    
    publishDir {"${meta.id}/toolsOutputRNA/isoseq/refine_cluster/"}, mode: 'copy', pattern: '*.bam'

    input:
    tuple val(meta), path(data)
    
    output:
    tuple val(meta),
        path("${meta.id}.${meta.npnRNA}.flnc.clustered.bam"),
        path("${meta.id}.${meta.npnRNA}.flnc.clustered.bam.pbi"),  emit: isoseq_flnc_clustered
    
    tuple val(meta), 
        path("${meta.id}.${meta.npnRNA}.flnc.clustered.bam"),       emit: bam_for_fofn
    
    tuple val(meta), 
        path("${meta.id}.${meta.npnRNA}.flnc.refined.bam"),
        path("${meta.id}.${meta.npnRNA}.flnc.refined.bam.pbi"),     emit: isoseq_flnc_refined

    tuple val(meta), 
        path("${meta.id}.${meta.npnRNA}.flnc.refined.*.report.json"), emit: refine_report_json

    script:
    """
    isoseq refine \
    ${data[0]} \
    ${params.kinnex_barcodes12plex} \
    ${meta.id}.${meta.npnRNA}.flnc.refined.bam

    isoseq cluster2 \
    ${meta.id}.${meta.npnRNA}.flnc.refined.bam \
    ${meta.id}.${meta.npnRNA}.flnc.clustered.bam

    pbindex ${meta.id}.${meta.npnRNA}.flnc.clustered.bam
    pbindex ${meta.id}.${meta.npnRNA}.flnc.refined.bam 
    """
}

process pbmm2_align_clusteredFLNC {
    label "high"
    tag "$meta.id"
    conda "${params.pbmm2}"

    publishDir {"${meta.id}/alignments/"}, mode: 'copy', pattern: '*.pbmm2.*'
    //publishDir {"${meta.id}/RNA/pbmm2/clustered/"}, mode: 'copy', pattern: '*.pbmm2.*'
    input:
    tuple val(meta), path(bam), path(pbi)
    
    output:
    tuple val(meta), 
        path("${meta.prefixRNA}.clustered.pbmm2.bam"),
        path("${meta.prefixRNA}.clustered.pbmm2*bai"),  emit: bam
 
    script:
    """
    pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.npn} \
    ${params.genome_mmi} \
    ${bam} \
    ${meta.prefixRNA}.clustered.pbmm2.bam
    """
}

process pbmm2_align_refinedFLNC {
    label "high"
    tag "$meta.id"
    conda "${params.pbmm2}"

    publishDir {"${meta.id}/alignments/"}, mode: 'copy', pattern: '*.pbmm2.*'
    
    input:
    tuple val(meta), path(bam), path(pbi)
    
    output:
    tuple val(meta), 
        path("${meta.prefixRNA}.refined.pbmm2.bam"), 
        path("${meta.prefixRNA}.refined.pbmm2*bai"),  emit: bam
 
    script:
    """
    pbmm2 align \
    --preset ISOSEQ \
    --sort \
    --num-threads ${task.cpus} \
    --bam-index BAI \
    --sample ${meta.npn} \
    ${params.genome_mmi} \
    ${bam} \
    ${meta.prefixRNA}.refined.pbmm2.bam
    """
}


process isoseq_collapse {
    label "high"
    tag "$meta.id"
    conda "${params.isoseq_pbtk}"
    
    publishDir {"${meta.id}/toolsOutputRNA/isoseq/collapsed"}, mode: 'copy'

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.collapsed.*"), emit: all
    
    tuple val(meta),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.gff"),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.flnc_count.txt"),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.abundance.txt"), emit: collapsed_list
    
    tuple val(meta),
        path("${meta.prefixRNA}.refinedFLNC.collapsed.read_stat.txt"), emit: read_stat


    script:
    """
    isoseq collapse \
    --do-not-collapse-extra-5exons \
    ${data.refinedPbmm2BAM} \
    ${data.refinedBAM}\
    ${meta.prefixRNA}.refinedFLNC.collapsed.gff
    """
}

process pigeon_classify {
    label "high"
    tag "$meta.id"
    conda "${params.pigeon}"

    publishDir {"${meta.id}/toolsOutputRNA/isoseq/pigeon"}, mode: 'copy'
    
    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), 
        path("*.pigeon*"),                           emit: pigeon
    
    tuple val(meta), 
        path("*.filtered_lite_classification.txt"),  emit: classification

    tuple val(meta), 
        path("*.pigeon_classification.txt"),         emit: classification_unfiltered

    tuple val(meta),
        path("*.collapsed.sorted.gff"),              emit: sortedGFF
    
    tuple val(meta),
        path("*.filtered.report.json"),
        path("*.pigeon.report.json"),                emit: pigeon_reports_json
   
    script:
    """
    cp ${data.collapsedGFF} ${meta.prefixRNA}.pbmm2.collapsed.gff 
    pigeon prepare ${meta.prefixRNA}.pbmm2.collapsed.gff

    pigeon classify \
    ${meta.prefixRNA}.pbmm2.collapsed.sorted.gff \
    ${params.pigeon_gtf} \
    ${params.genome_fasta} \
    --fl ${data.flncCounts} \
    --cage-peak ${params.pigeon_refTSS} \
    --poly-a ${params.pigeon_polyA} \
    -o ${meta.prefixRNA}.pigeon

    pigeon filter \
    ${meta.prefixRNA}.pigeon_classification.txt \
    ${meta.prefixRNA}.pbmm2.collapsed.sorted.gff

    pigeon report \
    --exclude-singletons \
    ${meta.prefixRNA}.pigeon_classification.filtered_lite_classification.txt \
    ${meta.prefixRNA}.pigeon.report.txt

    """
}

process sqanti3_QC {
    label "high"
    tag "$meta.id"
    conda "${params.sqanti3}"

    publishDir {"${meta.id}/toolsOutputRNA/sqanti3/"}, mode: 'copy'
    
    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("sqanti3_results/*"),  emit: sqanti3QC
 
    script:
   // def (refined_bam,refined_pbi,pbmm2_bam,pbmm2_bai, collapsed_gff,flnccounts,abundance) = data
    """
    singularity run -B ${params.s_bind} ${params.simgpath}/sqanti3.sif sqanti3_qc.py \
    --isoforms ${data.collapsedGFF} \
    --refGTF ${params.gencode_gtf} \
    --refFasta ${params.genome_fasta} \
    --CAGE_peak ${params.pigeon_refTSS} \
    --polyA_motif_list ${params.pigeon_polyA} \
    --fl_count ${data.flncAbundance} \
    --saturation \
    --tusco human \
    -d sqanti3_results \
    -n 4
    """
}


process pbfusion {
    label "medium"
    tag "$meta.id"
    conda "${params.pbfusion}"

    publishDir {"${meta.id}/toolsOutputRNA/pbfusion"}, mode: 'copy'
   // publishDir {"${meta.id}/TUMORBOARDFILES/"}, mode: 'copy',pattern: "*.INHOUSE.*"
   

    input:
    tuple val(meta), val(data)
    
    output:
    tuple val(meta), path("*.{pdf,bed,txt,vcf,idx}"),  emit: fusion
    tuple val(meta), path("${meta.prefixRNA}.refinedBAM.PBfusion.INHOUSE.txt"),  emit: inhouse_fusion 
    script:
    //def (refined_bam,refined_pbi,pbmm2_bam,pbmm2_bai) = data
    """
    pbfusion discover \
    -b ${data.refinedBAM} \
    --threads ${task.cpus} \
    --gtf ${params.gencode_gtf} \
    --min-coverage 2 \
    -o  ${meta.prefixRNA}.refinedBAM.fusion

    cat  ${meta.prefixRNA}.refinedBAM.fusion.breakpoints.groups.bed| grep -w -f ${params.inhouse_fusionGenelist} > ${meta.prefixRNA}.refinedBAM.PBfusion.INHOUSE.txt

    """
}


process isocallProfile {

    label "medium"
    tag "$meta.id"

    publishDir {"${meta.id}/toolsOutputRNA/isoCall/profile/"}, mode: 'copy'
    publishDir {"${params.lrsStorageBase}/RNA/isoCallProfiles/"}, mode: 'copy'
    input:
    tuple val(meta), path(bam),path(bai)
    
    output:
    tuple val(meta), path("${meta.prefixRNA}.isoCallProfile.gz"),  emit: profile
 
    script:
    //def (refined_bam,refined_pbi,pbmm2_bam,pbmm2_bai) = data
    """
    ${params.isocall} profile \
    --reads ${bam} \
    --output ${meta.prefixRNA}.isoCallProfile.gz
    """

}

process isocallMerge {

    label "medium"
    tag "$meta.id"

    publishDir {"${meta.id}/toolsOutputRNA/isoCall/"}, mode: 'copy'
    
    input:
    tuple val(meta), path(bam),path(bai)
    
    output:
    tuple val(meta), path("*.{pdf,bed,txt,vcf,idx}"),  emit: fusion
 
    script:
    //def (refined_bam,refined_pbi,pbmm2_bam,pbmm2_bai) = data
    """
    ${params.isocall} merge \
    --reads ${bam} \
    -o TODO
    """

}

process isocallCall {

    label "medium"
    tag "$meta.id"

    publishDir {"${meta.id}/toolsOutputRNA/isoCall/"}, mode: 'copy'
    
    input:
    tuple val(meta), path(profile)
    
    output:
    tuple val(meta), path("*.{gz,txt}"),  emit: fusion
 
    script:
    //def (refined_bam,refined_pbi,pbmm2_bam,pbmm2_bai) = data
    """
    ${params.isocall} call \
    --merged-profile ${profile} \
    --known-isoforms ${params.isocall_gtf} \
    --reference ${params.genome_fasta} \
    --output-prefix ${meta.prefixRNA}.isoCall  
    """

}

// --------------- RNA EXPRESSION ONLY -------------------

process oarFish {
    label "medium"
    tag "$meta.id"
    conda "${params.oarfish}"

    publishDir {"${meta.id}/toolsOutputRNA/oarFish/"}, mode: 'copy', pattern: '*.oarFish.*'

    input:
    tuple val(meta), val(data)
   //     tuple val(meta), path(bam),path(bai)
    output:
    tuple val(meta), path("${meta.prefixRNA}.oarFish.quant"),path("${meta.prefixRNA}.oarFish.ambig_info.tsv"),  emit: quant
    tuple val(meta), path("${meta.prefixRNA}.oarFish.meta_info.json"), emit: meta_json

    script:
    """
    samtools collate -@ ${task.cpus} ${data.refinedPbmm2BAM} -o ${meta.prefixRNA}.oarFish.collated.bam
    
    oarfish \
        --genome-alignments ${meta.prefixRNA}.oarFish.collated.bam \
        --genome-fasta ${params.genome_fasta} \
        --annotation ${params.oarfish_gtf} \
        --threads ${task.cpus} \
        --output ${meta.prefixRNA}.oarFish
    """
}


// -------------- INTEGRATION PROCESSES ------------------ 
process whatshap_haplotag {
    label 'medium'
    tag "$meta.id"
    conda "${params.whatshap}"                    // NB: env must also provide samtools

    publishDir "${meta.id}/toolsOutputRNA/haplotag/",    mode: 'copy', pattern: "*.haplotag.tsv"
    publishDir "${meta.id}/alignments/",  mode: 'copy', pattern: "*.haplotagged.ba*"

    input:
    tuple val(meta), val(data)   // data: bam, phasedVcf   (phasedVcf .tbi sits beside it)

    output:
    tuple val(meta), path("${meta.prefixRNA}.refined.pbmm2.haplotagged.bam"),
                     path("${meta.prefixRNA}.refined.pbmm2.haplotagged.bam.bai"), emit: bam
    tuple val(meta), path("${meta.prefixRNA}.refined.pbmm2.haplotag.tsv"),        emit: haplotag_list

    script:
    """
    whatshap haplotag \
        --reference ${params.genome_fasta} \
        --ignore-read-groups \
        --output-haplotag-list ${meta.prefixRNA}.refined.pbmm2.haplotag.tsv \
        --output ${meta.prefixRNA}.refined.pbmm2.haplotagged.bam \
        ${data.phasedVcf} \
        ${data.bam}

    samtools index ${meta.prefixRNA}.refined.pbmm2.haplotagged.bam
    """

}

process collect_clinical_summaryRNA {
    label "low"
    tag "$meta.id"
    conda "${params.somaticSummaryEnv}"  // needs pyyaml, pandas, python-calamine

    publishDir "${meta.id}/TUMORBOARDFILES/", mode: 'copy', pattern: "*.clinical_summaryRNA.html"
    publishDir "${meta.id}/summaryFiles/", mode: 'copy', pattern: "*.clinical_summaryRNA.*"
   
    publishDir "${params.lrsStorageBase}/clinicalSummaries/rna/json/", mode: 'copy', pattern: "*.clinical_summaryRNA.json"
    publishDir "${params.lrsStorageBase}/clinicalSummaries/rna/yaml/", mode: 'copy', pattern: "*.clinical_summaryRNA.yaml"

    input:
    tuple val(meta), val(data)

    output:
    tuple val(meta), path("${meta.prefixTN}.clinical_summaryRNA.yaml"), emit: yaml
    tuple val(meta), path("${meta.prefixTN}.clinical_summaryRNA.json"), emit: json
    tuple val(meta), path("${meta.prefixTN}.clinical_summaryRNA.html"), emit: html
    script:
    """
    python3 ${params.clinical_summaryRNA_py} \
        --case-id            ${meta.id} \
        --npn-rna            ${meta.npn} \
        --genome-version     ${params.genome_version} \
        --refine-json        ${data.refineReportJSON} \
        --pigeon-raw         ${data.pigeonRawJSON} \
        --pigeon-filtered    ${data.pigeonFilteredJSON} \
        --pbfusion           ${data.fusionInhouse} \
        --html-template      ${params.clinical_summaryRNA_html} \
        --output             ${meta.prefixTN}.clinical_summaryRNA
    """
}


//  --pcgr-code          ${meta.pcgr} \

/* ----------------------------------------------------------------------------
 * ase_readcounter — allele-specific expression on the tumour RNA.
 *   run 1: germline biallelic het SNPs from the NORMAL (constitutional baseline,
 *          avoids tumour CNV/LOH artefacts) -> genome-wide allelic imbalance
 *   run 2: somatic SNVs (DeepSomatic PASS, tumour-only) -> is the mutant allele
 *          transcribed / preferentially expressed / NMD-degraded
 * Self-contained (GATK ASEReadCounter + bcftools). No external R.
 * -------------------------------------------------------------------------- */
process ase_readcounter {
    label 'medium'
    tag "$meta.id"
    conda "${params.aseGatkEnv}"                  // gatk4 + bcftools + htslib(tabix)

    publishDir "${meta.id}/toolsOutputRNA/ASE/", mode: 'copy', pattern: "*.ASE.*.tsv"

    input:
    tuple val(meta), val(data)
    // data: bam, germlineVcf, somaticVcf  (somaticVcf may be an empty list -> skipped)

    output:
    tuple val(meta), path("${meta.prefixRNA}.ASE.germlineHet.tsv"), emit: germline_ase
    tuple val(meta), path("${meta.prefixRNA}.ASE.somatic.tsv"),     emit: somatic_ase, optional: true

    script:
    def somaticVcf = (data.somaticVcf instanceof List) ? '' : "${data.somaticVcf}"
    """
    # --- germline biallelic het SNPs from the NORMAL ---
    bcftools view -f PASS -m2 -M2 -v snps -g het \
        ${data.germlineVcf} -Oz -o germline.hets.vcf.gz
    tabix -p vcf germline.hets.vcf.gz

    gatk ASEReadCounter \
        -R ${params.genome_fasta} \
        -I ${data.bam} \
        -V germline.hets.vcf.gz \
        --min-mapping-quality ${params.ase_minMapQ} \
        --min-base-quality ${params.ase_minBaseQ} \
        -O ${meta.prefixRNA}.ASE.germlineHet.tsv

    # --- somatic SNVs (tumour-only PASS) ---
    if [ -n "${somaticVcf}" ] && [ -s "${somaticVcf}" ]; then
        bcftools view -m2 -M2 -v snps \
            ${somaticVcf} -Oz -o somatic.snvs.vcf.gz
        tabix -p vcf somatic.snvs.vcf.gz

        gatk ASEReadCounter \
            -R ${params.genome_fasta} \
            -I ${data.bam} \
            -V somatic.snvs.vcf.gz \
            --min-mapping-quality ${params.ase_minMapQ} \
            --min-base-quality ${params.ase_minBaseQ} \
            -O ${meta.prefixRNA}.ASE.somatic.tsv || true
    fi
    """
}


/* ----------------------------------------------------------------------------
 * expression_outlier — n=1 "DE": gene-level FL counts vs a background cohort.
 * WRAPPER around your existing expression_outlier.R (from the somatic-RNA chat).
 * Emits <sample>.geneCounts.tsv to a shared cohort dir so the background grows.
 *
 * [INFER] CLI below is reconstructed from the handoff. Reconcile with the
 *         script's actual argparse before first run.
 * -------------------------------------------------------------------------- */
process expression_outlier {
    label 'low'
    tag "$meta.id"
    conda "${params.rnaOUTRIDER}"                        // R: OUTRIDER (+ deps)

    publishDir "${meta.id}/toolsOutputRNA/expression/",              mode: 'copy', pattern: "*.outliers.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/",         mode: 'copy', pattern: "*.outliers.tsv"
    publishDir "${params.exprCohortDir}/",                mode: 'copy', pattern: "*.geneCounts.tsv"

    input:
    tuple val(meta), val(data)         // pigeon filtered_lite_classification.txt

    output:
    tuple val(meta), path("${meta.prefixRNA}.outliers.tsv"),   emit: outliers,   optional: true
    tuple val(meta), path("${meta.prefixRNA}.geneCounts.tsv"), emit: geneCounts, optional: true

    script:
    """
    Rscript ${params.expression_outlier_R} \
    --sample         ${meta.prefixRNA} \
    --classification ${data.classification} \
    --abundance      ${data.abundance} \
    --background     ${params.exprBackground} \
    --min-cohort     ${params.expr_minCohort} \
    --out            ${meta.prefixRNA}
    """
}


/* ----------------------------------------------------------------------------
 * splicing_isoformswitch — exon skipping etc. via IsoformSwitchAnalyzeR in
 * ANNOTATION-ONLY mode (n=1: no switch test). WRAPPER around your existing
 * isoformSwitch_singleSample.R (from the somatic-RNA chat).
 *
 * [INFER] CLI below is reconstructed from the handoff. Reconcile with the
 *         script's actual argparse before first run.
 * -------------------------------------------------------------------------- */
process splicing_isoformswitch {
    label 'medium'
    tag "$meta.id"
    conda "${params.rnaISO}"                        // R: IsoformSwitchAnalyzeR (+ BSgenome.Hsapiens.UCSC.hg38)

    publishDir "${meta.id}/toolsOutputRNA/splicing/",        mode: 'copy'
    //publishDir "${meta.id}/TUMORBOARDFILES/RNA/", mode: 'copy', pattern: "*.exonSkipping.panel.tsv"

    input:
    tuple val(meta), val(data)                    // data: gff (pigeon sorted), classification

    output:
    tuple val(meta), path("${meta.prefixRNA}.exonSkipping.panel.tsv"), emit: panel, optional: true
    tuple val(meta), path("${meta.prefixRNA}.splicing*"),              emit: all,   optional: true

    script:
    """
    Rscript ${params.isoformSwitch_R} \
        --sample         ${meta.prefixRNA} \
        --gff            ${data.gff} \
        --abundance      ${data.abundance} \
        --classification ${data.classification} \
        --gtf            ${params.pigeon_gtf} \
        --genelist       ${params.inhouse_splicing_genelist} \
        --min-if         ${params.splicing_minIF} \
        --out            ${meta.prefixRNA}
    """
}


