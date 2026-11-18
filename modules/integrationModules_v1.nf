#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
 integrationModules.nf   (cross-arm DNA <-> RNA integration)
 House convention: file inputs travel inside a val(data) map as absolute paths
================================================================================
*/


// TO DO:
// channel to main
/*
hiPhase.out.hiphase_bam_normal
    | mix(hiPhase.out.hiphase_bam_tumor)
    | map { meta, bam, bai -> [meta, bam.name, bai.name] }
    | set { align_links_ch }

alignmentLinks_tumorboard(align_links_ch)
//symlink process:
process alignmentLinks_tumorboard {
    label 'low'
    tag  "$meta.id"

    input:
    tuple val(meta), val(bamName), val(baiName)

    output:
    tuple val(meta), val(bamName), emit: linked   // token; keeps it in the DAG

    script:
    def realDir = "${launchDir}/${meta.id}/alignments"
    def tbDir   = "${launchDir}/${meta.id}/${params.tumorboardAlignSubdir}"
    """
    mkdir -p '${tbDir}'
    for f in '${bamName}' '${baiName}'; do
        rel=\$(realpath -ms --relative-to='${tbDir}' "${realDir}/\$f")
        ln -sf "\$rel" "${tbDir}/\$f"
    done
    """
}

*/


process methylation_expression {
    label 'low'
    tag "$meta.id"
    conda "${params.condaEnvs.integration}"              // python3 + pandas + numpy

    publishDir "${meta.id}/toolsOutputRNA/integration/methylExpr/", mode: 'copy', pattern: "*.methylExpr*.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/",        mode: 'copy', pattern: "*.methylExpr.panel.tsv"

    input:
    tuple val(meta), val(data)                    // data: geneCounts, normalBed, tumorBed

    output:
    tuple val(meta), path("${meta.prefixRNA}.methylExpr.tsv"),       emit: table
    tuple val(meta), path("${meta.prefixRNA}.methylExpr.panel.tsv"), emit: panel, optional: true

    script:
    def genelist = params.methylExprGenelist ? "--genelist ${params.methylExprGenelist}" : ""
    """
    python3 ${params.methylExprPy} \
        --gtf         ${params.gencodeV49Gtf} \
        --tumor-meth  ${data.tumorBed} \
        --normal-meth ${data.normalBed} \
        --expression  ${data.geneCounts} \
        --sample      ${meta.prefixRNA} \
        --out         ${meta.prefixRNA} \
        --upstream    ${params.promoterUp} \
        --downstream  ${params.promoterDown} \
        --meth-col    ${params.methCol} \
        --cov-col     ${params.methCovCol} \
        --min-cov     ${params.methMinCov} \
        --hyper       ${params.methHyper} \
        --delta       ${params.methDelta} \
        ${genelist}
    """
}


process fusion_sv_concordance {
    label 'low'
    tag "$meta.id"
    conda "${params.condaEnvs.integration}"              // python3 + numpy

    publishDir "${meta.id}/toolsOutputRNA/fusionSV_integration/", mode: 'copy', pattern: "*.fusionSV.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/",      mode: 'copy', pattern: "*.fusionSV.tsv"

    input:
    tuple val(meta), val(data)                    // data: fusionBed, sv (severus vcf)

    output:
    tuple val(meta), path("${meta.prefixRNA}.fusionSV.tsv"), emit: table, optional: true

    script:
    """
    python3 ${params.fusionSvPy} \
        --fusion ${data.fusionBed} \
        --sv     ${data.sv} \
        --sample ${meta.prefixRNA} \
        --out    ${meta.prefixRNA} \
        --window ${params.fusionSvWindow}
    """
}
