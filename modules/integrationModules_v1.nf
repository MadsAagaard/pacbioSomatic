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
    def tbDir   = "${launchDir}/${meta.id}/${params.tumorboard_align_subdir}"
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
    conda "${params.integrationEnv}"              // python3 + pandas + numpy

    publishDir "${meta.id}/toolsOutputRNA/integration/methylExpr/", mode: 'copy', pattern: "*.methylExpr*.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/",        mode: 'copy', pattern: "*.methylExpr.panel.tsv"

    input:
    tuple val(meta), val(data)                    // data: geneCounts, normalBed, tumorBed

    output:
    tuple val(meta), path("${meta.prefixRNA}.methylExpr.tsv"),       emit: table
    tuple val(meta), path("${meta.prefixRNA}.methylExpr.panel.tsv"), emit: panel, optional: true

    script:
    def genelist = params.methylExpr_genelist ? "--genelist ${params.methylExpr_genelist}" : ""
    """
    python3 ${params.methyl_expr_py} \
        --gtf         ${params.gencode_gtf} \
        --tumor-meth  ${data.tumorBed} \
        --normal-meth ${data.normalBed} \
        --expression  ${data.geneCounts} \
        --sample      ${meta.prefixRNA} \
        --out         ${meta.prefixRNA} \
        --upstream    ${params.promoter_up} \
        --downstream  ${params.promoter_down} \
        --meth-col    ${params.meth_col} \
        --cov-col     ${params.meth_cov_col} \
        --min-cov     ${params.meth_minCov} \
        --hyper       ${params.meth_hyper} \
        --delta       ${params.meth_delta} \
        ${genelist}
    """
}


process fusion_sv_concordance {
    label 'low'
    tag "$meta.id"
    conda "${params.integrationEnv}"              // python3 + numpy

    publishDir "${meta.id}/toolsOutputRNA/fusionSV_integration/", mode: 'copy', pattern: "*.fusionSV.tsv"
    publishDir "${meta.id}/TUMORBOARDFILES/RNA/",      mode: 'copy', pattern: "*.fusionSV.tsv"

    input:
    tuple val(meta), val(data)                    // data: fusionBed, sv (severus vcf)

    output:
    tuple val(meta), path("${meta.prefixRNA}.fusionSV.tsv"), emit: table, optional: true

    script:
    """
    python3 ${params.fusion_sv_py} \
        --fusion ${data.fusionBed} \
        --sv     ${data.sv} \
        --sample ${meta.prefixRNA} \
        --out    ${meta.prefixRNA} \
        --window ${params.fusion_sv_window}
    """
}
