#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
================================================================================
 KG Vejle — PacBio LRS   DNA (somatic T/N WGS) + RNA (Kinnex Iso-Seq)
 Unified, case-centric pipeline.
--------------------------------------------------------------------------------
 Samplesheet (tab-sep, >=6 cols; header optional):
   caseID  DNAnormalSampleID  DNAtumorSampleID  RNAtumorSampleID  gender  pcgrCode
   RNAtumorSampleID may be empty / NA / - / .   -> that case is DNA-only.
   --skipRNA disables the RNA arm globally regardless of the samplesheet.

 Entry points (-entry):
   FULL   (default)  DNA + RNA + INTEGRATION      [RNA/INTEGRATION wired next]
   DNA               DNA arm only
   RNA               RNA arm only                 [wired next]
================================================================================
*/

date  = new Date().format('yyMMdd')
date2 = new Date().format('yyMMdd HH:mm:ss')
user  = "$USER"
runID = "${date}.${user}"
 
// -------------------------- naming (single source) ---------------------------
//def samplePrefix = { npn, sampletype -> "${npn}.${sampletype}.${params.genome_version}.${params.readSet}" }

// ---- shared pcgr tag helper (define BEFORE the prefix closures) ----
def pcgrTag = { pcgr ->
    def code = pcgr?.toString()?.trim()
    (code && !(code.toLowerCase() in ['', 'na', '-', '.', 'null'])) ? ".pcgr_${code}" : ''
}

def samplePrefix = { npn, sampletype, pcgr = null ->
    def tag = (sampletype in ['tumor', 'rna']) ? pcgrTag(pcgr) : ''
    "${npn}.${sampletype}${tag}.${params.genome_version}.${params.readSet}"
}

def tnPrefix = { id, pcgr = null ->
    "${id}${pcgrTag(pcgr)}.${params.genome_version}.${params.readSet}"
}
// 
def isMissing = { v -> (v == null) || (v.toString().trim().toLowerCase() in ['', 'na', '-', '.', 'null']) }


/*

def samplePrefix = { npn, sampletype, pcgr = null ->
    def code    = pcgr?.toString()?.trim()
    def useCode = code && !(code.toLowerCase() in ['', 'na', '-', '.', 'null'])
    def tag     = (sampletype in ['tumor', 'rna'] && useCode) ? ".pcgr_${code}" : ''
    "${npn}.${sampletype}${tag}.${params.genome_version}.${params.readSet}"
}

def tnPrefix     = { id  -> "${id}.${params.genome_version}.${params.readSet}" }

//def rnaPrefix    = { id, npn            -> "${npn}.${params.genome_version}.rna" }

def isMissing = { v -> (v == null) || (v.toString().trim().toLowerCase() in ['', 'na', '-', '.', 'null']) }
*/

/* =============================================================================
 *  RAW INPUT GLOBS  (two arms, two filename grammars)
 * ========================================================================== */

// ---- DNA ubams: <archive>/**/<npn>...hifi_reads[.allReads].bam ----
def dnaInputBam
if (params.inputDNA) {
    dnaInputBam = params.allReads ? "${params.inputDNA}/*.bam" : "${params.inputDNA}/*.hifi_reads.*.bam"
} else {
    dnaInputBam = params.allReads ? "${params.dataArchiveDNA}/**/*.bam" : "${params.dataArchiveDNA}/**/*.hifi_reads.*.bam"
}

// ---- RNA demux ubams: <archive>/<run>/rna_demultiplex/<npn>.*.lima.*.SMfixed.bam ----
def rnaInputBam = params.inputRNA ? "${params.inputRNA}/**/*.lima.*.bam"
                                  : "${params.dataArchiveRNA}/**/*.lima.*.bam"


/* =============================================================================
 *  UBAM CHANNELS  (keyed by sample npn, grouped over multiple movies)
 * ========================================================================== */

Channel.fromPath(dnaInputBam, followLinks: true)
    | map { bam ->
        def base = bam.baseName
        def (samplenameFull, pacbioID, readset, barcode) = base.tokenize('.')
        def (samplename, material, testlist, gender)      = samplenameFull.tokenize('_')
        def meta = [ npn: samplename, material: material, testlist: testlist, gender: gender, pacbioID: pacbioID ]
        tuple(meta.npn, bam)
    }
    | groupTuple(sort: true)
    | set { dnaUbam_ch }

Channel.fromPath(rnaInputBam, followLinks: true)
    | map { bam -> tuple(bam.baseName, bam) }
    | map { id, bam ->
        // e.g. 113125099265.m84328_..._s4.hifi_reads.lima.IsoSeqX_bc02_5p.SMfixed
        def npn = id.tokenize('.')[0]
        tuple(npn, bam)
    }
    | groupTuple(sort: true)
    | set { rnaUbam_ch }


/* =============================================================================
 *  SAMPLESHEET  ->  one enriched meta per case  ->  keyed views
 * ========================================================================== */

if (!params.samplesheet) { exit 1, "ERROR: --samplesheet is required (caseID, DNAnormal, DNAtumor, RNAtumor, gender, pcgr)" }

Channel.fromPath(params.samplesheet)
    | splitCsv(sep: '\t')
    | filter { row -> row && row[0] && !row[0].toString().startsWith('#') && row[0].toString().toLowerCase() != 'caseid' }
    | map { row ->
        def (caseID, npnNormal, npnTumor, npnRNA, gender, pcgr) = row
        def hasRNA = !params.skipRNA && !isMissing(npnRNA)
        def meta = [
            id:        caseID,
            npnNormal: npnNormal,
            npnTumor:  npnTumor,
            npnRNA:    hasRNA ? npnRNA.toString().trim() : null,
            gender:    gender,
            pcgr:      pcgr,
            hasRNA:    hasRNA
        ]
        // attach canonical prefixes once
        meta + [
            prefixNormal: samplePrefix(meta.npnNormal, 'normal'),
            prefixTumor:  samplePrefix(meta.npnTumor,  'tumor',meta.pcgr),
            prefixRNA:    hasRNA ? samplePrefix(meta.npnRNA, 'rna',meta.pcgr) : null,
            prefixTN:     tnPrefix(meta.id,meta.pcgr)
        ]
    }
    | set { cases_ch }

// keyed views
cases_ch | map { meta -> [meta.npnNormal, meta] } | set { ss_normal_ch }
cases_ch | map { meta -> [meta.npnTumor,  meta] } | set { ss_tumor_ch  }

// DNA per-sample channel (normal + tumor rows), joined to ubams
ss_normal_ch.join(dnaUbam_ch)
    | map { npn, meta, bams -> [ meta + [sampletype: 'normal', npn: npn, prefix: samplePrefix(meta.npnNormal, 'normal')], bams ] }
    | set { normal_ch }

ss_tumor_ch.join(dnaUbam_ch)
    | map { npn, meta, bams -> [ meta + [sampletype: 'tumor', npn: npn, prefix: samplePrefix(meta.npnTumor, 'tumor',meta.pcgr)], bams ] }
    | set { tumor_ch }

normal_ch.concat(tumor_ch)
    | set { per_sample_ch }

// RNA per-case channel (only hasRNA cases; empty when --skipRNA)
cases_ch
    | filter { meta -> meta.hasRNA }
    | map    { meta -> [meta.npnRNA, meta] }
    | set    { ss_rna_ch }

ss_rna_ch.join(rnaUbam_ch)
    | map { npn, meta, bams -> [ meta + [sampletype: 'rnaTumor', npn: npn, prefix: samplePrefix(meta.npnRNA, 'rna',meta.pcgr)], bams ] }
    | set { rna_ch }


/* =============================================================================
 *  MODULE IMPORTS  
 * ========================================================================== */
include {
    create_fofn;
    inputFiles_symlinks_ubam;
    pbmm2_align_mergedData;
    mosdepthROI;
    multiQC;
    cramino;
    deepvariant; 
    sawFish2; 
    svdb_SawFish;
    hiPhase;
    pbCPGtools;
    methBat;
    methBatNEW_pileup;
    methBatNEW_profile_single;
    deepSomatic; 
    deepSomatic_edits;
    severus;
    severus_edits;
    cobalt;
    amber;
    purple;
    owl_msi;
    //chord_hrd;
    hrd_scores;
    scarhrd_purple;
    //scarhrd_wakhan;
    scarhrd_wakhan_bed;
    pcgr_v212_deepSomatic;
    collect_clinical_summary;
    purple_genome_view;
    wakhan;
    alignmentLinks_tumorboard;
} from './modules/dnaModules.nf'


include {
    inputFiles_symlinks_ubamRNA;
    merge_ubams;
    isoseq_refine_cluster;
    pbmm2_align_clust;
    pbmm2_align_refined_forIsocall;
    isoseq_collapse;
    pigeon_classify;
    sqanti3_QC;
    pbfusion;
    isocallProfile;
    isocallCall;
    oarFish;
    whatshap_haplotag;
    collect_clinical_summaryRNA;
    ase_readcounter;
    expression_outlier;
    splicing_isoformswitch;
} from './modules/rnaModules.nf'


// ---- cross-arm integration (methylation x expression, fusion x SV)
include {
    methylation_expression;
    fusion_sv_concordance;
} from './modules/integrationModules_v1.nf'


/* =============================================================================
 *  DNA SUBWORKFLOWS
 * ========================================================================== */

workflow DNA_PREPROCESS {
    take: per_sample_ch
    main:
        inputFiles_symlinks_ubam(per_sample_ch)
        create_fofn(per_sample_ch)
        pbmm2_align_mergedData(create_fofn.out)

        pbmm2_align_mergedData.out.bam
            | map { meta, bam, bai -> tuple(meta, [bam, bai]) }
            | set { aligned_ch }

        pbmm2_align_mergedData.out.bam
            | map { meta, bam, bai -> [meta.id, meta, bam, bai] }
            | groupTuple(by: 0)
            | map { id, metas, bams, bais ->
                def n = metas.findIndexOf { it.sampletype == 'normal' }
                def t = metas.findIndexOf { it.sampletype == 'tumor'  }
                def meta = [
                    id: id,
                    npnNormal: metas[n].npn,
                    npnTumor:  metas[t].npn,
                    gender:    metas[n].gender,
                    pcgr:      metas[n].pcgr,
                    npnRNA:    metas[n].npnRNA,      // carried through for RNA integration
                    hasRNA:    metas[n].hasRNA
                ]
                meta = meta + [
                    prefixNormal: samplePrefix(meta.npnNormal, 'normal'),
                    prefixTumor:  samplePrefix(meta.npnTumor,  'tumor',meta.pcgr),
                ]
                def data = [ bamNormal: bams[n], baiNormal: bais[n], bamTumor: bams[t], baiTumor: bais[t] ]
                tuple(meta, data)
            }
            | set { tn_aligned_ch }
    emit:
        aligned              = pbmm2_align_mergedData.out.bam
        aligned_per_sample   = aligned_ch
        aligned_tn_paired    = tn_aligned_ch
}

workflow DNA_PREPHASE {
    take: tn_aligned_ch
    main:
        deepvariant(tn_aligned_ch)
        sawFish2(tn_aligned_ch)

        deepvariant.out.dv_vcf | map { meta, vcf, idx -> tuple(meta, [vcf, idx]) } | set { dv_vcf_ch }
        sawFish2.out.sv_vcf    | map { meta, vcf, idx -> tuple(meta, [vcf, idx]) } | set { sawfish_vcf_ch }

        tn_aligned_ch.join(dv_vcf_ch).join(sawfish_vcf_ch) | set { hiphase_input_ch }
    emit:
        dv_vcf                   = dv_vcf_ch
        dv_gvcf                  = deepvariant.out.dv_gvcf
        sawfish_vcf              = sawfish_vcf_ch
        hiphaseInput             = hiphase_input_ch
        sawfish_supporting_reads = sawFish2.out.sv_supporting_reads
}

workflow DNA_PHASE {
    take:
        hiphase_input_ch
        sawfish_supporting_reads
    main:
        hiPhase(hiphase_input_ch)

        hiPhase.out.hiphase_bam_normal
            .join(hiPhase.out.hiphase_bam_tumor)
            .join(hiPhase.out.hiphase_dv_vcf)
            .join(hiPhase.out.hiphase_sv_vcf)
            .join(sawfish_supporting_reads)
            | map { meta, bamN, baiN, bamT, baiT, dv_vcf, dv_idx, sv_vcf, sv_idx, sv_jsonReads ->
                tuple(
                    meta + [
                        prefixNormal: samplePrefix(meta.npnNormal, 'normal'),
                        prefixTumor:  samplePrefix(meta.npnTumor,  'tumor',meta.pcgr),
                        prefixTN:     tnPrefix(meta.id,meta.pcgr)
                    ],
                    [
                        bamNormal: bamN, baiNormal: baiN,
                        bamTumor:  bamT, baiTumor:  baiT,
                        dv_vcf:    dv_vcf, dv_idx: dv_idx,
                        sawfish_vcf: sv_vcf, sawfish_idx: sv_idx,
                        sawfish_reads: sv_jsonReads
                    ]
                )
            }
            | set { phasedAll_ch }

         hiPhase.out.hiphase_bam_normal
         |mix(hiPhase.out.hiphase_bam_tumor)
         | map { meta, bam, bai -> [meta, bam.name, bai.name] }
         | set { align_links_ch }

    alignmentLinks_tumorboard(align_links_ch)

    emit:
        phasedAll = phasedAll_ch
}

workflow DNA_SOMATIC {
    take: phasedAll
    main:
        // per-sample fan-out (normal + tumor) from the paired map
        phasedAll
            | flatMap { meta, data ->
                [
                    [ meta + [sampletype: 'normal', npn: meta.npnNormal, prefix: samplePrefix(meta.npnNormal, 'normal')],
                      [bam: data.bamNormal, bai: data.baiNormal, dv_vcf: data.dv_vcf, dv_idx: data.dv_idx] ],
                    [ meta + [sampletype: 'tumor',  npn: meta.npnTumor,  prefix: samplePrefix(meta.npnTumor,  'tumor',meta.pcgr)],
                      [bam: data.bamTumor,  bai: data.baiTumor,  dv_vcf: data.dv_vcf, dv_idx: data.dv_idx] ]
                ]
            }
            | set { phasedPerSample }

        svdb_SawFish(phasedAll)
        pbCPGtools(phasedPerSample)
        methBat(pbCPGtools.out)

        methBatNEW_pileup(phasedPerSample)
        methBatNEW_profile_single(methBatNEW_pileup.out.met5mC)

        // per-case T/N combined bedMethyl (for INTEGRATION methylation x expression)
        pbCPGtools.out
            | map { meta, files ->
                def fl  = (files instanceof List) ? files : [files]
                def bed = fl.find { it.name.contains('combined.bed') }
                [meta.id, meta.sampletype, bed]
            }
            | filter { it[2] != null }
            | groupTuple(by: 0)
            | map { id, sts, beds -> [id, beds[sts.indexOf('normal')], beds[sts.indexOf('tumor')]] }
            | set { methyl_bed_by_case }
        cramino(phasedPerSample)
        mosdepthROI(phasedPerSample)

        deepSomatic(phasedAll)
        deepSomatic_edits(deepSomatic.out.vcf)
        owl_msi(phasedPerSample)
        cobalt(phasedAll)
        amber(phasedAll)
        severus(phasedAll)
        severus_edits(severus.out.vcf)

        amber.out.amberDir.join(cobalt.out).join(severus_edits.out.vcf).join(deepSomatic_edits.out.vcf)
            | set { purple_pass_input }

        phasedAll.join(severus_edits.out.vcf)
            | map { meta, data, severusVCF -> tuple(meta, data + [severusVCF: severusVCF]) }
            | set { phasedAll_with_severusVCF }



        wakhan(phasedAll_with_severusVCF)
        purple(purple_pass_input)

        // HRD
        
        deepSomatic_edits.out.vcf.join(severus_edits.out.vcf)
            | map { meta, dsVCF, svVCF -> tuple(meta, [dsVCF, svVCF]) }
            | set { chord_input }
        //chord_hrd(chord_input)

        deepSomatic_edits.out.vcf.join(severus_edits.out.vcf).join(purple.out.purple_pass_for_hrd)
            | map { meta, dsVCF, svVCF, purpleCNV -> tuple(meta, [dsVCF, svVCF, purpleCNV]) }
            | set { hrd_input }
        hrd_scores(hrd_input)

/*
        wakhan.out.vcf
        | map { meta, vcf ->
            def m = (vcf.name =~ /_([\d.]+)_([\d.]+)_([\d.]+)_wakhan_cna_integers/)
            def ploidy = m ? m[0][1] : null
            def purity = m ? m[0][2] : null
            tuple(meta + [wakhanPloidy: ploidy, wakhanPurity: purity], vcf)
        }
        | set { wakhan_for_scarHRD }
*/
        wakhan.out.cnBed
        | map { meta, hp1, hp2 ->
            // <genome_name>_<ploidy>_<purity>_<conf>_copynumbers_segments_HP_1.bed
            def m = (hp1.name =~ /_([\d.]+)_([\d.]+)_([\d.]+)_copynumbers_segments_HP_1\.bed$/)
            return tuple(meta + [wakhanPloidy: m ? m[0][1] : null,
                                 wakhanPurity: m ? m[0][2] : null,
                                 wakhanConf:   m ? m[0][3] : null], hp1, hp2)
        }
        | set { wakhan_bed_for_scarHRD }


        scarhrd_purple(purple.out.purple_pass_for_hrd)

       // scarhrd_wakhan(wakhan_for_scarHRD)
        
        scarhrd_wakhan_bed(wakhan_bed_for_scarHRD)

        // PCGR
        deepSomatic.out.pcgr_vcf.join(purple.out.cna_for_pcgr)
            | map { meta, pcgr_vcf, pcgr_idx, cna -> tuple(meta, [pcgr_vcf, pcgr_idx, cna]) }
            | set { pcgr_input }
        pcgr_v212_deepSomatic(pcgr_input)

        // ---------------- summary assembly (unchanged logic) ----------------
        cramino.out.for_yaml_summary
            | map { meta, txt -> [meta.id, meta.sampletype, txt] } | groupTuple(by: 0)
            | map { id, st, txts -> [id, txts[st.indexOf('normal')], txts[st.indexOf('tumor')]] }
            | set { cramino_for_yaml_ch }

        owl_msi.out.for_yaml_summary
            | map { meta, txt -> [meta.id, meta.sampletype, txt] } | groupTuple(by: 0)
            | map { id, st, txts -> [id, txts[st.indexOf('normal')], txts[st.indexOf('tumor')]] }
            | set { owl_for_yaml_ch }

        methBat.out.for_yaml_summary
            | map { meta, j -> [meta.id, meta.sampletype, j] } | groupTuple(by: 0)
            | map { id, st, js -> [id, js[st.indexOf('normal')], js[st.indexOf('tumor')]] }
            | set { methbat_for_yaml_ch }
/*
        amber.out.for_yaml_summary
            | map { meta, amberQC, amberBAF -> [meta.id, meta, amberQC, amberBAF] }
            | join( purple.out.for_yaml_summary | map { meta, pur, dr, cnv -> [meta.id, pur, dr, cnv] } )
            | join( chord_hrd.out.for_yaml_summary | map { meta, f -> [meta.id, f] } )
            | join( scarhrd.out.for_yaml_summary | map { meta, f -> [meta.id, f] } )
            | join( pcgr_v212_deepSomatic.out.for_yaml_summary | map { meta, f -> [meta.id, f] } )
            | join( wakhan.out.wakhanTSV | map { meta, f -> [meta.id, f] } )
            | join( methbat_for_yaml_ch )
            | join( owl_for_yaml_ch )
            | join( cramino_for_yaml_ch )
            | map { id, meta, amberQC, amberBAF, pur, driver, cnv, chord, scar, pcgr, wak, mb_n, mb_t, owl_n, owl_t, cr_n, cr_t ->
                tuple(meta, [
                    amberQC: amberQC, amberBAF: amberBAF,
                    purple_purity: pur, purple_driver: driver, purple_cnv: cnv,
                    chord: chord, scarhrd: scar, pcgr: pcgr, wakhan: wak,
                    methbat_n: mb_n, methbat_t: mb_t, owl_n: owl_n, owl_t: owl_t,
                    cramino_n: cr_n, cramino_t: cr_t
                ])
            }
            | set { for_summary_final_ch }
*/
        amber.out.for_yaml_summary
            | map { meta, amberQC, amberBAF -> [meta.id, meta, amberQC, amberBAF] }
            | join( purple.out.for_yaml_summary | map { meta, pur, dr, cnv -> [meta.id, pur, dr, cnv] } )
            | join( scarhrd_purple.out.for_yaml_summary | map { meta, f -> [meta.id, f] } )
           // | join( scarhrd_wakhan.out.for_yaml_summary | map { meta, f -> [meta.id, f] } )
            | join( scarhrd_wakhan_bed.out.for_yaml_summary | map { meta, txt, json -> [meta.id, txt, json] } )
            | join( pcgr_v212_deepSomatic.out.for_yaml_summary | map { meta, f -> [meta.id, f] } )
            | join( wakhan.out.wakhanTSV | map { meta, f -> [meta.id, f] } )
            | join( methbat_for_yaml_ch )
            | join( owl_for_yaml_ch )
            | join( cramino_for_yaml_ch )
            | map { id, meta, amberQC, amberBAF, pur, driver, cnv, scar_purple,scarhrd_wakhan_txt,scarhrd_wakhan_json, pcgr, wak, mb_n, mb_t, owl_n, owl_t, cr_n, cr_t ->
                tuple(meta, [
                    amberQC: amberQC,
                    amberBAF: amberBAF,
                    purple_purity: pur,
                    purple_driver: driver,
                    purple_cnv: cnv,
                    scarhrd: scar_purple,
                    scarhrd_wakhan_txt:scarhrd_wakhan_txt,scarhrd_wakhan_json:scarhrd_wakhan_json,
                    pcgr: pcgr,
                    wakhan: wak,
                    methbat_n: mb_n,
                    methbat_t: mb_t,
                    owl_n: owl_n,
                    owl_t: owl_t,
                    cramino_n: cr_n,
                    cramino_t: cr_t
                ])
            }
            | set { for_summary_final_ch }


        purple_genome_view(for_summary_final_ch)

        for_summary_final_ch
            .join(purple_genome_view.out.png)
            .join(hrd_scores.out.hrdetect_json)
            .join(hrd_scores.out.chord_json)
            .join(hrd_scores.out.snv2020_json)
            .join(hrd_scores.out.snv2015_json)
            .join(hrd_scores.out.indel_json)
            .join(hrd_scores.out.dbs_json)
            | map { meta, data, cnv_plot, hrdetect, chord, snv2020,snv2015,indel,dbs->
                tuple(meta, data + [cnv_plot: cnv_plot, hrdetectJson: hrdetect, chordJson: chord,snv2020Json:snv2020,snv2015Json:snv2015,indelJson:indel,dbsJson:dbs])
            }
            | set { for_summary_final_with_cnv_plot_ch }

        collect_clinical_summary(for_summary_final_with_cnv_plot_ch)





    emit:
        phasedAll        = phasedAll
        somaticVcfPass   = deepSomatic.out.pcgr_vcf     // (meta, PASS vcf, idx) — tumour-only, for RNA ASE
        methylBedByCase  = methyl_bed_by_case           // (caseID, normalBed, tumorBed) — for INTEGRATION
        severusSV        = severus_edits.out.vcf         // (meta, somatic SV vcf) — for fusion x SV
        clinicalSummary  = collect_clinical_summary.out.json
}


/* =============================================================================
 *  DNA ARM (composed) — reusable by FULL and DNA entry points
 * ========================================================================== */
workflow DNA_ARM {
    take: per_sample_ch
    main:
        DNA_PREPROCESS(per_sample_ch)
        DNA_PREPHASE(DNA_PREPROCESS.out.aligned_tn_paired)
        DNA_PHASE(DNA_PREPHASE.out.hiphaseInput, DNA_PREPHASE.out.sawfish_supporting_reads)
        DNA_SOMATIC(DNA_PHASE.out.phasedAll)
    emit:
        phasedAll       = DNA_PHASE.out.phasedAll
        somaticVcfPass  = DNA_SOMATIC.out.somaticVcfPass
        methylBedByCase = DNA_SOMATIC.out.methylBedByCase
        severusSV       = DNA_SOMATIC.out.severusSV
}


/* =============================================================================
 *  RNA SUBWORKFLOWS  (Kinnex Iso-Seq)
 * ========================================================================== */

workflow RNA_PREPROCESS {
    take: rna_ch                         // (meta, bams)  bams = grouped list (>=1 movie)
    main:
        // one demux BAM per RNA sample is typical; merge only when >1 movie
        rna_ch
            | branch { meta, bams ->
                multi:  (bams instanceof List) && (bams.size() > 1)
                single: true
              }
            | set { rna_branched }

        merge_ubams(rna_branched.multi)

        rna_branched.single
            | map { meta, bams -> tuple(meta, (bams instanceof List) ? bams[0] : bams) }
            | mix(merge_ubams.out)
            | set { ubam_final }         // (meta, singleBam)

        isoseq_refine_cluster(ubam_final)

        // per-molecule (refined FLNC) alignment — the BAM haplotag + ASE consume
        pbmm2_align_refined_forIsocall(isoseq_refine_cluster.out.isoseq_bam_refined)
        isocallProfile(pbmm2_align_refined_forIsocall.out.bam)
        //isocallCall(isocallProfile.out.profile)

        // clustered alignment — feeds collapse / pigeon / fusion
        pbmm2_align_clust(isoseq_refine_cluster.out.isoseq_bam_clustered)

        isoseq_refine_cluster.out.isoseq_bam_refined
            .join(pbmm2_align_clust.out.bam)
            | map { meta, refinedBam, refinedPbi, clusteredBam, clusteredPbi -> tuple(meta, [refinedBam, refinedPbi, clusteredBam, clusteredPbi]) }
            | set { isoseq_pbmm2_joined }

        isoseq_collapse(isoseq_pbmm2_joined)

        // assemble the canonical preprocess map (keys consumed downstream)
        isoseq_refine_cluster.out.isoseq_bam_refined
            .join(isoseq_refine_cluster.out.refine_report_json)
            .join(isoseq_refine_cluster.out.isoseq_bam_clustered)
            .join(pbmm2_align_clust.out.bam)
            .join(isoseq_collapse.out.collapsed_gff)
            .join(pbmm2_align_refined_forIsocall.out.bam)
            | map { meta, rb, rp, refine_json, cb, cp, pb, pi, gff, counts, abund,pbbm2_r, pbbm2_i ->
                tuple(meta, [
                    refinedBAM:   rb, refinedPBI:   rp,
                    clusteredBAM: cb, clusteredPBI: cp,
                    pbmm2BAM:     pb, pbmm2BAI:     pi,
                    collapsedGFF: gff, flncCounts:  counts, flncAbundance: abund,
                    refinedPbmm2BAM: pbbm2_r, refinedPbmm2BAI: pbbm2_i,
                    refineReportJSON: refine_json
                ])
            }
            | set { preprocess_all_joined }
    emit:
        preprocessFullOutput = preprocess_all_joined
        // (meta, bam, bai) aligned refined FLNC BAM — per-molecule, for RNA_HAPLOTAG/ASE
        isoseqForSummary    = isoseq_refine_cluster.out.refine_report_json
        refinedAlignedBam    = pbmm2_align_refined_forIsocall.out.bam
}

workflow RNA_TRANSCRIPTOME {
    take: preprocessFullOutput
    main:
        pigeon_classify(preprocessFullOutput)
        sqanti3_QC(preprocessFullOutput)
        oarFish(preprocessFullOutput)

    emit:
        classification = pigeon_classify.out.classification   // filtered_lite classification
        sortedGFF      = pigeon_classify.out.sortedGFF        // for ISA importGTF (splicing)
        pigeon         = pigeon_classify.out.pigeon
        sqanti3        = sqanti3_QC.out.sqanti3QC
        pigeonForSummary    = pigeon_classify.out.pigeon_reports_json
}

workflow RNA_FUSION {
    take: preprocessFullOutput
    main:
        pbfusion(preprocessFullOutput)
    emit:
        fusion   = pbfusion.out.fusion
        fusionForSummary = pbfusion.out.inhouse_fusion
}

workflow RNA_SUMMARY {
    take:
        isoseqForSummary
        pigeonForSummary
        fusionForSummary
    main:
        isoseqForSummary.join(pigeonForSummary).join(fusionForSummary)
            | map { meta, refine_json, pigeonFiltered_json, pigeonRaw_json, fusions ->  
                tuple(meta, [
                    refineReportJSON: refine_json,
                    pigeonFilteredJSON: pigeonFiltered_json,
                    pigeonRawJSON: pigeonRaw_json,
                    fusionInhouse: fusions
                ])
            }
            | set { rna_summary_joined }
        collect_clinical_summaryRNA(rna_summary_joined)


}


/* RNA arm (composed) — reusable by FULL and RNA entry points */
workflow RNA_ARM {
    take: rna_ch
    main:
        RNA_PREPROCESS(rna_ch)
        RNA_TRANSCRIPTOME(RNA_PREPROCESS.out.preprocessFullOutput)
        RNA_FUSION(RNA_PREPROCESS.out.preprocessFullOutput)
        RNA_SUMMARY(
            RNA_PREPROCESS.out.isoseqForSummary,
            RNA_TRANSCRIPTOME.out.pigeonForSummary,
            RNA_FUSION.out.fusionForSummary
        )
    emit:
        preprocessFullOutput = RNA_PREPROCESS.out.preprocessFullOutput
        refinedAlignedBam    = RNA_PREPROCESS.out.refinedAlignedBam   // -> RNA_HAPLOTAG (next)
        classification       = RNA_TRANSCRIPTOME.out.classification
        sortedGFF            = RNA_TRANSCRIPTOME.out.sortedGFF
        fusion               = RNA_FUSION.out.fusion
}


/* =============================================================================
 *  CROSS-ARM SUBWORKFLOWS  (RNA consumes DNA germline/somatic; join key = caseID)
 * ========================================================================== */

workflow RNA_HAPLOTAG {
    take:
        rna_refined_bam       // (meta, bam, bai)     meta.id = caseID
        germline_by_case      // (caseID, phasedVcf)  hiPhase-phased NORMAL DeepVariant
    main:
        rna_refined_bam
            | map { meta, bam, bai -> [meta.id, meta, bam] }
            | join(germline_by_case)                       // [caseID, meta, bam, phasedVcf]
            | map { id, meta, bam, vcf -> tuple(meta, [bam: bam, phasedVcf: vcf]) }
            | set { haplotag_in }
        whatshap_haplotag(haplotag_in)
    emit:
        haplotaggedBam = whatshap_haplotag.out.bam            // (meta, bam, bai)
        haplotagList   = whatshap_haplotag.out.haplotag_list
}

workflow RNA_SOMATIC {
    take:
        rna_refined_bam       // (meta, bam, bai)   ASE runs on the pre-haplotag refined BAM
        germline_by_case      // (caseID, vcf)      NORMAL phased DeepVariant
        somatic_by_case       // (caseID, vcf)      DeepSomatic PASS (tumour-only)
        preprocessFullOutput  // (meta, dataMap)    -> flncAbundance
        classification        // (meta, filtered_lite_classification.txt)
        sortedGFF             // (meta, collapsed.sorted.gff)
    main:
        // ---- ASE: bam + germline (required) + somatic (attach if present) ----
        rna_refined_bam
            | map { meta, bam, bai -> [meta.id, meta, bam] }
            | join(germline_by_case)                          // [caseID, meta, bam, gVcf]
            | join(somatic_by_case, remainder: true)          // + sVcf (null if absent)
            | filter { it[1] != null }                        // drop somatic-only rows
            | map { id, meta, bam, gVcf, sVcf ->
                tuple(meta, [ bam: bam, germlineVcf: gVcf, somaticVcf: (sVcf ?: []) ])
            }
            | set { ase_in }
        ase_readcounter(ase_in)

        // abundance (isoseq collapse *.abundance.txt) keyed by caseID
        preprocessFullOutput
            | map { meta, d -> [meta.id, d.flncAbundance] }
            | set { abundance_by_case }

        // ---- expression outliers: classification + abundance ----
        classification
            | map { meta, cls -> [meta.id, meta, cls] }
            | join(abundance_by_case)                         // + abundance
            | map { id, meta, cls, ab -> tuple(meta, [classification: cls, abundance: ab]) }
            | set { expr_in }
        expression_outlier(expr_in)

        // ---- splicing / exon-skipping: gff + abundance + classification ----
        classification
            | map { meta, cls -> [meta.id, meta, cls] }
            | join( sortedGFF | map { meta, gff -> [meta.id, gff] } )
            | join(abundance_by_case)                         // + abundance
            | map { id, meta, cls, gff, ab -> tuple(meta, [classification: cls, gff: gff, abundance: ab]) }
            | set { splice_in }
        splicing_isoformswitch(splice_in)
    emit:
        germlineAse = ase_readcounter.out.germline_ase
        somaticAse  = ase_readcounter.out.somatic_ase
        geneCounts  = expression_outlier.out.geneCounts       // -> INTEGRATION (methylation x expression)
        outliers    = expression_outlier.out.outliers
        splicing    = splicing_isoformswitch.out.panel
        //asEvents    = splicing_isoformswitch.out.as_events    // carries NMD column (NMD x ASE, later)
}


/* =============================================================================
 *  INTEGRATION  (DNA <-> RNA cross-arm readouts; join key = caseID)
 * ========================================================================== */
workflow INTEGRATION {
    take:
        methyl_bed_by_case    // (caseID, normalBed, tumorBed)  pb-CpG-tools combined bedMethyl
        geneCounts            // (meta, geneCounts.tsv)         RNA gene FL counts
        fusion_out            // (meta, [fusion files])         pbfusion_v2 outputs
        severus_sv            // (meta, severus somatic SV vcf)
    main:
        // ---- methylation x expression ----
        geneCounts
            | map { meta, gc -> [meta.id, meta, gc] }
            | join(methyl_bed_by_case)                        // [caseID, meta, gc, nBed, tBed]
            | map { id, meta, gc, nBed, tBed ->
                tuple(meta, [geneCounts: gc, normalBed: nBed, tumorBed: tBed])
            }
            | set { methExpr_in }
        methylation_expression(methExpr_in)

        // ---- fusion x SV (extract the pbfusion breakpoints bed from the fusion outputs) ----
        fusion_out
            | map { meta, files ->
                def fl  = (files instanceof List) ? files : [files]
                def bed = fl.find { it.name.endsWith('.groups.bed') } ?: fl.find { it.name.contains('breakpoints') && it.name.endsWith('.bed') }
                [meta.id, meta, bed]
            }
            | filter { it[2] != null }
            | join( severus_sv | map { meta, vcf -> [meta.id, vcf] } )
            | map { id, meta, bed, vcf -> tuple(meta, [fusionBed: bed, sv: vcf]) }
            | set { fusSv_in }
        fusion_sv_concordance(fusSv_in)
    emit:
        methylExpr = methylation_expression.out.table
        fusionSV   = fusion_sv_concordance.out.table
}


/* =============================================================================
 *  BANNER
 * ========================================================================== */
log.info """\
======================================================
KG Vejle: PacBio LRS  DNA(somatic) + RNA(Kinnex)  v2
======================================================
Genome        : ${params.genome}
Genome FASTA  : ${params.genome_fasta}
Genome ver.   : ${params.genome_version}
ROI           : ${params.ROI}
GATK image    : ${params.gatk_image}
skipRNA       : ${params.skipRNA}
Samplesheet   : ${params.samplesheet}
RunID         : ${runID}
workDir       : ${workflow.workDir}
"""


/* =============================================================================
 *  ENTRY POINTS
 * ========================================================================== */

// default = FULL
workflow {
    DNA_ARM(per_sample_ch)

    if (!params.skipRNA) {
        RNA_ARM(rna_ch)

        // ---- cross-arm keyed channels (join key = caseID) ----
        DNA_ARM.out.phasedAll
            | map { meta, data -> [meta.id, data.dv_vcf] }        // NORMAL phased DeepVariant
            | set { germline_by_case }

        DNA_ARM.out.somaticVcfPass
            | map { meta, vcf, idx -> [meta.id, vcf] }            // DeepSomatic PASS (tumour-only)
            | set { somatic_by_case }

        RNA_HAPLOTAG(RNA_ARM.out.refinedAlignedBam, germline_by_case)
        

        
        if (params.somaticRNA) {
        RNA_SOMATIC(
            RNA_ARM.out.refinedAlignedBam,
            germline_by_case,
            somatic_by_case,
            RNA_ARM.out.preprocessFullOutput,
            RNA_ARM.out.classification,
            RNA_ARM.out.sortedGFF
        )

        // ---- cross-arm integration (methylation x expression, fusion x SV) ----
        INTEGRATION(
            DNA_ARM.out.methylBedByCase,
            RNA_SOMATIC.out.geneCounts,
            RNA_ARM.out.fusion,
            DNA_ARM.out.severusSV
        )
        }
    }
}

// DNA-only
workflow DNA {
    DNA_ARM(per_sample_ch)
}

// RNA-only
workflow RNA {
    if (params.skipRNA) { exit 1, "RNA entry called with --skipRNA set." }
    RNA_ARM(rna_ch)
}
