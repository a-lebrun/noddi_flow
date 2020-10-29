#!/usr/bin/env nextflow

params.input = false
params.load_kernels = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    engine = new groovy.text.SimpleTemplateEngine()

    bindings = ["noddi_nb_threads":"$params.noddi_nb_threads"]

    template = engine.createTemplate(usage.text).make(bindings)

    print template.toString()
    return
}

log.info "Nextflow NODDI pipeline"
log.info "======================="
log.info ""
log.info "Start time: $workflow.start"
log.info ""

log.debug "[Command-line]"
log.debug "$workflow.commandLine"
log.debug ""

log.info "[Git Info]"
log.info "$workflow.repository - $workflow.revision [$workflow.commitId]"
log.info ""

workflow.onComplete {
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    log.info "Execution duration: $workflow.duration"
}

if (params.input){
    log.info "Input: $params.input"
    input = file(params.input)
    in_data = Channel
        .fromFilePairs("$input/**/*{ad.nii.gz,brain_mask.nii.gz,bval,bvec,dwi.nii.gz,fa.nii.gz,md.nii.gz}",
                       size: 7,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}

(data_for_priors, data_for_noddi) = in_data
    .map{sid, ad, brain_mask, bvals, bvecs, dwi, fa, md -> [tuple(sid, ad, fa, md),
                                                            tuple(sid, dwi, bvals, bvecs, brain_mask)]}
    .separate(2)

process Compute_NODDI_Priors {
    cpus 1

    input:
    set sid, file(ad), file(fa), file(md) from data_for_priors

    output:
    file "${sid}__average_ad_1fiber.txt"
    file "${sid}__average_md_vent.txt"
    file "${sid}__average_ad_1fiber.json" into all_ad_to_collect
    file "${sid}__average_md_vent.json" into all_md_to_collect

    script:
    """
    scil_compute_NODDI_priors.py $fa $ad $md\
        --out_txt_1fiber ${sid}__average_ad_1fiber.txt\
        --out_txt_ventricles ${sid}__average_md_vent.txt

    jq '[.]' < ${sid}__average_ad_1fiber.txt >\
        ${sid}__average_ad_1fiber.json
    jq '[.]' < ${sid}__average_md_vent.txt >\
        ${sid}__average_md_vent.json

    """
}

all_ad_to_collect
    .collect()
    .set{all_ad_stats_for_mean_ad}

process Compute_Mean_AD {
    tag = {"Mean_AD"}
    publishDir = params.MeanStatsPublishDir

    input:
    file(all_ad_stats) from all_ad_stats_for_mean_ad

    output:
    file "mean_ad.json" into mean_ad

    script:
    """
    cat $all_ad_stats | jq -s add | jq add/length > mean_ad.json
    """
}

all_md_to_collect
    .collect()
    .set{all_md_stats_for_mean_md}

process Compute_Mean_MD {
    tag = {"Mean_MD"}
    publishDir = params.MeanStatsPublishDir

    input:
    file(all_md_stats) from all_md_stats_for_mean_md

    output:
    file "mean_md.json" into mean_md

    script:
    """
    cat $all_md_stats | jq -s add | jq add/length > mean_md.json
    """
}

process Compute_NODDI {
    cpus params.noddi_nb_threads

    input:
    set sid, file(dwi), file(bvals), file(bvecs), file(brain_mask) from data_for_noddi
    file(mean_ad) from mean_ad.first()
    file(mean_md) from mean_md

    output:
    file "${sid}__FIT_dir.nii.gz"
    file "${sid}__FIT_ISOVF.nii.gz"
    file "${sid}__FIT_ICVF.nii.gz"
    file "${sid}__FIT_OD.nii.gz"

    script:
    option_noddi=""
    if (params.load_kernels) {
      option_noddi="--load_kernels $params.load_kernels"
    }
    """
    mean_ad=\$(<$mean_ad)
    mean_md=\$(<$mean_md)

    scil_compute_NODDI.py $dwi $bvals $bvecs --mask $brain_mask\
        --para_diff \$mean_ad\
        --iso_diff \$mean_md\
        --processes $params.noddi_nb_threads \
        $option_noddi

    mv results/FIT_dir.nii.gz ${sid}__FIT_dir.nii.gz
    mv results/FIT_ICVF.nii.gz ${sid}__FIT_ICVF.nii.gz
    mv results/FIT_ISOVF.nii.gz ${sid}__FIT_ISOVF.nii.gz
    mv results/FIT_OD.nii.gz ${sid}__FIT_OD.nii.gz

    rm -rf results
    """
}
