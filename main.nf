#!/usr/bin/env nextflow

params.input = false
params.help = false

if(params.help) {
    usage = file("$baseDir/USAGE")

    engine = new groovy.text.SimpleTemplateEngine()

    bindings = ["nb_threads":"$params.nb_threads",
                "memory_limit":"$params.memory_limit",
                "para_diff": "$params.para_diff",
                "iso_diff": "$params.iso_diff",
                "output_dir":"$params.output_dir",
                "b_thr":"$params.b_thr"]

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

log.info "[Inputs]"
log.info "Input: $params.input"
log.info "Output directory: $params.output_dir"
log.info ""

log.info "Options"
log.info "======="
log.info ""
log.info "[NODDI fitting]"
log.info "Parallel diff: $params.para_diff"
log.info "Iso diff: $params.iso_diff"
log.info "b-threshold: $params.b_thr"
log.info ""

log.info "Number of processes per tasks"
log.info "============================="
log.info "NODDI fitting: $params.nb_threads"
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
        .fromFilePairs("$input/**/*{brain_mask.nii.gz,bval,bvec,dwi.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}

(all_data_for_kernels, data_for_noddi) = in_data
    .map{sid, mask, bval, bvec, dwi ->
        [tuple(sid, mask, bval, bvec, dwi),
         tuple(sid, mask, bval, bvec, dwi)]}
    .separate(2)

all_data_for_kernels.first().set{unique_data_for_kernels}

process Compute_Kernel {
  cpus 1
  publishDir = "${params.output_dir}/Compute_Kernel"

  input:
    set sid, file(brain_mask), file(bval), file(bvec), file(dwi) from unique_data_for_kernels

  output:
    file("kernels/") into kernel_for_noddi

  script:
    """
    scil_compute_NODDI.py $dwi $bval $bvec --mask $brain_mask\
      --para_diff $params.para_diff\
      --iso_diff $params.iso_diff\
      --processes $params.nb_threads\
      --b_thr $params.b_thr\
      --save_kernels kernels/ \
      --compute_only
    """
}

data_for_noddi
  .combine(kernel_for_noddi)
  .set{data_with_kernel_for_noddi}

process Compute_NODDI {
    cpus params.nb_threads
    memory params.memory_limit

    input:
      set sid, file(brain_mask), file(bval), file(bvec), file(dwi), file(kernels) from data_with_kernel_for_noddi

    output:
      file "${sid}__FIT_dir.nii.gz"
      file "${sid}__FIT_ISOVF.nii.gz"
      file "${sid}__FIT_ICVF.nii.gz"
      file "${sid}__FIT_OD.nii.gz"

    script:
      """
      scil_compute_NODDI.py $dwi $bval $bvec --mask $brain_mask\
        --para_diff $params.para_diff\
        --iso_diff $params.iso_diff\
        --processes $params.nb_threads \
        --load_kernels $kernels

      mv results/FIT_dir.nii.gz ${sid}__FIT_dir.nii.gz
      mv results/FIT_ICVF.nii.gz ${sid}__FIT_ICVF.nii.gz
      mv results/FIT_ISOVF.nii.gz ${sid}__FIT_ISOVF.nii.gz
      mv results/FIT_OD.nii.gz ${sid}__FIT_OD.nii.gz

      rm -rf results
      """
}
