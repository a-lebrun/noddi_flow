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
                "lambda1": "$params.lambda1",
                "lambda2": "$params.lambda2",
                "output_dir":"$params.output_dir",
                "b_thr":"$params.b_thr",
                "run_priors_only":"$params.run_priors_only",
                "nb_subjects_for_priors":"$params.nb_subjects_for_priors",
                "fa_min":"$params.fa_min",
                "fa_max":"$params.fa_max",
                "md_min":"$params.md_min",
                "roi_radius":"$params.roi_radius"]

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
log.info "Lambda 1: $params.lambda1"
log.info "Lambda 2: $params.lambda2"
log.info "b-threshold: $params.b_thr"
log.info ""
log.info "[NODDI priors]"
log.info "Run priors only: $params.run_priors_only"
log.info "Nb subjects for priors: $params.nb_subjects_for_priors"
log.info "FA min: $params.fa_min"
log.info "FA max: $params.fa_max"
log.info "MD min: $params.md_min"
log.info "ROI radius: $params.roi_radius"
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

    in_data_priors = Channel
        .fromFilePairs("$input/**/*{ad.nii.gz,fa.nii.gz,md.nii.gz,rd.nii.gz}",
                       size: 4,
                       maxDepth:1,
                       flat: true) {it.parent.name}
}

(all_data_for_kernels, data_for_noddi) = in_data
    .map{sid, mask, bval, bvec, dwi ->
        [tuple(sid, mask, bval, bvec, dwi),
         tuple(sid, mask, bval, bvec, dwi)]}
    .separate(2)

in_data_priors.take(params.nb_subjects_for_priors).set{sub_in_data_priors}

process Compute_Priors {
  cpus 1

  input:
    set sid, file(ad), file(fa), file(rd), file(md) from sub_in_data_priors

  output:
    set val("Priors"), "${sid}__para_diff.txt", "${sid}__iso_diff.txt", "${sid}__mask_1fiber.nii.gz", "${sid}__mask_ventricles.nii.gz" into priors_for_mean

  when:
    params.run_priors_only

  script:
    """
    scil_NODDI_priors.py $fa $ad $rd $md\
      --fa_min $params.fa_min\
      --fa_max $params.fa_max\
      --md_min $params.md_min\
      --roi_radius $params.roi_radius\
      --out_txt_1fiber_para ${sid}__para_diff.txt\
      --out_txt_ventricles ${sid}__iso_diff.txt\
      --out_mask_1fiber ${sid}__mask_1fiber.nii.gz\
      --out_mask_ventricles ${sid}__mask_ventricles.nii.gz\
    """
}

priors_for_mean
    .groupTuple()
    .set{all_priors_for_mean}

process Mean_Priors {
  cpus 1
  publishDir = "${params.output_dir}/Mean_Priors"

  input:
    set sid, file(para_diff), file(iso_diff) from all_priors_for_mean

  output:
    file "para_diff.txt"
    file "iso_diff.txt"

  script:
    """
    cat $para_diff > all_para_diff.txt
    awk '{ total += \$1; count++ } END { print total/count }' all_para_diff.txt > para_diff.txt
    cat $iso_diff > all_iso_diff.txt
    awk '{ total += \$1; count++ } END { print total/count }' all_iso_diff.txt > iso_diff.txt
    """
}

all_data_for_kernels.first().set{unique_data_for_kernels}

process Compute_Kernel {
  cpus params.nb_threads
  publishDir = "${params.output_dir}/Compute_Kernel"

  input:
    set sid, file(brain_mask), file(bval), file(bvec), file(dwi) from unique_data_for_kernels

  output:
    file("kernels/") into kernel_for_noddi

  when:
    !params.run_priors_only

  script:
    """
    export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=
    export OMP_NUM_THREADS=
    export OPENBLAS_NUM_THREADS=1
    
    scil_NODDI_maps.py $dwi $bval $bvec --mask $brain_mask\
      --para_diff $params.para_diff\
      --iso_diff $params.iso_diff\
      --lambda1 $params.lambda1\
      --lambda2 $params.lambda2\
      --processes $params.nb_threads\
      --tolerance $params.b_thr\
      --save_kernels kernels/ \
      --compute_only
    """
}

data_for_noddi
  .combine(kernel_for_noddi)
  .set{data_with_kernel_for_noddi}

process Compute_NODDI {
    cpus params.nb_threads
    memory { params.memory_limit * task.attempt }

    input:
      set sid, file(brain_mask), file(bval), file(bvec), file(dwi), file(kernels) from data_with_kernel_for_noddi

    output:
      file "${sid}__fit_dir.nii.gz"
      file "${sid}__fit_FWF.nii.gz" // FreeWater partie ISO
      file "${sid}__fit_NDI.nii.gz" // Neurite density index - Intra cerebral volume fraction 
      file "${sid}__fit_ODI.nii.gz" // nufo like - number of orientation (float 0-1) - correlé à nufo (1,2,3,4,5)

    script:
      """
      export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=
      export OMP_NUM_THREADS=
      export OPENBLAS_NUM_THREADS=1
      
      scil_NODDI_maps.py $dwi $bval $bvec --mask $brain_mask\
        --para_diff $params.para_diff\
        --iso_diff $params.iso_diff\
        --processes $params.nb_threads \
        --load_kernels $kernels

      mv results/fit_dir.nii.gz ${sid}__fit_dir.nii.gz
      mv results/fit_FWF.nii.gz ${sid}__fit_FWF.nii.gz
      mv results/fit_NDI.nii.gz ${sid}__fit_NDI.nii.gz
      mv results/fit_ODI.nii.gz ${sid}__fit_ODI.nii.gz

      rm -rf results
      """
}
