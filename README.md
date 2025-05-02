Nextflow NODDI pipeline
=======================

Run the NODDI priors estimation and NODDI modelling pipeline.

If you use this pipeline, please cite:

```
Kurtzer GM, Sochat V, Bauer MW (2017)
Singularity: Scientific containers for mobility of compute.
PLoS ONE 12(5): e0177459. 

P. Di Tommaso, et al. (2017) Nextflow enables reproducible computational workflows.
Nature Biotechnology 35, 316–319

H. Zhang, H., et al. (2012). NODDI: Practical in vivo neurite orientation
dispersion and density imaging of the human brain. NeuroImage, 61(4), 1000–1016. 

A. Daducci, et al. (2015). Accelerated Microstructure Imaging via Convex Optimization (AMICO)
from diffusion MRI data. NeuroImage, 105, 32–44. 
```

Requirements
------------

- [Nextflow](https://www.nextflow.io)
- Scilpy (https://github.com/scilus/scilpy)
- AMICO (https://github.com/daducci/AMICO)

Singularity/Docker
-----------
If you are on Linux, we recommend using the Singularity to run tractometry_flow pipeline.
If you have Apptainer (Singularity), launch your Nextflow command with:
`-with-singularity ABSOLUTE_PATH/scilus-2.1.0.sif`

Image is available [here](http://scil.dinf.usherbrooke.ca/en/containers_list/scilus-2.1.0.sif)

If you are on MacOS or Windows, we recommend using the Docker container to run tractometry_flow pipeline.
Launch your Nextflow command with:
`-with-docker scilus/scilus:2.1.0`

:warning: WARNING :warning:
---------
The official release 2.1.0 is **NOT** available now.

Please, either build the singularity container using this command:

`singularity build scilus_latest.sif docker://scilus/scilus:latest` 

and then launch your Nextflow command with:
`-with-singularity ABSOLUTE_PATH/scilus_latest.sif`

Or launch your Nextflow command with docker:
`-with-docker scilus/scilus:latest`

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`