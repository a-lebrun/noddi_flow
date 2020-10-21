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

Usage
-----

See *USAGE* or run `nextflow run main.nf --help`