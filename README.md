# vespa.net - Signaling reconstruction module for the VESPA R-package

## 1. Installation
This module requires the [Snakemake](https://snakemake.readthedocs.io/en/stable/) workflow manager and Docker to be installed and ready to use. It has only been tested on a CentOS 7 cluster, but it might work on other systems too.

### Download repository
Clone the repository to your working directory.

```
git clone git@github.com:califano-lab/vespa.net.git
```

### Singularity
If you are using [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) on a HPC environment, build Singularity images for ARACNe and VESPA first:

```
singularity build aracne.simg docker://ghcr.io/califano-lab/vespa.aracne:latest
singularity build vespa.simg docker://ghcr.io/califano-lab/vespa:latest
```

## 2. Using vespa.net
### Proteomic data
Add one or several RDS phosphopeptide abundance files (e.g. `CPTAC_S045_COAD_phospho.rds`) with different sample tags (e.g. `CPTAC_S045_COAD`) to the root directory. Add corresponding whole proteome files (`CPTAC_S045_COAD_proteo.rds`) with paired tags. Note, that they only differ by file ending (`_phospho.rds` / `_proteo.rds`). If no whole proteome files are available, make a duplicate copy of the `_phospho.rds` file and rename it to `_proteo.rds`, the algorithm will then skip this step.

### Proteomic reference data for optimization
The ``vespa.net`` module optimizes signalons using the ``metaVIPER`` algorithm for a target phosphoproteome. A single RDS phosphopeptide abundance file (i.e. `reference.rds`) needs to be added to the root directory. This file can be identical to one of the learning RDS phosphopeptide abundance files (e.g. `CPTAC_S045_COAD_phospho.rds`).

### FASTA library
Place the FASTA library (without decoys) used for the MS analysis in file ``library.fasta`` in the root directory.

## 3. Run analysis

Submit the Snakemake job. On a local computer, this could for example be:

```
snakemake --snakefile Snakefile -j 64 --restart-times 2

```

In a HPC environment with SLURM, the following command might be adapted:

```
sbatch --qos=1day --time=1-00:00:00 --mem-per-cpu=8192 snakemake --snakefile Snakefile --use-singularity -j 64 --restart-times 2 --cluster-config envs/res.json --cluster "sbatch --ntasks {cluster.nCPUs} --mem-per-cpu {cluster.memory} --qos=6hours --time=6:00:0"
```
