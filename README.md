# ChIPseq Workflow

This workflow performs ChIPseq analysis using [BWA-MEM](http://bio-bwa.sourceforge.net/), [MACS2](https://github.com/taoliu/MACS), and [deepTools](https://deeptools.readthedocs.io/en/develop/)

## Authors
* Dean Pettinga (@deanpettinga), https://github.com/deanpettinga

## Usage
**NOTE** this workflow is optimized for HPC3 @ Van Andel Institute.

### Step 1: Installation

make sure you are running [![Snakemake](https://img.shields.io/badge/snakemake-≥5.4.4-green.svg)](https://snakemake.bitbucket.io)

The following recipe provides established best practices for running and extending this workflow in a reproducible way.

1. [Fork](https://help.github.com/en/articles/fork-a-repo) the repo to a project directory on /secondary
2. [Clone](https://help.github.com/en/articles/cloning-a-repository) the fork to the desired working directory for your project.
3. [Create a new branch](https://git-scm.com/docs/gittutorial#_managing_branches) (the project-branch) within the clone and switch to it. The branch will contain any project-specific modifications (e.g. to configuration, but also to code).

### Step 2: Configure the workflow
* Modify the config, and any necessary sheets, e.g.:
  * src/samples.txt
  * src/config.yaml
  * src/cluster.yaml
* Move your sequencing reads to `raw_reads/`
