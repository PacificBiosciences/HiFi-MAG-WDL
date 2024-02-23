Either `miniwdl` or `Cromwell` can be used to run workflows on the HPC.

# Installing and configuring `miniwdl`

## Requirements

- [`miniwdl`](https://github.com/chanzuckerberg/miniwdl) >= 1.9.0
- [`miniwdl-slurm`](https://github.com/miniwdl-ext/miniwdl-slurm)

## Configuring

An [example miniwdl.cfg file](miniwdl.cfg) is provided here. This should be placed at `~/.config/miniwdl.cfg` and edited to match your slurm configuration. This allows running workflows using a basic SLURM setup.

# Installing and configuring `Cromwell`

Cromwell supports a number of different HPC backends; see [Cromwell's documentation](https://cromwell.readthedocs.io/en/stable/backends/HPC/) for more information on configuring each of the backends.

# Configuring and running the workflow

## Filling out workflow inputs

Fill out any information missing in [the inputs file](inputs.hpc.json). Once you have downloaded the reference data bundle, ensure that you have replaced the `<local_path_prefix>` in the input template file with the local path to the reference datasets on your HPC.

See [the inputs section of the main README](../../README.md#workflow-inputs) for more information on the structure of the inputs.json file.

## Running the workflow

### Running via miniwdl

`miniwdl run workflows/main.wdl -i <inputs_json_file>`

### Running via Cromwell

`cromwell run workflows/main.wdl -i <inputs_json_file>`

# Reference data

The pipeline requires two publicly licensed reference databases, [UniRef100](https://www.uniprot.org/help/uniref) and the [Genome Taxonomy Database Toolkit (GTDB-Tk)](https://ecogenomics.github.io/GTDBTk/). For convenience we've provided the links below. The GTDB-tk DB can stay bundled as a `*.tar.gz` file, but UniRef100 needs to be unzipped.

```bash

# make a directory for and download the reference databases
mkdir dataset; cd dataset
wget https://data.gtdb.ecogenomic.org/releases/release207/207.0/auxillary_files/gtdbtk_r207_v2_data.tar.gz
wget https://zenodo.org/records/4626519/files/uniref100.KO.v1.dmnd.gz
gunzip uniref100.KO.v1.dmnd.gz

```