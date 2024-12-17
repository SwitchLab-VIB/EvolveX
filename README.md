# EvolveX

This repository contains the code of EvolveX, a *de novo* antibody computational design pipeline introduced in [link to paper](). Specifically, it corresponds to the computational pipeline that generates antibody designs given an initial set of antibody-antigen docks and a set of positions to mutate and explore, which is driven by the FoldX force field to optimize the binding affinity while maintaining the thermodynamic stability of the designed antibodies.

# Installation

- **Python version >= 3.9**
- Strongly recommended to create a virtual environment (run "python -m venv venv", then "source venv/bin/activate")
- Download and unzip the github repository.
- Run "pip install ." from the directory that contains the "setup.py file". This will create an "evolvex" command that you can run from the command line.
- Running EvolveX additionally requires FoldX, which can be obtained [here](https://foldxsuite.crg.eu/licensing-and-services).

The code has been tested on Linux and MacOS operating systems.

# How to run EvolveX

The evolvex command only takes a single input, which is a YAML configuration file, for which an example can be found in "evolvex_config_example.yaml". 

**If you wish to run EvolveX on the human Vsig4 example** showcased in our publication, simply set the number of CPU cores and the path to your FoldX folder in the pre-filled configuration file "evolvex_config_Vsig4.yaml" and run "evolvex evolvex_config_Vsig4.yaml". It will use the pre-generated antibody-antigen docks and input files in the "Vsig4_example" folder. Note that the search parameters have been set to a reduced version as it only runs 10 iterations, 2 models per dock and performs recombination every 5 iterations, which should take 1-2 hours to run on a personal laptop with 10 CPU cores. To run the same search as we did using 500 iterations, a population of 50 models per dock and recombination every 50 iterations you would need to run it on a lab cluster or HPC (see the [additional details section](#additional-details) to run it on SLURM-based HPCs), as it would take weeks to run on a personal laptop.

For details about the configuration file and additional input files needed to run EvolveX, read below.

## Additional details
The YAML configuration parameters are the following:

- The "antibody_chains" and "antigen_chains" parameters are self-explanatory. The antibody can be a single chain (i.e nanobody) or a standard double chain Fv.

- In the "Required paths" section:
  - The working_dir is where evolvex will write all the files it generates. If you are running on an HPC, make sure this folder points to a directory with enough space and which can perform fast writes.
  - The foldx_path should contain an executable file named "foldx" and a "rotabase.txt" file.
  - The Backbones_dir should contain PDB files corresponding to the initial set of antibody-antigen docks.
  - The PositionsToExplore_file_path should be a TSV-formated file containing the positions to mutate for each PDB file in Backbones_dir (see PositionsToExplore_example.tsv for format). "AA_Allowed" can either be a string of pre-selected amino acids in single letter format (e.g KRHDE), or set to "AUTO" to let evolvex test each individual mutation and determine which ones are worth trying during the MC-GA search. "MakeAla" can be set to "Y", in which case that position will be mutated to Alanine before the MC-GA search, or to "N", in which case the wildtype amino acid will be kept as the starting amino acid at that position before the MC-GA search.

- In the "Search algorithm settings" section, the population_size corresponds to the number of models that will be generated and explored PER DOCK. So if you have 100 PDBs and set this to 100, 10000 models will be generated and explored over the number of iterations you have selected. By default, we run 500 iterations, 50 models per dock and do a recombination step every 50 iterations.

- In the "Compute settings" section:
  - The "compute_env" can be set to "local" (default) or "SLURM" if running on a SLURM-based HPC.
  - The "n_cores" sets the number of parallel CPU cores to use, both for the local of SLURM compute environments.
  - When running on SLURM, additional parameters are required:
    - The number of CPU cores are split across "max_SLURM_jobs". Most SLURM HPCs limit the number of jobs a user can have in the queue at any time, so you should set the "max_SLURM_jobs" accordingly. For example, if n_cores=250 and max_SLURM_jobs=25, then 25 jobs with 10 CPU cores each will be submitted.
    - The "walltime" corresponds to the maximum time the jobs will run for.
    - Set the "account_name", "cluster_name" and "cluster_partition" according to your HPC.
    - Adapt the "SLURM_job_prologue" to your HPC, making sure the Python version you load is the same as the one used to create the virtual environment used to install Evolvex.

Once the YAML file is ready, run "evolvex evolvex_config.yaml". While running, you will find a "generated_models_info.csv" file in the working_dir, which contains information about each of the models generated at every iteration for every initial PDB dock. This is the file that should be used to then select the antibody designs. The corresponding PDB file of each model is saved in the "model_PDB_files" directory. 

**NOTE: When running on a SLURM environment, if a "tcp connection error" or any other similar error that suggests that the head process has lost communication with the workers arises when launching EvolveX from a login node, then try launching the evolvex command through a SLURM script so that the head process runs on a compute node instead (see the "evolvex_slurm_head_example.sbatch").**