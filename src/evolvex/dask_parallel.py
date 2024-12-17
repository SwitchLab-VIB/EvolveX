import math

from dask.distributed import Client, LocalCluster, wait
from dask_jobqueue import SLURMCluster

def setup_dask_parallel_executor(GLOBALS):
    compute_env = GLOBALS.compute_env
    n_cores = GLOBALS.n_cores

    if compute_env == 'local':
        cluster = LocalCluster(
            n_workers = n_cores,
            threads_per_worker = 1,
            processes = True,
        )

    elif compute_env == 'SLURM':
        account = GLOBALS.account_name
        clusters = GLOBALS.cluster_name
        SLURM_job_prologue = GLOBALS.SLURM_job_prologue
        max_SLURM_jobs = GLOBALS.max_SLURM_jobs
        partition = GLOBALS.cluster_partition
        walltime = GLOBALS.walltime
        
        # SLURMCluster determines which resources are mobilized by each worker, which are then scaled according to the number of CPUs by the cluster.adapt function.
        SLURM_info_files_dir = GLOBALS.working_dir / 'slurm_info_files'; SLURM_info_files_dir.mkdir()

        n_cores_per_job = math.ceil(n_cores / max_SLURM_jobs)
        print(f'{max_SLURM_jobs} jobs, each running {n_cores_per_job} single core workers.', flush=True)

        cluster = SLURMCluster(
            processes = n_cores_per_job,
            cores = n_cores_per_job,
            memory = f'{3 * n_cores_per_job} GB', # Could turn this into a parameter for the config file
            job_extra_directives = [f'--{account=}', f'--{clusters=}', f'--{partition=}', f'--output={str(SLURM_info_files_dir)}/slurm-%j.out'],
            walltime = walltime,
            job_script_prologue = SLURM_job_prologue, # Used to setup each job and associated workers with the necessary modules and virtual environment
            death_timeout=300
        )
        cluster.scale(jobs = max_SLURM_jobs)
        #cluster.adapt(minimum_jobs = 1, maximum_jobs = max_SLURM_jobs, interval = '4s') # NOTE: Adaptive scaling does not work, workers randomly die which kills the whole run. Setting minimum_jobs = max_SLURM_jobs does not solve the problem.

    else:
        raise ValueError("compute_env must be either 'local' or 'SLURM'.")

    parallel_executor = Client(cluster) 

    print(f'Dask dashboard link: ', parallel_executor.dashboard_link, flush=True) # To monitor how well the parallelization is going

    return parallel_executor

def wait_and_remove(parallel_executor, futures):
    wait(futures)
    if futures:
        parallel_executor.cancel(futures)
    return