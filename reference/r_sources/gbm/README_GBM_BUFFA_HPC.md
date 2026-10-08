# GBM Buffa 2010 signature derivation on HPC

## Folder layout on HPC

```text
/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/
  CODE/
    run_buffa2010_hpc.R
    environment.yml
    create_conda_env.sh
  SLURM/
    run_buffa2010_smoke_test.sbatch
    run_buffa2010_final.sbatch
  DATA/
    gencode.v36.annotation.gtf.gene.probemap
    NORMALIZED/
      IT_GBM_RNA_Deseq2_vst_Combat.csv
      FR_GBM_RNA_Deseq2_vst_Combat.csv
      TCGA_GBM_RNA_Deseq2_vst.csv
  RESULTS/
  logs/
  envs/
```

## 1. Copy scripts to HPC

From your local machine, copy the generated folders to the HPC project folder. Example:

```bash
scp -r CODE SLURM README MoscatoG@<cluster>:/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/
```

## 2. Create the conda environment

On the HPC login node:

```bash
cd /home/MoscatoG/GBM_HYPOXIA_SIGNATURE
bash CODE/create_conda_env.sh
```

This creates:

```text
/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/envs/gbm_buffa_r
```

## 3. Run a smoke test

```bash
cd /home/MoscatoG/GBM_HYPOXIA_SIGNATURE
sbatch SLURM/run_buffa2010_smoke_test.sbatch
```

Monitor:

```bash
squeue -u $USER
sacct -j <JOBID> --format=JobID,JobName,State,Elapsed,MaxRSS,AllocCPUS
less logs/gbm_buffa_test_<JOBID>.out
less logs/gbm_buffa_test_<JOBID>.err
```

## 4. Run the final analysis

Default final parameters are:

```text
B_BOOT = 1000
N_MC   = 5000
CPUs   = 32
Memory = 128G
Time   = 48h
```

Submit:

```bash
cd /home/MoscatoG/GBM_HYPOXIA_SIGNATURE
sbatch SLURM/run_buffa2010_final.sbatch
```

Override parameters at submission time if needed:

```bash
sbatch --export=ALL,B_BOOT=500,N_MC=2000 SLURM/run_buffa2010_final.sbatch
```

## 5. Resume behavior

The final script runs with `--resume TRUE`. It writes checkpoint chunks under:

```text
RESULTS/BUFFA2010_EXACT_HPC/latest/run_path.txt
RESULTS/BUFFA2010_EXACT_HPC/run_<timestamp>/checkpoints/
```

If the job stops, resubmit the same sbatch. Existing bootstrap/Monte Carlo chunks will be skipped.

## 6. Main outputs

Inside the run directory:

```text
csv/Buffa2010_exact_ITxFR_full_ranking.csv
csv/GBM_Hypoxia_Buffa2010_exact_ITxFR_signature.csv
csv/GBM_Hypoxia_Buffa2010_exact_ITxFR_signature_top51.csv
csv/Buffa2010_exact_ITxFR_HS_all_signature_IT_FR_TCGA.csv
csv/Buffa2010_exact_ITxFR_HS_top51_IT_FR_TCGA.csv
rds/Buffa2010_exact_ITxFR_full_ranking.rds
rds/GBM_Hypoxia_Buffa2010_exact_ITxFR_signature.rds
rds/sessionInfo.rds
logs/run_buffa2010_hpc.log
plots/*.pdf
```

## 7. Copy results back locally

```bash
scp -r MoscatoG@<cluster>:/home/MoscatoG/GBM_HYPOXIA_SIGNATURE/RESULTS/BUFFA2010_EXACT_HPC/latest .
```

Note: `latest/run_path.txt` points to the real timestamped run folder.
