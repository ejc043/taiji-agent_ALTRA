# taiji-run-altra

Faithfully re-implements the Taiji input construction from `Taiji_ALTRA/scripts/prepare_input_yaml.R`. Builds a `taiji_input.xlsx` (RNA + ATAC + HiC per cluster), generates per-sample TSVs/configs, and submits a SLURM array using the Taiji binary.

## Three modalities are required

| Modality | Tag | Source | Why |
|---|---|---|---|
| RNA-seq | `GeneQuant` | `{proj}_C{i}_rna.tsv` from pseudobulk-altra | TF expression |
| ATAC-seq | `PairedEnd` | `{proj}_C{i}_atac.bed.gz` from pseudobulk-altra | Taiji calls peaks internally via macs2 |
| HiC | `ChromosomeLoop` | `epitensor/hg38/epitensor_loop_top10p_87311.txt` (shared across all clusters) | TF→target edge evidence |

**Do not pre-call peaks.** The `PairedEnd` tag tells Taiji to run macs2 internally. Providing pre-called narrowPeak files requires a different tag and workflow.

## Inputs

| File | Description |
|---|---|
| `pseudobulk/` dir | Contains `{proj}_C{i}_rna.tsv` and `{proj}_C{i}_atac.bed.gz` for eligible clusters |
| `epitensor/hg38/epitensor_loop_top10p_87311.txt` | EpiTensor predicted chromatin loops (hg38 PBMC) |
| `dependencies_data/hg38/genome.fa` | hg38 genome FASTA |
| `dependencies_data/hg38/genes.gtf` | GENCODE v45 GTF |
| `cisbp_database/cisbp_human_2.meme` | CIS-BP motif file |
| `samples.csv` | `group,cohort` — one row per eligible cluster |
| `taiji_config.template.yml` | Three-placeholder template (input, output_dir, genome) |

## samples.csv format

```csv
group,cohort
GSM8554115_C2,CD14_Mono
GSM8554115_C3,NK
...
```

`group` = cluster name matching pseudobulk file prefix. `cohort` = dominant cell type from cluster_identity.csv (spaces replaced with underscores).

## taiji_config.template.yml

```yaml
input: "[insert_input_filepath_here]"
output_dir: "[insert_output_directory_here]"
tmp_dir: "${REPO_ROOT}/data/<dataset>/taiji_run/tmp"
picard: "/path/to/picard.jar"
genome: "[insert_genome_filepath_here]"
annotation: "${REPO_ROOT}/dependencies_data/hg38/genes.gtf"
motif_file: "${REPO_ROOT}/cisbp_database/cisbp_human_2.meme"
resource:
    ATAC_Align: { memory: 40 }
    RNA_Align: { memory: 40 }
    RNA_Make_Expr_Table: { memory: 25 }
    ATAC_Find_TFBS_Union: { memory: 10 }
    ATAC_Get_TFBS: { memory: 10 }
```

The three `[insert_*]` placeholders are filled per-sample by `taiji-runner`.

## SLURM array script requirements

The generated `taiji_array.sbatch` must have these tools in PATH **before** the Taiji binary runs:

```bash
# Required: macs2 for ATAC peak calling
export PATH="/path/to/macs2/bin:${PATH}"
# Required: bedGraphToBigWig for ATAC track visualization
export PATH="/path/to/ucsc-tools/bin:${PATH}"
```

On this HPC:
- `macs2`: `/stg3/data1/eunice/bin/miniconda3/envs/circlehunter/bin/macs2`
- `bedGraphToBigWig`: `/stg3/data1/eunice/bin/bedGraphToBigWig`

## Taiji binary

Use a binary compatible with the HPC OS. The CentOS binary (`taiji-CentOS-x86_64`) requires `libtinfo.so.5` which is absent on Debian/Ubuntu nodes — use the locally compiled Taiji 1.2.0 binary instead:

```
/stg3/data1/eunice/bin/Taiji.1.2.0/taiji
```

## SLURM resources per sample

```
--cpus-per-task=4
--mem=30G
--time=24:00:00
```

Genome index build (~3 GB) takes ~10 min on first run per sample. Subsequent stages complete in 10–40 min depending on cluster size.

## Expected outputs per sample

```
Output/Partial/<group>_output/
├── GeneRanks.tsv           ← per-TF PageRank scores (primary output)
├── GeneRanks_PValues.tsv
├── sciflow.db              ← workflow state (checkpoint; delete to force full rerun)
├── RNASeq/
│   ├── expression_profile.tsv
│   └── RNA-<group>_rep1_gene_quant.tsv
├── ATACSeq/
│   ├── Bed/                ← merged fragment BED
│   ├── GeneQuant/          ← gene accessibility scores
│   └── openChromatin.bed.gz ← called peaks
└── GENOME/
    └── genome.index        ← bowtie2 index (~3 GB)
```

`GeneRanks.tsv` has ~1100 rows (one per TF in cisbp_human_2.meme).

## Adaptation for a new dataset

1. Replace `samples.csv` with your cluster list (group = pseudobulk prefix, cohort = cell type)
2. Replace `epitensor/hg38/epitensor_loop_top10p_87311.txt` with the appropriate genome's EpiTensor file:
   - hg19: check `Annotation/Epitensor/hg19/`
   - mm10: use `Annotation/Epitensor/mm10/loops.txt`
3. Update `dependencies_data/<genome>/` paths in the config template
4. Adjust `--mem` if clusters are very large (>5000 cells → consider 50G)

## Reference scripts

- `data/ALTRA/build_taiji_input_p1.py` — builds the xlsx
- `data/ALTRA/run_taiji_p1.sh` — orchestration wrapper
- `data/ALTRA/taiji_run/taiji_array.sbatch` — SLURM array (edit PATH before reuse)
- `data/ALTRA/taiji_run/code/samples.csv` — P1 example
- `data/ALTRA/taiji_run/code/taiji_config.template.yml` — P1 example
