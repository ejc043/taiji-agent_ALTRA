# Environment management

This repo uses conda + lockfiles + (optionally) `conda-pack` to ship a
reproducible, "Docker-like" environment to SLURM nodes that don't allow
container runtimes. This doc walks through:

1. **Profiles** — install only what your dataset needs (the headline change as of 2026-04-26).
2. The one-command install (Tier 1).
3. Locking the env for byte-exact reproducibility (Tier 2).
4. Packing the env into a portable tarball (Tier 3).

## Profiles

Not every dataset needs the full R + Seurat/Signac stack. The env is split
into composable profiles so users only install what they'll use:

| Profile | Contents | Enables | Disk | Time |
|---------|----------|---------|------|------|
| `base` (default) | Python + xlsx I/O + MACS3 + local taiji-agent package | 5 of 6 skills (everything except pseudobulk-construct) | ~500 MB | ~5 min |
| `sc` (additive) | R-base + r-seurat + bioconductor-signac + supporting Bioconductor + r-remotes + GitHub R package (MuDataSeurat for `.h5mu`) | pseudobulk-construct | +3-4 GB | +15-30 min |
| `dev` (orthogonal) | pytest + pytest-cov + ruff + mypy + ipython | author tooling for editing the skills themselves | +500 MB | +3 min |
| `full` | base + sc + dev | all skills + dev tooling | ~5 GB | ~25 min |

**Profiles compose additively**: passing `--profile sc` always installs `base`
first, then layers SC on top. A user who installed `base` and later wants SC
can run `bash bin/install.sh --profile sc` and only the new R packages get
added — no rebuild.

### Per-skill membership

Each skill's `dependencies.yml` declares which profile it belongs to:

```yaml
# skills/detect-dataset-type/dependencies.yml
profile: base

# skills/pseudobulk-construct/dependencies.yml
profile: sc
```

`bin/doctor.sh --profile <name>` filters its verification table to skills
that actually need the requested profile, so you don't see false-failures
for SC skills when you only installed base.

### Picking a profile

| If your dataset is... | Use profile |
|------------------------|-------------|
| Bulk RNA-seq / ATAC-seq / HiC (TSV / narrowPeak / bedpe) | `base` |
| Single-cell (.rds / .h5ad / .h5mu) requiring pseudobulk | `sc` |
| Mixed cohort with both bulk + SC samples | `sc` (covers both branches) |
| Hacking on the skills themselves | `dev` (or `full` if also processing data) |
| Unsure | `base` — upgrade with `--profile sc` later if needed |

### Auto-detected install

`bin/auto-install.sh` runs `detect-dataset-type` on your data directory and
picks the profile for you:

```bash
bash bin/auto-install.sh --data-dir data/<dataset> --system macos
# → classifies bulk → installs --profile base
# → classifies single-cell → installs --profile sc
# → mixed → installs --profile sc

bash bin/auto-install.sh --data-dir data/<dataset> --system macos --dry-run
# → prints decision without installing
```

Useful when onboarding a new collaborator who's not yet sure what's in their
dataset.

## Tier 1 — One-command install

The end-to-end install:

```bash
bash bin/install.sh                          # base profile, default
```

That single command:

1. Creates the conda env `taiji-agent` from `environment.base.yml` using
   `micromamba` (5-10× faster than `conda`; falls back to `mamba`/`conda`
   if you pass `--solver`).
2. If profile includes `sc`: runs `bin/postinstall.R` inside the env to
   install MuDataSeurat from GitHub (it isn't on bioconda; needed for
   `.h5mu` input). SeuratDisk (formerly auto-installed for `.h5ad` input)
   was removed — install manually if you need it.
3. Runs `bin/install-taiji.sh` to download the Taiji binary for the
   current OS into `binaries/taiji`.

Common variants:

```bash
bash bin/install.sh --system centos                    # base + CentOS Taiji binary
bash bin/install.sh --profile sc                       # base + sc (R stack)
bash bin/install.sh --profile full                     # base + sc + dev
bash bin/install.sh --profile dev                      # base + dev
bash bin/install.sh --skip-taiji                       # env-only, no Taiji binary
bash bin/install.sh --skip-r                           # force-skip postinstall.R
bash bin/install.sh --solver mamba                     # if micromamba isn't installed
bash bin/install.sh --use-lockfile                     # install full env from conda-lock
```

After the install finishes:

```bash
micromamba activate taiji-agent
bash bin/doctor.sh
```

`doctor.sh` walks every `skills/*/dependencies.yml`, runs `which` for binaries
and `Rscript -e 'library(X)'` for R packages, and prints a per-skill PASS/FAIL
table. It exits 0 only if every required dep for every skill is present.
Useful flags: `--skill <name>`, `--no-r`, `--quiet`, `--json`.

### Installing micromamba (one-time, no root)

```bash
"${SHELL}" <(curl -L micro.mamba.pm/install.sh)
# follow prompts, then restart your shell
```

That drops a static binary into `~/.local/bin/micromamba`. No conda
installation, no system packages, no root.

## Tier 2 — Lock the env for byte-exact reproducibility

`environment.yml` declares packages with version *constraints* (e.g.
`r-seurat>=5.0`). The conda solver picks specific versions that satisfy those
constraints at install time. To freeze the result so the exact same versions
land on every machine:

```bash
# One-time install of conda-lock (drop into your base env or any env)
micromamba install -n base -c conda-forge conda-lock

# Generate the lockfile (commit the result to the repo)
conda-lock --file environment.yml --platform linux-64
```

This produces `conda-lock.linux-64.yml` with every transitive package +
version + hash + URL pinned. Subsequent installs that use the lockfile
will produce literally identical envs:

```bash
bash bin/install.sh --use-lockfile
```

When `environment.yml` changes, regenerate the lockfile:

```bash
conda-lock --file environment.yml --platform linux-64 --update <package>
# or wholesale:
conda-lock --file environment.yml --platform linux-64
git add conda-lock.linux-64.yml && git commit -m "lockfile: refresh"
```

For multi-platform support (e.g. linux-64 for SLURM + osx-arm64 for local
dev), pass `--platform` multiple times — `conda-lock` will produce one
lockfile per platform.

## Tier 3 — Pre-baked tarball (the closest thing to a Docker image without containers)

Once `taiji-agent` is built and verified on a machine, you can pack the env
into a portable tarball that other compatible Linux nodes can untar and use
in under 2 minutes — no network, no resolver, no compilation.

### Build the tarball

```bash
# Install conda-pack into the env (one-time)
micromamba install -n taiji-agent -c conda-forge conda-pack

# Pack the env. ~3-4 GB compressed for the full Python + R + Seurat stack.
micromamba run -n taiji-agent conda-pack -n taiji-agent \
  -o taiji-agent-env.tar.gz
```

Stage the tarball wherever your lab keeps shared artifacts, e.g.
`/stg3/data1/eunice/envs/`.

### Restore on a new node

```bash
mkdir -p ~/envs/taiji-agent
tar -xzf /stg3/data1/eunice/envs/taiji-agent-env.tar.gz -C ~/envs/taiji-agent

# Activate
source ~/envs/taiji-agent/bin/activate

# Rewrite the internal absolute paths so binaries find their libs.
# This is the conda-pack "unpack" step. Run it ONCE per untar.
~/envs/taiji-agent/bin/conda-unpack

# Verify
bash <repo>/bin/doctor.sh
```

That's it. The Taiji binary stays separate (it's one executable, not worth
packing) and lives in `binaries/taiji`.

### When to rebuild the tarball

- After any `environment.yml` or lockfile change.
- After any `bin/postinstall.R` change (e.g. updating MuDataSeurat).
- After a new R package gets added to the env.

The build is the slow step (~10-20 min); the restore is the fast one (<2 min),
which is what makes this approach scale across many nodes.

## Comparison vs. the alternatives

| Approach           | Root needed?  | Install time on target node | Reproducibility       | Works on SLURM? |
|--------------------|---------------|-----------------------------|-----------------------|-----------------|
| Docker             | yes (or rootless setup) | <1 min            | image hash            | usually no      |
| Singularity .sif   | yes to build  | <1 min                      | image hash            | depends on site |
| Plain conda        | no            | 10-30 min (network solve)   | loose                 | yes             |
| micromamba + lockfile | no         | 5-10 min (network)          | byte-pinned           | yes             |
| **conda-pack tarball** | **no**     | **<2 min (no network)**     | **byte-identical**    | **yes**         |

For your situation (UCSD SLURM, no Singularity/Docker), the conda-pack
tarball gives you Docker-image-equivalent immutability and portability
without any container runtime.

## Limitations

- **conda-pack is platform-specific.** A linux-64 tarball won't run on
  macOS. If you need both, build two tarballs from two source machines.
- **GitHub-only R packages are baked in at install time.** Updating
  MuDataSeurat means rebuilding the tarball. That's the right tradeoff
  for reproducibility.
- **Reference data (FASTA, GTF, MEME files) does NOT go in the tarball.**
  Genome files are too big and shared across many projects — they belong
  in `/stg3/data1/eunice/database/` (or your lab's equivalent) with
  `--genome` paths in the per-skill configs pointing at them.
- **Tarball size: 3-4 GB compressed** for Python + R + Seurat + Signac +
  Bioconductor deps. Manageable, but worth knowing before you try to email it.

## Per-skill dependency declarations

Each skill ships a `skills/<name>/dependencies.yml` that declares its
runtime needs (binaries, Python packages, R packages, optional data
files). These are documentation-as-code: they're what `bin/doctor.sh`
reads to verify the install, and they're easy to scan when extracting
a single skill from this repo.

Schema:

```yaml
skill: pseudobulk-construct
purpose: One-line description.

binaries:
  - {name: Rscript,  min_version: "4.2",  purpose: "..."}
  - {name: macs3,    min_version: "3.0",  alternatives: [macs2]}
  - {name: tabix,    optional: true}

python_packages:
  - {name: python,   min_version: "3.10"}
  - {name: pyyaml,   import_name: yaml,   min_version: "6.0"}  # pip-name vs import-name

r_packages:
  - {name: Seurat,        min_version: "5.0",  source: "bioconda"}
  - {name: MuDataSeurat,  source: "github",   repo: "PMBio/MuDataSeurat"}

data:
  - description: "fragments.tsv.gz"
    required_for: ATAC peak calling
    note: "supplied via --fragments at runtime"

notes: |
  Free-text additional context for the doctor's report or for
  human readers extracting this skill.
```

Field semantics:

- **`alternatives`**: a binary is satisfied if either `name` or any entry in
  `alternatives` is on PATH (e.g., macs3 OR macs2 satisfies the peak-caller
  requirement).
- **`optional`**: missing optional deps don't fail the doctor; they're
  reported as `[opt]` and ignored for the overall PASS/FAIL.
- **`import_name`**: for Python packages whose pip name differs from the
  import name (e.g. `pyyaml` -> `import yaml`), specify the import name here.
- **`source: github` + `repo`**: tells the doctor (and humans) that this
  package isn't installed by conda — it comes from `bin/postinstall.R`.

When you add a new dependency to a skill, update its `dependencies.yml`
entry too. A small CI check that greps each script for `subprocess.run`
calls and cross-checks against the YAML closes the "I forgot to update the
manifest" gap.
