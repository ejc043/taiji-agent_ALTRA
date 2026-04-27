# Environment management

This repo uses conda + lockfiles + (optionally) `conda-pack` to ship a
reproducible, "Docker-like" environment to SLURM nodes that don't allow
container runtimes. This doc walks through:

1. The one-command install (Tier 1).
2. Locking the env for byte-exact reproducibility (Tier 2).
3. Packing the env into a portable tarball (Tier 3).

## Tier 1 — One-command install

The end-to-end install:

```bash
bash bin/install.sh
```

That single command:

1. Creates the conda env `taiji-agent` from `environment.yml` using
   `micromamba` (5-10× faster than `conda`; falls back to `mamba`/`conda`
   if you pass `--solver`).
2. Runs `bin/postinstall.R` inside the env to install SeuratDisk and
   MuDataSeurat from GitHub (these are not on bioconda).
3. Runs `bin/install-taiji.sh` to download the Taiji binary for the
   current OS into `binaries/taiji`.

Common variants:

```bash
bash bin/install.sh --system centos          # explicit Taiji target
bash bin/install.sh --skip-taiji             # env only
bash bin/install.sh --skip-r                 # env + Taiji, no GitHub R packages
bash bin/install.sh --solver mamba           # if micromamba isn't installed
bash bin/install.sh --use-lockfile           # install from conda-lock.linux-64.yml
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
- After any `bin/postinstall.R` change (e.g. updating SeuratDisk).
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
- **The two GitHub-only R packages are baked in at install time.** Updating
  SeuratDisk or MuDataSeurat means rebuilding the tarball. That's the right
  tradeoff for reproducibility.
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
  - {name: SeuratDisk,    source: "github",   repo: "mojaveazure/seurat-disk"}

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
