# sh4oqs

Standalone **surface-hopping** reference dynamics for the spin-boson model: FSSH (MASH retained
but disabled) with a classical bath sampled per mode, for the harmonic/anharmonic and
linear/quadratic coupling cases.  The companion repositories are the private memory-kernel code
(`mkct4oqs`) and the public `dmrg4oqs` (HEOM / MPS-tDMRG; the shared
`source/config/` and `source/library/{qutils,fsutils,pltools,extended_deom,mre,corrmat}.py`).

Scientific background and methodology are described in
[`.opencode/skills/SKILL.md`](./.opencode/skills/SKILL.md); load the `surface-hopping` skill.

## Layout

All code is under `source/`, the working directory for the drivers (paths such as `./config` and
`./output` are relative to it):

- `source/sh.py` — driver: runs the selected methods/cases for a config and writes the
  basis correlation matrix `output/data/<prefix>/corrmat_<method>.npz (the observables are contracted from it at plot
  time by `library/corrmat.observables`, so only the correlation matrix is stored)` (`save_group`,
  key = case, `branch = component = ("00","01","10","11")`, `C(0) = I`).  The five reported
  observables are contractions of it (`library/corrmat.py`).
  (seed files averaged first); labelled sweep runs stay per case.
- `source/library/surface_hopping/` — `fssh.py` (fewest-switches) and `mash.py` (mapping approach),
  both enabled in `METHODS`, plus `validation.py` (regression suite).
- `source/sh_error.py` — accuracy of the surface-hopping runs against HEOM, with the relative
  Frobenius metric of `library/mre.py` (`||C - C'||_F / ||C'||_F`); the
  reference and the MKCT matrices come from the canonical `corrmat_heom.npz` /
  `corrmat_mkct.npz`; writes `output/sh_error.{json,md}`.
- `source/environment.yaml`, `source/requirements.txt` — the `sh` conda environment.

Shared pieces (`source/reference.py`, `source/config/`, and
`source/library/{qutils,fsutils,pltools,extended_deom,mre,corrmat}.py`) are **copies** of the canonical files in the private
`mkct4oqs` repository, so a public clone and the CI work without a token.

> **Sync reminder.**  After changing any of them in `mkct4oqs`, run
> `bash tools/sync_from_mkct.sh` here, then commit and push (same for `dmrg4oqs`).

## Conventions

- Run everything from `source/` with the `sh` environment active.
- Four cases: `QHM`, `QAH`, `CHM`, `CAH`; the anharmonic coefficient is `a = anharmon * sqrt(7/6)`.
- The classical cases use the classical per-mode variance (below the quantum zero-point).
- Linear coupling `V = sigma_z`; quadratic (`K-linear_quad`) `V = |1><1|`.
- Selection by `METHODS`/`CASES` or `-method`; trajectory budgets `N_TRAJ` (FSSH 20000,
  MASH 40000); flags `-n_traj -N -wmax_factor -dt -seed`.  Both methods take the inner pool,
  so `-nproc`/`-njobs` apply to either.
- Parallelism: `-nproc` runs the cases as separate processes (default 4), `-njobs` splits each
  case's ensemble over inner workers (default -1: share the cores).  The parent process writes
  every archive, so no two workers touch the same file.
- **Step size.** Harmonic runs use `-dt` (default 0.005).  The quartic term stiffens the fastest
  mode, so an anharmonic run caps the step at `dt <= 0.15 / w_max`; otherwise a rare
  high-energy trajectory near the Verlet stability limit can diverge (`w_max_factor = 10`
  gives `w_max = 120` and `dt = 0.00125` for the fast-bath scenario).
- Classical nuclei: a cross-check, reliable mainly for slow baths and high temperature.
- **Parameter sweep** (A-base_line, `mash`, documented labels, `dt = 0.01`): `base`
  (N=200, 40k trajectories, `wmax_factor` 10), `N400`/`N800`, `wf8`/`wf14`, `traj160k`
  (160k trajectories), and `all` (N=800, cut 14, 160k).  Each label is one `data.yaml`
  dispatch with `-label`; the per-case `corrmat_mash_<case>_<label>.npz` files stay unmerged
  and feed `sh_error`/`FIG7`.
- Coding style: see the maintainer's shared `repo-rules` skill.

## Workflows

```bash
cd source
python sh.py A-base_line -method fssh
python sh_merge.py --configs A-base_line --methods fssh
python library/surface_hopping/validation.py
```

## CI

- `.github/workflows/ci.yaml` — syntax check and the regression suite.
- `.github/workflows/data.yaml` — manual (`workflow_dispatch`) matrix over
  `<config> x <case> x <method>`; each job writes one combination and uploads
  `data/<config>/*.npz` as artifact `sh-<config>-<case>-<method>`, and a `combine` job packs the
  whole run into the `sh-data` tarball (both retained 30 days).  Inputs: `configs` (default
  `A-base_line` while the pipeline is still being tuned), `cases`, `methods` (default
  `fssh mash`), `n_traj`, `N`, `njobs`.  No automatic trigger: the configs live in the other
  repository, so a config change could not fire this workflow anyway.

```bash
gh workflow run data.yaml          # A-base_line x 4 cases x {fssh, mash} = 8 jobs
gh run watch <id> --exit-status
gh run download <id> -n sh-data -D /tmp
tar -xzf /tmp/sh-data/sh-data.tar.gz -C source/output/data/   # -> source/output/data/<config>/
```

- The repository is public, so runner minutes are free.  Keep the default config list at
  `A-base_line` until the pipeline is settled; pass `configs` to widen it (the whole `RUNLIST` is
  11 configs x 4 cases x 2 methods).

## Git style

See the maintainer's shared `repo-rules` skill (`type(scope): summary`; commit only when asked;
never commit secrets or `output/`).
