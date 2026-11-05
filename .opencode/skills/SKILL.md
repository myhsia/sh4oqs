---
name: surface-hopping
description: FSSH and MASH surface-hopping reference dynamics for the spin-boson model (harmonic and anharmonic bath, linear and quadratic coupling) — use when generating, scoring or interpreting the classical-nuclei reference data
license: MIT
compatibility: opencode
metadata:
  audience: maintainers
  workflow: github
---

# Surface-hopping reference dynamics

An independent, classical-nuclei reference for the same spin-boson model as the memory-kernel
calculation: the electronic two-level system is propagated on the adiabatic surfaces (FSSH) or
through the mapping variables (MASH) while the bath is sampled classically.

## Model

The four cases `QHM` (quantum harmonic), `QAH` (quantum anharmonic), `CHM` (classical harmonic)
and `CAH` (classical anharmonic), with the Ohmic density
`J(omega) = pi*lambd/wC * omega * exp(-omega/wC)` discretized by the shared Gauss-Legendre modes.
The linear coupling is `V = sigma_z`; the quadratic case (`K-linear_quad`) uses `V = |1><1|`.
The anharmonic quartic coefficient is `a = anharmon * sqrt(7/6)`.

The initial bath is sampled per mode: a **quantum** case draws from the Wigner/Gaussian
distribution with the quantum thermal variance, a **classical** case from the classical
(Boltzmann) variance, which lies below the quantum zero-point.  The dynamics is the same for both;
only the initial ensemble differs.  `-sampling {auto,quantum,classical}` overrides the per-case
choice.

## Methods

- `library/surface_hopping/fssh.py` — fewest-switches surface hopping, and its mean-field variant
  (`-method ehrenfest`).
- `library/surface_hopping/mash.py` — the mapping approach to surface hopping, and its mean-field
  variant (`-method lsc`).  Both `fssh` and `mash` are enabled (`METHODS = ["fssh", "mash"]` in
  `sh.py`); the two mean-field variants exist as comparison methods and share the same machinery.
- The drivers use numba-jitted kernels; the trajectory loop is the expensive part.  Ensembles are
  split into independent slices with deterministic per-slice seeds (`n_jobs`), and the slice
  workers return unnormalized sums that the parent combines, so results do not depend on the number
  of workers.

`sh.py` runs the selected methods and cases for each config in its `RUNLIST` and writes one basis
correlation matrix archive per method, `output/data/<prefix>/corrmat_<method>[_<label>].npz`
(`save_group`, one `(nt, 4, 4)` matrix per case, `C(0) = I`, `<method>` = `fssh`, `mash`, `lsc` or
`ehrenfest`).  The five observables are contractions of `C` (`library/corrmat.py`), computed at plot
time only; `sh_error.py` scores `C` against HEOM by the relative Frobenius error
`||C - C'||_F / ||C'||_F` over the calculation window `[0, tf]`, in percent.  The synced
`FIG0@corrmat.py` / `FIG4@robustness.py` draw this repository's figures, e.g.
`python FIG0@corrmat.py -ref ../../mkct4oqs/source/output/data` and
`python FIG4@robustness.py -ref ../../mkct4oqs/source/output/data`
(`source/output/figures/`), overlaying the reference when the private checkout is present.

## Conventions

- Run from `./source` with the `sh` conda environment active.
- Selection is by the `METHODS` / `CASES` lists in `sh.py` (comment out entries) or by the
  `-method` CLI flag; `-n_traj` defaults to the per-method `N_TRAJ` budget
  (`mash` 40000, `fssh` 20000, `lsc` 40000, `ehrenfest` 20000).  The **production** setting is
  larger and passed explicitly: `N = 3200`, `wmax_factor = 10`, `dt = 0.01`, MASH 480 000,
  FSSH 120 000.  Flags: `-n_traj -N -wmax_factor -dt -seed -sampling -label -nproc -njobs`.
- **Step size:** canonical `-dt` is **0.01**; an anharmonic run caps it at `dt <= 0.4 / w_max`
  (the velocity-Verlet bound), e.g. `0.0033` for the fast-bath scenario, where `w_max = 120`.
  Note the `argparse` default is still `0.005` although the help text, `run()` and the workflow
  description say `0.01`; pass `-dt 0.01` until the parser is corrected.
- Parallelism: `-nproc` (default 4) runs the cases as separate processes; `-njobs` (default -1)
  splits each case's ensemble over inner workers, sharing the cores over the cases.  The parent
  process writes every archive, so no two processes write the same file.  With one case per process
  (`nproc >= #cases`) the inner pool is skipped, since nested pools only add overhead.
- Archive naming: one file per `(config, method)`, all four cases inside, same basis and keys as
  the canonical `corrmat_mkct.npz` / `corrmat_heom.npz`, so the two codes are compared without any
  per-observable intermediate file.
- **Caveat:** classical nuclei, so this is a cross-check rather than an exact reference; it is
  reliable mainly for slow baths and high temperature, and its deviation from HEOM grows as
  `beta` rises.  That residual is the method's regime, not a convergence parameter.
- The shared pieces (`reference.py`, `config/`, `FIG0@corrmat.py`, `FIG4@robustness.py`,
  `qutils.py`, `fsutils.py`, `pltools.py`, `extended_deom.py`, `mre.py`, `corrmat.py`) are copies of
  the canonical files in `mkct4oqs`; after a change there, run `tools/sync_from_mkct.sh` and
  commit.
- System-wide conventions (repo etiquette, Python style, git/PR, LaTeX, tooling) are in the
  maintainer's global `repo-rules` skill.

## Workflows and validation

```bash
cd source
python sh.py A-base_line -method mash                    # one config, all four cases
python sh.py A-base_line -method mash -n_traj 480000 -N 3200 -dt 0.01   # production point
python sh_error.py --configs A-base_line                 # Frobenius table -> output/sh_error.{json,md}
python library/surface_hopping/validation.py             # regression suite
```

The suite compares against numerically exact references and exits non-zero on failure:

| check | tolerance |
|---|---|
| bare two-level, MASH estimator + `sigma_y` sign | 0.05 |
| bare two-level, FSSH density matrix | 0.05 |
| one harmonic mode (NAC sign, hop rescaling) | 0.13 |
| one quartic mode (anharmonic bath) | 0.06 |
| one QSBM mode (quadratic coupling) | 0.08 |

The bare-system error falls as `1/sqrt(n_traj)`, so a rising floor there indicates an estimator
bias rather than Monte-Carlo noise; the coupled single-mode checks bound the method error instead.
CI runs the suite on every push.

CI: `.github/workflows/ci.yaml` (syntax + environment + regression suite) and
`.github/workflows/data.yaml` (manual, matrix over `<config> x <method>`; each job runs the four
cases internally and uploads `source/output/data/` as `sh-<config>-<method>`; a `combine` job packs
the run into the `sh-data` tarball, retained 30 days).  The matrix has no case dimension — the
workflow's `cases` input is currently unused.  Defaults: `configs=A-base_line`, `methods="fssh
mash"`.  The repository is public, so the runner minutes are free and the production data are
generated on CI rather than locally:

```bash
gh workflow run data.yaml                       # A-base_line x {fssh, mash} = 2 jobs
gh run watch <id> --exit-status
PROXY=$(cat ~/.config/opencode/proxy)
https_proxy="$PROXY" http_proxy="$PROXY" gh run download <id> -n sh-data -D /tmp
tar -xzf /tmp/sh-data.tar.gz -C source/output/data/   # -> source/output/data/<config>/
```

A job has 350 minutes; at `N = 3200` MASH fits 480 000 trajectories while FSSH (about nine initial
states per sample) fits only about 120 000.  There is no automatic trigger, because the configs are
owned by the other repository.
