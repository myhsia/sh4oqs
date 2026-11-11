# Guard Rules

The following directories, files, and patterns MUST NOT be read or accessed under any circumstances:

- `.git/`
- `**/__pycache__/`

Any tool that reads files (`read`, `glob`, `grep`, `shell` with file-reading commands) must respect these restrictions.

# sh4oqs

Standalone **surface-hopping** reference dynamics for the spin-boson model: fewest-switches surface
hopping (FSSH) and the mapping approach to surface hopping (MASH), **both enabled**, propagated with
classical bath nuclei, for the quantum/classical and harmonic/anharmonic cases and the linear and
linear-quadratic couplings.  The companion repositories are the private memory-kernel code
(`mkct4oqs`) and the public `dmrg4oqs` (HEOM / MPS-tDMRG; the shared `source/config/` and the shared
library modules).

Scientific background and methodology are described in
[`.opencode/skills/SKILL.md`](./.opencode/skills/SKILL.md); load the `surface-hopping` skill.

## Shared rules

Repository etiquette (preserve uncommitted changes and staging, commit by pathspec, keep the guides
in sync), Python coding style, the git/PR workflow, LaTeX writing and source conventions,
tooling/upkeep, and local/CI parity live in the global `repo-rules` skill
(`~/.config/opencode/skills/repo-rules/SKILL.md`).

## Environment

1. A conda environment `sh` is created from `./source/environment.yaml`: the pinned interpreter
   `python=3.15` (channels `defaults`, `conda-pypi`) plus `./source/requirements.txt`, which stays
   unpinned (`numpy`, `numba`, `matplotlib`, `pyyaml`).
2. All commands assume the working directory is `./source` and the `sh` env is active; the drivers
   use paths relative to `./source` (e.g. `./config`, `./output`).
3. Entry points: `sh.py` (generate the data), `sh_error.py` (the accuracy table),
   `FIG0@corrmat.py` / `FIG4@robustness.py` (figures), and
   `library/surface_hopping/validation.py` (the regression suite).

## Project layout

- `source/sh.py` — driver: for one config it runs the selected methods over the four cases and
  writes **one archive per method**, `output/data/<prefix>/corrmat_<method>[_<label>].npz`
  (`save_group`: `t_<case>` plus the `(nt, 4, 4)` matrix under the case name,
  `branch = component = ("00","01","10","11")`, `C(0) = I`).  The four cases run as parallel
  processes and all of them are collected before the single write, so no post-processing merge is
  needed.
- `source/sh_error.py` — accuracy of the surface-hopping runs against HEOM: the relative Frobenius
  error `||C - C'||_F / ||C'||_F` over `[0, tf]` from `library/mre.py`, per config and method;
  `mkct` and the HEOM reference are read from the canonical checkout (`SH_REF_DATA`, default
  `../../mkct4oqs/source/output/data`); writes `output/sh_error.{json,md}`.
- `source/library/surface_hopping/fssh.py` — fewest-switches hopping; also the mean-field
  (Ehrenfest) variant, selected with `-method ehrenfest`.
- `source/library/surface_hopping/mash.py` — the mapping approach; also its mean-field variant,
  spin-LSC, selected with `-method lsc`.
- `source/library/surface_hopping/validation.py` — regression suite (exact single-mode and bare
  two-level references).
- `source/library/{qutils,fsutils,pltools,extended_deom,mre,corrmat}.py` — **synced** shared
  modules: quadrature and the Pauli basis, file IO (`save_group`/`load_group`), the window parser,
  the extended-DEOM reference, the error metric, the correlation-matrix basis/contractions, and
  the label pinning (`_pin_ylabel`/`_pin_text`/`_pin_xlabel`).
  The five reported observables are contractions of `C` (`corrmat.py`), computed at plot time only.
- `source/reference.py` — **synced** reference driver (HEOM / extended DEOM, writes
  `corrmat_heom.npz`).
- `source/config/` — **synced** configs (`A-base_line` … `K-linear_quad`, plus the retired
  `A-base_line_tf6` test config).  `mkct4oqs` owns them.
- `source/FIG0@corrmat.py`, `source/FIG4@robustness.py` — **synced** canonical figure scripts:
  `FIG0` draws the 4x4 `C_ij(t)` pages of every method found here and `FIG4` the
  `<method>@<prefix>.pdf` comparison figures, overlaying the reference `corrmat_heom.npz` that
  `tools/link_heom.sh` places in `source/output/data/<prefix>/`.  Both read the data root
  `./output/data` by default (`-ref DIR` for another root; without a reference the `qhm` pages are
  drawn alone).  `FIG0` also takes `-label TAG` (a tagged run) and `-tmax T` (axis window, file
  suffixed `_t<T>`); `FIG4` draws a group that holds only some of the four cases from the cases
  present, so a single-case method still yields its reduced pages.  Both flags are defined in
  `mkct4oqs/AGENTS.md`.
- `source/environment.yaml`, `source/requirements.txt` — the `sh` environment.
- `.github/workflows/ci.yaml`, `.github/workflows/data.yaml`, `.github/actions/sh-env/action.yml` —
  CI (see below).
- `tools/sync_from_mkct.sh` — copies the shared pieces from `../mkct4oqs/source`.
- `source/tools/link_heom.sh` — links the reference `corrmat_heom.npz` (generated in `dmrg4oqs`)
  into `source/output/data/<prefix>/`; run it once per checkout after the data land.
- `README.md`, `LICENSE`.

> **Sync reminder.**  After changing any of the shared pieces in `mkct4oqs`, run
> `bash tools/sync_from_mkct.sh` here, then commit and push (same for `dmrg4oqs`).

## Conventions and pitfalls

- Set `sys.pycache_prefix = os.path.abspath("./output/pycache")` **after** third-party imports and
  **before** the first project import in every run script.
- Four cases: `QHM`, `QAH`, `CHM`, `CAH` (`_CASE_MATRIX` in `sh.py`: quantum sampling × anharmonic).
  The anharmonic quartic coefficient is `a = anharmon * sqrt(7/6)`.
- The quantum cases draw the Wigner (quantum thermal) variance per mode, the classical cases the
  Boltzmann variance, which lies below the quantum zero-point; `-sampling {auto,quantum,classical}`
  overrides the per-case choice (`auto` is the default).
- Linear coupling `V = sigma_z`; quadratic (`K-linear_quad`) `V = |1><1|`.
- Selection: the `METHODS` / `CASES` lists in `sh.py` (comment out entries) or `-method`
  (`fssh`, `mash`, `lsc`, `ehrenfest`, `all`).  `RUNLIST` holds only `A-base_line` for now; the rest
  of the demos are commented out.
- Flags: `-n_traj -N -wmax_factor -dt -seed -sampling -label -nproc -njobs`.
- Trajectory budgets: `N_TRAJ` (the defaults, used when `-n_traj` is absent) is `mash` 40000,
  `fssh` 20000, `lsc` 40000, `ehrenfest` 20000.  The **production** setting is larger and is passed
  explicitly: `N = 3200`, `wmax_factor = 10`, `dt = 0.01`, **MASH 480 000**, **FSSH 120 000**.
- **Step size.** `-dt` is **0.01** by default, in the parser, its help text, `run()`'s signature and
  the `data.yaml` input description alike.  An anharmonic run additionally caps the step at
  `dt <= 0.4 / w_max`; without that, a rare high-energy trajectory diverges at the velocity-Verlet
  stability limit.
- Parallelism: `-nproc` runs the cases as separate processes (default 4), `-njobs` splits each
  case's ensemble over inner workers (default -1: share the cores over the cases).  Workers return
  results and the parent writes every archive, so no two processes write the same file.  With one
  case per process (`nproc >= #cases`) the inner pool is skipped, since nested pools only add
  overhead.
- `numba` is a declared dependency: the trajectory kernels import `njit` directly, with no
  `except ImportError` fallback decorator, so a missing numba fails at import instead of silently
  running the kernels as pure Python.
- Classical nuclei: a cross-check, not an exact reference; reliable mainly for slow baths and high
  temperature.
- Archive schema is shared with MKCT: `corrmat_<method>[_<label>].npz`, one matrix per case, same
  basis — so `sh_error.py` and the synced figure scripts need no per-observable files.

## Workflows

```bash
cd source
python sh.py A-base_line -method mash                    # four cases, one corrmat archive
python sh.py A-base_line -method mash -n_traj 480000 -N 3200 -dt 0.01   # production point
python sh.py A-base_line -method mash -label tst -n_traj 500 -N 60      # labelled smoke test
python sh_error.py --configs A-base_line                 # accuracy table
bash tools/link_heom.sh                                 # reference links (once)
python FIG0@corrmat.py
python FIG4@robustness.py
python library/surface_hopping/validation.py
```

## Validation

`python library/surface_hopping/validation.py` runs five checks against numerically exact
references and exits non-zero if any fails:

| check | tolerance |
|---|---|
| bare two-level, MASH estimator + `sigma_y` sign | 0.05 |
| bare two-level, FSSH density matrix | 0.05 |
| one harmonic mode (NAC sign, hop rescaling) | 0.13 |
| one quartic mode (anharmonic bath) | 0.06 |
| one QSBM mode (quadratic coupling) | 0.08 |

Each line prints `[PASS]`/`[FAIL]` with the achieved maximum error.  The bare-system error falls as
`1/sqrt(n_traj)`, so a rising floor there means an estimator bias, not noise.  The suite runs in CI
on every push (see below).

## CI
**Dependencies: the bot proposes, CI decides** (the policy and its reasoning are in the
machine-global `repo-rules` skill, "Local and CI parity").  `.github/dependabot.yml` watches
`github-actions` daily and, for the `conda` ecosystem, ignores **patch** releases of `python` -- no
functional change, and the class most likely to break the compiled extensions.  The pin in
`source/environment.yaml` moves deliberately, in step with the sibling repositories, and only once
`build_env` resolves the environment and the suite is green.  `requirements.txt` is intentionally
unpinned, so the libraries float to latest at install time.

- `.github/workflows/ci.yaml` — `Syntax (compileall)`, `Build environment (conda)`, and the
  `Surface-hopping suite`; triggers are push/PR (Markdown ignored), the daily 09:41 UTC schedule,
  and manual dispatch.
- `.github/workflows/data.yaml` — **manual only** (`workflow_dispatch`; the configs live in the
  other repository, so a config change could not fire it).  Shape: `plan` → `build_env` →
  a matrix over **`<config> × <method>`** (`max-parallel: 20`, `timeout-minutes: 350`) → `combine`.
  Each job runs the four cases inside itself as parallel processes and uploads
  `source/output/data/` as artifact `sh-<config>-<method>`; `combine` packs every artifact of the
  run into the `sh-data` tarball (30-day retention for both).  `run-name` puts the config and
  method in the run title.
- `.github/workflows/collect.yaml` — **manual only** (`workflow_dispatch`): packs the `sh-data`
  artifacts of several Data runs (input `runs`, space-separated ids) into one tarball, so a batch
  dispatched as per-config runs is a single download.
- Inputs: `configs` (default `A-base_line`), `methods` (default `fssh mash`), `n_traj` (empty: the
  `N_TRAJ` budget), `N`, `njobs`, `dt` (empty: the driver's default), `wmax`, `label`, `sampling`.
  The `cases` input still exists and is exported by `plan`, but it is **unused**: the matrix has no
  case dimension and `sh.run` always runs all four `CASES`.
- The repository is public, so runner minutes are free; the production data are generated on CI
  rather than locally.  A method that needs more than the 350-min job budget must be split by
  `label` or have `n_traj` reduced (see Validity).

```bash
gh workflow run data.yaml                        # A-base_line × {fssh, mash} = 2 jobs
gh run watch <id> --exit-status
PROXY=$(cat ~/.config/opencode/proxy)            # artifacts download through the proxy
https_proxy="$PROXY" http_proxy="$PROXY" gh run download <id> -n sh-data -D /tmp
tar -xzf /tmp/sh-data.tar.gz -C source/output/data/   # -> source/output/data/<config>/
```

Passing `configs="<all eleven>"` gives 11 × 2 = 22 jobs, one per `(config, method)`.

## Git style

See the maintainer's shared `repo-rules` skill (`type(scope): summary`; commit only when asked;
never commit secrets or `output/`).

## Validity (measured, not assumed)

Labels: **measured** = a run cited by its artifact/command; **derived** = a stated formula;
**recorded** = a number kept in a committed or ignored output file, re-read rather than re-run.

- **Step size `dt = 0.01`** — derived: the anharmonic guard keeps `w_max * dt <= 0.4`, inside the
  velocity-Verlet stability bound; the cap was `0.15 / w_max` before and was relaxed because the
  step only smoothed the curves while dominating wall time.  For `E-fast_bath` (`wC = 12`,
  `w_max = 120`) the guard gives `dt = 0.0033`.
- **Case → sampling mapping** — derived, not fitted: `qhm/qah` use the quantum marginal, `chm/cah`
  the classical one, so the four cases separate the two approximations (sampling and quartic term).
- **Production budgets MASH 480 000 / FSSH 120 000** — measured on CI: MASH at `N = 3200`,
  `dt = 0.01`, 480 000 fits the 350-min job; FSSH costs about nine initial states per sample, so
  480 000 exceeded the cap and 120 000 fits.
- **`N = 3200`, `wmax_factor = 10`** — measured: a ten-point sweep of `N` (200 → 800), the cutoff
  (8 → 14), `dt` (5e-3 → 2.5e-3) and the ensemble (40 000 → 160 000), individually and combined,
  left the error flat, and later runs reached `N = 3200` and cutoff 30 without moving it.  That
  sweep was scored with the earlier observable-contraction metric, so it is quoted qualitatively;
  the largest setting is used for production to put all the studied axes at their limit.
- **Accuracy, `A-base_line`** — recorded in `output/sh_error.md` under the Frobenius metric:
  MASH 8.3 %, MKCT 17.4 %, FSSH 53.9 %.  MASH being better than MKCT here is a property of that
  demo, not a general ordering.
- **Classical-nuclei residual** — measured: the deviation from HEOM shrinks as `beta` falls and
  grows as it rises, i.e. it is the classical-nuclei approximation, not a convergence parameter.
  No numerical setting removes it; report it as the method's regime rather than tuning it away.

