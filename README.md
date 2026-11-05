# sh4oqs

Reference **surface-hopping** dynamics for the spin-boson model used by the MKCT manuscript:
fewest-switches surface hopping (FSSH) and the mapping approach to surface hopping (MASH), both
enabled, with classical bath nuclei, covering the quantum/classical and harmonic/anharmonic cases
with linear and linear-quadratic coupling.

## Install

```bash
conda env create -f source/environment.yaml
conda activate sh
```

## Run

From `source/`:

```bash
python sh.py A-base_line -method mash                    # all four cases, one archive
python sh.py A-base_line -method mash -n_traj 480000 -N 3200 -dt 0.01   # production point
```

The four cases run as parallel processes (`-nproc`, default 4) and, when there are spare cores,
each case splits its ensemble over inner workers (`-njobs`, default -1: share the cores).

Results are written to `output/data/<prefix>/corrmat_<method>[_<label>].npz` — one basis
correlation matrix per case (`t_<case>` plus a `(nt, 4, 4)` matrix), the same schema as the MKCT
archives.  The five observables are contractions of that matrix, computed at plot time.

## CI data generation

The full production run is cheaper on GitHub Actions: dispatch the `Data` workflow and download
the results (public repository, so the runner minutes are free).

```bash
gh workflow run data.yaml                     # A-base_line x {fssh, mash} = 2 jobs
gh run watch <id> --exit-status
PROXY=$(cat ~/.config/opencode/proxy)         # artifacts download through the proxy
https_proxy="$PROXY" http_proxy="$PROXY" gh run download <id> -n sh-data -D /tmp
tar -xzf /tmp/sh-data.tar.gz -C source/output/data/   # -> source/output/data/<config>/
```

Each job handles one `<config> x <method>`, runs the four cases internally, and uploads
`data/<config>/*.npz`; the `combine` job packs every job into the single `sh-data` tarball
(retained 30 days).  Both methods (FSSH and MASH) are enabled; the other scenarios are commented
out in `RUNLIST` for now.  A job has 350 minutes, so MASH fits 480 000 trajectories at `N = 3200`
while FSSH (about nine initial states per sample) fits about 120 000.

## Regression suite

```bash
python library/surface_hopping/validation.py
```

## Shared pieces

`source/reference.py`, `source/config/`, `source/FIG0@corrmat.py`, `source/FIG4@robustness.py` and
`source/library/{qutils,fsutils,pltools,extended_deom,mre,corrmat}.py` are copies of
the canonical files in the private `mkct4oqs` repository.  After a change there, run
`bash tools/sync_from_mkct.sh` and commit (the public
[dmrg4oqs](https://github.com/myhsia/dmrg4oqs) mirror does the same).

## Environment

`conda env create -f source/environment.yaml` builds the `sh` environment (pinned interpreter,
unpinned `requirements.txt`).
