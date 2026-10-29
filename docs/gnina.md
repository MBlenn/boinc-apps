# gnina (GPU docking test app)

[gnina](https://github.com/gnina/gnina) docks small molecules into protein structures and
rescores the poses with convolutional neural networks (CNN) on the GPU. This app is a **test
app**: it runs gnina on a handful of well-known systems from gnina's own test data, so that
results can be checked automatically. It exercises the full chain (GPU in the container, input
files, output upload, a real validator) before production docking campaigns.

## Job mix

`boinc-app submit gnina N` creates N jobs (`apps/gnina/jobs.py`):

| Share | Kind | What it does | Validation |
|---|---|---|---|
| ~60 % | `redock` | native ligand into its own receptor, default model, random seed | plausibility checks **and** RMSD to the crystal pose below the per-system limit |
| ~30 % | `cross` | ligand of system A into receptor B (A may equal B), random CNN model, scoring mode and seed | plausibility checks only |
| ~10 % | `broken` | truncated ligand file | must fail on the client (exit code 2); one attempt only |

Workunit names encode the kind and the systems: `gnina_<kind>_<receptor>_<ligand>_<stamp>_<n>`.
The estimated run time (`rsc_fpops_est`) depends on the scoring mode: `rescore` 1×,
`metrorescore` 3×, `refinement` 5×, `metrorefine` 6×, `all` 20×.

Test systems (`apps/gnina/testdata/systems.txt`): 184l, 10gs, 3rod, grid, noelem. The
redocking limits (2.0 to 3.5 Å) come from a measured baseline (gnina 1.3.3 on the CPU, several
seeds): 184l 0.5–0.9 Å, 3rod 0.5 Å, grid 0.6 Å, 10gs 0.8–2.3 Å, noelem 2.2–2.6 Å.

## Validator

`script_validator` calls `apps/gnina/server/gnina_validate.py` (copied to the project's `bin/`)
for every result. It checks:

- `docked.sdf` contains at least one pose
- every pose has a numeric `minimizedAffinity` (−50 to 100) and a `CNNscore` (0 to 1)
- every pose has as many heavy atoms as the input ligand
- the pose centroids lie within 12 Å of the box ligand's centroid
- redocking: at least one pose within the system's RMSD limit of the crystal ligand
  (atoms of the same element are matched optimally, no superposition)

Exit code 0 = valid, 1 = invalid, 3 = database not reachable (BOINC retries later). One line per
result goes to `validator_gnina.log` in the project's `log_<host>/` directory:

```bash
tail -f /srv/boinc/<project>/boinc/<project>/log_<project>-daemons/validator_gnina.log
```

Manual check of a result file (inside the daemons container, in the project directory):

```bash
python3 bin/gnina_validate.py --wu gnina_redock_184l_184l_x_0 path/to/docked.sdf
```

`GNINA_RMSD_MAX` overrides the per-system limit, `GNINA_TESTDATA` the test data directory
(default: `boinc_apps/gnina/testdata` in the project, created by `boinc-app install`).

## Client requirements

- Linux x86_64, BOINC client 8.2 or newer
- NVIDIA GPU with a driver that supports CUDA 12.8 (driver 570 or newer recommended)
- Podman or Docker usable by the client user, NVIDIA Container Toolkit (CDI for Podman), see
  [boinc-stack docs/docker-apps.md](https://github.com/MBlenn/boinc-stack/blob/main/docs/docker-apps.md)
- several GB of disk space for the image (pulled once per host and gnina version)

Plan class `docker_nvidia_cpu1`: one GPU and one CPU core per task.

## App image

`app-images/gnina/Containerfile` builds `ghcr.io/mblenn/boinc-app-gnina:<gnina version>`:
gnina's static release binary plus the CUDA 12 / cuDNN 9 runtime libraries from NVIDIA's PyPI
wheels (pinned versions; the official gnina image is far larger). The driver library itself is
injected by the NVIDIA Container Toolkit. The CI job runs a CPU redocking of 184l and validates
the result.

After the first build, set the GHCR package `boinc-app-gnina` to public.

The image is tagged `<gnina version>-<IMAGE_REV>` (e.g. `1.3.3-2`) and `<gnina version>`.
Clients build the app's Dockerfile `FROM` a fixed tag and never pull a tag again that they
already have, so every image change needs a new tag. To update the image or gnina:

1. change `IMAGE_REV` (and `GNINA_VERSION` for a new gnina release) in the workflow
2. use the new tag in `apps/gnina/files/Dockerfile`
3. raise `APP_VERSION` in `apps/gnina/app.env` (new app version, new local image name)
4. run `boinc-app install gnina --from ...` again on the server

PyTorch inside gnina loads `libnvrtc` only at run time for CNN scoring on the GPU. The image
build checks that it is present, because CPU runs (including the CI smoke test) do not need it.
