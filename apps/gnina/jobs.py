"""Job generator for the gnina test apps, used by "boinc-app submit <app> N".

The same file serves both apps (apps/gnina_cpu/jobs.py is a link to it); the variant
follows from the directory it is loaded from:
  gnina      NVIDIA GPU (plan class docker_nvidia_cpu1): CNN-heavy scoring modes, so the
             GPU is actually used (rescore only touches the GPU for a fraction of a second)
  gnina_cpu  CPU only (plan class docker_cpu1, no GPU reserved): --no_gpu, CNN rescoring

Mix (see docs/gnina.md):
  ~60 % redock  native ligand into its own receptor, default model, random seed.
                Validated against the crystal pose (RMSD, see server/gnina_validate.py).
  ~30 % cross   ligand of system A into receptor B (also A == B), random CNN model and seed
                (GPU: random scoring mode). Checked for plausibility only.
  ~10 % broken  truncated ligand file. Must end as a client error, never as a valid result.

make_jobs() returns dicts that boinc-app turns into create_work calls:
  name     unique workunit name, <app>_<kind>_<receptor>_<ligand>_...
  inputs   test data files (receptor, ligand, box) in template order
  args     gnina options, passed to the container as $ARGS
  options  extra create_work options
"""
import os

APP = os.path.basename(os.path.dirname(os.path.abspath(__file__)))
CPU_ONLY = APP.endswith("_cpu")

# built-in CNN models of gnina 1.3 (gnina --help); None = default ensemble
MODELS = [None, "crossdock_default2018", "general_default2018_1", "redock_default2018_1_3",
          "all_default_to_default_1_3_1", "crossdock_default2018_1_3_2"]

# (scoring mode, weight, run time factor). The CNN modes beyond "rescore" take minutes
# per job on a CPU core (10gs: rescore 2.5 min, refinement > 10 min), so the CPU app
# only rescores.
if CPU_ONLY:
    MODES = [("rescore", 1, 1)]
    REDOCK_MODE = ("rescore", 1)
    EXTRA_ARGS = " --no_gpu"
    # rsc_fpops_est of a factor-1 job and the abort limit relative to it.
    # CPU: ~1 min on a typical core (BOINC measures ~5e10 FLOPS per core).
    BASE_FPOPS, BOUND_FACTOR = 3e12, 1000
else:
    MODES = [("refinement", 40, 1), ("metrorescore", 20, 0.6), ("metrorefine", 25, 1.2),
             ("all", 15, 4)]
    REDOCK_MODE = ("refinement", 1)
    EXTRA_ARGS = ""
    # GPU: BOINC divides by the GPU's peak FLOPS (tens of TFLOPS), so the estimate has to
    # be large; a too small bound aborts tasks with "exceeded elapsed time limit".
    BASE_FPOPS, BOUND_FACTOR = 1e16, 200


def _options(fpops_factor=1, max_errors=3, max_total=5):
    est = BASE_FPOPS * fpops_factor
    return ["--rsc_fpops_est", "%g" % est, "--rsc_fpops_bound", "%g" % (est * BOUND_FACTOR),
            "--max_error_results", str(max_errors), "--max_total_results", str(max_total)]


def load_systems(testdata):
    systems = {}
    for line in open(os.path.join(testdata, "systems.txt"), encoding="utf-8"):
        f = line.split("#")[0].split()
        if len(f) == 4:
            systems[f[0]] = {"rec": f[1], "lig": f[2]}
    return systems


def _pick_mode(rng):
    r = rng.uniform(0, sum(w for _, w, _ in MODES))
    for mode, w, factor in MODES:
        if r < w:
            return mode, factor
        r -= w
    return MODES[0][0], MODES[0][2]


def make_jobs(count, rng, stamp, testdata):
    systems = load_systems(testdata)
    codes = sorted(systems)
    jobs = []
    for i in range(count):
        seed = rng.randint(1, 2**31 - 1)
        base = f"--seed {seed} --num_modes 9 --cpu 1{EXTRA_ARGS}"
        r = rng.random()
        if r < 0.6:
            code = rng.choice(codes)
            s = systems[code]
            mode, factor = REDOCK_MODE
            jobs.append({
                "name": f"{APP}_redock_{code}_{code}_{stamp}_{i}",
                "inputs": [s["rec"], s["lig"], s["lig"]],
                "args": f"{base} --cnn_scoring {mode}",
                "options": _options(factor),
            })
        elif r < 0.9:
            rec, lig = rng.choice(codes), rng.choice(codes)
            model = rng.choice(MODELS)
            mode, factor = _pick_mode(rng)
            args = f"{base} --cnn_scoring {mode}"
            if model:
                args += f" --cnn {model}"
            jobs.append({
                "name": f"{APP}_cross_{rec}_{lig}_{stamp}_{i}",
                "inputs": [systems[rec]["rec"], systems[lig]["lig"], systems[rec]["lig"]],
                "args": args,
                "options": _options(factor),
            })
        else:
            jobs.append({
                "name": f"{APP}_broken_10gs_broken_{stamp}_{i}",
                "inputs": [systems["10gs"]["rec"], "broken_lig.sdf", systems["10gs"]["lig"]],
                "args": base,
                # one attempt is enough: the error is expected
                "options": _options(max_errors=1, max_total=1),
            })
    return jobs
