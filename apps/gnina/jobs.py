"""Job generator for the gnina test app, used by "boinc-app submit gnina N".

Mix (see docs/gnina.md):
  ~60 % redock  native ligand into its own receptor, default model and scoring, random seed.
                Validated against the crystal pose (RMSD, see server/gnina_validate.py).
  ~30 % cross   ligand of system A into receptor B (also A == B), random CNN model, scoring
                mode and seed. Variety and different run times; checked for plausibility only.
  ~10 % broken  truncated ligand file. Must end as a client error, never as a valid result.

make_jobs() returns dicts that boinc-app turns into create_work calls:
  name     unique workunit name, gnina_<kind>_<receptor>_<ligand>_...
  inputs   test data files (receptor, ligand, box) in template order
  args     gnina options, passed to the container as $ARGS
  options  extra create_work options
"""
import os

# built-in CNN models of gnina 1.3 (gnina --help); None = default ensemble
MODELS = [None, "crossdock_default2018", "general_default2018_1", "redock_default2018_1_3",
          "all_default_to_default_1_3_1", "crossdock_default2018_1_3_2"]
# (scoring mode, weight, run time factor relative to rescore)
MODES = [("rescore", 50, 1), ("refinement", 25, 5), ("metrorescore", 10, 3),
         ("metrorefine", 10, 6), ("all", 5, 20)]
BASE_FPOPS = 3e12   # rsc_fpops_est of a default redocking job


def _options(fpops_factor=1, max_errors=3, max_total=5):
    return ["--rsc_fpops_est", "%g" % (BASE_FPOPS * fpops_factor),
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
        r = rng.random()
        if r < 0.6:
            code = rng.choice(codes)
            s = systems[code]
            jobs.append({
                "name": f"gnina_redock_{code}_{code}_{stamp}_{i}",
                "inputs": [s["rec"], s["lig"], s["lig"]],
                "args": f"--seed {seed} --num_modes 9 --cpu 1",
                "options": _options(),
            })
        elif r < 0.9:
            rec, lig = rng.choice(codes), rng.choice(codes)
            model = rng.choice(MODELS)
            mode, factor = _pick_mode(rng)
            args = f"--seed {seed} --num_modes 9 --cpu 1 --cnn_scoring {mode}"
            if model:
                args += f" --cnn {model}"
            jobs.append({
                "name": f"gnina_cross_{rec}_{lig}_{stamp}_{i}",
                "inputs": [systems[rec]["rec"], systems[lig]["lig"], systems[rec]["lig"]],
                "args": args,
                "options": _options(factor),
            })
        else:
            jobs.append({
                "name": f"gnina_broken_10gs_broken_{stamp}_{i}",
                "inputs": [systems["10gs"]["rec"], "broken_lig.sdf", systems["10gs"]["lig"]],
                "args": f"--seed {seed} --num_modes 9 --cpu 1",
                # one attempt is enough: the error is expected
                "options": _options(max_errors=1, max_total=1),
            })
    return jobs
