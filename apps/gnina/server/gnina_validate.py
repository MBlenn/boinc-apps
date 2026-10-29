#!/usr/bin/env python3
"""Validator for gnina results (called by BOINC's script_validator as init script).

    script_validator --app gnina --init_script gnina_validate.py

script_validator passes the paths of the result's output files. Only a bare script name
works: BOINC's split() keeps the delimiter in each part ("files " != "files"), so the
argument keywords in --init_script "script files result_id" are never recognized and
the script is called without arguments (BOINC server 1.6.2, still in master 10/2026).
The workunit name is taken from the output file names instead:
  <workunit>_<instance>_r<random>_<file index>

Exit code 0 = valid, 1 = invalid, 3 = cannot decide yet (server-side problem: BOINC retries later).
Prints one line per result to the validator log.

Checks for every job:
  - docked.sdf contains at least one pose
  - every pose has numeric minimizedAffinity and CNNscore (0..1)
  - every pose has as many heavy atoms as the input ligand
  - the poses lie inside the docking box (centroid near the box ligand)
Additionally for redocking jobs (gnina_redock_*):
  - at least one pose within rmsd_max (per system, systems.txt) of the crystal ligand
    (symmetry-tolerant RMSD: atoms of the same element are matched optimally)

The job type and the systems come from the workunit name, created by boinc-app submit:
  gnina_<kind>_<receptor>_<ligand>_<...>   kind = redock | cross | broken
"""
import math
import os
import re
import sys

# boinc-app install copies this script to <project>/bin/ and the app definition
# (including testdata/) to <project>/boinc_apps/gnina/
_PROJECT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TESTDATA = os.environ.get("GNINA_TESTDATA", os.path.join(_PROJECT, "boinc_apps", "gnina", "testdata"))
RMSD_MAX = os.environ.get("GNINA_RMSD_MAX")  # overrides the per-system limit in systems.txt
CENTROID_MAX = 12.0  # A, distance of a pose centroid from the box ligand centroid


# ---------------------------------------------------------------- parsing

def element(sym):
    sym = sym.strip()
    return sym[:1].upper() + sym[1:].lower() if sym else "*"


def parse_sdf(text):
    """Return list of (atoms, props); atoms = [(element, x, y, z)]."""
    mols = []
    for block in text.split("$$$$"):
        lines = block.strip("\n").splitlines()
        if len(lines) < 4:
            continue
        # counts line is the first line ending in V2000 (header may have empty lines)
        ci = next((i for i, l in enumerate(lines) if l.rstrip().endswith("V2000")), None)
        if ci is None:
            continue
        n = int(lines[ci][:3])
        atoms = []
        for l in lines[ci + 1:ci + 1 + n]:
            atoms.append((element(l[31:34]), float(l[0:10]), float(l[10:20]), float(l[20:30])))
        props = dict(re.findall(r">\s*<([^>]+)>[^\n]*\n([^\n]*)", block))
        mols.append((atoms, props))
    return mols


def parse_pdb(text):
    atoms = []
    for l in text.splitlines():
        if l.startswith(("ATOM", "HETATM")):
            el = l[76:78].strip() or re.sub(r"[^A-Za-z]", "", l[12:16])[:1]
            atoms.append((element(el), float(l[30:38]), float(l[38:46]), float(l[46:54])))
    return atoms


def read_ligand(path):
    text = open(path, encoding="utf-8", errors="replace").read()
    if path.endswith(".pdb"):
        return parse_pdb(text)
    mols = parse_sdf(text)
    return mols[0][0] if mols else []


def heavy(atoms):
    return [a for a in atoms if a[0] != "H"]


# ---------------------------------------------------------------- geometry

def centroid(atoms):
    n = len(atoms)
    return tuple(sum(a[i] for a in atoms) / n for i in (1, 2, 3))


def hungarian(cost):
    """Minimal assignment for a square cost matrix (list of lists), O(n^3)."""
    n = len(cost)
    INF = float("inf")
    u, v = [0.0] * (n + 1), [0.0] * (n + 1)
    p, way = [0] * (n + 1), [0] * (n + 1)
    for i in range(1, n + 1):
        p[0] = i
        j0 = 0
        minv = [INF] * (n + 1)
        used = [False] * (n + 1)
        while True:
            used[j0] = True
            i0, delta, j1 = p[j0], INF, 0
            for j in range(1, n + 1):
                if not used[j]:
                    cur = cost[i0 - 1][j - 1] - u[i0] - v[j]
                    if cur < minv[j]:
                        minv[j], way[j] = cur, j0
                    if minv[j] < delta:
                        delta, j1 = minv[j], j
            for j in range(n + 1):
                if used[j]:
                    u[p[j]] += delta
                    v[j] -= delta
                else:
                    minv[j] -= delta
            j0 = j1
            if p[j0] == 0:
                break
        while True:
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            if j0 == 0:
                break
    return sum(cost[p[j] - 1][j - 1] for j in range(1, n + 1))


def rmsd_matched(pose, ref):
    """RMSD without superposition; atoms of the same element are matched optimally."""
    total, n = 0.0, 0
    for el in set(a[0] for a in ref):
        A = [a for a in pose if a[0] == el]
        B = [b for b in ref if b[0] == el]
        if len(A) != len(B):
            return float("inf")
        cost = [[(a[1] - b[1]) ** 2 + (a[2] - b[2]) ** 2 + (a[3] - b[3]) ** 2 for b in B] for a in A]
        total += hungarian(cost)
        n += len(A)
    return math.sqrt(total / n) if n else float("inf")


# ---------------------------------------------------------------- job info

def systems():
    """{code: (receptor file, ligand file, rmsd_max)} from testdata/systems.txt"""
    out = {}
    for line in open(os.path.join(TESTDATA, "systems.txt"), encoding="utf-8"):
        line = line.split("#")[0].split()
        if len(line) == 4:
            out[line[0]] = (os.path.join(TESTDATA, line[1]), os.path.join(TESTDATA, line[2]),
                            float(line[3]))
    return out



# ---------------------------------------------------------------- validation

def validate(files, wu_name):
    m = re.match(r"gnina_(redock|cross|broken)_([a-z0-9]+)_([a-z0-9]+)_", wu_name or "")
    if not m:
        return False, "unknown workunit name %r" % wu_name
    kind, rec_code, lig_code = m.groups()
    if kind == "broken":
        return False, "broken-input job reported success (should have failed on the client)"
    sysmap = systems()
    if rec_code not in sysmap or lig_code not in sysmap:
        return False, "unknown system in %s" % wu_name
    box_lig = heavy(read_ligand(sysmap[rec_code][1]))
    in_lig = heavy(read_ligand(sysmap[lig_code][1]))

    # output files arrive in template order: OUTFILE_0 = docked.sdf, OUTFILE_1 = gnina.log
    sdf = files[0] if files else None
    if not sdf or not os.path.exists(sdf):
        return False, "docked.sdf missing"
    poses = parse_sdf(open(sdf, encoding="utf-8", errors="replace").read())
    if not poses:
        return False, "no poses"

    box_c = centroid(box_lig)
    best = float("inf")
    for i, (atoms, props) in enumerate(poses, 1):
        try:
            aff = float(props["minimizedAffinity"])
            cnn = float(props["CNNscore"])
        except (KeyError, ValueError):
            return False, "pose %d: scores missing" % i
        if not (math.isfinite(aff) and -50 < aff < 100 and 0.0 <= cnn <= 1.0):
            return False, "pose %d: implausible scores (affinity %s, CNNscore %s)" % (i, aff, cnn)
        h = heavy(atoms)
        if len(h) != len(in_lig):
            return False, "pose %d: %d heavy atoms, input ligand has %d" % (i, len(h), len(in_lig))
        if math.dist(centroid(h), box_c) > CENTROID_MAX:
            return False, "pose %d: outside the docking box" % i
        if kind == "redock":
            best = min(best, rmsd_matched(h, in_lig))

    top = poses[0][1]
    info = "%d poses, top affinity %s, top CNNscore %s" % (
        len(poses), top.get("minimizedAffinity"), top.get("CNNscore"))
    if kind == "redock":
        limit = float(RMSD_MAX) if RMSD_MAX else sysmap[lig_code][2]
        info += ", best RMSD %.2f A (limit %.1f)" % (best, limit)
        if best > limit:
            return False, info
    return True, info


def wu_name_from_path(path):
    """upload/xyz/<workunit>_<instance>_r<random>_<index> -> <workunit>"""
    m = re.match(r"(.+)_\d+_r\d+_\d+$", os.path.basename(path))
    return m.group(1) if m else None


def main():
    args = sys.argv[1:]
    if len(args) >= 2 and args[0] == "--wu":           # manual test: --wu NAME files...
        wu_name, files = args[1], args[2:]
    else:                                               # script_validator: output files
        files = args
        wu_name = wu_name_from_path(files[0]) if files else None
        if not wu_name:
            print("gnina_validate: cannot derive the workunit name from %r" % (files,))
            sys.exit(3)  # VAL_RESULT_TRANSIENT_ERROR: configuration problem, retry later
    if not os.path.isfile(os.path.join(TESTDATA, "systems.txt")):
        # server-side problem (app not installed with boinc-app?): never invalidate results for it
        print("gnina_validate: test data not found in %s (set GNINA_TESTDATA)" % TESTDATA)
        sys.exit(3)  # VAL_RESULT_TRANSIENT_ERROR
    ok, info = validate(files, wu_name)
    print("gnina_validate: %s %s: %s" % (wu_name, "VALID" if ok else "INVALID", info))
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
