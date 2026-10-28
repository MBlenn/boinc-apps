#!/bin/sh
# gnina job: dock "ligand" into "receptor", box defined by "box" (crystal ligand of the receptor).
# Extra gnina options (seed, CNN model, scoring mode) come from the workunit's command line,
# which docker_wrapper passes in $ARGS (BOINC also adds "--device N" for the assigned GPU).
# Outputs: docked.sdf (poses with scores), gnina.log
set -u

# gnina picks the file format from the extension: detect PDB vs. SDF by content
fmt() { if grep -qE '^(ATOM|HETATM)' "$1"; then echo pdb; else echo sdf; fi; }
cp receptor rec.pdb
cp ligand "lig.$(fmt ligand)"
cp box "box.$(fmt box)"

echo "gnina job, args: ${ARGS:-}" > gnina.log
# shellcheck disable=SC2086
gnina -r rec.pdb -l lig.* --autobox_ligand box.* -o docked.sdf ${ARGS:-} >> gnina.log 2>&1
rc=$?
cat gnina.log >&2

if [ $rc -ne 0 ]; then
    echo "gnina failed with exit code $rc" >&2
    exit $rc
fi
# gnina exits 0 even if it could not read the ligand: treat "no poses" as an error
if ! grep -q '\$\$\$\$' docked.sdf 2>/dev/null; then
    echo "gnina produced no poses (unreadable input?)" >&2
    exit 2
fi
exit 0
