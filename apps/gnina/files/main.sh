#!/bin/sh
# gnina job: dock "ligand" into "receptor", box defined by "box" (crystal ligand of the receptor).
# Extra gnina options (seed, CNN model, scoring mode) come from the workunit's command line,
# which docker_wrapper passes in $ARGS. The GPU assigned by the BOINC client is the only one
# CUDA sees (CUDA_VISIBLE_DEVICES, set by the docker_wrapper), so gnina's default device 0 is
# the right one; do not pass --device.
# Outputs: docked.sdf (poses with scores, gzip-compressed), gnina.log
set -u
# timestamps for the steps in the task's stderr; UTC, as the container has no time zone
# (the docker_wrapper's own lines before and after use the host's local time)
ts() { date -u '+%Y-%m-%d %H:%M:%S UTC'; }
echo "$(ts) main.sh @APP_VERSION@ started" >&2

# gnina picks the file format from the extension: detect PDB vs. SDF by content
fmt() { if grep -qE '^(ATOM|HETATM)' "$1"; then echo pdb; else echo sdf; fi; }
cp receptor rec.pdb
cp ligand "lig.$(fmt ligand)"
cp box "box.$(fmt box)"

echo "$(ts) gnina started, args: ${ARGS:-}" >&2
echo "gnina job, args: ${ARGS:-}" > gnina.log
# shellcheck disable=SC2086
gnina -r rec.pdb -l lig.* --autobox_ligand box.* -o docked.sdf ${ARGS:-} >> gnina.log 2>&1
rc=$?
echo "$(ts) gnina finished, exit code $rc; gnina.log follows" >&2
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
# compress the poses for the upload (SDF text shrinks to a third or less); the file keeps its
# name, the validator and "boinc-app results" recognize gzip by its magic bytes
if gzip -9 -n -c docked.sdf > docked.sdf.gz; then
    mv docked.sdf.gz docked.sdf
else
    rm -f docked.sdf.gz
    echo "gzip failed, uploading docked.sdf uncompressed" >&2
fi
echo "$(ts) main.sh finished" >&2
exit 0
