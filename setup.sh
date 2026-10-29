#!/usr/bin/env bash
# One-command bootstrap for a fresh checkout on any machine.
#
# A Python venv is NOT portable (baked-in absolute interpreter path,
# compiled C-extension wheels tied to this machine's OS/arch) -- copying
# .venv itself across machines does not work. This script instead
# recreates it from scratch, which is fast and always works. Re-run it
# after every fresh `git clone`/copy, not just once.
set -euo pipefail
cd "$(dirname "$0")"

rm -rf .venv
python3 -m venv .venv
.venv/bin/pip install --upgrade pip
.venv/bin/pip install -r requirements.txt
.venv/bin/pip install -e .

# Registers this venv as a named Jupyter kernel (~/.local/share/jupyter/
# kernels/ion_track_mc/kernel.json) so the notebooks -- whose kernelspec
# metadata names "ion_track_mc" -- find the right interpreter regardless
# of which Jupyter install (this venv's, or a system-wide one) opens them.
.venv/bin/python -m ipykernel install --user --name=ion_track_mc \
    --display-name "ion_track_mc (.venv)"

echo
echo "Setup complete."
echo "  Activate:       source .venv/bin/activate"
echo "  Run tests:      .venv/bin/pytest -q"
echo "  Open notebooks: .venv/bin/jupyter lab"
