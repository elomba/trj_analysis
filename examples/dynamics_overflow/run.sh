#!/bin/bash
set -euo pipefail

example_dir=$(cd "$(dirname "$0")" && pwd)
binary=${1:-"$example_dir/../../bin/trj_analysis"}
work_dir=${2:-"$example_dir/work"}

python3 "$example_dir/generate.py"
mkdir -p "$work_dir"
ln -sfn "$example_dir/run.nc" "$work_dir/run.nc"
cp "$example_dir/input.nml" "$work_dir/input.nml"
(cd "$work_dir" && "$binary" input.nml > run.log 2>&1)
python3 "$example_dir/check.py" "$work_dir" --require-overflow
