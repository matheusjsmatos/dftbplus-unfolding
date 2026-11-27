#!/usr/bin/env bash
# Exercise the Fortran examples shipped with the manual, including their plotting command.
set -euo pipefail
if [[ $# -ne 3 ]]; then
  echo "Usage: bash Testes/check_howto.sh BUILD_APP MIO_ROOT OUTPUT" >&2
  exit 2
fi

repo=$(cd -- "$(dirname -- "$0")/.." && pwd)
BIN=$(cd -- "$1" && pwd)
SKDIR=$(cd -- "$2/skfiles" && pwd)
export BIN SKDIR
mkdir -p -- "$3"
out=$(cd -- "$3" && pwd)
work=$(mktemp -d "$out/run.XXXXXX")
echo "HOWTO results: $work"

diamond="$repo/claude-fortran-unfold/app_unfold/examples/diamond"
plot="$repo/claude-fortran-unfold/tools/plot_unfolded.py"
cp -R "$diamond" "$work/diamond"
bash "$work/diamond/run_all.sh"

cp -R "$repo/validation/graphene" "$work/graphene"
cp "$diamond/wfc_C.hsd" "$work/graphene/"
(cd "$work/graphene" && python gen.py && bash mk.sh)

for example in diamond graphene; do
  map="$work/$example/run2/unfold.refoldedBands"
  test -s "$map"
  figure="$work/$example/refolded.png"
  python "$plot" "$map" -o "$figure"
  test -s "$figure"
done
