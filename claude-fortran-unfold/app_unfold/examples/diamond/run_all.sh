#!/bin/bash
# End-to-end example: an 8-atom cubic diamond cell is unfolded onto the fcc primitive cell.
#
# Usage:  SKDIR=/path/to/mio-1-1 BIN=/path/to/build/app ./run_all.sh
#
# The Slater-Koster files and the STO basis (wfc_C.hsd) must belong to the same parameter set
# (here mio-1-1). Using the STOs of another set makes the sum rule check fail (see README).
set -e
: "${SKDIR:?set SKDIR to the directory with the mio-1-1 SK-files}"
: "${BIN:?set BIN to the build/app directory of DFTB+}"
cd "$(dirname "$0")"
rm -rf run1 run2; mkdir run1 run2
sed "s#@SKDIR@#${SKDIR%/}#" dftb_in.hsd.template > dftb_in.hsd

# 1) SCC calculation of the supercell (charges) -------------------------------------------------
cp wfc_C.hsd dftb_in.hsd run1/
cp unfold_list.hsd run1/unfold_in.hsd
( cd run1 && $BIN/dftb+/dftb+ > scc.out 2>&1 )

# 2) which k-points are needed for the requested q-path? ----------------------------------------
( cd run1 && $BIN/unfold/unfold > list.out 2>&1 )

# 3) non-SCC run at exactly these k-points, with eigenvectors and band.out ----------------------
cp run1/charges.bin run1/unfold.kpoints.hsd run2/
sed -e '/KPointsAndWeights = SupercellFolding {/,/^  }/d' dftb_in.hsd | python3 -c "
import sys
s=sys.stdin.read()
s=s.replace('  SCCTolerance = 1e-8\n','  SCCTolerance = 1e-8\n  MaxSCCIterations = 1\n  ReadInitialCharges = Yes\n  <<+ \"unfold.kpoints.hsd\"\n')
print(s)" > run2/dftb_in.hsd
cp wfc_C.hsd run2/
( cd run2 && $BIN/dftb+/dftb+ > scc2.out 2>&1 )

# 4) unfolding ------------------------------------------------------------------------------------
EF=$(grep "Fermi level" run1/detailed.out | awk '{print $5}')
echo "Fermi level (eV) of the SCC run: $EF"
sed "s/@FERMI@/$EF/" unfold_run.hsd > run2/unfold_in.hsd
( cd run2 && $BIN/unfold/unfold > unfold.out 2>&1 )
tail -8 run2/unfold.out
echo "Results: run2/unfold.unfoldedBands and run2/unfold.refoldedBands"
