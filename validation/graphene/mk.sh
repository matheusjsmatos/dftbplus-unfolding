#!/bin/bash
set -e
B="${BIN:?defina BIN com o diretorio build/app do DFTB+}"
cd "$(dirname "$0")"
rm -rf run1 run2 prim; mkdir run1 run2 prim prim/bands
SK="${SKDIR:?defina SKDIR com o diretorio dos SK mio-1-1}/"
hsd() { # genfile kpoints-block extra
cat <<EOT
Geometry = GenFormat {
  <<< "$1"
}
Hamiltonian = DFTB {
  SCC = Yes
  SCCTolerance = 1e-8
  MaxAngularMomentum { C = "p" }
  SlaterKosterFiles = Type2FileNames { Prefix = "$SK"
    Separator = "-"
    Suffix = ".skf" }
  $3
  KPointsAndWeights = $2
}
Options { WriteChargesAsText = No
  WriteDetailedXML = Yes }
Analysis { WriteEigenvectors = Yes }
ParserOptions { ParserVersion = 14 }
EOT
}
cp sc.gen wfc_C.hsd run1/; hsd sc.gen "SupercellFolding {
 4 0 0
 0 4 0
 0 0 1
 0.5 0.5 0.0 }" "" > run1/dftb_in.hsd
cp unfold_list.hsd run1/unfold_in.hsd
( cd run1 && $B/dftb+/dftb+ > scc.out 2>&1 && $B/unfold/unfold > list.out 2>&1 )
cp run1/charges.bin run1/unfold.kpoints.hsd run2/; cp sc.gen wfc_C.hsd run2/
hsd sc.gen "{}" "MaxSCCIterations = 1
  ReadInitialCharges = Yes
  <<+ \"unfold.kpoints.hsd\"" | sed '/KPointsAndWeights = {}/d' > run2/dftb_in.hsd
( cd run2 && $B/dftb+/dftb+ > scc2.out 2>&1 )
EF=$(grep "Fermi level" run1/detailed.out | awk '{print $5}'); echo "Fermi $EF"
sed "s/@FERMI@/$EF/" unfold_run.hsd > run2/unfold_in.hsd
( cd run2 && $B/unfold/unfold > unfold.out 2>&1 ) || true
# primitive reference
cp prim.gen prim/; hsd prim.gen "SupercellFolding {
 12 0 0
 0 12 0
 0 0 1
 0.5 0.5 0.0 }" "" > prim/dftb_in.hsd
( cd prim && $B/dftb+/dftb+ > scc.out 2>&1 ) 
cp prim/charges.bin prim/prim.gen prim/bands/
hsd prim.gen "Klines { 1 0.3333333333 0.3333333333 0 
 30 0 0 0
 20 0.5 0 0
 25 0.3333333333 0.3333333333 0 }" "MaxSCCIterations = 1
  ReadInitialCharges = Yes" > prim/bands/dftb_in.hsd
( cd prim/bands && $B/dftb+/dftb+ > b.out 2>&1 )
tail -9 run2/unfold.out
echo "primitive Fermi: $(grep 'Fermi level' prim/detailed.out | awk '{print $5}')"
