# Validação

- `diamond/compare.py`: compara `unfold.refoldedBands` (Fortran) com `band.out` da célula primitiva (uso: `python3 compare.py <nivel_de_fermi>`).
  O fluxo completo que gera esses arquivos é `claude-fortran-unfold/app_unfold/examples/diamond/run_all.sh`.
- `diamond/dftb_in_writeRealHS.hsd` e `diamond/unfold.json`: entradas usadas para rodar a versão Python (ChatGPT) no mesmo diamante.
- `graphene/`: supercélula 3×3 de grafeno (`gen.py` e `mk.sh` geram a geometria e rodam). Defina `BIN` (diretório `build/app` do DFTB+) e
  `SKDIR` (SK mio-1-1) antes de rodar `mk.sh`.
- Os arquivos Slater–Koster (mio-1-1) e `wfc` correspondentes não são redistribuídos aqui.
