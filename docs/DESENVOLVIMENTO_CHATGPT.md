# Implementação Python (ChatGPT)

Código recebido em `DFTBplus_Unfolding.zip` e incluído sem modificações em `chatgpt-python-unfold/` (tag **chatgpt**).

## O que faz
Abordagem diferente da Fortran: usa as matrizes de espaço real exportadas pelo DFTB+ (`WriteRealHS = Yes`: `hamreal1.dat`,
`overreal.dat`). Reconstrói H(K) e S(K) em qualquer K, resolve o problema generalizado (SciPy/LAPACK, ou Cholesky em NumPy),
calcula as mesmas amplitudes de Fourier e escreve um CSV (energia, peso completo e peso refolded por banda e spin).

Módulos: `dftb_unfold.py` (fluxo), `unfold_io.py` (geometria, matrizes, base), `unfold_fourier.py` (vetores g, Bessel,
harmônicas, integral radial Gauss–Legendre), `spectrum.py` (alargamento gaussiano, NPZ, PDF opcional), `test_unfold.py` (8 testes).

## O que foi verificado aqui
- Os 8 testes unitários passam.
- Execução no diamante de 8 átomos (`validation/diamond/unfold.json`): 113 vetores g, 45 K distintos, 1472 linhas no CSV.
- Comparação com a versão Fortran: `COMPARACAO.md`.
- Os exemplos `examples/graphene` exigem dftb+, arquivos SK e wfc; não foram reexecutados.

## Vantagem e custo
Aceita K arbitrário sem listar k-pontos antes, mas exige `WriteRealHS` e diagonalização densa em Python (mais lento para sistemas grandes).
