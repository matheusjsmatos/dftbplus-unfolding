# Implementação Python (ChatGPT)

Código recebido em `DFTBplus_Unfolding.zip` e importado originalmente em `chatgpt-python-unfold/` (tag **chatgpt**).

## O que faz
Abordagem diferente da Fortran: usa as matrizes de espaço real exportadas pelo DFTB+ (`WriteRealHS = Yes`: `hamreal1.dat`,
`overreal.dat`). Reconstrói H(K) e S(K) em qualquer K, resolve o problema generalizado (SciPy/LAPACK, ou Cholesky em NumPy),
calcula as mesmas amplitudes de Fourier e escreve um CSV (energia, peso completo e peso refolded por banda e spin).

Módulos: `dftb_unfold.py` (fluxo), `unfold_io.py` (geometria, matrizes, base), `unfold_fourier.py` (vetores g, Bessel,
harmônicas, integral radial Gauss–Legendre), `spectrum.py` (depósito linear ou gaussiano, NPZ e PDF opcional). A suíte do projeto está em [Testes/](../Testes/README.md).

## Validação da revisão

Os 11 testes Python e as comparações nativas passaram no grafeno 2×2 e no diamante de oito átomos com DFTB+ 26.1 e mio-1-1. O [relatório](../chatgpt-python-unfold/docs/DEPURACAO.md) apresenta as diferenças de energia e peso espectral, e o [HOWTO](HOWTO.md) mostra como repetir os cálculos.

## Resultados da importação original

- Os 8 testes unitários passam.
- Execução no diamante de 8 átomos (`validation/diamond/unfold.json`): 113 vetores g, 45 K distintos, 1472 linhas no CSV.
- Comparação com a versão Fortran: `COMPARACAO.md`.
- Os exemplos `examples/graphene` exigem dftb+, arquivos SK e wfc; não foram reexecutados.

## Vantagem e custo
Aceita K arbitrário sem listar k-pontos antes, mas exige `WriteRealHS` e diagonalização densa em Python (mais lento para sistemas grandes).
