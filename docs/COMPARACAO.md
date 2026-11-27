# Comparação entre as duas implementações

A suíte em [Testes/](../Testes/README.md) reproduz a comparação com o DFTB+ 26.1 e os parâmetros mio-1-1. No diamante, os pesos refoldados integrados de −25 a +15 eV, com depósito linear, foram **323,17491 no Fortran** e **323,17275 no Python**. O excesso de 1,7–3,7% observado no ensaio original não se reproduziu com as entradas atuais. A origem da diferença histórica segue em aberto. [Resultados por material e por alargamento](../chatgpt-python-unfold/docs/DEPURACAO.md)

## Ensaio original do diamante

Mesmo sistema: diamante, supercélula de 8 átomos (matriz [[-1,1,1],[1,-1,1],[1,1,-1]]), mio-1-1, caminho L–Γ–X em 46 pontos,
corte de 20 Ry (113 vetores g nos dois casos), energia relativa a −1,0731 eV.

| Quantidade | Fortran | Python | Diferença |
|---|---|---|---|
| Peso desdobrado (g = 0), total na janela | 39,4585 | 39,4595 | 5e-4 (relativa, L2) |
| Peso refolded, total na janela de −25 a 15 eV | ≈ 329 | ≈ 323 | Fortran 1,7–3,7% maior, dependendo de q |

Os dois concordam no unfolding completo. No refolded há uma discrepância de poucos por cento, **cuja causa não foi
identificada**. Os dois usam o mesmo conjunto de g e a mesma janela. Hipóteses não testadas:
tratamento de caudas em |q+g| grande da transformada radial (truncamento no raio de corte) e arredondamento no binning.
Ambos validam contra a célula primitiva no caso Fortran; o Python não foi comparado diretamente com bandas primitivas aqui.
Antes de usar números absolutos do refolded em publicação, investigue esse ponto.
