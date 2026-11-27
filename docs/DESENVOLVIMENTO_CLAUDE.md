# Desenvolvimento: implementação Fortran (Claude)

## Decisão de arquitetura
Pós-processador autônomo ("caminho A"): não altera o núcleo do DFTB+. Lê saídas existentes e fica em `app/unfold`,
compilado pelo CMake do DFTB+ (reusa o parser HSD, a leitura de geometria e as constantes).

## Procedimento
1. **Análise** do repositório DFTB+ (commit `fedabbc`) e do artigo: o que o DFTB+ já escreve (`eigenvec.bin`, `band.out`,
   `detailed.xml`) e como o Waveplot define a base STO e as harmônicas reais. Conclusão: viável sem mudar o núcleo.
2. **Módulo radial** (`radial.F90`): STO a partir do `wfc`, transformada F_l(q) tabelada (dq = 0,01, Simpson em 6001 pontos),
   interpolação linear, Bessel esféricas até l = 3, harmônicas reais na ordem do Waveplot.
3. **Entrada** (`initunfold.F90`): arquivo `unfold_in.hsd` com ações `Unfold`, `ListKPoints` e `CheckBasis`.
   Gera o conjunto de vetores g e lê o caminho q.
4. **Programa principal** (`unfold.F90`): mapeia cada q+g para o k-ponto do arquivo (módulo vetor recíproco), lê autovetores
   por acesso direto, calcula os pesos, deposita em malha de energia e escreve os arquivos no formato do SIESTA.
5. **Validação** (diamante com supercélula de 8 átomos e grafeno 3×3), comparando com a célula primitiva.
6. **Ajustes dos testes:** gen em fracionárias (`F`); sinal de K verificado em k-pontos genéricos; mesma base STO e SK.

## Fluxo de uso
1. `ListKPoints` escreve `<Output>.kpoints.hsd` com os k-pontos necessários para o DFTB+.
2. Rodar o DFTB+ (não-SCC, cargas lidas) com esses k-pontos, gravando `eigenvec.bin`.
3. `Unfold` escreve `.unfoldedBands` (g = 0) e `.refoldedBands` (soma sobre g).

## Resultados de validação
- Diamante (supercélula 8 átomos, mio-1-1): 100% do peso refolded a menos de 0,15 eV das bandas primitivas; desvio de centróide 0,7 meV.
- Grafeno 3×3: 100% do peso nas bandas primitivas; desvio < 0,1 meV; regra de soma 0,95–1,00.
- Figura: `docs/figures/diamond_unfold.png`. Scripts: `validation/`.

## Armadilhas encontradas
`getChild(..., modifier=)` para entradas com unidade; `invert33` sem determinante não inicializado; `<<<` em linha própria no HSD;
conjunto wfc igual ao dos SK.
