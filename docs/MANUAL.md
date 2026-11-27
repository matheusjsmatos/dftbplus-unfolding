# Manual de referência

Este manual descreve os dois pós-processadores de bandas eletrônicas do DFTB+ em base LCAO. Para executar e plotar um exemplo, consulte o [HOWTO](HOWTO.md). As equações de projeção e as convenções de fase estão em [Método](METODO.md).

| | **unfold** (Fortran, Claude) | **dftb_unfold.py** (Python, ChatGPT) |
|---|---|---|
| Pasta | `claude-fortran-unfold/` | `chatgpt-python-unfold/` |
| O que lê do DFTB+ | `eigenvec.bin`, `band.out`, `detailed.xml` | `hamreal1.dat`, `overreal.dat` (e `hamreal2.dat`), geometria |
| Pontos de Bloch disponíveis | Lê os autovetores dos pontos calculados pelo DFTB+. `ListKPoints` gera a lista necessária para o caminho escolhido. | Reconstrói e diagonaliza $H(\mathbf K)$ e $S(\mathbf K)$ nos pontos necessários a partir dos blocos em espaço real. |
| Configuração | `unfold_in.hsd` (HSD) | `unfold.json` |
| Saída | `.unfoldedBands`, `.refoldedBands` (formato SIESTA) | CSV, e `spectrum.py` gera NPZ/PDF |
| Custo | uma execução do DFTB+ nos k-pontos listados, depois leitura | uma diagonalização densa por (K, spin) |

Use em ambos os programas os orbitais radiais `wfc` correspondentes à parametrização Slater–Koster do cálculo. As matrizes SK e os coeficientes radiais desempenham papéis diferentes na projeção. [Método](METODO.md) · [Diagnóstico de base no HOWTO](HOWTO.md)

## 1. Requisitos

**Fortran.** A compilação e as comparações no Linux foram verificadas com o código-fonte do DFTB+ **26.1**. A integração original também foi testada no commit `fedabbcb993ef26967141f275bb7d31a984e4bdc`. Copie `claude-fortran-unfold/app_unfold` para `app/unfold` na árvore do DFTB+ e aplique `claude-fortran-unfold/unfold-cmake.patch`. A compilação usa o compilador Fortran, o CMake, o pré-processador `fypp` e as bibliotecas BLAS/LAPACK do DFTB+. O executável fica em `<build>/app/unfold/unfold`. A versão 25.1 tem assinaturas incompatíveis nas chamadas `printDftbHeader` e `warnUnprocessedNodes`. [Comandos de compilação testados](../.github/workflows/validar-unfolding.yml)

**Python.** Python 3.9 ou superior e `numpy>=1.23`. O SciPy, quando instalado, fornece a diagonalização generalizada por LAPACK. O caminho com NumPy usa a fatoração de Cholesky. `spectrum.py --pdf` requer Pillow e ReportLab. O script de figuras da implementação Fortran requer NumPy e Matplotlib.

**Entradas comuns.**

- Cálculo DFTB+ da supercélula com a parametrização SK escolhida.
- Arquivo `wfc.<conjunto>.hsd` do mesmo conjunto SK, em formato Waveplot (Slater, com coeficientes e expoentes e `Cutoff`). Até $l=3$.

## 2. Fortran: `unfold`

### 2.1 Execução

```bash
unfold            # lê ./unfold_in.hsd
```

O programa lê sempre `unfold_in.hsd` no diretório corrente e escreve o HSD processado em `unfold_pin.hsd`.

### 2.2 Ações (`Action`)

| Valor | O que faz |
|---|---|
| `ListKPoints` | Lê `detailed.xml` e o caminho `QPoints` e escreve `<Output>.kpoints.hsd`, com os k-pontos que o DFTB+ precisa calcular. |
| `Unfold` (padrão) | Lê `detailed.xml`, `eigenvec.bin`, `band.out` e a base, e escreve `<Output>.unfoldedBands` e (se há célula de refolding) `<Output>.refoldedBands`. |
| `CheckBasis` | Verifica a regra de soma nos primeiros 12 pontos k disponíveis, ou em todos se houver menos, usando os autovetores do primeiro canal de spin. Imprime mínimo, média e máximo sobre os estados. |

### 2.3 Chaves de `unfold_in.hsd`

| Chave | Padrão | Significado |
|---|---|---|
| `Action` | `Unfold` | Veja 2.2. |
| `DetailedXML` | `detailed.xml` | Geometria, pontos k, dimensões e identificador da execução. Para `Unfold`, use o XML da mesma execução de bandas que produziu `eigenvec.bin` e `band.out` (`Options { WriteDetailedXML = Yes }`). |
| `EigenvecBin` | `eigenvec.bin` | Autovetores (acesso *stream*, padrão do DFTB+). Precisa de `WriteEigenvectors = Yes` em `Analysis`. |
| `BandOut` | `band.out` | Autovalores, em eV. |
| `Output` | `unfold` | Prefixo dos arquivos de saída. |
| `Basis` | — | STOs no formato Waveplot, por espécie. O bloco pode ser incluído com `Basis { <<+ "wfc_C.hsd" }`. `ListKPoints` dispensa esse bloco. |
| `NormalizeOrbitals` | `Yes` | Renormaliza as radiais na projeção. O comparador nativo usa `No` para aplicar os mesmos coeficientes Waveplot que o Python. |
| `QPoints` | — | Caminho em $\mathbf q$ no estilo `Klines`: `1 q1 q2 q3` para o ponto inicial, `n q1 q2 q3` para $n$ pontos equidistantes terminando no ponto dado. Coordenadas fracionárias na base recíproca da célula de refolding (ou da supercélula, se não houver). Vários caminhos: blocos `Path { ... }`, cada um gera arquivo próprio. Fora da primeira zona dá o espectro na zona estendida. |
| `RefoldingLatticeVectors` | — | Vetores da célula de refolding (linhas), com unidade, por exemplo `[Angstrom]`. Sem esta chave só se calcula o espectro desdobrado completo. |
| `RefoldingGCutoff` | 10 Hartree | Inclui vetores com $\lvert\mathbf g\rvert^2/2\le c$ em unidades atômicas. `RefoldingGCutoff [Hartree] = 12.5` corresponde a 25 Ry. |
| `Energy { Window }` | −20 a 10 eV | Limites inferior e superior de $E-E_{\mathrm{ref}}$. Para valores em eV, escreva `Window [eV] = { -20 10 }`. |
| `Energy { NrOfIntervals }` | 300 | Intervalos da malha de energia (a saída tem N+1 pontos). |
| `Energy { FermiLevel }` | 0 | Referência subtraída dos autovalores. Use `FermiLevel [eV]` com o Fermi do SCC para representar $E-E_F$. |
| `Energy { Broadening }` | 0 | `Broadening [eV]` é o desvio-padrão σ da gaussiana. O valor zero reparte o peso linearmente entre dois pontos vizinhos da malha. |
| `KPointSign` | −1 | Sinal entre o $\mathbf k$ do arquivo e o $\mathbf K$ físico (seção 4.2 de `METODO.md`). Não altere sem motivo. |
| `SumRuleCheck` | `Yes` | Em `Unfold`, verifica o primeiro ponto k utilizado, no primeiro canal de spin. `CheckBasis` amplia a amostra de k. |
| `SumRuleQMax` | 8 bohr⁻¹ | Faixa de $q$ na soma da regra. |

### 2.4 Formato de saída

O cabeçalho traz `nq  ne+1  Emin  Emax  EFermi`, com energias em eV. Para cada $\mathbf q$ há uma linha com três componentes cartesianas em bohr⁻¹ e o índice do segmento. Seguem-se `ne+1` valores, um por linha, para $E_j=E_{\min}+j\,\Delta E$, com $\Delta E=(E_{\max}-E_{\min})/\texttt{ne}$. As energias do eixo já incluem a subtração de `FermiLevel`.

O mapa Fortran armazena **peso por ponto da malha de energia**. Divida por $\Delta E$ em eV para compará-lo à densidade espectral por eV do mapa Python. A soma das projeções em ondas planas se aproxima de um por estado conforme convergem o corte e a representação radial da base. Cada peso individual e cada integral em uma janela finita dependem da projeção e da janela escolhidas.

Cálculos com dois canais geram sufixos `.spinUp` e `.spinDn`. Se houver vários caminhos, cada um recebe também `.path1`, `.path2` etc. Por exemplo, `unfold.refoldedBands.spinUp.path1` contém o primeiro canal e o primeiro caminho.

A seção 4 descreve a leitura desses arquivos para visualização.

### 2.5 Mensagens e verificações

A saída padrão, que pode ser salva com `unfold > unfold.out`, contém o número de pontos $\mathbf q$, o número de vetores $\mathbf g$, os pontos k no arquivo, o volume, as normas das STOs e, quando solicitada, a regra de soma. `Missing k-points` indica que faltam autovetores para os pontos pedidos. Use `ListKPoints` e execute novamente o cálculo de bandas. `The eigenvector file has not the expected size` indica um `eigenvec.bin` incompatível com as dimensões ou o formato *stream* esperado. Verifique a origem comum dos três arquivos e `WriteEigenvectors = Yes`.

## 3. Python: `dftb_unfold.py` e `spectrum.py`

### 3.1 Linha de comando

```
python dftb_unfold.py CONFIG.json [-o saida.csv] [--plan]
python spectrum.py saida.csv --range EMIN EMAX [--sigma 0.06] [--step 0.01]
                   [--spin 1] [--weight full|refold] [-o mapa.npz] [--pdf mapa.pdf]
```

`--plan` conta pontos $\mathbf q$, vetores $\mathbf g$, K distintos e se a referência é comensurável, sem diagonalizar. Rode-o antes de uma execução grande.

### 3.2 Chaves de `unfold.json`

`geometry`, `matrix_dir` e os nomes em `wfc_files` são relativos ao diretório do JSON. O caminho passado a `-o` na linha de comando é relativo ao diretório de execução.

| Chave | Obrigatória | Significado |
|---|---|---|
| `geometry` | sim | `detailed.xml` ou GEN periódico (`S` cartesiano em Å, ou `F` fracionário). Mesma ordem de átomos da matriz. |
| `matrix_dir` | não (`.`) | Pasta com `hamreal1.dat`, `overreal.dat` e, para spin colinear, `hamreal2.dat`. |
| `wfc_files` | sim | Lista de arquivos Waveplot (um por espécie ou combinados). |
| `supercell_matrix` | uma das duas | Matriz inteira $M$ com $A_{\rm SC}=M A_{\rm ref}$, vetores nas linhas. |
| `reference_cell_angstrom` | uma das duas | Três vetores da célula de referência, em Å, explícitos. |
| `path_frac` | sim | Vértices do caminho, em coordenadas fracionárias recíprocas da **referência**. |
| `intervals` | sim | Intervalos de cada segmento. Um por segmento, vértices compartilhados aparecem uma vez. |
| `labels` | não | Rótulo de cada vértice. O leitor verifica a contagem. O PDF de `spectrum.py` exibe a distância acumulada em Å⁻¹ no eixo horizontal. |
| `gcut_ry` | não (0) | Corte nos vetores de refolding, $\lvert\mathbf g\rvert^2<c_{\rm Ry}$ em unidades atômicas. Com corte zero, só entra $\mathbf g=0$ e o CSV grava `weight_refold = weight_full`. |
| `reference_ev` | não (0) | Energia de referência subtraída dos autovalores. Use o Fermi do SCC, em eV, para representar $E-E_F$. Com zero, o grafeno arquivado conserva o zero de energia das saídas DFTB+. |
| `window_ev` | não | Faixa `[min, max]` relativa a `reference_ev`. Restringe os estados gravados no CSV. Para comparar espectros com alargamento, grave todas as bandas e aplique a janela em `spectrum.py`. |
| `q_batch` | não (64) | Quantos vetores são projetados de uma vez. Reduz memória. |

### 3.3 Colunas do CSV

| Coluna | Significado |
|---|---|
| `point`, `distance_inv_ang`, `q1,q2,q3` | Índice (a partir de 1), distância acumulada em Å⁻¹ e coordenadas fracionárias da referência. |
| `K1_sc,K2_sc,K3_sc` | $\mathbf K$ físico da supercélula, em coordenadas recíprocas fracionárias reduzidas a $[-0{,}5,0{,}5)$. |
| `spin`, `band` | Índices a partir de 1. |
| `energy_ev`, `energy_rel_ev` | Autovalor e autovalor menos `reference_ev`, em eV. |
| `weight_full` | Peso em $\mathbf q$ (unfolding completo, $\mathbf g=0$). |
| `weight_refold` | Soma dos pesos em $\mathbf q+\mathbf g$ para $\lvert\mathbf g\rvert^2<$ `gcut_ry`. |

### 3.4 `spectrum.py`

Com `--sigma` positivo, o programa deposita cada linha em uma gaussiana de área unitária e desvio-padrão em eV. Com `--sigma 0`, usa depósito linear. O NPZ guarda `distance_inv_ang`, `energy_rel_ev` e a matriz `spectral_weight_per_ev`, com dimensões (pontos do caminho, energias). A integral em uma janela discreta é a soma da densidade multiplicada pelo espaçamento real do eixo. `--pdf` desenha o mapa. Para bandas degeneradas, compare a soma no subespaço.

Se a referência for incomensurável, vários $\mathbf K$ podem contribuir ao mesmo ponto do caminho. O CSV conserva uma linha por estado e grupo de contribuições com o mesmo $\mathbf K$. Some os pesos na energia de cada linha, como faz `spectrum.py`.

### 3.5 Limites do código Python

A implementação lê blocos reais $H(R)$ e $S(R)$ do DFTB+ periódico. As matrizes de Bloch podem ser complexas fora de Γ. A base tem uma camada radial por $l$, e o tratamento de spin abrange o caso escalar e o colinear. SOC, spin não colinear, contatos e fônons ficam fora desse escopo. A diagonalização custa $\mathcal O(N_{\rm orb}^3)$ por par (K, spin), com memória $\mathcal O(N_{\rm orb}^2)$ para as matrizes densas.

## 4. Visualização

**Saída Fortran.** `python claude-fortran-unfold/tools/plot_unfolded.py ARQUIVO [-o fig.png] [--vmax FRAC] [--ref band.out FERMI]` desenha o peso em escala de cinza. O eixo horizontal é o índice do ponto no caminho. `--ref` superpõe em vermelho as bandas da primitiva, após subtrair `FERMI` em eV. Use a mesma trajetória, amostragem e referência energética. Para spin colinear, a referência deve conter apenas o canal desejado, pois o script lê todos os blocos `KPT` na ordem do arquivo.

**Saída Python.** `spectrum.py --pdf` usa a distância acumulada no caminho. Para ver as bandas π do grafeno plano, plote `--weight refold` e use um corte que inclua vetores com componente fora do plano. No espectro completo ao longo de $q_z=0$, a paridade ímpar dos orbitais $p_z$ pode anular a projeção.

## 5. Unidades e convenções

| Grandeza | Fortran | Python |
|---|---|---|
| Energias | Entradas HSD com `[eV]` para valores em eV. Sem modificador, a unidade de entrada é Hartree. Saída em eV. | Configuração e saídas de energia em eV. |
| $\mathbf q$ na entrada | fracionário (base recíproca da célula de refolding) | fracionário (idem) |
| $\mathbf q$ na saída | cartesiano, bohr⁻¹ | fracionário no CSV, distância em Å⁻¹ |
| Corte em $\mathbf g$ | Hartree, $\lvert\mathbf g\rvert^2/2\le c_{\rm Ha}$ | Ry, $\lvert\mathbf g\rvert^2<c_{\rm Ry}$, com $c_{\rm Ry}=2c_{\rm Ha}$. Uma casca exatamente no limiar pode ser tratada de modo diferente. |
| Peso espectral | Peso adimensional por ponto da malha. | Peso adimensional por estado no CSV. Densidade em eV⁻¹ no NPZ. |

## 6. Solução de problemas

| Sintoma | Causa provável | O que fazer |
|---|---|---|
| Regra de soma longe de 1 (por exemplo 0,4–3,8) | `wfc` de outro conjunto que os SK, ordem de orbitais errada, ou `KPointSign` errado | HOWTO 6 |
| A regra de soma se aproxima de 1 lentamente | Corte insuficiente ou diferença entre a sobreposição das STOs e a matriz $S$ do conjunto SK | Aumente `SumRuleQMax` e acompanhe a convergência antes de atribuir o desvio ao ajuste radial. |
| `Missing k-points` (Fortran) | `eigenvec.bin` não calculado nos k-pontos de `ListKPoints` | Refaça o passo 3 do HOWTO 1 |
| `geometry and matrix have different numbers of atoms` (Python) | Geometria e matrizes de cálculos diferentes | Use o `detailed.xml` do mesmo cálculo |
| `commensurate: False` (Python, `--plan`) | A relação entre os vetores da célula simulada e os da referência não é uma matriz inteira | Confira os vetores e as unidades. Se a referência incomensurável for intencional, use o número de K distintos mostrado por `--plan` para estimar o custo. |
| Pesos fora das bandas do cristal perfeito | Defeito ou desordem, ou erro de referência | HOWTO 7 |
| Bandas π ausentes do mapa `full` em grafeno plano | A projeção em $q_z=0$ pode se anular pela paridade de $p_z$ | Veja a seção 4 e o [HOWTO](HOWTO.md) para gerar o mapa refoldado. |

## 7. Testes de referência

Na raiz do repositório, execute `python Testes/run.py` para a suíte offline com NumPy. Ela reproduz o grafeno 2×2 a partir das matrizes e das bandas DFTB+ arquivadas. A modalidade `native` usa os executáveis `dftb+` e `unfold` e os parâmetros mio-1-1 para recalcular grafeno e diamante e comparar os dois pós-processadores. [Como executar os testes](../Testes/README.md) · [Critérios, parâmetros e resultados numéricos](../chatgpt-python-unfold/docs/DEPURACAO.md)
