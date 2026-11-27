# Manual de referência

Este manual descreve as entradas, as saídas e as opções dos dois programas do repositório. Para tarefas passo a passo, veja `HOWTO.md`. A teoria está em `METODO.md`.

| | **unfold** (Fortran, Claude) | **dftb_unfold.py** (Python, ChatGPT) |
|---|---|---|
| Pasta | `claude-fortran-unfold/` | `chatgpt-python-unfold/` |
| O que lê do DFTB+ | `eigenvec.bin`, `band.out`, `detailed.xml` | `hamreal1.dat`, `overreal.dat` (e `hamreal2.dat`), geometria |
| K arbitrário | Não. Só os k-pontos que o DFTB+ calculou | Sim. Diagonaliza $H(\mathbf K)$ e $S(\mathbf K)$ |
| Configuração | `unfold_in.hsd` (HSD) | `unfold.json` |
| Saída | `.unfoldedBands`, `.refoldedBands` (formato SIESTA) | CSV, e `spectrum.py` gera NPZ/PDF |
| Custo | uma execução do DFTB+ nos k-pontos listados, depois leitura | uma diagonalização densa por (K, spin) |

Em ambos, **o conjunto de STOs do `wfc` deve ser o mesmo dos arquivos Slater–Koster usados no cálculo** (veja a seção 5 de `METODO.md` e o HOWTO 6).

## 1. Requisitos

**Fortran.** Árvore do DFTB+ no commit `fedabbcb993ef26967141f275bb7d31a984e4bdc` (testado), com `app/unfold` copiada e `add_subdirectory(unfold)` em `app/CMakeLists.txt` (o arquivo `unfold-cmake.patch` faz isso). Compilador Fortran 2008, CMake e `fypp`, como para compilar o DFTB+. O executável fica em `<build>/app/unfold/unfold`.

**Python.** Python 3.9 ou superior e `numpy>=1.23`. O SciPy é opcional, e usa LAPACK na diagonalização generalizada. Sem SciPy o código usa Cholesky com NumPy. Para PDF, `pillow` e `reportlab`.

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
| `CheckBasis` | Só roda o teste da regra de soma e imprime min/média/máx. |

### 2.3 Chaves de `unfold_in.hsd`

| Chave | Padrão | Significado |
|---|---|---|
| `Action` | `Unfold` | Veja 2.2. |
| `DetailedXML` | `detailed.xml` | Geometria, k-pontos, número de orbitais. Obrigatório que o DFTB+ tenha escrito (`WriteDetailedXML = Yes`). |
| `EigenvecBin` | `eigenvec.bin` | Autovetores (acesso *stream*, padrão do DFTB+). Precisa de `WriteEigenvectors = Yes` em `Analysis`. |
| `BandOut` | `band.out` | Autovalores, em eV. |
| `Output` | `unfold` | Prefixo dos arquivos de saída. |
| `Basis` | — | STOs no formato Waveplot, por espécie. Pode incluir: `Basis { <<+ "wfc_C.hsd" }`. |
| `NormalizeOrbitals` | `Yes` | Normaliza as funções radiais. |
| `QPoints` | — | Caminho em $\mathbf q$ no estilo `Klines`: `1 q1 q2 q3` para o ponto inicial, `n q1 q2 q3` para $n$ pontos equidistantes terminando no ponto dado. Coordenadas fracionárias na base recíproca da célula de refolding (ou da supercélula, se não houver). Vários caminhos: blocos `Path { ... }`, cada um gera arquivo próprio. Fora da primeira zona dá o espectro na zona estendida. |
| `RefoldingLatticeVectors` | — | Vetores da célula de refolding (linhas), com unidade, por exemplo `[Angstrom]`. Sem esta chave só se calcula o espectro desdobrado completo. |
| `RefoldingGCutoff` | 10 Hartree | Corte dos vetores $\mathbf g$: $\lvert\mathbf g\rvert^2/2<$ corte. |
| `Energy { Window }` | −20 a 10 eV | Janela, relativa ao nível de Fermi. |
| `Energy { NrOfIntervals }` | 300 | Intervalos da malha de energia (a saída tem N+1 pontos). |
| `Energy { FermiLevel }` | 0 | Zero de energia, em eV. Dê o nível de Fermi do cálculo SCC. |
| `Energy { Broadening }` | 0 | Largura gaussiana. 0 deposita linearmente em dois canais, como o SIESTA. |
| `KPointSign` | −1 | Sinal entre o $\mathbf k$ do arquivo e o $\mathbf K$ físico (seção 4.2 de `METODO.md`). Não altere sem motivo. |
| `SumRuleCheck` | `Yes` | Imprime a regra de soma para o(s) primeiro(s) k-ponto(s). |
| `SumRuleQMax` | 8 bohr⁻¹ | Faixa de $q$ na soma da regra. |

### 2.4 Formato de saída

Primeira linha: `nq  ne+1  Emin  Emax  EFermi` (energias em eV). Depois, para cada $\mathbf q$: uma linha com `q1 q2 q3` cartesiano em bohr⁻¹ e o índice do segmento do caminho, seguida de `ne+1` valores de peso, um por linha. As energias são relativas ao nível de Fermi. Os pesos são normalizados de modo que cada estado contribua com 1 no total sobre todos os $\mathbf G$. Cada intervalo de energia contém um número (peso) de estados, não uma densidade.

Escritores de figuras: `claude-fortran-unfold/tools/plot_unfolded.py` (seção 4).

### 2.5 Mensagens e verificações

O programa imprime em `unfold.out` (redirecione a saída): número de pontos $\mathbf q$, número de vetores $\mathbf g$, número de k-pontos no arquivo, volume, normas das STOs (devem ser ≈ 1) e a regra de soma. A mensagem `Missing k-points` (com "N of M needed (q+g) combinations have no matching k-point") indica que o `eigenvec.bin` não foi calculado nos pontos pedidos. Use `ListKPoints` e rode o DFTB+ com a lista. A mensagem `The eigenvector file has not the expected size` indica arquivo sem `WriteEigenvectors = Yes`, ou escrito em formato que não é *stream*.

## 3. Python: `dftb_unfold.py` e `spectrum.py`

### 3.1 Linha de comando

```
python dftb_unfold.py CONFIG.json [-o saida.csv] [--plan]
python spectrum.py saida.csv --range EMIN EMAX [--sigma 0.06] [--step 0.01]
                   [--spin 1] [--weight full|refold] [-o mapa.npz] [--pdf mapa.pdf]
```

`--plan` conta pontos $\mathbf q$, vetores $\mathbf g$, K distintos e se a referência é comensurável, sem diagonalizar. Rode-o antes de uma execução grande.

### 3.2 Chaves de `unfold.json`

Caminhos são relativos ao próprio JSON.

| Chave | Obrigatória | Significado |
|---|---|---|
| `geometry` | sim | `detailed.xml` ou GEN periódico (`S` cartesiano em Å, ou `F` fracionário). Mesma ordem de átomos da matriz. |
| `matrix_dir` | não (`.`) | Pasta com `hamreal1.dat`, `overreal.dat` e, para spin colinear, `hamreal2.dat`. |
| `wfc_files` | sim | Lista de arquivos Waveplot (um por espécie ou combinados). |
| `supercell_matrix` | uma das duas | Matriz inteira $M$ com $A_{\rm SC}=M A_{\rm ref}$, vetores nas linhas. |
| `reference_cell_angstrom` | uma das duas | Três vetores da célula de referência, em Å, explícitos. |
| `path_frac` | sim | Vértices do caminho, em coordenadas fracionárias recíprocas da **referência**. |
| `intervals` | sim | Intervalos de cada segmento. Um por segmento, vértices compartilhados aparecem uma vez. |
| `labels` | não | Rótulos dos vértices, usados nos gráficos. |
| `gcut_ry` | não (0) | Corte $\lvert\mathbf g\rvert^2<$ `gcut_ry` em unidades atômicas (Ry). Corte 0 dá só `weight_full`. Equivale a 2 × o corte em Hartree do Fortran. |
| `reference_ev` | não (0) | Energia subtraída (nível de Fermi do SCC). |
| `window_ev` | não | `[min, max]` relativo a `reference_ev`. Limita as bandas gravadas. |
| `q_batch` | não (64) | Quantos vetores são projetados de uma vez. Reduz memória. |

### 3.3 Colunas do CSV

| Coluna | Significado |
|---|---|
| `point`, `distance_inv_ang`, `q1,q2,q3` | Índice (a partir de 1), distância acumulada em Å⁻¹ e coordenadas fracionárias da referência. |
| `K1_sc,K2_sc,K3_sc` | $\mathbf K$ físico reduzido à SBZ, fracionário. |
| `spin`, `band` | Índices a partir de 1. |
| `energy_ev`, `energy_rel_ev` | Autovalor e autovalor menos `reference_ev`, em eV. |
| `weight_full` | Peso em $\mathbf q$ (unfolding completo, $\mathbf g=0$). |
| `weight_refold` | Soma dos pesos em $\mathbf q+\mathbf g$ para $\lvert\mathbf g\rvert^2<$ `gcut_ry`. |

### 3.4 `spectrum.py`

Alarga cada linha com gaussiana de área unitária e largura `--sigma` (eV) e grava um NPZ com densidade por ponto, eixo de energia e distância. `--pdf` desenha uma figura (requer pillow e reportlab). Os pesos de bandas degeneradas se somam: compare sempre o subespaço degenerado.

### 3.5 Limites do código Python

DFTB periódico com Hamiltoniano real, base com uma camada radial por $l$, spin escalar ou colinear. Sem SOC, spin não colinear, contatos ou fônons. Custo $\mathcal O(N_{\rm orb}^3)$ por (K, spin), memória $\mathcal O(N_{\rm orb}^2)$.

## 4. Visualização

**Saída Fortran.** `claude-fortran-unfold/tools/plot_unfolded.py ARQUIVO [-o fig.png] [--vmax FRAC] [--ref band.out FERMI]` desenha o mapa em tons de cinza e, com `--ref`, sobrepõe em vermelho as bandas de um cálculo primitivo de referência. O formato é o mesmo do SIESTA, e os scripts de `Util/Unfolding/Src` devem funcionar, mas isso não foi testado aqui.

**Saída Python.** `spectrum.py` e seu `--pdf`.

## 5. Unidades e convenções

| Grandeza | Fortran | Python |
|---|---|---|
| Energias | eV | eV |
| $\mathbf q$ na entrada | fracionário (base recíproca da célula de refolding) | fracionário (idem) |
| $\mathbf q$ na saída | cartesiano, bohr⁻¹ | fracionário no CSV, distância em Å⁻¹ |
| Corte em $\mathbf g$ | Hartree: $\lvert\mathbf g\rvert^2/2<c$ | Ry: $\lvert\mathbf g\rvert^2<c_{\rm Ry}$, com $c_{\rm Ry}=2c_{\rm Ha}$ |
| Peso | adimensional, soma 1 por estado sobre $\mathbf G$ | adimensional |

## 6. Solução de problemas

| Sintoma | Causa provável | O que fazer |
|---|---|---|
| Regra de soma longe de 1 (por exemplo 0,4–3,8) | `wfc` de outro conjunto que os SK, ordem de orbitais errada, ou `KPointSign` errado | HOWTO 6 |
| Regra de soma 0,9 a 0,97 | `SumRuleQMax` pequeno, ou ajuste das STOs imperfeito | Aumente `SumRuleQMax`. Pequenos desvios são do ajuste STO. |
| `Missing k-points` (Fortran) | `eigenvec.bin` não calculado nos k-pontos de `ListKPoints` | Refaça o passo 3 do HOWTO 1 |
| `geometry and matrix have different numbers of atoms` (Python) | Geometria e matrizes de cálculos diferentes | Use o `detailed.xml` do mesmo cálculo |
| `commensurate: False` (Python, `--plan`) | Referência não é submúltiplo da SC | Normal para amorfos e defeitos deformados, mas pesa em custo |
| Pesos fora das bandas do cristal perfeito | Defeito ou desordem, ou erro de referência | HOWTO 7 |
| Refolded sem peso em grafeno 2D | $p_z$ ímpar em $z$ e $q_z=0$ | HOWTO 8 |
