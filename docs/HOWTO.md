# HOWTO: receitas por tarefa

Receitas práticas para os dois programas. As opções estão detalhadas em `MANUAL.md` e a teoria em `METODO.md`. Os exemplos usam o diamante cúbico de 8 átomos desdobrado na célula primitiva fcc (parâmetros mio-1-1), que é o caso validado.

Conteúdo:

1. Tutorial completo no diamante, com o código Fortran
2. Tutorial completo no diamante, com o código Python
3. Definir a supercélula e a célula de referência
4. Escolher caminho em k, janela de energia e nível de Fermi
5. Escolher o corte dos vetores g
6. Diagnosticar uma regra de soma ruim
7. Desdobrar um sistema com defeito e compará-lo ao cristal perfeito
8. Grafeno e sistemas 2D
9. Sistemas grandes: custo e escolha do código
10. Comparar com as bandas da célula primitiva
11. Spin

---

## 1. Tutorial completo no diamante (Fortran)

**Pré-requisito:** DFTB+ compilado com `app/unfold` (veja o `README.md` da raiz), os arquivos SK **mio-1-1** e os executáveis em `<build>/app`.

O script do exemplo faz tudo:

```bash
cd claude-fortran-unfold/app_unfold/examples/diamond
SKDIR=/caminho/para/mio-1-1 BIN=/caminho/para/build/app ./run_all.sh
```

O que ele faz, passo a passo (útil para adaptar a outros sistemas):

1. **Cálculo SCC da supercélula** (`run1/`): converge as cargas e grava `charges.bin`.
2. **`Action = ListKPoints`** (`unfold_list.hsd`): `unfold` lê o caminho `QPoints` e escreve `unfold.kpoints.hsd` com os k-pontos necessários.
3. **Cálculo não-SCC nesses k-pontos** (`run2/`): o `dftb_in.hsd` ganha

   ```hsd
   Hamiltonian = DFTB {
     SCC = Yes
     MaxSCCIterations = 1
     ReadInitialCharges = Yes
     <<+ "unfold.kpoints.hsd"
   }
   Analysis { WriteEigenvectors = Yes }
   Options { WriteDetailedXML = Yes }
   ```

4. **`Action = Unfold`** (`unfold_run.hsd`): o `@FERMI@` é substituído pelo nível de Fermi de `run1/detailed.out`. Saída: `run2/unfold.unfoldedBands` e `run2/unfold.refoldedBands`.

Para ver o resultado:

```bash
python3 claude-fortran-unfold/tools/plot_unfolded.py run2/unfold.refoldedBands -o refolded.png
```

Esperado no diamante: o mapa segue as bandas da célula primitiva (veja `docs/figures/diamond_unfold.png`) e a regra de soma no `unfold.out` fica em 0,976–0,999.

## 2. Tutorial completo no diamante (Python)

1. **Cálculo SCC da supercélula**, para obter `charges.bin`, igual ao passo 1 da receita anterior. Ative em `Options`: `WriteDetailedXML = Yes`.
2. **Exporte as matrizes.** Em uma pasta nova, copie `charges.bin` e use o arquivo `validation/diamond/dftb_in_writeRealHS.hsd` como modelo. Ele traz `ReadInitialCharges = Yes`, `RecomputeAfterDensity = No` e `WriteRealHS = Yes` (veja o próprio arquivo, incluindo o bloco `SlaterKosterFiles`, cujo `Prefix` você precisa ajustar ao seu caminho). O DFTB+ grava `hamreal1.dat` e `overreal.dat` e termina com "Hamilton/Overlap written, exiting program.".
3. **Copie** para a mesma pasta o `detailed.xml` do passo 1 e o `wfc_C.hsd` do conjunto mio-1-1.
4. **Configure** `validation/diamond/unfold.json`:

   ```json
   {"geometry": "detailed.xml", "wfc_files": ["wfc_C.hsd"],
    "supercell_matrix": [[-1,1,1],[1,-1,1],[1,1,-1]],
    "path_frac": [[0.5,0.5,0.5],[0,0,0],[0,0.5,0.5]],
    "intervals": [20, 25], "gcut_ry": 20.0, "reference_ev": -1.0731}
   ```

   `reference_ev` é o nível de Fermi do passo 1.

5. **Planeje e rode:**

   ```bash
   python dftb_unfold.py unfold.json --plan
   python dftb_unfold.py unfold.json -o saida.csv
   python spectrum.py saida.csv --range -25 15 --sigma 0.1 --weight refold -o mapa.npz --pdf mapa.pdf
   ```

   O `--plan` deve informar 8 átomos, 46 pontos, 113 vetores de refolding e 45 K distintos, e `commensurate: True`. O CSV tem 1472 linhas.

## 3. Definir a supercélula e a célula de referência

- **Fortran:** informe a célula de refolding em `RefoldingLatticeVectors` (linhas, com unidade). A supercélula vem do `detailed.xml`.
- **Python:** `supercell_matrix` é a matriz inteira $M$ tal que $A_{\rm SC}=M\,A_{\rm ref}$, com vetores nas **linhas**. Se a supercélula não for múltipla inteira da referência, use `reference_cell_angstrom` com os três vetores.
- **Conferência do $M$:** para o diamante, a célula cúbica de 8 átomos contém 4 células primitivas fcc, e a matriz acima tem determinante 4. O determinante de $M$ é o número de células de referência na supercélula.
- Células incomensuráveis (amorfos, defeitos com relaxação forte) funcionam no Python (veja `--plan`). Cada $\mathbf g$ pode então exigir uma diagonalização própria, e o custo cresce com o corte.

## 4. Escolher caminho em k, janela de energia e nível de Fermi

- **Caminho.** Os dois programas usam coordenadas fracionárias da base recíproca da **referência**. Para L–Γ–X do fcc: `0.5 0.5 0.5` → `0 0 0` → `0 0.5 0.5`.
  - Fortran: `QPoints = { 1 0.5 0.5 0.5  20 0.0 0.0 0.0  25 0.0 0.5 0.5 }`.
  - Python: `"path_frac": [[0.5,0.5,0.5],[0,0,0],[0,0.5,0.5]]`, `"intervals": [20, 25]`.
- **Fortran:** depois de mudar o caminho, **rode de novo** `ListKPoints` e o cálculo no DFTB+ (receita 1, passos 2 e 3), pois os k-pontos calculados mudam. O Python não precisa disso.
- **Nível de Fermi.** Pegue de `detailed.out` do cálculo SCC (linha "Fermi level"). Fortran: `Energy { FermiLevel [eV] = ... }`. Python: `reference_ev`.
- **Janela.** Fortran: `Energy { Window [eV] = {-25 15} }`. Python: `window_ev` no JSON, ou `--range` no `spectrum.py`. Em semicondutores a banda de valência inicia perto de −20 eV em relação ao Fermi (diamante), então uma janela de −25 a 15 eV cobre todas as bandas ocupadas.
- **Alargamento.** Fortran: `Energy { Broadening }` (0 = sem alargamento). Python: `--sigma` em `spectrum.py` (padrão 0,06 eV).

## 5. Escolher o corte dos vetores g

O corte controla a precisão do refolding. O peso em $\mathbf q+\mathbf g$ decai com $\lvert\mathbf g\rvert$, mas devagar, por causa do truncamento brusco da STO em $r_c$.

1. Rode com um corte moderado (Fortran: 10 Hartree, o padrão; Python: 20 Ry).
2. Rode de novo com o dobro.
3. Compare o peso total refolded na janela de energia. Se mudar pouco, o corte basta.

**Conversão:** `gcut_ry = 2 × corte_Hartree`. As duas escolhas acima são equivalentes e dão 113 vetores no diamante.

Cuidado: um corte maior aumenta o custo linearmente no número de vetores. No Fortran, as diagonalizações já foram feitas pelo DFTB+, então o corte só afeta o pós-processamento.

## 6. Diagnosticar uma regra de soma ruim

A regra de soma é o teste mais importante (veja a seção 5 de `METODO.md`). Ela é impressa em `unfold.out` pelo Fortran. O valor esperado é perto de 1 (0,97–1,00).

**Passo 1: confira o conjunto.** O `wfc` e os arquivos `.skf` devem ser do **mesmo conjunto** (por exemplo, os dois mio-1-1). Com SK 3ob e `wfc` do mio, o resultado foi 0,4–3,8.

**Passo 2: rode só o teste.** No Fortran, `Action = CheckBasis`. Ele examina k-pontos genéricos do arquivo, e é rápido.

**Passo 3: confira o sinal de K.** O TRIM (por exemplo Γ) não distingue `KPointSign = ±1`. Se os valores estão ruins só em k-pontos genéricos, teste `KPointSign = 1`. O padrão −1 foi verificado no diamante e no grafeno.

**Passo 4: confira a faixa de $q$.** Aumente `SumRuleQMax` (padrão 8 bohr⁻¹). Se a soma sobe e satura, era a faixa.

**Passo 5: confira a ordem dos orbitais e $l$.** As harmônicas seguem a ordem do Waveplot (p: y, z, x). Bases com $l>3$ não são suportadas.

No Python não há impressão da regra de soma. Para checá-la, rode primeiro o Fortran no mesmo sistema, ou use o `test_unfold.py`, que verifica $C^\dagger S C=I$ em casos de teste.

## 7. Desdobrar um sistema com defeito e compará-lo ao cristal perfeito

Receita genérica (a do artigo para a vacância em Si [1]):

1. Calcule o cristal perfeito na mesma supercélula (ou na primitiva) e o sistema com defeito, com a mesma parametrização e a mesma célula de refolding.
2. Desdobre os dois com o mesmo caminho, a mesma janela e o mesmo `Fermi`/`reference_ev` (cada um com o **seu** nível de Fermi).
3. Compare os mapas, ou a diferença. Estados novos aparecem como peso dentro do gap. O topo da banda de valência perde peso e pode descer em energia.
4. Atenção à comparação direta: o nível de Fermi do defeito difere do perfeito. Se desenhar os dois com origem em $E_F$, a diferença de alinhamento aparece como deslocamento. Comente isso nos resultados.

Mantenha os parâmetros numéricos (malha de k, tolerância SCC) consistentes entre os dois cálculos, para que as diferenças venham do defeito e não do método numérico.

## 8. Grafeno e sistemas 2D

- Mantenha o **mesmo vácuo** na supercélula e na célula de referência.
- A referência 3D para uma folha tem o terceiro vetor igual ao da supercélula.
- Orbitais $p_z$ são ímpares em $z$, então sua transformada se anula em $q_z=0$. O peso total no plano $q_z=0$ pode desaparecer nas bandas π só no `weight_full`. O refolding com $g_z\neq0$ recupera o peso. Para ver as bandas π, use o `weight_refold` (Python) ou `RefoldingLatticeVectors` com o vetor vertical (Fortran).
- O exemplo validado é a supercélula 3×3 (18 átomos), em `validation/graphene/` (`gen.py`, `mk.sh`). Os caminhos nos scripts apontam para o diretório de compilação local: ajuste-os.
- Para o Python, o exemplo `examples/graphene/validate_graphene.py` exige `dftb+`, os arquivos SK e o `wfc` (veja o `README.md` da pasta). Não foi reexecutado neste projeto.

## 9. Sistemas grandes: custo e escolha do código

| | Fortran | Python |
|---|---|---|
| Quantas diagonalizações? | Nenhuma extra. Usa as do DFTB+ nos k-pontos de `ListKPoints` | Uma por (K, spin), densa em Python |
| Para supercélula comensurável | Poucos k-pontos distintos, porque muitos $\mathbf q+\mathbf g$ dobram no mesmo $\mathbf K$ | Idem. Veja `--plan` |
| Gargalo | O cálculo do DFTB+ | Memória $\mathcal O(N_{\rm orb}^2)$ e tempo $\mathcal O(N_{\rm orb}^3)$ por K |

Recomendações:

- Para muitos átomos e caminho fixo, prefira o **Fortran**.
- Para explorar caminhos diferentes sem repetir cálculos no DFTB+, ou para K arbitrário, prefira o **Python**.
- No Python, rode sempre `--plan` antes. Se faltar memória, reduza `q_batch` e use `window_ev` para limitar as bandas.

## 10. Comparar com as bandas da célula primitiva

1. Rode um cálculo DFTB+ de bandas da célula primitiva, com o mesmo caminho. Isso gera um `band.out`.
2. Fortran: use `plot_unfolded.py` com `--ref`:

   ```bash
   python3 claude-fortran-unfold/tools/plot_unfolded.py run2/unfold.refoldedBands \
       --ref prim/bands/band.out -1.0731 -o comparacao.png
   ```

   O último número é o nível de Fermi usado como origem das bandas de referência.

3. Para uma comparação numérica (fração do peso a menos de 0,15 eV das bandas primitivas e desvio do centróide), veja `validation/diamond/compare.py`. Uso: `python3 compare.py <nivel_de_fermi>`. Ele lê `run2/unfold.refoldedBands` e `prim/bands/band.out`, então rode-o na pasta do cálculo.

## 11. Spin

- **Sem polarização de spin:** caso validado.
- **Colinear polarizado:** o Fortran grava um arquivo por canal (segundo o README do `app_unfold`), o Python usa `hamreal1.dat` e `hamreal2.dat` e produz as colunas `spin` (1 ou 2). Esse caso **não foi validado** neste projeto.
- **Não colinear e SOC:** sem suporte nos dois.
