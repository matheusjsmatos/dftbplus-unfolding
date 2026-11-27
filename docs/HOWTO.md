# HOWTO de unfolding de bandas no DFTB+

Os comandos partem da **raiz deste repositório**. O exemplo de grafeno arquivado roda em Python no Windows ou no Linux. As receitas que usam `bash`, `sed` ou variáveis como `$BIN` usam Bash em Linux. Consulte o [Manual](MANUAL.md) para a referência das entradas e o [Método](METODO.md) para as equações.

## Rodar e plotar o grafeno 2×2 já calculado

O pacote inclui matrizes H/S e bandas de referência calculadas no DFTB+ para o grafeno puro. Para reproduzir o unfolding com a base mio-1-1 e gerar a figura em PDF, execute os comandos abaixo. A instalação de Pillow e ReportLab serve para a figura.

```bash
python -m pip install -r chatgpt-python-unfold/requirements.txt
python -m pip install Pillow reportlab
python chatgpt-python-unfold/tools/check_graphene.py --output results/graphene
python chatgpt-python-unfold/spectrum.py results/graphene/unfolded.csv --weight refold --range -12 6 --step 0.02 --sigma 0.08 -o results/graphene/mapa.npz --pdf results/graphene/bandas.pdf
```

Abra **`results/graphene/bandas.pdf`**. A intensidade representa o peso espectral refoldado ao longo de Γ–K–M–Γ. A abscissa é a distância acumulada em Å⁻¹. A janela vertical vai de −12 a +6 eV em relação ao zero do arquivo DFTB+ (`reference_ev = 0` no exemplo). Para apresentar $E-E_F$ em um cálculo próprio, use em `reference_ev` o Fermi do seu SCC e gere novamente o CSV. Para representar apenas a projeção em $\mathbf q$, troque `--weight refold` por `--weight full`. A seção 8 explica a diferença para as bandas π.

A execução gera **2.272 estados no CSV** (71 pontos × 32 bandas). Em `results/graphene/validation.json`, procure `"passed": true`. Esse teste compara os autovalores com os da supercélula e a posição dos estados com peso com as bandas da célula primitiva. Os arquivos de entrada e a configuração estão em [referência do grafeno](../chatgpt-python-unfold/examples/graphene/reference/).

## Executar os testes do projeto

Com as dependências Python instaladas, rode a suíte offline a partir da raiz.

```bash
python Testes/run.py --output results/testes
```

Para gerar novos cálculos de grafeno **e** diamante e comparar os pesos absolutos de Python e Fortran, use os executáveis `dftb+` e `unfold` compilados com **DFTB+ 26.1**. O workflow contém os comandos para obter o pacote-fonte 26.1, incluir o diretório `app/unfold` e compilar ambos os alvos. A opção `--mio` recebe a raiz do repositório de parâmetros, contendo `skfiles/C-C.skf` e `extras/wfc.hsd`. A suíte verifica o SHA-256 desses dois arquivos. [Compilação testada](../.github/workflows/validar-unfolding.yml) · [Versão de mio usada](https://github.com/dftbparams/mio/tree/263e52419ce15222c864df9b3f6fbf58bfbdfa0f)

```bash
python Testes/run.py native --dftb /caminho/build/native/app/dftb+/dftb+ --unfold /caminho/build/native/app/unfold/unfold --mio /caminho/mio --output results/native
```

A execução cria uma pasta nova por material. Em `results/native`, procure `validation.json` para a comparação com as bandas primitivas, `states.json` para a verificação dos autovetores na métrica S e `comparison.json` para os mapas. O teste compara separadamente depósito linear e alargamento gaussiano com σ = 0,08 eV. Nesse último caso, a diferença L1 relativa no mapa refoldado foi **0,0236% para o grafeno** e **0,0255% para o diamante**. [Critérios, resultados e referências](../chatgpt-python-unfold/docs/DEPURACAO.md)

## 1. Tutorial completo no diamante (Fortran)

O exemplo usa uma célula cúbica de oito átomos de carbono e a desdobra sobre a primitiva fcc. Ajuste os dois caminhos absolutos abaixo. `SKDIR` aponta para a pasta que contém `C-C.skf`, e `BIN` para a pasta comum dos dois executáveis. A cópia do exemplo em uma pasta nova preserva os resultados anteriores, pois `run_all.sh` recria `run1/` e `run2/`.

```bash
export SKDIR=/caminho/mio/skfiles
export BIN=/caminho/build/native/app
mkdir -p results
DIAM=$(mktemp -d "$PWD/results/diamante.XXXXXX")
cp -R claude-fortran-unfold/app_unfold/examples/diamond/. "$DIAM/"
bash "$DIAM/run_all.sh"
python -m pip install matplotlib
python claude-fortran-unfold/tools/plot_unfolded.py "$DIAM/run2/unfold.refoldedBands" -o "$DIAM/refolded.png"
```

1. O cálculo SCC em `run1/` converge a densidade e grava `charges.bin`.
2. Com `Action = ListKPoints`, o `unfold` lê o caminho da referência em `unfold_list.hsd` e escreve `run1/unfold.kpoints.hsd`.
3. O exemplo copia as cargas e a lista de k-pontos para `run2/`. O DFTB+ lê as cargas com `ReadInitialCharges = Yes`, executa uma iteração com `MaxSCCIterations = 1` e grava `band.out`, `detailed.xml` e `eigenvec.bin`. Mantenha juntos os três arquivos dessa mesma execução.
4. Com `Action = Unfold`, o script substitui `@FERMI@` pelo Fermi em eV do SCC e produz `run2/unfold.unfoldedBands` e `run2/unfold.refoldedBands`. Abra `refolded.png` para ver o mapa. A figura usa o índice do ponto de caminho no eixo horizontal. A regra de soma fica em `run2/unfold.out`. [Entradas do exemplo](../claude-fortran-unfold/app_unfold/examples/diamond/)

## 2. Tutorial completo no diamante (Python)

Os comandos a seguir reutilizam o SCC da receita 1, na mesma sessão Bash, com `$DIAM`, `$BIN` e `$SKDIR` definidos. Se você já tem um cálculo próprio, substitua o modelo pela **mesma geometria, parametrização, ordem de átomos e orbitais** usadas para convergir as cargas.

1. Exporte os blocos reais de H e S em uma pasta separada. O modelo recebe as cargas convergidas e ativa `ReadInitialCharges`, `WriteRealHS` e `ShowFoldedCoords`, com `RecomputeAfterDensity = No`.

   ```bash
   EXP="$DIAM/export"
   mkdir -p "$EXP"
   cp "$DIAM/run1/charges.bin" "$EXP/"
   sed "s|/caminho/para/mio-1-1/|${SKDIR%/}/|" validation/diamond/dftb_in_writeRealHS.hsd > "$EXP/dftb_in.hsd"
   (cd "$EXP" && "$BIN/dftb+/dftb+" > matrizes.log 2>&1)
   grep 'Hamilton/Overlap written' "$EXP/matrizes.log"
   ls -l "$EXP/hamreal1.dat" "$EXP/overreal.dat"
   ```

   O DFTB+ pode encerrar a exportação com código de saída diferente de zero. Confira a mensagem no log **e** os dois arquivos preenchidos antes de prosseguir.

2. Copie a geometria XML do SCC, os orbitais de carbono e a configuração do exemplo. Os caminhos dentro do JSON são relativos à pasta do próprio JSON.

   ```bash
   cp "$DIAM/run1/detailed.xml" "$DIAM/wfc_C.hsd" "$EXP/"
   cp validation/diamond/unfold.json "$EXP/"
   ```

   Confira `reference_ev` em `$EXP/unfold.json`. O valor fornecido, −1,0731 eV, é a referência usada no teste do diamante. Para plotar $E-E_F$ do seu cálculo, use o Fermi do SCC em eV.

3. Planeje e calcule a partir da raiz do repositório.

   ```bash
   python chatgpt-python-unfold/dftb_unfold.py "$EXP/unfold.json" --plan
   python chatgpt-python-unfold/dftb_unfold.py "$EXP/unfold.json" -o "$EXP/saida.csv"
   python chatgpt-python-unfold/spectrum.py "$EXP/saida.csv" --range -25 15 --sigma 0.1 --weight refold -o "$EXP/mapa.npz" --pdf "$EXP/mapa.pdf"
   ```

   O plano informa 8 átomos, 46 pontos, 113 vetores de refolding, 45 K distintos e `commensurate: True`. O CSV contém 1.472 estados além do cabeçalho. A figura usa σ = 0,1 eV. Para reproduzir a comparação gaussiana da suíte nativa no diamante, use `--sigma 0.08 --step 0.1`.

## 3. Definir a supercélula e a célula de referência

- No Fortran, informe a célula de refolding em `RefoldingLatticeVectors [Angstrom]`, um vetor por linha. A supercélula vem de `detailed.xml`.
- No Python, `supercell_matrix` é a matriz inteira $M$ tal que $A_{\rm SC}=M A_{\rm ref}$, com vetores nas linhas. Para o diamante, $M$ tem linhas $(-1,1,1)$, $(1,-1,1)$ e $(1,1,-1)$. Seu determinante em módulo é 4, o número de primitivas na célula cúbica. [JSON do diamante](../validation/diamond/unfold.json)
- Para uma referência incomensurável, use `reference_cell_angstrom` com os três vetores. A comensurabilidade depende dos **vetores de rede**. Uma vacância ou o relaxamento interno dos átomos pode conservar esses vetores. Use `--plan` para saber se diferentes $\mathbf q+\mathbf g$ exigirão diferentes $\mathbf K$.

## 4. Escolher caminho, janela de energia e nível de Fermi

- Os dois programas recebem o caminho em coordenadas fracionárias da base recíproca da **referência**. Para L–Γ–X do diamante, use `[0.5,0.5,0.5]`, `[0,0,0]` e `[0,0.5,0.5]`. No Python, `"intervals": [20,25]` produz $1+20+25=46$ pontos. No Fortran, veja a entrada multilinha `QPoints`. [Exemplo HSD](../claude-fortran-unfold/app_unfold/examples/diamond/unfold_run.hsd)
- Ao mudar caminho, célula de referência ou corte no Fortran, rode `ListKPoints` para saber quais estados são necessários. Refaça o cálculo de bandas quando faltarem k-pontos no arquivo de autovetores. O Python pode reutilizar as matrizes H/S para outro caminho.
- Para usar $E-E_F$, leia o valor em **eV** da linha `Fermi level` do `detailed.out` do SCC. Ajuste `Energy { FermiLevel [eV] = ... }` no Fortran e `reference_ev` no JSON. A mesma origem deve valer para o mapa e para as bandas sobrepostas.
- No Fortran, `Energy { Window [eV] = {-25 15} }` define a malha. No Python, `window_ev` filtra os estados **antes** de gravar o CSV. Para preservar as caudas das linhas próximas às bordas, guarde o CSV completo e limite a figura com `--range`. A janela do teste do diamante é de −25 a +15 eV.
- `Broadening [eV]` e `--sigma` são o desvio padrão da gaussiana. Com zero, ambos distribuem linearmente o peso pelos dois pontos de energia vizinhos. Para comparar intensidades absolutas, use o mesmo passo e a mesma largura. O Fortran grava peso por intervalo, enquanto `spectrum.py` fornece densidade em eV⁻¹. Divida os valores do mapa Fortran por ΔE ou use o comparador da suíte.

## 5. Escolher o corte dos vetores g

Comece com `RefoldingGCutoff [Hartree] = 10` no Fortran ou `"gcut_ry": 20` no Python. Essa escolha produz 113 vetores no diamante. Dobre o corte e compare o peso integrado e as bandas de interesse mantendo janela, passo e caminho fixos. A conversão é $E_{\rm Ry}=2 E_{\rm Ha}$. A quantidade de vetores depende da célula de referência.

Cortes maiores aumentam a soma sobre $\mathbf g$. Em uma supercélula comensurável, $\mathbf q+\mathbf g$ pertence ao mesmo $\mathbf K$ dobrado para cada $\mathbf q$. Para uma referência incomensurável, reveja a lista de k-pontos conforme a seção 4.

## 6. Diagnosticar a regra de soma

O Fortran estima a soma de pesos em ondas planas para autovetores normalizados. No exemplo original do diamante, o intervalo relatado foi 0,976–0,999. O resultado depende da cobertura em momento e da consistência dos radiais Waveplot com a sobreposição DFTB+. [Método e regra de soma](METODO.md)

1. Confira se os radiais `wfc` e os arquivos `.skf` pertencem ao mesmo conjunto de parâmetros. O exemplo usa mio-1-1 nos dois.
2. Rode `Action = CheckBasis` com os arquivos nativos da execução de bandas. Inclua k-pontos genéricos. Só Γ é insuficiente para testar a convenção de fase.
3. No DFTB+ validado, use `KPointSign = -1`. Confirme que os autovetores, o XML e as energias vieram da mesma execução antes de investigar outra convenção.
4. Aumente `SumRuleQMax` (padrão de 8 bohr⁻¹) e examine a convergência. Confira também a ordem dos orbitais do Waveplot, por exemplo $p_y,p_z,p_x$. As implementações cobrem $s,p,d,f$.
5. Para testar a etapa anterior à projeção, consulte o `states.json` do modo nativo. A métrica $C^\dagger S C=I$ verifica a normalização dos autovetores na base LCAO. A soma em ondas planas testa também os radiais e o corte.

## 7. Comparar um defeito com o cristal perfeito

1. Calcule o cristal perfeito e o sistema com defeito com parametrização, tolerância SCC e densidades de malha k comparáveis.
2. Desdobre os dois sobre a mesma célula de referência e use o mesmo caminho, janela, passo, corte e alargamento.
3. Para estudar os estados em torno do potencial químico, represente cada mapa em relação ao **seu próprio** $E_F$. Para medir um deslocamento entre os dois sistemas, alinhe antes as energias por uma referência física comum e registre essa escolha.
4. Compare a posição e a intensidade das bandas. Procure perda de peso de uma banda ideal ou estados dentro do gap. Para comparar amplitudes ou subtrair os mapas, mantenha a mesma escala de intensidade. [Mayo, Yndurain e Soler, exemplo da vacância em Si](https://doi.org/10.1088/1361-648X/ab6e8e)

## 8. Grafeno e sistemas 2D

- Mantenha o mesmo vetor de vácuo na supercélula e na referência para comparar pesos com a mesma normalização de volume.
- No grafeno plano, a simetria dos orbitais $p_z$ pode suprimir a projeção direta em $q_z=0$. O mapa refoldado soma intensidades em $\mathbf q+\mathbf g$, incluindo termos com $g_z\ne0$ quando o corte os alcança. Assim as bandas π podem reaparecer. No Python use `--weight refold`. No Fortran forneça os três `RefoldingLatticeVectors` e examine `.refoldedBands`.
- O exemplo arquivado e a comparação nativa automatizada usam grafeno ideal **2×2, 8 C**, 25 Ry e o caminho Γ–K–M–Γ. [Geometria e parâmetros](../chatgpt-python-unfold/examples/graphene/)
- Há uma receita adicional **3×3, 18 C** em [validation/graphene](../validation/graphene/). `gen.py` grava a primitiva e a supercélula no diretório corrente. Defina `BIN` e `SKDIR` como na receita 1 e copie o bloco radial do carbono mio-1-1 para a pasta da execução.

  ```bash
  mkdir -p results
  G3=$(mktemp -d "$PWD/results/grafeno3.XXXXXX")
  cp -R validation/graphene/. "$G3/"
  cp claude-fortran-unfold/app_unfold/examples/diamond/wfc_C.hsd "$G3/"
  (cd "$G3" && python gen.py && bash mk.sh)
  ```

  Confira a existência de `$G3/run2/unfold.refoldedBands` e o log `$G3/run2/unfold.out`, pois o script prossegue para o cálculo primitivo mesmo se o unfolding falhar. O caminho fornecido começa em (⅓,⅓,0). Para usar o mesmo vértice K = (⅔,⅓,0) do teste 2×2, ajuste os dois arquivos `unfold_*.hsd` e o bloco `Klines` da primitiva em `mk.sh`, mantendo a mesma ordem de pontos.

## 9. Estimar o custo em sistemas grandes

| Operação | Fortran | Python |
|---|---|---|
| Autovetores | Lê os que o DFTB+ calculou nos k-pontos pedidos por `ListKPoints`. | Reconstrói H e S e diagonaliza uma vez por par (K, spin). |
| Custos a estimar | Cálculo DFTB+, tamanho de `eigenvec.bin` e soma sobre $\mathbf g$. | Matrizes densas de ordem $N_{\rm orb}$, diagonalização cúbica por K e projeções por lote. |

Use `--plan` para contar pontos e vetores antes de um cálculo grande. `q_batch` controla a memória temporária da projeção e `window_ev` limita os estados projetados e gravados. A diagonalização densa ainda exige memória da ordem de $N_{\rm orb}^2$ por K. O Python permite testar vários caminhos com o mesmo conjunto H/S.

## 10. Sobrepor as bandas da célula primitiva

1. Gere um `band.out` da primitiva com o **mesmo caminho físico e os mesmos pontos, na mesma ordem**. O `run_all.sh` do diamante produz o mapa da supercélula. O modo `native` em `Testes/run.py` gera também as bandas primitivas, em `pc_bands/band.out`.
2. O gráfico Fortran usa o índice do ponto ao longo do caminho. Defina `MAPA` como o arquivo `.refoldedBands` e `PRIMITIVA` como o `band.out` correspondente. O exemplo usa a origem −1,0731 eV do teste de diamante.

   ```bash
   python claude-fortran-unfold/tools/plot_unfolded.py "$MAPA" --ref "$PRIMITIVA" -1.0731 -o comparacao.png
   ```

3. Para medir a fração de peso a menos de 0,15 eV das bandas primitivas e o desvio dos centróides, use [compare.py](../validation/diamond/compare.py). O script lê os nomes fixos `run2/unfold.refoldedBands` e `prim/bands/band.out` relativamente à pasta do cálculo. A suíte em `Testes/` compara também os pesos absolutos dos dois programas.

## 11. Spin

Os exemplos nativos de grafeno e diamante usam um canal escalar. Para spin colinear, o Fortran grava os sufixos `.spinUp` e `.spinDn`. No Python, inclua `hamreal1.dat` e `hamreal2.dat` com a sobreposição comum e escolha `--spin 1` ou `--spin 2` em `spectrum.py`. Há teste unitário com uma base colinear sintética. A comparação nativa publicada usa carbono sem polarização. Os dois programas requerem Hamiltonianos escalares ou colineares e não implementam SOC.
