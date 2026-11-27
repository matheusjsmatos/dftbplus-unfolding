# Depuração e equivalência das duas implementações

O alvo são **bandas eletrônicas calculadas com DFTB+** em bases LCAO não ortogonais. A implementação Python do ChatGPT e o aplicativo Fortran identificado como Claude pelo mantenedor calculam o mesmo observável quando se igualam os autovetores, a geometria, os orbitais radiais e as condições de projeção. A comparação numérica depende também de como os estados são depositados na malha de energia. [Resultados anteriores do repositório](../../docs/COMPARACAO.md)

## O que cada implementação faz

O Python usa `WriteRealHS` para reconstruir as matrizes de Bloch dos blocos reais H e S e resolve o problema generalizado $HC=SCE$, com $C^\dagger SC=I$. As integrais radiais são calculadas por Gauss–Legendre. O CSV preserva a energia e os dois pesos de cada estado antes do alargamento. O Fortran lê os coeficientes já diagonalizados de `eigenvec.bin`, os pontos k e a geometria de `detailed.xml` e as energias de `band.out`. Tabula as radiais por Simpson em 6.001 pontos, com interpolação linear em passos de 0,01 bohr⁻¹. A versão publicada renormaliza as radiais por padrão. Para cotejar os coeficientes Waveplot tais como distribuídos, o comparador escreve `NormalizeOrbitals = No`. A opção `--normalize-native-orbitals` reproduz o ajuste radial do programa original.

Para cada momento $q$ escolhe-se o estado em $K=q\pmod{G_s}$, na zona da supercélula. Na convenção do Waveplot empregada aqui, o k nativo do arquivo de autovetores é $-K$. Com volume $V_s$, centros $\tau_a$ e harmônicos reais na ordem do DFTB+, a amplitude é

$$
A_n(q)=\frac{4\pi}{\sqrt{V_s}}\sum_{a,l,m}C_{alm,n}(K)e^{-iq\cdot\tau_a}(-i)^lY_{lm}(\hat q)
\int_0^{r_c}r^2 R_{al}(r)j_l(|q|r)\,dr.
$$

A interferência ocorre antes do módulo quadrado. O espectro totalmente desdobrado usa $|A_n(q)|^2$ na energia $E_{nK}$. O refolding soma as **intensidades** nos momentos $k+g$ sobre a rede recíproca da célula de referência. Cada termo usa o $K$ compatível. A soma de $g$ com componente perpendicular ao plano recupera peso das bandas π do grafeno que a simetria de $p_z$ pode suprimir na projeção direta com $q_z=0$. Em degenerescências, comparam-se pesos somados no subespaço ou na malha de energia. [Artigo de Mayo, Yndurain e Soler, Eqs. 6, 9 e 10–14](https://doi.org/10.1088/1361-648X/ab6e8e) · [Utilitário SIESTA](https://gitlab.com/siesta-project/siesta/-/tree/dev/Util/Unfolding) · [Waveplot](https://dftbplus-recipes.readthedocs.io/en/latest/basics/waveplot.html)

## O que os números permitem concluir

| Ensaio | Resultado |
|---|---|
| Testes Python | Os 11 testes passaram. Oito cobrem o núcleo do unfolding, incluindo cadeia duplicada, interferência complexa s/p, spin colinear e matrizes não ortogonais. Três verificam o comparador e o depósito espectral na malha de energia. |
| Grafeno puro 2×2, Python contra `band.out` DFTB+ da supercélula | 71 pontos, 32 bandas e 2.272 linhas. Maior diferença de energia de **5,183 × 10⁻⁵ eV**. |
| Grafeno 2×2, comparação com as bandas nativas da primitiva | Maior peso refoldado fora das bandas primitivas de **2,12 × 10⁻²⁰**, usando 2 × 10⁻⁴ eV de tolerância. |
| Reimplementação independente das equações radiais do anexo em NumPy, mesmos autovetores H/S | Maior erro, dividido pelo maior peso Python, de **3,66 × 10⁻⁵** no full e **2,87 × 10⁻⁵** no refold, somando degenerescências. Esse ensaio verifica o projetor e não executa o leitor binário Fortran. |
| Diamante no relatório anterior do repositório | O relatório registra diferença relativa L2 de **5 × 10⁻⁴** no full e excesso de **1,7–3,7%** no refolded Fortran, dependendo de q. A página original não traz as matrizes e os dois mapas usados nesse número. O novo ensaio recria as duas saídas com o mesmo Hamiltoniano e registra o resultado por modo de depósito. |

No carbono mio-1-1, mediram-se normas radiais 0,9999993382 (s) e 0,9999991246 (p). São próximas da unidade. Nessa base, a mudança de normalização é pequena demais para explicar, isoladamente, o excesso percentual relatado no diamante. A comparação nativa usa a execução de bandas e a exportação H/S a partir das mesmas cargas convergidas, as mesmas imagens atômicas e os mesmos parâmetros. Os três arquivos nativos de bandas devem vir da mesma execução. O relatório `states.json` compara os autovetores nativos e reconstruídos na métrica S, verifica a diagonalidade de H e soma o peso fora de cada subespaço degenerado. Isso permite localizar uma divergência antes da projeção radial.

O depósito linear com $\sigma=0$ reparte a linha entre os dois pontos de energia adjacentes, inclusive na borda da janela. Com $\sigma>0$, ambos usam o mesmo desvio-padrão em eV. O Fortran publicado corta a gaussiana a 5σ e o Python a 8σ. A saída Fortran é peso por intervalo. A comparação divide por ΔE antes de calcular as diferenças. Relatam-se o maior erro absoluto em peso/eV, sua razão com o pico Python e a norma L1 relativa, sem renormalizar os mapas. O JSON também registra os dois pesos integrados por bins, o peso do CSV com janela rígida e os extremos da variação relativa por ponto q com sinal significativo. Esses totais permitem testar o efeito das bordas no número originalmente publicado. [Convenções do repositório](../../docs/METODO.md)

## Comparação com o executável Fortran

O GitHub Actions compilou o DFTB+ 26.1 junto com o `unfold` de `claude-fortran-unfold/` e rodou a diagonalização, a projeção e o refolding para os dois cristais. Os 11 testes Python e as quatro comparações nativas passaram. A tabela mede a diferença L1 relativa em todo o mapa, $100\sum_{q,E}|I_{\rm Fortran}-I_{\rm Python}|/\sum_{q,E}|I_{\rm Python}|$. [Execução e relatórios no GitHub Actions](https://github.com/matheusjsmatos/dftbplus-unfolding/actions/runs/45932713924)

| Material | Depósito na malha | L1 do unfolding completo | L1 do refolding |
|---|---|---:|---:|
| Grafeno 2×2 | Linear, ΔE = 0,05 eV | 0,0877% | 0,0956% |
| Grafeno 2×2 | Gaussiano, σ = 0,08 eV | 0,0222% | 0,0236% |
| Diamante, 8 átomos | Linear, ΔE = 0,10 eV | 0,0534% | 0,0515% |
| Diamante, 8 átomos | Gaussiano, σ = 0,08 eV | 0,0270% | 0,0255% |

Os autovalores reconstruídos a partir de H/S apresentaram diferenças máximas de 5,18 × 10⁻⁵ eV no grafeno e 5,07 × 10⁻⁵ eV no diamante em relação ao `band.out`. A maior componente de $C^\dagger S C-I$ em módulo foi 2,90 × 10⁻¹³ nos dois ensaios considerados em conjunto.

No diamante, a soma do peso refoldado nos 46 pontos, com depósito linear na janela de −25 a +15 eV, foi **323,17491** no Fortran e **323,17275** no Python. A diferença relativa integrada por ponto variou de −0,00224% a +0,00137%. O excesso de 1,7–3,7% do ensaio anterior não se repetiu. Para atribuir uma causa à diferença histórica seriam necessários os arquivos daquele cálculo. Esta validação cobre carbono com base s,p do mio-1-1 em grafeno ideal e diamante.

## Reproduzir a validação nativa

O aplicativo Fortran deste repositório usa as assinaturas do DFTB+ 26.1 para `printDftbHeader` e `warnUnprocessedNodes`. A tentativa de compilação com 25.1 identificou uma incompatibilidade nessas duas chamadas. O teste usa o pacote-fonte 26.1 e verifica seu SHA-256 antes de compilar. [Release 26.1](https://github.com/dftbplus/dftbplus/releases/tag/26.1)

Na pasta `chatgpt-python-unfold`, com `dftb+` e o `unfold` deste repositório compilados e o repositório mio-1-1 descompactado, execute

```bash
python tools/check_native.py graphene --dftb /caminho/dftb+ --unfold /caminho/unfold --mio /caminho/mio --output results/native
python tools/check_native.py diamond  --dftb /caminho/dftb+ --unfold /caminho/unfold --mio /caminho/mio --output results/native
```

O script verifica os hashes de `skfiles/C-C.skf` e `extras/wfc.hsd`, cria diretórios novos e faz um cálculo com depósito linear e outro com σ=0,08 eV para cada material. Para comparar uma saída já disponível, use `tools/compare_backends.py --help`. O fluxo [GitHub Actions](../../.github/workflows/validar-unfolding.yml) contém os comandos de compilação do DFTB+ 26.1 e publica os relatórios `comparison.json`, as entradas e os logs. Os critérios-padrão do ensaio são 0,003 para o erro máximo relativo ao pico e 0,002 para o L1 relativo. Uma janela vazia é registrada como erro. Os limites são critérios de reprodução numérica para estas malhas e parâmetros.

O conjunto **mio-1-1** fornece os parâmetros C–C e os radiais de carbono. O teste arquivado usou grafeno com $a=2,46$ Å, altura de 18 Å, 100 K e malhas SCC 12×12×1 e 6×6×1. O corte de 25 Ry (12,5 Ha) inclui 1.375 vetores de refolding. No diamante, a célula convencional tem 8 C, $a=3,567$ Å, caminho L–Γ–X com 46 pontos e corte de 20 Ry. A comparação contempla os casos eletrônicos escalar e colinear suportados pelo software. A validação numérica aqui usa carbono s,p. [DFTB+, Hourahine et al.](https://doi.org/10.1063/1.5143190) · [SCC-DFTB, Elstner et al.](https://doi.org/10.1103/PhysRevB.58.7260) · [Parâmetros usados](https://github.com/dftbparams/mio/tree/263e52419ce15222c864df9b3f6fbf58bfbdfa0f)
