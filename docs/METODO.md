# Fundamentação teórica do método de unfolding

Este documento deduz as equações usadas nas duas implementações do repositório (Fortran e Python) e mostra como elas se relacionam com a implementação original do SIESTA. As equações seguem Mayo, Yndurain e Soler [1], citado aqui como "o artigo". Os números entre colchetes remetem à lista de referências no fim, que é a mesma de `REFERENCIAS.md`. As equações do artigo são citadas como (A-n).

Os trechos marcados como **"nossa dedução"** ou **"verificado numericamente"** indicam o que foi feito neste projeto e não vem diretamente do artigo.

## 1. Problema

Um cálculo de supercélula (SC), com vetores $\mathbf A_i$ e volume $V$, produz estados de Bloch $\psi_{\mathbf K,i}$ com $\mathbf K$ na zona de Brillouin da supercélula (SBZ). Como a SBZ é pequena, as bandas ficam "dobradas" e difíceis de comparar com as do cristal primitivo. O unfolding recupera a informação perdida em dois passos [1]:

1. **Unfolding completo:** distribuir o peso de cada estado em momentos $\mathbf q$ de todo o espaço recíproco.
2. **Refolding:** somar esses pesos na zona de Brillouin de uma célula de referência (RBZ), obtendo as "bandas desdobradas" convencionais.

## 2. Unfolding completo (qLDOS)

### 2.1 Densidade de estados no espaço dos momentos

A densidade de estados na SBZ é (A-1):

$$ n_{\rm SBZ}(\mathbf K,E)=\sum_i \delta(E-\varepsilon_{\mathbf K,i}). $$

Cada estado de Bloch, $\psi_{\mathbf K,i}(\mathbf r)=u_{\mathbf K,i}(\mathbf r)\,e^{i\mathbf K\cdot\mathbf r}$ com $u$ periódica na SC, é expandido em ondas planas da rede recíproca da SC ($\mathbf G$):

$$ u_{\mathbf K,i}(\mathbf r)=\sum_{\mathbf G}\tilde u_{\mathbf K,i,\mathbf G}\,e^{i\mathbf G\cdot\mathbf r},\qquad
\tilde u_{\mathbf K,i,\mathbf G}=\frac1V\int_{\rm SC}u_{\mathbf K,i}(\mathbf r)\,e^{-i\mathbf G\cdot\mathbf r}\,d^3r \quad\text{(A-5)}. $$

Logo $\psi_{\mathbf K,i}=\sum_{\mathbf G}\tilde u_{\mathbf K,i,\mathbf G}\,e^{i(\mathbf K+\mathbf G)\cdot\mathbf r}$ só tem componentes de Fourier em $\mathbf q=\mathbf K+\mathbf G$. Com a transformada unitária $\tilde\psi(\mathbf q)=(2\pi)^{-3/2}\int\psi\,e^{-i\mathbf q\cdot\mathbf r}d^3r$ e $\int e^{i(\mathbf q'-\mathbf q)\cdot\mathbf r}d^3r=(2\pi)^3\delta(\mathbf q'-\mathbf q)$ obtém-se (A-4):

$$ \tilde\psi_{\mathbf K,i}(\mathbf q)=(2\pi)^{3/2}\sum_{\mathbf G}\delta(\mathbf K+\mathbf G-\mathbf q)\,\tilde u_{\mathbf K,i,\mathbf G}. $$

$|\tilde\psi(\mathbf q)|^2$ é a probabilidade de medir o momento $\mathbf q$. Somando sobre estados e energias define-se o **espectro desdobrado completo**, equivalente à densidade local de estados no espaço recíproco (qLDOS) (A-2, A-6):

$$ n(\mathbf q,E)=\sum_i\int_{\rm SBZ}d\mathbf K\,|\tilde\psi_{\mathbf K,i}(\mathbf q)|^2\,\delta(E-\varepsilon_{\mathbf K,i})
=\sum_i|\tilde u_{\mathbf K_q,i,\mathbf G_q}|^2\,\delta(E-\varepsilon_{\mathbf K_q,i}), $$

onde $\mathbf K_q$ e $\mathbf G_q$ são os únicos vetores tais que $\mathbf K_q\in{\rm SBZ}$, $\mathbf G_q$ é vetor recíproco da SC e $\mathbf K_q+\mathbf G_q=\mathbf q$. Observações:

- $n(\mathbf q,E)$ **não é periódico** em $\mathbf q$. As energias em que ele é não nulo são periódicas, mas os pesos mudam de zona para zona.
- Pelo teorema de Parseval, $\sum_{\mathbf G}|\tilde u_{\mathbf K,i,\mathbf G}|^2=\frac1V\int_{\rm SC}|u|^2=1$ (com a normalização do artigo). Assim $\sum_{\mathbf G}n(\mathbf K+\mathbf G,E)=n_{\rm SBZ}(\mathbf K,E)$: dobrar o espectro completo de volta para a SBZ recupera as bandas originais.

### 2.2 Refolding

Dada uma célula de referência com vetores recíprocos $\mathbf g$ (A-9):

$$ n_{\rm RBZ}(\mathbf k,E)=\sum_{\mathbf g}n(\mathbf k+\mathbf g,E). $$

**Caso ideal (verificação analítica).** Se o sistema é um cristal perfeito e a SC é um múltiplo da célula primitiva, um estado com vetor de onda primitivo $\mathbf k_0$ contém apenas componentes $\mathbf q=\mathbf k_0+\mathbf g$. Então $\tilde u_{\mathbf G}\neq0$ só para esses $\mathbf q$, e o refolding com a célula primitiva devolve peso total 1 exatamente em $\mathbf k=\mathbf k_0$. Defeitos e desordem espalham o peso, e o refolding mostra esse espalhamento como alargamento das bandas.

**Células incomensuráveis.** A SC não precisa ser múltipla da célula de referência (amorfos, defeitos com forte deformação). Nesse caso o par $(\mathbf K_q,\mathbf G_q)$ depende de $\mathbf g$ e é preciso uma diagonalização por $\mathbf k+\mathbf g$ [1, nota 2]. É o que o SIESTA faz com a variável `notSuperCell`. Nossas duas implementações tratam o caso comensurável. O código Python também detecta o caso incomensurável (`commensurate: False`).

## 3. Base LCAO: dedução de (A-11) e (A-14)

### 3.1 Estados na base atômica

Com orbitais $\phi_\mu$ centrados em $\mathbf R+\mathbf r_\mu$ (A-10):

$$ \psi_{\mathbf K,i}(\mathbf r)=\sum_{\mathbf R}\sum_\mu c_{\mathbf K,i,\mu}\,\phi_\mu(\mathbf r-\mathbf R-\mathbf r_\mu)\,e^{i\mathbf K\cdot(\mathbf R+\mathbf r_\mu)}. $$

Escreva $\int_{\rm all}\psi\,e^{-i\mathbf q\cdot\mathbf r}d^3r$ de duas formas. Pela expansão acima, só os termos com $\mathbf q=\mathbf K+\mathbf G$ sobrevivem na soma sobre $\mathbf R$, que contribui com um fator $N_{\rm SC}$ que se cancela na normalização por SC. Em seguida desloca-se a origem do integrando por $\mathbf R+\mathbf r_\mu$:

$$ \int\phi_\mu(\mathbf r-\mathbf R-\mathbf r_\mu)\,e^{-i\mathbf q\cdot\mathbf r}d^3r=(2\pi)^{3/2}\,\tilde\phi_\mu(\mathbf q)\,e^{-i\mathbf q\cdot(\mathbf R+\mathbf r_\mu)}. $$

Os fatores $e^{i\mathbf K\cdot\mathbf R}e^{-i\mathbf q\cdot\mathbf R}=e^{-i\mathbf G\cdot\mathbf R}=1$ e, para a posição do átomo, $e^{i\mathbf K\cdot\mathbf r_\mu}e^{-i\mathbf q\cdot\mathbf r_\mu}=e^{-i\mathbf G\cdot\mathbf r_\mu}$. O resultado é (A-11):

$$ \tilde u_{\mathbf K,i,\mathbf G}=\frac{(2\pi)^{3/2}}{V}\sum_\mu c_{\mathbf K,i,\mu}\,\tilde\phi_\mu(\mathbf K+\mathbf G)\,e^{-i\mathbf G\cdot\mathbf r_\mu}. $$

### 3.2 Transformada de um orbital atômico (nossa dedução de A-14)

Seja $\phi_\mu(\mathbf r)=R_\mu(r)\,Y_{l m}(\hat r)$ com harmônicas esféricas **reais** (as do DFTB+). A expansão da onda plana é

$$ e^{-i\mathbf q\cdot\mathbf r}=4\pi\sum_{l,m}(-i)^l\,j_l(qr)\,Y_{lm}(\hat q)\,Y_{lm}(\hat r), $$

e a ortogonalidade das $Y_{lm}$ reais dá

$$ \tilde\phi_\mu(\mathbf q)=\frac{1}{(2\pi)^{3/2}}\int\phi_\mu e^{-i\mathbf q\cdot\mathbf r}d^3r
=(-i)^{l}\,Y_{lm}(\hat q)\,F_l(q),\qquad F_l(q)=\sqrt{\frac2\pi}\int_0^\infty r^2R_\mu(r)\,j_l(qr)\,dr. $$

O fator $\sqrt{2/\pi}$ vem de $4\pi/(2\pi)^{3/2}$. É a (A-14) do artigo, que usa $(-i)^l$, e o SIESTA também (`phi = (-ii)**l * ylm * phi` em `unfold.F90`).

### 3.3 Peso espectral final

Combinando 3.1 e 3.2, com $\mathbf q=\mathbf K+\mathbf G$, e multiplicando por $V$ para que a regra de soma valha com $c^\dagger Sc=1$ (ver 5):

$$ \boxed{\,w_i(\mathbf q)=V\,|\tilde u_{\mathbf K,i,\mathbf G}|^2=\frac{(2\pi)^3}{V}\Big|\sum_\mu c_{\mathbf K,i,\mu}\,(-i)^{l_\mu}\,Y_{l_\mu m_\mu}(\hat q)\,F_{l_\mu}(q)\,e^{-i\mathbf G\cdot\mathbf r_\mu}\Big|^2\,} $$

$$ n(\mathbf q,E)=\sum_i w_i(\mathbf q)\,\delta(E-\varepsilon_{\mathbf K_q,i}). $$

No SIESTA, $c_0=(2\pi)^{3/2}/V$ e `ddos = cdos*vol*abs(ukg)**2` com `cdos = 1`, ou seja, exatamente $w=V|c_0\sum\ldots|^2$.

## 4. Convenções do DFTB+ (nossa dedução e verificação)

### 4.1 Gauge de célula

O SIESTA usa a fase $e^{i\mathbf K\cdot\mathbf r_\mu}$ nos coeficientes (A-10), com as fases do vetor interatômico em $H(\mathbf K)$. O DFTB+ monta $H(\mathbf k)$ com fase apenas nas translações de célula $\mathbf R$ (gauge de célula). Os coeficientes dos dois gauges se relacionam por

$$ \tilde c_\mu=c_\mu\,e^{i\mathbf K\cdot\mathbf r_\mu}. $$

Substituindo em 3.1, a fase $e^{-i\mathbf G\cdot\mathbf r_\mu}$ vira $e^{-i\mathbf q\cdot\mathbf r_\mu}$ com $\mathbf q=\mathbf K+\mathbf G$ completo:

$$ w_i(\mathbf q)=\frac{(2\pi)^3}{V}\Big|\sum_\mu\tilde c_{\mathbf K,i,\mu}\,(-i)^{l_\mu}Y_{l_\mu m_\mu}(\hat q)\,F_{l_\mu}(q)\,e^{-i\mathbf q\cdot\mathbf r_\mu}\Big|^2. $$

É esta forma que as duas implementações usam, porque os autovetores do DFTB+ estão no gauge de célula.

### 4.2 Sinal do vetor de onda

O DFTB+ usa o fator de fase $e^{+i\mathbf k\cdot\mathbf R}$ com a convenção de $\mathbf R$ própria do código. O vetor de Bloch físico, no sentido da seção 3, é $\mathbf K=-\mathbf k_{\rm arquivo}$. Isso foi **verificado numericamente**: em k-pontos genéricos, a regra de soma só dá ≈ 1 com `KPointSign = -1` (ação `CheckBasis`), e o TRIM, por simetria, não distingue o sinal. O programa Fortran tem `KPointSign` configurável; no Python a convenção está dentro de `RealMatrix.bloch`, e foi verificada pelos testes unitários e pela concordância do peso desdobrado com o Fortran.

### 4.3 Base radial e harmônicas

- A parte radial vem de `wfc.*.hsd` (formato do Waveplot): $R(r)=r^{l}\sum_{p,a}c_{pa}\,r^{p-1}e^{-\alpha_a r}$ até o raio de corte $r_c$.
- Harmônicas reais na ordem do DFTB+/Waveplot: s; p = (y, z, x); d = (xy, yz, 3z²−r², xz, x²−y²); até $l=3$.
- **A base do `wfc` deve ser a mesma usada para gerar os arquivos Slater–Koster.** Caso contrário, os orbitais da transformada não são os orbitais cujo overlap está em $S$, e a regra de soma falha (ver 6).

## 5. Regra de soma

Pela relação de Parseval, para estados com $c^\dagger S_{\rm STO}\,c=1$:

$$ \sum_{\mathbf G}w_i(\mathbf K+\mathbf G)=c_i^\dagger S_{\rm STO}\,c_i, $$

onde $S_{\rm STO}$ é o overlap calculado com as STOs do `wfc`. Se o overlap do cálculo for o mesmo (mesmo conjunto de parâmetros SK), a soma dá ≈ 1. Desvios medem (i) base inconsistente ou (ii) faixa de $q$ insuficiente. Usamos isso como teste (`SumRuleCheck`) e como critério de aceitação: 0,976–0,999 com mio-1-1 nos dois lados, 0,4–3,8 com SK 3ob e STO mio.

## 6. Discretização numérica

1. **Conjunto de $\mathbf g$:** os vetores recíprocos da célula de referência com $|\mathbf g|^2<g_{\rm cut}$. O SIESTA usa `RefoldingGcutoff` em Ry (padrão 50 Ry) com $G^2<$ corte. No DFTB+ (Hartree), $|\mathbf g|^2/2\le g_{\rm cut}$, o que é a mesma condição física. Escolha igual nos dois códigos: 113 vetores para o diamante a 20 Ry.
2. **Mapeamento $\mathbf q\to\mathbf K$:** em coordenadas fracionárias da SC, $x=\mathbf q\cdot\mathbf A/2\pi$, $\mathbf K=\mathbf B\,(x-\mathrm{nint}\,x)$, $\mathbf G=\mathbf q-\mathbf K$.
3. **$F_l(q)$:** o SIESTA tabela por `radfft` em malha com $q_c=50\ \mathrm{bohr}^{-1}$ (`dr = π/qc`). No Fortran daqui, tabela com $dq=0{,}01$, Simpson em 6001 pontos radiais e interpolação linear. No Python, quadratura de Gauss–Legendre de 12 pontos por painel, sem tabela. O truncamento brusco em $r_c$ faz $F_l(q)$ decair lentamente em $q$ grande, e isso limita a precisão dos pesos em $|\mathbf q+\mathbf g|$ alto.
4. **Energia:** deposição linear em dois canais vizinhos da malha em energia (igual ao SIESTA, com os pesos `we` e `1-we`). O Fortran também tem alargamento gaussiano opcional.
5. **Duas saídas:** `unfoldedBands` (apenas $\mathbf g=0$) e `refoldedBands` (soma sobre $\mathbf g$), no formato de arquivos do SIESTA.

## 7. A implementação original do SIESTA (Mayo & Soler)

O utilitário `Util/Unfolding` do SIESTA [3, 4] é a referência. Foi escrito por S. G. Mayo e J. M. Soler (2018) e depois paralelizado por P. Ordejón (2021). A estrutura do `unfold.F90` é:

1. Lê o arquivo `HSX` com $H$, $S$ e os vetores interatômicos $\mathbf x_{ij}$. O `SaveHS` do SIESTA o escreve.
2. Tabela os orbitais numéricos em $r$ e calcula a transformada $\tilde\phi(q)$ por `radfft` ($q_c=50\ \mathrm{bohr}^{-1}$).
3. Para cada ponto $\mathbf q$ do caminho e cada $\mathbf g$ de `RefoldingLatticeVectors`:
   - Obtém $\mathbf K_q$ e $\mathbf G_q$ da forma de 6.2. Se a SC é supercélula da referência, diagonaliza apenas para $\mathbf g=0$.
   - Monta $H(\mathbf K_q)=\sum e^{-i\mathbf K\cdot\mathbf x_{ij}}H_{ij}$ e idem $S$, e resolve o problema generalizado (`cdiag`).
   - Calcula $u_{\mathbf K G}=c_0\sum_\mu c_\mu\,\tilde\phi_\mu(|\mathbf q+\mathbf g|)\,(-i)^l\,Y_{lm}(\widehat{q+g})\,e^{-i\mathbf G_q\cdot\mathbf r_\mu}$ (funções `rlylm` para $Y_{lm}$).
   - Acumula `ddos = vol·|u|²` na malha de energia: em `udos` apenas quando $|\mathbf g|\approx0$, em `rdos` sempre.
4. Escreve `.unfoldedBands` e `.refoldedBands`. O formato é o mesmo que usamos.

Diferença de gauge: o SIESTA usa $e^{-i\mathbf K\cdot\mathbf x_{ij}}$ com $\mathbf x_{ij}$ interatômico, o que corresponde ao gauge atômico de (A-10), e por isso só a fase $e^{-i\mathbf G\cdot\mathbf r_\mu}$ aparece. No DFTB+ o gauge é de célula (seção 4.1).

## 8. Como as duas implementações deste repositório se encaixam

| Etapa | SIESTA [3, 4] | Fortran (Claude) | Python (ChatGPT) |
|---|---|---|---|
| Dados do cálculo | `HSX` ($H$, $S$ em espaço real) | `eigenvec.bin`, `band.out`, `detailed.xml` | `hamreal1.dat`, `overreal.dat` (`WriteRealHS`) |
| Autovetores | diagonaliza em cada $\mathbf K_q$ | usa os do DFTB+ nos k-pontos listados por `ListKPoints` | diagonaliza $H(\mathbf K)$, $S(\mathbf K)$ em $\mathbf K$ arbitrário (Cholesky/SciPy) |
| Gauge | atômico | célula, com $\mathbf K=-\mathbf k$ | célula |
| $\tilde\phi(q)$ | `radfft`, malha | tabela + Simpson + interpolação | Gauss–Legendre, sem tabela |
| Base | orbitais numéricos | STO do `wfc` | STO do `wfc` |
| Fases | $e^{-i\mathbf G\cdot\mathbf r_\mu}$ | $e^{-i\mathbf q\cdot\mathbf r_\mu}$ | $e^{-i\mathbf q\cdot\mathbf r_\mu}$ |
| Saída | `.unfoldedBands`, `.refoldedBands` | mesmo formato | CSV (peso completo e refolded por banda) |

As três implementações calculam a mesma expressão (seção 3.3). As diferenças são de entrada, gauge, quadratura radial e fluxo de trabalho. A comparação numérica entre as duas implementações deste repositório está em `COMPARACAO.md`.

## 9. Limites do que foi validado

- Testado apenas em cálculos sem polarização de spin. Spin colinear polarizado e não colinear não foram validados.
- Orbitais até $l=3$ (limite das harmônicas implementadas nas duas leituras).
- Validado em diamante (supercélula de 8 átomos) e grafeno 3×3, com SK mio-1-1 e STO do mesmo conjunto.
- Células incomensuráveis: o Python as detecta, mas o fluxo completo nelas não foi validado aqui.
- O peso refolded das duas implementações difere em poucos por cento; a causa não foi identificada (`COMPARACAO.md`).

## Referências

[1] S. G. Mayo, F. Yndurain, J. M. Soler, "Band unfolding made simple", J. Phys.: Condens. Matter **32**, 205902 (2020). doi:10.1088/1361-648X/ab6e8e.
[2] S. G. Mayo, F. Yndurain, J. M. Soler, arXiv:1812.03925 (versão em preprint do artigo [1]).
[3] J. M. Soler, E. Artacho, J. D. Gale, A. García, J. Junquera, P. Ordejón, D. Sánchez-Portal, "The SIESTA method for ab initio order-N materials simulation", J. Phys.: Condens. Matter **14**, 2745 (2002).
[4] SIESTA, `Util/Unfolding` (README e `Src/unfold.F90`), S. G. Mayo e J. M. Soler (2018), com modificações de P. Ordejón (2021). Código e README consultados na cópia fornecida pelo autor deste repositório.
[5] T. Kosugi, H. Nishi, Y. Kato, Y. Matsushita, J. Phys. Soc. Jpn. **86**, 124717 (2017), citado em [1] como origem do espectro desdobrado em ondas planas.
[6] B. Hourahine et al., "DFTB+, a software package for efficient approximate density functional theory based atomistic simulations", J. Chem. Phys. **152**, 124101 (2020).
[7] B. Aradi, B. Hourahine, Th. Frauenheim, "DFTB+, a sparse matrix-based implementation of the DFTB method", J. Phys. Chem. A **111**, 5678 (2007).
[8] M. Elstner et al., "Self-consistent-charge density-functional tight-binding method for simulations of complex materials properties", Phys. Rev. B **58**, 7260 (1998). Parâmetros mio-1-1.
