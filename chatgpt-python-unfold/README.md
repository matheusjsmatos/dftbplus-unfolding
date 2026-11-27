# Unfolding de bandas eletrônicas do DFTB+

Este pós-processador calcula o espectro eletrônico no espaço dos momentos e o refolding na zona de Brillouin de uma célula de referência pelo método de Mayo, Yndurain e Soler. Lê as matrizes em espaço real exportadas pelo DFTB+ e as funções radiais do Waveplot. A saída é um CSV de energias e pesos por banda e spin. A execução independe de recompilar o DFTB+.

Os comandos abaixo reproduzem a cadeia não ortogonal a partir da pasta `dftb_unfold`.

```bash
python -m pip install -r requirements.txt
python dftb_unfold.py examples/chain/unfold.json -o cadeia.csv
python spectrum.py cadeia.csv --range -11.6 -3.9 --sigma 0.07 -o cadeia.npz
python -m unittest -v test_unfold
```

Para produzir também uma figura, instale `pillow reportlab` e acrescente `--pdf cadeia.pdf` ao comando `spectrum.py`. A figura desse exemplo já está em `examples/chain/spectrum.pdf`. O núcleo requer Python 3.9 ou superior e NumPy. Se estiver instalado, SciPy fornece a diagonalização generalizada por LAPACK. O caminho alternativo em NumPy usa Cholesky e diagonalização hermitiana.

## Viabilidade e correspondência entre os códigos

O artigo anexo define a densidade espectral em momento na Eq. (6), o refolding na Eq. (9) e a projeção LCAO nas Eqs. (10)–(14) da versão publicada. A inspeção da rotina `unfold.F90` do SIESTA mostra onde entram as matrizes HSX, os orbitais numéricos, as fases e a soma de intensidades. A mesma construção é viável no DFTB+ quando a parametrização disponibiliza as funções radiais das espécies usadas. [Artigo publicado](https://doi.org/10.1088/1361-648X/ab6e8e) · [Util/Unfolding do SIESTA](https://gitlab.com/siesta-project/siesta/-/tree/dev/Util/Unfolding?ref_type=heads) · [Tutorial oficial do SIESTA](https://docs.siesta-project.org/projects/siesta/en/5.0/tutorials/advanced/unfolding/index.html)

| Etapa | SIESTA | DFTB+ e implementação entregue |
|---|---|---|
| Recuperar o Hamiltoniano e o overlap | `unfold.F90` lê os blocos do HSX | `writeHSAndStop` e `writeSparse` exportam os blocos por par de átomos e translação. `RealMatrix` reconstrói o triângulo hermitiano de acordo com `unpackHS_cmplx_kpts`. |
| Obter os coeficientes LCAO | Diagonalização generalizada para cada momento necessário | Resolve $H(K)C=S(K)CE$ com $C^\dagger S C=I$. Para spin colinear, a exportação do DFTB+ converte os canais para ↑ e ↓. |
| Transformar cada orbital | `radfft` atua sobre a base numérica do SIESTA | `FourierBasis` integra os radiais de Slater com $j_l(qr)$ até o corte do Waveplot e usa os harmônicos reais na ordem do DFTB+. |
| Projetar e somar | A rotina combina a fase, o fator angular e os coeficientes, toma o módulo ao quadrado e acumula no refolding | O CSV preserva o peso no momento especificado e o peso somado sobre vetores recíprocos da referência. O mapa aplica alargamento em energia posteriormente. |

A revisão analisada do DFTB+ é `fedabbcb993ef26967141f275bb7d31a984e4bdc`. No SIESTA, a rotina consultada pertence à árvore `dev` indicada pelo usuário. [Exportação DFTB+](https://github.com/dftbplus/dftbplus/blob/fedabbcb993ef26967141f275bb7d31a984e4bdc/src/dftbp/dftbplus/mainio.F90) · [Formato dos blocos](https://github.com/dftbplus/dftbplus/blob/fedabbcb993ef26967141f275bb7d31a984e4bdc/src/dftbp/io/formatout.F90) · [Montagem de Bloch](https://github.com/dftbplus/dftbplus/blob/fedabbcb993ef26967141f275bb7d31a984e4bdc/src/dftbp/dftb/sparse2dense.F90) · [Rotina de referência SIESTA](https://gitlab.com/siesta-project/siesta/-/blob/dffb40802f327353576f222ab8a10f9a8c223c4b/Util/Unfolding/Src/unfold.F90)

## Como preparar um cálculo DFTB+

1. Convirja as cargas SCC da estrutura na geometria final, com uma malha de integração adequada. No cálculo eletrônico estático que gera `charges.bin`, habilite as opções abaixo. O `detailed.xml` guarda as coordenadas das imagens centrais efetivamente selecionadas pelo DFTB+. Essa escolha importa para as fases dos átomos próximos às faces da supercélula.

```hsd
Options {
  WriteCharges = Yes
  WriteDetailedXML = Yes
  ShowFoldedCoords = Yes
}
```

2. Em outro diretório, copie `charges.bin`. Use a mesma geometria, ordenação dos átomos, base, parâmetros SK, spin e demais termos do Hamiltoniano. No bloco `Hamiltonian = DFTB { ... }` mantenha `SCC = Yes` e defina `ReadInitialCharges = Yes` e `RecomputeAfterDensity = No`. No bloco `Options` defina `WriteRealHS = Yes`. Se o SCC gravou cargas em texto com `Options { WriteChargesAsText = Yes }`, copie `charges.dat` e habilite `ReadChargesAsText = Yes` em `Options` na segunda execução. O DFTB+ grava as matrizes e encerra antes de iterar o SCC. O marcador de conclusão é `Hamilton/Overlap written, exiting program.`. Para cálculos originalmente sem SCC, conserve `SCC = No` e exporte as matrizes diretamente. [Manual DFTB+](https://github.com/dftbplus/dftbplus/blob/fedabbcb993ef26967141f275bb7d31a984e4bdc/doc/dftb%2B/manual/dftbp.tex) · [Receita de bandas com cargas convergidas](https://dftbplus-recipes.readthedocs.io/en/latest/basics/bandstruct.html)

3. Reúna os arquivos abaixo. As funções radiais devem corresponder à parametrização Slater–Koster usada nas matrizes. A compatibilidade do número de orbitais é verificada automaticamente. A identidade física da parametrização deve ser confirmada na documentação da base. Para grafeno com B substitucional, por exemplo, são necessários os integrais de pares compatíveis B–B, B–C, C–B e C–C utilizados pelo DFTB+ e os radiais de B e C. Os arquivos `.skf` fornecem integrais entre orbitais. A transformada de Fourier individual de cada orbital requer também sua forma radial. [Documentação do Waveplot](https://dftbplus-recipes.readthedocs.io/en/latest/basics/waveplot.html)

| Entrada | Conteúdo e unidades |
|---|---|
| `hamreal1.dat` e `overreal.dat` | Blocos esparsos do DFTB+, Hamiltoniano em hartree e overlap adimensional. Para dois canais colineares, inclua também `hamreal2.dat` da mesma execução. |
| `detailed.xml` ou geometria `.gen` | Supercélula com átomos na ordem da matriz. Prefira o XML do cálculo estático com `ShowFoldedCoords = Yes`. GEN periódico aceita `S`, cartesiano em Å, ou `F`, fracionário. Nos exemplos GEN, a origem é zero. |
| `wfc.*.hsd` | Definições numéricas literais de orbitais do Waveplot, nas unidades do arquivo original. Aceita arquivos separados por espécie e a sequência de camadas $s,p,d,f$ truncada no momento angular usado pelo DFTB+. Passe diretamente os arquivos de coeficientes. |
| `unfold.json` | Célula de referência, caminho recíproco, corte para refolding, referência de energia e localização dos arquivos. Os caminhos dos insumos são relativos ao JSON. |

4. Ajuste uma cópia de `examples/graphene/unfold.json`. Use `"geometry": "scf/detailed.xml"` e `"matrix_dir": "export"`, por exemplo, quando esses diretórios estiverem ao lado do JSON. Para uma supercélula construída com matriz inteira $M$, a convenção é $A_{\rm SC}=M A_{\rm ref}$ com **vetores nas linhas**. `supercell_matrix` aceita $M$ não diagonal. A alternativa `reference_cell_angstrom` especifica os três vetores de referência explicitamente. `path_frac` contém vértices em coordenadas fracionárias recíprocas da **referência**. `intervals` contém o número de intervalos de cada segmento e cada vértice compartilhado aparece uma única vez.

```bash
python dftb_unfold.py unfold.json --plan
python dftb_unfold.py unfold.json -o bandas.csv
python spectrum.py bandas.csv --range -4 4 --sigma 0.05 --weight refold --spin 1 -o mapa.npz --pdf mapa.pdf
```

`reference_ev` é subtraído das energias em eV. Para representar $E-E_F$, introduza o Fermi do SCC convergido e use a mesma referência nas comparações. `window_ev` é opcional e limita as bandas gravadas, em relação a essa origem. Para produzir mapas com todos os pontos do caminho, pode-se omitir `window_ev` e cortar apenas em `spectrum.py`. O código calcula o espectro para os estados ocupados e vazios, com pesos por canal, sem fator de ocupação nem fator 2 de spin.

## O que significam os pesos

Para $\mathbf q=\mathbf K+\mathbf G_{\rm SC}$, usando ondas planas normalizadas no volume $V_{\rm SC}$, a amplitude empregada é

$$
a_{n\mathbf K}(\mathbf q)=\frac{4\pi}{\sqrt{V_{\rm SC}}}
\sum_\mu C_{\mu n}(\mathbf K)\,e^{-i\mathbf q\cdot\boldsymbol\tau_\mu}
(-i)^{l_\mu}Y^{\rm real}_{l_\mu m_\mu}(\widehat{\mathbf q})
\int_0^{r_{c,\mu}}r^2 j_{l_\mu}(qr)R_\mu(r)\,dr .
$$

O artigo e o utilitário SIESTA incluem a posição do centro orbital dentro da célula na escolha de fase dos coeficientes de Bloch. O Waveplot do DFTB+ usa uma soma de translações com $\exp(-i\mathbf K_{\rm DFTB+}\cdot\mathbf R)$. No pós-processador o índice CSV é o momento físico, definido por $\psi(\mathbf r+\mathbf R)=e^{+i\mathbf K\cdot\mathbf R}\psi(\mathbf r)$. Com a transposição dos blocos do exportador, isso equivale a chamar a montagem nativa em $-\mathbf K$ e usar $e^{-i\mathbf q\cdot\boldsymbol\tau_\mu}$ na projeção. Para $p$, por exemplo, a ordem interna é $p_y,p_z,p_x$. Essas convenções alteram a interferência entre sítios e orbitais. [Soma de Bloch do Waveplot](https://github.com/dftbplus/dftbplus/blob/fedabbcb993ef26967141f275bb7d31a984e4bdc/app/waveplot/molorb.F90) · [Radiais e harmônicos reais](https://github.com/dftbplus/dftbplus/blob/fedabbcb993ef26967141f275bb7d31a984e4bdc/app/waveplot/slater.F90)

| Campo | Interpretação |
|---|---|
| `point`, `distance_inv_ang`, `q1,q2,q3` | Índice do ponto do caminho a partir de 1, distância acumulada em Å⁻¹ e coordenadas recíprocas fracionárias da referência. |
| `K1_sc,K2_sc,K3_sc`, `spin`, `band` | $\mathbf K$ físico reduzido módulo a rede recíproca da supercélula, em coordenadas fracionárias centradas em zero. Índices de spin e banda começam em 1. |
| `energy_ev`, `energy_rel_ev` | Autovalor em eV e autovalor menos `reference_ev`. |
| `weight_full` | $|a_{n\mathbf K}(\mathbf q)|^2$, peso adimensional da onda plana no momento indicado. |
| `weight_refold` | Soma **incoerente** dos pesos em $\mathbf q+\mathbf g_{\rm ref}$ para $|\mathbf g_{\rm ref}|^2<\texttt{gcut_ry}$ em unidades atômicas. O termo $\mathbf g=0$ é sempre incluído. Corte zero produz apenas `weight_full`. |

Numa referência incomensurável, diferentes $\mathbf q+\mathbf g_{\rm ref}$ podem corresponder a diferentes $\mathbf K$. Nessa situação há múltiplos blocos de bandas para o mesmo `point` e cada linha carrega o seu próprio $\mathbf K$, energia e contribuição. Somar os pesos **a cada energia** produz o espectro do artigo. `spectrum.py` faz essa operação usando gaussianas de área unitária e largura `--sigma` (desvio padrão, em eV). O NPZ guarda a densidade em eV⁻¹ por ponto do caminho, o eixo de energia e a distância. Em degenerescências, compare a soma dos pesos do subespaço degenerado.

**Grafeno e sistemas 2D.** Para grafeno plano, o orbital $p_z$ é ímpar em $z$. Sua transformada se anula em $q_z=0$. As bandas $\pi$ podem, portanto, desaparecer na coluna `weight_full` ao longo de $\Gamma$–K–M no plano. O refolding tridimensional inclui termos com $g_z\ne0$ e recupera esse peso mesmo quando o caminho tem $q_z=0$. Mantenha o mesmo vácuo nas células de referência e da supercélula. Verifique a estabilidade das intensidades aumentando `gcut_ry` e, para comparar sistemas, mantenha a mesma convenção de normalização. A magnitude desses pesos representa uma densidade espectral de momento. Uma intensidade fotoemissiva experimental exigiria também os elementos de matriz da medida. [Artigo](https://doi.org/10.1088/1361-648X/ab6e8e)

## Validação e limites

| Verificação | Resultado |
|---|---|
| Cadeia periódica homogênea dobrada, com overlap de primeiros vizinhos $s=0{,}07097404418511578$ | Em 81 pontos e a 25 Ry, erro máximo de $3{,}6\times10^{-15}$ eV na dispersão permitida. Peso máximo na réplica energeticamente distinta de $5{,}1\times10^{-30}$. Peso refoldado da banda física entre 0,999752 e 0,999930, somando degenerescências. O modelo analítico é inteiramente incluído. |
| Cadeia com perturbação em um único sítio | O teste detecta peso nas duas bandas, verifica $C^\dagger S C=I$, invariância de fase global e a soma de peso entre momentos parceiros. |
| Matrizes de H₂O exportadas pelo DFTB+ no conjunto público de testes do `sisl` | Leitura independente dos blocos com seis orbitais. Na indexação Python, $H[0,4]=-0{,}2885462383598692$ hartree e $S[0,4]=0{,}2439705710593631$. O maior desvio elemento a elemento de $C^\dagger S C-I$ foi $1{,}1\times10^{-15}$. Este teste verifica a leitura e a diagonalização molecular. [Matrizes originais](https://github.com/zerothi/sisl-files/tree/ed8c733c3b71c12b37aba258c169cdeeb665d7e0/tests/dftb/h2o) |
| Projeção e formatos | Os testes incluídos verificam a transformada de um estado com interferência complexa $s/p$ por integração direta, Bessel e harmônicos até $f$, supercélula oblíqua, referência incomensurável, dois canais colineares e leitura da geometria periódica. |
| Grafeno cristalino $2\times2$ | O caso de integração abaixo está preparado. O DFTB+ executável não estava disponível no ambiente de desenvolvimento, de modo que a comparação periódica com `band.out` permanece pendente. |

Na cadeia, a solução de referência é $E(k)=[-0{,}30-0{,}16\cos(2\pi k)]/[1+2s\cos(2\pi k)]$ hartree. A sobreposição decorre do mesmo orbital compacto usado na projeção. A soma de pesos a corte finito mede também a convergência da expansão em ondas planas. Em bases realistas, há ainda a qualidade do ajuste das funções radiais do Waveplot em relação às integrais usadas no SK. O programa preserva os pesos absolutos calculados.

### Comparação reproduzível no grafeno

O comando abaixo requer `dftb+` no PATH, um conjunto contendo `C-C.skf` e o arquivo radial correspondente.

```bash
python examples/graphene/validate_graphene.py --sk-dir /caminho/para/sk --wfc /caminho/para/wfc.hsd --workdir resultados
```

O script cria uma pasta nova, converge cargas em malhas equivalentes $12\times12\times1$ e $6\times6\times1$, executa bandas nativas de referência, exporta $H$ e $S$ e percorre $\Gamma$–K–M–$\Gamma$ em 71 pontos. Grava os `dftb_in.hsd`, logs, `unfolded.csv` e `validation.json`. A comparação em energia usa $2\times10^{-4}$ eV de tolerância, compatível com as quatro casas decimais do `band.out`. Também verifica se cada energia da primitiva encontra peso $>10^{-6}$ e se o peso fora das bandas primitivas fica abaixo de $10^{-6}$. Use o cristal puro para esse critério. Quebrar a simetria translacional por dopagem ou defeitos pode gerar peso físico nas réplicas. Os arquivos de parametrização para o grafeno devem ser fornecidos pelo usuário.

O núcleo atual cobre DFTB periódico com Hamiltoniano real em espaço real, base mínima de uma camada radial por $l$, escalar ou spin colinear. SOC, spin não colinear, transporte com contatos, células helicoidais e fônons ficam fora do escopo desta versão. Modelos com termos dependentes de $\mathbf K$ não recuperáveis desses blocos exigem outra interface. Cada par $(\mathbf K,\text{spin})$ requer uma diagonalização densa, com custo de ordem $N_{\rm orb}^3$ e memória de ordem $N_{\rm orb}^2$. `--plan` informa antecipadamente quantos momentos distintos e quantos vetores de refolding serão usados. `q_batch` limita o número de vetores na projeção simultânea. Para uma referência incomensurável o número de diagonalizações pode crescer com o corte. A janela em energia reduz a projeção e o tamanho do CSV. As rotinas de diagonalização desta versão calculam o espectro completo.
