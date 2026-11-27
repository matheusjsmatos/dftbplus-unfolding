# Testes de unfolding no DFTB+

Na raiz do repositório, instale o NumPy no ambiente Python e rode a suíte com os dados arquivados do grafeno.

```bash
python -m pip install -r chatgpt-python-unfold/requirements.txt
python Testes/run.py
```

A execução funciona sem acesso à rede depois da instalação e usa somente o Python e as entradas do repositório. O comando termina com código zero quando todas as etapas passam.

| Etapa | O que se verifica |
|---|---|
| Oito testes em `test_band_unfolding.py` | Fases de Bloch, interferência entre orbitais, autovetores na métrica S, cadeia duplicada, imagens periódicas e entradas inválidas. |
| Três testes em `test_spectral_maps.py` | Normalização por intervalo de energia, depósito linear em dois pontos da malha e detecção de mapas incompatíveis. |
| Grafeno puro 2×2 | Compara energias da supercélula com `band.out` e verifica o peso refoldado nas bandas da célula primitiva ao longo de Γ–K–M–Γ. |
| Projeção radial do carbono | Recalcula em NumPy a quadratura de Simpson e os pesos absolutos com os mesmos autovetores H/S. A comparação com o binário Fortran está no modo `native`. |

A pasta de saída padrão é `results/testes/`, relativa ao diretório de execução. Nela, `unfolded.csv` contém os pesos por estado, `validation.json` traz a comparação com as bandas do DFTB+ e `projection.json` traz os erros da projeção radial. Para guardar outra execução, use `python Testes/run.py offline --output results/repeticao`. Para executar somente os testes unitários, use `python -m unittest discover -s Testes -v`.

## Comparação com o executável Fortran

Com os executáveis `dftb+` e `unfold` compilados com DFTB+ 26.1 e o diretório do mio-1-1 usado no projeto, execute a partir da raiz. Substitua os três caminhos pelos locais da instalação.

```bash
python Testes/run.py native --dftb /caminho/dftb+ --unfold /caminho/unfold --mio /caminho/mio --output results/native
```

O programa verifica os hashes de `skfiles/C-C.skf` e `extras/wfc.hsd` do mio-1-1 e calcula o grafeno 2×2 e o diamante de oito átomos em subpastas novas. Para cada material, o teste confronta as bandas primitivas e os autovetores nativos com a reconstrução H/S. Também compara os dois mapas espectrais com depósito linear e com alargamento gaussiano de 0,08 eV. Uma falha em qualquer material resulta em código de saída diferente de zero.

Em `results/native/`, procure `validation.json` para as bandas, `states.json` para os autovetores e `comparison.json` nas subpastas `linear/` e `gaussian/` para os pesos. Os limites por métrica, os parâmetros físicos e os valores medidos estão no [relatório de depuração](../chatgpt-python-unfold/docs/DEPURACAO.md). A [configuração do GitHub Actions](../.github/workflows/validar-unfolding.yml) mostra como obter o mio e compilar os dois executáveis.

## Receitas do HOWTO no Linux

Para conferir também os exemplos Fortran distribuídos com o projeto, instale Matplotlib e rode os dois scripts de exemplo pelo comando abaixo. Passe a pasta `app` da compilação, a raiz do mio e o destino dos resultados.

```bash
python -m pip install matplotlib
bash Testes/check_howto.sh /caminho/build/native/app /caminho/mio results/howto
```

Em uma subpasta nova, o teste roda o diamante e o grafeno 3×3 e gera um PNG para cada mapa refoldado. O modo `native` acima faz a comparação quantitativa no grafeno 2×2 e no diamante.

Os testes Python reunidos nesta pasta acompanham a implementação identificada como ChatGPT e sua [licença MIT](../chatgpt-python-unfold/LICENSE). Os arquivos de referência do grafeno indicam a procedência das entradas e dos parâmetros no [README dos dados](../chatgpt-python-unfold/examples/graphene/reference/README.md).
