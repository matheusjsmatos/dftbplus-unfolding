# Unfolding de bandas para o DFTB+

Duas implementações independentes de *band unfolding* para o **DFTB+**, ambas baseadas no método de
Mayo, Yndurain e Soler, J. Phys.: Condens. Matter **32**, 205902 (2020), "Band unfolding made simple".
O utilitário `Util/Unfolding` do SIESTA serviu de modelo de referência.

| Pasta | Origem | Linguagem | Entrada | Observação |
|---|---|---|---|---|
| `claude-fortran-unfold/` | Claude (Anthropic) | Fortran 2008 (dentro da árvore do DFTB+) | `eigenvec.bin` + `band.out` + `detailed.xml` | Programa `unfold`, compilado junto com o DFTB+ |
| `chatgpt-python-unfold/` | **ChatGPT (OpenAI)** | Python / NumPy | `hamreal1.dat` + `overreal.dat` (WriteRealHS) | Diagonaliza em K arbitrário, saída CSV |

Para reproduzir o grafeno 2×2, gerar o gráfico ou rodar o diamante a partir do DFTB+ 26.1, siga o [HOWTO](docs/HOWTO.md). O [Manual](docs/MANUAL.md) descreve os arquivos, as unidades e os parâmetros de cada implementação.

## Rodar os testes

Na raiz do repositório, com Python e NumPy, execute

```bash
python -m pip install -r chatgpt-python-unfold/requirements.txt
python Testes/run.py
```

A suíte executa os testes unitários, recompõe as bandas de grafeno a partir das matrizes arquivadas e confere os pesos da projeção. Para repetir a comparação com o executável Fortran em grafeno e diamante, veja [Testes](Testes/README.md). Os exemplos e as entradas auxiliares estão em [validation](validation/README.md).

## Método e desenvolvimento

- [Método e convenções do DFTB+](docs/METODO.md).
- [Desenvolvimento da implementação Fortran](docs/DESENVOLVIMENTO_CLAUDE.md).
- [Desenvolvimento da implementação Python](docs/DESENVOLVIMENTO_CHATGPT.md).
- [Comparação numérica](docs/COMPARACAO.md) e [relatório da execução nativa](chatgpt-python-unfold/docs/DEPURACAO.md).
- [Referências e parametrização](docs/REFERENCIAS.md).

## Software de destino

Os dois programas foram testados com o **DFTB+ 26.1** e usam uma base LCAO
de orbitais tipo Slater lidos de arquivos `wfc.*.hsd` (formato do Waveplot). Não funcionam com outros códigos.

## Uso rápido (Fortran)

Compile o `unfold` junto com o DFTB+ 26.1. O [HOWTO](docs/HOWTO.md) mostra como rodar o exemplo e a [receita usada no GitHub Actions](.github/workflows/validar-unfolding.yml) traz o pacote-fonte, as dependências e os comandos de compilação em Linux. O [guia do aplicativo Fortran](claude-fortran-unfold/app_unfold/README.rst) detalha o HSD.

## Uso rápido (Python)

```bash
cd chatgpt-python-unfold
python -m pip install -r requirements.txt
python dftb_unfold.py examples/chain/unfold.json -o cadeia.csv
```

## Aviso importante

O conjunto de STOs do `wfc` **deve ser o mesmo** dos arquivos Slater-Koster usados no cálculo (por exemplo, mio-1-1
nos dois). Com 3ob nos SK e mio no wfc, a regra de soma falhou (0,4–3,8). Com mio-1-1 nos dois, ficou em 0,976–0,999.

## Licença

O repositório contém o [Fortran sob LGPL-3.0](claude-fortran-unfold/LICENSE) e o [Python sob MIT](chatgpt-python-unfold/LICENSE). [Licenças por diretório](LICENSES.md)
