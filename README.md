# Unfolding de bandas para o DFTB+

Duas implementações independentes de *band unfolding* para o **DFTB+**, ambas baseadas no método de
Mayo, Yndurain e Soler, J. Phys.: Condens. Matter **32**, 205902 (2020), "Band unfolding made simple".
O utilitário `Util/Unfolding` do SIESTA serviu de modelo de referência.

| Pasta | Origem | Linguagem | Entrada | Observação |
|---|---|---|---|---|
| `claude-fortran-unfold/` | Claude (Anthropic), com validação | Fortran 2008 (dentro da árvore do DFTB+) | `eigenvec.bin` + `band.out` + `detailed.xml` | Programa `unfold`, compilado junto com o DFTB+ |
| `chatgpt-python-unfold/` | **ChatGPT (OpenAI)**, sem modificações | Python / NumPy | `hamreal1.dat` + `overreal.dat` (WriteRealHS) | Diagonaliza em K arbitrário, saída CSV |

Documentação em `docs/`:

- `METODO.md`: como o unfolding é feito (equações e convenções do DFTB+).
- `DESENVOLVIMENTO_CLAUDE.md`: procedimento de desenvolvimento do código Fortran.
- `DESENVOLVIMENTO_CHATGPT.md`: o que o código do ChatGPT faz e como foi analisado.
- `COMPARACAO.md`: validação cruzada entre as duas.
- `REFERENCIAS.md`: referências usadas como parâmetros.

## Software de destino

Ambos são feitos para o **DFTB+** (commit analisado `fedabbcb993ef26967141f275bb7d31a984e4bdc`), com base LCAO
de orbitais tipo Slater lidos de arquivos `wfc.*.hsd` (formato do Waveplot). Não funcionam com outros códigos.

## Uso rápido (Fortran)

```bash
git clone https://github.com/dftbplus/dftbplus && cd dftbplus
git checkout fedabbcb993ef26967141f275bb7d31a984e4bdc
cp -r /caminho/deste/repo/claude-fortran-unfold/app_unfold app/unfold
git apply /caminho/deste/repo/claude-fortran-unfold/unfold-cmake.patch   # adiciona add_subdirectory(unfold)
# compile o DFTB+ como de costume; o executável "unfold" aparece em build/app/unfold
```

Fluxo e chaves de entrada: `claude-fortran-unfold/app_unfold/README.rst`. Exemplo: `app_unfold/examples/diamond/run_all.sh`.

## Uso rápido (Python)

```bash
cd chatgpt-python-unfold
python -m pip install -r requirements.txt
python dftb_unfold.py examples/chain/unfold.json -o cadeia.csv
python -m unittest -v test_unfold
```

## Aviso importante

O conjunto de STOs do `wfc` **deve ser o mesmo** dos arquivos Slater-Koster usados no cálculo (por exemplo, mio-1-1
nos dois). Com 3ob nos SK e mio no wfc, a regra de soma falhou (0,4–3,8). Com mio-1-1 nos dois, ficou em 0,976–0,999.

## Licença

Duas licenças, detalhadas em `LICENSES.md`: código Fortran em LGPL-3.0 (compatível com o DFTB+) e código Python do ChatGPT em MIT.
