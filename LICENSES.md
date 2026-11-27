# Licenças

Este repositório combina duas partes com licenças diferentes:

| Pasta | Licença | Arquivo |
|---|---|---|
| `claude-fortran-unfold/` | LGPL-3.0 | `claude-fortran-unfold/LICENSE` |
| `chatgpt-python-unfold/` | MIT | `chatgpt-python-unfold/LICENSE` |

O código Fortran é feito para compilar dentro do DFTB+ e reutiliza partes dele (cabeçalhos, parser HSD e a convenção
das harmônicas reais do Waveplot). Por isso segue a licença do DFTB+, a LGPL-3.0. O código Python é independente do código-fonte do DFTB+.
A documentação (`docs/`) e os scripts de `validation/` acompanham a licença da pasta de código a que se referem.
