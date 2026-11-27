# Referências

A numeração é a mesma de `METODO.md`.

**Método de unfolding**

1. S. G. Mayo, F. Yndurain, J. M. Soler, "Band unfolding made simple", J. Phys.: Condens. Matter **32**, 205902 (2020), doi:10.1088/1361-648X/ab6e8e. Equações (A-1) a (A-14): densidade espectral no espaço dos momentos, refolding e projeção LCAO.
2. S. G. Mayo, F. Yndurain, J. M. Soler, arXiv:1812.03925 (preprint de [1]).
3. J. M. Soler et al., "The SIESTA method for ab initio order-N materials simulation", J. Phys.: Condens. Matter **14**, 2745 (2002).
4. SIESTA, `Util/Unfolding` (README e `Src/unfold.F90`), S. G. Mayo e J. M. Soler (2018), com modificações de P. Ordejón (2021): implementação de referência.
5. T. Kosugi, H. Nishi, Y. Kato, Y. Matsushita, J. Phys. Soc. Jpn. **86**, 124717 (2017): espectro desdobrado em ondas planas, citado em [1].

Outros trabalhos sobre unfolding citados em [1], para contexto: T. B. Boykin, G. Klimeck, Phys. Rev. B **71**, 115215 (2005); W. Ku, T. Berlijn, C.-C. Lee, Phys. Rev. Lett. **104**, 216401 (2010); V. Popescu, A. Zunger, Phys. Rev. Lett. **104**, 236403 (2010) e Phys. Rev. B **85**, 085201 (2012); P. B. Allen, T. Berlijn, D. A. Casavant, J. M. Soler, Phys. Rev. B **87**, 085322 (2013).

**DFTB+ e parâmetros**

6. B. Hourahine et al., J. Chem. Phys. **152**, 124101 (2020).
7. B. Aradi, B. Hourahine, Th. Frauenheim, J. Phys. Chem. A **111**, 5678 (2007).
8. M. Elstner et al., Phys. Rev. B **58**, 7260 (1998): parâmetros mio-1-1, usados na validação.
9. Código do DFTB+ analisado: commit `fedabbcb993ef26967141f275bb7d31a984e4bdc` (`waveplot`, `sparse2dense`, `mainio`, `formatout`).

**Origem do código**

10. Código Python: gerado com ChatGPT (OpenAI), incluído sem modificações em `chatgpt-python-unfold/`.
11. Código Fortran: gerado com Claude (Anthropic) e validado neste projeto, em `claude-fortran-unfold/`.
