# Referências usadas como parâmetros

1. S. G. Mayo, F. Yndurain, J. M. Soler, "Band unfolding made simple", J. Phys.: Condens. Matter **32**, 205902 (2020),
   doi:10.1088/1361-648X/ab6e8e. Equações (6), (9) e (10)–(14): densidade espectral em momento, refolding e projeção LCAO.
2. SIESTA, `Util/Unfolding` (`unfold.F90` e README): rotina de referência. Mayo & Soler, arXiv:1812.03925.
3. DFTB+: B. Aradi, B. Hourahine, Th. Frauenheim, J. Phys. Chem. A **111**, 5678 (2007) e B. Hourahine et al.,
   J. Chem. Phys. **152**, 124101 (2020). Código analisado no commit `fedabbcb993ef26967141f275bb7d31a984e4bdc`
   (`waveplot`, `sparse2dense`, `mainio`, `formatout`).
4. M. Elstner et al., Phys. Rev. B **58**, 7260 (1998): parâmetros mio-1-1, usados na validação (SK e STO no mesmo conjunto).
5. Lista completa de referências do artigo (1)–(30+): ver o próprio artigo. As citadas diretamente aqui
   são as de unfolding em TB/LCAO, como Boykin & Klimeck 2005, Ku et al. 2010, Popescu & Zunger 2010/2012 e Allen et al. 2013.
