# Como o unfolding é feito

## Ideia (Mayo, Yndurain, Soler 2020)

Um estado de Bloch da supercélula, com vetor K na zona de Brillouin da supercélula, é decomposto em ondas planas
nos momentos q = K + G (G da rede recíproca da supercélula). O peso espectral em q é o módulo ao quadrado
da amplitude de Fourier.

Para uma base LCAO, com orbitais φ_μ = R(r) Y_lm centrados em r_μ, e coeficientes c_μ(K):

    ũ(q) = (2π)^{3/2} / V_sc · Σ_μ c_μ(K) φ̃_μ(q) e^{-i q·r_μ}
    φ̃_μ(q) = (−i)^l Y_lm(q̂) F_l(|q|),   F_l(q) = √(2/π) ∫ r² R(r) j_l(qr) dr
    w(q) = (2π)³ / V_sc · |Σ_μ c_μ φ̃_μ e^{-i q·r_μ}|²

- **Unfolding completo:** n(q,E) = Σ_i w_i(q) δ(E − ε_i). Cada estado é distribuído no espaço recíproco inteiro.
- **Refolding** na zona de Brillouin da célula de referência: n_RBZ(k,E) = Σ_g n(k+g,E), com g vetores recíprocos
  da célula de referência e |g|² < corte (`RefoldingGCutoff` no Fortran, `gcut_ry` no Python).
- **Regra de soma:** Σ_G w = c†S_STO c. Vale ≈ 1 se a base STO usada na transformada é a do cálculo SK.
  Serve como teste de consistência da base.

## Passos computados

1. Ler geometria, k-pontos, autovalores e autovetores LCAO do DFTB+.
2. Tabelar F_l(q) para cada camada de STO (integração radial até o raio de corte do `wfc`).
3. Para cada q e cada estado, somar sobre orbitais com a fase e^{-iq·r_μ} e a harmônica real do orbital.
4. Acumular pesos em malha de energia (binning linear, alargamento gaussiano opcional).
5. Para o refolding, somar sobre g.

## Convenções específicas do DFTB+

- H(k) usa a fase exp(+ik·R) só nas translações de célula, e não nas posições atômicas. Logo o vetor
  físico de Bloch é K = −k_arquivo (`KPointSign = −1`, verificado numericamente). As fases atômicas entram via e^{-iq·r_μ}.
- Harmônicas reais na ordem do DFTB+/Waveplot: s; p = y,z,x; d = xy,yz,z²,xz,x²−y²; até l = 3.
- Base radial lida de `wfc.*.hsd`: soma de STOs com expoentes e coeficientes polinomiais e raio de corte.
- `eigenvec.bin`: arquivo *stream* com runId (4 bytes) e registros de nOrb por estado (complexo, ou real em Γ real).
- Banda não-SCC: `MaxSCCIterations = 1` com `ReadInitialCharges = Yes`.
