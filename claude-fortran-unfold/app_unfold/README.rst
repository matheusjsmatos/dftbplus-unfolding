======
unfold
======

Band unfolding of DFTB+ supercell calculations, following S. G. Mayo, F. Yndurain and J. M. Soler,
*Band unfolding made simple*, J. Phys.: Condens. Matter **32**, 205902 (2020) (the method of the
``Util/Unfolding`` tool of SIESTA).

The Bloch states of the supercell calculation are Fourier transformed. With the atomic-orbital
expansion of DFTB+ (valence orbitals on the atoms) the amplitude at the wave vector
``q = K + G`` (``K`` in the supercell Brillouin zone, ``G`` a reciprocal supercell vector) is::

  u(q) = (2 pi)^(3/2) / V_sc * sum_mu c_mu(K) * phi~_mu(q) * exp(-i q.r_mu)

and the *fully unfolded bands* are ``n(q, E) = sum_i |u_i(q)|^2 delta(E - e_i(K))``. Summing over
the reciprocal vectors ``g`` of a refolding cell (e.g. the primitive cell of the pristine crystal)
gives the *refolded bands*, ``n_RBZ(k, E) = sum_g n(k + g, E)``. The tool needs no symmetry and no
commensurability between supercell and refolding cell, it can be used for defects, alloys,
amorphous systems, etc.

``unfold`` is a post-processing tool (like ``waveplot``): it reads files written by DFTB+ and does
**not** modify the DFTB+ program.


Workflow
========

1. Converged SCC calculation of the supercell (any k-point set). Keep ``charges.bin``.

2. ``Action = ListKPoints``: ``unfold`` reads ``detailed.xml`` and the q-path and writes
   ``<Output>.kpoints.hsd`` with the k-points the DFTB+ calculation must contain
   (``KPointsAndWeights = { ... }``, in units of the reciprocal supercell vectors).

3. Second, non-SCC DFTB+ run at these k-points using the charges of step 1::

     Hamiltonian = DFTB {
       SCC = Yes
       MaxSCCIterations = 1
       ReadInitialCharges = Yes
       <<+ "unfold.kpoints.hsd"
       ...
     }
     Analysis { WriteEigenvectors = Yes }
     Options { WriteDetailedXML = Yes }

   (``WriteBandOut`` is on by default.) For a commensurate supercell the number of different
   k-points is small, because many ``q + g`` fold to the same ``K``.

4. ``Action = Unfold`` (default): reads ``detailed.xml``, ``eigenvec.bin``, ``band.out`` and the
   STO basis and writes ``<Output>.unfoldedBands`` and ``<Output>.refoldedBands`` (the latter
   only if ``RefoldingLatticeVectors`` is given). The format is the one of the SIESTA tool, so
   its plot scripts can be used.

A complete example (cubic 8-atom diamond cell refolded onto the fcc primitive cell, including the
comparison with the primitive-cell band structure) is in ``examples/diamond``.


Input (``unfold_in.hsd``)
=========================

``Action``
  ``Unfold`` (default), ``ListKPoints`` or ``CheckBasis`` (only runs the sum rule test below).

``DetailedXML``, ``EigenvecBin``, ``BandOut``, ``Output``
  Names of the DFTB+ files (defaults ``detailed.xml``, ``eigenvec.bin``, ``band.out``) and the
  prefix of the output files (default ``unfold``). The eigenvector file must be written with
  stream access (the DFTB+ default).

``Basis``
  The Slater-type orbitals, in the format used by ``waveplot`` (``wfc.<set>.hsd`` of the
  Slater-Koster parameter set). The file can be included: ``Basis { <<+ "wfc_C.hsd" }``.
  **The STOs must belong to the same parameter set as the SK-files.** Otherwise the overlap of the
  STOs is inconsistent with the one DFTB+ used to normalise the eigenvectors (see the sum rule).
  Angular momenta up to f are supported. ``NormalizeOrbitals`` (default ``Yes``) normalises the
  radial functions.

``QPoints``
  Lines as in ``KPointsAndWeights = Klines``: ``1 q1 q2 q3`` for the first point and
  ``n q1 q2 q3`` for ``n`` equidistant points ending at the given one. Coordinates are relative to
  the reciprocal vectors of the refolding cell (the supercell, if none is given). Several paths
  can be given as ``Path { ... }`` blocks; each is written to a separate file. Coordinates
  outside the first zone give the fully unfolded bands in the extended zone.

``RefoldingLatticeVectors``
  Lattice vectors of the refolding cell (rows, with unit modifier like ``[Angstrom]``).
  Without it only the fully unfolded bands are calculated.

``RefoldingGCutoff``
  Cutoff for the refolding vectors, ``|g|^2/2 < cutoff`` (default 10 Hartree). This determines
  accuracy and cost of the refolding (the weight of ``q + g`` decays quickly with ``|g|``).

``Energy``
  ``Window`` (default -20..10 eV, relative to the Fermi level), ``NrOfIntervals`` (300),
  ``FermiLevel`` (energy zero; default 0, so give the Fermi level of the SCC run), and
  ``Broadening`` (Gaussian width; default 0 = linear binning of the weights like SIESTA).

``SumRuleCheck``, ``SumRuleQMax``
  See below.


Output format
=============

First line: ``nq ne+1 Emin Emax EFermi``. Then for every q: a line with the Cartesian ``q`` (1/Bohr)
and the path-segment index, followed by the ``ne+1`` weights on the energy mesh. Energies are
relative to the Fermi level. The weights are normalised so that every state contributes in
total 1 (summed over all ``G``); a bin therefore contains the number of states (weight) in it,
not a density.


Conventions and checks
======================

* DFTB+ assembles ``H(k)`` with the phase ``exp(+i k.R)`` of the *cell* translation only. The
  Bloch vector of the state in the sense of ``sum_R exp(iK.R) phi(r-R-r_mu)`` is therefore
  ``K = -k`` (the default, ``KPointSign = -1``), and the phase of the atoms enters via
  ``exp(-i q.r_mu)`` with the positions of ``detailed.xml``. Time reversal makes the final
  spectra insensitive to the sign, but not the sum rule at generic k-points, which was used to
  verify it.

* **Sum rule.** For every state ``sum_G (2 pi)^3/V |sum_mu c_mu phi~_mu(K+G) exp(-i(K+G).r_mu)|^2``
  equals ``c^+ S_STO(K) c``, which is 1 if the STO basis is consistent with the DFTB+ overlap.
  ``unfold`` prints min/mean/max of this sum for the first k-point(s) (``SumRuleCheck``,
  ``SumRuleQMax``, default 8/Bohr). Typical values for a matching set are 0.97-1.00; values like
  0.4-3.8 indicate a wrong ``wfc`` set (this happened when using the 3ob SK-files with the
  STOs of mio-1-1), a wrong orbital order or a wrong ``KPointSign``.

* The STOs are approximations of the confined atomic orbitals and the SK integrals are two-center
  approximations, so the unfolded spectral weights are consistent only up to the quality of that
  fit (the sum rule quantifies it). Non-collinear spin is not supported. Spin-polarised
  collinear calculations write one file per spin channel.

* The radial transforms are computed by numerical quadrature (``unfold_radial``) on a q-grid of
  0.01 / Bohr with linear interpolation.


Validation
==========

* 8-atom cubic diamond (mio-1-1) unfolded on the fcc primitive cell: all refolded weight lies
  within 0.15 eV of the bands of a primitive-cell DFTB+ calculation; the weighted band centroids
  agree on average to 1 meV; sum rule 0.976-0.999.
* 3x3 graphene supercell (18 atoms, hexagonal, 2D) unfolded on the primitive cell: refolded
  weight on the primitive bands; mean centroid deviation below 0.1 meV; sum rule 0.95-1.00.
