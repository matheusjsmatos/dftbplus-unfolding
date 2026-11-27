!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Band unfolding of supercell calculations with a linear combination of atomic orbitals basis.
!!
!! Implementation of the method of S. G. Mayo, F. Yndurain and J. M. Soler, J. Phys.: Condens.
!! Matter 32, 205902 (2020) for DFTB+ eigenvectors. The Bloch states of the supercell calculation
!! are Fourier transformed (using the Slater type orbitals of the wfc.*.hsd files) to get the
!! spectral weight n(q, E) at arbitrary q ("fully unfolded bands"), which is then optionally summed
!! over the reciprocal vectors g of a refolding cell ("refolded bands").
program unfold
  use unfold_initunfold, only : actionCheck, actionListK, actionUnfold, TProgramVariables,&
      & TProgramVariables_init
  use unfold_radial, only : getRadialFT, realTessY, TRadialTable, TRadialTable_init
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : pi
  use dftbp_common_globalenv, only : destructGlobalEnv, initGlobalEnv
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : determinant33
  use, intrinsic :: iso_fortran_env, only : int64, output_unit
  implicit none (type, external)

  type(TProgramVariables) :: prog

  !! Orbital bookkeeping: atom, species-shell, l and m of every orbital of the supercell
  integer, allocatable :: orbAtom(:), orbShell(:), orbL(:), orbM(:)
  integer :: nOrbTot

  !! One radial table per species
  type(TRadialTable), allocatable :: tables(:)

  integer :: nq, ng, nK, nSpin, nState, nOrb
  integer, allocatable :: kIdx(:,:)
  real(dp), allocatable :: kNeeded(:,:)
  integer :: nNeeded, nMissing
  real(dp) :: qMaxNeeded, volume

  call initGlobalEnv()
  call TProgramVariables_init(prog, output_unit)

  nq = size(prog%qPoint, dim=2)
  ng = size(prog%gVec, dim=2)
  nK = prog%input%nKPoint
  nSpin = prog%input%nSpin
  nState = prog%input%nState
  nOrb = prog%input%nOrb
  volume = abs(determinant33(prog%input%geo%latVecs))

  write(output_unit, "(A,I0,A,I0,A,I0)") "Number of q-points: ", nq, ", refolding vectors g: ",&
      & ng, ", k-points in file: ", nK
  write(output_unit, "(A,F12.4,A)") "Supercell volume: ", volume, " Bohr^3"

  call findKPoints(prog, kIdx, kNeeded, nNeeded, nMissing, qMaxNeeded)

  select case (prog%action)
  case (actionListK)
    call writeKPointList(prog, kNeeded(:, 1:nNeeded))
  case (actionCheck)
    call checkInputFiles(prog)
    call setupOrbitals()
    call checkBasis()
  case (actionUnfold)
    if (nMissing > 0) then
      write(output_unit, "(A,I0,A,I0,A)") "ERROR: ", nMissing, " of ", nq * ng,&
          & " needed (q+g) combinations have no matching k-point in the eigenvector file."
      write(output_unit, "(A,I0,A)") "Run unfold with 'Action = ListKPoints' to get the ",&
          & nNeeded, " k-points needed and use them in a second DFTB+ calculation."
      call error("Missing k-points")
    end if
    call checkInputFiles(prog)
    call setupOrbitals()
    call computeBands()
  end select

  call destructGlobalEnv()


contains


  !> Finds for every (q, g) the k-point of the DFTB+ calculation, which contains the Bloch state
  !! with the wave vector K = q + g - G (G being a reciprocal supercell vector).
  !!
  !! DFTB+ builds H(k) with the phase exp(+i k.R) for the matrix element between an orbital at the
  !! origin and one in the cell R, which corresponds to the physical Bloch vector K = -k, if the
  !! Bloch sum is defined as sum_R exp(i K.R) phi(r - R - r_mu).
  subroutine findKPoints(prog, kIdx, kNeeded, nNeeded, nMissing, qMaxNeeded)

    !> Program variables
    type(TProgramVariables), intent(in) :: prog

    !> Index of the file k-point for every (q, g), zero if not found
    integer, allocatable, intent(out) :: kIdx(:,:)

    !> List of the distinct k-points needed (as in the DFTB+ convention, reduced to [-0.5, 0.5))
    real(dp), allocatable, intent(out) :: kNeeded(:,:)

    !> Number of entries in kNeeded
    integer, intent(out) :: nNeeded

    !> Number of (q, g) pairs without k-point in the file
    integer, intent(out) :: nMissing

    !> Largest |q + g|
    real(dp), intent(out) :: qMaxNeeded

    integer :: iq, ig, ik, ii, jj
    real(dp) :: qg(3), frac(3), kNeed(3), dd(3)
    logical :: tFound, tHave

    allocate(kIdx(nq, ng), source=0)
    allocate(kNeeded(3, nq * ng))
    nNeeded = 0
    nMissing = 0
    qMaxNeeded = 0.0_dp
    do ig = 1, ng
      do iq = 1, nq
        qg(:) = prog%qPoint(:, iq) + prog%gVec(:, ig)
        qMaxNeeded = max(qMaxNeeded, norm2(qg))
        do ii = 1, 3
          frac(ii) = dot_product(prog%input%geo%latVecs(:, ii), qg) / (2.0_dp * pi)
        end do
        kNeed(:) = real(prog%kSign, dp) * frac
        tFound = .false.
        if (prog%action == actionUnfold) then
          do ik = 1, nK
            dd(:) = kNeed - prog%input%kPoint(:, ik)
            if (all(abs(dd - nint(dd)) < prog%kTol)) then
              kIdx(iq, ig) = ik
              tFound = .true.
              exit
            end if
          end do
          if (.not. tFound) nMissing = nMissing + 1
        end if
        if (tFound) cycle
        ! register in the list of needed k-points
        kNeed(:) = kNeed - nint(kNeed)
        tHave = .false.
        do jj = 1, nNeeded
          dd(:) = kNeed - kNeeded(:, jj)
          if (all(abs(dd - nint(dd)) < prog%kTol)) then
            tHave = .true.
            exit
          end if
        end do
        if (.not. tHave) then
          nNeeded = nNeeded + 1
          kNeeded(:, nNeeded) = kNeed
        end if
      end do
    end do

  end subroutine findKPoints


  !> Writes the needed k-points as input for DFTB+.
  subroutine writeKPointList(prog, kList)

    !> Program variables
    type(TProgramVariables), intent(in) :: prog

    !> k-points (relative to the reciprocal supercell vectors)
    real(dp), intent(in) :: kList(:,:)

    integer :: fd, ii
    character(len=:), allocatable :: fileName

    fileName = trim(prog%outPrefix) // ".kpoints.hsd"
    open(newunit=fd, file=fileName, action="write", status="replace")
    write(fd, "(A)") "KPointsAndWeights = {"
    do ii = 1, size(kList, dim=2)
      write(fd, "(3(F20.14,1X),ES14.6)") kList(:, ii), 1.0_dp / real(size(kList, dim=2), dp)
    end do
    write(fd, "(A)") "}"
    close(fd)
    write(output_unit, "(A,I0,A)") "Wrote ", size(kList, dim=2), " k-points to '" // fileName // "'"
    if (size(kList, dim=2) == 0) write(output_unit, "(A)") "(Nothing to do: all needed k-points&
        & are already present.)"

  end subroutine writeKPointList


  !> Consistency checks of the DFTB+ output against the input.
  subroutine checkInputFiles(prog)

    !> Program variables
    type(TProgramVariables), intent(in) :: prog

    integer :: fd, runId, ii
    integer(int64) :: fileSize, expected, elemBytes

    if (prog%input%nState /= prog%input%nOrb) then
      call error("Unsupported: number of states differs from the number of orbitals&
          & (non-collinear spin?)")
    end if
    open(newunit=fd, file=trim(prog%eigVecBin), access="stream", form="unformatted",&
        & action="read", status="old", iostat=ii)
    if (ii /= 0) call error("Cannot open eigenvector file '" // trim(prog%eigVecBin) // "'")
    inquire(unit=fd, size=fileSize)
    read(fd) runId
    close(fd)
    if (runId /= prog%input%identity) then
      call warning(output_unit, "Run id in '" // trim(prog%eigVecBin) // "' differs from the one in&
          & detailed.xml")
    end if
    elemBytes = 16_int64
    if (prog%input%tRealHam) elemBytes = 8_int64
    expected = 4_int64 + int(nK, int64) * nSpin * nState * nOrb * elemBytes
    if (fileSize /= expected) then
      write(output_unit, "(A,I0,A,I0)") "Size of eigenvector file: ", fileSize, ", expected: ",&
          & expected
      call error("The eigenvector file has not the expected size. Binary files must be written&
          & with stream access (default) and WriteEigenvectors = Yes.")
    end if

  end subroutine checkInputFiles


  !> Creates the orbital tables of the supercell.
  subroutine setupOrbitals()

    integer :: iAt, iSp, iSh, ll, mm, ind
    integer :: nTot
    real(dp) :: qTab

    nTot = 0
    do iAt = 1, prog%input%geo%nAtom
      iSp = prog%input%geo%species(iAt)
      do iSh = 1, size(prog%basis(iSp)%shells)
        nTot = nTot + 2 * prog%basis(iSp)%shells(iSh)%ll + 1
      end do
    end do
    if (nTot /= nOrb) then
      write(output_unit, "(A,I0,A,I0)") "Orbitals from basis file: ", nTot, ", in DFTB+ run: ",&
          & nOrb
      call error("Basis file and DFTB+ calculation have different numbers of orbitals.")
    end if
    nOrbTot = nTot
    allocate(orbAtom(nTot), orbShell(nTot), orbL(nTot), orbM(nTot))
    ind = 0
    do iAt = 1, prog%input%geo%nAtom
      iSp = prog%input%geo%species(iAt)
      do iSh = 1, size(prog%basis(iSp)%shells)
        ll = prog%basis(iSp)%shells(iSh)%ll
        do mm = -ll, ll
          ind = ind + 1
          orbAtom(ind) = iAt
          orbShell(ind) = iSh
          orbL(ind) = ll
          orbM(ind) = mm
        end do
      end do
    end do

    ! Fourier transformed radial functions
    qTab = max(qMaxNeeded * 1.01_dp + 0.1_dp, prog%sumRuleQMax + 0.1_dp)
    write(output_unit, "(A,F8.2,A)") "Tabulating radial Fourier transforms up to q = ", qTab,&
        & " 1/Bohr"
    allocate(tables(size(prog%basis)))
    do iSp = 1, size(prog%basis)
      call TRadialTable_init(tables(iSp), prog%basis(iSp)%shells, qTab, 0.01_dp)
      do iSh = 1, size(prog%basis(iSp)%shells)
        write(output_unit, "(A,A,A,I0,A,I0,A,F9.5)") "  species ",&
            & trim(prog%input%geo%speciesNames(iSp)), ", shell ", iSh, ", l = ",&
            & prog%basis(iSp)%shells(iSh)%ll, ", norm of STO = ", prog%basis(iSp)%shells(iSh)%norm
      end do
    end do

  end subroutine setupOrbitals


  !> Amplitudes A_mu(q) = (-i)^l Y_lm(q^) F_l(|q|) exp(-i q.r_mu) of all orbitals.
  subroutine orbitalAmplitudes(qq, amp)

    !> Wave vector (Cartesian, 1/Bohr)
    real(dp), intent(in) :: qq(3)

    !> Amplitudes of all orbitals
    complex(dp), intent(out) :: amp(:)

    complex(dp), parameter :: iUnit = (0.0_dp, 1.0_dp)
    complex(dp) :: pref(0:3)
    real(dp) :: qMod, radial, phi
    integer :: iMu, iSp, iAt

    pref(0) = (1.0_dp, 0.0_dp)
    pref(1) = -iUnit
    pref(2) = (-1.0_dp, 0.0_dp)
    pref(3) = iUnit

    qMod = norm2(qq)
    do iMu = 1, nOrbTot
      iAt = orbAtom(iMu)
      iSp = prog%input%geo%species(iAt)
      radial = getRadialFT(tables(iSp), orbShell(iMu), qMod)
      phi = -dot_product(qq, prog%input%geo%coords(:, iAt))
      amp(iMu) = pref(orbL(iMu)) * realTessY(orbL(iMu), orbM(iMu), qq) * radial&
          & * cmplx(cos(phi), sin(phi), dp)
    end do

  end subroutine orbitalAmplitudes


  !> Reads one block of eigenvectors (all states of one k-point and spin).
  subroutine readEigenvectors(fd, iKS, block)

    !> Unit of the stream file
    integer, intent(in) :: fd

    !> Index of the (k-point, spin) combination (k-point running fastest)
    integer, intent(in) :: iKS

    !> Eigenvectors (orbital, state)
    complex(dp), intent(out) :: block(:,:)

    real(dp), allocatable :: realBlock(:,:)
    integer(int64) :: pos

    if (prog%input%tRealHam) then
      allocate(realBlock(nOrb, nState))
      pos = 5_int64 + int(iKS - 1, int64) * nState * nOrb * 8_int64
      read(fd, pos=pos) realBlock
      block(:,:) = cmplx(realBlock, 0.0_dp, dp)
    else
      pos = 5_int64 + int(iKS - 1, int64) * nState * nOrb * 16_int64
      read(fd, pos=pos) block
    end if

  end subroutine readEigenvectors


  !> Reads the eigenvalues from band.out (in eV).
  subroutine readBandOut(fileName, eig)

    !> File name
    character(len=*), intent(in) :: fileName

    !> Eigenvalues (state, k-point, spin) in eV
    real(dp), intent(out) :: eig(:,:,:)

    integer :: fd, iErr, iK, iS, ii, idx
    character(len=256) :: line, dummy1, dummy2, dummy3
    real(dp) :: weight, occ

    open(newunit=fd, file=fileName, action="read", status="old", iostat=iErr)
    if (iErr /= 0) call error("Cannot open '" // fileName // "'. Is WriteBandOut = Yes?")
    do iS = 1, size(eig, dim=3)
      do iK = 1, size(eig, dim=2)
        read(fd, "(A)", iostat=iErr) line
        if (iErr /= 0) call error("Unexpected end of '" // fileName // "'")
        read(line, *, iostat=iErr) dummy1, ii, dummy2, idx, dummy3, weight
        if (iErr /= 0 .or. ii /= iK .or. idx /= iS) then
          call error("Unexpected header in '" // fileName // "': " // trim(line))
        end if
        do ii = 1, size(eig, dim=1)
          read(fd, *, iostat=iErr) idx, eig(ii, iK, iS), occ
          if (iErr /= 0) call error("Cannot read eigenvalues from '" // fileName // "'")
        end do
        read(fd, "(A)", iostat=iErr) line
      end do
    end do
    close(fd)

  end subroutine readBandOut


  !> Adds the weight of one state to an energy histogram, with optional Gaussian broadening.
  subroutine depositWeight(hist, ee, ww, eMin, dE, sigma)

    !> Histogram (0:nE)
    real(dp), intent(inout) :: hist(0:)

    !> Energy of the state, relative to the Fermi level (eV)
    real(dp), intent(in) :: ee

    !> Weight
    real(dp), intent(in) :: ww

    !> Lowest energy of the mesh and mesh spacing
    real(dp), intent(in) :: eMin, dE

    !> Gaussian width (eV), zero for linear binning
    real(dp), intent(in) :: sigma

    integer :: je, j0, j1, jj, nE
    real(dp) :: we, xx

    nE = ubound(hist, dim=1)
    if (sigma <= 0.0_dp) then
      je = floor((ee - eMin) / dE)
      we = (eMin + real(je + 1, dp) * dE - ee) / dE
      if (je >= 0 .and. je <= nE) hist(je) = hist(je) + ww * we
      if (je >= -1 .and. je < nE) hist(je + 1) = hist(je + 1) + ww * (1.0_dp - we)
    else
      j0 = max(0, ceiling((ee - 5.0_dp * sigma - eMin) / dE))
      j1 = min(nE, floor((ee + 5.0_dp * sigma - eMin) / dE))
      do jj = j0, j1
        xx = (eMin + real(jj, dp) * dE - ee) / sigma
        hist(jj) = hist(jj) + ww * dE * exp(-0.5_dp * xx**2) / (sigma * sqrt(2.0_dp * pi))
      end do
    end if

  end subroutine depositWeight


  !> Main calculation of the fully unfolded and refolded spectral weights.
  subroutine computeBands()

    complex(dp), allocatable :: block(:,:), amps(:), proj(:)
    real(dp), allocatable :: eig(:,:,:), unf(:,:,:), refold(:,:,:)
    real(dp) :: qg(3), dE, pref, ee
    integer :: fd, iS, iK, iq, ig, iState, nE, iErr, iKS
    logical :: tDone

    nE = prog%nEnergy
    dE = (prog%eMax - prog%eMin) / real(nE, dp)
    allocate(eig(nState, nK, nSpin))
    call readBandOut(trim(prog%bandOut), eig)
    eig(:,:,:) = eig - prog%fermiLevel

    allocate(unf(0:nE, nq, nSpin), source=0.0_dp)
    allocate(refold(0:nE, nq, nSpin), source=0.0_dp)
    allocate(block(nOrb, nState), amps(nOrb), proj(nState))
    pref = (2.0_dp * pi)**3 / volume

    open(newunit=fd, file=trim(prog%eigVecBin), access="stream", form="unformatted",&
        & action="read", status="old")

    tDone = .false.
    do iS = 1, nSpin
      do iK = 1, nK
        if (.not. any(kIdx == iK)) cycle
        iKS = (iS - 1) * nK + iK
        call readEigenvectors(fd, iKS, block)
        if (.not. tDone .and. prog%sumRuleCheck) then
          call sumRule(iK, block)
          tDone = .true.
        end if
        do ig = 1, ng
          do iq = 1, nq
            if (kIdx(iq, ig) /= iK) cycle
            qg(:) = prog%qPoint(:, iq) + prog%gVec(:, ig)
            call orbitalAmplitudes(qg, amps)
            proj(:) = matmul(amps, block)
            do iState = 1, nState
              ee = eig(iState, iK, iS)
              if (ee < prog%eMin - 5.0_dp * prog%broadening - dE&
                  & .or. ee > prog%eMax + 5.0_dp * prog%broadening + dE) cycle
              if (ig == 1) call depositWeight(unf(:, iq, iS), ee, pref * abs(proj(iState))**2,&
                  & prog%eMin, dE, prog%broadening)
              call depositWeight(refold(:, iq, iS), ee, pref * abs(proj(iState))**2, prog%eMin,&
                  & dE, prog%broadening)
            end do
          end do
        end do
      end do
    end do
    close(fd)

    call writeResults(unf, refold)

  end subroutine computeBands


  !> Runs the sum rule check for the first few k-points of the eigenvector file (all spins).
  subroutine checkBasis()

    complex(dp), allocatable :: block(:,:)
    integer :: fd, iK

    allocate(block(nOrb, nState))
    open(newunit=fd, file=trim(prog%eigVecBin), access="stream", form="unformatted",&
        & action="read", status="old")
    do iK = 1, min(nK, 12)
      call readEigenvectors(fd, iK, block)
      call sumRule(iK, block)
    end do
    close(fd)

  end subroutine checkBasis


  !> Checks the sum rule sum_G |u_G|^2 = 1 for the states of one k-point. It tests the transforms
  !! of the radial functions, the orbital ordering, and the consistency of the overlap matrix and
  !! the STO functions.
  subroutine sumRule(iK, block)

    !> k-point
    integer, intent(in) :: iK

    !> Eigenvectors
    complex(dp), intent(in) :: block(:,:)

    complex(dp), allocatable :: amps(:), proj(:)
    real(dp), allocatable :: total(:)
    real(dp) :: qq(3), qLim, pref, kFrac(3)
    integer :: nMax(3), i1, i2, i3, ii, nVec
    real(dp) :: nEst

    qLim = prog%sumRuleQMax
    do ii = 1, 3
      nMax(ii) = ceiling(qLim * norm2(prog%input%geo%latVecs(:, ii)) / (2.0_dp * pi)) + 1
    end do
    nEst = 4.0_dp / 3.0_dp * pi * qLim**3 * volume / (2.0_dp * pi)**3
    if (nEst * nOrb * nState > 5.0e10_dp) then
      write(output_unit, "(A)") "Sum rule check skipped (system too large). Use SumRuleCheck = No&
          & to suppress this message."
      return
    end if
    allocate(amps(nOrb), proj(nState))
    allocate(total(nState), source=0.0_dp)
    pref = (2.0_dp * pi)**3 / volume
    kFrac(:) = real(prog%kSign, dp) * prog%input%kPoint(:, iK)
    nVec = 0
    do i3 = -nMax(3), nMax(3)
      do i2 = -nMax(2), nMax(2)
        do i1 = -nMax(1), nMax(1)
          qq(:) = matmul(prog%supRecVec, kFrac + real([i1, i2, i3], dp))
          if (norm2(qq) > qLim) cycle
          nVec = nVec + 1
          call orbitalAmplitudes(qq, amps)
          proj(:) = matmul(amps, block)
          total(:) = total + pref * abs(proj)**2
        end do
      end do
    end do
    write(output_unit, "(A,I0,A,F6.2,A)") "Sum rule check at file k-point ", iK, " (|q| <= ",&
        & qLim, " 1/Bohr):"
    write(output_unit, "(A,I0,A,3F10.6)") "  ", nVec, " plane waves; sum over G of the weights per&
        & state (min, mean, max): ", minval(total), sum(total) / real(nState, dp), maxval(total)
    write(output_unit, "(A)") "  (should be close to 1; deviations indicate an inconsistent&
        & basis/overlap or a too small q range)"

  end subroutine sumRule


  !> Writes the unfolded and refolded bands (format of the SIESTA unfold utility).
  subroutine writeResults(unf, refold)

    !> Fully unfolded weight (0:nE, nq, nSpin)
    real(dp), intent(in) :: unf(0:,:,:)

    !> Refolded weight
    real(dp), intent(in) :: refold(0:,:,:)

    integer :: iPath, iS
    character(len=1024) :: fileName

    do iPath = 1, prog%nPath
      do iS = 1, nSpin
        call outName("unfoldedBands", iPath, iS, fileName)
        call writeOne(fileName, iPath, unf(:, :, iS))
        if (prog%tRefold) then
          call outName("refoldedBands", iPath, iS, fileName)
          call writeOne(fileName, iPath, refold(:, :, iS))
        end if
      end do
    end do

  end subroutine writeResults


  !> Constructs the name of an output file.
  subroutine outName(suffix, iPath, iS, name)

    !> Kind of the file
    character(len=*), intent(in) :: suffix

    !> Path and spin index
    integer, intent(in) :: iPath, iS

    !> Resulting name
    character(len=*), intent(out) :: name

    character(len=32) :: tmp

    name = trim(prog%outPrefix) // "." // suffix
    if (nSpin == 2) then
      if (iS == 1) then
        name = trim(name) // ".spinUp"
      else
        name = trim(name) // ".spinDn"
      end if
    end if
    if (prog%nPath > 1) then
      write(tmp, "(I0)") iPath
      name = trim(name) // ".path" // trim(tmp)
    end if

  end subroutine outName


  !> Writes the spectral weights of one path.
  subroutine writeOne(name, iPath, data)

    !> File name
    character(len=*), intent(in) :: name

    !> Path index
    integer, intent(in) :: iPath

    !> Weights (0:nE, nq)
    real(dp), intent(in) :: data(0:,:)

    integer :: fd, iq, jj, nE

    nE = prog%nEnergy
    open(newunit=fd, file=trim(name), action="write", status="replace")
    write(fd, "(2(I6),2(F10.3),F13.6)") prog%lastQ(iPath) - prog%lastQ(iPath - 1), nE + 1,&
        & prog%eMin, prog%eMax, prog%fermiLevel
    do iq = prog%lastQ(iPath - 1) + 1, prog%lastQ(iPath)
      write(fd, "(3(F13.8),I5)") prog%qPoint(:, iq), prog%iLine(iq)
      do jj = 0, nE
        write(fd, "(E15.6)") data(jj, iq)
      end do
    end do
    close(fd)
    write(output_unit, "(A)") "Written '" // trim(name) // "'"

  end subroutine writeOne

end program unfold
