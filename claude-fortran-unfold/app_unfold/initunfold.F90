!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Input handling for the unfold tool (band unfolding of supercell calculations).
module unfold_initunfold
  use unfold_radial, only : TStoShell, TStoShell_init
  use dftbp_common_accuracy, only : dp
  use dftbp_common_constants, only : eV__Hartree, Hartree__eV, pi
  use dftbp_common_unitconversion, only : energyUnits, inverseLengthUnits, lengthUnits
  use dftbp_extlibs_xmlf90, only : char, destroyNode, fnode, fNodeList, getItem1, getLength, string
  use dftbp_io_charmanip, only : i2c, unquote
  use dftbp_io_formatout, only : printDftbHeader
  use dftbp_io_hsdparser, only : dumpHSD, parseHSD
  use dftbp_io_hsdutils, only : detailedError, detailedWarning, getChild, getChildren,&
      & getChildValue
  use dftbp_io_hsdutils2, only : convertUnitHsd, readHSDAsXML, warnUnprocessedNodes
  use dftbp_io_message, only : error, warning
  use dftbp_math_simplealgebra, only : determinant33, invert33
  use dftbp_type_linkedlist, only : asArray, asVector, destruct, init, len, TListIntR1,&
      & TListReal, TListRealR1
  use dftbp_type_typegeometry, only : TGeometry
  use dftbp_type_typegeometryhsd, only : readTGeometryHSD
  implicit none (type, external)

  private
  public :: TProgramVariables, TProgramVariables_init, TSpeciesShells
  public :: actionUnfold, actionListK, actionCheck


  !> Action: compute unfolded/refolded bands
  integer, parameter :: actionUnfold = 1

  !> Action: write the list of k-points, which are needed for the requested q-points
  integer, parameter :: actionListK = 2

  !> Action: only check the sum rule of the Fourier weights (consistency of basis and overlap)
  integer, parameter :: actionCheck = 3


  !> Radial shells of one species.
  type TSpeciesShells

    !> Atomic number as given in the basis file
    integer :: atomicNumber = 0

    !> Shells
    type(TStoShell), allocatable :: shells(:)

  end type TSpeciesShells


  !> Data from detailed.xml and the eigenvector file.
  type TInput

    !> Geometry of the (super)cell used in the DFTB+ calculation
    type(TGeometry) :: geo

    !> Identity of the DFTB+ run
    integer :: identity = 0

    !> Real Hamiltonian (Gamma point only)?
    logical :: tRealHam = .false.

    !> Number of k-points, spin channels, states, orbitals
    integer :: nKPoint = 0, nSpin = 0, nState = 0, nOrb = 0

    !> k-points in units of the reciprocal supercell vectors, as in the DFTB+ run (3, nKPoint)
    real(dp), allocatable :: kPoint(:,:)

  end type TInput


  !> Container for all program settings.
  type TProgramVariables

    !> What to do
    integer :: action = actionUnfold

    !> Input data from the DFTB+ run
    type(TInput) :: input

    !> File names
    character(len=1024) :: eigVecBin = "eigenvec.bin", bandOut = "band.out", outPrefix = "unfold"

    !> Radial basis for every species
    type(TSpeciesShells), allocatable :: basis(:)

    !> Energy window (eV, relative to the Fermi level), number of intervals and Fermi level (eV)
    real(dp) :: eMin = -20.0_dp, eMax = 10.0_dp, fermiLevel = 0.0_dp, broadening = 0.0_dp
    integer :: nEnergy = 300

    !> Refolding cell (columns are lattice vectors, Bohr); equals supercell if not given
    real(dp) :: refCell(3,3) = 0.0_dp

    !> Reciprocal vectors of the refolding cell (columns, 1/Bohr, including 2 pi)
    real(dp) :: refRecVec(3,3) = 0.0_dp

    !> Reciprocal vectors of the supercell (columns, 1/Bohr, including 2 pi)
    real(dp) :: supRecVec(3,3) = 0.0_dp

    !> Was an explicit refolding cell given?
    logical :: tRefold = .false.

    !> Cutoff for the refolding vectors: |g|^2 / 2 < gCutoff (Hartree)
    real(dp) :: gCutoff = 10.0_dp

    !> Refolding vectors (Cartesian, 1/Bohr)
    real(dp), allocatable :: gVec(:,:)

    !> q-points in Cartesian coordinates (1/Bohr), and the line they belong to
    real(dp), allocatable :: qPoint(:,:)
    integer, allocatable :: iLine(:)

    !> Number of separate q paths and index of the last q-point of each path
    integer :: nPath = 0
    integer, allocatable :: lastQ(:)

    !> Normalise the radial functions?
    logical :: tNormalise = .true.

    !> Tolerance (in fractional coordinates) for identifying k-points
    real(dp) :: kTol = 1.0e-5_dp

    !> Relation between the physical Bloch vector and the k-points of the file: K = kSign * k
    integer :: kSign = -1

    !> Check the sum rule of the Fourier weights for one k-point
    logical :: sumRuleCheck = .true.

    !> Largest |q| used for the sum rule check (1/Bohr)
    real(dp) :: sumRuleQMax = 8.0_dp

  end type TProgramVariables


  character(len=*), parameter :: version = "0.1"

  character(len=*), parameter :: rootTag = "unfold"

  character(len=*), parameter :: hsdInput = "unfold_in.hsd"

  character(len=*), parameter :: hsdParsedInput = "unfold_pin.hsd"

  integer, parameter :: parserVersion = 1

contains


  !> Reads the input file and all data from the DFTB+ run, which is independent of eigenvectors.
  subroutine TProgramVariables_init(this, stdOut)

    !> Instance
    type(TProgramVariables), intent(out) :: this

    !> Unit of standard output
    integer, intent(in) :: stdOut

    type(fnode), pointer :: hsdTree, root, tmp, detailed, child, node
    type(string) :: strBuffer, modifier
    integer :: inputVersion, ii
    real(dp) :: tmpWindow(2), gCutRaw
    real(dp) :: latVec(3,3)
    logical :: tAnyBasis

    call printDftbHeader(stdOut, '(UNFOLD ' // version // ')', 2025)

    call parseHSD(rootTag, hsdInput, hsdTree)
    call getChild(hsdTree, rootTag, root)
    write(stdOut, "(A)") "Interpreting input file '" // hsdInput // "'"

    call getChildValue(root, "InputVersion", inputVersion, parserVersion)
    if (inputVersion /= parserVersion) then
      call error("Version of input (" // i2c(inputVersion) // ") and parser ("&
          & // i2c(parserVersion) // ") do not match")
    end if

    call getChildValue(root, "Action", strBuffer, "Unfold")
    select case (trim(adjustl(char(strBuffer))))
    case ("Unfold", "unfold")
      this%action = actionUnfold
    case ("ListKPoints", "listkpoints")
      this%action = actionListK
    case ("CheckBasis", "checkbasis")
      this%action = actionCheck
    case default
      call detailedError(root, "Unknown action '" // trim(char(strBuffer)) // "'. Valid choices&
          & are Unfold, ListKPoints and CheckBasis.")
    end select

    ! Data from the preceeding DFTB+ run
    call getChildValue(root, "DetailedXML", strBuffer, "detailed.xml")
    call readHSDAsXML(unquote(char(strBuffer)), tmp)
    call getChild(tmp, "detailedout", detailed)
    call readDetailed(this%input, detailed)
    call destroyNode(tmp)
    if (.not. this%input%geo%tPeriodic) call error("Unfolding needs a periodic system.")

    call getChildValue(root, "EigenvecBin", strBuffer, "eigenvec.bin")
    this%eigVecBin = unquote(char(strBuffer))
    call getChildValue(root, "BandOut", strBuffer, "band.out")
    this%bandOut = unquote(char(strBuffer))
    call getChildValue(root, "Output", strBuffer, "unfold")
    this%outPrefix = unquote(char(strBuffer))

    ! Reciprocal vectors of the supercell
    latVec(:,:) = this%input%geo%latVecs
    call invert33(this%supRecVec, latVec)
    this%supRecVec(:,:) = 2.0_dp * pi * transpose(this%supRecVec)

    ! Energy mesh
    call getChild(root, "Energy", node, requested=.false.)
    if (associated(node)) then
      call getChildValue(node, "Window", tmpWindow, [-20.0_dp, 10.0_dp] * eV__Hartree,&
          & modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), energyUnits, child, tmpWindow)
      this%eMin = tmpWindow(1) * Hartree__eV
      this%eMax = tmpWindow(2) * Hartree__eV
      call getChildValue(node, "NrOfIntervals", this%nEnergy, 300)
      call getChildValue(node, "FermiLevel", this%fermiLevel, 0.0_dp, modifier=modifier,&
          & child=child)
      call convertUnitHsd(char(modifier), energyUnits, child, this%fermiLevel)
      this%fermiLevel = this%fermiLevel * Hartree__eV
      call getChildValue(node, "Broadening", this%broadening, 0.0_dp, modifier=modifier,&
          & child=child)
      call convertUnitHsd(char(modifier), energyUnits, child, this%broadening)
      this%broadening = this%broadening * Hartree__eV
    end if
    if (this%eMax <= this%eMin) call error("Energy window is empty.")
    if (this%nEnergy < 1) call error("NrOfIntervals must be positive.")

    ! Refolding cell
    call getChild(root, "RefoldingLatticeVectors", node, requested=.false., modifier=modifier)
    if (associated(node)) then
      this%tRefold = .true.
      call getChildValue(root, "RefoldingLatticeVectors", latVec, modifier=modifier, child=child)
      call convertUnitHsd(char(modifier), lengthUnits, child, latVec)
      this%refCell(:,:) = latVec
      if (abs(determinant33(this%refCell)) < 1.0e-8_dp) then
        call detailedError(child, "Singular refolding cell.")
      end if
    else
      this%refCell(:,:) = this%input%geo%latVecs
    end if
    call invert33(this%refRecVec, this%refCell)
    this%refRecVec(:,:) = 2.0_dp * pi * transpose(this%refRecVec)

    call getChildValue(root, "RefoldingGCutoff", gCutRaw, 10.0_dp, modifier=modifier,&
        & child=child)
    call convertUnitHsd(char(modifier), energyUnits, child, gCutRaw)
    this%gCutoff = gCutRaw
    call makeRefoldingVectors(this)

    call readQLines(this, root)

    call getChildValue(root, "NormalizeOrbitals", this%tNormalise, .true.)
    call getChildValue(root, "KPointSign", this%kSign, -1)
    call getChildValue(root, "SumRuleCheck", this%sumRuleCheck, .true.)
    call getChildValue(root, "SumRuleQMax", this%sumRuleQMax, 8.0_dp)

    tAnyBasis = (this%action /= actionListK)
    if (tAnyBasis) then
      call getChild(root, "Basis", tmp)
      call readBasis(this, tmp)
    end if

    call warnUnprocessedNodes(stdOut, root, .true.)
    call dumpHSD(hsdTree, hsdParsedInput)
    write(stdOut, "(A)") "Processed input written as HSD to '" // hsdParsedInput // "'"
    call destroyNode(hsdTree)

  end subroutine TProgramVariables_init


  !> Reads the relevant parts of detailed.xml.
  subroutine readDetailed(input, detailed)

    !> Resulting data
    type(TInput), intent(inout) :: input

    !> Root node of detailed.xml
    type(fnode), pointer :: detailed

    type(fnode), pointer :: tmp
    real(dp), allocatable :: kw(:,:)

    call getChildValue(detailed, "Identity", input%identity)
    call getChild(detailed, "Geometry", tmp)
    call readTGeometryHSD(tmp, input%geo)
    call getChildValue(detailed, "Real", input%tRealHam)
    call getChildValue(detailed, "NrOfKPoints", input%nKPoint)
    call getChildValue(detailed, "NrOfSpins", input%nSpin)
    call getChildValue(detailed, "NrOfStates", input%nState)
    call getChildValue(detailed, "NrOfOrbitals", input%nOrb)
    allocate(kw(4, input%nKPoint))
    call getChildValue(detailed, "KPointsAndWeights", kw)
    allocate(input%kPoint(3, input%nKPoint))
    input%kPoint(:,:) = kw(1:3, :)

  end subroutine readDetailed


  !> Generates the refolding vectors g (including g = 0 as first element).
  subroutine makeRefoldingVectors(this)

    !> Instance
    type(TProgramVariables), intent(inout) :: this

    real(dp) :: gMax, gg(3)
    real(dp), allocatable :: buf(:,:)
    integer :: nMax(3), i1, i2, i3, ng

    if (.not. this%tRefold) then
      allocate(this%gVec(3, 1))
      this%gVec(:,:) = 0.0_dp
      return
    end if

    gMax = sqrt(2.0_dp * this%gCutoff)
    ! Number of reciprocal vectors needed along each direction: |b_i . (cell^-1 row)|
    ! n_i <= gMax * |a_i| / (2 pi)
    do i1 = 1, 3
      nMax(i1) = ceiling(gMax * norm2(this%refCell(:, i1)) / (2.0_dp * pi))
    end do
    allocate(buf(3, product(2 * nMax + 1)))
    ng = 1
    buf(:, 1) = 0.0_dp
    do i3 = -nMax(3), nMax(3)
      do i2 = -nMax(2), nMax(2)
        do i1 = -nMax(1), nMax(1)
          if (i1 == 0 .and. i2 == 0 .and. i3 == 0) cycle
          gg(:) = this%refRecVec(:, 1) * i1 + this%refRecVec(:, 2) * i2 + this%refRecVec(:, 3) * i3
          if (norm2(gg) <= gMax) then
            ng = ng + 1
            buf(:, ng) = gg
          end if
        end do
      end do
    end do
    allocate(this%gVec(3, ng))
    this%gVec(:,:) = buf(:, 1:ng)

  end subroutine makeRefoldingVectors


  !> Reads the q-point lines in the same format as the Klines option of DFTB+. The coordinates are
  !! relative to the reciprocal vectors of the refolding cell (or the supercell, if none is given).
  !! Several paths can be given as separate "Path" blocks.
  subroutine readQLines(this, root)

    !> Instance
    type(TProgramVariables), intent(inout) :: this

    !> Root of the input
    type(fnode), pointer :: root

    type(fnode), pointer :: qNode, pathNode
    type(fNodeList), pointer :: paths
    integer :: iPath, nPaths

    call getChild(root, "QPoints", qNode)
    call getChildren(qNode, "Path", paths)
    nPaths = getLength(paths)
    if (nPaths == 0) then
      ! A single path directly given
      this%nPath = 1
      allocate(this%lastQ(0:1))
      this%lastQ(0) = 0
      call readOnePath(this, qNode, 0)
      this%lastQ(1) = size(this%qPoint, dim=2)
    else
      this%nPath = nPaths
      allocate(this%lastQ(0:nPaths))
      this%lastQ(0) = 0
      do iPath = 1, nPaths
        call getItem1(paths, iPath, pathNode)
        call readOnePath(this, pathNode, this%lastQ(iPath - 1))
        this%lastQ(iPath) = size(this%qPoint, dim=2)
      end do
    end if

  end subroutine readQLines


  !> Reads one path and appends the q-points.
  subroutine readOnePath(this, node, nOld)

    !> Instance
    type(TProgramVariables), intent(inout) :: this

    !> Node with the path
    type(fnode), pointer :: node

    !> Number of q-points already stored
    integer, intent(in) :: nOld

    type(TListIntR1) :: li
    type(TListRealR1) :: lr
    integer, allocatable :: nStep(:), tmpI(:,:)
    real(dp), allocatable :: pts(:,:), qNew(:,:)
    integer, allocatable :: lNew(:)
    real(dp) :: dd(3), qRel(3), qq(3)
    integer :: ii, jj, kk, nTot, ind

    call init(li)
    call init(lr)
    call getChildValue(node, "", 1, li, 3, lr)
    if (len(li) < 1) call detailedError(node, "At least one point must be given.")
    allocate(tmpI(1, len(li)))
    allocate(nStep(len(li)))
    call asArray(li, tmpI)
    nStep(:) = tmpI(1, :)
    allocate(pts(3, len(lr)))
    call asArray(lr, pts)
    call destruct(li)
    call destruct(lr)
    if (nStep(1) /= 1) call detailedError(node, "The first line of a path must have 1 point.")
    if (any(nStep < 1)) call detailedError(node, "Number of points must be positive.")

    nTot = sum(nStep)
    allocate(qNew(3, nOld + nTot))
    allocate(lNew(nOld + nTot))
    if (nOld > 0) then
      qNew(:, 1:nOld) = this%qPoint
      lNew(1:nOld) = this%iLine
      deallocate(this%qPoint, this%iLine)
    end if
    ind = nOld
    do ii = 1, size(nStep)
      do jj = 1, nStep(ii)
        if (ii == 1) then
          qRel(:) = pts(:, 1)
        else
          dd(:) = (pts(:, ii) - pts(:, ii - 1)) / real(nStep(ii), dp)
          qRel(:) = pts(:, ii - 1) + real(jj, dp) * dd
        end if
        ind = ind + 1
        qNew(:, ind) = matmul(this%refRecVec, qRel)
        lNew(ind) = ii
      end do
    end do
    call move_alloc(qNew, this%qPoint)
    call move_alloc(lNew, this%iLine)

  end subroutine readOnePath


  !> Reads the STO basis for all species.
  subroutine readBasis(this, node)

    !> Instance
    type(TProgramVariables), intent(inout) :: this

    !> Node containing the basis (same format as for waveplot)
    type(fnode), pointer :: node

    type(fnode), pointer :: spNode, orbNode
    type(fNodeList), pointer :: children
    real(dp), allocatable :: exps(:), coeffs(:)
    real(dp) :: cutoff, res
    integer :: iSp, iOrb, nOrb, ll
    character(len=:), allocatable :: spName

    call getChildValue(node, "Resolution", res, 0.01_dp)
    allocate(this%basis(this%input%geo%nSpecies))
    do iSp = 1, this%input%geo%nSpecies
      spName = trim(this%input%geo%speciesNames(iSp))
      call getChild(node, spName, spNode)
      call getChildValue(spNode, "AtomicNumber", this%basis(iSp)%atomicNumber)
      call getChildren(spNode, "Orbital", children)
      nOrb = getLength(children)
      if (nOrb < 1) call detailedError(spNode, "Missing orbital definitions")
      allocate(this%basis(iSp)%shells(nOrb))
      do iOrb = 1, nOrb
        call getItem1(children, iOrb, orbNode)
        call getChildValue(orbNode, "AngularMomentum", ll)
        if (ll < 0 .or. ll > 3) call detailedError(orbNode, "Angular momentum must be 0..3")
        call getChildValue(orbNode, "Cutoff", cutoff)
        call readRealList(orbNode, "Exponents", exps)
        call readRealList(orbNode, "Coefficients", coeffs)
        if (mod(size(coeffs), size(exps)) /= 0) then
          call detailedError(orbNode, "Number of coefficients incompatible with number of&
              & exponents")
        end if
        call TStoShell_init(this%basis(iSp)%shells(iOrb), ll, cutoff, exps, coeffs,&
            & this%tNormalise)
      end do
    end do

  end subroutine readBasis


  !> Reads a free-length list of reals from a node.
  subroutine readRealList(node, name, vals)

    !> Parent node
    type(fnode), pointer :: node

    !> Name of the child
    character(len=*), intent(in) :: name

    !> Values on return
    real(dp), allocatable, intent(out) :: vals(:)

    type(TListReal) :: buffer
    type(fnode), pointer :: child

    call init(buffer)
    call getChildValue(node, name, buffer, child=child)
    if (len(buffer) == 0) call detailedError(child, "Missing values for " // name)
    allocate(vals(len(buffer)))
    call asArray(buffer, vals)
    call destruct(buffer)

  end subroutine readRealList

end module unfold_initunfold
