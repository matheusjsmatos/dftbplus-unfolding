!--------------------------------------------------------------------------------------------------!
!  DFTB+: general package for performing fast atomistic simulations                                !
!  Copyright (C) 2006 - 2025  DFTB+ developers group                                               !
!                                                                                                  !
!  See the LICENSE file for terms of usage and distribution.                                       !
!--------------------------------------------------------------------------------------------------!

#:include 'common.fypp'

!> Radial parts of the (Slater-type) basis functions, their spherical Bessel transforms and the
!! real tesseral harmonics needed to Fourier transform the DFTB atomic orbitals.
!!
!! A basis function of an atom is phi(r) = R(r) Y_lm(r^), where R(r) is the STO expansion stored in
!! the wfc.*.hsd files (the same data the waveplot tool uses):
!!
!!   R(r) = sum_j exp(-alpha_j r) sum_i a_ij r^(l+i-1)
!!
!! The 3D Fourier transform with the unitary convention
!!
!!   phi~(q) = (2 pi)^(-3/2) int phi(r) exp(-i q.r) d^3r = (-i)^l Y_lm(q^) F_l(|q|)
!!
!! contains the radial transform  F_l(q) = sqrt(2/pi) int r^2 R(r) j_l(q r) dr.
module unfold_radial
  use dftbp_common_accuracy, only : dp
  implicit none (type, external)

  private
  public :: TStoShell, TRadialTable, TRadialTable_init, TStoShell_init, evalSto, getRadialFT
  public :: realTessY, sphericalBessel


  !> One radial shell (one angular momentum) of a species.
  type TStoShell

    !> Angular momentum
    integer :: ll = 0

    !> Cutoff radius, beyond which the function is taken as zero (Bohr)
    real(dp) :: cutoff = 0.0_dp

    !> Positive exponents of the exponentials (1/Bohr)
    real(dp), allocatable :: alpha(:)

    !> Polynomial coefficients (nPow, nAlpha)
    real(dp), allocatable :: aa(:,:)

    !> Factor, which normalises int R^2 r^2 dr to one (1 if normalisation is switched off)
    real(dp) :: scale = 1.0_dp

    !> Norm int R^2 r^2 dr of the unscaled function (for diagnostics)
    real(dp) :: norm = 1.0_dp

  end type TStoShell


  !> Tabulated radial Fourier transforms F_l(q) for all shells of one species.
  type TRadialTable

    !> Spacing of the q grid (1/Bohr)
    real(dp) :: dq = 0.0_dp

    !> Number of q grid intervals
    integer :: nq = 0

    !> Radial transforms (0:nq, nShell)
    real(dp), allocatable :: val(:,:)

  end type TRadialTable


  !> Number of radial quadrature points (must be odd for Simpson's rule)
  integer, parameter :: nRadial = 6001

contains


  !> Initialises a shell and computes its norm.
  subroutine TStoShell_init(this, ll, cutoff, alpha, coeffs, normalise)

    !> Instance
    type(TStoShell), intent(out) :: this

    !> Angular momentum
    integer, intent(in) :: ll

    !> Cutoff radius (Bohr)
    real(dp), intent(in) :: cutoff

    !> Exponents (positive values)
    real(dp), intent(in) :: alpha(:)

    !> Flat list of coefficients; first the nPow coefficients of alpha(1), then those of alpha(2)...
    real(dp), intent(in) :: coeffs(:)

    !> Whether the radial function should be normalised to int R^2 r^2 dr = 1
    logical, intent(in) :: normalise

    real(dp), allocatable :: rr(:), integrand(:)
    integer :: nPow, ii

    nPow = size(coeffs) / size(alpha)
    this%ll = ll
    this%cutoff = cutoff
    allocate(this%alpha(size(alpha)))
    this%alpha(:) = alpha
    allocate(this%aa(nPow, size(alpha)))
    this%aa(:,:) = reshape(coeffs, [nPow, size(alpha)])
    this%scale = 1.0_dp

    allocate(rr(nRadial), integrand(nRadial))
    call radialGrid(cutoff, rr)
    do ii = 1, nRadial
      integrand(ii) = (evalSto(this, rr(ii)) * rr(ii))**2
    end do
    this%norm = simpson(integrand, rr(2) - rr(1))
    if (normalise .and. this%norm > 0.0_dp) this%scale = 1.0_dp / sqrt(this%norm)

  end subroutine TStoShell_init


  !> Evaluates R(r) of a shell (including the normalisation factor) at a given distance.
  function evalSto(shell, rr) result(val)

    !> Shell
    type(TStoShell), intent(in) :: shell

    !> Distance (Bohr)
    real(dp), intent(in) :: rr

    !> Value R(r)
    real(dp) :: val

    real(dp) :: rPow, poly
    integer :: ii, jj, nPow

    val = 0.0_dp
    if (rr > shell%cutoff) return
    nPow = size(shell%aa, dim=1)
    do jj = 1, size(shell%alpha)
      poly = 0.0_dp
      if (shell%ll == 0) then
        rPow = 1.0_dp
      else
        rPow = rr**shell%ll
      end if
      do ii = 1, nPow
        poly = poly + shell%aa(ii, jj) * rPow
        rPow = rPow * rr
      end do
      val = val + poly * exp(-shell%alpha(jj) * rr)
    end do
    val = val * shell%scale

  end function evalSto


  !> Tabulates the radial Fourier transforms F_l(q) of all shells of a species on a uniform q grid.
  subroutine TRadialTable_init(this, shells, qMax, dq)

    !> Instance
    type(TRadialTable), intent(out) :: this

    !> Shells of the species
    type(TStoShell), intent(in) :: shells(:)

    !> Largest q needed (1/Bohr)
    real(dp), intent(in) :: qMax

    !> q spacing (1/Bohr)
    real(dp), intent(in) :: dq

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), allocatable :: rr(:), integrand(:), rad(:)
    real(dp) :: qq
    integer :: iShell, iq, ir

    this%dq = dq
    this%nq = ceiling(qMax / dq) + 2
    allocate(this%val(0:this%nq, size(shells)))
    allocate(rr(nRadial), integrand(nRadial), rad(nRadial))

    do iShell = 1, size(shells)
      call radialGrid(shells(iShell)%cutoff, rr)
      do ir = 1, nRadial
        rad(ir) = rr(ir)**2 * evalSto(shells(iShell), rr(ir))
      end do
      do iq = 0, this%nq
        qq = real(iq, dp) * dq
        do ir = 1, nRadial
          integrand(ir) = rad(ir) * sphericalBessel(shells(iShell)%ll, qq * rr(ir))
        end do
        this%val(iq, iShell) = sqrt(2.0_dp / pi) * simpson(integrand, rr(2) - rr(1))
      end do
    end do

  end subroutine TRadialTable_init


  !> Returns F_l(q) of a shell by linear interpolation in the table.
  function getRadialFT(this, iShell, qq) result(val)

    !> Table
    type(TRadialTable), intent(in) :: this

    !> Shell index
    integer, intent(in) :: iShell

    !> Wave vector length (1/Bohr)
    real(dp), intent(in) :: qq

    !> Interpolated transform
    real(dp) :: val

    integer :: ind
    real(dp) :: frac

    ind = int(qq / this%dq)
    if (ind >= this%nq) then
      val = 0.0_dp
      return
    end if
    frac = qq / this%dq - real(ind, dp)
    val = (1.0_dp - frac) * this%val(ind, iShell) + frac * this%val(ind + 1, iShell)

  end function getRadialFT


  !> Spherical Bessel function j_l(x) for l = 0..3.
  function sphericalBessel(ll, xx) result(val)

    !> Order
    integer, intent(in) :: ll

    !> Argument (non-negative)
    real(dp), intent(in) :: xx

    !> Value
    real(dp) :: val

    real(dp) :: sx, cx, term, x2
    integer :: kk, ii

    if (xx < 2.0_dp) then
      ! Power series  j_l(x) = x^l/(2l+1)!! sum_k (-x^2/2)^k / (k! (2l+3)(2l+5)...(2l+2k+1))
      x2 = -0.5_dp * xx**2
      term = 1.0_dp
      do ii = 1, ll
        term = term * xx / real(2 * ii + 1, dp)
      end do
      val = term
      do kk = 1, 25
        term = term * x2 / (real(kk, dp) * real(2 * ll + 2 * kk + 1, dp))
        val = val + term
        if (abs(term) < 1.0e-17_dp * abs(val)) exit
      end do
    else
      sx = sin(xx)
      cx = cos(xx)
      select case (ll)
      case (0)
        val = sx / xx
      case (1)
        val = sx / xx**2 - cx / xx
      case (2)
        val = (3.0_dp / xx**3 - 1.0_dp / xx) * sx - 3.0_dp * cx / xx**2
      case (3)
        val = (15.0_dp / xx**4 - 6.0_dp / xx**2) * sx - (15.0_dp / xx**3 - 1.0_dp / xx) * cx
      case default
        val = 0.0_dp
      end select
    end if

  end function sphericalBessel


  !> Real tesseral harmonics for l = 0..3 in the m-ordering and sign convention used by DFTB+
  !! (identical to the one of the waveplot tool). The vector needs not to be normalised.
  function realTessY(ll, mm, vec) result(rty)

    !> Angular momentum (0 <= l <= 3)
    integer, intent(in) :: ll

    !> Magnetic quantum number (-l..l)
    integer, intent(in) :: mm

    !> Direction
    real(dp), intent(in) :: vec(3)

    !> Value of the harmonic
    real(dp) :: rty

    real(dp) :: rr, xx, yy, zz

    rr = norm2(vec)
    xx = vec(1)
    yy = vec(2)
    zz = vec(3)
    rty = 0.0_dp
    if (ll == 0) then
      rty = 0.2820947917738782_dp
      return
    end if
    if (rr < epsilon(1.0_dp)) return

    select case (ll)
    case(1)
      select case(mm)
      case(-1)
        rty = 0.4886025119029198_dp * yy / rr
      case(0)
        rty = 0.4886025119029198_dp * zz / rr
      case(1)
        rty = 0.4886025119029198_dp * xx / rr
      end select
    case(2)
      select case(mm)
      case(-2)
        rty = 1.092548430592079_dp * xx * yy / rr**2
      case(-1)
        rty = 1.092548430592079_dp * yy * zz / rr**2
      case(0)
        rty = -0.3153915652525200_dp * (-2.0_dp * zz**2 + xx**2 + yy**2) / rr**2
      case(1)
        rty = 1.092548430592079_dp * xx * zz / rr**2
      case(2)
        rty = 0.5462742152960395_dp * (xx**2 - yy**2) / rr**2
      end select
    case(3)
      select case (mm)
      case(-3)
        rty = 0.5900435899266435_dp * yy * (3.0_dp * xx**2 - yy**2) / rr**3
      case(-2)
        rty = 2.890611442640554_dp * xx * yy * zz / rr**3
      case(-1)
        rty = -0.4570457994644658_dp * (-4.0_dp * zz**2 + xx**2 + yy**2) * yy / rr**3
      case(0)
        rty = -0.3731763325901155_dp * zz * (-2.0_dp * zz**2 + 3.0_dp * xx**2 + 3.0_dp * yy**2)&
            & / rr**3
      case(1)
        rty = -0.4570457994644658_dp * (-4.0_dp * zz**2 + xx**2 + yy**2) * xx / rr**3
      case(2)
        rty = 1.445305721320277_dp * zz * (xx**2 - yy**2) / rr**3
      case(3)
        rty = 0.5900435899266435_dp * xx * (xx**2 - 3.0_dp * yy**2) / rr**3
      end select
    end select

  end function realTessY


  !> Uniform radial grid from 0 to the cutoff.
  subroutine radialGrid(cutoff, rr)

    !> Cutoff radius
    real(dp), intent(in) :: cutoff

    !> Grid points (size nRadial)
    real(dp), intent(out) :: rr(:)

    integer :: ii

    do ii = 1, size(rr)
      rr(ii) = cutoff * real(ii - 1, dp) / real(size(rr) - 1, dp)
    end do

  end subroutine radialGrid


  !> Simpson integration of equidistant data (odd number of points).
  function simpson(ff, hh) result(res)

    !> Function values
    real(dp), intent(in) :: ff(:)

    !> Grid spacing
    real(dp), intent(in) :: hh

    !> Integral
    real(dp) :: res

    integer :: nn

    nn = size(ff)
    res = ff(1) + ff(nn) + 4.0_dp * sum(ff(2:nn-1:2)) + 2.0_dp * sum(ff(3:nn-2:2))
    res = res * hh / 3.0_dp

  end function simpson

end module unfold_radial
