! ============================================================
! Module: hf_kernels
! Builds the one- and two-electron integrals needed for the
! Hartree-Fock SCF procedure for closed-shell systems with
! contracted s-type Gaussian basis functions.
! ============================================================

MODULE hf_kernels

    IMPLICIT NONE

CONTAINS

FUNCTION distance(r1, r2) RESULT(dist)
    
    ! Compute the Euclidean distance between two 3D points

    double precision, dimension(3), intent(in) :: r1, r2
    double precision :: dist

    dist = sqrt((r1(1) -r2(1))**2.0d0 + (r1(2) - r2(2))**2.0d0 + (r1(3) - r2(3))**2.0d0)

END FUNCTION distance

FUNCTION nuclear_repulsion(natoms, Z, atom_coords) RESULT(enuc)
    
    ! Compute the nuclear-nuclear repulsion energy from the atomic numbers and coordinates

    integer, intent(in) :: natoms
    integer, dimension(:), intent(in) :: Z
    double precision, dimension(:,:), intent(in) :: atom_coords
    double precision :: enuc, rij
    integer :: i, j

    enuc = 0.0d0

    ! Loop over all unique pairs of nuclei to compute the repulsion energy

    do i = 1, natoms
        do j = i + 1, natoms
            rij = distance(atom_coords(i, :), atom_coords(j, :))
            enuc = enuc + (dble(Z(i)) * dble(Z(j))) / rij
        end do
    end do

END FUNCTION nuclear_repulsion

SUBROUTINE fill_symmetry(nbasis, ERI, mu, nu, lambda, sigma, val)

    ! Fill all equivalent permutations of a two-electron integral

    integer, intent(in) :: nbasis, mu, nu, lambda, sigma
    double precision, intent(in) :: val
    double precision, dimension(nbasis, nbasis, nbasis, nbasis), intent(inout) :: ERI

    ERI(mu, nu, lambda, sigma) = val
    ERI(mu, nu, sigma, lambda) = val
    ERI(nu, mu, lambda, sigma) = val
    ERI(nu, mu, sigma, lambda) = val
    ERI(lambda, sigma, mu, nu) = val
    ERI(lambda, sigma, nu, mu) = val
    ERI(sigma, lambda, mu, nu) = val
    ERI(sigma, lambda, nu, mu) = val

END SUBROUTINE fill_symmetry

SUBROUTINE boys_func(arg, f0, f1)

    ! Evaluate Boys functions F0 and F1, including the small-argument limit

    double precision, intent(in) :: arg
    double precision, intent(out) :: f0, f1
    double precision, parameter :: pi = acos(-1.0d0)

    ! For very small arguments, use the limiting values F0(0)=1 and F1(0)=1/3
    
    if (arg < 1.0d-10) then
        F0 = 1.0d0
        F1 = 1.0d0 / 3.0d0
    else
        F0 = 0.5d0 * sqrt(pi / arg) * derf(sqrt(arg))
        F1 = (F0 - exp(-arg)) / (2.0d0 * arg)
    end if

END SUBROUTINE boys_func

SUBROUTINE compute_overlap_matrix(nbasis, n_primitives, exponents, coefficients, atom_coords, atom_map, S)

    ! Build the overlap matrix S for contracted s-type Gaussian basis functions.
    ! Only unique basis-function pairs are evaluated explicitly and matrix symmetry is used.

    integer, intent(in) :: nbasis
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:,:), intent(in) :: exponents, coefficients, atom_coords
    integer, dimension(:), intent(in) :: atom_map
    double precision, dimension(nbasis, nbasis), intent(out) :: S

    integer :: mu, nu, k, l
    double precision :: alpha, beta, zeta, xi, dist_sq, S_val
    double precision, parameter :: pi = acos(-1.0d0)

    ! Initialize the overlap matrix to zero

    S = 0.0d0
    
    ! Loop over unique basis-function pairs (mu,nu)
    
    do mu = 1, nbasis
        
        do nu = 1, mu
            
            ! Squared distance between the centers of basis functions mu and nu

            dist_sq = distance(atom_coords(atom_map(mu), :), atom_coords(atom_map(nu), :))**2.0d0

            S_val = 0.0d0

                ! Sum over primitive Gaussians belonging to basis functions mu and nu
                
                do k = 1, n_primitives(mu)
                    
                    alpha = exponents(mu, k)

                    do l = 1, n_primitives(nu)
                        
                        ! Composite exponents and prefactor for the primitive Gaussian pair

                        beta = exponents(nu, l)
                        zeta = alpha + beta
                        xi = (alpha * beta) / zeta
                        S_val = S_val + coefficients(mu, k) * coefficients(nu, l) * exp(-xi * dist_sq) * ((pi / zeta) ** 1.5d0)
                    
                    end do
                end do

                ! Use symmetry S(mu,nu) = S(nu,mu)

                S(mu, nu) = S_val
                S(nu, mu) = S_val

            end do
        end do

END SUBROUTINE compute_overlap_matrix

SUBROUTINE compute_kinetic_matrix(nbasis, n_primitives, exponents, coefficients, atom_coords, atom_map, T)

    ! Build the kinetic-energy matrix T for contracted s-type Gaussian basis functions.
    ! Only unique basis-function pairs are evaluated explicitly and matrix symmetry is used.

    integer, intent(in) :: nbasis
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:,:), intent(in) :: exponents, coefficients, atom_coords
    integer, dimension(:), intent(in) :: atom_map
    double precision, dimension(nbasis, nbasis), intent(out) :: T

    integer :: mu, nu, k, l
    double precision :: alpha, beta, zeta, xi, dist_sq, T_val, S_prim, T_prim
    double precision, parameter :: pi = acos(-1.0d0)

    ! Initialize the kinetic-energy matrix

    T = 0.0d0

    ! Loop over unique basis-function pairs (mu,nu)

    do mu = 1, nbasis
        
        do nu = 1, mu
            
            ! Squared distance between the centers of basis functions mu and nu

            dist_sq = distance(atom_coords(atom_map(mu), :), atom_coords(atom_map(nu), :))**2.0d0

            T_val = 0.0d0

                ! Sum over primitive Gaussians belonging to basis functions mu and nu
    
                do k = 1, n_primitives(mu)
                    
                    alpha = exponents(mu, k)

                    do l = 1, n_primitives(nu)
                        
                        ! Composite exponents for the primitive Gaussian pair

                        beta = exponents(nu, l)
                        zeta = alpha + beta
                        xi = (alpha * beta) / zeta

                        ! Primitive overlap contribution reused in the kinetic-energy expression

                        S_prim = exp(-xi * dist_sq) * ((pi / zeta) ** 1.5d0)
                        T_prim = xi * (3.0d0 - 2.0d0 * xi * dist_sq) * S_prim
                        T_val = T_val + coefficients(mu, k) * coefficients(nu, l) * T_prim

                    end do
                end do

                ! Use symmetry T(mu,nu) = T(nu,mu)

                T(mu, nu) = T_val
                T(nu, mu) = T_val

            end do
        end do

END SUBROUTINE compute_kinetic_matrix

SUBROUTINE compute_nuclear_atraction_matrix(nbasis, n_primitives, exponents, coefficients, atom_coords, atom_map, Z, V )

    ! Build the electron-nucleus attraction matrix V for contracted s-type Gaussian basis functions.
    ! Contributions from all nuclei are accumulated for each basis-function pair.

    integer, intent(in) :: nbasis
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:,:), intent(in) :: exponents, coefficients, atom_coords
    integer, dimension(:), intent(in) :: atom_map, Z
    double precision, dimension(nbasis, nbasis), intent(out) :: V

    integer :: mu, nu, k, l, A
    double precision :: alpha, beta, zeta, xi, dist_sq, center_P(3), V_val, V_prim, S_prim, rIP_sq, arg, f0, f1
    double precision, parameter :: pi = acos(-1.0d0)

    ! Initialize the nuclear-attraction matrix

    V = 0.0d0

    ! Loop over unique basis-function pairs (mu,nu)

    do mu = 1, nbasis
        
        do nu = 1, mu
            
            ! Squared distance between the centers of basis functions mu and nu

            dist_sq = distance(atom_coords(atom_map(mu), :), atom_coords(atom_map(nu), :))**2.0d0
            V_val = 0.0d0

            ! Sum over primitive Gaussians belonging to basis functions mu and nu

            do k = 1, n_primitives(mu)
                
                alpha = exponents(mu, k)

                do l = 1, n_primitives(nu)
                    
                    beta = exponents(nu, l)
                    zeta = alpha + beta
                    xi = (alpha * beta) / zeta

                    ! Gaussian product center and primitive overlap contribution

                    center_P = (alpha * atom_coords(atom_map(mu), :) + beta * atom_coords(atom_map(nu), :)) / zeta
                    S_prim = exp(-xi * dist_sq) * ((pi / zeta) ** 1.5d0)
                    V_prim = 0.0d0

                    ! Sum the attraction contribution from all nuclei

                    do A = 1, size(Z)
                        
                        rIP_sq = distance(atom_coords(A, :), center_P)**2.0d0
                        arg = zeta * rIP_sq

                        ! Boys function value for the electron-nucleus attraction term

                        call boys_func(arg, f0, f1)

                        V_prim = V_prim - 2.0d0 * dble(Z(A)) * sqrt(zeta / pi) * S_prim * f0
                    
                    end do
                
                V_val = V_val + coefficients(mu, k) * coefficients(nu, l) * V_prim
                
                end do
            end do

        ! Use symmetry V(mu,nu) = V(nu,mu)

        V(mu, nu) = V_val
        V(nu, mu) = V_val
        
        end do
    end do

END SUBROUTINE compute_nuclear_atraction_matrix
                        
SUBROUTINE compute_eri(nbasis, n_primitives, exponents, coefficients, atom_coords, atom_map, eri)

    ! Build the two-electron integral tensor for contracted s-type Gaussian basis functions.
    ! Only unique ERIs are computed explicitly; all symmetry-equivalent permutations are filled afterward.

    integer, intent(in) :: nbasis
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:,:), intent(in) :: exponents, coefficients, atom_coords
    integer, dimension(:), intent(in) :: atom_map
    double precision, dimension(nbasis, nbasis, nbasis, nbasis), intent(out) :: eri

    integer :: mu, nu, lambda, sigma, k, l, m, n, munu, lasi
    double precision :: k1, k2, zeta1, zeta2, xi1, xi2, rho, center_P(3), center_Q(3), center_PQ, eri_prim, eri_val, f0, f1
    double precision, parameter :: pi = acos(-1.0d0)

    ! Initialize the two-electron integral tensor

    eri = 0.0d0

    ! Loop over unique basis-function quartets (mu,nu,lambda,sigma)
    do mu = 1, nbasis
        
        do nu = 1, mu
            
            ! Composite index used to keep only unique ERIs

            munu = mu * (mu -1) / 2 + nu

            do lambda = 1, nbasis
                
                do sigma = 1, lambda

                    lasi = lambda * (lambda -1) / 2 + sigma

                    ! Skip quartets related by ERI symmetry

                    if (munu < lasi) cycle

                    eri_val = 0.0d0

                    ! Sum over primitive Gaussian quartets

                    do k = 1, n_primitives(mu)
                        
                        do l = 1, n_primitives(nu)

                            ! Parameters for the first Gaussian pair (mu,nu)

                            zeta1 = exponents(mu, k) + exponents(nu, l)
                            xi1 = exponents(mu, k) * exponents(nu, l) / zeta1
                            center_P = (exponents(mu, k) * atom_coords(atom_map(mu), :) + exponents(nu, l) * atom_coords(atom_map(nu), :)) / zeta1
                            k1 = sqrt(2.0d0) * (pi**1.25d0 / zeta1) * exp(-xi1 * distance(atom_coords(atom_map(mu), :), atom_coords(atom_map(nu), :))**2.0d0)

                            do m = 1, n_primitives(lambda)
                              
                                do n = 1, n_primitives(sigma)

                                    ! Parameters for the second Gaussian pair (lambda,sigma)

                                    zeta2 = exponents(lambda, m) + exponents(sigma, n)
                                    xi2 = exponents(lambda, m) * exponents(sigma, n) / zeta2
                                    center_Q = (exponents(lambda, m) * atom_coords(atom_map(lambda), :) + exponents(sigma, n) * atom_coords(atom_map(sigma), :)) / zeta2
                                    k2 = sqrt(2.0d0) * (pi**1.25d0 / zeta2) * exp(-xi2 * distance(atom_coords(atom_map(lambda), :), atom_coords(atom_map(sigma), :)) **2.0d0)

                                    ! Distance and Boys-function argument for the primitive ERI

                                    rho = (zeta1 * zeta2) / (zeta1 + zeta2)
                                    center_PQ = distance(center_P, center_Q) ** 2

                                    call boys_func(rho * center_PQ, f0, f1)
                                    eri_prim = (1.0d0 / sqrt(zeta1 + zeta2)) * k1 * k2 * f0

                                    eri_val = eri_val + coefficients(mu, k) * coefficients(nu, l) * coefficients(lambda, m) * coefficients(sigma, n) * eri_prim

                                end do
                            end do
                        end do
                    end do
                    
                    ! Fill all symmetry-related permutations of the current ERI

                    call fill_symmetry(nbasis, eri, mu, nu, lambda, sigma, eri_val)
                
                end do
            end do
        end do
    end do

END SUBROUTINE compute_eri

END MODULE hf_kernels