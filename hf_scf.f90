! ============================================================
! hf_scf
! Routines for the Hartree-Fock SCF procedure:
! orthogonalization of the AO basis, Fock-matrix construction,
! solution of the Roothaan equations, density-matrix building,
! and electronic-energy evaluation.
! ============================================================

MODULE hf_scf

    IMPLICIT NONE

CONTAINS

SUBROUTINE orthogonalizer(nbasis, S, X)

    ! Build the symmetric orthogonalization matrix X = S^(-1/2)
    ! by diagonalizing the overlap matrix S.

    integer, intent(in) :: nbasis
    double precision, dimension(:,:), intent(in) :: S
    double precision, dimension(nbasis, nbasis), intent(out) :: X

    integer :: i, j, k, info
    double precision, dimension(nbasis) :: eigenvalues
    double precision, dimension(nbasis, nbasis) :: eigenvectors
    integer :: lwork
    double precision, dimension(3 * nbasis - 1):: work

    ! Copy S because the LAPACK diagonalizer overwrites its input

    lwork = 3 * nbasis - 1
    eigenvectors = S

    ! Diagonalize the overlap matrix: S = L * diag(eigenvalues) * L^T

    call dsyev('V', 'U', nbasis, eigenvectors, nbasis, eigenvalues, work, lwork, info)

    if (info /= 0) then
        print *, 'Error in dsyev: Diagonalization failed with info = ', info
        stop
    end if

    ! Reconstruct X = S^(-1/2) from the eigenvectors and eigenvalues of S

    X = 0.0d0
    do i = 1, nbasis
        
        do j = 1, nbasis
        
            do k = 1, nbasis
        
                ! Ignore very small eigenvalues to avoid numerical instability

                if (eigenvalues(k) > 1.0d-12) then        
        
                    X(i, j) = X(i, j) + eigenvectors(i, k) * (1.0d0 / sqrt(eigenvalues(k))) * eigenvectors(j, k)
        
                end if
        
            end do
        end do
    end do

END SUBROUTINE orthogonalizer

SUBROUTINE compute_g_matrix(nbasis, P, ERI, G)

    ! Build the two-electron contribution G to the Fock matrix from
    ! the current density matrix P and the two-electron integrals.

    integer, intent(in) :: nbasis
    double precision, dimension(:,:), intent(in) :: P
    double precision, dimension(:,:,:,:), intent(in) :: ERI
    double precision, dimension(nbasis, nbasis), intent(out) :: G

    integer :: i, j, k, l

    ! Initialize the two-electron part of the Fock matrix

    G = 0.0d0

    ! Contract the density matrix with Coulomb and exchange terms

    do i = 1, nbasis
        
        do j = 1, nbasis
        
            do k = 1, nbasis
        
                do l = 1, nbasis
        
                    G(i, j) = G(i, j) + P(k, l) * (ERI(i, j, k, l) - 0.5d0 * ERI(i, k, j, l))
        
                end do
            end do
        end do
    end do   

END SUBROUTINE compute_g_matrix

SUBROUTINE solve_fock_equations(nbasis, F, X, C, epsilon)

    ! Solve the Roothaan equations by transforming the Fock matrix
    ! to the orthonormal basis, diagonalizing it, and back-transforming
    ! the molecular-orbital coefficients.

    INTEGER, intent(in) :: nbasis
    double precision, dimension(:,:), intent(in) :: F, X
    double precision, dimension(nbasis, nbasis), intent(out) :: C
    double precision, dimension(nbasis), intent(out) :: epsilon

    double precision, dimension(nbasis, nbasis) :: F_prime, C_prime, temp
    double precision, dimension(3*nbasis-1) :: work
    INTEGER :: lwork, info

    ! Transform Fock Matrix: F' = X^T * F * X
    ! First: temp = F * X
    temp = matmul(F, X)
    ! Second: F' = X^T * temp
    F_prime = matmul(transpose(X), temp)

    ! Diagonalize F' in the orthonormal basis using LAPACK

    C_prime = F_prime
    lwork = 3*nbasis - 1
    call dsyev('V', 'U', nbasis, C_prime, nbasis, epsilon, work, lwork, info)

    if (info /= 0) then
        print *, "Error: Fock diagonalization failed."
        stop
    end if

    ! Back-transform coefficients: C = X * C'

    C = matmul(X, C_prime)

END SUBROUTINE solve_fock_equations

SUBROUTINE compute_density_matrix(nbasis, nocc, C, P)
    
    ! Build the closed-shell density matrix from the occupied
    ! molecular-orbital coefficients.

    integer, intent(in) :: nbasis, nocc
    double precision, dimension(:,:), intent(in) :: C
    double precision, dimension(nbasis, nbasis), intent(out) :: P

    integer :: mu, nu, i

    P = 0.0d0
    
    ! Sum only over occupied orbitals

    do mu = 1, nbasis
        
        do nu = 1, nbasis
        
            do i = 1, nocc
        
                P(mu, nu) = P(mu, nu) + 2.0d0 * C(mu, i) * C(nu, i)
        
            end do
        end do
    end do
END SUBROUTINE compute_density_matrix

SUBROUTINE calculate_electronic_energy(nbasis, P, Hcore, F, energy)

    ! evaluate the hartree-fock electronic energy from the density,
    ! core hamiltonian, and fock matrix.

    integer, intent(in) :: nbasis
    double precision, dimension(:,:), intent(in) :: P, Hcore, F
    double precision, intent(out) :: energy
    
    integer :: mu, nu
    
    energy = 0.0d0
    
    do mu = 1, nbasis
    
        do nu = 1, nbasis
    
            ! e = 0.5 * sum( P * (Hcore + F) )
    
            energy = energy + 0.5d0 * P(nu, mu) * (Hcore(mu, nu) + F(mu, nu))
    
        end do
    end do

END SUBROUTINE calculate_electronic_energy

END MODULE hf_scf