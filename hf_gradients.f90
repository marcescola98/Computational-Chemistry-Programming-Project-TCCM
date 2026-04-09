! ============================================================
! hf_gradients
! Routines to evaluate the different contributions to the
! Hartree-Fock energy gradient for closed-shell systems with
! contracted s-type Gaussian basis functions.
! ============================================================

MODULE hf_gradients

    IMPLICIT NONE

CONTAINS

SUBROUTINE boys_func(arg, f0, f1)

    ! Evaluate Boys functions F0 and F1.
    ! For very small arguments, use the limiting values F0(0)=1 and F1(0)=1/3.

    double precision, intent(in) :: arg
    double precision, intent(out) :: f0, f1
    double precision, parameter :: pi = acos(-1.0d0)

    ! Special treatment avoids division by very small numbers

    if (arg < 1.0d-10) then
        f0 = 1.0d0
        f1 = 1.0d0 / 3.0d0
    else
        f0 = 0.5d0 * sqrt(pi / arg) * derf(sqrt(arg))
        f1 = (f0 - exp(-arg)) / (2.0d0 * arg)
    end if

END SUBROUTINE boys_func


SUBROUTINE compute_vnn_gradient(natoms, Z, atom_coords, vnn_grad)

    ! compute the nuclear-nuclear repulsion contribution to the gradient.

    integer, intent(in) :: natoms
    integer, dimension(:), intent(in) :: z
    double precision, dimension(:,:), intent(in) :: atom_coords
    double precision, dimension(natoms, 3), intent(out) :: vnn_grad

    integer :: i, j
    double precision :: rij, dx, dy, dz, factor

    ! accumulate the derivative with respect to each atomic coordinate

    vnn_grad = 0.0d0

    ! loop over all atom pairs i /= j

    do i = 1, natoms
        do j = 1, natoms
            if (i == j) cycle
            
            dx = atom_coords(i, 1) - atom_coords(j, 1)
            dy = atom_coords(i, 2) - atom_coords(j, 2)
            dz = atom_coords(i, 3) - atom_coords(j, 3)
            rij = sqrt(dx**2 + dy**2 + dz**2)

            ! derivative of the pairwise nuclear repulsion term
                        
            factor = - dble(z(i)) * dble(z(j)) / (rij**3)

            vnn_grad(i, 1) = vnn_grad(i, 1) + factor * dx
            vnn_grad(i, 2) = vnn_grad(i, 2) + factor * dy
            vnn_grad(i, 3) = vnn_grad(i, 3) + factor * dz
        
        end do
    end do

END SUBROUTINE compute_vnn_gradient


SUBROUTINE compute_w_matrix(nbasis, nocc, eigenvalues, C, W)

    ! build the energy-weighted density matrix w used in the overlap
    ! (lagrangian) contribution to the gradient.

    integer, intent(in) :: nbasis, nocc
    double precision, dimension(:), intent(in) :: eigenvalues
    double precision, dimension(:,:), intent(in) :: c
    double precision, dimension(nbasis, nbasis), intent(out) :: w
    
    integer :: mu, nu, i

    w = 0.0d0
    
    do mu = 1, nbasis
    
        do nu = 1, nbasis
    
            do i = 1, nocc
    
                w(mu, nu) = w(mu, nu) + 2.0d0 * eigenvalues(i) * c(mu, i) * c(nu, i)
    
            end do
        end do
    end do

END SUBROUTINE compute_w_matrix


SUBROUTINE compute_overlap_grad_contrib(nbasis, natoms, n_primitives, exponents, coefficients, atom_coords, atom_map, Q, grad_S)

    ! compute the overlap (lagrangian) contribution to the gradient
    ! using the energy-weighted density matrix.

    integer, intent(in) :: nbasis, natoms
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:, :), intent(in) :: exponents, coefficients
    double precision, dimension(:, :), intent(in) :: atom_coords
    integer, dimension(:), intent(in) :: atom_map
    double precision, dimension(:, :), intent(in) :: q
    double precision, dimension(natoms, 3), intent(out) :: grad_s

    integer :: mu, nu, k, l, at_mu, at_nu
    double precision :: alpha, beta, zeta, xi, dist_sq, sab
    double precision :: ra(3), rb(3), rp(3), ds(3), pref

    double precision, parameter :: pi = acos(-1.0d0)

    ! initialize the overlap contribution to the gradient

    grad_s = 0.0d0

    ! loop over unique basis-function pairs and use translational symmetry

    do mu = 1, nbasis
        
        at_mu = atom_map(mu)
        ra = atom_coords(at_mu, :)

        do nu = 1, mu
        
            at_nu = atom_map(nu)
            rb = atom_coords(at_nu, :)

            if (at_mu == at_nu) cycle

            dist_sq = sum((ra - rb)**2)

            ! include the symmetric partner when mu /= nu

            pref = q(mu,nu)

            if (mu /= nu) pref = pref + q(nu,mu)

            ! sum over primitive gaussian pairs belonging to basis functions mu and nu

            do k = 1, n_primitives(mu)
                
                alpha = exponents(mu,k)

                do l = 1, n_primitives(nu)
                   
                    beta = exponents(nu,l)
                    zeta = alpha + beta
                    xi   = alpha * beta / zeta

                    sab = coefficients(mu,k) * coefficients(nu,l) * exp(-xi * dist_sq) * (pi / zeta)**1.5d0

                    ! gaussian product center and derivative of the primitive overlap term

                    rp = (alpha * ra + beta * rb) / zeta
                    ds = 2.0d0 * alpha * (rp - ra) * sab

                    grad_s(at_mu,:) = grad_s(at_mu,:) - pref * ds(:)
                    grad_s(at_nu,:) = grad_s(at_nu,:) + pref * ds(:)
                
                end do
            end do
        end do
    end do
END SUBROUTINE compute_overlap_grad_contrib

SUBROUTINE compute_kinetic_grad_contrib(nbasis, natoms, n_primitives, exponents, coefficients, atom_coords, atom_map, P, grad_T)
    
    ! compute the kinetic-energy contribution to the gradient.
    
    integer, intent(in) :: nbasis, natoms
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:, :), intent(in) :: exponents, coefficients
    double precision, dimension(:, :), intent(in) :: atom_coords
    integer, dimension(:), intent(in) :: atom_map
    double precision, dimension(:, :), intent(in) :: p
    double precision, dimension(natoms, 3), intent(out) :: grad_t

    integer :: mu, nu, k, l, at_mu, at_nu
    double precision :: alpha, beta, zeta, xi, dist_sq, sab, tab
    double precision :: ra(3), rb(3), rp(3), dt(3), pref

    double precision, parameter :: pi = acos(-1.0d0)

    ! initialize the kinetic contribution to the gradient

    grad_t = 0.0d0

    ! loop over unique basis-function pairs and use translational symmetry

    do mu = 1, nbasis
        
        at_mu = atom_map(mu)
        ra = atom_coords(at_mu, :)

        do nu = 1, mu
        
            at_nu = atom_map(nu)
            rb = atom_coords(at_nu, :)

            if (at_mu == at_nu) cycle

            dist_sq = sum((ra - rb)**2)

            ! include the symmetric partner when mu /= nu

            pref = p(mu,nu)

            if (mu /= nu) pref = pref + p(nu,mu)

            ! sum over primitive gaussian pairs belonging to basis functions mu and nu

            do k = 1, n_primitives(mu)
                
                alpha = exponents(mu,k)

                do l = 1, n_primitives(nu)
                
                    beta = exponents(nu,l)
                    zeta = alpha + beta
                    xi   = alpha * beta / zeta

                    sab = coefficients(mu,k) * coefficients(nu,l) * exp(-xi * dist_sq) * (pi / zeta)**1.5d0

                    ! primitive kinetic-energy term and its derivative

                    tab = xi * (3.0d0 - 2.0d0 * xi * dist_sq) * sab
                    rp  = (alpha * ra + beta * rb) / zeta

                    dt = 2.0d0 * alpha * (rp - ra) * (tab + 2.0d0 * xi * sab)

                    grad_t(at_mu,:) = grad_t(at_mu,:) + pref * dt(:)
                    grad_t(at_nu,:) = grad_t(at_nu,:) - pref * dt(:)
                
                end do
            end do
        end do
    end do

END SUBROUTINE compute_kinetic_grad_contrib

SUBROUTINE compute_nuclear_attraction_grad_contrib(nbasis, natoms, n_primitives, exponents, coefficients, atom_coords, atom_map, Z, P, grad_V)
    
    ! compute the electron-nucleus attraction contribution to the gradient.
    ! contributions from the basis-function centers and from the nuclei are
    ! accumulated explicitly.    

    integer, intent(in) :: nbasis, natoms
    integer, dimension(:), intent(in) :: n_primitives
    integer, dimension(:), intent(in) :: atom_map
    integer, dimension(:), intent(in) :: z
    double precision, dimension(:, :), intent(in) :: exponents, coefficients
    double precision, dimension(:, :), intent(in) :: atom_coords
    double precision, dimension(:, :), intent(in) :: p
    double precision, dimension(natoms, 3), intent(out) :: grad_v

    integer :: mu, nu, k, l, i
    integer :: at_mu, at_nu
    double precision :: alpha, beta, zeta, xi
    double precision :: rab2, rip2
    double precision :: f0, f1, pref
    double precision :: sab, v0, v1
    double precision :: ra(3), rb(3), rivec(3), rp(3)
    double precision :: dvdmu(3), dvdnu(3), dvdi(3)

    double precision, parameter :: pi = acos(-1.0d0)

    ! initialize the nuclear-attraction contribution to the gradient

    grad_v = 0.0d0

    ! loop over unique basis-function pairs

    do mu = 1, nbasis
        
        at_mu = atom_map(mu)
        ra = atom_coords(at_mu, :)

        do nu = 1, mu
        
            at_nu = atom_map(nu)
            rb = atom_coords(at_nu, :)

            ! include the symmetric partner when mu /= nu

            pref = p(mu, nu)
        
            if (mu /= nu) pref = pref + p(nu, mu)

            rab2 = sum((ra - rb)**2)

            ! sum over all nuclei contributing to the attraction term

            do i = 1, natoms
                
                rivec = atom_coords(i, :)

                ! sum over primitive gaussian pairs

                do k = 1, n_primitives(mu)
                
                    alpha = exponents(mu, k)

                    do l = 1, n_primitives(nu)
                
                        beta = exponents(nu, l)
                        zeta = alpha + beta
                        xi   = alpha * beta / zeta
                        rp   = (alpha * ra + beta * rb) / zeta

                        sab = coefficients(mu, k) * coefficients(nu, l) * exp(-xi * rab2) * (pi / zeta)**1.5d0

                        rip2 = sum((rivec - rp)**2)

                        ! boys-function values for the primitive attraction integrals

                        call boys_func(zeta * rip2, f0, f1)

                        v0 = -2.0d0 * dble(z(i)) * sqrt(zeta / pi) * sab * f0
                        v1 = -2.0d0 * dble(z(i)) * sqrt(zeta / pi) * sab * f1

                        ! derivatives with respect to the two basis-function centers and the nucleus

                        dvdmu = 2.0d0 * alpha * ( (rp - ra) * v0 - (rp - rivec) * v1 )
                        dvdnu = 2.0d0 * beta  * ( (rp - rb) * v0 - (rp - rivec) * v1 )
                        dvdi  = 2.0d0 * zeta  * (rp - rivec) * v1

                        grad_v(at_mu, :) = grad_v(at_mu, :) + pref * dvdmu(:)
                        grad_v(at_nu, :) = grad_v(at_nu, :) + pref * dvdnu(:)
                        grad_v(i,     :) = grad_v(i,     :) + pref * dvdi(:)
                    end do
                end do
            end do
        end do
    end do

end SUBROUTINE compute_nuclear_attraction_grad_contrib

SUBROUTINE compute_eri_grad_contrib(nbasis, natoms, n_primitives, exponents, coefficients, atom_coords, atom_map, P, grad_eri)
    
    ! Compute the two-electron contribution to the gradient from the
    ! derivatives of the electron-repulsion integrals.

    integer, intent(in) :: nbasis, natoms
    integer, dimension(:), intent(in) :: n_primitives
    double precision, dimension(:, :), intent(in) :: exponents, coefficients
    double precision, dimension(:, :), intent(in) :: atom_coords
    integer, dimension(:), intent(in) :: atom_map
    double precision, dimension(:, :), intent(in) :: p
    double precision, dimension(natoms, 3), intent(out) :: grad_eri

    integer :: mu, nu, lambda, sigma, k, l, m, n
    integer :: at_mu, at_nu, at_lambda, at_sigma
    double precision :: alpha, beta, gamma, delta
    double precision :: zeta, zetap, xi, xip, rho
    double precision :: rab2, rcd2, rpq2
    double precision :: f0, f1, eri0, eri1, pref
    double precision :: ra(3), rb(3), rc(3), rd(3)
    double precision :: rp(3), rq(3), rw(3)
    double precision :: dmu(3), dnu(3), dlambda(3), dsigma(3)
    double precision :: kab, kcd

    double precision, parameter :: pi = acos(-1.0d0)

    ! initialize the two-electron contribution to the gradient

    grad_eri = 0.0d0

    ! loop over all basis-function quartets

    do mu = 1, nbasis
        
        at_mu = atom_map(mu)
        ra = atom_coords(at_mu, :)

        do nu = 1, nbasis
        
            at_nu = atom_map(nu)
            rb = atom_coords(at_nu, :)

            rab2 = sum((ra - rb)**2)

            do lambda = 1, nbasis
           
                at_lambda = atom_map(lambda)
                rc = atom_coords(at_lambda, :)

                do sigma = 1, nbasis
           
                    at_sigma = atom_map(sigma)
                    rd = atom_coords(at_sigma, :)

                    ! combined coulomb and exchange prefactor for this quartet

                    pref = 0.5d0 * p(mu,nu) * p(lambda,sigma) - 0.25d0 * p(mu,lambda) * p(nu,sigma)

                    ! skip quartets whose density prefactor is negligibly small

                    if (abs(pref) < 1.0d-14) cycle

                    rcd2 = sum((rc - rd)**2)

                    ! sum over primitive gaussian quartets

                    do k = 1, n_primitives(mu)
                        alpha = exponents(mu,k)

                        do l = 1, n_primitives(nu)

                            beta = exponents(nu,l)
                            zeta = alpha + beta
                            xi = alpha * beta / zeta
                            rp = (alpha * ra + beta * rb) / zeta

                            ! prefactor for the first gaussian pair (mu,nu)

                            kab = coefficients(mu,k) * coefficients(nu,l) * sqrt(2.0d0) * pi**1.25d0 / zeta * exp(-xi * rab2)

                            do m = 1, n_primitives(lambda)

                                gamma = exponents(lambda,m)

                                do n = 1, n_primitives(sigma)

                                    delta = exponents(sigma,n)
                                    zetap = gamma + delta
                                    xip = gamma * delta / zetap
                                    rq = (gamma * rc + delta * rd) / zetap

                                    ! prefactor for the second gaussian pair (lambda,sigma)

                                    kcd = coefficients(lambda,m) * coefficients(sigma,n) * sqrt(2.0d0) * pi**1.25d0 / zetap * exp(-xip * rcd2)

                                    rho = zeta * zetap / (zeta + zetap)
                                    rw = (zeta * rp + zetap * rq) / (zeta + zetap)
                                    rpq2 = sum((rp - rq)**2)

                                    ! boys-function values for the primitive eri derivative

                                    call boys_func(rho * rpq2, f0, f1)

                                    eri0 = kab * kcd * f0 / sqrt(zeta + zetap)
                                    eri1 = kab * kcd * f1 / sqrt(zeta + zetap)  

                                    ! derivatives with respect to the four centers;
                                    ! the last one is obtained from translational invariance

                                    dmu = 2.0d0 * alpha * ((rp - ra) * eri0 + (rw - rp) * eri1)
                                    dnu = 2.0d0 * beta * ((rp - rb) * eri0 + (rw - rp) * eri1)
                                    dlambda = 2.0d0 * gamma * ((rq - rc) * eri0 + (rw - rq) * eri1)
                                    dsigma = -(dmu + dnu + dlambda)

                                    grad_eri(at_mu,:) = grad_eri(at_mu,:) + pref * dmu(:)
                                    grad_eri(at_nu,:) = grad_eri(at_nu,:) + pref * dnu(:)
                                    grad_eri(at_lambda,:) = grad_eri(at_lambda,:) + pref * dlambda(:)
                                    grad_eri(at_sigma,:) = grad_eri(at_sigma,:) + pref * dsigma(:)
                                
                                end do
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end do
end subroutine compute_eri_grad_contrib

END MODULE hf_gradients