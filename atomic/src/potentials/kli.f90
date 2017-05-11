module kli
use kinds, ONLY : dp
use utils
use ld1inc, only: psi,& ! all electron wavefuntions
                  rho, &   ! all electron density rho(:,1)(up), rho(:,2) (down)
                  num_wave_functions => nwf, &
                  get_spin => isw , &
                  grid, & ! radial grid
                  orbital_angular_momentum =>  ll, &
                  shell_occupancy => oc

use radia_grids, only  : ndmx

implicit none
private 
public compute_kli_potential

real(dp), allocatable :: mat_m_up(:,:),mat_m_down(:,:)
real(dp) :: ux_kli(mdmx, num_wave_functions)

integer :: num_up_wavefunctions, num_down_wavefunctions
integer :: wavefunc_index(2,num_wave_functions)
contains



    subroutine init_module()
        implicit none
        integer :: i
        num_up_wavefunctions = 0
        num_down_wavefunctions = 0
        do i=1,num_wave_functions
            if(get_spin(i) == 1) then
             num_up_wavefunctions = num_up_wavefunctions + 1
             wavefunc_index(1,num_up_wavefunctions) = i
            enddo   
            if(get_spin(i) == 2) then
                num_down_wavefunctions = num_down_wavefunctions + 1
                wavefunc_index(2,num_down_wavefunctions) = i
            enddo
        enddo
        call assert( (num_up_wavefunctions + num_down_wavefunctions) == num_wave_functions)

        allocate(mat_m_up(num_up_wavefunctions,num_up_wavefunctions))
        allocate(mat_m_down(num_down_wavefunctions, num_down_wavefunctions)))

        allocate(matrix_k_up(num_wave_functions, num_wave_functions, ndmx))

    end subroutine init_module()


    subroutine compute_mat_m(N, wavefuntions, rho,  mat):
        integer,intent(in) :: N
        real,intent(in) :: wavefuntions(:,:)
        real,intent(in) :: rho(:)
        real, intent(out) :: mat(N,N)

        ! aux variables
        integer :: i_idx, j_idx
        integer :: grid_size 
        grid_size = size(rho)

        do i =1, N
            do j =1, i
                M(i,j) = integrate(abs( psi(:,i)) * abs(psi(:,j))/rho)
            enddo
        enddo

        !M_ij = \int |\phi_i|^2 |\phi_j|^2 / rho


    end subroutine compute_mat_m

    subroutine compute_ux_kli()
        ! compute the value of ux_kli given by 
        ! equation (118) of the GKK paper
        integer :: i, j
        integer :: spin_i, spin_j
        integer :: l1, l2, lambda
        integer :: lmin, lmax
        real(dp) :: pseudodensity(ndmx)
        real(dp) :: pseudopot(ndmx)
        


        do i = 1, num_wave_functions
            spin_i = get_spin(i)
            l1 =  orbital_angular_momentum(i)
            ux_kli(:, i) = 0.0
            do j = 1, num_down_wavefunctions
                spin_j = get_spin(j)
                l2 =  orbital_angular_momentum(j)
                if ( spin_i == spin_j) then
                    ! only wavefunctions with the same spin
                    lmin =  abs(l1 - l2)
                    lmax =  l1 + l2
                    pseudodensity =  psi(i) * psi(j)
                    do lambda = lmin, lmax
                        3l_symbol_squared = sl3(l1,l2,lambda) * 0.5_dp
                        ! compute the radial part 
                        
                        ux_kli(:,i) = ux_kli(:, i) +  3l_symbol_squared * pseudopot * shell_occupancy(i)
                    enddo
                enddo
            enddo
        enddo


    end subroutine compute_ux_kli

    subroutine compute_radial_part(l1,l2, k, Rnl1, Rnl2, solution)
        integer,intent(in) :: l1,l2, k
        real(dp),intent(in) :: Rnl1(ndmx), Rnl2(ndmx)
        real(dp),intent(out) :: solution(ndmx)
        real(dp) :: f(ndmx), zeta(ndmx)
        real(dp) :: dx, fact
        integer :: i
        
        dx = grid%dx
        f =  grid%mesh * Rnl1 * Rnl2
        zeta(1) = f(1)/(l1 + l2 + k + 3)
        zeta(2) = f(2)/(l1 + l2 + k + 3)
        
        !integrate aux z function outwards
        do i = 3, ndmx
            zeta(i) = z(i-2)*exp(- 2 *dx * k) + (dx/3)*(f(i) + 4* f(i-1) * exp(-dx*k) + exp(-2*dx*k)*f(i-2))
        enddo

        solution(ndmx) =  zeta(ndmx)
        solution(ndmx-1) = zeta(ndmx - 1) 
        do i = ndmx - 2, 1
            fact = zeta(i) + 4 * exp(-dx * (k + 1)) * z(i + 1) + exp(-2 * dx *(k + 1)) * z(i + 2))
            solution(i) = solution(i + 2) * exp(-2*dx*(k+1)) + (2 * k + 1)*(dx/3)*fact
        enddo
        

    end subroutine compute_radial_part


    function integrate(N,func) result(retval)
        integer,intent(in) :: N
        real(dp),intent(in) :: func(N)
        real(dp),intent(out) :: retval
    end function integrate


    

    subroutine compute_kli_potential(ndmx,  exchange_potential)
        implicit none
        integer intent(in) ::  ndmx  ! number of grid points
        real(dp) :: exchange_potential(ndmx, 2) 

        ! access the wavefunctions


    end subroutine compute_kli_potential
end module kli