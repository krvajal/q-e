module kli
use kinds, ONLY : dp
use ld1_parameters, only:  nwfx ! max number of wavefunctions
use ld1inc, only: psi,& ! all electron wavefuntions
                  rho, &   ! all electron density rho(:,1)(up), rho(:,2) (down)
                  num_wave_functions => nwf, &
                  get_spin => isw , &
                  grid, & ! radial grid
                  orbital_angular_momentum =>  ll, &
                  shell_occupancy => oc, &
                  sl3

use radial_grids, only  : ndmx, hartree

implicit none
private 
public compute_kli_potential

real(dp), allocatable :: mat_m_up(:,:),mat_m_down(:,:)
real(dp) :: ux_kli(ndmx,nwfx) ! ndmx x nwf 
real(dp) :: potential_s(nwfx) ! ndmx x nwf
real(dp) :: average_ux_kli(nwfx), average_kli_potential(nwfx) 
integer :: num_up_wavefunctions, num_down_wavefunctions
real(dp) :: mat_m(nwfx, nwfx)
! linear problem variables A x = y (with x = V^{KLI}_X)
real(dp) :: A(nwfx, nwfx), y(nwfx), AA(nwfx, nwfx)

contains

    subroutine init_module()
        implicit none
        integer :: i
        num_up_wavefunctions = 0
        num_down_wavefunctions = 0
      
       
    end subroutine init_module



    subroutine compute_mat_m(N,grid_size, rho)
        integer,intent(in) :: N
        
        real(dp),intent(in) :: rho(ndmx)
        
        real(dp) :: int_0_inf_dr
        
        ! aux variables
        integer :: i, j, k
        integer :: grid_size 
        real (dp) :: func(ndmx) , fact1,fact2
        
        print *,"N", N
        print *, "shape", grid_size
        print *, ndmx
        do i = 1, N
            do j = 1, i
                ! psi (i) and psi(j) must have the same spin
                print *, "(i,j)" , (/ i,j /)
                func = abs(  psi(:,1,i))**2 * abs(psi(:,1,j))**2
                do k = 1, ndmx
                    if (abs(rho(k)) > tiny(1.0_dp)) then
                        func(k)  = func(k) / rho(k)
                    else if(abs(func(k)) > tiny(1.0_dp)) then
                        print *, "Density to small for k", k
                        stop 
                    end if
                enddo
                call print_vec(grid_size, func)
                ! call print_vec(size(func),func)
                fact1 = shell_occupancy(i)
                fact2 = shell_occupancy(j)
                mat_m(i,j) = fact1 * fact2 * int_0_inf_dr(func, &
                                        grid, &
                                        grid_size, &
                                        0) ! this is the asyntotic behavior for r -> 0 of the integrand  
            enddo
        enddo
        print *, mat_m(1,1:num_wave_functions)
        print *, mat_m(2,1: num_wave_functions)
        print *, mat_m(3,1:num_wave_functions)
    
        !M_ij = \int |\phi_i|^2 |\phi_j|^2 / rho


    end subroutine compute_mat_m

    subroutine print_vec(N, v)
        integer:: N,i
        real(dp) :: v(:)
        do i  = 1, N
            print *, "v[",i,"]=", v(i)
        enddo
    end subroutine print_vec


    subroutine compute_ux_kli(grid_size)
        ! compute the value of ux_kli given by 
        ! equation (118) of the GKK paper
        ! this is the so called orbital dependent potential
        integer, intent(in) :: grid_size

        integer :: i, j
        integer :: spin_i, spin_j
        integer :: l1, l2, lambda
        integer :: lmin, lmax
        real(dp) :: pseudodensity(ndmx)
        real(dp) :: pseudopot(ndmx)
        real(dp) :: l3_symbol_squared
        
        print *, "num_wave_functions", num_wave_functions
        
        do i = 1, num_wave_functions
            spin_i = get_spin(i)
            
            l1 =  orbital_angular_momentum(i)
            print *, "info_i", i, spin_i, l1
            ux_kli(:, i) = 0.0
            do j = 1, num_wave_functions
                spin_j = get_spin(j)
                
                l2 =  orbital_angular_momentum(j)
                print *, "info_j", j, spin_j, l2
                if ( spin_i == spin_j) then
                    ! only wavefunctions with the same spin
                    lmin =  abs(l1 - l2)
                    lmax =  l1 + l2
                    print *, "lmin=", lmin
                    print *, "lmax=", lmax, l1 + l2
                    print *,"------------------"
                    do lambda = lmin, lmax
                        l3_symbol_squared = sl3(l1,l2,lambda) * 0.5_dp
                        pseudodensity = psi(:,1,i) * psi(:,1,j)
                        print  *,"density",pseudodensity(1:10)
                        print *,"l3_symbol_squared", l3_symbol_squared
                        ! compute the radial part 
                        call hartree(lambda, l1 + l2 + 2, grid_size, grid, pseudodensity, pseudopot)
                        
                        ux_kli(:,i) = ux_kli(:, i) +  l3_symbol_squared * pseudopot * shell_occupancy(i)
                    enddo
                endif
            enddo
            call print_vec(grid_size, ux_kli(:,i))

            ! stop
        enddo

    end subroutine compute_ux_kli




    subroutine  compute_potential_s(grid_size)
        integer, intent(in) :: grid_size

        integer :: i,j
        real(dp) :: work(ndmx)
        real(dp) :: fact
        real(dp) :: int_0_inf_dr
        work = 0.0_dp
        do j = 1, num_wave_functions
            
            do i = 1, num_wave_functions
                fact = shell_occupancy(i)
                work = work + abs(psi(:, 1,i))**2 * ( ux_kli(:,i)  + ux_kli(:,i)) * 0.5_dp  * fact
            enddo
            work = work * ( shell_occupancy(j) * abs( psi(:, 1, j ))**2 )
            potential_s(j) = int_0_inf_dr(work, grid, grid_size, 0)
            print *,j, potential_s(j)
        enddo
    end subroutine compute_potential_s

    subroutine compute_average_ux_kli(grid_size)
        ! computes  <i| u_xi| i> 
        integer, intent(in) :: grid_size
        integer :: i, fact
        real(dp) :: work(ndmx)
        real(dp) :: int_0_inf_dr
        work = 0.0_dp
        do i  = 1, num_wave_functions
            work =  psi(:,1,i) * ux_kli(:,i) * psi(:,1,i)
            average_ux_kli(i) = shell_occupancy(i) * int_0_inf_dr(work, grid, grid_size, 0)
        enddo
        print *,average_ux_kli(1: num_wave_functions)
    end subroutine compute_average_ux_kli


    subroutine compute_radial_part(l1, l2, k, grid_size, Rnl1, Rnl2, solution)
        integer,intent(in) :: l1,l2, k, grid_size
        real(dp),intent(in) :: Rnl1(ndmx), Rnl2(ndmx)
        real(dp),intent(out) :: solution(ndmx)
        real(dp) :: f(ndmx), zeta(ndmx)
        real(dp) :: dx, fact
        integer :: i
        
        dx = grid%dx

        print *, "dx", dx
        
        f =  grid%mesh * Rnl1 * Rnl2
        print *, f(1:10)
        print *, "l1,l2,k", (/l1, l2, k /)
        
        zeta(1) = f(1)/(l1 + l2 + k + 3)
        zeta(2) = f(2)/(l1 + l2 + k + 3)
        
        !integrate aux z function outwards
        do i = 3, grid_size
            zeta(i) = zeta(i-2)*exp(- 2 *dx * k) + (dx/3)*(f(i) + 4* f(i-1) * exp(-dx*k) + exp(-2*dx*k)*f(i-2))
        enddo
        print *,"zeta", zeta(1:10)

        solution(grid_size) =  zeta(grid_size)
        solution(grid_size - 1) = zeta(grid_size - 1) 
        do i = grid_size - 2, 1
            fact = zeta(i) + 4 * exp(-dx * (k + 1)) * zeta(i + 1) + exp(-2 * dx *(k + 1)) * zeta(i + 2)
            solution(i) = solution(i + 2) * exp(-2*dx*(k+1)) + (2 * k + 1)*(dx/3)*fact
        enddo
        

    end subroutine compute_radial_part


    subroutine solve_linear_problem(num_wave_functions, average_kli_potential)
        integer,intent(in) :: num_wave_functions
        integer, intent(out) :: average_kli_potential(nwfx)
        integer :: i, info
        integer :: ipivot(num_wave_functions)
        integer :: lda = nwfx

        print *, "Solving linear problem"
        y = 0 
        A = 0
        A = - mat_m
        ! A = A + I
        do i = 1, num_wave_functions
            A(i,i) = A(i,i) + 1.0_dp
            y(i) = potential_s(i) - dot_product(mat_m(1:num_wave_functions,i), average_ux_kli(1:num_wave_functions))
            print *,"y", y(i)
        enddo
        print *, A(1,1:num_wave_functions)
        print *, A(2,1:num_wave_functions)
        print *, A(3,1:num_wave_functions)

       ! dim(A) = num_wave_functions * num_wave_functions
       ! solve real matrix Ax = b using blas with double precision
       AA = A ! store original
       call DGETRF(num_wave_functions, num_wave_functions, A,lda,ipivot, info)
       call DGETRS('N', num_wave_functions, 1, A, lda, ipivot, y, num_wave_functions, info)
       print *, "info", info
       print *, y(1:num_wave_functions)
       average_kli_potential = y  ! save the solution
       stop
    end subroutine solve_linear_problem



    subroutine compute_kli_potential(grid_size,  exchange_potential)
        
        integer,  intent(in) ::  grid_size  ! number of grid points
        real(dp), intent(out) :: exchange_potential(ndmx, 2) 
      

        ! access the wavefunctions
        print *, "compute kli potential"
        print *, num_wave_functions
        
       
        call compute_mat_m(num_wave_functions, grid_size,rho)
        call compute_ux_kli(grid_size)
        print *, "computing s potential"
        call compute_potential_s(grid_size)
        print *, "computing averagSe orbital dependent potential"
        call compute_average_ux_kli(grid_size)
        call solve_linear_problem(num_wave_functions, average_kli_potential)


        do i = 1, num_wave_functions
        enddo
        


    end subroutine compute_kli_potential
end module kli