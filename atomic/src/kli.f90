module kli
use kinds, ONLY : dp
use ld1_parameters, only:  nwfx ! max number of wavefunctions
use ld1inc, only: psi,& ! all electron wavefuntions
                  rho, &   ! all electron density rho(:,1)(up), rho(:,2) (down)
                  num_wave_functions => nwf, &
                  get_spin => isw , &
                  grid, & ! radial grid
                  orbital_angular_momentum =>  ll, &
                  oc, &
                  nspin,&
                  sl3

use radial_grids, only  : ndmx, hartree

implicit none
private 
public compute_kli_potential

real(dp), allocatable :: mat_m_up(:,:),mat_m_down(:,:)
real(dp) :: ux_kli(ndmx,nwfx) ! ndmx x nwf 
real(dp) :: v_x_hf(ndmx,nwfx) ! ndmx x nwf 
real(dp) :: potential_s(nwfx) ! ndmx x nwf
real(dp) :: average_ux_kli(nwfx), average_kli_potential(nwfx) 
integer :: num_up_wavefunctions, num_down_wavefunctions
real(dp) :: mat_m(nwfx, nwfx)
! linear problem variables A x = y (with x = V^{KLI}_X)
real(dp) :: A(nwfx, nwfx), y(nwfx), AA(nwfx, nwfx)

contains

    function shell_occupancy(i) result(occup)
        integer :: i
        real(dp) :: occup
        occup = oc(i)  !* 0.5_dp * nspin
    end function shell_occupancy


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
                
                ! call print_vec(size(func),func)
                fact1 = shell_occupancy(i) 
                fact2 = shell_occupancy(j) 
                mat_m(i,j) = fact1 * fact2 * int_0_inf_dr(func, &
                                        grid, &
                                        grid_size, &
                                        2*orbital_angular_momentum(i) + 2) ! this is the asyntotic behavior for r -> 0 of the integrand  
            enddo
        enddo
        if( mat_m(1,1) /= mat_m(1,1)) then
            print *, "We got a problem!!!!"
            print *, psi(1:10,1,1)
            stop
        endif

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
        
        v_x_hf  = 0.0_dp
        ux_kli = 0.0_dp
        do i = 1, num_wave_functions
            call dvex(i,v_x_hf(:, i))
            ! spin_i = get_spin(i)
            
            ! l1 =  orbital_angular_momentum(i)
            ! print *, "info_i", i, spin_i, l1
            ! ux_kli(:, i) = 0.0
            ! do j = 1, num_wave_functions
            !     spin_j = get_spin(j)
                
            !     l2 =  orbital_angular_momentum(j)
            !     print *, "info_j", j, spin_j, l2
            !     if ( spin_i == spin_j) then
            !         ! only wavefunctions with the same spin
            !         lmin =  abs(l1 - l2)
            !         lmax =  l1 + l2
            !         print *, "lmin=", lmin
            !         print *, "lmax=", lmax, l1 + l2
            !         print *,"------------------"
            !         do lambda = lmin, lmax
            !             l3_symbol_squared = sl3(l1,l2,lambda) * 0.5_dp
            !             pseudodensity = psi(:,1,i) * psi(:,1,j)
            !             print  *,"density",pseudodensity(1:10)
            !             print *,"l3_symbol_squared", l3_symbol_squared
            !             ! compute the radial part 
            !             call hartree(lambda, l1 + l2 + 2, grid_size, grid, pseudodensity, pseudopot)
                        
            !             ux_kli(:,i) = ux_kli(:, i) +  l3_symbol_squared * pseudopot * shell_occupancy(i)
            !         enddo
            !     endif
                
            ! enddo
            
            ux_kli(:,i) = v_x_hf(:,i) / (psi(:,1,i)* shell_occupancy(i))
            
            ! call print_vec(grid_size, ux_kli(:,i))
            
        enddo
        

    end subroutine compute_ux_kli


    subroutine  compute_potential_s(grid_size)
        integer, intent(in) :: grid_size

        integer :: i,j
        real(dp) :: work(ndmx)
        real(dp) :: fact
        real(dp) :: int_0_inf_dr
        integer :: nst
        
        
        do j = 1, num_wave_functions
            work = 0.0_dp    
            do i = 1, num_wave_functions
                fact = shell_occupancy(i) 
                work = work + abs(psi(:, 1,i)) * ( v_x_hf(:,i))  * fact
                ! work = work + psi(:, 1,i) * ( ux_kli(:,i)  + ux_kli(:,i)) * 0.5_dp  * fact
            enddo
            nst = 2 * orbital_angular_momentum(j) + 2
            
            work = work * ( shell_occupancy(j) * abs( psi(:, 1, j ))**2 ) / rho(:,1)
            ! work = work * ( shell_occupancy(j) * psi(:, 1, j )) / rho(:,1)
            potential_s(j) = int_0_inf_dr(work, grid, grid_size, nst) 
            print *,j, potential_s(j)
        enddo
    end subroutine compute_potential_s

    subroutine compute_average_ux_kli(grid_size)
        ! computes  <i| u_xi| i> 
        integer, intent(in) :: grid_size
        integer :: i, fact
        real(dp) :: work(ndmx)
        real(dp) :: int_0_inf_dr
        integer :: nst
        
        work = 0.0_dp
        do i  = 1, num_wave_functions
            work =  psi(:,1,i)**2 * ux_kli(:,i)
            fact = shell_occupancy(i)
            nst = 2 * orbital_angular_momentum(i) + 2
            average_ux_kli(i) = fact * int_0_inf_dr(work, grid, grid_size, nst) 
        enddo
        
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
        

        solution(grid_size) =  zeta(grid_size)
        solution(grid_size - 1) = zeta(grid_size - 1) 
        do i = grid_size - 2, 1
            fact = zeta(i) + 4 * exp(-dx * (k + 1)) * zeta(i + 1) + exp(-2 * dx *(k + 1)) * zeta(i + 2)
            solution(i) = solution(i + 2) * exp(-2*dx*(k+1)) + (2 * k + 1)*(dx/3)*fact
        enddo
        

    end subroutine compute_radial_part


    subroutine solve_linear_problem(num_wave_functions, average_kli_potential)

        integer,intent(in) :: num_wave_functions
        real(dp), intent(out) :: average_kli_potential(nwfx)

        integer :: i, info, N
        integer :: ipivot(num_wave_functions)
        integer :: lda = nwfx
        
        print *,"=============="
        print *,"linear system"

        N = num_wave_functions -1
        y = 0 
        A = 0
        A = - mat_m
        ! A = A + I
        do i = 1, N
            A(i,i) = A(i,i) + 1.0_dp
            y(i) = potential_s(i) - dot_product(mat_m(1:N,i), average_ux_kli(1:N))
            print *,"y", y(i)
        enddo

        print *,"A", A(1,1:N)
        print *,"M", mat_m(1,1)
        print *, "V^{S}", potential_s(1)
        print *, "U_X", average_ux_kli(1)

        ! print *, A(2,1:num_wave_functions)
        ! print *, A(3,1:num_wave_functions)
        if (N > 0) then
            ! dim(A) = num_wave_functions * num_wave_functions
            ! solve real matrix Ax = b using blas with double precision
            AA = A ! store original
            call DGETRF(N, N, A,lda,ipivot, info)
            call DGETRS('N', N, 1, A, lda, ipivot, y, N, info)
            !    print *, "info", info
            !    print *, y(1:num_wave_functions)
            average_kli_potential = y  ! save the solution
        endif

        average_kli_potential(N + 1) = average_ux_kli(N + 1) !last term
        print *,"res", average_kli_potential(1)
    !    stop
    end subroutine solve_linear_problem

    subroutine compute_kli_potential(grid_size,  exchange_potential)
        
        integer,  intent(in) ::  grid_size  ! number of grid points
        real(dp), intent(out) :: exchange_potential(ndmx, 2) 
        real(dp) :: work(ndmx), fact(ndmx), fact1
        real(dp) :: slater_potential(ndmx)
        integer :: i, j

        do i =  1, num_wave_functions
            do j = 1, grid_size
                if(psi(j,1,i) /= psi(j,1,i)) then
                    print *, "We got a problem!!!"
                    print *, "Invalid wavefunction passed"
                    stop    
                endif
            enddo
        enddo
        

        ! access the wavefunctions
        print *, "compute kli potential"
        print *, num_wave_functions
        
        mat_m = 0.0_dp
        call compute_mat_m(num_wave_functions, grid_size,rho)
        print *,"======="
        print *, "Matrix M"
        do i = 1, num_wave_functions
            print *,mat_m(i,1:num_wave_functions)
        enddo
        
        call compute_ux_kli(grid_size)
        print *,"KLI UX"
        print *, ux_kli(1:10,1)
        print *,"!!!!!!!!!"
        ! stop
        print *, "computing s potential"
        call compute_potential_s(grid_size)
        print *, "computing average orbital dependent potential"
        call compute_average_ux_kli(grid_size)
        call solve_linear_problem(num_wave_functions, average_kli_potential)

        ! print *, "here", average_ux_kli(1) - average_kli_potential(1)
        ! print *, average_ux_kli(1)
        ! stop
        slater_potential = 0

        do i = 1, num_wave_functions
            slater_potential = slater_potential +  psi(:,i,1) * v_x_hf(:,i) * shell_occupancy(i) 
        enddo

        slater_potential =  slater_potential / rho(:,1)

        work = 0
        do i = 1, num_wave_functions - 1
            fact =  (average_kli_potential(i) - average_ux_kli(i))
            ! work =  work + abs(psi(:,1,i))**2 * fact * shell_occupancy(i)
             work =  work + psi(:,1,i)**2 * fact * shell_occupancy(i)/2
            ! work =  work + psi(:,1,i)  * fact
        enddo

        work =  work / rho(:,1)
        print *, "average kli", average_kli_potential(1:num_wave_functions)
        print *, "average kli", average_ux_kli(1:num_wave_functions)

        exchange_potential(:,1) = slater_potential + work
        
        
    end subroutine compute_kli_potential

    subroutine savetxtv2(filename,x,y)
    ! Saves a x,y array into a textfile.
    !
    ! Arguments
    ! ---------
    !
    character(len=*), intent(in) :: filename  ! File to save the array to
    real(dp), intent(in) :: x(:)           ! The x data points
    real(dp), intent(in) :: y(:)           ! The y data points
    !
    ! Example
    ! -------
    !
    ! real(dp) :: data(3, 2)
    ! call savetxt("log.txt", data)

    integer :: s, i
    open(newunit=s, file=filename, status="replace")
    do i = 1, size(x, 1)
        write(s, *) x(i), y(i)
    end do
    close(s)

    end subroutine savetxtv2
end module kli