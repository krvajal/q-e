!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!


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
                  sl3, &
                  title

use radial_grids, only  : ndmx, hartree

implicit none
private 
public compute_kli_potential



real(dp) :: v_x_hf(ndmx,nwfx) !  V^{HF}_i =  \sum \phi_j(r) \int d^3r' \phi_i(r')\phi_j(r')/ | r - r'|
real(dp) :: ux_kli(ndmx,nwfx) !  u_{xi\sigma} = V^{HF}_i / \psi_i
real(dp) :: potential_s(nwfx)! < i| V^{AFA}_\sigma| i> 
real(dp) :: average_ux_kli(nwfx), &  ! < i| u_{xi\sigma}| i> 
            average_kli_potential(nwfx) ! < i| Vx^{KLI}| i> 
integer :: num_up_wavefunctions, num_down_wavefunctions
real(dp) :: mat_m(nwfx, nwfx, 2)  ! M_ij = \int |\psi_i|^2 |\psi_j|^2 / rho
integer :: num_wf(2) ! number of up and down wavefunctions
real(dp) :: oc2(2)
integer :: idx(nwfx,2) ! mapping between the vector of all wawefunctions and the values depending on the spin
! linear problem variables A x = y (with x = V^{KLI}_X)
real(dp) :: A(nwfx, nwfx), y(nwfx),  ysol(nwfx), AA(nwfx, nwfx)
real(dp) :: slater_potential(ndmx,2) ! V^{AFA}_\sigma \sum_i^{N_\sigma} |\phi_{i\sigma}|^2/\rho_sigma * u_{xi\sigma}
contains

    subroutine compute_num_wf(num_wave_functions, num_wf)
    ! determine the number of up and down wavefunctions
        integer , intent(in) :: num_wave_functions
        integer , intent(out) :: num_wf(2) 


        integer :: i, s ! aux variables
        integer :: nu
        num_wf = 0
       
        do i =  1, num_wave_functions
            s = get_spin(i)
            if (oc(i) > 0) then
                num_wf(s) = num_wf(s) + 1
                idx(num_wf(s) , s) = i !index of the jth wf with spin s in the global array
            endif
        enddo   
        if (nspin == 1) then
            if (num_wf(2) /= 0) stop "error"
        endif
    end subroutine compute_num_wf

    function shell_occupancy(i) result(occup)
        integer :: i
        real(dp) :: occup
        
        occup = oc(i)  !* (nspin / 2.0_dp)
    end function shell_occupancy


    subroutine init_module()
        implicit none
        integer :: i
        num_up_wavefunctions = 0
        num_down_wavefunctions = 0
      
       
    end subroutine init_module


    subroutine compute_mat_m(N,grid_size, rho)

        integer,intent(in) :: N
        
        real(dp),intent(in) :: rho(ndmx,2)
        
        real(dp) :: int_0_inf_dr
        
        ! aux variables
        integer :: i, j, k, ii
        integer :: grid_size 
        integer ::  s_i, s_j, idx_i, idx_j
        real (dp) :: func(ndmx) , fact1,fact2
        real (dp) :: retval
        integer :: s ! spin interation index
        integer :: nst ! leading behavior for r -> 0 
        real(dp) :: half
        !M_ij = \int |\phi_i|^2 |\phi_j|^2 / rho

        do s = 1,nspin
            do i = 1, num_wf(s)
                idx_i = idx(i,s)
                s_i = get_spin(idx_i) 
                ! s_i should be equal s
                
                do j = 1, num_wf(s)
                    ! psi (i) and psi(j) must have the same spin
                    idx_j = idx(j,s)
                    s_j =  get_spin(idx_j)
                    if( s_i == s_j) then
                        
                        func = abs(  psi(:,1,idx_i))**2 * abs(psi(:,1,idx_j))**2
                        do k = 1, ndmx
                            ! if (abs(rho(k,s)) > tiny(1.0_dp)) then
                              func(k)  = func(k) /( rho(k,s))
                            ! else if(abs(func(k)) > tiny(1.0_dp)) then
                                ! print *, "Density to small for k", k
                                ! stop 
                            ! end if
                        enddo
                        ! call print_vec(size(func),func)
                        fact1 = (2 * orbital_angular_momentum(idx_i)  + 1)
                        fact2 =  shell_occupancy(idx_j)
                        
                        nst =   2 * orbital_angular_momentum(idx_i) + 2
                        retval = fact1 * fact2 * int_0_inf_dr(func, &
                                            grid, &
                                            grid_size, &
                                            nst) ! this is the asyntotic behavior for r -> 0 of the integrand  

                        ! retval goes to the i,j entry of matrix with spin s
                         mat_m(i,j,s) = retval 
                        print *, i,j, mat_m(i,j,s)
#ifdef DEBUG
                         if (retval < tiny(1.0_dp)) then    
                            print *,"Too small"
                            print *,psi(1:10,1,idx_i)
                            print *,psi(1:10,1,idx_j)
                            print *, fact1
                            print *, fact2
                            stop
                         endif
#endif
                    else
                        
                        stop "We got into trouble here"
                    endif ! the elements does not have the same spin               
                enddo ! loop over j
            enddo ! lopp over i
            
#ifdef DEBUG      
            ! small check for non nan
            if( mat_m(1,1,s) /= mat_m(1,1,s)) then
                print *, "We got a problem!!!!"
                print *, psi(1:10,1,1)
                stop
            endif

#endif
            
        enddo ! loop over s
        
        

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
        
        
        v_x_hf  = 0.0_dp
        ux_kli = 0.0_dp
        do i = 1, num_wave_functions
            call dvex(i,v_x_hf(:, i))
            ! ux_kli(:,i) = v_x_hf(:,i) / (psi(:,1,i))
            ! call print_vec(grid_size, ux_kli(:,i))
        enddo
        if(v_x_hf(1,1) == 0) then 
            print *, psi(1:10,1,2)
            stop "here"
        endif

    end subroutine compute_ux_kli


    subroutine  compute_potential_s(grid_size)
        integer, intent(in) :: grid_size

        integer :: i,j,s
        real(dp) :: work(ndmx)
        real(dp) :: fact
        real(dp) :: int_0_inf_dr
        integer :: nst, spin_i, spin_j
        
        slater_potential = 0.0 ! initialize to zero
        do s = 1, nspin
            do i = 1, num_wf(s)
                fact =  shell_occupancy(idx(i,s))
                slater_potential(:,s) = slater_potential(:,s) +  psi(:, 1, idx(i,s)) *  v_x_hf(:, idx(i,s)) * fact 
            enddo
            slater_potential(:,s) = slater_potential(:,s)/ (rho(:,s))
        enddo

        do j = 1, num_wave_functions
            spin_j = get_spin(j)
            nst = 2 * orbital_angular_momentum(j) + 2
            fact = 2 * orbital_angular_momentum(j) + 1
            work = psi(:, 1 ,  j ) * psi(:, 1, j) * slater_potential(:,spin_j)
            ! work = work * ( shell_occupancy(j) * psi(:, 1, j )) / rho(:,1)
            potential_s(j) = int_0_inf_dr(work, grid, grid_size, nst)  * shell_occupancy(j)
#ifdef DEBUG
            print *, j, potential_s(j)
#endif

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
            work =  psi(:,1,i) *  v_x_hf(:,i) * shell_occupancy(i)
            ! fact = shell_occupancy(i)* (2 * orbital_angular_momentum(i) + 1)
            nst = 2 * orbital_angular_momentum(i) + 2
            average_ux_kli(i) =  int_0_inf_dr(work, grid, grid_size, nst) 

            print *, i , average_ux_kli(i), get_spin(i)
            
            if(average_ux_kli(i) /= average_ux_kli(i)) then 
                print *, "We got a problem" ! NaN value
                print *, "Invalid average kli"
                print *, ux_kli(:,i)
                stop
            endif
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

        integer :: i, info, N, s,l
        integer :: ipivot(num_wave_functions)
        integer :: lda = nwfx
        
        average_kli_potential = 0
        do s  = 1, nspin
            N = num_wf(s) - 1
            
            if (N > 0) then
                y = 0 
                A = 0
                A = - mat_m(:,:,s)
                ! A = A + I
                do i = 1, N 
                    if(shell_occupancy(idx(i,s)) == 0) then
                        print *, "error"
                        stop
                    endif
                    A(i,i) = A(i,i) + 1.0_dp 
                    y(i) = ( potential_s(idx(i,s)) -  average_ux_kli(idx(i,s))) 
                    ! print *,"y", i, y(i)
                enddo

                ! print *,"A", A(1,1:N)
                ! print *,"M", mat_m(1,1,s)
                ! print *, "V^{S}", potential_s(1)
                ! print *, "U_X", average_ux_kli(1)

            ! print *, A(2,1:num_wave_functions)
            ! print *, A(3,1:num_wave_functions)
        
                ! dim(A) = num_wave_functions * num_wave_functions
                ! solve real matrix Ax = b using blas with double precision
                AA = A ! store original
                call DGETRF(N, N, A,lda,ipivot, info)
                if(info /= 0) then
                    print *, "Failed to factorize matrix"
                    stop
                endif
                call DGETRS('N', N, 1, A, lda, ipivot, y, N, info)
                if (info /= 0) then
                    print *, "Failed to solve linear system"
                    stop
                endif
                !    print *, "info", info
                !    print *, y(1:num_wave_functions)
                do l = 1, N

                   if( shell_occupancy(idx(l,s)) == 0) stop 
                   ysol(idx(l,s)) = y(l) /shell_occupancy(idx(l,s))
                enddo
                ! print *, "y =", y(1:N)
            endif
            
        enddo
        
       
    end subroutine solve_linear_problem

    subroutine compute_kli_potential(grid_size,  exchange_potential)
    !----------------------------------------------------------------------------------
    ! calculation of the exchange potential in the KLI ( Kriegger - Li - Ifrate) 
    ! aproximation 
    ! Paper:  Krieger, J. B., Li, Y. & Iafrate, G. J. Phys. Rev. A 46, 5453–5458 (1992).
    ! ----------------------------------------------------------------------------------
        
    integer,  intent(in) ::  grid_size  ! number of grid points
    real(dp), intent(out) :: exchange_potential(ndmx, 2) 
    real(dp) :: work(ndmx,2), fact(ndmx), fact1
    real(dp) :: shift
    integer :: i, j, s, last_index

        ! check for valid input pased to the routine

        do j = 1, grid_size
            do i =  1, num_wave_functions
                if(psi(j,1,i) /= psi(j,1,i)) then
                    print *, "We got a problem!!!"
                    print *, "Invalid wavefunction passed"
                    stop    
                endif
                if( oc(i) > 0 .and. rho(j, get_spin(i)) < tiny(1.0_dp)) then
                        print *, "We got a problem here"
                        print *, "density is zero for spin=", s, i, grid%mesh
                        stop
                endif
            enddo
        enddo

        
        ! compute the N_up and N_down wawefunctions
        ! N_up =  num_wf(1)
        ! N_down =  num_wf(2)
        call compute_num_wf(num_wave_functions, num_wf)
        
        mat_m = 0.0_dp
        ! compute matrix M_sigma 
        ! M(:,:,1)  Matrix up
        ! M(:,:,2)  Matrix down
        call compute_mat_m(num_wave_functions, grid_size,rho)    

        call compute_ux_kli(grid_size)
        ! compute the average slater potential
        call compute_potential_s(grid_size)
        call compute_average_ux_kli(grid_size)
        ! solve proble A x = b with A =  I - M
        call solve_linear_problem(num_wave_functions, average_kli_potential)   
        ! fill the output array for the potential
        exchange_potential  = 0
        do s = 1,nspin
            work = 0
            do i = 1, num_wf(s) - 1
                fact = ysol(idx(i,s))  !average_kli_potential(idx(i,s)) - average_ux_kli(idx(i,s))
                ! work =  work + abs(psi(:,1,i))**2 * fact * shell_occupancy(i)
                work(:,s) =  work(:,s) + psi(:,1,idx(i,s))**2 * fact * shell_occupancy(idx(i,s))
            ! work =  work + psi(:,1,i)  * fact
            enddo   
            ! last_index = idx(num_wf(s),s)
            ! shift =  psi(:,1,idx(i,s))**2 * ysol(last_index) * (shell_occupancy(last_index) - 1)
            shift = 0
            if(num_wf(s) > 0) then ! s can be zero as in the hidrogen case
                work(:,s) =  (work(:,s) +  shift)/rho(:,s)
                exchange_potential(:,s) = (slater_potential(:,s)  +  work(:,s))
            endif
        enddo
        
        ! save exchange potential to file
        call savetxtv2("kli_pot.up", grid%r, exchange_potential(:,1))
        if (nspin > 1) call savetxtv2("kli_pot.dw", grid%r, exchange_potential(:,2))
        
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