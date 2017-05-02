module kli
use kinds, ONLY : dp
use utils
use ld1inc, only: psi,& ! all electron wavefuntions
                  rho, &   ! all electron density rho(:,1)(up), rho(:,2) (down)
                  num_wave_functions => nwf,
                  get_spin => isw 


implicit none
private 
public compute_kli_potential
real(dp), allocatable :: mat_m_up(:,:),mat_m_down(:,:)
integer :: num_up_wavefunctions, num_down_wavefunctions

contains



    subroutine init_module()
        implicit none
        integer :: i
        num_up_wavefunctions = 0
        num_down_wavefunctions = 0
        do i=1,num_wave_functions
            if(get_spin(i) == 1) num_up_wavefunctions = num_up_wavefunctions + 1
            if(get_spin(i) == 2) num_down_wavefunctions = num_down_wavefunctions + 1
        enddo
        call assert( (num_up_wavefunctions + num_down_wavefunctions) == num_wave_functions)

        allocate(mat_m_up(num_up_wavefunctions,num_up_wavefunctions))
        allocate(mat_m_down(num_down_wavefunctions, num_down_wavefunctions)))
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