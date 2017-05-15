!
! Copyright (C) 2004 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------
subroutine ascheq(nn,lam,energy,mesh,grid,vpot,ze2,thresh0, solution, nstop)
  !---------------------------------------------------------------
  !
  !  numerical integration of the radial schroedinger equation for
  !  bound states in a local potential.
  !  thresh determines the absolute accuracy for the eigenvalue
  !
  use kinds, only : DP
  use radial_grids, only: radial_grid_type, series
  implicit none

  type(radial_grid_type), intent(in) :: grid

  integer, intent(in) :: mesh,&
                        lam ! angular momentum value l
  integer, intent(in) :: nn   ! main quantum number
  real(dp), intent(in) :: vpot(mesh) 
  real(dp), intent(inout) :: energy  ! energy eigenvalue
  real(dp), intent(out) :: solution(mesh) ! eigenfunction

  integer:: nstop, &
            iter, &
            l1,& 
            i, &
            ik, &
            ncross, &
            n, &
            nstart, &
            ns, &
            ndcr

  real(DP) :: ze2, &
             ddx12,&
             eup, &
             elw, &
             b0e, &
             ymx, &
             rap, &
             rstart, &
             expn,  &
             c1, &
             c2, &
             fe, &
             sum0, &
             f2, &
             sum, &
             sqlhf,&
             f0,&
             f1, &
             dfe, &
             de,&
             eps, &
             yln, &
             xp, &
             xl1, &
             x4l6

  real(DP),allocatable:: c(:), el(:), f(:)
  real(DP):: b(0:3), thresh0, thresh

  integer, parameter ::  maxter =  50
  !
  !  set up constants and initialize
  !
  call allocate_vars()
 
  thresh = thresh0
  if (energy < -5.e+2) thresh = thresh0 * 10.0_DP
  
  ddx12 = grid%dx * grid%dx / 12.0_dp
  l1 = lam + 1
  sqlhf = (DBLE(lam) + 0.5_dp)**2
  ndcr = nn - lam - 1

  !
  !  set initial lower and upper bounds to the eigenvalue
  !
  eup = vpot(mesh) + sqlhf / grid%r2(mesh)
  elw = eup

  do i = 1, mesh
     elw = min(elw, vpot(i) + sqlhf/grid%r2(i))
  enddo

  if(eup .eq. elw) call stop_error(200)

  call adjust_energy(elw, eup, energy)

  !
  !  series development of the potential near the origin
  !
  do i = 1, 4
     solution(i) = vpot(i) - ze2 / grid%r(i)
  enddo

  call series(solution,grid%r,grid%r2,b)
  
  
  do iter = 1, maxter
      
      !
      !  set up the f-function and determine the position of its last
      !  change of sign
      !  f < 0 (approximatively) means classically allowed   region
      !  f > 0         "           "        "      forbidden   "
      !
      f(1) = ddx12*(grid%r2(1)*(vpot(1)-energy)+sqlhf)

      do i = 2,mesh
        f(i) = ddx12 *(grid%r2(i)*(vpot(i) - energy) + sqlhf)
        !  print  '(5f20.6)', f(i), vpot(i), sqlhf, grid%r(i), sqlhf
        if( f(i) * f(i-1) < 0 ) ik = i
      enddo

      !---------------------------------
      ! check for error and get out
      ! if(ik .ge. mesh-2) call stop_error(302)  
      !-------------------------------

      do i = 1, mesh
        f(i)=1.0_dp-f(i)
      enddo
      
      !
      solution(:) = 0.0_dp

      !
      !  determination of the wave-function in the first two points by
      !  series development
      !
      xl1 = lam + 1.0_DP
      x4l6= 4.0_dp *lam + 6.0_dp
      b0e = b(0)-energy
      c1 =0.5_dp*ze2/xl1
      c2 = (c1*ze2+b0e)/x4l6
      call start_scheq( lam, energy, b, grid, ze2, solution)
      !
      !  start outward integration and count number of crossings
      !
      ncross=0
      ymx=0.0_dp
      do n = 2, ik - 1
        solution(n+1)=((12.0_dp-10.0_dp*f(n))*solution(n)-f(n-1)*solution(n-1))/f(n+1)
        if ( solution(n) * solution(n+1) < 0 ) ncross=ncross+1
        ymx = max (ymx, abs(solution(n+1)))
      end do
      !
      !  matching radius has been reached going out. if ncross is not
      !  equal to ndcr, modify the trial eigenvalue.
      !
      if(ndcr < ncross) then
        !
        !  too many crossings. e is an upper bound to the true eigen-
        !  value. increase abs(e)
        !
        eup = energy
        rap = (DBLE(ncross+l1)/DBLE(nn))**2
        energy = ( energy - vpot(mesh))*rap+vpot(mesh)

        ! should apply only the lower check,
        ! my modification is applying both
        ! i dont know if this introduces any error
        call adjust_energy(elw, eup, energy)
        cycle ! skip rest of the loop body and start again
      else if (ndcr > ncross) then
        !
        !  too few crossings. e is a lower bound to the true eigen-
        !  value. decrease abs(e)
        !
        elw = energy
        rap=(DBLE(ncross+l1)/DBLE(nn))**2
        energy =(energy - vpot(mesh))*rap+vpot(mesh)
        if(energy  .gt. eup) energy=0.9_dp*eup+0.1_dp*elw
        cycle ! skip rest of the loop body and start again
      end if

      call inward_integration_setup()

      call integrate_inwards()
    

      !
      !  if necessary, improve the trial eigenvalue by the cooley's
      !  procedure. jw cooley math of comp 15,363(1961)
      !
      fe=(12.0_dp-10.0_dp*f(ik))*solution(ik)-f(ik-1)*solution(ik-1)-f(ik+1)*solution(ik+1)
      !
      !  calculate the normalization
      !
      if(ymx.ge.1.0e10_dp) solution = solution/ymx

      call compute_sum0()

      f2 = grid%r2(1  )* solution(1)*solution(1)
      sum = grid%r(1)*f2/DBLE(2*l1+1)

      do n = 1, nstart-2 ,2
        f0 = f2
        f1 = grid%r2(n+1)*solution(n+1)*solution(n+1)
        f2 = grid%r2(n+2)*solution(n+2)*solution(n+2)
        sum = sum + f0 + f2 + 4.0_DP * f1
      enddo

      sum = sum0 + grid%dx*sum/3.0_dp

      dfe=-solution(ik)*f(ik)/grid%dx/sum
      de=-fe*dfe
      eps=abs(de/energy)

      if(abs(de).lt.thresh) then
        ! convergence achieved
        exit ! exit cycle
      endif

      if(eps .gt. 0.25_dp) de=0.25_dp*de/eps

      if(de.gt.0.0_dp) elw= energy
      if(de.lt.0.0_dp) eup= energy

      energy =  energy + de
      call adjust_energy(elw, eup, energy)
      nstop = 50
  enddo ! main iteration loop

  if (iter ==  maxter) call stop_error(300)
  !
  !  normalize the eigenfunction and exit
  !

  do n = nstart, mesh - 1
     solution(n+1)=0.0_dp
     if(solution(n).eq.0.0_dp) cycle
     yln = log(abs(solution(n)))
     xp = -sqrt(12.0_dp*abs(1.0_dp-f(n)))
     expn = yln+xp
     if(expn.lt.-80.0_dp) cycle
     solution(n+1)=sign(exp(expn),solution(n))
  enddo

  call compute_sum()

  solution = grid%sqr * solution / sum
  
  ! if(nstop .lt. 100) then
  !   call stop_error(nstop)
  ! endif
  nstop = 0

  call deallocate_vars()

  contains 

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! internal subroutines
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine  adjust_energy(lower, upper, energy)
      real(DP),intent(in) :: lower, upper
      real(DP),intent(inout) :: energy
       if(energy .gt.upper) energy = 0.9_dp*upper+0.1_dp*lower
      if(energy.lt.lower) energy = 0.9_dp*lower+0.1_dp*upper
    end subroutine adjust_energy


    subroutine inward_integration_setup()
      !
      !  prepare inward integration
      !  charlotte froese can j phys 41,1895(1963)
      !
      !            start at  min( rmax, 10*rmatch )
      !
      implicit none
      real(DP) :: di  
      nstart=mesh
      ns = 10
      rstart = ns*grid%r(ik)

      if(rstart.lt.grid%r(mesh)) then
         do i=ik,mesh
            nstart=i
            if(grid%r(i).ge.rstart) exit
         enddo
         nstart=nstart/2
         nstart=2*nstart+1
      end if
      !
      !  set up a, l, and c vectors
      !
      n = ik + 1
      el(n)=10.0_dp*f(n)-12.0_dp
      c(n)=-f(ik)*solution(ik)
      
      do n= ik + 2, nstart
         di=10.0_dp*f(n)-12.0_dp
         el(n)=di-f(n)*f(n-1)/el(n-1)
         c(n)=-c(n-1)*f(n-1)/el(n-1)
      enddo
    end subroutine inward_integration_setup

    subroutine integrate_inwards()
      !
      !  start inward integration by the froese's tail procedure
      !

      expn=exp(-sqrt(12.0_dp*abs(1.0_dp-f(nstart-1))))

      solution(nstart-1)=c(nstart-1)/(el(nstart-1)+f(nstart)*expn)
      solution(nstart)=expn*solution(nstart-1)

      do n=nstart-2,ik+1,-1
        solution(n)=(c(n)-f(n+1)*solution(n+1))/el(n)
      end do

    end subroutine integrate_inwards

    subroutine allocate_vars()
      
      integer :: ierr
      allocate(c(mesh),stat=ierr)
      allocate(f(mesh),stat=ierr)
      allocate(el(mesh),stat=ierr)

    end subroutine allocate_vars


    subroutine deallocate_vars()
      if(allocated(el)) deallocate(el)
      if(allocated(f)) deallocate(f )
      if(allocated(c)) deallocate(c )
    end subroutine deallocate_vars

    subroutine compute_sum()
      
      real(DP) :: sum1
      sum1=0.0_dp

      do n=nstart,mesh-2,2
         f0=f2
         f1=grid%r2(n+1)*solution(n+1)*solution(n+1)
         f2=grid%r2(n+2)*solution(n+2)*solution(n+2)
         sum1= sum1 + f0 + f2+4.0_dp*f1
      enddo

      sum=sum+grid%dx*sum1/3.0_dp
      sum=sqrt(sum)

    end subroutine compute_sum

    subroutine compute_sum0()

      ! I don't have any idea what sum0 is
      real(dp) :: a0, a1, a2

      a0 = 1.0_dp/DBLE(2*lam+3)
      a1=c1/DBLE(lam+2)
      a2=(c1*c1+c2+c2)/DBLE(2*lam+5)
      sum0=(a0+grid%r(1)*(a1+grid%r(1)*a2))*grid%r(1)**(2*lam+3)

    end subroutine compute_sum0

    subroutine stop_error(code)
      !
      !  error exit
      !
      integer, intent(in) :: code
      call deallocate_vars()
      write(6,9000) code,nn,lam,elw,eup
      9000 format(5x,'error in ascheq: nstop =',i4,'. n l =',2i3,/ &
      & 5x,'elw =',f15.10,' eup =',f15.10)
      stop
    end subroutine stop_error
end subroutine ascheq
