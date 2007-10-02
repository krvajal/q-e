!
! Copyright (C) 2001-2007 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine read_ncpp (iunps, np, upf)
  !-----------------------------------------------------------------------
  !
  USE kinds, only: dp
  USE parameters, ONLY: nchix, lmaxx
  use atom,  only: numeric
  use pseud, only: cc, alpc, aps, alps, nlc, nnl, lmax, lloc, &
       a_nlcc, b_nlcc, alpha_nlcc
  use funct, only: set_dft_from_name, dft_is_meta, dft_is_hybrid
  USE pseudo_types

  implicit none
  !
  TYPE (pseudo_upf) :: upf
  integer :: iunps, np
  !
  real(DP) :: x, vll
  real(DP), allocatable:: vnl(:,:)
  real(DP), parameter :: rcut = 10.d0, e2 = 2.d0
  real(DP), external :: erf
  integer :: nb, ios, i, l, ir
  logical :: bhstype
  !
  !====================================================================
  ! read norm-conserving PPs
  !
  read (iunps, '(a2)', end=300, err=300, iostat=ios) upf%dft
  if (upf%dft(1:2) .eq.'**') upf%dft = 'PZ'
  read (iunps, *, err=300, iostat=ios) upf%psd, upf%zp, lmax(np), nlc(np), &
                                       nnl(np), upf%nlcc, lloc(np), bhstype
  if (nlc(np) > 2 .or. nnl(np) > 3) &
       call errore ('read_ncpp', 'Wrong nlc or nnl', np)
  if (nlc(np)*nnl(np) < 0) call errore ('read_ncpp', 'nlc*nnl < 0 ? ', np)
  if (upf%zp <= 0d0 .or. upf%zp > 100 ) &
       call errore ('read_ncpp', 'Wrong zp ', np)
  !
  !   In numeric pseudopotentials both nlc and nnl are zero.
  !
  numeric (np) = (nlc (np) <= 0) .and. (nnl (np) <= 0)
  !
  if (lloc (np) == -1000) lloc (np) = lmax (np)
  if (lloc (np) < 0 .or. lmax(np) < 0 .or. &
       .not.numeric(np) .and. (lloc(np) > min(lmax(np)+1,lmaxx+1) .or. &
       lmax(np) > max(lmaxx,lloc(np))) .or. &
       numeric(np) .and. (lloc(np) > lmax(np) .or. lmax(np) > lmaxx) ) &
       call errore ('read_ncpp', 'wrong lmax and/or lloc', np)
  if (.not.numeric (np) ) then
     !
     !   read here pseudopotentials in analytic form
     !
     read (iunps, *, err=300, iostat=ios) &
          (alpc(i,np), i=1,2), (cc(i,np), i=1,2)
     if (abs (cc(1,np)+cc(2,np)-1.d0) > 1.0d-6) &
          call errore ('read_ncpp', 'wrong pseudopotential coefficients', 1)
     do l = 0, lmax (np)
        read (iunps, *, err=300, iostat=ios) (alps(i,l,np), i=1,3), &
                                             (aps(i,l,np),  i=1,6)
     enddo
     if (upf%nlcc ) then
        read (iunps, *, err=300, iostat=ios) &
             a_nlcc(np), b_nlcc(np), alpha_nlcc(np)
        if (alpha_nlcc(np) <= 0.d0) call errore('read_ncpp','alpha_nlcc=0',np)
     endif
  endif
  read (iunps, *, err=300, iostat=ios) upf%zmesh, upf%xmin, upf%dx, &
                                       upf%mesh, upf%nwfc 
  if ( upf%mesh <= 0) &
       call errore ('read_ncpp', 'wrong nuymber of mesh points', np)
  if ( upf%nwfc > nchix .or. &
       (upf%nwfc < lmax(np)   .and. lloc(np) == lmax(np)) .or. & 
       (upf%nwfc < lmax(np)+1 .and. lloc(np) /= lmax(np)) ) &
       call errore ('read_ncpp', 'wrong no. of wfcts', np)
  !
  !  Here pseudopotentials in numeric form are read
  !
  ALLOCATE ( upf%chi(upf%mesh,upf%nwfc), upf%rho_atc(upf%mesh) )
  ALLOCATE ( upf%lchi(upf%nwfc), upf%oc(upf%nwfc) )
  allocate (vnl(upf%mesh, 0:lmax(np)))
  if (numeric (np) ) then
     do l = 0, lmax (np)
        read (iunps, '(a)', err=300, iostat=ios)
        read (iunps, *, err=300, iostat=ios) (vnl(ir,l), ir=1,upf%mesh )
     enddo
     if ( upf%nlcc ) then
        read (iunps, *, err=300, iostat=ios) (upf%rho_atc(ir), ir=1,upf%mesh)
     endif
  endif
  !
  !  Here pseudowavefunctions (in numeric form) are read
  !
  do nb = 1, upf%nwfc
     read (iunps, '(a)', err=300, iostat=ios)
     read (iunps, *, err=300, iostat=ios) upf%lchi(nb), upf%oc(nb)
     !
     !     Test lchi and occupation numbers
     !
     if (nb <= lmax(np) .and. upf%lchi(nb)+1 /= nb) &
          call errore ('read_ncpp', 'order of wavefunctions', 1)
     if (upf%lchi(nb) > lmaxx .or. upf%lchi(nb) < 0) &
                      call errore ('read_ncpp', 'wrong lchi', np)
     if (upf%oc(nb) < 0.d0 .or. upf%oc(nb) > 2.d0*(2*upf%lchi(nb)+1)) &
          call errore ('read_ncpp', 'wrong oc', np)
     read (iunps, *, err=300, iostat=ios) ( upf%chi(ir,nb), ir=1,upf%mesh )
  enddo
  !
  !====================================================================
  ! PP read: now setup 
  !
  call set_dft_from_name( upf%dft )
  !
#if defined (EXX)
#else
  IF ( dft_is_hybrid() ) &
    CALL errore( 'read_ncpp ', 'HYBRID XC not implemented in PWscf', 1 )
#endif
  !
  !    calculate the number of beta functions
  !
  upf%nbeta = 0
  do l = 0, lmax (np)
     if (l /= lloc (np) ) upf%nbeta = upf%nbeta + 1
  enddo
  ALLOCATE ( upf%lll(upf%nbeta) )
  nb = 0
  do l = 0, lmax (np)
     if (l /= lloc (np) ) then
        nb = nb + 1 
        upf%lll (nb) = l
     end if
  enddo
  !
  !    compute the radial mesh
  !
  ALLOCATE ( upf%r(upf%mesh), upf%rab(upf%mesh) )
  do ir = 1, upf%mesh
     x = upf%xmin + DBLE (ir - 1) * upf%dx
     upf%r(ir) = exp (x) / upf%zmesh 
     upf%rab(ir) = upf%dx * upf%r(ir)
  enddo
  ALLOCATE ( upf%kkbeta(upf%nbeta) )
  do ir = 1, upf%mesh
     if ( upf%r(ir) > rcut) then
        upf%kkbeta(:) = ir
        go to 5
     end if
  end do
  upf%kkbeta(:) = upf%mesh
  !
  ! ... force kkbeta to be odd for simpson integration (obsolete?)
  !
5 upf%kkbeta(:) = 2 * ( ( upf%kkbeta(:) + 1 ) / 2) - 1
  !
  ALLOCATE ( upf%vloc(upf%mesh) )
  upf%vloc (:) = 0.d0
  if (.not. numeric(np)) then
     !
     ! bring analytic potentials into numerical form
     !
     IF ( nlc(np) == 2 .AND. nnl(np) == 3 .AND. bhstype ) &
          CALL bachel( alps(1,0,np), aps(1,0,np), 1, lmax(np) )
     !
     do i = 1, nlc (np)
        do ir = 1, upf%kkbeta(1)
           upf%vloc (ir) = upf%vloc (ir) - upf%zp * e2 * cc (i, np) * &
               erf ( sqrt (alpc(i,np)) * upf%r(ir) ) / upf%r(ir)
        end do
     end do
     do l = 0, lmax (np)
        vnl (:, l) = upf%vloc (1:upf%mesh)
        do i = 1, nnl (np)
           vnl (:, l) = vnl (:, l) + e2 * (aps (i, l, np) + &
                   aps (i + 3, l, np) * upf%r (:) **2) * &
                   exp ( - upf%r(:) **2 * alps (i, l, np) )
        enddo
     enddo
     ! core corrections are still analytic!
     !!! numeric(np) =.true.
  end if
  !
  ! assume l=lloc as local part and subtract from the other channels
  !
  if (lloc (np) <= lmax (np) ) &
    upf%vloc (:) = vnl (:, lloc (np))
  ! lloc > lmax is allowed for PP in analytical form only
  ! it means that only the erf part is taken as local part 
  do l = 0, lmax (np)
     if (l /= lloc(np)) vnl (:, l) = vnl(:, l) - upf%vloc(:)
  enddo
  !
  !    compute the atomic charges
  !
  ALLOCATE ( upf%rho_at (upf%mesh) )
  upf%rho_at(:) = 0.d0
  do nb = 1, upf%nwfc
     if ( upf%oc(nb) > 0.d0) then
        do ir = 1, upf%mesh
           upf%rho_at(ir) = upf%rho_at(ir) + upf%oc(nb) * upf%chi(ir,nb)**2
        enddo
     endif
  enddo
  !====================================================================
  ! convert to separable (KB) form
  !
  ALLOCATE ( upf%beta (upf%mesh, upf%nbeta) ) 
  ALLOCATE ( upf%dion (upf%nbeta,upf%nbeta), upf%lll (upf%nbeta) ) 
  upf%dion (:,:) = 0.d0
  nb = 0
  do l = 0, lmax (np)
     if (l /= lloc (np) ) then
        nb = nb + 1
        ! betar is used here as work space
        do ir = 1, upf%kkbeta(1)
           upf%beta (ir, nb) = upf%chi(ir, l+1) **2 * vnl(ir, l)
        end do
        call simpson (upf%kkbeta (1), upf%beta (1, nb), upf%rab, vll )
        upf%dion (nb, nb) = 1.d0 / vll
        ! upf%beta stores projectors  |beta(r)> = |V_nl(r)phi(r)>
        do ir = 1, upf%kkbeta (1)
           upf%beta (ir, nb) = vnl (ir, l) * upf%chi (ir, l + 1)
        enddo
        upf%lll (nb) = l
     endif
  enddo
  deallocate (vnl)
  !
  ! for compatibility with USPP
  !
  upf%nqf = 0
  upf%nqlc= 0
  upf%tvanp =.false.
  upf%has_so=.false.
  !
  return

300 call errore ('read_ncpp', 'pseudo file is empty or wrong', abs (np) )
end subroutine read_ncpp

