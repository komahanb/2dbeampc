program problemPC

  use dimpce,only:probtype,id_proc

  implicit none
  !
  !     include the Ipopt return codes
  !
  include 'IpReturnCodes.inc'
  include 'mpif.h'
  !
  !     Size of the problem (number of variables and equality constraints)
  !
  integer     N,     M,     NELE_JAC,     NELE_HESS,      IDX_STY
  parameter  (N = 2, M = 3, NELE_JAC = 6, NELE_HESS = 3)
  parameter  (IDX_STY = 1 )
  !
  !     Space for multipliers and constraints
  !
  double precision LAM(M)
  double precision G(M)
  !
  !     Vector of variables
  !
  double precision X(N)
  !
  !     Vector of lower and upper bounds
  !
  double precision X_L(N), X_U(N), Z_L(N), Z_U(N)
  double precision G_L(M), G_U(M)
  !
  !     Private data for evaluation routines
  !     This could be used to pass double precision and integer arrays untouched
  !     to the evaluation subroutines EVAL_*
  !
  double precision DAT(2000)
  integer IDAT(2000)
  !
  !     Place for storing the Ipopt Problem Handle
  !
  integer*8 IPROBLEM
  integer*8 IPCREATE
  !
  integer IERR
  integer IPSOLVE, IPADDSTROPTION
  integer IPADDNUMOPTION, IPADDINTOPTION
  integer IPOPENOUTPUTFILE
  !
  double precision F,Fs,sigmax(N),pi
  integer i,kprob

  double precision  infbound
  parameter        (infbound = 1.d+20)
  !
  !     The following are the Fortran routines for computing the model
  !     functions and their derivatives - their code can be found further
  !     down in this file.
  !
  external EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS, ITER_CB


  call MPI_START

  pi=4.0*atan(1.0) ! constant for later use (visible globally)

!!$       !Store variables into DAT array to passthrough
!!$       ! Use this in your calling program
  dat(1000+1)=40.0d3 !Bending Moment M
  dat(1000+2)=150.0d3 ! Shear Force V
  dat(1000+3)=10.0d6 ! Max bending
  dat(1000+4)=2.0d6 ! max shear stress
  dat(1000+5)=1.0 !Factor of safety
  dat(1000+20)=77      ! filenum for PC

  !============
  !  DAT array
  !============
  
  ! 1 to N are used to store sigmax

  ! N+1--> fmeantmp
  ! N+2--> fvartmp

  ! N+2+1 to 2*N+2  --> fmeanprime(i)+fvarprime(i)
  ! 2N*+2+1 to 3N+2 --> x(i)



  !1000+ is used to pass data to PC 
  
  !==================!
  !  IDAT array      !
  !==================!
  !  IDAT(1)=kprob   !
  !  IDAT(2)=0       !
  !  IDAT(3)=probtype!
  !==================!

  !Other IPOPT params

  probtype(:)=1

  kprob=4

  ! SD for area design variables

  sigmax(1)=0.05
  sigmax(2)=0.05

  do i=i,n
     dat(i)=sigmax(i)
  end do

  IDAT(1)=kprob
  IDAT(2)=0
  IDAT(3:N+2)=probtype(1:N)

  ! Area design variables

  do i=1,N
     X(i)   = 1.0  
     X_L(i) = 0.25 
     X_U(i) = infbound 
  end do

  !
  !     Set bounds for the constraints
  !

  do i=1,M
     G_L(i)=-infbound
     G_U(i)=0.d0
  end do

  !
  !     First create a handle for the Ipopt problem (and read the options
  !     file)
  !


  IPROBLEM = IPCREATE(N, X_L, X_U, M, G_L, G_U, NELE_JAC, NELE_HESS,IDX_STY, EV_F, EV_G, EV_GRAD_F, EV_JAC_G, EV_HESS)
  if (IPROBLEM.eq.0) then
     write(*,*) 'Error creating an Ipopt Problem handle.'
     call stop_all
  endif
  !
  !     Open an output file
  !
  
if (id_proc.eq.0) open(unit=86,file='Opt.his',form='formatted',status='replace')


  IERR = IPOPENOUTPUTFILE(IPROBLEM, 'IPOPT.OUT', 5)
  if (IERR.ne.0 ) then
     write(*,*) 'Error opening the Ipopt output file.'
     goto 9000
  endif
  !

  !!
  !!     Set a callback function to give you control once per iteration.
  !!     You can use it if you want to generate some output, or to stop
  !!     the optimization early.
  !!

  call IPSETCALLBACK(IPROBLEM, ITER_CB)

  !
  !     Call optimization routine
  !

  if (id_proc.eq.0) then
     IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
     if (IERR.ne.0 ) goto 9990
  else
     IERR = IPADDINTOPTION(IPROBLEM, 'print_level', 0)
     if (IERR.ne.0 ) goto 9990
  end if
  


  IERR = IPSOLVE(IPROBLEM, X, G, F, LAM, Z_L, Z_U, IDAT, DAT)

  !
  !     Output:
  !
  if (id_proc.eq.0) then

     if( IERR.eq.IP_SOLVE_SUCCEEDED .or. IERR.eq.5) then
        write(*,*)
        write(*,*) 'The solution was found.'
        write(*,*)
     else
        write(*,*)
        write(*,*) 'An error occoured.'
        write(*,*) 'The error code is ',IERR
        write(*,*)
     endif

     write(*,*) 'The final value of the objective function is ',F
     write(*,*)
     write(*,*) 'The optimal values of X are:'
     write(*,*)
     do i = 1, N
        write(*,*) 'X  (',i,') = ',X(i),'m^2'
        end if
     enddo
     write(*,*)
     write(*,*) 'The multipliers for the equality constraints are:'
     write(*,*)
     do i = 1, M
        write(*,*) 'LAM(',i,') = ',LAM(i)
     enddo
     write(*,*)
     write(*,*) 'Weight and its variance:',DAT(N+1),DAT(N+2)

  end if
  !
9000 continue

  !
  !     Clean up
  !

  call IPFREE(IPROBLEM)

  if (id_proc.eq.0) close(86)

  call stop_all
  !
9990 continue
  write(*,*) 'Error setting an option'
  goto 9000

end program problemPC
!
! =============================================================================
!
!                    Computation of objective function
!
! =============================================================================
!

subroutine EV_F(N, X, NEW_X, F, IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc

  implicit none
  integer N, NEW_X,I
  double precision F, X(N),sigmax(N),fmeantmp,fvartmp,fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision DAT(*)
  integer IDAT(*),kprob,NMC
  integer IERR
  double precision fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
  double precision  rho, L, sigmay, pi, p, E, Fs 
  integer::myflag(10) 

  kprob=IDAT(1)
  probtype(1:N)=IDAT(3:N+2)

  do i=1,N
     sigmax(i)=DAT(i)
     Xsave(i)=X(i)
  end do

  !---- MEAN and VARIANCE OF worst OBJECTIVE FUNCTION

  !call  PCestimate(dim,xavgin,xstdin,fctin,fctindxin,DATIN,orderinitial,orderfinal,statin,probtypeIN,sampfac,fmeanout,fvarout,fmeanprimeout,fvarprimeout,fmeandbleprimeout,fvardbleprimeout)

  call  PCestimate(N,x,sigmax,10,0,DAT(1001:1020),1,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)


  if (IDAT(2).eq.1) then ! Deterministic with PC
     fvartmp=0.0d0
     fvarprimetmp=0.0d0
  end if

  dat(N+1)=fmeantmp
  dat(N+2)=fvartmp

  ! Store gradients for later use

  do i=1,N
     dat(N+2+i)=fmeanprimetmp(i)+fvarprimetmp(i)
  end do

  ! Store x
  do i=1,n
     DAT(2*N+2+i)=Xsave(i)
     X(i)=Xsave(i)
  end do

  if (id_proc.eq.0) then
     print*,''
     write(*,'(4x,a,3F13.4)') '>>Objective:',fmeantmp,fvartmp,fmeantmp+fvartmp
     print*,''
  end if

  !---- COMBINED OBJECTIVE FUNCTION

  F=fmeantmp+fvartmp

  IERR = 0
  return

end subroutine EV_F

!
! =============================================================================
!
!                     Computation of constraints
!
! =============================================================================
!

subroutine EV_G(N, X, NEW_X, M, G, IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc

  implicit none
  integer N, NEW_X, M
  double precision G(M), X(N), sigmax(N), cmean(M), cstd(M), fmeantmp, fvartmp
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n),dc(M,N)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  integer IDAT(*),kprob,NMC
  integer IERR, i, j, cnt
  double precision::fmin,fmax,gradmin(N-1),gradmax(N-1),gtol,low(N-1),up(N-1),Xsave(N)
  double precision :: rho, L, sigmay, pi, p, E, Fs
  integer::myflag(10) 

  kprob=IDAT(1)
  probtype(1:N)=IDAT(3:N+2)

  do i=1,N
     sigmax(i)=DAT(i)
     Xsave(i)=X(i)
  end do

  do i=1,M

     !---- MEAN OF INEQUALITY CONSTRAINT i
     !call  PCestimate(dim,xavgin,xstdin,fctin,fctindxin,DATIN,orderinitial,orderfinal,statin,probtypeIN,sampfac,fmeanout,fvarout,fmeanprimeout,fvarprimeout,fmeandbleprimeout,fvardbleprimeout)
!     if (id_proc.eq.0) print*,x
     call  PCestimate(N,x,sigmax,10,i,DAT(1001:1020),1,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
!     if (id_proc.eq.0) print*,fmeanprimetmp,fvarprimetmp

     cmean(i)=fmeantmp
     cstd(i)=sqrt(fvartmp)

     do j=1,N
        dc(i,j)=fmeanprimetmp(j)
        if (fvartmp.ne.0.0) then
           dc(i,j)=dc(i,j)+dble(kprob)*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
        endif
     end do

     if (IDAT(2).eq.1) then ! Deterministic with PC
        fvartmp=0.0d0
        fvarprimetmp=0.0d0
     end if

  end do

  G(1:M)=cmean(1:M)+dble(kprob)*cstd(1:M)

  !Just printing

  if (id_proc.eq.0) then
     print*,''
     write(*,'(4x,a)') '>>Normalized Constraint Values:'
     do i=1,M
        write(*,'(E13.2)'),g(i)
     end do
     print*,''
  end if

  do i=1,N
     DAT(3*N+2+i)=Xsave(i)
     X(i)=Xsave(i)
  end do

  cnt=0
  do i=1,M
     do j=1,N
        cnt=cnt+1
        DAT(4*N+2+cnt)=dc(i,j)
     end do
  end do


  IERR = 0
  return
end subroutine EV_G

!
! =============================================================================
!
!                Computation of gradient of objective function
!
! =============================================================================
!

subroutine EV_GRAD_F(N, X, NEW_X, GRAD, IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc

  implicit none
  integer N, NEW_X,i
  double precision GRAD(N), X(N), sigmax(N), fmeantmp, fvartmp
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  integer IDAT(*),kprob,NMC
  integer IERR
  double precision  rho, L, sigmay, pi, p, E, Fs 
  integer::myflag(10) 
  logical samex

  
  samex=.true.
  do i=1,N
     if (x(i).ne.DAT(2*N+2+i)) samex=.false. 
  end do
  
  if (samex) then

     if (id_proc.eq.0) print *,'Samex in obj'
     !---- TOTAL GRADIENT OF OBJECTIVE FUNCTION
     do i=1,n
        GRAD(i)=DAT(N+2+i)
     end do
     if (id_proc.eq.0) print *,'grad obj',grad

  else

     kprob=IDAT(1)
     probtype(1:N)=IDAT(3:N+2)

     do i=1,N
        sigmax(i)=DAT(i)
     end do
     
     !---- MEAN OF INEQUALITY CONSTRAINT i
     !call  PCestimate(dim,xavgin,xstdin,fctin,fctindxin,DATIN,orderinitial,orderfinal,statin,probtypeIN,sampfac,fmeanout,fvarout,fmeanprimeout,fvarprimeout,fmeandbleprimeout,fvardbleprimeout)
     
     call  PCestimate(N,x,sigmax,10,0,DAT(1001:1020),1,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)

     !---- GRADIENT OF OBJECTIVE FUNCTION

     if (id_proc.eq.0)then
        print *,' >> Gradient of obj:'
     end if

     do i=1,n

        if (IDAT(2).eq.1) then ! Deterministic with PC
           fvartmp=0.0d0
           fvarprimetmp=0.0d0
        end if

        if (id_proc.eq.0)  print*,fmeanprimetmp(i),fvarprimetmp(i)

        GRAD(i)=fmeanprimetmp(i)+fvarprimetmp(i)

     end do

  end if

  if (id_proc.eq.0)print *,''
  
  IERR = 0
  return
end subroutine EV_GRAD_F

!
! =============================================================================
!
!                Computation of Jacobian of constraints
!
! =============================================================================
!
subroutine EV_JAC_G(TASK, N, X, NEW_X, M, NZ, ACON, AVAR, A,IDAT, DAT, IERR)
  use dimpce,only:probtype,id_proc

  implicit none
  integer TASK, N, NEW_X, M, NZ
  double precision X(N), A(NZ),cgrad(M,N), sigmax(N), fmeantmp, fvartmp
  integer ACON(NZ), AVAR(NZ), I, J, K, cnt, NMC
  double precision DAT(*),fmeanprimetmp(n),fvarprimetmp(n)
  real*8 :: fmeandbleprimetmp(n,n),fvardbleprimetmp(n,n)
  double precision  rho, L, sigmay, pi, p, E, Fs
  integer IDAT(*)
  integer IERR, kprob
  logical samex
  integer::myflag(10) 
  
  
  if( TASK.eq.0 ) then 
     !
     !     structure of Jacobian:
     !
     ACON(1) = 1
     AVAR(1) = 1

     ACON(2) = 1
     AVAR(2) = 2

     ACON(3) = 2 
     AVAR(3) = 1

     ACON(4) = 2
     AVAR(4) = 2

     ACON(5) = 3
     AVAR(5) = 1

     ACON(6) = 3
     AVAR(6) = 2
  else
         
         
         samex=.true.
         do i=1,N
            if (x(i).ne.DAT(3*N+2+i)) samex=.false. 
         end do

         if (samex) then

            if (id_proc.eq.0) print *,'Samex in con'

            cnt=0
            do i=1,M
               do j=1,N
                  cnt=cnt+1
                  cgrad(i,j)=DAT(4*N+2+cnt)
               end do
            end do

         else !not samex

            !---- TOTAL GRADIENT OF CONSTRAINTS 

            kprob=IDAT(1)
            probtype(1:N)=IDAT(3:N+2)

            do i=1,N
               sigmax(i)=DAT(i)
            end do


            cgrad(:,:)=0.0

            do i=1,M

               !---- MEAN OF INEQUALITY CONSTRAINT i
               !call  PCestimate(dim,xavgin,xstdin,fctin,fctindxin,DATIN,orderinitial,orderfinal,statin,probtypeIN,sampfac,fmeanout,fvarout,fmeanprimeout,fvarprimeout,fmeandbleprimeout,fvardbleprimeout)

               call  PCestimate(N,x,sigmax,10,i,DAT(1001:1020),1,3,3,0,probtype,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)

               if (IDAT(2).eq.1) then ! Deterministic with PC
                  fvartmp=0.0d0
                  fvarprimetmp=0.0d0
               end if


               do j=1,N
                  cgrad(i,j)=fmeanprimetmp(j)
                  if (fvartmp.ne.0.0) then
                     cgrad(i,j)=cgrad(i,j)+dble(kprob)*fvarprimetmp(j)/(2.0*sqrt(fvartmp))
                  endif
               end do
            end do
         end if
            !if (id_proc.eq.0) print *,'Cons Gradients',jac(1:6)
         
        
         A(1)=cgrad(1,1)
         A(2)=cgrad(1,2)

         A(3)=cgrad(2,1)
         A(4)=cgrad(2,2)

         A(5)=cgrad(3,1)
         A(6)=cgrad(3,2)

      end if

         IERR = 0
  return
end subroutine EV_JAC_G
!
! =============================================================================
!
!                Computation of Hessian of Lagrangian
!
! =============================================================================
!
subroutine EV_HESS(TASK, N, X, NEW_X, OBJFACT, M, LAM, NEW_LAM,NNZH, IRNH, ICNH, HESS, IDAT, DAT, IERR)
  implicit none
  integer TASK, N, NEW_X, M, NEW_LAM, NNZH,  i,j,ii
  double precision X(N), OBJFACT, LAM(M), HESS(NNZH), sigmax(N)
  integer IRNH(NNZH), ICNH(NNZH)
  double precision::fmeantmp,fvartmp
  double precision OBJHESS(NNZH),CONHESS(M,NNZH)
  double precision DAT(*)
  integer IDAT(*), kprob
  integer IERR
  integer::myflag(10) 

  if( TASK.eq.0 ) then
     !
     !     structure of sparse Hessian (lower triangle):
     !
     
     IRNH(1) = 1
     ICNH(1) = 1

     IRNH(2) = 2
     ICNH(2) = 2

     IRNH(3) = 1
     ICNH(3) = 2

  else

!!$     
!!$     do ii=0,m
!!$
!!$        !      call PCestimate(dim,xavgin,xstdin,fctin,fctindxin,orderinitial,orderfinal,statin,fmeanout,fvarout,fmeanprimeout,fvarprimeout)
!!$        call  PCestimate(N,x,sigmax,11,ii,3,3,1,fmeantmp,fvartmp,fmeanprimetmp,fvarprimetmp,fmeandbleprimetmp,fvardbleprimetmp)
!!$
!!$        if (ii.eq.0) then
!!$           
!!$           cnt=0
!!$           do i=1,N
!!$              do j=1,N
!!$                 if (i.le.j) then
!!$                    cnt=cnt+1
!!$                    objhess(cnt)=fmeandbleprimetmp(i,j)+kprob*fvardbleprimetmp(i,j)
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$
!!$        else
!!$           
!!$           cnt=0
!!$           do i=1,N
!!$              do j=1,N
!!$                 if (i.le.j) then
!!$                    cnt=cnt+1
!!$                    conhess(ii,cnt)=fmeandbleprimetmp(i,j)+kprob*fvardbleprimetmp(i,j)
!!$                 end if
!!$              end do
!!$           end do
!!$
!!$        end if
!!$
!!$     end do
!!$
!!$     ! Assemble all into HESS
!!$     
!!$     HESS(:)=0.0
!!$     do i=1,NNZH
!!$        hesstmp=0.0
!!$        do j=1,m
!!$           hesstmp=hesstmp+lam(j)*conhess(j,i)
!!$        end do
!!$        hess(i)=hesstmp+objhess(i)
!!$     end do
     
     IERR = 0

  endif

  return
end subroutine EV_HESS











!
! =============================================================================
!
!                   Callback method called once per iteration
!
! =============================================================================
!
subroutine ITER_CB(ALG_MODE, ITER_COUNT,OBJVAL, INF_PR, INF_DU,MU, DNORM, REGU_SIZE, ALPHA_DU, ALPHA_PR, LS_TRIAL, IDAT,DAT, ISTOP)
  use dimpce,only:probtype,id_proc

  implicit none
  integer ALG_MODE, ITER_COUNT, LS_TRIAL
  double precision OBJVAL, INF_PR, INF_DU, MU, DNORM, REGU_SIZE
  double precision ALPHA_DU, ALPHA_PR
  double precision DAT(*)
  integer IDAT(*)
  integer ISTOP
  !
  !     You can put some output here
  !
  if (id_proc.eq.0) then

     if (ITER_COUNT .eq.0) then
        write(*,*) 
        write(*,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
     end if

     write(*,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU
     write(86,*) 'iter    objective      ||grad||        inf_pr          inf_du         lg(mu)'
     write(86,'(i5,5e15.7)') ITER_COUNT,OBJVAL,DNORM,INF_PR,INF_DU,MU

  end if
  !
  !     And set ISTOP to 1 if you want Ipopt to stop now.  Below is just a
  !     simple example.
  !
  
  if (ITER_COUNT .gt. 1 .and. DNORM.le.1.0D-02.and.inf_pr.le.5.0d-03) then

     ISTOP = 1

  end if

  return
end subroutine ITER_CB
