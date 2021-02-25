MODULE DEBUG
  contains
  SUBROUTINE FATAL(bool_condition,sub_name,error_name)
  ! This subroutine checks if bool_condition is true and if it is it stops the execution of the code.
  ! Additional optional inputs are the subroutine where the check is performed (sub_name) and an error message (error_name)
      implicit none
      LOGICAL,INTENT(in)                     :: bool_condition
      CHARACTER(len=*),INTENT(in),optional     :: sub_name, error_name

      if( bool_condition ) then
        write(0, '("-------------- ERROR ------------------------")')
        if (present(sub_name)) then
          write(0,*) TRIM(sub_name),"    :   "
        endif
        if (present(error_name)) then
  	  write(0,*) TRIM(error_name)
        end if
      STOP
      endif
      RETURN

  END SUBROUTINE FATAL

  SUBROUTINE CHECKPOINT(D_FLAG,opt_string,n_eigs,matrix_1,omega,mass,L,N,dx,t_step,egvs,GS,Time_Evol,norms)
  !This subroutines print variables on console to debug the code
  !For future applications more complex functions can be added

    implicit none
    LOGICAL,INTENT(in)                                        :: D_FLAG
    CHARACTER(len=*),INTENT(in),optional                      :: opt_string
    INTEGER*8,INTENT(in),optional                             :: n_eigs
    REAL*8,INTENT(in),optional                                :: omega,mass,dx,L,N,t_step
    REAL*8,ALLOCATABLE,INTENT(inout),optional                 :: matrix_1(:,:)
    REAL*8,ALLOCATABLE,INTENT(inout),optional                 :: egvs(:),norms(:)
    DOUBLE complex,allocatable,INTENT(inout),optional              :: GS(:)
    DOUBLE complex,allocatable,INTENT(inout),optional              :: Time_Evol(:,:)

    IF (D_FLAG) THEN
      print*,"-------------STARTING CHECKPOINT---------------"
      IF(PRESENT(opt_string)) THEN
       print*, "Called in"
       print*, opt_string
      END IF
      IF(PRESENT(n_eigs)) THEN
       print*, "n_eigs"
       print*, n_eigs
      END IF
      IF(PRESENT(matrix_1)) THEN
       print*, "matrix_1"
       print*, matrix_1
      END IF
      IF(PRESENT(omega)) THEN
       print*, "omega"
       print*, omega
      END IF
      IF(PRESENT(mass)) THEN
       print*, "mass"
       print*, mass
      END IF
      IF(PRESENT(N)) THEN
       print*, "N"
       print*, N
      END IF
      IF(PRESENT(dx)) THEN
       print*, "dx"
       print*, dx
      END IF
      IF(PRESENT(t_step)) THEN
       print*, "t_step"
       print*, t_step
      END IF
      IF(PRESENT(L)) THEN
       print*, "L"
       print*, L
      END IF
      IF(PRESENT(egvs)) THEN
       print*, "egvs"
       print*, egvs
      END IF
      IF(PRESENT(GS)) THEN
       print*, "GS"
       print*, gs
      END IF
      IF(PRESENT(Time_Evol)) THEN
       print*, "Time_Evol"
       print*, Time_Evol
      END IF
      IF(PRESENT(norms)) THEN
       print*, "norms"
       print*, norms
      END IF
      IF(PRESENT(opt_string)) THEN
       print*, "Called in"
       print*, opt_string
      END IF
      print*,"-------------ENDING CHECKPOINT---------------"
    END IF

  END SUBROUTINE
END MODULE


PROGRAM Exercise
! The main program computes the time evolution of the ground state of an harmonic oscillator with time dependent potential
  use DEBUG
  USE,INTRINSIC :: iso_c_binding
  implicit none
  include "FFTW/api/fftw3.f03"
  INTEGER*8                         :: i,j,n_eig,n_T,LWORK,INFO,local_dist,plan
  REAL*8                            :: omega,mass,dx,L,N,T,t_step,temp_norm,PI = 3.141592653589793238462643383279502
  INTEGER*8,DIMENSION(2)            :: shape_holder
  LOGICAL                           :: D_FLAG = .FALSE.
  REAL*8,DIMENSION(:,:),allocatable :: matrix_1
  REAL*8,DIMENSION(:),allocatable   :: egvals,errs,norms
  REAL*8,DIMENSION(:),allocatable   :: WORK
  character(len=100)                :: call_flags,input_file,output_rad
  DOUBLE complex,DIMENSION(:),allocatable:: GS,temp_mom_space
  DOUBLE complex,DIMENSION(:,:),allocatable:: Time_Evol

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialization and pre-checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL GET_COMMAND_ARGUMENT(1,input_file)
  open (unit = 5, file = input_file)
  CALL FATAL(  LEN_TRIM(input_file) == 0,      "main","Missing input file name" )

  CALL GET_COMMAND_ARGUMENT(2,output_rad)
  CALL FATAL(  LEN_TRIM(output_rad) == 0,      "main","Missing output file name" )
  open (unit = 10,Access = 'append', file = TRIM(output_rad)//"_egvals.txt")
  open (unit = 11,Access = 'append', file = TRIM(output_rad)//"_egvects.txt")
  open (unit = 12,Access = 'append', file = TRIM(output_rad)//"_err.txt")

  CALL GET_COMMAND_ARGUMENT(3,call_flags)
  IF (call_flags == "DEBUG")     D_FLAG = .TRUE.

  read (5,*) n_eig,mass,omega,L,n_T,T

  CALL FATAL(   n_eig  <= 0 ,      "main","Number of eigenvalues must be positive" )
  CALL FATAL(   MOD(n_eig,2)/= 0 , "main","Number of eigenvalues must be even" )
  CALL FATAL(   mass   <= 0 ,      "main","Mass must be positive" )
  CALL FATAL(   L      <= 0 ,      "main","Length must be positive" )
  CALL FATAL(   n_T    <= 0 ,      "main","Number of time intervals must be positive" )
  CALL FATAL(   T      <= 0 ,      "main","T must be positive" )

  N = n_eig
  n_eig = n_eig +1
  dx = 1.0*2*L/N
  t_step = 1.0*T/(n_T-1)

  LWORK = 3*n_eig-1
  ALLOCATE(matrix_1(n_eig,n_eig))
  ALLOCATE(egvals(n_eig))
  ALLOCATE(GS(n_eig))
  ALLOCATE(temp_mom_space(n_eig))
  ALLOCATE(Time_Evol(n_eig,n_T))
  ALLOCATE(norms(n_T))
  ALLOCATE(errs(n_eig))
  ALLOCATE(WORK(LWORK))

  CALL CHECKPOINT(D_FLAG,"After Allocation", n_eig,matrix_1,omega,mass,L,N,dx,t_step)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fill Hermitian matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1,n_eig
    matrix_1(i,i) = 2.0/dx/dx + (omega*(-1.0*L + (i-1)*dx))**2
    IF (i < n_eig) THEN
       matrix_1(i,i+1) = -1.0/dx/dx
       matrix_1(i+1,i) = -1.0/dx/dx
    END IF
  END DO
  CALL CHECKPOINT(D_FLAG,"After Initialization", n_eig,matrix_1,omega,mass,L,N,dx,t_step)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute eigenvalues and spacings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL DSYEV('V','L',n_eig,matrix_1,n_eig,egvals,WORK,LWORK,INFO)
  !CALL FATAL(INFO /= 0, "main", "Eigenvalues not found")

  write(10,'(E15.5)')egvals
  !write(11,'(E15.5)')matrix_1
  CALL CHECKPOINT(D_FLAG,"After Computation", n_eig,matrix_1,omega,mass,L,N,dx,t_step,egvals)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute Time evolution of groud state
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
temp_norm = 0

DO  j = 1,n_eig
  temp_norm = temp_norm + (matrix_1(j,1)**2)*dx
END DO
GS = cmplx(matrix_1(:,1)/SQRT(temp_norm))
Time_Evol(:,1) = GS


CALL CHECKPOINT(D_FLAG,"After GS", n_eig,matrix_1,omega,mass,L,N,dx,t_step,egvals,GS,Time_Evol)

DO  i = 2,n_T
  CALL dfftw_plan_dft_1d(plan,n_eig,Time_Evol(:,i),temp_mom_space,FFTW_FORWARD,FFTW_ESTIMATE)
  !Evolve in real space
  DO  j = 1,n_eig
    Time_Evol(j,i) = Time_Evol(j,i-1)* cexp( cmplx( 0.0 , -0.25*(omega**2*(-1.0*L + (j-1)*dx - t_step*(i-1)/T)**2)*t_step) )
  END DO

  temp_mom_space = 0

  !pass to momentum space
  CALL dfftw_execute_dft(plan, Time_Evol(:,i), temp_mom_space)
  CALL dfftw_destroy_plan(plan)
  temp_mom_space = 1.0*temp_mom_space/n_eig

  !evolve in momentum space
  DO  j = 1,n_eig
    IF (j <= 1.0*n_eig/2) THEN
      temp_mom_space(j) = temp_mom_space(j)* cexp( cmplx( 0.0 , -2.0*((PI*j/2.0/L)**2)/mass*t_step) )
    ELSE
      temp_mom_space(j) = temp_mom_space(j)* cexp( cmplx( 0.0 , -2.0*((PI*(j- 1.0*n_eig)/2.0/L)**2)/mass*t_step) )
    END IF
  END DO

  !back to real space
  CALL dfftw_plan_dft_1d(plan,n_eig,temp_mom_space,Time_Evol(:,i),FFTW_BACKWARD,FFTW_ESTIMATE)
  CALL dfftw_execute_dft(plan, temp_mom_space,Time_Evol(:,i))
  CALL dfftw_destroy_plan(plan)
  !Time_Evol(:,i) = Time_Evol(:,i)/ sqrt(DBLE(n_eig))

  DO  j = 1,n_eig
    Time_Evol(j,i) = Time_Evol(j,i)* cexp( cmplx( 0.0 , -0.25*(omega**2*(-1.0*L + (j-1)*dx - t_step*(i-1)/T)**2)*t_step) )
  END DO

END DO

DO  i = 1,n_T
  temp_norm = 0
  DO  j = 1,n_eig
    temp_norm = temp_norm + (Time_Evol(j,i)*CONJG(Time_Evol(j,i)))*dx
  END DO
  norms(i) = SQRT(temp_norm)
  Time_Evol(:,i) = Time_Evol(:,i)/norms(i)
END DO


CALL CHECKPOINT(D_FLAG,"After Evolution", n_eig,matrix_1,omega,mass,L,N,dx,t_step,egvals,GS,Time_Evol,norms)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute Error in respect to expected values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i = 1,n_eig
    errs(i) = ((i-0.5)*2*omega - egvals(i))/((i-0.5)*2*omega )
  END DO

  !write(12,'(E15.5)') ABS(errs)

  DO  i = 1,n_T
    DO j = 1,n_eig
       Time_Evol(j,i) = Time_Evol(j,i)*CONJG(Time_Evol(j,i))
    END do
  END DO
  DO i = 1,n_eig
          write(11,*)REAL(Time_Evol(i,:))
  END DO
  write(12,*) n_eig,mass,omega,L,SUM(ABS(errs(1:10)))
  CALL CHECKPOINT(D_FLAG,"Final", n_eig,matrix_1,omega,mass,L,N,dx,t_step,egvals,GS,Time_Evol,norms)



END PROGRAM Exercise
