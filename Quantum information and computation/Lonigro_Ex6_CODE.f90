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

  SUBROUTINE CHECKPOINT(D_FLAG,opt_string,n_eigs,matrix_1,omega,mass,L,N,dx,egvs)
  !This subroutines print variables on console to debug the code
  !For future applications more complex functions can be added

    implicit none
    LOGICAL,INTENT(in)                                        :: D_FLAG
    CHARACTER(len=*),INTENT(in),optional                      :: opt_string
    INTEGER*8,INTENT(in),optional                             :: n_eigs
    REAL*8,INTENT(in),optional                                  :: omega,mass,dx,L,N
    REAL*8,ALLOCATABLE,INTENT(inout),optional                :: matrix_1(:,:)
    REAL*8,ALLOCATABLE,INTENT(inout),optional                   :: egvs(:)

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
      IF(PRESENT(L)) THEN
       print*, "L"
       print*, L
      END IF
      IF(PRESENT(egvs)) THEN
       print*, "egvs"
       print*, egvs
      END IF
      print*,"-------------ENDING CHECKPOINT---------------"
    END IF

  END SUBROUTINE
END MODULE


PROGRAM Exercise
! The main program checks the different time needed to perform a matrix multiplication performing explicitely
! the product in different order and using the native function MATMUL in fortran
  use DEBUG
  implicit none
  INTEGER*8                         :: i,j,n_eig,LWORK,INFO,local_dist
  REAL*8                              :: omega,mass,dx,L,N
  INTEGER*8,DIMENSION(2)            :: shape_holder
  LOGICAL                           :: D_FLAG = .FALSE.
  REAL*8,DIMENSION(:,:),allocatable:: matrix_1
  REAL*8,DIMENSION(:),allocatable     :: egvals,errs
  REAL*8,DIMENSION(:),allocatable  :: WORK
  character(len=100)                :: call_flags,input_file,output_rad

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

  read (5,*) n_eig,mass,omega,L

  CALL FATAL(   n_eig  <= 0 ,      "main","Number of eigenvalues must be positive" )
  CALL FATAL(   mass   <= 0 ,      "main","Mass must be positive" )
  CALL FATAL(   L      <= 0 ,      "main","Length must be positive" )

  N = n_eig
  n_eig = n_eig +1
  dx = 1.0*2*L/N

  LWORK = 3*n_eig-1
  ALLOCATE(matrix_1(n_eig,n_eig))
  ALLOCATE(egvals(n_eig))
  ALLOCATE(errs(n_eig))
  ALLOCATE(WORK(LWORK))

  CALL CHECKPOINT(D_FLAG,"After Allocation", n_eig,matrix_1,omega,mass,L,N,dx)
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
  CALL CHECKPOINT(D_FLAG,"After Initialization", n_eig,matrix_1,omega,mass,L,N,dx)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute eigenvalues and spacings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL DSYEV('V','L',n_eig,matrix_1,n_eig,egvals,WORK,LWORK,INFO)
  CALL FATAL(INFO /= 0, "main", "Eigenvalues not found")

  write(10,'(E15.5)')egvals
  write(11,'(E15.5)')matrix_1
  CALL CHECKPOINT(D_FLAG,"After Computation", n_eig,matrix_1,omega,mass,L,N,dx,egvals)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute Error in respect to expected values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  DO i = 1,n_eig
    errs(i) = ((i-0.5)*2*omega - egvals(i))/((i-0.5)*2*omega )
  END DO

  !write(12,'(E15.5)') ABS(errs)
  write(12,*) n_eig,mass,omega,L,SUM(ABS(errs(1:10)))
  CALL CHECKPOINT(D_FLAG,"Final", n_eig,matrix_1,omega,mass,L,N,dx,egvals)
END PROGRAM Exercise
