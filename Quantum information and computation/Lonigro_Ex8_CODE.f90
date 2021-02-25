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

  SUBROUTINE CHECKPOINT(D_FLAG,opt_string,D,N, wavefunction , N_body,density,subden_1,subden_2)
  !This subroutines print variables on console to debug the code
  !For future applications more complex functions can be added

    implicit none
    LOGICAL,INTENT(in)                                                  :: D_FLAG
    CHARACTER(len=*),INTENT(in),optional                                :: opt_string
    INTEGER*16,INTENT(in),optional                                       :: D,N
    DOUBLE COMPLEX,ALLOCATABLE,INTENT(inout),optional                   :: wavefunction(:)
    DOUBLE COMPLEX,ALLOCATABLE,INTENT(inout),optional                   :: N_body(:)
    DOUBLE COMPLEX,allocatable,INTENT(inout),optional                   :: density(:,:),subden_1(:,:),subden_2(:,:)

    IF (D_FLAG) THEN
      print*,"-------------STARTING CHECKPOINT---------------"
      IF(PRESENT(opt_string)) THEN
       print*, "Called in"
       print*, opt_string
      END IF
      IF(PRESENT(D)) THEN
       print*, "D"
       print*, D
      END IF
      IF(PRESENT(N)) THEN
       print*, "N"
       print*, N
      END IF
      IF(PRESENT(wavefunction)) THEN
       print*, "wavefunction"
       print*, wavefunction
      END IF
      IF(PRESENT(N_body)) THEN
       print*, "N_body"
       print*, N_body
      END IF
      IF(PRESENT(density)) THEN
       print*, "density"
       print*, density
      END IF
      IF(PRESENT(subden_1)) THEN
       print*, "subden_1"
       print*, subden_1
      END IF
      IF(PRESENT(subden_2)) THEN
       print*, "subden_2"
       print*, subden_2
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
  INTEGER*16                                   :: D,N,i,j,k,ierror,ierror_2
  LOGICAL                                     :: D_FLAG = .FALSE.,Time_flag = .FALSE.,Generic_flag = .FALSE.
  DOUBLE COMPLEX,DIMENSION(:),allocatable     :: wavefunction
  DOUBLE COMPLEX,DIMENSION(:),allocatable     :: N_body
  DOUBLE COMPLEX,DIMENSION(:,:),allocatable   :: density,subden_1,subden_2
  character(len=100)                          :: call_flags,input_file,output_rad
  REAL*16                                     :: start, finish,start_2,finish_2,norm_wave,norm_n_body
  DOUBLE COMPLEX                              :: trace

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialization and pre-checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  CALL GET_COMMAND_ARGUMENT(1,input_file)
  open (unit = 5, file = input_file)
  CALL FATAL(  LEN_TRIM(input_file) == 0,      "main","Missing input file name" )

  CALL GET_COMMAND_ARGUMENT(2,output_rad)
  CALL FATAL(  LEN_TRIM(output_rad) == 0,      "main","Missing output file name" )
  open (unit = 10,Access = 'append', file = TRIM(output_rad)//".txt")

  CALL GET_COMMAND_ARGUMENT(3,call_flags)
  IF (call_flags == "DEBUG")     D_FLAG = .TRUE.
  IF (call_flags == "TIME")     Time_flag = .TRUE.
  IF (call_flags == "GENERIC")     Generic_flag = .TRUE.
  CALL GET_COMMAND_ARGUMENT(4,call_flags)
  IF (call_flags == "DEBUG")     D_FLAG = .TRUE.
  IF (call_flags == "TIME")     Time_flag = .TRUE.
  IF (call_flags == "GENERIC")     Generic_flag = .TRUE.


  IF (Time_flag .EQV. .FALSE.)  Generic_flag = .TRUE.
  read(5,'(I10)',iostat=ierror) D,N
  CALL FATAL( ierror > 0 ,      "main","Inputs are not integers" )

  CALL FATAL(  D <= 0 ,      "main","Dimension must be non-negative" )
  CALL FATAL(  N <= 0 ,      "main","Number of elements must be non-negative" )
  !read(*,'(2I10)',iostat=ierror) intval


  ALLOCATE(wavefunction(D*N))
  IF (Generic_flag) ALLOCATE(N_body(D**N))
  IF (Time_flag .EQV. .FALSE.) THEN
    ALLOCATE(density(D**2,D**2))
    ALLOCATE(subden_1(D,D))
    ALLOCATE(subden_2(D,D))
  END IF
  IF (Time_flag .EQV. .FALSE.) CALL CHECKPOINT(D_FLAG,"After Allocation", D,N,wavefunction,N_body,density)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fill with random numbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  IF (Time_flag)call cpu_time(start)
  DO i = 1,D*N
      wavefunction(i) =  COMPLEX(RAND(),RAND())
  END DO
  norm_wave = 0
  DO i = 1,D*N
      norm_wave = norm_wave + wavefunction(i) *CONJG(wavefunction(i))
  END DO
  wavefunction = wavefunction/SQRT(norm_wave)

  IF (Time_flag) call cpu_time (finish)
  !write(10,*)D,"  ", finish-start

  IF (Time_flag) call cpu_time(start_2)

  IF (Generic_flag) THEN
  DO i = 1,D**N
      N_body(i) = COMPLEX(RAND(),RAND())
  END DO
    norm_n_body = 0
    DO i = 1,D*N
      norm_n_body = norm_n_body + N_body(i) *CONJG(N_body(i))
    END DO
    N_body = N_body/SQRT(norm_n_body)
  ENDIF


  IF (Time_flag) call cpu_time (finish_2)
  IF (Time_flag) write(10,*)D,N, finish-start, finish_2 - start_2

  IF (Time_flag .EQV. .FALSE.) THEN
    DO i = 1,D**2
      DO j = 1,D**2
        density(i,j) =  N_body(i)*CONJG(N_body(j))
      END DO
    END DO
  END IF
  IF (Time_flag .EQV. .FALSE.) CALL CHECKPOINT(D_FLAG,"After Allocation", D,N,wavefunction,N_body,density,subden_1,subden_2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Partial trace
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF (Time_flag .EQV. .FALSE.) THEN
    DO i = 1,D
      DO j = 1,D
        subden_2(i,j) = 0
        DO k = 1,D
          subden_2(i,j) =  subden_2(i,j) + density(i+D*(k-1),j+D*(k-1))
        END DO
      END DO
    END DO
    DO i = 1,D
      DO j = 1,D
        subden_1(i,j) = 0
        DO k = 1,D
          subden_1(i,j) =  subden_1(i,j) + density(k+D*(i-1),k+D*(j-1))
        END DO
      END DO
    END DO
END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Post checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
norm_wave = 0
DO i = 1,D*N
    norm_wave = norm_wave + wavefunction(i) *CONJG(wavefunction(i))
END DO
norm_n_body = 0
DO i = 1,D*N
  norm_n_body = norm_n_body + N_body(i) *CONJG(N_body(i))
END DO

IF (Time_flag .EQV. .FALSE.) THEN
  write (10,*) "Separable state  ",wavefunction
  write (10,*) "Norm of separable state  ",SQRT(norm_wave)
  write (10,*) "Generic state  ",N_body
  write (10,*) "Norm of generic state  ",SQRT(norm_n_body)
  trace = 0
  DO i = 1,D**N
    trace = trace + density(i,i)
  END DO
  write (10,*) "Trace of density  ",trace
  write (10,*) "Density matrix    ",density
  trace = 0
  DO i = 1,D
    trace = trace + subden_1(i,i)
  END DO
  write (10,*) "Trace of density matrix of first system  ",trace
  write (10,*) "Density matrix of first system",subden_1
  trace = 0
  DO i = 1,D
    trace = trace + subden_2(i,i)
  END DO
  write (10,*) "Trace of density matrix of second system  ",trace
  write (10,*) "Density matrix of second system",subden_2
END IF


IF (Time_flag .EQV. .FALSE.) CALL CHECKPOINT(D_FLAG,"After Allocation", D,N,wavefunction,N_body,density,subden_1,subden_2)
END PROGRAM Exercise
