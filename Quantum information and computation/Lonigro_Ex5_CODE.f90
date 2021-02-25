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

  SUBROUTINE CHECKPOINT(D_FLAG,opt_string,dim_x_1,dim_y_1,matrix_1,egvs,spacings)
  !This subroutines print variables on console to debug the code
  !For future applications more complex functions can be added

    implicit none
    LOGICAL,INTENT(in)                                        :: D_FLAG
    CHARACTER(len=*),INTENT(in),optional                      :: opt_string
    INTEGER*8,INTENT(in),optional                             :: dim_x_1,dim_y_1
    COMPLEX,ALLOCATABLE,INTENT(inout),optional                   :: matrix_1(:,:)
    REAL,ALLOCATABLE,INTENT(inout),optional                   :: egvs(:),spacings(:)

    IF (D_FLAG) THEN
      print*,"-------------STARTING CHECKPOINT---------------"
      IF(PRESENT(opt_string)) THEN
       print*, "Called in"
       print*, opt_string
      END IF
      IF(PRESENT(dim_x_1)) THEN
       print*, "dim_x_1"
       print*, dim_x_1
      END IF
      IF(PRESENT(dim_y_1)) THEN
       print*, "dim_y_1"
       print*, dim_y_1
      END IF
      IF(PRESENT(matrix_1)) THEN
       print*, "matrix_1"
       print*, matrix_1
      END IF
      IF(PRESENT(egvs)) THEN
       print*, "egvs"
       print*, egvs
      END IF
      IF(PRESENT(spacings)) THEN
       print*, "spacings"
       print*, spacings
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
  INTEGER*8                         :: i,j,dim_x_1,dim_y_1,LWORK,INFO,local_dist
  INTEGER*8,DIMENSION(2)            :: shape_holder
  LOGICAL                           :: D_FLAG = .FALSE., diag_switch = .FALSE.
  COMPLEX,DIMENSION(:,:),allocatable   :: matrix_1
  REAL,DIMENSION(:),allocatable     :: egvs,RWORK,spacings, mean_spacing,ratios
  COMPLEX,DIMENSION(:),allocatable  :: WORK
  character(len=100)                :: call_flags,input_file,output_rad,local_char

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Initialization and pre-checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL GET_COMMAND_ARGUMENT(1,input_file)
  open (unit = 5, file = input_file)
  CALL FATAL(  LEN_TRIM(input_file) == 0,      "main","Missing input file name" )

  CALL GET_COMMAND_ARGUMENT(2,output_rad)
  CALL FATAL(  LEN_TRIM(output_rad) == 0,      "main","Missing output file name" )
  open (unit = 10,Access = 'append', file = TRIM(output_rad)//".txt")
  open (unit = 11,Access = 'append', file = TRIM(output_rad)//"_rat.txt")

  CALL GET_COMMAND_ARGUMENT(3,local_char)
  READ(local_char,*)local_dist
  CALL FATAL(  local_dist<= 0,      "main","Local distance must be positive integer " )

  CALL GET_COMMAND_ARGUMENT(4,call_flags)
  IF (call_flags == "DEBUG")     D_FLAG = .TRUE.
  IF (call_flags == "DIAGONAL")  diag_switch = .TRUE.

  CALL GET_COMMAND_ARGUMENT(5,call_flags)
  IF (call_flags == "DEBUG")     D_FLAG = .TRUE.
  IF (call_flags == "DIAGONAL")  diag_switch = .TRUE.

  read (5,*) dim_x_1
  dim_y_1 = dim_x_1

  CALL FATAL(  MOD(dim_x_1,local_dist) /= 0 ,      "main","Dimension must be multiple of local_dist" )
  CALL FATAL(  dim_x_1 <= 0 ,      "main","Dimension of first matrix must be non-negative" )
  CALL FATAL(  dim_y_1 <= 0 ,      "main","Dimension of first matrix must be non-negative" )

  LWORK = 2*dim_x_1-1
  ALLOCATE(matrix_1(dim_x_1,dim_y_1))
  ALLOCATE(egvs(dim_x_1))
  ALLOCATE(spacings(dim_x_1-2))
  ALLOCATE(ratios(dim_x_1-1))
  ALLOCATE(mean_spacing(dim_x_1-2))
  ALLOCATE(WORK(LWORK))
  ALLOCATE(RWORK(dim_x_1*3-2))

  CALL CHECKPOINT(D_FLAG,"After Allocation", dim_x_1,dim_y_1,matrix_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fill with random numbers
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  DO i = 1,dim_x_1
    matrix_1(i,i) = RAND()
  END DO
  IF(diag_switch .EQV. .FALSE.) THEN
    DO i = 2,dim_x_1
      DO j = 1,i-1
        matrix_1(i,j) = COMPLEX(RAND(),RAND())
      END DO
    END DO
  END IF
  CALL CHECKPOINT(D_FLAG,"After Initialization", dim_x_1,dim_y_1,matrix_1)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute eigenvalues and spacings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL CHEEV('N','L',dim_x_1,matrix_1,dim_x_1,egvs,WORK,LWORK,RWORK,INFO)
  CALL FATAL(INFO /= 0, "main", "Eigenvalues not found")
  DO i = 1,dim_x_1-2
    spacings(i) =   egvs(i+1) -  egvs(i+2)
  END DO
  DO i = 1,dim_x_1-3
    ratios(i) =   MIN(spacings(i),spacings(i+1))/MAX(spacings(i),spacings(i+1))
  END DO
  IF(local_dist /= 1) THEN
    DO i = 1,dim_x_1-2
      mean_spacing(i) = SUM(spacings(max(1,i-dim_x_1/local_dist):min(dim_x_1-1,i + dim_x_1/local_dist)))/ &
                    (min(dim_x_1-1,i+dim_x_1/local_dist) - max(1,i-dim_x_1/local_dist))
    END DO
    DO i = 1,dim_x_1-2
      spacings(i) = spacings(i)/mean_spacing(i)
    END DO
  ELSE
    mean_spacing(1) = SUM(spacings)/(dim_x_1-1)
    spacings = spacings/mean_spacing(1)
  END IF
  write(10,'(F15.8)')spacings
  write(11,'(F15.8)')SUM(ratios)/(dim_x_1-2)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Post checks
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  CALL CHECKPOINT(D_FLAG,"After Computation", dim_x_1,dim_y_1,matrix_1,egvs,spacings)
END PROGRAM Exercise
