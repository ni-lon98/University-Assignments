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

  SUBROUTINE CHECKPOINT(D_FLAG,opt_string,N,H,coupling_1,coupling_2,egvs)
  !This subroutines print variables on console to debug the code
  !For future applications more complex functions can be added

    implicit none
    LOGICAL,INTENT(in)                                        :: D_FLAG
    CHARACTER(len=*),INTENT(in),optional                      :: opt_string
    INTEGER,INTENT(in),optional                             :: N
    REAL,INTENT(in),optional                                :: coupling_1,coupling_2
    DOUBLE COMPLEX,ALLOCATABLE,INTENT(inout),optional       :: H(:,:)
    REAL*8,ALLOCATABLE,INTENT(inout),optional       :: egvs(:)

    IF (D_FLAG) THEN
      print*,"-------------STARTING CHECKPOINT---------------"
      IF(PRESENT(opt_string)) THEN
       print*, "Called in"
       print*, opt_string
      END IF
      IF(PRESENT(N)) THEN
       print*, "N"
       print*, N
      END IF
      IF(PRESENT(H)) THEN
       print*, "H"
       print*, H
      END IF
      IF(PRESENT(coupling_1)) THEN
       print*, "coupling_1"
       print*, coupling_1
      END IF
      IF(PRESENT(coupling_2)) THEN
       print*, "coupling_2"
       print*, coupling_2
      END IF
      IF(PRESENT(egvs)) THEN
       print*, "egvs"
       print*, egvs
      END IF
      print*,"-------------ENDING CHECKPOINT---------------"
    END IF

  END SUBROUTINE
END MODULE

MODULE HAMILTONIAN
  contains
  SUBROUTINE Def_Hamiltonian(N,H,coupling_1,coupling_2)
    implicit none
    INTEGER,INTENT(in)                          :: N
    REAL,INTENT(in),optional             :: coupling_1,coupling_2
    DOUBLE COMPLEX,ALLOCATABLE,INTENT(inout)       :: H(:,:)
    DOUBLE COMPLEX,ALLOCATABLE                     :: int_1(:,:),int_2(:,:)
    INTEGER                                        :: i,j,k,m
    INTEGER                                        :: swap

    ALLOCATE(int_1(2**N,2**N))
    ALLOCATE(int_2(2**N,2**N))
    !FIRST ORDER COUPLING
    IF (ABS(coupling_1) > EPSILON (coupling_1)) THEN
      H = 0
      DO k=1,N
        swap = -1
        DO i = 1,2**N
            IF (MOD(i-1,INT(2**N/2**K)) == 0) swap = -swap
            H(i,i) = H(i,i) + coupling_1*swap
        END DO
      END DO
    END IF
    !SECOND ORDER COUPLING
    IF(ABS(coupling_2) > EPSILON (coupling_2)) THEN
      DO k=1,N-1
        int_1 = 0
        int_2 = 0
        DO m=1,INT(2**N/2**(k))
          DO i = 1,2**k
            IF (i > 2**(k-1) ) THEN
              int_1(i+(m-1)*2**k,i -2**(k-1) +(m-1)*2**k) = int_1(i+(m-1)*2**k,i-2**(k-1)+(m-1)*2**k) + 1
            ELSE
              int_1(i+(m-1)*2**k,i +2**(k-1) +(m-1)*2**k) = int_1(i+(m-1)*2**k,i+2**(k-1)+(m-1)*2**k) + 1
            END IF
          END DO
        END DO
        DO m=1,INT(2**N/2**(k+1))
          DO i = 1,2**(k+1)
            IF (i > 2**(k) ) THEN
              int_2(i+(m-1)*2**(k+1),i -2**(k) +(m-1)*2**(k+1)) = int_2(i+(m-1)*2**(k+1),i -2**(k)+(m-1)*2**(k+1)) + 1
            ELSE
              int_2(i+(m-1)*2**(k+1),i +2**(k) +(m-1)*2**(k+1)) = int_2(i+(m-1)*2**(k+1),i+2**(k)+(m-1)*2**(k+1)) + 1
            END IF
          END DO
        END DO

      !print *, int_2
    !    print *, int_1
  !    print*,H
      H = H + coupling_2*MATMUL(int_2,int_1)
      END DO

    END IF

    DEALLOCATE(int_1)
    DEALLOCATE(int_2)

    !COULD ADD MORE COUPLING!

    RETURN
  END SUBROUTINE

  SUBROUTINE RG_iter(N,H,coupling_2,N_steps)
    use DEBUG
    implicit none
    INTEGER,INTENT(in)                          :: N,N_steps
    REAL,INTENT(inout),optional                    :: coupling_2
    DOUBLE COMPLEX,ALLOCATABLE,INTENT(inout)    :: H(:,:)
    DOUBLE COMPLEX,ALLOCATABLE                  :: H_copy(:,:),H_tot(:,:),int_1(:,:),int_2(:,:),P(:,:),P_dag(:,:), &
                                                    int_1_tot(:,:),int_2_tot(:,:)
    INTEGER                                     :: i,j,k,m,step,LWORK,INFO,dim
    REAL*8,DIMENSION(:),allocatable             :: egvals
    DOUBLE COMPLEX,DIMENSION(:),allocatable     :: WORK,RWORK

    ALLOCATE(H_tot(2**(2*N),2**(2*N)))
    ALLOCATE(H_copy(2**(2*N),2**(2*N)))
    ALLOCATE(int_1(2**N,2**N))
    ALLOCATE(int_2(2**N,2**N))
    ALLOCATE(int_1_tot(2**(2*N),2**(2*N)))
    ALLOCATE(int_2_tot(2**(2*N),2**(2*N)))
    ALLOCATE(P(N**2,N))
    ALLOCATE(P_dag(N**2,N))
    LWORK = 2*(2**(2*N))
    dim = 2**(2*N)
    ALLOCATE(egvals(2**(2*N)))
    ALLOCATE(WORK(LWORK))
    ALLOCATE(RWORK( 3*(2**(2*N))))
    H_tot = 0
    int_1 = 0
    int_2 = 0

    !Define Initial interaction terms
    DO i = 1,2**N
        IF (i > 2**(N-1) ) THEN
            int_1(i,i-2**(N-1)) = int_1(i,i-2**(N-1)) + 1
        ELSE
            int_1(i,i +2**(N-1)) = int_1(i,i+2**(N-1)) + 1
        END IF
    END DO

    DO m=1,2**(N-1)
      DO i = 1,2
        IF (i > 1 ) THEN
          int_2(i+(m-1)*2,i -1 +(m-1)*2) = int_2(i+(m-1)*2,i-1+(m-1)*2) + 1
        ELSE
          int_2(i+(m-1)*2,i +1 +(m-1)*2) = int_2(i+(m-1)*2,i+1+(m-1)*2) + 1
        END IF
      END DO
    END DO

    !RG iteration loop
    DO step = 1,N_steps
      H_tot = 0
      int_1_tot = 0
      int_2_tot = 0

      DO i = 1,2**N
        DO j = 1,2**N
          DO k = 1,2**N
            H_tot((i-1)*2**N+j,(i-1)*2**N+k) = H_tot((i-1)*2**N+j,(i-1)*2**N+k) + H(j,k)  ! Tensor Product of H and ID
            int_1_tot((i-1)*2**N+j,(i-1)*2**N+k) = int_1_tot((i-1)*2**N+j,(i-1)*2**N+k) + int_1(j,k)
          END DO
        END DO
      END DO

      DO i = 1,2**N
        DO j = 1,2**N
          DO k = 1,2**N
            H_tot((i-1)*2**N+k, (j-1)*2**N+k) = H_tot((i-1)*2**N+k, (j-1)*2**N+k) + H(i,j)  ! Tensor Product of ID and H
            int_2_tot((i-1)*2**N+j,(i-1)*2**N+k) = int_2_tot((i-1)*2**N+j,(i-1)*2**N+k) + int_2(j,k)
          END DO
        END DO
      END DO


      DO i=1, 2**N
        DO j=1, 2**N
          DO k=1, 2**N
            DO m=1, 2**N
              H_tot( (i-1)*2**N+k,(j-1)*2**N+m) =   H_tot((i-1)*2**N+k,(j-1)*2**N+m) + coupling_2*int_1(i,j)*int_2(k,m) !Interaction term
            END DO
          END DO
        END DO
      END DO

      H_copy = H_tot
      CALL zheev('V','U',dim,H_tot,dim,egvals,WORK,LWORK,RWORK,INFO)
      CALL FATAL(INFO /= 0, "main", "Eigenvalues not found")

      P = H_tot(:,1:2**(2*N-1))
      P_dag = TRANSPOSE(CONJG(P))
      coupling_2 = coupling_2*0.5
      H = 0.5*MATMUL(MATMUL(P_dag,H_copy),P)
      int_1 = MATMUL(P_dag,MATMUL(int_1_tot,P))
      int_2 = MATMUL(P_dag,MATMUL(int_2_tot,P))
    END DO

  END SUBROUTINE
END MODULE

PROGRAM Exercise
! The main program checks the different time needed to perform a matrix multiplication performing explicitely
! the product in different order and using the native function MATMUL in fortran
  use DEBUG
  use HAMILTONIAN
  implicit none
  INTEGER*8                      ::dim
  INTEGER                        :: N,i,j,LWORK,INFO,ierror,levels,N_steps
  REAL                            :: couple_1,couple_2
  LOGICAL                           :: D_FLAG = .FALSE.
  DOUBLE COMPLEX,DIMENSION(:,:),allocatable :: H
  REAL*8,DIMENSION(:),allocatable   :: egvals
  DOUBLE COMPLEX,DIMENSION(:),allocatable   :: WORK,RWORK
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
  !open (unit = 11,Access = 'append', file = TRIM(output_rad)//"_egvects.txt")

  CALL GET_COMMAND_ARGUMENT(3,call_flags)
  IF (call_flags == "DEBUG")     D_FLAG = .TRUE.


  read(5,'(I10)',iostat=ierror) N,levels,N_steps
  CALL FATAL( ierror > 0 ,      "main","Inputs are not integers" )

  read (5,*) couple_1,couple_2
  CALL FATAL(   N             <= 0 ,      "main","Number of q-bits" )
  CALL FATAL(   levels        <= 0 ,      "main","Number of q-bits" )
  CALL FATAL(   N_steps       <= 0 ,      "main","Number of q-bits" )

  LWORK = 2*(2**N)
  ALLOCATE(egvals(2**N))
  ALLOCATE(WORK(LWORK))
  ALLOCATE(RWORK( 3*(2**N)))
  ALLOCATE(H(2**N,2**N))

  CALL CHECKPOINT(D_FLAG,"After Allocation",N,H,couple_1,couple_2,egvals)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Fill Hermitian matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL Def_Hamiltonian(N,H,couple_1,couple_2)
  CALL CHECKPOINT(D_FLAG,"After Initialization",N,H,couple_1,couple_2,egvals)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Renormalization group loop
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  CALL RG_iter(N,H,couple_2,N_steps)
  CALL CHECKPOINT(D_FLAG,"After Initialization",N,H,couple_1,couple_2,egvals)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Compute final eigenvalues and spacings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  dim = 2**N
  CALL zheev('N','U',dim,H,dim,egvals,WORK,LWORK,RWORK,INFO)
  CALL FATAL(INFO /= 0, "main", "Eigenvalues not found")

  write(10,*)couple_1,couple_2,egvals(:levels)/N
  !write(11,'(E15.5)')H
  CALL CHECKPOINT(D_FLAG,"After Computation",N,H,couple_1,couple_2,egvals)

END PROGRAM Exercise
