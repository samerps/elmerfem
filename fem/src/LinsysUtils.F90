!/*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! *  This library is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU Lesser General Public
! *  License as published by the Free Software Foundation; either
! *  version 2.1 of the License, or (at your option) any later version.
! *
! *  This library is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! *  Lesser General Public License for more details.
! * 
! *  You should have received a copy of the GNU Lesser General Public
! *  License along with this library (in file ../LGPL-2.1); if not, write 
! *  to the Free Software Foundation, Inc., 51 Franklin Street, 
! *  Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Utilities for dealing with the linear system (matrix-vector, scaling, loads etc.)
! *
! ******************************************************************************
! *
! *  Authors: Juha Ruokolainen
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 28 Sep 1998
! *
! *****************************************************************************/

!> Basic utilities used by individual solvers. 
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE LinsysUtils

#include "../config.h"

   USE ElementUtils
   USE ParallelUtils
   USE ListMatrix
   USE CRSMatrix
   
   IMPLICIT NONE


CONTAINS


!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE MatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_MatrixVectorMultiply( A,u,v )

     CASE( MATRIX_BAND,MATRIX_SBAND )
       CALL Band_MatrixVectorMultiply( A,u,v )

     CASE( MATRIX_LIST )
       CALL Warn('MatrixVectorMultiply','Not implemented for List matrix type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MatrixVectorMultiply


!------------------------------------------------------------------------------
!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE MaskedMatrixVectorMultiply( A,u,v,ActiveRow,ActiveCol )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
     LOGICAL, DIMENSION(:) :: ActiveRow
     LOGICAL, DIMENSION(:) :: ActiveCol
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_MaskedMatrixVectorMultiply( A,u,v,ActiveRow, ActiveCol )

     CASE DEFAULT
       CALL Fatal('MaskedMatrixVectorMultiply','Not implemented for List matrix type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MaskedMatrixVectorMultiply
!------------------------------------------------------------------------------


!> Matrix vector multiplication of sparse matrices.
!------------------------------------------------------------------------------
   SUBROUTINE TransposeMatrixVectorMultiply( A,u,v )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n
     REAL(KIND=dp), DIMENSION(:) CONTIG :: u,v
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
     CASE( MATRIX_CRS )
       CALL CRS_TransposeMatrixVectorMultiply( A,u,v )

     CASE DEFAULT 
       CALL Fatal('TransposeMatrixVectorMultiply','Not implemented for other than CRS type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE TransposeMatrixVectorMultiply
!------------------------------------------------------------------------------


!> Create a copy of the linear system (Values,Rhs) to (BulkValues,BulkRhs).
!------------------------------------------------------------------------------
   SUBROUTINE CopyBulkMatrix( A, BulkMass, BulkDamp )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,n
     LOGICAL, OPTIONAL :: BulkMass, BulkDamp
     
     n = SIZE( A % Rhs )
     IF( ASSOCIATED( A % BulkRhs ) ) THEN
       IF( SIZE( A % BulkRhs ) /= n ) THEN
          DEALLOCATE( A % BulkRhs ) 
          A % BulkRHS => NULL()
       END IF
     END IF
     IF ( .NOT. ASSOCIATED( A % BulkRHS ) ) THEN
       ALLOCATE( A % BulkRHS( n ) )
     END IF
     DO i=1,n
       A % BulkRHS(i) = A % Rhs(i)
    END DO
     
     n = SIZE( A % Values )
     IF( ASSOCIATED( A % BulkValues ) ) THEN
       IF( SIZE( A % BulkValues ) /= n ) THEN
          DEALLOCATE( A % BulkValues ) 
          A % BulkValues => NULL()
       END IF
     END IF
     IF ( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
       ALLOCATE( A % BulkValues( n ) )
     END IF

     DO i=1,n
       A % BulkValues(i) = A % Values(i)
     END DO

     IF( PRESENT( BulkMass ) .AND. ASSOCIATED( A % MassValues) ) THEN
       IF( BulkMass ) THEN
         n = SIZE( A % MassValues )
         IF( ASSOCIATED( A % BulkMassValues ) ) THEN
           IF( SIZE( A % BulkMassValues ) /= n ) THEN
             DEALLOCATE( A % BulkMassValues ) 
             A % BulkMassValues => NULL()
           END IF
         END IF
         IF ( .NOT. ASSOCIATED( A % BulkMassValues ) ) THEN
           ALLOCATE( A % BulkMassValues( n ) )
         END IF

         DO i=1,n
           A % BulkMassValues(i) = A % MassValues(i)
         END DO
       END IF
     END IF

     IF( PRESENT( BulkDamp ) .AND. ASSOCIATED( A % DampValues) ) THEN
       IF( BulkDamp ) THEN
         n = SIZE( A % DampValues )
         IF( ASSOCIATED( A % BulkDampValues ) ) THEN
           IF( SIZE( A % BulkDampValues ) /= n ) THEN
             DEALLOCATE( A % BulkDampValues ) 
             A % BulkDampValues => NULL()
           END IF
         END IF
         IF ( .NOT. ASSOCIATED( A % BulkDampValues ) ) THEN
           ALLOCATE( A % BulkDampValues( n ) )
         END IF

         DO i=1,n
           A % BulkDampValues(i) = A % DampValues(i)
         END DO
       END IF
     END IF
     
   END SUBROUTINE CopyBulkMatrix
!------------------------------------------------------------------------------


!> Restores the saved bulk after
!------------------------------------------------------------------------------
   SUBROUTINE RestoreBulkMatrix( A )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,n
     
     IF( ASSOCIATED( A % BulkRhs ) ) THEN
       n = SIZE( A % Rhs )
       IF( SIZE( A % BulkRhs ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore rhs of different size!')
       END IF
       A % Rhs(1:n) = A % BulkRhs(1:n)
     END IF
     
     IF( ASSOCIATED( A % BulkValues ) ) THEN
       n = SIZE( A % Values )
       IF( SIZE( A % BulkValues ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore matrix of different size!')
       END IF
       DO i=1,n
         A % Values(i) = A % BulkValues(i)
       END DO
     END IF

     IF( ASSOCIATED( A % BulkMassValues ) ) THEN
       n = SIZE( A % MassValues )
       IF( SIZE( A % BulkMassValues ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore mass matrix of different size!')
       END IF
       DO i=1,n
         A % MassValues(i) = A % BulkMassValues(i)
       END DO
     END IF

     IF( ASSOCIATED( A % BulkDampValues ) ) THEN
       n = SIZE( A % DampValues )
       IF( SIZE( A % BulkDampValues ) /= n ) THEN
         CALL Fatal('RestoreBulkMatrix','Cannot restore damp matrix of different size!')
       END IF
       DO i=1,n
         A % DampValues(i) = A % BulkDampValues(i)
       END DO
     END IF
     
   END SUBROUTINE RestoreBulkMatrix
!------------------------------------------------------------------------------





  
!------------------------------------------------------------------------------
!>  Scale system Ax = b as:
!>  (DAD)y = Db, where D = 1/SQRT(Diag(A)), and y = D^-1 x
!------------------------------------------------------------------------------
  SUBROUTINE ScaleLinearSystem(Solver,A,b,x,DiagScaling, & 
          ApplyScaling,RhsScaling,ConstraintScaling)

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:),x(:)
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)
    LOGICAL, OPTIONAL :: ApplyScaling, RhsScaling,ConstraintScaling
    INTEGER :: n,i,j
    REAL(KIND=dp) :: bnorm,s
    COMPLEX(KIND=dp) :: DiagC
    LOGICAL :: ComplexMatrix, DoRHS, DoCM, Found
    REAL(KIND=dp), POINTER  :: Diag(:)

    TYPE(Matrix_t), POINTER :: CM

    n = A % NumberOfRows
    
    CALL Info('ScaleLinearSystem','Scaling diagonal entries to unity',Level=10)

    IF( PRESENT( DiagScaling ) ) THEN
      CALL Info('ScaleLinearSystem','Reusing existing > DiagScaling < vector',Level=12)
      Diag => DiagScaling 
    ELSE
      CALL Info('ScaleLinearSystem','Computing > DiagScaling < vector',Level=12)
      IF(.NOT. ASSOCIATED(A % DiagScaling)) THEN
        ALLOCATE( A % DiagScaling(n) ) 
      END IF
      Diag => A % DiagScaling
      Diag = 0._dp
    
      ComplexMatrix = Solver % Matrix % COMPLEX

      IF( ListGetLogical( Solver % Values,'Linear System Pseudo Complex',Found ) ) THEN
        ComplexMatrix = .TRUE.
      END IF
      
      IF ( ComplexMatrix ) THEN
        CALL Info('ScaleLinearSystem','Assuming complex matrix while scaling',Level=20)

        !$OMP PARALLEL DO &
        !$OMP SHARED(Diag, A, N) &
        !$OMP PRIVATE(i, j) &
        !$OMP DEFAULT(NONE)
        DO i=1,n,2
          j = A % Diag(i)
          IF(j>0) THEN
            Diag(i)   = A % Values(j)
            Diag(i+1) = A % Values(j+1)
          ELSE
            Diag(i) = 0._dp
            Diag(i+1) = 0._dp
          END IF
        END DO
        !$OMP END PARALLEL DO
      ELSE
        CALL Info('ScaleLinearSystem','Assuming real valued matrix while scaling',Level=25)

        !$OMP PARALLEL DO &
        !$OMP SHARED(Diag, A, N) &
        !$OMP PRIVATE(i, j) &
        !$OMP DEFAULT(NONE)
        DO i=1,n
          j = A % Diag(i)
          IF (j>0) Diag(i) = A % Values(j)
        END DO
        !$OMP END PARALLEL DO
      END IF
      
      IF ( ParEnv % PEs > 1 ) CALL ParallelSumVector(A, Diag)

      IF ( ComplexMatrix ) THEN
        !$OMP PARALLEL DO &
        !$OMP SHARED(Diag, A, N) &
        !$OMP PRIVATE(i, j, DiagC, s) &
        !$OMP DEFAULT(NONE)
        DO i=1,n,2
          DiagC = CMPLX(Diag(i),-Diag(i+1),KIND=dp)

          s = SQRT( ABS( DiagC ) )
          IF( s > TINY(s) ) THEN 
            Diag(i)   = 1.0_dp / s
            Diag(i+1) = 1.0_dp / s
          ELSE
            Diag(i)   = 1.0_dp
            Diag(i+1) = 1.0_dp
          END IF
        END DO
        !$OMP END PARALLEL DO
      ELSE
        s = 0.0_dp
        ! TODO: Add threading
        IF (ANY(ABS(Diag) <= TINY(bnorm))) s=1
        s = ParallelReduction(s,2) 

        IF(s > TINY(s) ) THEN 
          DO i=1,n
            IF ( ABS(Diag(i)) <= TINY(bnorm) ) THEN
              Diag(i) = SUM( ABS(A % Values(A % Rows(i):A % Rows(i+1)-1)) )
            ELSE
              j = A % Diag(i)
              IF (j>0) Diag(i) = A % Values(j)
            END IF
          END DO
          IF ( ParEnv % PEs > 1 ) CALL ParallelSumVector(A, Diag)
        END IF

        !$OMP PARALLEL DO &
        !$OMP SHARED(Diag, N, bnorm) &
        !$OMP PRIVATE(i) &
        !$OMP DEFAULT(NONE)
        DO i=1,n
          IF ( ABS(Diag(i)) > TINY(bnorm) ) THEN
            Diag(i) = 1.0_dp / SQRT(ABS(Diag(i)))
          ELSE
            Diag(i) = 1.0_dp
          END IF
        END DO
        !$OMP END PARALLEL DO
      END IF
    END IF

    
    ! Optionally we may just create the diag and leave the scaling undone
    !--------------------------------------------------------------------
    IF( PRESENT( ApplyScaling ) ) THEN
      IF(.NOT. ApplyScaling ) RETURN
    END IF

    CALL Info('ScaleLinearSystem','Scaling matrix values',Level=20)
    
    !$OMP PARALLEL &
    !$OMP SHARED(Diag, A, N) &
    !$OMP PRIVATE(i,j) &
    !$OMP DEFAULT(NONE)

    !$OMP DO
    DO i=1,n
      DO j = A % Rows(i), A % Rows(i+1)-1
        A % Values(j) = A % Values(j) * &
            ( Diag(i) * Diag(A % Cols(j)) )
      END DO
    END DO
    !$OMP END DO NOWAIT

    ! Dont know why this was temporarily commented off....
#if 1
    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN 
        CALL Info('ScaleLinearSystem','Scaling PrecValues',Level=20)
        !$OMP DO
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) * &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
        !$OMP END DO NOWAIT
      END IF
    END IF
#endif

    IF ( ASSOCIATED( A % MassValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
        CALL Info('ScaleLinearSystem','Scaling MassValues',Level=20)
        !$OMP DO
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % MassValues(j) = A % MassValues(j) * &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
        !$OMP END DO NOWAIT
      END IF
    END IF
    
    IF ( ASSOCIATED( A % DampValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
        CALL Info('ScaleLinearSystem','Scaling DampValues',Level=20)
        !$OMP DO
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % DampValues(j) = A % DampValues(j) * &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
        !$OMP END DO NOWAIT
      END IF
    END IF

    !$OMP END PARALLEL

    DoCM = .FALSE.
    IF(PRESENT(ConstraintScaling)) DoCm=ConstraintScaling

    IF(doCM) THEN
      CM => A % ConstraintMatrix
      IF (ASSOCIATED(CM)) THEN
        CALL Info('ScaleLinearSystem','Scaling Constraints',Level=20)
        !$OMP PARALLEL DO &
        !$OMP SHARED(Diag, CM) &
        !$OMP PRIVATE(i,j) &
        !$OMP DEFAULT(NONE)
        DO i=1,CM % NumberOFRows
          DO j=CM % Rows(i), CM % Rows(i+1)-1
            CM % Values(j) = CM % Values(j) * Diag(CM % Cols(j))
          END DO
        END DO
        !$OMP END PARALLEL DO
      END IF
    END IF

    ! Scale r.h.s. and initial guess
    !--------------------------------
    A % RhsScaling=1._dp
    ! TODO: Add threading
    IF( PRESENT( b ) ) THEN
      CALL Info('ScaleLinearSystem','Scaling Rhs vector',Level=20)
      
      b(1:n) = b(1:n) * Diag(1:n)
      DoRHS = .TRUE.
      IF (PRESENT(RhsScaling)) DoRHS = RhsScaling
      IF (DoRHS) THEN
        bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))

        IF( bnorm < SQRT( TINY( bnorm ) ) ) THEN
          CALL Info('ScaleLinearSystem','Rhs vector is almost zero, skipping rhs scaling!',Level=20)
          DoRhs = .FALSE.
          bnorm = 1.0_dp
        END IF
      ELSE
        bnorm = 1.0_dp
      END IF
      
      A % RhsScaling = bnorm

      IF( DoRhs ) THEN
        Diag(1:n) = Diag(1:n) * bnorm
        b(1:n) = b(1:n) / bnorm
      END IF
      
      IF( PRESENT( x) ) THEN
        x(1:n) = x(1:n) / Diag(1:n)
      END IF
    END IF

    
    !-----------------------------------------------------------------------------
  END SUBROUTINE ScaleLinearSystem
!-----------------------------------------------------------------------------


!-----------------------------------------------------------------------------
!>   Equilibrate the rows of the coefficient matrix A to
!>   minimize the condition number. The associated rhs vector f is also scaled.
!------------------------------------------------------------------------------
  SUBROUTINE RowEquilibration( A, f, Parallel )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: f(:)
    LOGICAL :: Parallel
!-----------------------------------------------------------------------------
    LOGICAL :: ComplexMatrix
    INTEGER :: i, j, n 
    REAL(kind=dp) :: norm, tmp
    INTEGER, POINTER :: Cols(:), Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-------------------------------------------------------------------------

    CALL Info('RowEquilibration','Scaling system such that abs rowsum is unity',Level=15)
        

    n = A % NumberOfRows
    ComplexMatrix = A % COMPLEX

    Rows   => A % Rows
    Cols   => A % Cols
    Values => A % Values

    IF( .NOT. ASSOCIATED(A % DiagScaling) ) THEN
      ALLOCATE( A % DiagScaling(n) ) 
    END IF
    Diag => A % DiagScaling    
    
    Diag = 0.0d0
    norm = 0.0d0

    !---------------------------------------------
    ! Compute 1-norm of each row
    !---------------------------------------------
    IF (ComplexMatrix) THEN
      DO i=1,n,2
        tmp = 0.0d0
        DO j=Rows(i),Rows(i+1)-1,2
          tmp = tmp + ABS( CMPLX( Values(j), -Values(j+1), kind=dp ) )
        END DO
        Diag(i) = tmp
        Diag(i+1) = tmp
      END DO
    ELSE
      DO i=1,n
        tmp = 0.0d0
        DO j=Rows(i),Rows(i+1)-1        
          tmp = tmp + ABS(Values(j))          
        END DO
        Diag(i) = tmp       
      END DO
    END IF

    IF (Parallel) THEN
      CALL ParallelSumVector(A, Diag)
    END IF
    norm = MAXVAL(Diag(1:n))
    IF( Parallel ) THEN
      norm = ParallelReduction(norm,2)
    END IF

    !--------------------------------------------------
    ! Now, define the scaling matrix by inversion and 
    ! perform the actual scaling of the linear system
    !--------------------------------------------------
    IF (ComplexMatrix) THEN    
      DO i=1,n,2
        IF (Diag(i) > TINY(norm) ) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
        Diag(i+1) = Diag(i)
      END DO
    ELSE
      DO i=1,n      
        IF (Diag(i) > TINY(norm)) THEN
          Diag(i) = 1.0_dp / Diag(i)
        ELSE
          Diag(i) = 1.0_dp
        END IF
      END DO
    END IF

    DO i=1,n    
      DO j=Rows(i),Rows(i+1)-1
        Values(j) = Values(j) * Diag(i)
      END DO
      f(i) = Diag(i) * f(i)
    END DO


    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) * Diag(i) 
          END DO
        END DO
      END IF
    END IF

    
    WRITE( Message, * ) 'Unscaled matrix norm: ', norm    
    CALL Info( 'OptimalMatrixScaling', Message, Level=5 )

!------------------------------------------------------------------------------
  END SUBROUTINE RowEquilibration
!------------------------------------------------------------------------------


  
!--------------------------------------------------------------
!>  Scale the system back to original.
!--------------------------------------------------------------
  SUBROUTINE BackScaleLinearSystem( Solver,A,b,x,DiagScaling,&
      ConstraintScaling, EigenScaling ) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t) :: A
    REAL(KIND=dp), OPTIONAL :: b(:),x(:)
    LOGICAL, OPTIONAL :: ConstraintScaling, EigenScaling
    REAL(KIND=dp), OPTIONAL, TARGET :: DiagScaling(:)

    REAL(KIND=dp), POINTER :: Diag(:)
    REAL(KIND=dp) :: bnorm
    INTEGER :: n,i,j
    LOGICAL :: doCM

    TYPE(Matrix_t), POINTER :: CM

    CALL Info('BackScaleLinearSystem','Scaling back to original scale',Level=14)

    
    n = A % NumberOfRows
    
    IF( PRESENT( DiagScaling ) ) THEN
      Diag => DiagScaling
    ELSE  
      Diag => A % DiagScaling
    END IF

    IF(.NOT. ASSOCIATED( Diag ) ) THEN
      CALL Warn('BackScaleLinearSystem','Diag not associated!')
      RETURN
    END IF
    IF( SIZE( Diag ) /= n ) THEN
      CALL Fatal('BackScaleLinearSystem','Diag of wrong size!')
    END IF 

    IF( PRESENT( b ) ) THEN
       ! TODO: Add threading
! 
!      Solve x:  INV(D)x = y, scale b back to orig
!      -------------------------------------------
      IF( PRESENT( x ) ) THEN
        x(1:n) = x(1:n) * Diag(1:n)
      END IF
      bnorm = A % RhsScaling
      Diag(1:n) = Diag(1:n) / bnorm
      b(1:n) = b(1:n) / Diag(1:n) * bnorm
    END IF
    
    IF( PRESENT( EigenScaling ) ) THEN
      IF( EigenScaling ) THEN
        ! TODO: Add threading
        DO i=1,Solver % NOFEigenValues
          !
          !           Solve x:  INV(D)x = y
          !           --------------------------
          IF ( Solver % Matrix % COMPLEX ) THEN
            Solver % Variable % EigenVectors(i,1:n/2) = &
                Solver % Variable % EigenVectors(i,1:n/2) * Diag(1:n:2)
          ELSE
            Solver % Variable % EigenVectors(i,1:n) = &
                Solver % Variable % EigenVectors(i,1:n) * Diag(1:n)
          END IF
        END DO
      END IF
    END IF
    
    !$OMP PARALLEL &
    !$OMP SHARED(Diag, A, N) &
    !$OMP PRIVATE(i, j) &
    !$OMP DEFAULT(NONE)
    
    !$OMP DO
    DO i=1,n
      DO j=A % Rows(i), A % Rows(i+1)-1
        A % Values(j) = A % Values(j) / (Diag(i) * Diag(A % Cols(j)))
      END DO
    END DO
    !$OMP END DO NOWAIT
    
#if 0
    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        !$OMP DO
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) / &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
        !$OMP END DO NOWAIT
      END IF
    END IF
#endif
    IF ( ASSOCIATED( A % MassValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % MassValues)) THEN
        !$OMP DO
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % MassValues(j) = A % MassValues(j) / &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
        !$OMP END DO NOWAIT
      END IF
    END IF
    
    IF ( ASSOCIATED( A % DampValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % DampValues)) THEN
        !$OMP DO
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % DampValues(j) = A % DampValues(j) / &
                ( Diag(i) * Diag(A % Cols(j)) )
          END DO
        END DO
        !$OMP END DO NOWAIT
      END IF
    END IF

    !$OMP END PARALLEL

    ! TODO: Add threading
    doCM=.FALSE.
    IF(PRESENT(ConstraintScaling)) doCM=ConstraintScaling
    IF(doCM) THEN
      CM => A % ConstraintMatrix
      IF (ASSOCIATED(CM)) THEN
        DO i=1,CM % NumberOFRows
          DO j=CM % Rows(i), CM % Rows(i+1)-1
            CM % Values(j) = CM % Values(j) / ( Diag(CM % Cols(j)) )
          END DO
        END DO
      END IF
    END IF

    A % RhsScaling=1._dp
    DEALLOCATE(A % DiagScaling); A % DiagScaling=>NULL()
    
  END SUBROUTINE BackScaleLinearSystem


!------------------------------------------------------------------------------
!> Scale the linear system back to original when the linear
!> system scaling has been done by row equilibration.
!------------------------------------------------------------------------------
  SUBROUTINE ReverseRowEquilibration( A, f )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    REAL(KIND=dp) :: f(:)
!-----------------------------------------------------------------------------
    INTEGER :: i, j, n
    INTEGER, POINTER :: Rows(:)
    REAL(KIND=dp), POINTER :: Values(:), Diag(:)
!-----------------------------------------------------------------------------
    n = A % NumberOfRows
    Diag => A % DiagScaling   
    Values => A % Values
    Rows => A % Rows

    IF(.NOT. ASSOCIATED( Diag ) ) THEN
      CALL Fatal('ReverseRowEquilibration','Diag not associated!')
    END IF
    IF( SIZE( Diag ) /= n ) THEN
      CALL Fatal('ReverseRowEquilibration','Diag of wrong size!')
    END IF 

    f(1:n) = f(1:n) / Diag(1:n)
    DO i=1,n    
      DO j = Rows(i), Rows(i+1)-1
        Values(j) = Values(j) / Diag(i)
      END DO
    END DO

    IF ( ASSOCIATED( A % PrecValues ) ) THEN
      IF (SIZE(A % Values) == SIZE(A % PrecValues)) THEN
        DO i=1,n
          DO j=A % Rows(i), A % Rows(i+1)-1
            A % PrecValues(j) = A % PrecValues(j) / Diag(i) 
          END DO
        END DO
      END IF
    END IF

    
    DEALLOCATE(A % DiagScaling)
    A % DiagScaling => NULL()

!------------------------------------------------------------------------------
  END SUBROUTINE ReverseRowEquilibration
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computing nodal weight may be good when one needs to transform nodal 
!> information back to continuous fields by dividing with the nodal weight. 
!> Active either for the permutation defined by the primary variable of the 
!> solver, or for a permutation vector defined by an optional flag that
!> is used as a mask to define the set of active nodes.
!------------------------------------------------------------------------------
  SUBROUTINE CalculateNodalWeights(Solver,WeightAtBoundary,&
      Perm,VarName,Var)
!------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(Solver_t) :: Solver
    LOGICAL :: WeightAtBoundary
    INTEGER, POINTER, OPTIONAL :: Perm(:)
    CHARACTER(*), OPTIONAL :: VarName
    TYPE(Variable_t), POINTER, OPTIONAL :: Var
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: IntVarName
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: WeightsVar
    TYPE(ValueList_t), POINTER :: ElemParams
    REAL(KIND=dp), POINTER :: Weights(:), Solution(:)    
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER ::k, e, t, n, ElemStart, ElemFin, Coordsys
    INTEGER, POINTER :: IntPerm(:), Indexes(:),LocalIndexes(:)
    REAL(KIND=dp) :: u,v,w,s,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    LOGICAL :: GotIt, stat, VariableOutput, UseMask, RequireLogical, Hit
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3)


    Mesh => Solver % Mesh
    CoordSys = CurrentCoordinateSystem()

    NULLIFY( WeightsVar ) 
    IF( PRESENT( VarName ) ) THEN
      IntVarName = VarName
    ELSE IF ( WeightAtBoundary ) THEN
      IntVarName = GetVarName(Solver % Variable) // ' Boundary Weights'
    ELSE
      IntVarName = GetVarName(Solver % Variable) // ' Weights'
    END IF
    WeightsVar => VariableGet( Mesh % Variables, IntVarName )

    IF( WeightAtBoundary ) THEN
      ElemStart = Mesh % NumberOfBulkElements + 1
      ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      UseMask = ListCheckPresentAnyBC( CurrentModel, IntVarName )
    ELSE
      ElemStart = 1
      ElemFin = Mesh % NumberOfBulkElements 
      UseMask = ListCheckPresentAnyBodyForce( CurrentModel, IntVarName )
    END IF

    RequireLogical = .FALSE.
    NULLIFY( IntPerm ) 
    IF ( .NOT. ASSOCIATED(WeightsVar) ) THEN
      IF( PRESENT( Perm ) ) THEN
        IntPerm => Perm 
      ELSE
        IntPerm => Solver % Variable % Perm
      END IF
      IF( ASSOCIATED( IntPerm ) ) THEN
	NULLIFY( Solution )
	n = MAXVAL( IntPerm ) 
        ALLOCATE( Solution(n))
        Solution = 0.0d0        
        CALL VariableAdd( Mesh % Variables, Mesh, Solver,&
            IntVarName, 1, Solution, IntPerm )
        NULLIFY( Solution )
      ELSE
        CALL Warn('CalculateNodalWeights','Permutation vector not present?')
        RETURN
      END IF
      WeightsVar => VariableGet( Mesh % Variables, IntVarName )
    END IF

    IF( .NOT. ASSOCIATED( WeightsVar ) ) THEN
      CALL Fatal('CalculateNodalWeights','Solution variable not present?')
    END IF
    Weights => WeightsVar % Values
    IntPerm => WeightsVar % Perm
    IF ( .NOT. ASSOCIATED(Weights) ) THEN
      CALL Warn('CalculateNodalWeights','Solution vector not present?')
      RETURN
    ELSE
      IF( PRESENT( Var) ) Var => WeightsVar
    END IF

    CALL Info('ComputeNodalWeights',&
        'Computing weights for solver to variable: '//TRIM(IntVarName))
    n = Mesh % MaxElementNodes

    ALLOCATE(Basis(n), ElementNodes % x(n), ElementNodes % y(n), &
        ElementNodes % z(n), LocalIndexes(n) )
    Weights = 0.0_dp

    DO e=ElemStart,ElemFin

      Element => Mesh % Elements( e )
      Indexes => Element % NodeIndexes

      n = Element % TYPE % NumberOfNodes
      LocalIndexes(1:n) = IntPerm( Indexes ) 
      IF( ANY( LocalIndexes(1:n) == 0 ) ) CYCLE

      IF( UseMask ) THEN
        Hit = .FALSE.
        IF( WeightAtBoundary ) THEN
          DO k=1,CurrentModel % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(k) % Tag ) THEN
              Hit = .TRUE.
              EXIT
            END IF
          END DO
          IF( .NOT. Hit ) CYCLE
          ElemParams => CurrentModel % BCs(k) % Values
        ELSE
          ElemParams => CurrentModel % Bodies(Element % BodyId) % Values
        END IF
        IF( RequireLogical ) THEN
          IF( .NOT. ListGetLogical( ElemParams, IntVarName, Stat ) ) CYCLE
        ELSE
          IF( .NOT. ListCheckPresent( ElemParams, IntVarName ) ) CYCLE
        END IF
      END IF

      n = Element % TYPE % NumberOfNodes
      ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes)

      IntegStuff = GaussPoints( Element )

      DO t=1,IntegStuff % n        
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)
        S = IntegStuff % s(t)

        stat = ElementInfo( Element, ElementNodes, U, V, W, detJ, Basis )

        IF ( CoordSys /= Cartesian ) THEN
          X = SUM( ElementNodes % X(1:n) * Basis(1:n) )
          Y = SUM( ElementNodes % Y(1:n) * Basis(1:n) )
          Z = SUM( ElementNodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * SqrtMetric
        END IF
        
        Weights( LocalIndexes(1:n) ) = &
            Weights( LocalIndexes(1:n) ) + s * detJ * Basis(1:n)
      END DO

    END DO

    DEALLOCATE(Basis, ElementNodes % x, ElementNodes % y, &
        ElementNodes % z, LocalIndexes )

    CALL Info('ComputeNodalWeights','All done')

  END SUBROUTINE CalculateNodalWeights
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Compute weights of entities i.e. their area and volume in the mesh.
!------------------------------------------------------------------------------
  SUBROUTINE CalculateEntityWeights(Model, Mesh)
!------------------------------------------------------------------------------
    IMPLICIT NONE 
    TYPE(Model_t) :: Model 
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element, Left, Right
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER ::i,j,k, e,t, n, Coordsys, TrueOwner, bc_id, bf_id, mat_id, body_id, &
        maxsize, ierr, PotOwner
    INTEGER :: NoBC, NoBodies, NoBF, NoMat
    INTEGER, POINTER :: Indexes(:)
    REAL(KIND=dp) :: u,v,w,s,detJ
    REAL(KIND=dp), ALLOCATABLE :: Basis(:)
    REAL(KIND=dp), POINTER :: bc_weights(:),body_weights(:),&
        mat_weights(:),bf_weights(:),tmp_weights(:)
    REAL(KIND=dp) :: x,y,z,Metric(3,3),SqrtMetric,Symb(3,3,3),dSymb(3,3,3,3), &
        Coeff
    LOGICAL :: Found, Stat, BodyElem


    CoordSys = CurrentCoordinateSystem()

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Warn('CalculateEntityWeights','Mesh not associated!')
      RETURN
    END IF

    CALL Info('ComputeNodalWeights','Computing weights for the mesh entities',Level=6)
    n = Mesh % MaxElementNodes

    NoBC = Model % NumberOfBCs
    NoBodies = Model % NumberOfBodies
    NoMat = Model % NumberOfMaterials
    NoBF = Model % NumberOfBodyForces
    

    ALLOCATE(Basis(n), &
        ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )

    IF( .NOT. Mesh % EntityWeightsComputed ) THEN
      IF( NoBC > 0 ) ALLOCATE( Mesh % BCWeight(NoBC) )
      IF( NoBodies > 0 ) ALLOCATE( Mesh % BodyWeight(NoBodies ) ) 
      IF( NoMat > 0 ) ALLOCATE( Mesh % MaterialWeight(NoMat) ) 
      IF( NoBF > 0 ) ALLOCATE( Mesh % BodyForceWeight(NoBF ) )       
    END IF

    IF( NoBC > 0 ) THEN
      bc_weights => Mesh % BCWeight
      bc_weights(1:NoBC ) = 0.0_dp
    END IF

    IF( NoBodies > 0 ) THEN
      body_weights => Mesh % BodyWeight
      body_weights(1:NoBodies ) = 0.0_dp
    END IF

    IF( NoMat > 0 ) THEN
      mat_weights => Mesh % MaterialWeight
      mat_weights(1:NoMat ) = 0.0_dp
    END IF

    IF( NoBF > 0 ) THEN
      bf_weights => Mesh % BodyForceWeight       
      bf_weights(1:NoBF ) = 0.0_dp
    END IF

    DO e=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      bf_id = 0
      mat_id = 0
      body_id = 0
      bc_id = 0
      Coeff = 1.0_dp

      BodyElem = ( e <= Mesh % NumberOfBulkElements ) 
      Element => Mesh % Elements( e )

      IF( BodyElem ) THEN
        body_id = Element % BodyId
        bf_id = ListGetInteger( Model % Bodies(body_id) % Values,&
            'Body Force',Found)
        mat_id = ListGetInteger( Model % Bodies(body_id) % Values,&
            'Material',Found)
      ELSE
        Found = .FALSE.
        DO bc_id = 1,Model % NumberOfBCs
          Found = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) 
          IF( Found ) EXIT
        END DO
        IF(.NOT. Found) CYCLE
      END IF

      Coeff = 1.0_dp
      
      ! In parallel compute the weight only at their true owners.
      ! Therefore cycle the halo elements. For shared BCs 
      ! take only half of the weight. 
      IF( ParEnv % PEs > 1 ) THEN
        IF( BodyElem ) THEN
          IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
        ELSE
          TrueOwner = 0
          PotOwner = 0
          Left => Element % BoundaryInfo % Left
          IF( ASSOCIATED( Left ) ) THEN
            PotOwner = PotOwner + 1
            IF( Left % PartIndex == ParEnv % MyPe ) TrueOwner = TrueOwner + 1
          END IF
          Right => Element % BoundaryInfo % Right
          IF( ASSOCIATED( Right ) ) THEN
            PotOwner = PotOwner + 1
            IF( Right % PartIndex == ParEnv % MyPe ) TrueOwner = TrueOwner + 1
          END IF
          IF( PotOwner > 0 ) THEN
            IF( TrueOwner == 0 ) CYCLE
            Coeff = 1.0_dp * TrueOwner / PotOwner
          END IF
        END IF
      END IF

      Indexes => Element % NodeIndexes

      n = Element % TYPE % NumberOfNodes
      ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes)

      IntegStuff = GaussPoints( Element )

      DO t=1,IntegStuff % n        
        U = IntegStuff % u(t)
        V = IntegStuff % v(t)
        W = IntegStuff % w(t)

        stat = ElementInfo( Element, ElementNodes, U, V, W, detJ, Basis )
        S = Coeff * DetJ * IntegStuff % s(t)

        IF ( CoordSys /= Cartesian ) THEN
          X = SUM( ElementNodes % X(1:n) * Basis(1:n) )
          Y = SUM( ElementNodes % Y(1:n) * Basis(1:n) )
          Z = SUM( ElementNodes % Z(1:n) * Basis(1:n) )
          CALL CoordinateSystemInfo( Metric,SqrtMetric,Symb,dSymb,X,Y,Z )
          s = s * 2 * PI * SqrtMetric
        END IF
                
        IF( bc_id > 0 .AND. bc_id <= NoBC) &
            bc_weights( bc_id ) = bc_weights( bc_id ) + s
        IF( body_id > 0 .AND. body_id <= NoBodies ) &
            body_weights( body_id ) = body_weights( body_id ) + s
        IF( mat_id > 0 .AND. mat_id <= NoMat ) &
            mat_weights( mat_id ) = mat_weights( mat_id ) + s
        IF( bf_id > 0 .AND. bf_id <= NoBF ) &
            bf_weights( bf_id ) = bf_weights( bf_id ) + s
      END DO

    END DO


    IF( ParEnv % PEs > 1 ) THEN
      maxsize = MAX( Model % NumberOfBCs, Model % NumberOfBodies ) 
      ALLOCATE( tmp_weights( maxsize ) ) 
      tmp_weights = 0.0_dp

      IF( NoBC > 0 ) THEN
        tmp_weights(1:NoBC ) = bc_weights
        CALL MPI_ALLREDUCE( tmp_weights, bc_weights, NoBC, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      IF( NoBF > 0 ) THEN
        tmp_weights(1:NoBF ) = bf_weights
        CALL MPI_ALLREDUCE( tmp_weights, bf_weights, NoBF, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      IF( NoBodies > 0 ) THEN
        tmp_weights(1:NoBodies ) = body_weights
        CALL MPI_ALLREDUCE( tmp_weights, body_weights, NoBodies, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      IF( NoMat > 0 ) THEN
        tmp_weights(1:NoMat ) = mat_weights
        CALL MPI_ALLREDUCE( tmp_weights, mat_weights, NoMat, &
            MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, ierr )
      END IF
      DEALLOCATE( tmp_weights )
    END IF

    IF( ParEnv % MyPe == 0 ) THEN
      DO i = 1, NoBC
        PRINT *,'BC weight:',i,bc_weights(i)
      END DO
      DO i = 1, NoBF
        PRINT *,'BF weight:',i,bf_weights(i)
      END DO
      DO i = 1, NoBodies
        PRINT *,'Body weight:',i,body_weights(i)
      END DO
      DO i = 1, NoMat
        PRINT *,'Mat weight:',i,mat_weights(i)
      END DO
    END IF

    DEALLOCATE(Basis, &
        ElementNodes % x, ElementNodes % y, ElementNodes % z )

    Mesh % EntityWeightsComputed = .TRUE.

    CALL Info('CalculateEntityWeights','All done',Level=10)

  END SUBROUTINE CalculateEntityWeights
!------------------------------------------------------------------------------





  SUBROUTINE CalculateLoads( Solver, Aaid, x, DOFs, UseBulkValues, NodalLoads ) 

    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER  :: Aaid
    REAL(KIND=dp) CONTIG :: x(:)
    INTEGER :: DOFs
    LOGICAL :: UseBulkValues
    TYPE(Variable_t), POINTER :: NodalLoads

    REAL(KIND=dp), POINTER :: LoadValues(:)
    INTEGER :: i,j,k,l,m,ii,This,DOF
    REAL(KIND=dp), POINTER :: TempRHS(:), TempVector(:), Rhs(:), TempX(:)
    REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
    REAL(KIND=dp) :: Energy, Energy_im
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: Found, Rotated


    REAL(KIND=dp), ALLOCATABLE :: BoundarySum(:), BufReal(:)
    INTEGER, ALLOCATABLE :: BoundaryShared(:),BoundaryActive(:),DofSummed(:),BufInteg(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: bc, ind, NoBoundaryActive, NoBCs, ierr
    LOGICAL :: OnlyGivenBCs


    IF( .NOT. ASSOCIATED(NodalLoads) ) RETURN
    ALLOCATE( TempVector(Aaid % NumberOfRows) )

    IF( UseBulkValues ) THEN
      SaveValues => Aaid % Values
      Aaid % Values => Aaid % BulkValues
      Rhs => Aaid % BulkRHS
    ELSE
      Rhs => Aaid % Rhs
    END IF


    IF ( ParEnv % PEs > 1 ) THEN
      ALLOCATE(TempRHS(SIZE(Rhs)))
      TempRHS = Rhs 
      CALL ParallelInitSolve( Aaid, x, TempRHS, Tempvector )
      CALL ParallelMatrixVector( Aaid, x, TempVector, .TRUE. )
    ELSE
      CALL MatrixVectorMultiply( Aaid, x, TempVector )
    END IF

    IF( ListGetLogical(Solver % Values, 'Calculate Energy Norm', Found) ) THEN
      Energy = 0._dp
      IF( ListGetLogical(Solver % Values, 'Linear System Complex', Found) ) THEN
        Energy_im = 0._dp
        DO i = 1, (Aaid % NumberOfRows / 2)
          IF ( ParEnv % Pes>1 ) THEN
            IF ( Aaid% ParMatrix % ParallelInfo % &
              NeighbourList(2*(i-1)+1) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
          END IF
          Energy    = Energy    + x(2*(i-1)+1) * TempVector(2*(i-1)+1) - x(2*(i-1)+2) * TempVector(2*(i-1)+2)
          Energy_im = Energy_im + x(2*(i-1)+1) * TempVector(2*(i-1)+2) + x(2*(i-1)+2) * TempVector(2*(i-1)+1) 
       END DO
       Energy    = ParallelReduction(Energy)
       Energy_im = ParallelReduction(Energy_im)

       CALL ListAddConstReal( Solver % Values, 'Energy norm', Energy)
       CALL ListAddConstReal( Solver % Values, 'Energy norm im', Energy_im)

       WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy Norm'
       CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

       WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy Norm im'
       CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy_im )

       WRITE( Message, * ) 'Energy Norm: ', Energy, Energy_im
       CALL Info( 'SolveLinearSystem', Message )
     ELSE 
       DO i=1,Aaid % NumberOfRows
         IF ( ParEnv % Pes>1 ) THEN
           IF ( Aaid % ParMatrix % ParallelInfo % &
                NeighbourList(i) % Neighbours(1) /= Parenv % MyPE ) CYCLE
         END IF
         Energy = Energy + x(i)*TempVector(i)
      END DO
      Energy = ParallelReduction(Energy)
      CALL ListAddConstReal( Solver % Values, 'Energy norm', Energy )

      WRITE( Message,'(A,A,A)') 'res: ',GetVarname(Solver % Variable),' Energy Norm'
      CALL ListAddConstReal( CurrentModel % Simulation, Message, Energy )

      WRITE( Message, * ) 'Energy Norm: ', Energy
      CALL Info( 'SolveLinearSystem', Message )
    END IF
  END IF

    IF ( ParEnv % PEs>1 ) THEN
      DO i=1,Aaid % NumberOfRows
        IF ( AAid % ParallelInfo % NeighbourList(i) % Neighbours(1) == ParEnv % Mype ) THEN
          TempVector(i) = TempVector(i) - TempRHS(i)
        ELSE
          TempVector(i) = 0
        END IF
      END DO
      CALL ParallelSumVector( AAid, Tempvector )
      DEALLOCATE( TempRhs ) 
    ELSE
      TempVector = TempVector - RHS
    END IF


    NoBCs = CurrentModel % NumberOfBCs
    DO This=1,NoBCs
      Projector => CurrentModel  % BCs(This) % PMatrix
      IF (ASSOCIATED(Projector))THEN
        DO DOF=1,DOFs
          DO i=1,Projector % NumberOfRows
            ii = Projector % InvPerm(i)
            IF( ii == 0 ) CYCLE
            k = Solver % Variable % Perm(ii)
            IF(k<=0) CYCLE
            k = DOFs * (k-1) + DOF
            TempVector(k)=0

            DO l = Projector % Rows(i), Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 ) CYCLE
              m = Solver % Variable % Perm( Projector % Cols(l) )
              IF ( m > 0 ) THEN
                m = DOFs * (m-1) + DOF
                TempVector(k) = TempVector(k) + Projector % Values(l)*TempVector(m)
              END IF
            END DO
          END DO
        END DO
      END IF
    END DO

    DO i=1,SIZE( NodalLoads % Perm )
      IF ( NodalLoads % Perm(i)>0 .AND. Solver % Variable % Perm(i)>0 ) THEN
        DO j=1,DOFs
          NodalLoads % Values(DOFs*(NodalLoads % Perm(i)-1)+j) =  &
              TempVector(DOFs*(Solver % Variable % Perm(i)-1)+j)
        END DO
      END IF
    END DO
    DEALLOCATE( TempVector )


    IF( ListGetLogical( Solver % Values,'Calculate Boundary Fluxes',Found ) ) THEN
      CALL Info('CalculateLoads','Computing boundary fluxes from nodal loads',Level=6)

      IF( Solver % Mesh % MaxEdgeDofs > 1 .OR. Solver % Mesh % MaxFaceDOFs > 1 ) THEN
        CALL Warn('CalculateLoads','Boundary flux computation implemented only for nodes for now!')
      END IF

      ALLOCATE( BoundarySum( NoBCs * DOFs ), &
          BoundaryActive( NoBCs ), &
          BoundaryShared( NoBCs ), &
          DofSummed( MAXVAL( NodalLoads % Perm ) ) )
      BoundarySum = 0.0_dp
      BoundaryActive = 0
      BoundaryShared = 0
      DofSummed = 0

      OnlyGivenBCs = ListCheckPresentAnyBC( CurrentModel,'Calculate Boundary Flux')
      
      k = Solver % Mesh % NumberOfBulkElements
      DO i = k+1,k + Solver % Mesh % NumberOfBoundaryElements
        Element => Solver % Mesh % Elements(i)
        bc = Element % BoundaryInfo % Constraint
           
        IF( bc == 0 ) CYCLE

        IF( OnlyGivenBCs ) THEN
          IF (.NOT. ListGetLogical( CurrentModel % BCs(bc) % Values,&
              'Calculate Boundary Flux',Found) ) CYCLE
        END IF

        DO j=1,Element % TYPE % NumberOfNodes
          ind = NodalLoads % Perm( Element % NodeIndexes(j) )
          IF( ind == 0 ) CYCLE

          ! In this partition sum up only the true owners
          IF ( ParEnv % PEs>1 ) THEN
            IF ( AAid % ParallelInfo % NeighbourList(ind) % Neighbours(1) &
                /= ParEnv % Mype ) CYCLE
          END IF

          ! Only sum each entry once. If there is a conflict we cannot 
          ! really resolve it with the chosen method so just warn. 
          IF( DofSummed(ind) == 0 ) THEN
            BoundarySum( DOFs*(bc-1)+1 :DOFs*bc ) = BoundarySum( DOFs*(bc-1)+ 1:DOFs*bc ) + &
                NodalLoads % Values( DOFs*(ind-1) + 1: DOFs * ind )
            DofSummed( ind ) = bc
            BoundaryActive( bc ) = 1
          ELSE IF( bc /= DofSummed(ind) ) THEN
            BoundaryShared(bc) = 1
            BoundaryShared(DofSummed(ind)) = 1
          END IF
        END DO
      END DO
      

      NoBoundaryActive = 0
      IF( ParEnv % PEs > 1 ) THEN
        ALLOCATE( BufInteg( NoBCs ), BufReal( NoBCs * DOFs ) )

        BufInteg = BoundaryActive
        CALL MPI_ALLREDUCE( BufInteg, BoundaryActive, NoBCs, MPI_INTEGER, &
            MPI_SUM, ParEnv % ActiveComm, ierr )
        
        BufInteg = BoundaryShared
        CALL MPI_ALLREDUCE( BufInteg, BoundaryShared, NoBCs, MPI_INTEGER, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        BufReal = BoundarySum 
        CALL MPI_ALLREDUCE( BufReal, BoundarySum, DOFs * NoBCs, MPI_DOUBLE_PRECISION, &
            MPI_SUM, ParEnv % ActiveComm, ierr )

        DEALLOCATE( BufInteg, BufReal ) 
      END IF


      DO i=1,CurrentModel % NumberOfBCs 
        IF( BoundaryActive(i) == 0 ) CYCLE
        IF( BoundaryShared(i) > 0) THEN
          CALL Warn('CalculateLoads','Boundary '//TRIM(I2S(i))//' includes inseparable dofs!')
        END IF
        NoBoundaryActive = NoBoundaryActive + 1

        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over BC '//TRIM(I2S(i))
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//TRIM(I2S(j))//' Flux over BC '//TRIM(I2S(i))
          END IF
          CALL ListAddConstReal( CurrentModel % Simulation, 'res: '//TRIM(Message), &
              BoundarySum(DOFs*(i-1)+j) )
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',BoundarySum(DOFs*(i-1)+j)
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END DO
      
      IF( NoBoundaryActive > 1 ) THEN
        DO j=1,DOFs
          IF( Dofs == 1 ) THEN
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' Flux over all BCs'
          ELSE
            WRITE( Message,'(A)') GetVarname(Solver % Variable)//&
                ' '//TRIM(I2S(j))//' Flux over all BCs'
          END IF
          WRITE( Message,'(A,ES12.5)') TRIM(Message)//': ',SUM(BoundarySum(j::DOFs))
          CALL Info('CalculateLoads',Message,Level=6)
        END DO
      END IF
      
      DEALLOCATE( DofSummed, BoundaryShared, BoundaryActive, BoundarySum )      
    END IF


    IF( UseBulkValues ) THEN
      Aaid % Values => SaveValues
    END IF

  END SUBROUTINE CalculateLoads

  
END MODULE LinsysUtils
