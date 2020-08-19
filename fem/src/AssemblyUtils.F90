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
! *  Utilities for *Solver - routines
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


MODULE AssemblyUtils

#include "../config.h"

   USE LoadMod
   USE ElementUtils
   USE TimeIntegrate
   USE ListMatrix
   USE CRSMatrix
   
   IMPLICIT NONE

   INTERFACE CondensateP
     MODULE PROCEDURE CondensatePR, CondensatePC
   END INTERFACE CondensateP


 CONTAINS


  !> Sets the matrix element to a desired value. 
!------------------------------------------------------------------------------
   SUBROUTINE SetMatrixElement( A, i, j, VALUE )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix
     INTEGER :: i                            !< Row index
     INTEGER :: j                            !< Column index
     REAL(KIND=dp) :: VALUE                  !< Value to be obtained
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_SetMatrixElement( A, i, j, VALUE )
         IF(A % FORMAT == MATRIX_LIST) THEN
           CALL List_toListMatrix(A)
           CALL List_SetMatrixElement( A % ListMatrix, i, j, VALUE )
         END IF

       CASE( MATRIX_LIST )
         CALL List_SetMatrixElement( A % ListMatrix, i, j, VALUE )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_SetMatrixElement( A, i, j, VALUE )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE SetMatrixElement
!------------------------------------------------------------------------------

!> Gets a matrix element. 
!------------------------------------------------------------------------------
   FUNCTION GetMatrixElement( A, i, j ) RESULT ( VALUE )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix
     INTEGER :: i                            !< Row index
     INTEGER :: j                            !< Column index
     REAL(KIND=dp) :: VALUE                  !< Value to be obtained
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         VALUE = CRS_GetMatrixElement( A, i, j )

      CASE( MATRIX_LIST )
         VALUE = List_GetMatrixElement( A % ListMatrix, i, j )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         VALUE = Band_GetMatrixElement( A, i, j )
     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION GetMatrixElement
!------------------------------------------------------------------------------

!> Changes the value of a given matrix element.
!------------------------------------------------------------------------------
   FUNCTION ChangeMatrixElement( A, i, j, NewValue ) RESULT ( OldValue )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,j
     REAL(KIND=dp) :: NewValue, OldValue
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         OldValue = CRS_ChangeMatrixElement( A, i, j, NewValue )

       CASE DEFAULT
         CALL Warn('ChangeMatrixElement','Not implemented for this type')

     END SELECT
!------------------------------------------------------------------------------
   END FUNCTION ChangeMatrixElement
!------------------------------------------------------------------------------


!> Adds to the value of a given matrix element.
!------------------------------------------------------------------------------
   SUBROUTINE AddToMatrixElement( A, i, j,VALUE )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i,j
     REAL(KIND=dp) :: VALUE
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_AddToMatrixElement( A, i, j, VALUE )
         IF(A % FORMAT == MATRIX_LIST) THEN
           CALL List_toListMatrix(A)
           CALL List_AddToMatrixElement( A % ListMatrix, i, j, VALUE )
         END IF

      CASE( MATRIX_LIST )
         CALL List_AddToMatrixElement( A % ListMatrix, i, j, VALUE )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_AddToMatrixElement( A, i, j, VALUE )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE AddToMatrixElement
!------------------------------------------------------------------------------

!> Adds CMPLX value to the value of a given CMPLX matrix element. -ettaka
!------------------------------------------------------------------------------
  SUBROUTINE AddToCmplxMatrixElement(CM, RowId, ColId, Re, Im)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: CM
    INTEGER :: RowId, ColId
    REAL(KIND=dp) :: Re, Im

    CALL AddToMatrixElement(CM, RowId, ColId, Re)
    CALL AddToMatrixElement(CM, RowId, ColId+1, -Im)
    CALL AddToMatrixElement(CM, RowId+1, ColId, Im)
    CALL AddToMatrixElement(CM, RowId+1, ColId+1, Re)

!------------------------------------------------------------------------------
  END SUBROUTINE AddToCmplxMatrixElement
!------------------------------------------------------------------------------

!> Moves a matrix element from one position adding it to the value of another one.
!------------------------------------------------------------------------------
   SUBROUTINE MoveMatrixElement( A, i1, j1, i2, j2 )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: i1,j1,i2,j2
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: VALUE

     VALUE = ChangeMatrixElement(A, i1, j1, 0.0_dp)
     CALL AddToMatrixElement(A, i2, j2, VALUE )
     
!------------------------------------------------------------------------------
   END SUBROUTINE MoveMatrixElement
!------------------------------------------------------------------------------


!> Zeros a row in matrix.
!------------------------------------------------------------------------------
   SUBROUTINE ZeroRow( A, n )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A  !< Structure holding the matrix 
      INTEGER :: n                           !< Row to be zerored.
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_ZeroRow( A,n )

       CASE( MATRIX_LIST )
         CALL List_ZeroRow( A % ListMatrix, n )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_ZeroRow( A,n )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE ZeroRow
!------------------------------------------------------------------------------

!> Moves a row and and sums it with the values of a second one, optionally 
!> multiplying with a constant.
!------------------------------------------------------------------------------
   SUBROUTINE MoveRow( A, n1, n2, Coeff, StayCoeff, MoveCoeff )
!------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     INTEGER :: n1, n2
     REAL(KIND=dp), OPTIONAL :: Coeff, StayCoeff, MoveCoeff
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_MoveRow( A,n1,n2,Coeff,StayCoeff )

         ! If entries are not found the format is changed on-the-fly
         IF( A % FORMAT == MATRIX_LIST ) THEN
           CALL CRS_MoveRow( A,n1,n2,Coeff,StayCoeff ) ! does this make sense?
         END IF
         
       CASE( MATRIX_LIST )
         CALL List_MoveRow( A % ListMatrix,n1,n2,Coeff,StayCoeff )

       CASE DEFAULT
         CALL Warn('MoveRow','Not implemented for this type')
         
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE MoveRow
!------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!> If we have antiperiodic DOFs in periodic system and want to do elimination
!> for conforming mesh, then we need to flip entries in stiffness/mass matrix.
!---------------------------------------------------------------------------
   SUBROUTINE FlipPeriodicLocalMatrix( Solver, n, Indexes, dofs, A )
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: n, dofs
     INTEGER :: Indexes(:)
     REAL(KIND=dp) :: A(:,:)

     LOGICAL, POINTER :: PerFlip(:)
     INTEGER :: i,j,k,l

     IF( .NOT. Solver % PeriodicFlipActive ) RETURN

     PerFlip => Solver % Mesh % PeriodicFlip           

     IF( .NOT. ANY( PerFlip( Indexes(1:n) ) ) ) RETURN
     
     IF( dofs == 1 ) THEN
       DO i=1,n
         DO j=1,n
           IF( XOR(PerFlip(Indexes(i)),PerFlip(Indexes(j))) ) THEN
             A(i,j) = -A(i,j)
           END IF
         END DO
       END DO
     ELSE
       DO i=1,n
         DO j=1,n
           IF( XOR(PerFlip(Indexes(i)),PerFlip(Indexes(j))) ) THEN
             DO k=1,dofs
               DO l=1,dofs
                 A(dofs*(i-1)+k,dofs*(j-1)+l) = -A(dofs*(i-1)+k,dofs*(j-1)+l)
               END DO
             END DO
           END IF
         END DO
       END DO       
     END IF
              
   END SUBROUTINE FlipPeriodicLocalMatrix


!---------------------------------------------------------------------------
!> If we have antiperiodic DOFs in periodic system and want to do elimination
!> for conforming mesh, then we need to flip entries in local force.
!---------------------------------------------------------------------------
   SUBROUTINE FlipPeriodicLocalForce( Solver, n, Indexes, dofs, F )
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: n, dofs
     INTEGER :: Indexes(:)
     REAL(KIND=dp) :: F(:)
     
     LOGICAL, POINTER :: PerFlip(:)
     INTEGER :: i,j

     IF( .NOT. Solver % PeriodicFlipActive ) RETURN

     PerFlip => Solver % Mesh % PeriodicFlip           
     
     IF( .NOT. ANY( PerFlip( Indexes(1:n) ) ) ) RETURN
     
     IF( dofs == 1 ) THEN
       DO i=1,n
         IF( PerFlip(Indexes(i))) F(i) = -F(i)
       END DO
     ELSE
       DO i=1,n
         IF( PerFlip(Indexes(i))) THEN
           DO j=1,dofs
             F(dofs*(i-1)+j) = -F(dofs*(i-1)+j)
           END DO
         END IF
       END DO
     END IF
          
   END SUBROUTINE FlipPeriodicLocalForce


!---------------------------------------------------------------------------
!> Check if there is something to flip.
!---------------------------------------------------------------------------
   FUNCTION AnyFlipPeriodic( Solver, n, Indexes ) RESULT ( DoFlip ) 
     TYPE(Solver_t), POINTER :: Solver
     INTEGER :: n
     INTEGER :: Indexes(:)
     LOGICAL :: DoFlip 
     
     LOGICAL, POINTER :: PerFlip(:)

     DoFlip = .FALSE.
     IF( .NOT. Solver % PeriodicFlipActive ) RETURN
    
     PerFlip => Solver % Mesh % PeriodicFlip                
     DoFlip = ANY( PerFlip(Indexes(1:n)))

   END FUNCTION AnyFlipPeriodic

   
   
!> Glues a local matrix to the global one.
!------------------------------------------------------------------------------
   SUBROUTINE GlueLocalSubMatrix( A,row0,col0,Nrow,Ncol,RowInds,ColInds,&
       RowDofs,ColDofs,LocalMatrix )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: LocalMatrix(:,:)
     TYPE(Matrix_t) :: A
     INTEGER :: Nrow,Ncol,RowDofs,ColDofs,Col0,Row0,RowInds(:),ColInds(:)
!------------------------------------------------------------------------------

     SELECT CASE( A % FORMAT )

       CASE( MATRIX_CRS )       
         CALL CRS_GlueLocalSubMatrix( A,row0,col0,Nrow,Ncol,RowInds,ColInds,&
             RowDofs,ColDofs,LocalMatrix )
      
       CASE( MATRIX_LIST )
         CALL List_GlueLocalSubMatrix( A % ListMatrix,row0,col0,Nrow,Ncol,RowInds,ColInds,&
             RowDofs,ColDofs,LocalMatrix )
        
       CASE DEFAULT
         CALL Warn('GlueLocalSubMatrix','Not implemented for this type')

     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE GlueLocalSubMatrix
!------------------------------------------------------------------------------


  !------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the local matrix containing other coefficients.
!------------------------------------------------------------------------------
   SUBROUTINE Add1stOrderTime( MassMatrix, StiffMatrix,  &
          Force, dt, n, DOFs, NodeIndexes, Solver, UElement )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: MassMatrix(:,:)   !< Local mass matrix.
     REAL(KIND=dp) :: StiffMatrix(:,:)  !< Local stiffness matrix.
     REAL(KIND=dp) :: Force(:)          !< Local right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     INTEGER :: n                       !< number of element nodes
     INTEGER :: DOFs                    !< variable degrees of freedom
     INTEGER :: NodeIndexes(:)          !< element nodes
     TYPE(Solver_t) :: Solver           !< Solver structure.
     TYPE(Element_t), TARGET, OPTIONAL :: UElement !< Element structure
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     INTEGER :: i,j,k,l,m,Order
     REAL(KIND=dp) :: s, t, zeta
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     REAL(KIND=dp) :: PrevSol(DOFs*n,Solver % Order), CurSol(DOFs*n), LForce(n*DOFs)
     TYPE(Variable_t), POINTER :: DtVar
     REAL(KIND=dp) :: Dts(Solver % Order)
     LOGICAL :: ConstantDt
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
     INTEGER :: PredCorrOrder       !< Order of predictor-corrector scheme

     IF ( PRESENT(UElement) ) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF

     IF ( Solver % Matrix % Lumped ) THEN
#ifndef OLD_LUMPING
       s = 0.d0
       t = 0.d0
       DO i=1,n*DOFs
         DO j=1,n*DOFs
           s = s + MassMatrix(i,j)
           IF (i /= j) THEN
             MassMatrix(i,j) = 0.d0
           END IF
         END DO
         t = t + MassMatrix(i,i)
       END DO
  
       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           L = DOFs * (NodeIndexes(i)-1) + j
           IF ( t /= 0.d0 ) THEN
             MassMatrix(K,K) = MassMatrix(K,K) * s / t
           END IF
         END DO
       END DO
#else
       DO i=1,n*DOFs
         s = 0.0d0
         DO j = 1,n*DOFs
           s = s + MassMatrix(i,j)
           MassMatrix(i,j) = 0.0d0
         END DO
         MassMatrix(i,i) = s
       END DO

       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           L = DOFs * (NodeIndexes(i)-1) + j
         END DO
       END DO
#endif
     END IF
!------------------------------------------------------------------------------
     Order = MIN(Solver % DoneTime, Solver % Order)

     DO i=1,n
       DO j=1,DOFs
         K = DOFs * (i-1) + j
         L = DOFs * (NodeIndexes(i)-1) + j
         DO m=1, Order
           PrevSol(K,m) = Solver % Variable % PrevValues(L,m)
         END DO
         CurSol(K) = Solver % Variable % Values(L)
       END DO
     END DO
     
     LForce(1:n*DOFs) = Force(1:n*DOFs)
     CALL UpdateGlobalForce( Solver % Matrix % Force(:,1), LForce, &
         n, DOFs, NodeIndexes, UElement=Element )
!------------------------------------------------------------------------------
!PrevSol(:,Order) needed for BDF
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )

     SELECT CASE( Method )
     CASE( 'fs' ) 
       CALL FractionalStep( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol(:,1), Solver % Beta, Solver )

     CASE('bdf')
       Dts(1) = Dt
       ConstantDt = .TRUE.
       IF(Order > 1) THEN
         DtVar => VariableGet( Solver % Mesh % Variables, 'Timestep size' )
         DO i=2,Order
           Dts(i) = DtVar % PrevValues(1,i-1)
           IF(ABS(Dts(i)-Dts(1)) > 1.0d-6 * Dts(1)) ConstantDt = .FALSE.
         END DO
       END IF
       
       IF(ConstantDt) THEN
         CALL BDFLocal( n*DOFs, dt, MassMatrix, StiffMatrix, Force, PrevSol, &
             Order )
       ELSE     
         CALL VBDFLocal( n*DOFs, dts, MassMatrix, StiffMatrix, Force, PrevSol, &
             Order )
       END IF
       
     CASE('runge-kutta')
       CALL RungeKutta( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol(:,1), CurSol )
       
     CASE('adams-bashforth')
       zeta = ListGetConstReal( Solver % Values, 'Adams Zeta', GotIt )
       IF ( .NOT. Gotit) zeta = 1.0_dp
       PredCorrOrder = ListGetInteger( Solver % Values, &
           'Predictor-Corrector Scheme Order', GotIt)
       IF (.NOT. GotIt) PredCorrOrder = 2
       PredCorrOrder = MIN(PredCorrOrder, Solver % DoneTime /2)       
       CALL AdamsBashforth( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol(:,1), zeta, PredCorrOrder)
       
     CASE('adams-moulton')
       PredCorrOrder = ListGetInteger( Solver % Values, &
           'Predictor-Corrector Scheme Order', GotIt)
       IF (.NOT. GotIt) PredCorrOrder = 2
       PredCorrOrder = MIN(PredCorrOrder, Solver % DoneTime /2)
       CALL AdamsMoulton( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol, PredCorrOrder )      
       
     CASE DEFAULT
       CALL NewmarkBeta( n*DOFs, dt, MassMatrix, StiffMatrix, Force, &
           PrevSol(:,1), Solver % Beta )
     END SELECT
     
!------------------------------------------------------------------------------
   END SUBROUTINE Add1stOrderTime
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the global matrix containing other coefficients.
!------------------------------------------------------------------------------
   SUBROUTINE Add1stOrderTime_CRS( Matrix, Force, dt, Solver )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: Matrix  !< Global matrix (including stiffness and mass)
     REAL(KIND=dp) :: Force(:)          !< Global right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     TYPE(Solver_t) :: Solver           !< Solver structure.
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     INTEGER :: i,j,k,l,m,n,Order
     REAL(KIND=dp) :: s, t, msum
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     TYPE(Variable_t), POINTER :: DtVar
     REAL(KIND=dp) :: Dts(Solver % Order)
     REAL(KIND=dp), POINTER :: PrevSol(:,:), ML(:), CurrSol(:)
     INTEGER, POINTER :: Rows(:), Cols(:)
     LOGICAL :: ConstantDt, Lumped, Found
!------------------------------------------------------------------------------

     CALL Info('Add1stOrderTime_CRS','Adding time discretization to CRS matrix',Level=20)

!------------------------------------------------------------------------------
     Order = MIN(Solver % DoneTime, Solver % Order)
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     CurrSol => Solver % Variable % Values
     PrevSol => Solver % Variable % PrevValues

     
     SELECT CASE( Method )
       
     CASE( 'fs' ) 
       CALL FractionalStep_CRS( dt, Matrix, Force, PrevSol(:,1), Solver )

     CASE('bdf')
       ConstantDt = .TRUE.
       IF(Order > 1) THEN
         Dts(1) = Dt
         DtVar => VariableGet( Solver % Mesh % Variables, 'Timestep size' )
         DO i=2,Order
           Dts(i) = DtVar % PrevValues(1,i-1)
           IF(ABS(Dts(i)-Dts(1)) > 1.0d-6 * Dts(1)) ConstantDt = .FALSE.
         END DO
       END IF

       IF(ConstantDt) THEN
         CALL BDF_CRS( dt, Matrix, Force, PrevSol, Order )
       ELSE     
         CALL VBDF_CRS( dts, Matrix, Force, PrevSol, Order )
       END IF

     CASE('runge-kutta')
       CALL RungeKutta_CRS( dt, Matrix, Force, PrevSol(:,1), CurrSol )

     CASE DEFAULT
       CALL NewmarkBeta_CRS( dt, Matrix, Force, PrevSol(:,1), &
             Solver % Beta )

     END SELECT

!------------------------------------------------------------------------------
   END SUBROUTINE Add1stOrderTime_CRS
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  For time dependent simulations add the time derivative coefficient terms
!>  to the matrix containing other coefficients.
!------------------------------------------------------------------------------
   SUBROUTINE Add2ndOrderTime( MassMatrix, DampMatrix, StiffMatrix,  &
         Force, dt, n, DOFs, NodeIndexes, Solver )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: MassMatrix(:,:)   !< Local mass matrix.
     REAL(KIND=dp) :: DampMatrix(:,:)   !< Local damping matrix.
     REAL(KIND=dp) :: StiffMatrix(:,:)  !< Local stiffness matrix.
     REAL(KIND=dp) :: Force(:)          !< Local right-hand-side vector.
     REAL(KIND=dp) :: dt                !< Simulation timestep size
     INTEGER :: n                       !< number of element nodes
     INTEGER :: DOFs                    !< variable degrees of freedom
     INTEGER :: NodeIndexes(:)          !< element nodes
     TYPE(Solver_t) :: Solver           !< Solver structure.
!------------------------------------------------------------------------------
     LOGICAL :: GotIt
     INTEGER :: i,j,k,l
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     REAL(KIND=dp) :: s,t
     REAL(KIND=dp) :: X(DOFs*n),V(DOFs*N),A(DOFs*N),LForce(n*DOFs)

!------------------------------------------------------------------------------

     IF ( Solver % Matrix % Lumped ) THEN
!------------------------------------------------------------------------------
#ifndef OLD_LUMPING
       s = 0.d0
       t = 0.d0
       DO i=1,n*DOFs
         DO j=1,n*DOFs
           s = s + MassMatrix(i,j)
           IF (i /= j) THEN
             MassMatrix(i,j) = 0.d0
           END IF
         END DO
         t = t + MassMatrix(i,i)
       END DO

       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           IF ( t /= 0.d0 ) THEN
             MassMatrix(K,K) = MassMatrix(K,K) * s / t
           END IF
         END DO
       END DO

       s = 0.d0
       t = 0.d0
       DO i=1,n*DOFs
         DO j=1,n*DOFs
           s = s + DampMatrix(i,j)
           IF (i /= j) THEN
             DampMatrix(i,j) = 0.d0
           END IF
         END DO
         t = t + DampMatrix(i,i)
       END DO

       DO i=1,n
         DO j=1,DOFs
           K = DOFs * (i-1) + j
           IF ( t /= 0.d0 ) THEN
             DampMatrix(K,K) = DampMatrix(K,K) * s / t
           END IF
         END DO
       END DO
#else
!------------------------------------------------------------------------------
!      Lump the second order time derivative terms ...
!------------------------------------------------------------------------------
       DO i=1,n*DOFs
         s = 0.0D0
         DO j=1,n*DOFs
           s = s + MassMatrix(i,j)
           MassMatrix(i,j) = 0.0d0
         END DO
         MassMatrix(i,i) = s
       END DO

!------------------------------------------------------------------------------
!      ... and the first order terms.
!------------------------------------------------------------------------------
       DO i=1,n*DOFs
         s = 0.0D0
         DO j=1,n*DOFs
           s = s + DampMatrix(i,j)
           DampMatrix(i,j) = 0.0d0
         END DO
         DampMatrix(i,i) = s
       END DO
#endif
!------------------------------------------------------------------------------
     END IF
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!    Get previous solution vectors and update current force
!-----------------------------------------------------------------------------
     DO i=1,n
       DO j=1,DOFs
         K = DOFs * (i-1) + j
         IF ( NodeIndexes(i) > 0 ) THEN
           L = DOFs * (NodeIndexes(i)-1) + j
           SELECT CASE(Method)
           CASE DEFAULT
             X(K) = Solver % Variable % PrevValues(L,3)
             V(K) = Solver % Variable % PrevValues(L,4)
             A(K) = Solver % Variable % PrevValues(L,5)
           END SELECT
         END IF
       END DO
     END DO

     LForce(1:n*DOFs) = Force(1:n*DOFs)
     CALL UpdateGlobalForce( Solver % Matrix % Force(:,1), LForce, &
                  n, DOFs, NodeIndexes )
!------------------------------------------------------------------------------
     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     SELECT CASE(Method)
     CASE DEFAULT
       CALL Bossak2ndOrder( n*DOFs, dt, MassMatrix, DampMatrix, StiffMatrix, &
                    Force, X, V, A, Solver % Alpha )
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE Add2ndOrderTime
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Update the right-hand-side of the global equation by adding the local entry. 
!------------------------------------------------------------------------------
   SUBROUTINE UpdateTimeForce( StiffMatrix, &
           ForceVector, LocalForce, n, NDOFs, NodeIndexes )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< Global stiffness matrix.
     REAL(KIND=dp) :: LocalForce(:)     !< Local right-hand-side vector.
     REAL(KIND=dp) :: ForceVector(:)    !< Global right-hand-side vector.
     INTEGER :: n                       !< number of element nodes
     INTEGER :: nDOFs                   !< variable degrees of freedom
     INTEGER :: NodeIndexes(:)          !< Element node to global node numbering mapping.
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
!------------------------------------------------------------------------------
     CALL UpdateGlobalForce( StiffMatrix % Force(:,1), LocalForce, &
                     n, NDOFs, NodeIndexes )
     LocalForce = 0.0d0
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateTimeForce
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Add element local matrices & vectors to global matrices and vectors.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateGlobalEquations( StiffMatrix, LocalStiffMatrix, &
      ForceVector, LocalForce, n, NDOFs, NodeIndexes, RotateNT, UElement, &
              GlobalValues )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< The global matrix
     REAL(KIND=dp) :: LocalStiffMatrix(:,:)  !< Local matrix to be added to the global matrix.
     REAL(KIND=dp) :: LocalForce(:)          !< Element local force vector.
     REAL(KIND=dp) :: ForceVector(:)         !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of element nodes. 
     INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
     REAL(KIND=dp), OPTIONAL :: GlobalValues(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,dim, Indexes(n)
     LOGICAL :: Rotate
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
!    Update global matrix and rhs vector....
!------------------------------------------------------------------------------
     IF (PRESENT(UElement)) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF
!------------------------------------------------------------------------------
!    Check first if this element has been defined passive
!------------------------------------------------------------------------------
     IF ( CheckPassiveElement(Element) )  RETURN

!------------------------------------------------------------------------------
     Rotate = .TRUE.
     IF ( PRESENT(RotateNT) ) Rotate = RotateNT

     dim = CoordinateSystemDimension()	
     IF ( Rotate .AND. NormalTangentialNOFNodes > 0 .AND. ndofs>=dim) THEN
       Indexes = 0
       Indexes(1:Element % TYPE % NumberOfNodes) = &
             BoundaryReorder(Element % NodeIndexes)
       CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
          Indexes, BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
     END IF
!------------------------------------------------------------------------------
     IF ( ASSOCIATED( StiffMatrix ) ) THEN
       SELECT CASE( StiffMatrix % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_GlueLocalMatrix( StiffMatrix,n,NDOFs, &
                      NodeIndexes, LocalStiffMatrix, GlobalValues )

       CASE( MATRIX_LIST )
         CALL List_GlueLocalMatrix( StiffMatrix % ListMatrix,n,NDOFs,NodeIndexes, &
                          LocalStiffMatrix )

       CASE( MATRIX_BAND,MATRIX_SBAND )
         CALL Band_GlueLocalMatrix( StiffMatrix,n,NDOFs,NodeIndexes, &
                          LocalStiffMatrix )
       END SELECT
     END IF

     DO i=1,n
       IF ( Nodeindexes(i) > 0 ) THEN
         DO j=1,NDOFs
           k = NDOFs * (NodeIndexes(i)-1) + j
!$omp atomic
           ForceVector(k) = ForceVector(k) + LocalForce(NDOFs*(i-1)+j)
         END DO
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalEquations
!------------------------------------------------------------------------------


!> Add element local matrices & vectors to global matrices and vectors.
!> Vectorized version, does not support normal or tangential boundary
!> conditions yet.
   SUBROUTINE UpdateGlobalEquationsVec( Gmtr, Lmtr, Gvec, Lvec, n, &
           NDOFs, NodeIndexes, RotateNT, UElement, MCAssembly )
     TYPE(Matrix_t), POINTER :: Gmtr         !< The global matrix
     REAL(KIND=dp) CONTIG :: Lmtr(:,:)              !< Local matrix to be added to the global matrix.
     REAL(KIND=dp) CONTIG :: Gvec(:)                !< Element local force vector.
     REAL(KIND=dp) CONTIG :: Lvec(:)                !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of degrees of free per node.
     INTEGER CONTIG :: NodeIndexes(:)               !< Element node to global node numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
     LOGICAL, OPTIONAL :: MCAssembly   !< Assembly process is multicoloured and guaranteed race condition free 

     ! Local variables
     INTEGER :: dim, i,j,k
     INTEGER :: Ind(n*NDOFs)
     REAL(KIND=dp) :: Vals(n*NDOFs)
!DIR$ ATTRIBUTES ALIGN:64::Ind, Vals

     TYPE(Element_t), POINTER :: Element
     LOGICAL :: Rotate
     LOGICAL :: ColouredAssembly, NeedMasking

     IF (PRESENT(UElement)) THEN
       Element => UElement
     ELSE
       Element => CurrentModel % CurrentElement
     END IF
     
     IF ( CheckPassiveElement(Element) )  RETURN
     Rotate = .TRUE.
     IF ( PRESENT(RotateNT) ) Rotate = RotateNT
     
     ColouredAssembly = .FALSE.
     IF ( PRESENT(MCAssembly) ) ColouredAssembly = MCAssembly

     dim = CoordinateSystemDimension()
     ! TEMP
     IF ( Rotate .AND. NormalTangentialNOFNodes > 0 .AND. ndofs>=dim) THEN
     
        DO i=1,Element % TYPE % NumberOfNodes
           Ind(i) = BoundaryReorder(Element % NodeIndexes(i))
        END DO

       ! TODO: See that RotateMatrix is vectorized
       CALL RotateMatrix( Lmtr, Lvec, n, dim, NDOFs, Ind, BoundaryNormals, &
                    BoundaryTangent1, BoundaryTangent2 )

       !IF ( Rotate .AND. NormalTangentialNOFNodes > 0 .AND. ndofs>=dim) THEN
       !  CALL Fatal('UpdateGlobalEquationsVec', &
       !          'Normal or tangential boundary conditions not supported yet!')
     END IF

     NeedMasking = .FALSE.
     DO i=1,n
       IF (NodeIndexes(i)<=0) THEN
         NeedMasking = .TRUE.
         EXIT
       END IF
     END DO
     
     IF ( ASSOCIATED( Gmtr ) ) THEN
       SELECT CASE( Gmtr % FORMAT )
       CASE( MATRIX_CRS )
         CALL CRS_GlueLocalMatrixVec(Gmtr, n, NDOFs, NodeIndexes, Lmtr, ColouredAssembly, NeedMasking)
       CASE DEFAULT
         CALL Fatal('UpdateGlobalEquationsVec','Not implemented for given matrix type')
       END SELECT
     END IF
     
     ! Check for multicolored assembly
     IF (ColouredAssembly) THEN
       IF (NeedMasking) THEN
         ! Vector masking needed, no ATOMIC needed
         !_ELMER_OMP_SIMD PRIVATE(j,k)
         DO i=1,n
           IF (NodeIndexes(i)>0) THEN
             DO j=1,NDOFs
               k = NDOFs*(NodeIndexes(i)-1) + j
               Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
             END DO
           END IF
         END DO
       ELSE
         ! No vector masking needed, no ATOMIC needed
         IF (NDOFS>1) THEN
           !_ELMER_OMP_SIMD PRIVATE(j,k)
           DO i=1,n
             DO j=1,NDOFs
               k = NDOFs*(NodeIndexes(i)-1) + j
               Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
             END DO
           END DO
         ELSE
           !_ELMER_OMP_SIMD
           DO i=1,n
             Gvec(NodeIndexes(i)) = Gvec(NodeIndexes(i)) + Lvec(i)
           END DO
         END IF
       END IF ! Vector masking
     ELSE
       IF (NeedMasking) THEN
         ! Vector masking needed, ATOMIC needed
         DO i=1,n
           IF (NodeIndexes(i)>0) THEN
!DIR$ IVDEP
             DO j=1,NDOFs
               k = NDOFs*(NodeIndexes(i)-1) + j
               !$OMP ATOMIC
               Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
             END DO
           END IF
         END DO
       ELSE
         ! No vector masking needed, ATOMIC needed
         DO i=1,n
!DIR$ IVDEP
           DO j=1,NDOFs
             k = NDOFs*(NodeIndexes(i)-1) + j
             !$OMP ATOMIC
             Gvec(k) = Gvec(k) + Lvec(NDOFs*(i-1)+j)
           END DO
         END DO
       END IF ! Vector masking
     END IF ! Coloured assembly
   END SUBROUTINE UpdateGlobalEquationsVec

!------------------------------------------------------------------------------
!> Update the global vector with the local vector entry.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateGlobalForce(ForceVector, LocalForce, n, &
             NDOFs, NodeIndexes, RotateNT, UElement )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: LocalForce(:)          !< Element local force vector.
     REAL(KIND=dp) :: ForceVector(:)         !< The global RHS vector.
     INTEGER :: n                            !< Number of nodes.
     INTEGER :: NDOFs                        !< Number of element nodes. 
     INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping.
     LOGICAL, OPTIONAL :: RotateNT           !< Should the global equation be done in local normal-tangential coordinates.
     TYPE(Element_t), OPTIONAL, TARGET :: UElement !< Element to be updated
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k, dim,indexes(n)
     LOGICAL :: Rotate
     REAL(KIND=dp) :: LocalStiffMatrix(n*NDOFs,n*NDOFs), LForce(n*NDOFs)
!------------------------------------------------------------------------------
!    Update global matrix and rhs vector....
!------------------------------------------------------------------------------
     IF (PRESENT(UElement)) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

     IF ( CheckPassiveElement( Element ) )  RETURN

     Rotate = .TRUE.
     IF ( PRESENT(RotateNT) ) Rotate=RotateNT

     IF ( Rotate .AND. NormalTangentialNOFNodes>0 ) THEN
       dim = CoordinateSystemDimension()
       Indexes = 0
       ! Element => CurrentModel % CurrentElement
       Indexes(1:Element % TYPE % NumberOfNodes) = &
             BoundaryReorder(Element % NodeIndexes)
       CALL RotateMatrix( LocalStiffMatrix, LocalForce, n, dim, NDOFs, &
          Indexes, BoundaryNormals, BoundaryTangent1, BoundaryTangent2 )
     END IF

     DO i=1,n
       IF ( NodeIndexes(i) > 0 ) THEN
         DO j=1,NDOFs
           k = NDOFs * (NodeIndexes(i)-1) + j
!$omp atomic
           ForceVector(k) = ForceVector(k) + LocalForce(NDOFs*(i-1)+j)
         END DO
       END IF
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateGlobalForce
!------------------------------------------------------------------------------


!> Updates the mass matrix only.
!------------------------------------------------------------------------------
   SUBROUTINE UpdateMassMatrix( StiffMatrix, LocalMassMatrix, &
              n, NDOFs, NodeIndexes, GlobalValues )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: StiffMatrix  !< The global matrix structure
     REAL(KIND=dp) :: LocalMassMatrix(:,:)   !< Local matrix to be added to the global matrix
     INTEGER :: n                            !<  number of nodes in element
     INTEGER :: NDOFs                        !< number of DOFs per node
     INTEGER :: NodeIndexes(:)               !< Element node to global node numbering mapping
     REAL(KIND=dp), OPTIONAL, TARGET :: GlobalValues(:)
!------------------------------------------------------------------------------
     INTEGER :: i,j,k
     REAL(KIND=dp) :: s,t
!------------------------------------------------------------------------------
!    Check first if this element has been defined passive
!------------------------------------------------------------------------------
     IF ( CheckPassiveElement() )  RETURN

!------------------------------------------------------------------------------
!    Update global matrix and rhs vector....
!------------------------------------------------------------------------------

     IF ( StiffMatrix % Lumped ) THEN
       s = 0.d0
       t = 0.d0
       DO i=1,n*NDOFs
          DO j=1,n*NDOFs
             s = s + LocalMassMatrix(i,j)
             IF (i /= j) LocalMassMatrix(i,j) = 0.0d0
          END DO
          t = t + LocalMassMatrix(i,i)
       END DO

        DO i=1,n*NDOFs
           LocalMassMatrix(i,i) = LocalMassMatrix(i,i) * s / t
        END DO
     END IF


     SELECT CASE( StiffMatrix % Format )
        CASE( MATRIX_CRS )
           CALL CRS_GlueLocalMatrix( StiffMatrix, &
                n, NDOFs, NodeIndexes, LocalMassMatrix, GlobalValues )

!       CASE( MATRIX_LIST )
!          CALL List_GlueLocalMatrix( StiffMatrix % ListMatrix, &
!               n, NDOFs, NodeIndexes, LocalMassMatrix )

!      CASE( MATRIX_BAND,MATRIX_SBAND )
!          CALL Band_GlueLocalMatrix( StiffMatrix, &
!               n, NDOFs, NodeIndexes, LocalMassMatrix )

        CASE DEFAULT
          CALL FATAL( 'UpdateMassMatrix', 'Unexpected matrix format')
     END SELECT
!------------------------------------------------------------------------------
   END SUBROUTINE UpdateMassMatrix
!------------------------------------------------------------------------------

 
!------------------------------------------------------------------------------
!> Eliminates bubble degrees of freedom from a local linear system.
!> This version is suitable for flow models with velocity and pressure as 
!> unknowns.
!------------------------------------------------------------------------------
SUBROUTINE NSCondensate( N, Nb, dim, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N, Nb, dim
    REAL(KIND=dp) :: K(:,:), F(:)
    REAL(KIND=dp), OPTIONAL :: F1(:)

    REAL(KIND=dp) :: Kbb(nb*dim,nb*dim)
    REAL(KIND=dp) :: Kbl(nb*dim,n*(dim+1)), Klb(n*(dim+1),nb*dim), Fb(nb*dim)

    INTEGER :: m, i, j, l, p, Cdofs((dim+1)*n), Bdofs(dim*nb)

    m = 0
    DO p = 1,n
      DO i = 1,dim+1
        m = m + 1
        Cdofs(m) = (dim+1)*(p-1) + i
      END DO
    END DO

    m = 0
    DO p = 1,nb
      DO i = 1,dim
        m = m + 1
        Bdofs(m) = (dim+1)*(p-1) + i + n*(dim+1)
      END DO
    END DO

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Cdofs)
    Klb = K(Cdofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb*dim )

    F(1:(dim+1)*n) = F(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    K(1:(dim+1)*n,1:(dim+1)*n) = &
    K(1:(dim+1)*n,1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb,Kbl ) )

    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:(dim+1)*n) = F1(1:(dim+1)*n) - MATMUL( Klb, MATMUL( Kbb, Fb ) )
    END IF
!------------------------------------------------------------------------------
END SUBROUTINE NSCondensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for the static condensation of element bubbles when there are
!> as many bubbles as DOFs left in the matrix (historically this convention
!> was used; now the count of elementwise bubble functions can be chosen
!> flexibly and then the subroutine CondensateP should be called instead).
!------------------------------------------------------------------------------
SUBROUTINE Condensate( N, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N
    REAL(KIND=dp) :: K(:,:),F(:)
    REAL(KIND=dp), OPTIONAL :: F1(:)
!------------------------------------------------------------------------------    
    IF ( PRESENT(F1) ) THEN
      CALL CondensateP( N, N, K, F, F1 )
    ELSE
      CALL CondensateP( N, N, K, F )
    END IF
!------------------------------------------------------------------------------
END SUBROUTINE Condensate
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for condensation of p element bubbles from linear problem.
!> Modifies given stiffness matrix and force vector(s) 
!------------------------------------------------------------------------------
SUBROUTINE CondensatePR( N, Nb, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N               !< Sum of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< Sum of internal (bubble) degrees of freedom.
    REAL(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    REAL(KIND=dp) :: F(:)      !< Local force vector.
    REAL(KIND=dp), OPTIONAL :: F1(:)  !< Local second force vector.
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Kbb(Nb,Nb), Kbl(Nb,N), Klb(N,Nb), Fb(Nb)
    INTEGER :: i, Ldofs(N), Bdofs(Nb)

    IF ( nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL InvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF

    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
END SUBROUTINE CondensatePR
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Subroutine for condensation of p element bubbles from complex-valued linear 
!> problem. Modifies given stiffness matrix and force vector(s) 
!------------------------------------------------------------------------------
SUBROUTINE CondensatePC( N, Nb, K, F, F1 )
!------------------------------------------------------------------------------
    USE LinearAlgebra
    INTEGER :: N               !< Sum of nodal, edge and face degrees of freedom.
    INTEGER :: Nb              !< Sum of internal (bubble) degrees of freedom.
    COMPLEX(KIND=dp) :: K(:,:)    !< Local stiffness matrix.
    COMPLEX(KIND=dp) :: F(:)      !< Local force vector.
    COMPLEX(KIND=dp), OPTIONAL :: F1(:)  !< Local second force vector.
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: Kbb(Nb,Nb), Kbl(Nb,N), Klb(N,Nb), Fb(Nb)
    INTEGER :: i, Ldofs(N), Bdofs(Nb)

    IF ( nb <= 0 ) RETURN

    Ldofs = (/ (i, i=1,n) /)
    Bdofs = (/ (i, i=n+1,n+nb) /)

    Kbb = K(Bdofs,Bdofs)
    Kbl = K(Bdofs,Ldofs)
    Klb = K(Ldofs,Bdofs)
    Fb  = F(Bdofs)

    CALL ComplexInvertMatrix( Kbb,nb )

    F(1:n) = F(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    IF (PRESENT(F1)) THEN
      Fb  = F1(Bdofs)
      F1(1:n) = F1(1:n) - MATMUL( Klb, MATMUL( Kbb, Fb  ) )
    END IF

    K(1:n,1:n) = K(1:n,1:n) - MATMUL( Klb, MATMUL( Kbb, Kbl ) )
!------------------------------------------------------------------------------
  END SUBROUTINE CondensatePC
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Rotate a vector to normal-tangential coordinate system.
!------------------------------------------------------------------------------
  SUBROUTINE RotateNTSystem( Vec, NodeNumber )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Vec(:)
     INTEGER :: NodeNumber
!------------------------------------------------------------------------------
     INTEGER :: i,j,k, dim
     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3)
!------------------------------------------------------------------------------

     IF ( NormalTangentialNOFNodes <= 0 ) RETURN

     dim = CoordinateSystemDimension()

     k = BoundaryReorder(NodeNumber)
     IF ( k <= 0 ) RETURN

     IF ( dim < 3 ) THEN
       Bu = Vec(1)
       Bv = Vec(2)
       Vec(1) =  BoundaryNormals(k,1)*Bu + BoundaryNormals(k,2)*Bv
       Vec(2) = -BoundaryNormals(k,2)*Bu + BoundaryNormals(k,1)*Bv
     ELSE
       Bu = Vec(1)
       Bv = Vec(2)
       Bw = Vec(3)

       RM(:,1) = BoundaryNormals(k,:)
       RM(:,2) = BoundaryTangent1(k,:)
       RM(:,3) = BoundaryTangent2(k,:)

       Vec(1) = RM(1,1)*Bu + RM(2,1)*Bv + RM(3,1)*Bw
       Vec(2) = RM(1,2)*Bu + RM(2,2)*Bv + RM(3,2)*Bw
       Vec(3) = RM(1,3)*Bu + RM(2,3)*Bv + RM(3,3)*Bw
     END IF
!------------------------------------------------------------------------------
  END SUBROUTINE RotateNTSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------------
!> Rotate all components of a solution vector to normal-tangential coordinate system
!------------------------------------------------------------------------------------
  SUBROUTINE RotateNTSystemAll( Solution, Perm, NDOFs )
!------------------------------------------------------------------------------
    REAL(KIND=dp) :: Solution(:)
    INTEGER :: Perm(:), NDOFs
!------------------------------------------------------------------------------
    INTEGER :: i,j,k, dim
    REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3)
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    IF ( NormalTangentialNOFNodes<=0.OR.ndofs<dim ) RETURN

    DO i=1,SIZE(BoundaryReorder)
       k = BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF ( dim < 3 ) THEN
          Bu = Solution(NDOFs*(j-1)+1)
          Bv = Solution(NDOFs*(j-1)+2)

          Solution(NDOFs*(j-1)+1) = BoundaryNormals(k,1)*Bu + BoundaryNormals(k,2)*Bv
          Solution(NDOFs*(j-1)+2) = -BoundaryNormals(k,2)*Bu + BoundaryNormals(k,1)*Bv

       ELSE
          Bu = Solution(NDOFs*(j-1)+1)
          Bv = Solution(NDOFs*(j-1)+2)
          Bw = Solution(NDOFs*(j-1)+3)
 
          RM(:,1) = BoundaryNormals(k,:)
          RM(:,2) = BoundaryTangent1(k,:)
          RM(:,3) = BoundaryTangent2(k,:)

          Solution(NDOFs*(j-1)+1) = RM(1,1)*Bu + RM(2,1)*Bv + RM(3,1)*Bw
          Solution(NDOFs*(j-1)+2) = RM(1,2)*Bu + RM(2,2)*Bv + RM(3,2)*Bw
          Solution(NDOFs*(j-1)+3) = RM(1,3)*Bu + RM(2,3)*Bv + RM(3,3)*Bw
       END IF
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE RotateNTSystemAll
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Backrotate a solution from normal-tangential coordinate system to cartesian one.
!------------------------------------------------------------------------------
  SUBROUTINE BackRotateNTSystem( Solution, Perm, NDOFs )
!------------------------------------------------------------------------------
     REAL(KIND=dp) :: Solution(:)
     INTEGER :: Perm(:), NDOFs
!------------------------------------------------------------------------------
     INTEGER :: i,j,k, dim
     REAL(KIND=dp) :: Bu,Bv,Bw,RM(3,3)
!------------------------------------------------------------------------------
     dim = CoordinateSystemDimension()

     IF ( NormalTangentialNOFNodes<=0.OR.ndofs<dim ) RETURN

     DO i=1,SIZE(BoundaryReorder)
       k = BoundaryReorder(i)
       IF ( k <= 0 ) CYCLE
       j = Perm(i)
       IF ( j <= 0 ) CYCLE

       IF ( dim < 3 ) THEN
         Bu = Solution(NDOFs*(j-1)+1)
         Bv = Solution(NDOFs*(j-1)+2)

         Solution(NDOFs*(j-1)+1) = BoundaryNormals(k,1) * Bu - &
                         BoundaryNormals(k,2) * Bv

         Solution(NDOFs*(j-1)+2) = BoundaryNormals(k,2) * Bu + &
                         BoundaryNormals(k,1) * Bv
       ELSE
         Bu = Solution(NDOFs*(j-1)+1)
         Bv = Solution(NDOFs*(j-1)+2)
         Bw = Solution(NDOFs*(j-1)+3)

         RM(1,:) = BoundaryNormals(k,:)
         RM(2,:) = BoundaryTangent1(k,:)
         RM(3,:) = BoundaryTangent2(k,:)

         Solution(NDOFs*(j-1)+1) = RM(1,1)*Bu + RM(2,1)*Bv + RM(3,1)*Bw
         Solution(NDOFs*(j-1)+2) = RM(1,2)*Bu + RM(2,2)*Bv + RM(3,2)*Bw
         Solution(NDOFs*(j-1)+3) = RM(1,3)*Bu + RM(2,3)*Bv + RM(3,3)*Bw
       END IF
     END DO 
!------------------------------------------------------------------------------
  END SUBROUTINE BackRotateNTSystem
!------------------------------------------------------------------------------

  
END MODULE AssemblyUtils
