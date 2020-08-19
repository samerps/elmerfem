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


MODULE SolverUtils

#include "../config.h"

   USE LoadMod
   USE DirectSolve
   USE Multigrid
   USE IterSolve
   USE ElementUtils
   USE ComponentUtils
   USE TimeIntegrate
   USE ModelDescription
   USE MeshUtils
   USE MortarProjector
   USE LinsysUtils
   USE AssemblyUtils
   USE ParallelUtils
   USE ParallelEigenSolve
   USE ListMatrix
   USE CRSMatrix
   
   IMPLICIT NONE

!   INTERFACE CondensateP
!     MODULE PROCEDURE CondensatePR, CondensatePC
!   END INTERFACE CondensateP

!   CHARACTER(LEN=MAX_NAME_LEN), PRIVATE :: NormalTangentialName
!   INTEGER, PRIVATE :: NormalTangentialNOFNodes
!   INTEGER, POINTER, PRIVATE :: NTelement(:,:)
!   LOGICAL, POINTER, PRIVATE :: NTzeroing_done(:,:)
!   INTEGER, POINTER, PRIVATE :: BoundaryReorder(:)
!   REAL(KIND=dp), POINTER, PRIVATE :: BoundaryNormals(:,:),  &
!                                      BoundaryTangent1(:,:), &
!                                      BoundaryTangent2(:,:)

!   SAVE BoundaryReorder, NormalTangentialNOFNodes, BoundaryNormals, &
!              BoundaryTangent1, BoundaryTangent2, NormalTangentialName

CONTAINS

!> Initialize matrix structure and vector to zero initial value.
!------------------------------------------------------------------------------
   SUBROUTINE InitializeToZero( A, ForceVector )
!------------------------------------------------------------------------------
     TYPE(Matrix_t), POINTER :: A  !< Matrix to be initialized
     REAL(KIND=dp) :: ForceVector(:)         !< vector to be initialized
!------------------------------------------------------------------------------
     INTEGER :: i,dim
     LOGICAL :: Found, AnyNT, AnyProj, DoDisplaceMesh
     TYPE(Solver_t), POINTER :: Solver
!------------------------------------------------------------------------------
     
     CALL Info('InitializeToZero','Initializing the linear system to zero',Level=12)
     
     IF ( ASSOCIATED( A ) ) THEN
       SELECT CASE( A % FORMAT )
         CASE( MATRIX_CRS )
           CALL CRS_ZeroMatrix( A )

         CASE( MATRIX_BAND,MATRIX_SBAND )
           CALL Band_ZeroMatrix( A )
       END SELECT

       IF ( ASSOCIATED(A % PrecValues) ) THEN
         A % PrecValues(:) = 0._dp 
       END IF

       IF ( ASSOCIATED( A % MassValues ) ) THEN
         A % MassValues(:) = 0.d0
       END IF

       IF ( ASSOCIATED( A % DampValues ) ) THEN
         A % DampValues(:) = 0.d0
       END IF

       IF ( ASSOCIATED( A % Force ) ) THEN
         A % Force(:,1) = 0.0d0
       END IF

       IF ( ASSOCIATED( A % RHS_im ) )  THEN
         A % RHS_im(:) = 0.0d0
       END IF
     END IF

     ForceVector = 0.0d0
     Solver => CurrentModel % Solver

     NormalTangentialNOFNodes = 0
     IF ( Solver % Variable % DOFs <= 1 ) RETURN

     NormalTangentialName = 'Normal-Tangential'
     IF ( SEQL(Solver % Variable % Name, 'flow solution') ) THEN
       NormalTangentialName = TRIM(NormalTangentialName) // ' Velocity'
     ELSE
       NormalTangentialName = TRIM(NormalTangentialName) // ' ' // &
                   GetVarName(Solver % Variable)
     END IF

     AnyNT = ListGetLogicalAnyBC( CurrentModel, NormalTangentialName ) 
     AnyProj =  ListGetLogicalAnyBC( CurrentModel, 'Mortar BC Nonlinear')
     IF( .NOT. (AnyNT .OR. AnyProj ) ) RETURN

     DoDisplaceMesh = ListGetLogical( Solver % Values,'Displace Mesh At Init',Found )
     IF( DoDisplaceMesh ) THEN
       CALL Info('InitializeToZero','Displacing mesh for nonlinear projectors',Level=8)
       CALL DisplaceMesh( Solver % Mesh, Solver % variable % Values, 1, &
           Solver % Variable % Perm, Solver % variable % Dofs )
     END IF

     IF( AnyNT ) THEN
       dim = CoordinateSystemDimension()
       CALL CheckNormalTangentialBoundary( CurrentModel, NormalTangentialName, &
           NormalTangentialNOFNodes, BoundaryReorder, &
           BoundaryNormals, BoundaryTangent1, BoundaryTangent2, dim )
       
       CALL AverageBoundaryNormals( CurrentModel, NormalTangentialName, &
           NormalTangentialNOFNodes, BoundaryReorder, &
           BoundaryNormals, BoundaryTangent1, BoundaryTangent2, &
           dim )
     END IF

     IF( AnyProj ) THEN
       CALL GenerateProjectors(CurrentModel,Solver,Nonlinear = .TRUE. )
     END IF

     IF( DoDisplaceMesh ) THEN
       CALL DisplaceMesh( Solver % Mesh, Solver % variable % Values, -1, &
           Solver % Variable % Perm, Solver % variable % Dofs )
     END IF
!------------------------------------------------------------------------------
   END SUBROUTINE InitializeToZero
!------------------------------------------------------------------------------

   

!> Create a child matrix of same toopology but optioanally different size than the
!> parent matrix.
!------------------------------------------------------------------------------

   FUNCTION CreateChildMatrix( ParentMat, ParentDofs, Dofs, ColDofs, CreateRhs, &
       NoReuse, Diagonal ) RESULT ( ChildMat )
     TYPE(Matrix_t) :: ParentMat
     INTEGER :: ParentDofs
     INTEGER :: Dofs
     TYPE(Matrix_t), POINTER :: ChildMat
     INTEGER, OPTIONAL :: ColDofs
     LOGICAL, OPTIONAL :: CreateRhs
     LOGICAL, OPTIONAL :: NoReuse
     LOGICAL, OPTIONAL :: Diagonal
     INTEGER :: i,j,ii,jj,k,l,m,n,nn,Cdofs
     LOGICAL :: ReuseMatrix

     IF( ParentMat % FORMAT /= MATRIX_CRS ) THEN
       CALL Fatal('CreateChildMatrix','Only available for CRS matrix format!')
     END IF

     ChildMat => AllocateMatrix()

     CALL CRS_CreateChildMatrix( ParentMat, ParentDofs, ChildMat, Dofs, ColDofs, CreateRhs, &
         NoReuse, Diagonal )

   END FUNCTION CreateChildMatrix

   

!> Search faces between passive / non-passive domains; add to boundary
!> elements with given bc-id.
!------------------------------------------------------------------------------
  SUBROUTINE GetPassiveBoundary(Model,Mesh,BcId)
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    INTEGER :: BcId
    TYPE(Mesh_t) :: Mesh 

    INTEGER, ALLOCATABLE :: arr(:)
    INTEGER :: i,j,n,cnt,ind, sz
    LOGICAL :: L1,L2
    TYPE(Element_t), POINTER :: Faces(:), Telems(:), Face, P1, P2

    CALL FindMeshEdges(Mesh,.FALSE.)
    SELECT CASE(Mesh % MeshDim)
    CASE(2)
      Faces => Mesh % Edges
      n = Mesh % NumberOfEdges
    CASE(3)
      Faces => Mesh % Faces
      n = Mesh % NumberOfFaces
    END SELECT

    ALLOCATE(arr(n)); cnt=0
    DO i=1,n
      P1 => Faces(i) % BoundaryInfo % Right
      P2 => Faces(i) % BoundaryInfo % Left
      IF ( .NOT. ASSOCIATED(P1) .OR. .NOT. ASSOCIATED(P2) ) CYCLE

      L1 = CheckPassiveElement(P1)
      L2 = CheckPassiveElement(P2)

      IF ( L1.NEQV.L2) THEN
        cnt = cnt+1
        arr(cnt) = i
      END IF
    END DO

    sz = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements - &
             Mesh % PassBCcnt
    IF ( sz+cnt>SIZE(Mesh % Elements) ) THEN
      Telems => Mesh % Elements
      ALLOCATE(Mesh % Elements(sz+cnt))
      IF ( ASSOCIATED(Model % Elements,Telems) ) &
        Model % Elements => Mesh % Elements

      Mesh % Elements(1:sz) = Telems

      ! fix boundary element parent pointers to use new array ...
      ! --------------------------------------------------------
      DO i=1,Mesh % NumberOfBoundaryElements-Mesh % PassBCcnt
        ind = i+Mesh % NumberOfBulkElements
        Face => Mesh % Elements(ind)
        IF ( ASSOCIATED(Face % BoundaryInfo % Left) ) &
          Face % BoundaryInfo % Left  => &
             Mesh % Elements(Face % BoundaryInfo % Left % ElementIndex)
        IF ( ASSOCIATED(Face % BoundaryInfo % Right ) ) &
          Face % BoundaryInfo % Right => &
             Mesh % Elements(Face % BoundaryInfo % Right % ElementIndex)
      END DO

      ! ...likewise for  faces (edges).
      ! -------------------------------
      DO i=1,n
        Face => Faces(i)
        IF ( ASSOCIATED(Face % BoundaryInfo % Left) ) &
          Face % BoundaryInfo % Left  => &
             Mesh % Elements(Face % BoundaryInfo % Left % ElementIndex)
        IF ( ASSOCIATED(Face % BoundaryInfo % Right ) ) &
          Face % BoundaryInfo % Right => &
             Mesh % Elements(Face % BoundaryInfo % Right % ElementIndex)
      END DO

      DEALLOCATE(Telems)
    END IF

    DO i=1,cnt
      sz = sz+1
      Mesh % Elements(sz) = Faces(arr(i))
      Mesh % Elements(sz) % Copy = .TRUE.
      Mesh % Elements(sz) % ElementIndex = sz
      Mesh % Elements(sz) % BoundaryInfo % Constraint = BcId
    END DO
    Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements - &
                Mesh % PassBCcnt + cnt
    Mesh % PassBCcnt = cnt
    IF ( ASSOCIATED(Model % Elements,Mesh % Elements) ) &
      Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
  END SUBROUTINE GetPassiveBoundary
!------------------------------------------------------------------------------


  

!------------------------------------------------------------------------------
!> This subroutine seeks for nodes which are adjacent to the given target node
!> and then creates a couple which corresponds to a given torque. If the 
!> optional definition of the director vector d is given, the torque arm should 
!> ideally be parallel to d and the couple created does not have a d-component. 
!> This version may be more convenient when the torque comes from a dimensionally
!> reduced model over a thin body. Without specifying the director, this 
!> subroutine expects a 3-D geometry.
!
! TO DO: - The target nodes can now be defined only by their indices
!        - Add a way to find the director from the specification of a shell model.
!------------------------------------------------------------------------------
   SUBROUTINE SetCoupleLoads(Model, Perm, A, F, Dofs)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     TYPE(Model_t) :: Model                     !< The current model structure
     INTEGER, POINTER, INTENT(IN) :: Perm(:)    !< The permutation of the associated variable
     TYPE(Matrix_t), INTENT(INOUT) :: A         !< The coefficient matrix of the problem
     REAL(KIND=dp), POINTER, INTENT(INOUT) :: F(:) !< The RHS vector of the problem
     INTEGER, INTENT(IN) :: Dofs                !< The DOF count of the associated variable
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(ValueList_t), POINTER :: ValueList

     LOGICAL :: WithDirector
     LOGICAL :: Found, NoUpperNode, NoLowerNode
     
     INTEGER, ALLOCATABLE :: NearNodes(:) 
     INTEGER, POINTER :: NodeIndexes(:)
     INTEGER, POINTER :: Cols(:), Rows(:), Diag(:)
     INTEGER :: Row, TargetNode, TargetInd, BC, TargetCount
     INTEGER :: i, j, k, l, n, p
     INTEGER :: jx, lx, jy, ly, jz, lz 
     INTEGER :: intarray(1)

     REAL(KIND=dp), ALLOCATABLE :: NearCoordinates(:,:), AllDirectors(:,:), Work(:,:)
     REAL(KIND=dp) :: E(3,3)
     REAL(KIND=dp) :: Torque(3)  ! The torque vector with respect to the global frame
     REAL(KIND=dp) :: d(3)       ! Director at a solid-shell/plate interface    
     REAL(KIND=dp) :: ex(3), ey(3), ez(3)
     REAL(KIND=dp) :: e1(3), e2(3), e3(3)
     REAL(KIND=dp) :: T(3), Force(3), v(3)
     REAL(KIND=dp) :: M1, M2, F1, F2, F3
     REAL(KIND=dp) :: res_x, maxres_x, minres_x
     REAL(KIND=dp) :: res_y, maxres_y, minres_y
     REAL(KIND=dp) :: res_z, maxres_z, minres_z
     REAL(KIND=dp) :: rlower, rupper, FVal, MVal
!------------------------------------------------------------------------------
     IF (.NOT. ListCheckPrefixAnyBC(Model, 'Torque')) RETURN

     Mesh => Model % Solver % Mesh

     IF (.NOT. ASSOCIATED(A % InvPerm)) THEN
       ALLOCATE(A % InvPerm(A % NumberOfRows))
       DO i = 1,SIZE(Perm)
         IF (Perm(i) > 0) THEN
           A % InvPerm(Perm(i)) = i
         END IF
       END DO
     END IF

     ex = [1.0d0, 0.0d0, 0.0d0]
     ey = [0.0d0, 1.0d0, 0.0d0]
     ez = [0.0d0, 0.0d0, 1.0d0]
     E(:,1) = ex
     E(:,2) = ey
     E(:,3) = ez

     Diag   => A % Diag
     Rows   => A % Rows
     Cols   => A % Cols

     DO BC=1,Model % NumberOfBCs
       ValueList => Model % BCs(BC) % Values
       IF (.NOT.ListCheckPresent(ValueList, 'Torque 1') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Torque 2') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Torque 3')) CYCLE
       NodeIndexes => ListGetIntegerArray(ValueList, 'Target Nodes', UnfoundFatal=.TRUE.)

       TargetCount = SIZE(NodeIndexes)
       ALLOCATE(Work(3,TargetCount))
       Work(1,1:TargetCount) = ListGetReal(ValueList, 'Torque 1', TargetCount, NodeIndexes, Found)
       Work(2,1:TargetCount) = ListGetReal(ValueList, 'Torque 2', TargetCount, NodeIndexes, Found)
       Work(3,1:TargetCount) = ListGetReal(ValueList, 'Torque 3', TargetCount, NodeIndexes, Found)

       !
       ! Check whether the torque arm is given by the director vector. This option
       ! is not finalized yet. Here the director definition is sought from the BC
       ! definition, while the director might already be available from the specification 
       ! of a shell model.
       !
       IF (.NOT.ListCheckPresent(ValueList, 'Director 1') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Director 2') .AND. &
           .NOT.ListCheckPresent(ValueList, 'Director 3')) THEN
         WithDirector = .FALSE.
       ELSE
         WithDirector = .TRUE.
         ALLOCATE(AllDirectors(3,TargetCount))
         AllDirectors(1,1:TargetCount) = ListGetReal(ValueList, 'Director 1', TargetCount, NodeIndexes, Found)
         AllDirectors(2,1:TargetCount) = ListGetReal(ValueList, 'Director 2', TargetCount, NodeIndexes, Found)
         AllDirectors(3,1:TargetCount) = ListGetReal(ValueList, 'Director 3', TargetCount, NodeIndexes, Found)
       END IF

       DO p=1,TargetCount
         TargetNode = NodeIndexes(p)
         TargetInd = Perm(NodeIndexes(p))
         IF (TargetInd == 0) CYCLE

         !------------------------------------------------------------------------------
         ! Find nodes which can potentially be used to make a representation of couple:
         !------------------------------------------------------------------------------
         Row = TargetInd * Dofs
         n = (Rows(Row+1)-1 - Rows(Row)-Dofs+1)/DOFs + 1
         ALLOCATE(NearNodes(n), NearCoordinates(3,n))

         k = 0
         DO i = Rows(Row)+Dofs-1, Rows(Row+1)-1, Dofs
           j = Cols(i)/Dofs
           k = k + 1
           NearNodes(k) = A % InvPerm(j)
         END DO
         ! PRINT *, 'POTENTIAL NODE CONNECTIONS:'
         ! print *, 'Nodes near target=', NearNodes(1:k)

         !
         ! The position vectors for the potential nodes where forces may be applied:
         !
         NearCoordinates(1,1:n) = Mesh % Nodes % x(NearNodes(1:n)) - Mesh % Nodes % x(TargetNode)
         NearCoordinates(2,1:n) = Mesh % Nodes % y(NearNodes(1:n)) - Mesh % Nodes % y(TargetNode)
         NearCoordinates(3,1:n) = Mesh % Nodes % z(NearNodes(1:n)) - Mesh % Nodes % z(TargetNode)


         IF (WithDirector) THEN
           !
           ! In this case the torque arm should ideally be parallel to the director vector d.
           ! Construct an orthonormal basis, with d giving the third basis vector.
           !
           d = AllDirectors(:,p)
           e3 = d/SQRT(DOT_PRODUCT(d,d))
           v(1:3) = ABS([DOT_PRODUCT(ex,e3), DOT_PRODUCT(ey,e3), DOT_PRODUCT(ez,e3)]) 
           intarray = MINLOC(v)
           k = intarray(1)
           v(1:3) = E(1:3,k)
           e1 = v - DOT_PRODUCT(v,e3)*e3
           e1 = e1/SQRT(DOT_PRODUCT(e1,e1))
           e2 = CrossProduct(e3,e1)
           !
           ! The torque is supposed to have no component in the direction of d, so remove it
           ! and also find the representation of the altered torque with respect to the local basis:
           !
           Torque = Work(:,p)
           v = DOT_PRODUCT(Torque,e3)*e3
           T = Torque - v
           M1 = DOT_PRODUCT(T,e1)
           M2 = DOT_PRODUCT(T,e2)

           !------------------------------------------------------------------------------
           ! Seek torque arms which are closest to be parallel to d:
           !------------------------------------------------------------------------------
           maxres_z = 0.0d0
           minres_z = 0.0d0
           jz = 0
           lz = 0
           DO i=1,n
             IF (NearNodes(i) == TargetNode) CYCLE
             res_z = DOT_PRODUCT(e3(:), NearCoordinates(:,i)) / &
                 SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
             IF (res_z > 0.0d0) THEN
               !
               ! A near node is on +d side
               !
               IF (res_z > maxres_z) THEN
                 jz = NearNodes(i)
                 maxres_z = res_z
               END IF
             ELSE
               !
               ! A near node is on -d side
               !
               IF (res_z < minres_z) THEN
                 lz = NearNodes(i)
                 minres_z = res_z
               END IF
             END IF
           END DO

           !
           ! Calculate arm lengths with respect to the coordinate axis parallel to d:
           !
           NoUpperNode = .FALSE.
           NoLowerNode = .FALSE.
           IF (jz == 0 .OR. ABS(maxres_z) < AEPS) THEN
             NoUpperNode = .TRUE.
           ELSE
             rupper = DOT_PRODUCT(e3(:), [ Mesh % Nodes % x(jz) - Mesh % Nodes % x(TargetNode), &
                 Mesh % Nodes % y(jz) - Mesh % Nodes % y(TargetNode), &
                 Mesh % Nodes % z(jz) - Mesh % Nodes % z(TargetNode) ])
             ! print *, 'THE NODE ON +d SIDE = ', JZ
             ! print *, 'TORQUE ARM = ', rupper
           END IF

           IF (lz == 0 .OR. ABS(minres_z) < AEPS) THEN
             NoLowerNode = .TRUE.
           ELSE
             rlower = DOT_PRODUCT(-e3(:), [ Mesh % Nodes % x(lz) - Mesh % Nodes % x(TargetNode), &
                 Mesh % Nodes % y(lz) - Mesh % Nodes % y(TargetNode), &
                 Mesh % Nodes % z(lz) - Mesh % Nodes % z(TargetNode) ])
             ! print *, 'THE NODE ON -d SIDE = ', LZ
             ! print *, 'TORQUE ARM = ', rlower
           END IF

           IF (NoUpperNode .OR. NoLowerNode) THEN
             CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite sides')
           ELSE
             !
             ! The torque generated from point loads as M1 * e1 + M2 * e2 = (r e3) x (f1 * e1 - f2 * e2) = 
             ! (r*f2)* e1 + (r*f1)* e2
             !
             F2 = M1/(rupper + rlower)
             F1 = M2/(rupper + rlower)
             Force = F1 * e1 - F2 * e2
             !
             ! Finally compute the components of force with respect to the global frame and
             ! add to the RHS: 
             !
             F1 = DOT_PRODUCT(Force,ex)
             F2 = DOT_PRODUCT(Force,ey)
             F3 = DOT_PRODUCT(Force,ez)

             k = Perm(jz)
             F((k-1)*Dofs+1) = F((k-1)*Dofs+1) + F1
             F((k-1)*Dofs+2) = F((k-1)*Dofs+2) + F2
             IF (Dofs > 2) F((k-1)*Dofs+3) = F((k-1)*Dofs+3) + F3
             k = Perm(lz)
             F((k-1)*Dofs+1) = F((k-1)*Dofs+1) - F1
             F((k-1)*Dofs+2) = F((k-1)*Dofs+2) - F2
             IF (Dofs > 2) F((k-1)*Dofs+3) = F((k-1)*Dofs+3) - F3
           END IF

         ELSE
           !------------------------------------------------------------------------------
           ! Seek torque arms which are closest to be parallel to the global coordinate
           ! axes: 
           !------------------------------------------------------------------------------
           maxres_x = 0.0d0
           minres_x = 0.0d0
           maxres_y = 0.0d0
           minres_y = 0.0d0
           maxres_z = 0.0d0
           minres_z = 0.0d0
           jx = 0
           lx = 0
           jy = 0
           ly = 0
           jz = 0
           lz = 0
           DO i=1,n
             IF (NearNodes(i) == TargetNode) CYCLE

             IF (ABS(Torque(3)) > AEPS) THEN
               res_x = DOT_PRODUCT(ex(:), NearCoordinates(:,i)) / &
                   SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
               IF (res_x > 0.0d0) THEN
                 !
                 ! A near node is on +E_X side
                 !
                 IF (res_x > maxres_x) THEN
                   jx = NearNodes(i)
                   maxres_x = res_x
                 END IF
               ELSE
                 !
                 ! A near node is on -E_X side
                 !
                 IF (res_x < minres_x) THEN
                   lx = NearNodes(i)
                   minres_x = res_x
                 END IF
               END IF
             END IF

             IF (ABS(Torque(1)) > AEPS) THEN
               res_y = DOT_PRODUCT(ey(:), NearCoordinates(:,i)) / &
                   SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
               IF (res_y > 0.0d0) THEN
                 !
                 ! A near node is on +E_Y side
                 !
                 IF (res_y > maxres_y) THEN
                   jy = NearNodes(i)
                   maxres_y = res_y
                 END IF
               ELSE
                 !
                 ! A near node is on -E_Y side
                 !
                 IF (res_y < minres_y) THEN
                   ly = NearNodes(i)
                   minres_y = res_y
                 END IF
               END IF
             END IF

             IF (ABS(Torque(2)) > AEPS) THEN
               res_z = DOT_PRODUCT(ez(:), NearCoordinates(:,i)) / &
                   SQRT(DOT_PRODUCT(NearCoordinates(:,i), NearCoordinates(:,i)))
               IF (res_z > 0.0d0) THEN
                 !
                 ! A near node is on +E_Z side
                 !
                 IF (res_z > maxres_z) THEN
                   jz = NearNodes(i)
                   maxres_z = res_z
                 END IF
               ELSE
                 !
                 ! A near node is on -E_Z side
                 !
                 IF (res_z < minres_z) THEN
                   lz = NearNodes(i)
                   minres_z = res_z
                 END IF
               END IF
             END IF
           END DO

           IF (ABS(Torque(1)) > AEPS) THEN
             !------------------------------------------------------------------------------
             ! Calculate arm lengths with respect to the Y-axis:
             !------------------------------------------------------------------------------
             NoUpperNode = .FALSE.
             NoLowerNode = .FALSE.
             IF (jy == 0) THEN
               NoUpperNode = .TRUE.
             ELSE
               rupper = DOT_PRODUCT(ey(:), [ Mesh % Nodes % x(jy) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(jy) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(jy) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (ly == 0) THEN
               NoLowerNode = .TRUE.
             ELSE
               rlower = DOT_PRODUCT(-ey(:), [ Mesh % Nodes % x(ly) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(ly) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(ly) - Mesh % Nodes % z(TargetNode) ])
             END IF

             !------------------------------------------------------------------------------
             ! Finally, create a couple which tends to cause rotation about the X-axis 
             ! provided nodes on both sides have been identified
             !------------------------------------------------------------------------------
             IF (NoUpperNode .OR. NoLowerNode) THEN
               CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite Y-sides')
             ELSE
               !
               ! The torque M_X E_X = (r E_Y) x (f E_Z), with the force f>0 applied on +E_Y side:
               !
               MVal = Torque(1)
               FVal = Mval/(rupper + rlower)
               k = Perm(jy)
               F((k-1)*Dofs+3) = F((k-1)*Dofs+3) + Fval
               k = Perm(ly)
               F((k-1)*Dofs+3) = F((k-1)*Dofs+3) - Fval
             END IF
           END IF

           IF (ABS(Torque(2)) > AEPS) THEN
             !
             ! Calculate arm lengths with respect to the Z-axis:
             !
             NoUpperNode = .FALSE.
             NoLowerNode = .FALSE.
             IF (jz == 0) THEN
               NoUpperNode = .TRUE.
             ELSE
               rupper = DOT_PRODUCT(ez(:), [ Mesh % Nodes % x(jz) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(jz) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(jz) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (lz == 0) THEN
               NoLowerNode = .TRUE.
             ELSE
               rlower = DOT_PRODUCT(-ez(:), [ Mesh % Nodes % x(lz) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(lz) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(lz) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (NoUpperNode .OR. NoLowerNode) THEN
               CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite Z-sides')
             ELSE
               !
               ! The torque M_Y E_Y = (r E_Z) x (f E_X), with the force f>0 applied on +E_Z side:
               !
               MVal = Torque(2)
               FVal = Mval/(rupper + rlower)
               k = Perm(jz)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) + Fval
               k = Perm(lz)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) - Fval
             END IF
           END IF

           IF (ABS(Torque(3)) > AEPS) THEN
             !
             ! Calculate arm lengths with respect to the X-axis:
             !
             NoUpperNode = .FALSE.
             NoLowerNode = .FALSE.
             IF (jx == 0) THEN
               NoUpperNode = .TRUE.
             ELSE
               rupper = DOT_PRODUCT(ex(:), [ Mesh % Nodes % x(jx) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(jx) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(jx) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (lx == 0) THEN
               NoLowerNode = .TRUE.
             ELSE
               rlower = DOT_PRODUCT(-ex(:), [ Mesh % Nodes % x(lx) - Mesh % Nodes % x(TargetNode), &
                   Mesh % Nodes % y(lx) - Mesh % Nodes % y(TargetNode), &
                   Mesh % Nodes % z(lx) - Mesh % Nodes % z(TargetNode) ])
             END IF

             IF (NoUpperNode .OR. NoLowerNode) THEN
               CALL Warn('SetCoupleLoads', 'A couple BC would need two nodes on opposite Y-sides')
             ELSE
               !
               ! The torque M_Z E_Z = (r E_X) x (f E_Y), with the force f>0 applied on +E_X side:
               !
               MVal = Torque(3)
               FVal = Mval/(rupper + rlower)
               k = Perm(jx)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) - Fval
               k = Perm(lx)
               F((k-1)*Dofs+1) = F((k-1)*Dofs+1) + Fval
             END IF
           END IF
         END IF

         DEALLOCATE(NearNodes, NearCoordinates)
       END DO
       DEALLOCATE(Work)
       IF (WithDirector) DEALLOCATE(AllDirectors)
     END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetCoupleLoads
!------------------------------------------------------------------------------

 

  !-------------------------------------------------------------------------------
  !> Communicate logical tag related to linear system.
  !> This could related to setting Neumann BCs to zero, for example.
  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateLinearSystemTag(A,ZeroDof)
  !-------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A
     LOGICAL, POINTER :: ZeroDof(:)
         
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
     INTEGER :: i,j,k,l,n,nn,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)
     INTEGER :: NewZeros
     
     IF( ParEnv % PEs<=1 ) RETURN

     ALLOCATE( fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )

     nn = 0
     ineigh = 0
     DO i=0, ParEnv % PEs-1
       k = i+1
       IF(.NOT.ParEnv % Active(k) ) CYCLE
       IF(i==ParEnv % myPE) CYCLE
       IF(.NOT.ParEnv % IsNeighbour(k) ) CYCLE
       nn = nn + 1
       fneigh(nn) = k
       ineigh(k) = nn
     END DO

     n = COUNT(ZeroDof .AND. A % ParallelInfo % Interface)
     ALLOCATE( s_e(n, nn ), r_e(n) )

     CALL CheckBuffer( nn*3*n )

     ii = 0
     DO i=1, A % NumberOfRows
       IF(ZeroDof(i) .AND. A % ParallelInfo % Interface(i) ) THEN
          DO j=1,SIZE(A % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = A % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              s_e(ii(k),k) = A % ParallelInfo % GlobalDOFs(i)
            END IF
          END DO
       END IF
     END DO

     DO i=1, nn
       j = fneigh(i) 
       CALL MPI_BSEND( ii(i),1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD,ierr )
       IF( ii(i) > 0 ) THEN
         CALL MPI_BSEND( s_e(1:ii(i),i),ii(i),MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,ierr )
       END IF
     END DO

     NewZeros = 0
     
     DO i=1, nn
       j = fneigh(i)
       CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
       IF ( n>0 ) THEN
         IF( n>SIZE(r_e)) THEN
           DEALLOCATE(r_e)
           ALLOCATE(r_e(n))
         END IF

         CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
         DO j=1,n
           k = SearchNode( A % ParallelInfo, r_e(j), Order=A % ParallelInfo % Gorder )
           IF ( k>0 ) THEN
             IF(.NOT. ZeroDof(k)) THEN
               ZeroDof(k) = .TRUE.
               NewZeros = NewZeros + 1
             END IF
           END IF
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e )
     
     !PRINT *,'New Zeros:',ParEnv % MyPe, NewZeros
     
  !-------------------------------------------------------------------------------
   END SUBROUTINE CommunicateLinearSystemTag
  !-------------------------------------------------------------------------------

#if 0    
!------------------------------------------------------------------------------
  FUNCTION sGetElementDOFs( Indexes, UElement, USolver )  RESULT(NB)
!------------------------------------------------------------------------------
     TYPE(Element_t), OPTIONAL, TARGET :: UElement
     TYPE(Solver_t),  OPTIONAL, TARGET :: USolver
     INTEGER :: Indexes(:)

     TYPE(Solver_t),  POINTER :: Solver
     TYPE(Element_t), POINTER :: Element, Parent

     LOGICAL :: Found, GB
     INTEGER :: nb,i,j,EDOFs, FDOFs, BDOFs,FaceDOFs, EdgeDOFs, BubbleDOFs

     IF ( PRESENT( UElement ) ) THEN
        Element => UElement
     ELSE
        Element => CurrentModel % CurrentElement
     END IF

     IF ( PRESENT( USolver ) ) THEN
        Solver => USolver
     ELSE
        Solver => CurrentModel % Solver
     END IF

     NB = 0

     IF ( Solver % DG ) THEN
        DO i=1,Element % DGDOFs
           NB = NB + 1
           Indexes(NB) = Element % DGIndexes(i)
        END DO

        IF ( ASSOCIATED( Element % BoundaryInfo ) ) THEN
           IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
              DO i=1,Element % BoundaryInfo % Left % DGDOFs
                 NB = NB + 1
                 Indexes(NB) = Element % BoundaryInfo % Left % DGIndexes(i)
              END DO
           END IF
           IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              DO i=1,Element % BoundaryInfo % Right % DGDOFs
                 NB = NB + 1
                 Indexes(NB) = Element % BoundaryInfo % Right % DGIndexes(i)
              END DO
           END IF
        END IF

        IF ( NB > 0 ) RETURN
     END IF

     DO i=1,Element % NDOFs
        NB = NB + 1
        Indexes(NB) = Element % NodeIndexes(i)
     END DO

     FaceDOFs   = Solver % Mesh % MaxFaceDOFs
     EdgeDOFs   = Solver % Mesh % MaxEdgeDOFs
     BubbleDOFs = Solver % Mesh % MaxBDOFs

     IF ( ASSOCIATED( Element % EdgeIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFEdges
          EDOFs = Solver % Mesh % Edges( Element % EdgeIndexes(j) ) % BDOFs
          DO i=1,EDOFs
             NB = NB + 1
             Indexes(NB) = EdgeDOFs*(Element % EdgeIndexes(j)-1) + &
                      i + Solver % Mesh % NumberOfNodes
          END DO
        END DO
     END IF

     IF ( ASSOCIATED( Element % FaceIndexes ) ) THEN
        DO j=1,Element % TYPE % NumberOFFaces
           FDOFs = Solver % Mesh % Faces( Element % FaceIndexes(j) ) % BDOFs
           DO i=1,FDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*(Element % FaceIndexes(j)-1) + i + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
           END DO
        END DO
     END IF

     GB = ListGetLogical( Solver % Values, 'Bubbles in Global System', Found )
     IF (.NOT. Found) GB = .TRUE.

     IF ( ASSOCIATED(Element % BoundaryInfo) ) THEN
       IF (.NOT. isActivePElement(Element) ) RETURN

       Parent => Element % BoundaryInfo % Left
       IF (.NOT.ASSOCIATED(Parent) ) &
         Parent => Element % BoundaryInfo % Right
       IF (.NOT.ASSOCIATED(Parent) ) RETURN

       IF ( ASSOCIATED( Parent % EdgeIndexes ) ) THEN
         EDOFs = Element % BDOFs
         DO i=1,EDOFs
           NB = NB + 1
           Indexes(NB) = EdgeDOFs*(Parent % EdgeIndexes(Element % PDefs % LocalNumber)-1) + &
                    i + Solver % Mesh % NumberOfNodes
         END DO
       END IF

       IF ( ASSOCIATED( Parent % FaceIndexes ) ) THEN
         FDOFs = Element % BDOFs
         DO i=1,FDOFs
           NB = NB + 1
           Indexes(NB) = FaceDOFs*(Parent % FaceIndexes(Element % PDefs % LocalNumber)-1) + i + &
              Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges
         END DO
       END IF
     ELSE IF ( GB ) THEN
        IF ( ASSOCIATED( Element % BubbleIndexes ) ) THEN
           DO i=1,Element % BDOFs
              NB = NB + 1
              Indexes(NB) = FaceDOFs*Solver % Mesh % NumberOfFaces + &
                 Solver % Mesh % NumberOfNodes + EdgeDOFs*Solver % Mesh % NumberOfEdges + &
                   Element % BubbleIndexes(i)
           END DO
        END IF
     END IF
!------------------------------------------------------------------------------
  END FUNCTION SgetElementDOFs
!------------------------------------------------------------------------------
#endif
  
!------------------------------------------------------------------------------
!> Check if Normal / Tangential vector boundary conditions present and
!> allocate space for normals, and if in 3D for two tangent direction
!> vectors.
!------------------------------------------------------------------------------
   SUBROUTINE CheckNormalTangentialBoundary( Model, VariableName, &
     NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals,     &
        BoundaryTangent1, BoundaryTangent2, dim )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model

    CHARACTER(LEN=*) :: VariableName

    INTEGER, POINTER :: BoundaryReorder(:)
    INTEGER :: NumberOfBoundaryNodes,dim

    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
                       BoundaryTangent2(:,:)
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: CurrentElement
    INTEGER :: i,j,k,n,t,ierr,iter, proc
    LOGICAL :: GotIt, Found, Conditional
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)

    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE, TARGET :: n_index(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:)
!------------------------------------------------------------------------------

    ! need an early initialization to average normals across partitions:
    !-------------------------------------------------------------------
    IF ( Parenv  % PEs >1 ) THEN
      IF (.NOT. ASSOCIATED(Model % Solver % Matrix % ParMatrix) ) &
         CALL ParallelInitMatrix( Model % Solver, Model % Solver % Matrix )
    END IF

    NumberOfBoundaryNodes = 0

    Found = .FALSE.
    DO i=1,Model % NumberOfBCs
      IF ( ListGetLogical(Model % BCs(i) % Values, VariableName, Gotit) ) THEN
        Found = ListGetLogical( Model % BCs(i) % Values, &
           TRIM(VariableName) // ' Rotate',Gotit )
        IF (.NOT. Gotit ) Found = .TRUE.
        IF ( Found ) EXIT
      END IF
    END DO
    IF ( .NOT. Found ) RETURN

    Mesh => Model % Mesh
    n = Mesh % NumberOFNodes

    IF ( .NOT. ASSOCIATED( BoundaryReorder ) ) THEN
      ALLOCATE( BoundaryReorder(n) )
    ELSE IF ( SIZE(BoundaryReorder)<n ) THEN
      DEALLOCATE( BoundaryReorder )
      ALLOCATE( BoundaryReorder(n) )
    END IF
    BoundaryReorder = 0

!------------------------------------------------------------------------------
    DO t=Mesh % NumberOfBulkElements + 1, Mesh % NumberOfBulkElements + &
                  Mesh % NumberOfBoundaryElements

      CurrentElement => Model % Elements(t)
      IF ( CurrentElement % TYPE % ElementCode == 101 )  CYCLE

      n = CurrentElement % TYPE % NumberOfNodes
      NodeIndexes => CurrentElement % NodeIndexes
      ALLOCATE( Condition(n)  )
      DO i=1,Model % NumberOfBCs
        IF ( CurrentElement % BoundaryInfo % Constraint == &
                  Model % BCs(i) % Tag ) THEN
          IF ( ListGetLogical( Model % BCs(i) % Values,VariableName, gotIt) ) THEN
            Found = ListGetLogical( Model % BCs(i) % Values, &
                 TRIM(VariableName) // ' Rotate',gotIt)
            IF ( Found .OR. .NOT. GotIt ) THEN
              Condition(1:n) = ListGetReal( Model % BCs(i) % Values, &
                 TRIM(VariableName) // ' Condition', n, NodeIndexes, Conditional )

              DO j=1,n
                IF ( Conditional .AND. Condition(j)<0._dp ) CYCLE

                k = NodeIndexes(j)
                IF ( BoundaryReorder(k)==0 ) THEN
                  NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                  BoundaryReorder(k) = NumberOfBoundaryNodes
                END IF
              END DO
            END IF
          END IF
        END IF
      END DO
      DEALLOCATE( Condition )
    END DO

    IF (ParEnv % PEs>1 )  THEN
!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
      ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
      n_count = 0
      IF ( NumberOfBoundaryNodes>0 ) THEN
        DO i=1,Mesh % NumberOfNodes
          IF (BoundaryReorder(i)<=0 ) CYCLE
          IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

          nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
          DO j=1,SIZE(nlist)
            k = nlist(j)+1
            IF ( k-1 == ParEnv % myPE ) CYCLE
            n_count(k) = n_count(k)+1
          END DO
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) &
            ALLOCATE( n_index(i) % buff(n_count(i)) )
        END DO
        n_count = 0
        DO i=1,Mesh % NumberOfNodes
          IF (BoundaryReorder(i)<=0 ) CYCLE
          IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

          nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
          DO j=1,SIZE(nlist)
            k = nlist(j)+1
            IF ( k == ParEnv % myPE+1 ) CYCLE
            n_count(k) = n_count(k)+1
            n_index(k) % buff(n_count(k)) = Mesh % Parallelinfo % &
                 GlobalDOFs(i)
          END DO
        END DO
      END IF

      DO i=1,ParEnv % PEs
        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                800, ELMER_COMM_WORLD, ierr )
           IF ( n_count(i)>0 ) &
             CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                 801, ELMER_COMM_WORLD, ierr )
        END IF
      END DO

      DO i=1,ParEnv % PEs
        IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff)

        IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
           CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                800, ELMER_COMM_WORLD, status, ierr )
           IF ( n>0 ) THEN
             ALLOCATE( gbuff(n) )
             proc = status(MPI_SOURCE)
             CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                 801, ELMER_COMM_WORLD, status, ierr )

             DO j=1,n
               k = SearchNodeL( Mesh % ParallelInfo, gbuff(j), Mesh % NumberOfNodes )
               IF ( k>0 ) THEN
                 IF ( BoundaryReorder(k)<= 0 ) THEN
                   NumberOfBoundaryNodes = NumberOfBoundaryNodes + 1
                   BoundaryReorder(k) = NumberOfBoundaryNodes
                 END IF
               END IF
             END DO
             DEALLOCATE(gbuff)
           END IF
        END IF
      END DO
      DEALLOCATE( n_index, n_count )
    END IF

!------------------------------------------------------------------------------

    IF ( NumberOfBoundaryNodes == 0 ) THEN
!     DEALLOCATE( BoundaryReorder )
!     NULLIFY( BoundaryReorder, BoundaryNormals,BoundaryTangent1, &
!                        BoundaryTangent2)
    ELSE
      IF ( ASSOCIATED(BoundaryNormals) ) THEN
        DEALLOCATE( BoundaryNormals, BoundaryTangent1, &
                    BoundaryTangent2, NTelement, NTzeroing_done)
      END IF

      ALLOCATE( NTelement(NumberOfBoundaryNodes,3) )
      ALLOCATE( NTzeroing_done(NumberOfBoundaryNodes,3) )
      ALLOCATE( BoundaryNormals(NumberOfBoundaryNodes,3)  )
      ALLOCATE( BoundaryTangent1(NumberOfBoundaryNodes,3) )
      ALLOCATE( BoundaryTangent2(NumberOfBoundaryNodes,3) )

      BoundaryNormals  = 0.0d0
      BoundaryTangent1 = 0.0d0
      BoundaryTangent2 = 0.0d0
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE CheckNormalTangentialBoundary
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Average boundary normals for nodes. The average boundary normals
!> may be beneficial as they provide more continuous definition of normal
!> over curved boundaries. 
!------------------------------------------------------------------------------
   SUBROUTINE AverageBoundaryNormals( Model, VariableName,    &
       NumberOfBoundaryNodes, BoundaryReorder, BoundaryNormals, &
       BoundaryTangent1, BoundaryTangent2, dim )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model

    INTEGER, POINTER :: BoundaryReorder(:)
    INTEGER :: NumberOfBoundaryNodes,DIM

    REAL(KIND=dp), POINTER :: BoundaryNormals(:,:),BoundaryTangent1(:,:), &
                       BoundaryTangent2(:,:)

    CHARACTER(LEN=*) :: VariableName
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: i,j,k,l,m,n,t, iBC, ierr, proc
    LOGICAL :: GotIt, Found, PeriodicNormals, Conditional
    REAL(KIND=dp) :: s,Bu,Bv,Nrm(3),Basis(32),DetJ
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Matrix_t), POINTER :: Projector
    REAL(KIND=dp), ALLOCATABLE :: Condition(:)

    TYPE(Variable_t), POINTER :: NrmVar, Tan1Var, Tan2Var

    LOGICAL, ALLOCATABLE :: Done(:), NtMasterBC(:), NtSlaveBC(:)
  
    REAL(KIND=dp), POINTER :: SetNormal(:,:), Rot(:,:)

    REAL(KIND=dp), TARGET :: x(Model % MaxElementNodes)
    REAL(KIND=dp), TARGET :: y(Model % MaxElementNodes)
    REAL(KIND=dp), TARGET :: z(Model % MaxElementNodes)

    TYPE buff_t
      INTEGER, ALLOCATABLE :: buff(:)
      REAL(KIND=dp), ALLOCATABLE :: normals(:)
    END TYPE buff_t
    INTEGER, DIMENSION(MPI_STATUS_SIZE) :: status
    INTEGER, POINTER :: nlist(:)
    TYPE(Buff_t), ALLOCATABLE :: n_index(:)
    REAL(KIND=dp), ALLOCATABLE :: nbuff(:)
    INTEGER, ALLOCATABLE :: n_count(:), gbuff(:), n_comp(:)

    LOGICAL :: MassConsistent, LhsSystem, RotationalNormals
    LOGICAL, ALLOCATABLE :: LhsTangent(:),RhsTangent(:)
    INTEGER :: LhsConflicts

    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Origin(3),Axis(3)
    REAL(KIND=dp), POINTER :: Pwrk(:,:)
    LOGICAL :: GotOrigin,GotAxis
    CHARACTER(*), PARAMETER :: Caller = 'AverageBoundaryNormals'

    !------------------------------------------------------------------------------

    ElementNodes % x => x
    ElementNodes % y => y
    ElementNodes % z => z

    Mesh => Model % Mesh
    NrmVar => VariableGet( Mesh % Variables, 'Normals' )

    
    IF ( ASSOCIATED(NrmVar) ) THEN

      IF ( NumberOfBoundaryNodes >0 ) THEN
        BoundaryNormals = 0._dp
        DO i=1,Model % NumberOfNodes
           k = BoundaryReorder(i)
           IF (k>0 ) THEN
             DO l=1,NrmVar % DOFs
                BoundaryNormals(k,l) = NrmVar % Values( NrmVar % DOFs* &
                             (NrmVar % Perm(i)-1)+l)
             END DO
           END IF
         END DO
      END IF

    ELSE

!------------------------------------------------------------------------------
!   Compute sum of elementwise normals for nodes on boundaries
!------------------------------------------------------------------------------
      ALLOCATE( n_comp(Model % NumberOfNodes) )
      n_comp = 0

      IF ( NumberOfBoundaryNodes>0 ) THEN
        BoundaryNormals = 0._dp

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
                      Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE

          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes

          ElementNodes % x(1:n) = Model % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Model % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Model % Nodes % z(NodeIndexes)

          ALLOCATE(Condition(n))

          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              BC => Model % BCs(i) % Values

              IF ( ListGetLogical( BC, VariableName, gotIt) ) THEN
                Found = ListGetLogical( BC, TRIM(VariableName) // ' Rotate',gotIt)
                IF ( Found .OR. .NOT. Gotit ) THEN
                  MassConsistent = ListGetLogical( BC,'Mass Consistent Normals',gotIt)
                  RotationalNormals = ListGetLogical(BC,'Rotational Normals',gotIt)

                  IF( RotationalNormals ) THEN
                    Pwrk => ListGetConstRealArray(BC,'Normals Origin',GotOrigin )
                    IF( GotOrigin ) THEN
                      IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
                        CALL Fatal(Caller,'Size of > Normals Origin < should be 3!')
                      END IF
                      Origin = Pwrk(1:3,1)
                    END IF
                    Pwrk => ListGetConstRealArray(BC,'Normals Axis',GotAxis )
                    IF( GotAxis ) THEN
                      IF( SIZE(Pwrk,1) /= 3 .OR. SIZE(Pwrk,2) /= 1 ) THEN
                        CALL Fatal(Caller,'Size of > Normals Axis < should be 3!')
                      END IF
                      Axis = Pwrk(1:3,1)
                      ! Normalize axis is it should just be used for the direction
                      Axis = Axis / SQRT( SUM( Axis*Axis ) )
                    END IF
                  END IF
                  
                  Condition(1:n) = ListGetReal( BC,&
                       TRIM(VariableName) // ' Condition', n, NodeIndexes, Conditional )

                  DO j=1,n
                    IF ( Conditional .AND. Condition(j) < 0._dp ) CYCLE

                    k = BoundaryReorder( NodeIndexes(j) )
                    IF (k>0) THEN
                      nrm = 0._dp
                      IF (MassConsistent) THEN
                        CALL IntegMassConsistent(j,n,nrm)
                      ELSE IF( RotationalNormals ) THEN
                        nrm(1) = ElementNodes % x(j)
                        nrm(2) = ElementNodes % y(j)
                        nrm(3) = ElementNodes % z(j)

                        IF( GotOrigin ) nrm = nrm - Origin
                        IF( GotAxis ) THEN
                          nrm = nrm - SUM( nrm * Axis ) * Axis
                        ELSE ! Default axis is (0,0,1)
                          nrm(3) = 0.0_dp
                        END IF

                        nrm = nrm / SQRT( SUM( nrm * nrm ) )
                      ELSE
                        Bu = Element % TYPE % NodeU(j)
                        Bv = Element % TYPE % NodeV(j)
                        nrm = NormalVector(Element,ElementNodes,Bu,Bv,.TRUE.)
                      END IF
                      n_comp(NodeIndexes(j)) = 1
                      BoundaryNormals(k,:) = BoundaryNormals(k,:) + nrm
                    END IF
                  END DO
                END IF
              END IF
            END IF
          END DO
          DEALLOCATE(Condition)
        END DO

        DO iBC=1,Model % NumberOfBCs
          Projector => Model % BCs(iBC) % PMatrix
          IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE

          !
          ! TODO: consistent normals, if rotations given:
          ! ---------------------------------------------
          BC => Model % BCs(iBC) % Values
          Rot => ListGetConstRealArray(BC,'Periodic BC Rotate', Found )
          IF ( Found .AND. ASSOCIATED(Rot) ) THEN
            IF ( ANY(Rot/=0) ) THEN
              ALLOCATE( Done(SIZE(BoundaryNormals,1)) )
              Done=.FALSE.
              DO i=1,Projector % NumberOfRows
                 k = BoundaryReorder(Projector % InvPerm(i))
                 IF ( k <= 0 ) CYCLE
                 DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                   IF ( Projector % Cols(l) <= 0 ) CYCLE
                   m = BoundaryReorder(Projector % Cols(l))
                   IF ( m>0 ) THEN
                     IF ( .NOT.Done(m) ) THEN
                       Done(m) = .TRUE.
                       BoundaryNormals(m,:) = -BoundaryNormals(m,:)
                     END IF
                   END IF
                 END DO
              END DO
              DEALLOCATE(Done)
              CYCLE
            END IF
          END IF

          DO i=1,Projector % NumberOfRows
            k = BoundaryReorder(Projector % InvPerm(i))
            IF ( k <= 0 ) CYCLE
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 ) CYCLE
              m = BoundaryReorder(Projector % Cols(l))
              IF ( m>0 ) BoundaryNormals(m,:) = 0._dp
            END DO
          END DO
        END DO

        DO iBC=1,Model % NumberOfBCs
           Projector => Model % BCs(iBC) % PMatrix
           IF ( .NOT. ASSOCIATED( Projector ) ) CYCLE

           !
           ! TODO: consistent normals, if rotations given:
           ! ---------------------------------------------
           BC => Model % BCs(iBC) % Values
           Rot => ListGetConstRealArray(BC,'Periodic BC Rotate', Found )
           IF ( Found .AND. ASSOCIATED(Rot) ) THEN
             IF ( ANY(Rot/=0) ) CYCLE
           END IF

           DO i=1,Projector % NumberOfRows
              k = BoundaryReorder(Projector % InvPerm(i))
              IF ( k <= 0 ) CYCLE
              DO l=Projector % Rows(i),Projector % Rows(i+1)-1
                IF ( Projector % Cols(l) <= 0 ) CYCLE
                m = BoundaryReorder(Projector % Cols(l))
                IF ( m > 0 ) &
                   BoundaryNormals(m,:) = BoundaryNormals(m,:) + &
                     Projector % Values(l) * BoundaryNormals(k,:)
              END DO
           END DO
        END DO
      END IF

      IF (ParEnv % PEs>1 ) THEN
        ALLOCATE( n_count(ParEnv% PEs),n_index(ParEnv % PEs) )
        n_count = 0

        IF ( NumberOfBoundaryNodes>0 ) THEN
          DO i=1,Mesh % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE
  
            nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
              k = nlist(j)+1
              IF ( k-1 == ParEnv % myPE ) CYCLE
              n_count(k) = n_count(k)+1
            END DO
          END DO
          DO i=1,ParEnv % PEs
            IF ( n_count(i)>0 ) &
                ALLOCATE( n_index(i) % buff(n_count(i)), &
                        n_index(i) % normals(3*n_count(i)) )
          END DO

          n_count = 0
          DO i=1,Model % NumberOfNodes
            IF (BoundaryReorder(i)<=0 .OR. n_comp(i)<=0 ) CYCLE
            IF (.NOT.Mesh % ParallelInfo % INTERFACE(i) ) CYCLE

            nlist => Mesh % ParallelInfo % NeighbourList(i) % Neighbours
            DO j=1,SIZE(nlist)
              k = nlist(j)+1
              IF ( k-1 == ParEnv % myPE ) CYCLE
              n_count(k) = n_count(k)+1
              n_index(k) % buff(n_count(k)) = Mesh % Parallelinfo % &
                 GlobalDOFs(i)
              l = BoundaryReorder(i)
              n_index(k) % normals(3*n_count(k)-2)=BoundaryNormals(l,1)
              n_index(k) % normals(3*n_count(k)-1)=BoundaryNormals(l,2)
              n_index(k) % normals(3*n_count(k)-0)=BoundaryNormals(l,3)
            END DO
          END DO
        END IF

        DO i=1,ParEnv % PEs
          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
            CALL MPI_BSEND( n_count(i), 1, MPI_INTEGER, i-1, &
                900, ELMER_COMM_WORLD, ierr )
            IF ( n_count(i)>0 ) THEN
              CALL MPI_BSEND( n_index(i) % buff, n_count(i), MPI_INTEGER, i-1, &
                  901, ELMER_COMM_WORLD, ierr )
              CALL MPI_BSEND( n_index(i) % normals, 3*n_count(i), MPI_DOUBLE_PRECISION, &
                    i-1,  902, ELMER_COMM_WORLD, ierr )
            END IF
          END IF
        END DO
        DO i=1,ParEnv % PEs
          IF ( n_count(i)>0 ) DEALLOCATE( n_index(i) % Buff, n_index(i) % Normals)

          IF ( ParEnv % Active(i) .AND. ParEnv % IsNeighbour(i) ) THEN
             CALL MPI_RECV( n, 1, MPI_INTEGER, MPI_ANY_SOURCE, &
                    900, ELMER_COMM_WORLD, status, ierr )
             IF ( n>0 ) THEN
               proc = status(MPI_SOURCE)
               ALLOCATE( gbuff(n), nbuff(3*n) )
               CALL MPI_RECV( gbuff, n, MPI_INTEGER, proc, &
                   901, ELMER_COMM_WORLD, status, ierr )

               CALL MPI_RECV( nbuff, 3*n, MPI_DOUBLE_PRECISION, proc, &
                    902, ELMER_COMM_WORLD, status, ierr )

               DO j=1,n
                 k = SearchNodeL( Mesh % ParallelInfo, gbuff(j), Mesh % NumberOfNodes )
                 IF ( k>0 ) THEN
                   n_comp(k) = n_comp(k)+1
                   l = BoundaryReorder(k)
                   IF ( l>0 ) THEN
                     BoundaryNormals(l,1)=BoundaryNormals(l,1)+nbuff(3*j-2)
                     BoundaryNormals(l,2)=BoundaryNormals(l,2)+nbuff(3*j-1)
                     BoundaryNormals(l,3)=BoundaryNormals(l,3)+nbuff(3*j-0)
                   END IF
                 END IF
               END DO
               DEALLOCATE(gbuff, nbuff)
             END IF
          END IF
        END DO
        DEALLOCATE( n_index, n_count )
      END IF

      DEALLOCATE(n_comp)
    END IF

!------------------------------------------------------------------------------
!   normalize 
!------------------------------------------------------------------------------
    IF ( NumberOfBoundaryNodes>0 ) THEN

      LhsSystem = ListGetLogical(Model % Simulation,'Use Lhs System',Found) 
      IF(.NOT. Found ) LhsSystem = ( dim == 3 )

      IF( LhsSystem ) THEN
        ALLOCATE( NtMasterBC( Model % NumberOfBCs ), NtSlaveBC( Model % NumberOfBCs ) )
        NtMasterBC = .FALSE.; NtSlaveBC = .FALSE.

        DO i = 1, Model % NumberOfBcs
          IF( .NOT. ListCheckPrefix( Model % BCs(i) % Values,'Normal-Tangential') ) CYCLE
          
          j = ListGetInteger( Model % BCs(i) % Values,'Mortar BC',Found )
          IF( .NOT. Found ) THEN
            j = ListGetInteger( Model % BCs(i) % Values,'Contact BC',Found )
          END IF
          IF( j == 0 .OR. j > Model % NumberOfBCs ) CYCLE

          NtSlaveBC( i ) = .TRUE.
          NtMasterBC( j ) = .TRUE.
        END DO
        LhsSystem = ANY( NtMasterBC )
      END IF

      IF( LhsSystem ) THEN
        DO i = 1, Model % NumberOfBcs
          IF( NtSlaveBC( i ) .AND. NtMasterBC( i ) ) THEN
            CALL Warn('AverageBoundaryNormals','BC '//TRIM(I2S(i))//' is both N-T master and slave!')
          END IF
        END DO

        ALLOCATE( LhsTangent( Model % NumberOfNodes ) )
        LhsTangent = .FALSE.

        ALLOCATE( RhsTangent( Model % NumberOfNodes ) )
        RhsTangent = .FALSE. 

        DO t=Model % NumberOfBulkElements + 1, Model % NumberOfBulkElements + &
            Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          IF ( Element % TYPE  % ElementCode < 200 ) CYCLE
          
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          
          DO i=1,Model % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
              IF( NtMasterBC(i) ) LhsTangent( NodeIndexes ) = .TRUE.
              IF( NtSlaveBC(i) ) RhsTangent( NodeIndexes ) = .TRUE.
              EXIT
            END IF
          END DO
        END DO

        LhsConflicts = COUNT( LhsTangent .AND. RhsTangent )
        IF( LhsConflicts > 0 ) THEN
          CALL Warn('AverageBoundaryNormals',&
              'There are '//TRIM(I2S(LhsConflicts))//' nodes that could be both rhs and lhs!')
        END IF
      END IF


      DO i=1,Model % NumberOfNodes
        k = BoundaryReorder(i) 
        IF ( k > 0 ) THEN
          s = SQRT( SUM( BoundaryNormals(k,:)**2 ) )
          IF ( s /= 0.0d0 ) &
            BoundaryNormals(k,:) = BoundaryNormals(k,:) / s
          IF ( dim > 2 ) THEN
            CALL TangentDirections( BoundaryNormals(k,:),  &
                BoundaryTangent1(k,:), BoundaryTangent2(k,:) )
            IF( LhsSystem ) THEN
              IF( LhsTangent(i) ) THEN
                BoundaryTangent2(k,:) = -BoundaryTangent2(k,:)
              END IF
            END IF
          END IF
        END IF
      END DO
      
      IF( ListGetLogical( Model % Simulation,'Save Averaged Normals',Found ) ) THEN
        CALL Info('AverageBoundaryNormals','Saving averaged boundary normals to variable: Averaged Normals')
        NrmVar => VariableGet( Mesh % Variables, 'Averaged Normals' )
        
        IF(.NOT. ASSOCIATED( NrmVar ) ) THEN
          CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,'Averaged Normals',3,&
              Perm = BoundaryReorder )
          NrmVar => VariableGet( Mesh % Variables, 'Averaged Normals' )
        END IF
            
        DO i=1,Model % NumberOfNodes
          k = BoundaryReorder(i)
          IF (k>0 ) THEN
            DO l=1,NrmVar % DOFs
              NrmVar % Values( NrmVar % DOFs* &
                  (NrmVar % Perm(i)-1)+l)  = BoundaryNormals(k,l)
            END DO
          END IF
        END DO

        IF( dim > 2 .AND. ListGetLogical( Model % Simulation,'Save Averaged Tangents',Found ) ) THEN
          Tan1Var => VariableGet( Mesh % Variables, 'Averaged First Tangent' )
          Tan2Var => VariableGet( Mesh % Variables, 'Averaged Second Tangent' )

          IF(.NOT. ASSOCIATED( Tan1Var ) ) THEN
            CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,&
                'Averaged First Tangent',3, Perm = BoundaryReorder )
            Tan1Var => VariableGet( Mesh % Variables, 'Averaged First Tangent' )
            CALL VariableAddVector( Mesh % Variables, Mesh, Model % Solver,&
                'Averaged Second Tangent',3, Perm = BoundaryReorder )
            Tan2Var => VariableGet( Mesh % Variables, 'Averaged Second Tangent' )
          END IF
          
          DO i=1,Model % NumberOfNodes
            k = BoundaryReorder(i)
            IF (k>0 ) THEN
              DO l=1,Tan1Var % DOFs
                Tan1Var % Values( Tan1Var % DOFs* &
                    (Tan1Var % Perm(i)-1)+l)  = BoundaryTangent1(k,l)
                Tan2Var % Values( Tan2Var % DOFs* &
                    (Tan2Var % Perm(i)-1)+l)  = BoundaryTangent2(k,l)
              END DO
            END IF
          END DO
        END IF
      END IF
    END IF



 CONTAINS

    SUBROUTINE IntegMassConsistent(j,n,nrm)
      INTEGER :: t,j,n
      LOGICAL :: stat
      REAL(KIND=dp) :: detJ,Basis(n),nrm(:),lnrm(3)

      TYPE(GaussIntegrationPoints_t) :: IP

      !----------------------
      IP = GaussPoints(Element)
      DO t=1,IP % n
        stat = ElementInfo(Element, ElementNodes, IP % U(t), &
               IP % v(t), IP % W(t), detJ, Basis)

        lnrm = NormalVector(Element,ElementNodes, &
              IP % U(t),IP % v(t),.TRUE.)

        nrm = nrm + IP % s(t) * lnrm * detJ * Basis(j)
      END DO
    END SUBROUTINE IntegMassConsistent

!------------------------------------------------------------------------------
  END SUBROUTINE AverageBoundaryNormals
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Search an element QueriedNode from an ordered set Nodes and return
!> Index to Nodes structure. Return value -1 means QueriedNode was
!> not found.
!------------------------------------------------------------------------------
FUNCTION SearchNodeL( ParallelInfo, QueriedNode,n ) RESULT(Indx)

  USE Types
  IMPLICIT NONE

  TYPE (ParallelInfo_t) :: ParallelInfo
  INTEGER :: QueriedNode, Indx,n

  ! Local variables

  INTEGER :: Lower, Upper, Lou, i

!------------------------------------------------------------------------------

  Indx = -1
  Upper = n
  Lower = 1

  ! Handle the special case

  IF ( Upper == 0 ) RETURN

10 CONTINUE
  IF ( ParallelInfo % GlobalDOFs(Lower) == QueriedNode ) THEN
     Indx = Lower
     RETURN
  ELSE IF ( ParallelInfo % GlobalDOFs(Upper) == QueriedNode ) THEN
     Indx = Upper
     RETURN
  END IF

  IF ( (Upper - Lower) > 1 ) THEN
     Lou = ISHFT((Upper + Lower), -1)
     IF ( ParallelInfo % GlobalDOFs(Lou) < QueriedNode ) THEN
        Lower = Lou
        GOTO 10
     ELSE
        Upper = Lou
        GOTO 10
     END IF
  END IF

  RETURN
!------------------------------------------------------------------------------
END FUNCTION SearchNodeL
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Initialize solver for next timestep.
!------------------------------------------------------------------------------
  SUBROUTINE InitializeTimestep( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver  !< Solver to be initialized.
!------------------------------------------------------------------------------
     CHARACTER(LEN=MAX_NAME_LEN) :: Method
     LOGICAL :: GotIt
     INTEGER :: i, Order,ndofs
     REAL(KIND=dp), POINTER CONTIG :: SaveValues(:)
     TYPE(Matrix_t), POINTER :: A
     TYPE(Variable_t), POINTER :: Var
     
!------------------------------------------------------------------------------
     Solver % DoneTime = Solver % DoneTime + 1
!------------------------------------------------------------------------------

     IF ( .NOT. ASSOCIATED( Solver % Matrix ) .OR. &
          .NOT. ASSOCIATED( Solver % Variable % Values ) ) RETURN

     IF ( Solver % TimeOrder <= 0 ) RETURN
!------------------------------------------------------------------------------

     Method = ListGetString( Solver % Values, 'Timestepping Method', GotIt )
     IF ( Method == 'none' ) RETURN
    
     IF ( .NOT.GotIt ) THEN

       Solver % Beta = ListGetConstReal( Solver % Values, 'Newmark Beta', GotIt )
       IF ( .NOT. GotIt ) THEN
         Solver % Beta = ListGetConstReal( CurrentModel % Simulation, 'Newmark Beta', GotIt )
       END IF

       IF ( .NOT.GotIt ) THEN
         IF (Solver % TimeOrder > 1) THEN
           Method = 'bossak'
           Solver % Beta = 1.0d0
         ELSE
           CALL Warn( 'InitializeTimestep', &
               'Timestepping method defaulted to IMPLICIT EULER' )

           Solver % Beta = 1.0D0
           Method = 'implicit euler'
         END IF
       END IF

     ELSE

       Solver % Beta = 1._dp
       SELECT CASE( Method )
         CASE('implicit euler')
           Solver % Beta = 1.0d0

         CASE('explicit euler')
           Solver % Beta = 0.0d0

         CASE('runge-kutta')
           Solver % Beta = 0.0d0

         CASE('crank-nicolson')
           Solver % Beta = 0.5d0

         CASE('fs')
           Solver % Beta = 0.5d0

         CASE('adams-bashforth')
           Solver % Beta = 0.0d0

         CASE('adams-moulton')
           Solver % Beta = 1.0d0

         CASE('newmark')
           Solver % Beta = ListGetConstReal( Solver % Values, 'Newmark Beta', GotIt )
           IF ( .NOT. GotIt ) THEN
              Solver % Beta = ListGetConstReal( CurrentModel % Simulation, &
                              'Newmark Beta', GotIt )
           END IF

           IF ( Solver % Beta<0 .OR. Solver % Beta>1 ) THEN
             WRITE( Message, * ) 'Invalid value of Beta ', Solver % Beta
             CALL Warn( 'InitializeTimestep', Message )
           END IF

         CASE('bdf')
           IF ( Solver % Order < 1 .OR. Solver % Order > 5  ) THEN
             WRITE( Message, * ) 'Invalid order BDF ',  Solver % Order
             CALL Fatal( 'InitializeTimestep', Message )
           END IF

         CASE('bossak')
           Solver % Beta = 1.0d0

         CASE DEFAULT 
           WRITE( Message, * ) 'Unknown timestepping method: ',Method
           CALL Fatal( 'InitializeTimestep', Message )
       END SELECT

     END IF

     ndofs = Solver % Matrix % NumberOfRows
     Var => Solver % Variable
     
     IF ( Method /= 'bdf' .OR. Solver % TimeOrder > 1 ) THEN

       IF ( Solver % DoneTime == 1 .AND. Solver % Beta /= 0.0d0 ) THEN
         Solver % Beta = 1.0d0
       END IF
       IF( Solver % TimeOrder == 2 ) THEN         
         Solver % Alpha = ListGetConstReal( Solver % Values, &
             'Bossak Alpha', GotIt )
         IF ( .NOT. GotIt ) THEN
           Solver % Alpha = ListGetConstReal( CurrentModel % Simulation, &
               'Bossak Alpha', GotIt )
         END IF
         IF ( .NOT. GotIt ) Solver % Alpha = -0.05d0
       END IF
       
       SELECT CASE( Solver % TimeOrder )
         
       CASE(1)
         Order = MIN(Solver % DoneTime, Solver % Order)
         DO i=Order, 2, -1
           Var % PrevValues(:,i) = Var % PrevValues(:,i-1)
         END DO
         Var % PrevValues(:,1) = Var % Values
         Solver % Matrix % Force(:,2) = Solver % Matrix % Force(:,1)
         
       CASE(2)
         Var % PrevValues(:,3) = Var % Values
         Var % PrevValues(:,4) = Var % PrevValues(:,1)
         Var % PrevValues(:,5) = Var % PrevValues(:,2)
       END SELECT
     ELSE
       Order = MIN(Solver % DoneTime, Solver % Order)
       DO i=Order, 2, -1
         Var % PrevValues(:,i) = Var % PrevValues(:,i-1)
       END DO
       Var % PrevValues(:,1) = Var % Values
     END IF


     IF( ListGetLogical( Solver % Values,'Nonlinear Timestepping', GotIt ) ) THEN
       IF( Solver % DoneTime > 1 ) THEN
         A => Solver % Matrix
         CALL Info('InitializeTimestep','Saving previous linear system for timestepping',Level=12)
         IF( .NOT. ASSOCIATED( A % BulkValues ) ) THEN
           CALL Fatal('InitializeTimestep','BulkValues should be associated!')
         END IF
         
         IF( .NOT. ASSOCIATED( A % BulkResidual ) ) THEN
           ALLOCATE( A % BulkResidual( SIZE( A % BulkRhs ) ) )
         END IF
         
         SaveValues => A % Values
         A % Values => A % BulkValues
         CALL MatrixVectorMultiply( A, Var % Values, A % BulkResidual )
         A % Values => SaveValues
         A % BulkResidual = A % BulkResidual - A % BulkRhs
       END IF
     END IF


     ! Advance also the exported variables if they happen to be time-dependent
     ! They only have normal prevvalues, when writing this always 2. 
     BLOCK
       INTEGER :: VarNo,n
       CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name
       LOGICAL :: Found
       
       VarNo =0      
       DO WHILE( .TRUE. )
         VarNo = VarNo + 1         
         str = ComponentName( 'exported variable', VarNo )    
         
         var_name = ListGetString( Solver % Values, str, Found )    
         IF(.NOT. Found) EXIT
         
         CALL VariableNameParser( var_name ) 
         
         Var => VariableGet( Solver % Mesh % Variables, Var_name )
         IF( .NOT. ASSOCIATED(Var)) CYCLE
         IF( .NOT. ASSOCIATED(Var % PrevValues) ) CYCLE
         
         n = SIZE( Var % PrevValues,2 )
         DO i=n,2,-1
           Var % PrevValues(:,i) = Var % PrevValues(:,i-1)
         END DO
         Var % PrevValues(:,1) = Var % Values
       END DO
     END BLOCK
     
!------------------------------------------------------------------------------
  END SUBROUTINE InitializeTimestep
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Update force vector AFTER ALL OTHER ASSEMBLY STEPS BUT BEFORE SETTING
!> DIRICHLET CONDITIONS. Required only for time dependent simulations..
!------------------------------------------------------------------------------
  SUBROUTINE FinishAssembly( Solver, ForceVector )
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: ForceVector(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Method, Simulation
    INTEGER :: Order
    LOGICAL :: Found
!------------------------------------------------------------------------------

    IF ( Solver % Matrix % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(Solver % Matrix)
    END IF

    Simulation = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
    IF ( Simulation == 'transient' ) THEN
      Method = ListGetString( Solver % Values, 'Timestepping Method' )
      Order = MIN(Solver % DoneTime, Solver % Order)

      IF ( Order <= 0 .OR. Solver % TimeOrder /= 1 .OR. Method=='bdf' ) RETURN

      IF ( Solver % Beta /= 0.0d0 ) THEN
        ForceVector = ForceVector + ( Solver % Beta - 1 ) * &
            Solver % Matrix % Force(:,1) + &
                ( 1 - Solver % Beta ) * Solver % Matrix % Force(:,2)
      END IF
    END IF

!------------------------------------------------------------------------------
  END SUBROUTINE FinishAssembly
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE InvalidateVariable( TopMesh,PrimaryMesh,Name )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Name
    TYPE(Mesh_t),  POINTER :: TopMesh,PrimaryMesh
!------------------------------------------------------------------------------
    CHARACTER(LEN=MAX_NAME_LEN) :: tmpname
    INTEGER :: i
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Variable_t), POINTER :: Var,Var1, PrimVar
!------------------------------------------------------------------------------
    Mesh => TopMesh

    PrimVar => VariableGet( PrimaryMesh % Variables, Name, ThisOnly=.TRUE.)
    IF ( .NOT.ASSOCIATED( PrimVar) ) RETURN

    DO WHILE( ASSOCIATED(Mesh) )
      ! Make the same variable invalid in all other meshes.
      IF ( .NOT.ASSOCIATED( PrimaryMesh, Mesh) ) THEN
        Var => VariableGet( Mesh % Variables, Name, ThisOnly=.TRUE.)
        IF ( ASSOCIATED( Var ) ) THEN
          Var % Valid = .FALSE.
          Var % PrimaryMesh => PrimaryMesh
        END IF

        IF ( PrimVar % DOFs > 1 ) THEN
          DO i=1,PrimVar % DOFs
            tmpname = ComponentName( Name, i )
            Var1 => VariableGet( Mesh % Variables, tmpname, .TRUE. )
            IF ( ASSOCIATED( Var1 ) ) THEN
              Var1 % Valid = .FALSE.
              Var1 % PrimaryMesh => PrimaryMesh
            END IF
          END DO
        END IF
      END IF
      Mesh => Mesh % Next
    END DO 

    ! Tell that values have changed in the primary mesh.
    ! Interpolation can then be activated if we request the same variable in the
    ! other meshes. 
    PrimVar % ValuesChanged = .TRUE.
    IF ( PrimVar % DOFs > 1 ) THEN
      DO i=1,PrimVar % DOFs
        tmpname = ComponentName( Name, i )
        Var => VariableGet( PrimaryMesh % Variables, tmpname, .TRUE. )
        IF ( ASSOCIATED(Var) ) Var % ValuesChanged = .TRUE.
      END DO
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE InvalidateVariable
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Computes the norm related to a solution vector of the Solver.
!------------------------------------------------------------------------------
  FUNCTION ComputeNorm(Solver, nin, values) RESULT (Norm)
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t), TARGET :: Solver
    INTEGER :: nin
    REAL(KIND=dp), TARGET, OPTIONAL :: values(:)
    
    INTEGER :: NormDim, NormDofs, Dofs,i,j,k,n,totn,PermStart
    INTEGER, POINTER :: NormComponents(:)
    INTEGER, ALLOCATABLE :: iPerm(:)
    REAL(KIND=dp) :: Norm, nscale, val
    LOGICAL :: Stat, ComponentsAllocated, ConsistentNorm
    REAL(KIND=dp), POINTER :: x(:)
    REAL(KIND=dp), ALLOCATABLE, TARGET :: y(:)

    CALL Info('ComputeNorm','Computing norm of solution',Level=10)

    IF(PRESENT(values)) THEN
      x => values
    ELSE
      x => Solver % Variable % Values
    END IF
    

    NormDim = ListGetInteger(Solver % Values,'Nonlinear System Norm Degree',Stat)
    IF(.NOT. Stat) NormDim = 2

    Dofs = Solver % Variable % Dofs

    ComponentsAllocated = .FALSE.
    NormComponents => ListGetIntegerArray(Solver % Values,&
        'Nonlinear System Norm Components',Stat)
    IF(Stat) THEN
      NormDofs = SIZE( NormComponents ) 
    ELSE
      NormDofs = ListGetInteger(Solver % Values,'Nonlinear System Norm Dofs',Stat)
      IF(Stat) THEN
        ALLOCATE(NormComponents(NormDofs))
        ComponentsAllocated = .TRUE.
        DO i=1,NormDofs
          NormComponents(i) = i
        END DO
      ELSE
        NormDofs = Dofs        
      END IF
    END IF
 
    n = nin
    totn = 0

    IF( ParEnv % PEs > 1 ) THEN
      ConsistentNorm = ListGetLogical(Solver % Values,'Nonlinear System Consistent Norm',Stat)
      IF (ConsistentNorm) CALL Info('ComputeNorm','Using consistent norm in parallel',Level=10)
    ELSE
      ConsistentNorm = .FALSE.
    END IF


    PermStart = ListGetInteger(Solver % Values,'Norm Permutation',Stat)
    IF ( Stat ) THEN
      ALLOCATE(iPerm(SIZE(Solver % Variable % Perm))); iPerm=0
      n = 0
      DO i=PermStart,SIZE(iPerm)
        IF ( Solver % Variable % Perm(i)>0 ) THEN
          n = n + 1
          iPerm(n) = Solver % Variable % Perm(i)
        END IF
      END DO
      ALLOCATE(y(n))
      y = x(iPerm(1:n))
      x => y
      DEALLOCATE(iPerm)
    END IF


    IF( NormDofs < Dofs ) THEN
      IF( ConsistentNorm ) THEN
        CALL Warn('ComputeNorm','Consistent norm not implemented for selective norm')
      END IF

      totn = NINT( ParallelReduction(1._dp*n) )
      nscale = NormDOFs*totn/(1._dp*DOFs)
      Norm = 0.0_dp

      SELECT CASE(NormDim)
      CASE(0)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = MAX(Norm, MAXVAL( ABS(x(j::Dofs))) )
        END DO
        Norm = ParallelReduction(Norm,2)
      CASE(1)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( ABS(x(j::Dofs)) )
        END DO
        Norm = ParallelReduction(Norm)/nscale
      CASE(2)
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( x(j::Dofs)**2 )
        END DO
        Norm = SQRT(ParallelReduction(Norm)/nscale)
      CASE DEFAULT
        DO i=1,NormDofs
          j = NormComponents(i)
          Norm = Norm + SUM( x(j::Dofs)**NormDim )
        END DO
        Norm = (ParallelReduction(Norm)/nscale)**(1.0d0/NormDim)
      END SELECT
    ELSE IF( ConsistentNorm ) THEN
      ! In consistent norm we have to skip the dofs not owned by the partition in order
      ! to count each dof only once. 

      Norm = 0.0_dp
      totn = 0
      DO j=1,n
        IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
            == ParEnv % MyPE ) totn = totn + 1
      END DO        

      totn = NINT( ParallelReduction(1._dp*totn) )
      nscale = 1.0_dp * totn

      SELECT CASE(NormDim)
        
      CASE(0) 
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = MAX( Norm, ABS( val ) )
        END DO
        
      CASE(1)
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + ABS(val)
        END DO
        
      CASE(2)          
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          
          Norm = Norm + val**2 
        END DO
        
      CASE DEFAULT
        DO j=1,n
          IF( Solver % Matrix % ParallelInfo % NeighbourList(j) % Neighbours(1) &
              /= ParEnv % MyPE ) CYCLE
          val = x(j)
          Norm = Norm + val**NormDim 
        END DO
      END SELECT

      SELECT CASE(NormDim)
      CASE(0)
        Norm = ParallelReduction(Norm,2)
      CASE(1)
        Norm = ParallelReduction(Norm) / nscale
      CASE(2)
        Norm = SQRT(ParallelReduction(Norm)/nscale)
      CASE DEFAULT
        Norm = (ParallelReduction(Norm)/nscale)**(1.0d0/NormDim)
      END SELECT
      
    ELSE      
      val = ParallelReduction(1.0_dp*n)
      totn = NINT( val )
      IF (totn == 0) THEN
         CALL Warn('ComputeNorm','Requested norm of a variable with no Dofs')
         Norm = 0.0_dp
      ELSE
         nscale = 1.0_dp * totn
         
         val = 0.0_dp
         SELECT CASE(NormDim)
         CASE(0)
            IF (n>0) val = MAXVAL(ABS(x(1:n)))
            Norm = ParallelReduction(val,2)
         CASE(1)
            IF (n>0) val = SUM(ABS(x(1:n)))
            Norm = ParallelReduction(val)/nscale
         CASE(2)
            IF (n>0) val = SUM(x(1:n)**2)
            Norm = SQRT(ParallelReduction(val)/nscale)
         CASE DEFAULT
            IF (n>0) val = SUM(x(1:n)**NormDim)
            Norm = (ParallelReduction(val)/nscale)**(1.0d0/NormDim)
         END SELECT
      END IF
    END IF

!   PRINT *,'ComputedNorm:',Norm, NormDIm
    
    IF( ComponentsAllocated ) THEN
      DEALLOCATE( NormComponents ) 
    END IF
!------------------------------------------------------------------------------
  END FUNCTION ComputeNorm
!------------------------------------------------------------------------------


  SUBROUTINE UpdateDependentObjects( Solver, SteadyState )

    TYPE(Solver_t), TARGET :: Solver
    LOGICAL :: SteadyState

    TYPE(ValueList_t), POINTER :: SolverParams
    LOGICAL :: Found, DoIt
    REAL(KIND=dp) :: dt
    TYPE(Variable_t), POINTER :: dtVar, VeloVar
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER, POINTER :: UpdateComponents(:)
    CHARACTER(*), PARAMETER :: Caller = 'UpdateDependentObjects'
  
    SolverParams => Solver % Values
    
    IF( SteadyState ) THEN
      CALL Info(Caller,'Updating objects depending on primary field in steady state',Level=20)
    ELSE
      CALL Info(Caller,'Updating objects depending on primary field in nonlinear system',Level=20)
    END IF
    
    
    ! The update of exported variables on nonlinear or steady state level.
    ! In nonlinear level the nonlinear iteration may depend on the updated values.
    ! Steady-state level is often sufficient if the dependendence is on some other solver.
    !-----------------------------------------------------------------------------------------    
    IF( SteadyState ) THEN
      DoIt = ListGetLogical( SolverParams,&
          'Update Exported Variables', Found )
    ELSE        
      DoIt = ListGetLogical( SolverParams,&
          'Nonlinear Update Exported Variables',Found )
    END IF
    IF( DoIt ) THEN
      CALL Info(Caller,'Updating exported variables',Level=20)
      CALL UpdateExportedVariables( Solver )	
    END IF
       
    ! Update components that depende on the solution of the solver.
    ! Nonlinear level allows some nonlinear couplings within the solver. 
    !-----------------------------------------------------------------------------------------
    IF( SteadyState ) THEN
      UpdateComponents => ListGetIntegerArray( SolverParams, &
          'Update Components', DoIt )
    ELSE
      UpdateComponents => ListGetIntegerArray( SolverParams, &
          'Nonlinear Update Components', DoIt )
    END IF
    IF( DoIt ) THEN
      CALL Info(Caller,'Updating components',Level=20)
      CALL UpdateDependentComponents( UpdateComponents )	
    END IF

    ! Compute derivative of solution with time i.e. velocity 
    ! For 2nd order schemes there is direct pointer to the velocity component
    ! Thus only 1st order schemes need to be computed.
    !-----------------------------------------------------------------------------------------
    DoIt = .FALSE.
    IF( SteadyState ) THEN
      DoIt = ListGetLogical( SolverParams,'Calculate Velocity',Found )
    ELSE
      DoIt = ListGetLogical( SolverParams,'Nonlinear Calculate Velocity',Found )
    END IF
    
    IF( DoIt ) THEN
      CALL Info(Caller,'Updating variable velocity')
      IF( .NOT. ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        CALL Warn(Caller,'Cannot calculate velocity without previous values!')
      ELSE IF( Solver % TimeOrder == 1) THEN
        dtVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
        dt = dtVar % Values(1) 
        str = TRIM( Solver % Variable % Name ) // ' Velocity'
        VeloVar => VariableGet( Solver % Mesh % Variables, str )        
        VeloVar % Values = (Solver % Variable % Values - Solver % Variable % PrevValues(:,1)) / dt
      END IF
    END IF
    
    ! Finally compute potentially velocities related to exported variables.
    ! Do this on nonlinear level only when 'Nonlinear Calculate Velocity' is set true.
    !-----------------------------------------------------------------------------------------
    IF( SteadyState .OR. DoIt ) THEN          
      CALL DerivateExportedVariables( Solver )  
    END IF       
    
  END SUBROUTINE UpdateDependentObjects


  
!------------------------------------------------------------------------------
!> When a new field has been computed compare it to the previous one.
!> Different convergence measures may be used. 
!> Also performs relaxation if a non-unity relaxation factor is given.
!------------------------------------------------------------------------------
  SUBROUTINE ComputeChange(Solver,SteadyState,nsize,values,values0,Matrix,RHS)
!------------------------------------------------------------------------------
    TYPE(Solver_t), TARGET :: Solver
    LOGICAL :: SteadyState
    TYPE(Matrix_t), OPTIONAL, TARGET :: Matrix
    INTEGER, OPTIONAL :: nsize
    REAL(KIND=dp), OPTIONAL, TARGET :: values(:), values0(:), RHS(:)
!------------------------------------------------------------------------------
    INTEGER :: i, n, nn, RelaxAfter, IterNo, MinIter, MaxIter, dofs
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: b(:), x(:), r(:)
    REAL(KIND=dp), POINTER :: x0(:)
    REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Change, PrevChange, Relaxation, tmp(1),dt, &
        Tolerance, MaxNorm, eps, Ctarget, Poffset, nsum, dpsum
    CHARACTER(LEN=MAX_NAME_LEN) :: ConvergenceType
    INTEGER, TARGET  ::  Dnodes(1)
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: iterVar, VeloVar, dtVar, WeightVar
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, str
    LOGICAL :: Stat, ConvergenceAbsolute, Relax, RelaxBefore, DoIt, Skip, &
        SkipConstraints, ResidualMode, RelativeP
    TYPE(Matrix_t), POINTER :: MMatrix
    REAL(KIND=dp), POINTER CONTIG :: Mx(:), Mb(:), Mr(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec
    INTEGER :: ipar(1)
    TYPE(ValueList_t), POINTER :: SolverParams
    CHARACTER(*), PARAMETER :: Caller = 'ComputeChange'

    
    SolverParams => Solver % Values
    RelativeP = .FALSE.
    
    IF(SteadyState) THEN	
      Skip = ListGetLogical( SolverParams,'Skip Compute Steady State Change',Stat)
      IF( Skip ) THEN
        CALL Info(Caller,'Skipping the computation of steady state change',Level=15)
        RETURN
      END IF
        
      ! No residual mode for steady state analysis
      ResidualMode = .FALSE.

      ConvergenceType = ListGetString(SolverParams,&
          'Steady State Convergence Measure',Stat)
      IF(.NOT. Stat) ConvergenceType = 'norm' 

      ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Steady State Convergence Absolute',Stat)
      IF(.NOT. Stat) ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Use Absolute Norm for Convergence',Stat)

      Relaxation = ListGetCReal( SolverParams, &
          'Steady State Relaxation Factor', Relax )
      Relax = Relax .AND. ABS(Relaxation-1.0_dp) > EPSILON(Relaxation)

      iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
      IterNo = NINT( iterVar % Values(1) )
      IF( Relax ) THEN
        RelaxAfter = ListGetInteger(SolverParams,'Steady State Relaxation After',Stat)
        IF( Stat .AND. RelaxAfter >= IterNo ) Relax = .FALSE.
      END IF	

      RelaxBefore = .TRUE.
      IF(Relax) THEN
        RelaxBefore = ListGetLogical( SolverParams, &
            'Steady State Relaxation Before', Stat )      
        IF (.NOT. Stat ) RelaxBefore = .TRUE.
      END IF

      ! Steady state system has never any constraints
      SkipConstraints = .FALSE.
      
    ELSE
      iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      IterNo = NINT( iterVar % Values(1) )
      Solver % Variable % NonlinIter = IterNo

      Skip = ListGetLogical( SolverParams,'Skip Advance Nonlinear iter',Stat)
      IF( .NOT. Skip )  iterVar % Values(1) = IterNo + 1 

      IF( .NOT. Solver % NewtonActive ) THEN
        i = ListGetInteger( SolverParams, 'Nonlinear System Newton After Iterations',Stat )
        IF( Stat .AND. i <= IterNo ) Solver % NewtonActive = .TRUE.
      END IF     
      
      Skip = ListGetLogical( SolverParams,'Skip Compute Nonlinear Change',Stat)

      IF(Skip) THEN
        CALL Info(Caller,'Skipping the computation of nonlinear change',Level=15)
        RETURN
      END IF
        
      ResidualMode = ListGetLogical( SolverParams,'Linear System Residual Mode',Stat)

      ConvergenceType = ListGetString(SolverParams,&
          'Nonlinear System Convergence Measure',Stat)
      IF(.NOT. stat) ConvergenceType = 'norm' 

      ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Nonlinear System Convergence Absolute',Stat)
      IF(.NOT. Stat) ConvergenceAbsolute = &
          ListGetLogical(SolverParams,'Use Absolute Norm for Convergence',Stat)
              
      Relaxation = ListGetCReal( SolverParams, &
          'Nonlinear System Relaxation Factor', Relax )
      Relax = Relax .AND. ( ABS( Relaxation - 1.0_dp) > EPSILON( Relaxation ) )
      IF( Relax ) THEN
        RelaxAfter = ListGetInteger(SolverParams,'Nonlinear System Relaxation After',Stat)
        IF( Stat .AND. RelaxAfter >= Solver % Variable % NonlinIter ) Relax = .FALSE.

        RelativeP = ListGetLogical( SolverParams,'Relative Pressure Relaxation',Stat) 
        IF( RelativeP) CALL Info(Caller,'Using relative pressure relaxation',Level=10)
      END IF
      
      SkipConstraints = ListGetLogical(SolverParams,&
          'Nonlinear System Convergence Without Constraints',Stat) 

      RelaxBefore = .TRUE.
      IF(Relax) THEN
        RelaxBefore = ListGetLogical( SolverParams, &
            'Nonlinear System Relaxation Before', Stat )
        IF (.NOT. Stat ) RelaxBefore = .TRUE.
      END IF
    END IF

    
    IF(PRESENT(values)) THEN
      x => values
    ELSE
      x => Solver % Variable % Values      
    END IF

    IF ( .NOT. ASSOCIATED(x) ) THEN
      Solver % Variable % Norm = 0.0d0 
      IF(SteadyState) THEN
        Solver % Variable % SteadyChange = 0.0d0
      ELSE
        Solver % Variable % NonlinChange = 0.0d0
      END IF
      RETURN
    END IF

    
    IF(PRESENT(nsize)) THEN
      n = nsize 
    ELSE 
      n = SIZE( x )
    END IF

    IF( SkipConstraints ) n = MIN( n, Solver % Matrix % NumberOfRows )

    Stat = .FALSE.
    x0 => NULL()
    IF(PRESENT(values0)) THEN
      x0 => values0
      Stat = .TRUE.
    ELSE IF(SteadyState) THEN
      IF( ASSOCIATED(Solver % Variable % SteadyValues) ) THEN
        x0 => Solver % Variable % SteadyValues
        Stat = .TRUE.
      END IF
    ELSE 
      IF( ASSOCIATED(Solver % Variable % NonlinValues)) THEN
        x0 => Solver % Variable % NonlinValues
        Stat = .TRUE.
      ELSE
        x0 => Solver % Variable % Values
        Stat = .TRUE.
      END IF
    END IF
    
    IF(Stat .AND. .NOT. SkipConstraints ) THEN
      IF (SIZE(x0) /= SIZE(x)) CALL Info(Caller,'WARNING: Possible mismatch in length of vectors!',Level=10)
    END IF

    ! This ensures that the relaxation does not affect the mean of the pressure
    IF( RelativeP ) THEN
      dofs = Solver % Variable % Dofs

      dpsum = SUM(x(dofs:n:dofs)) - SUM(x0(dofs:n:dofs)) 
      nsum = 1.0_dp * n / dofs

      dpsum = ParallelReduction( dpsum ) 
      nsum = ParallelReduction( nsum )

      Poffset = (1-Relaxation) * dpsum / nsum
    END IF

    
    IF( ResidualMode ) THEN
      IF(Relax .AND. RelaxBefore) THEN
        x(1:n) = x0(1:n) + Relaxation*x(1:n)
      ELSE
        x(1:n) = x0(1:n) + x(1:n)
      END IF
    ELSE 
      IF(Relax .AND. RelaxBefore) THEN
        x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
        IF( RelativeP ) x(dofs:n:dofs) = x(dofs:n:dofs) + Poffset
      END IF
    END IF

    IF(SteadyState) THEN
      PrevNorm = Solver % Variable % PrevNorm
    ELSE
      PrevNorm = Solver % Variable % Norm
    END IF

    Norm = ComputeNorm(Solver, n, x)
    Solver % Variable % Norm = Norm
    
    !--------------------------------------------------------------------------
    ! The norm should be bounded in order to reach convergence
    !--------------------------------------------------------------------------
    IF( Norm /= Norm ) THEN
      CALL NumericalError(Caller,'Norm of solution appears to be NaN')
    END IF

    IF( SteadyState ) THEN
      MaxNorm = ListGetCReal( SolverParams, &
          'Steady State Max Norm', Stat )
    ELSE
      MaxNorm = ListGetCReal( SolverParams, &
          'Nonlinear System Max Norm', Stat )
    END IF    

    IF( Stat ) THEN
      CALL Info(Caller,Message)
      CALL NumericalError(Caller,'Norm of solution exceeded given bounds')
    END IF
      
    SELECT CASE( ConvergenceType )
        
    CASE('residual')
      !--------------------------------------------------------------------------
      ! x is solution of A(x0)x=b(x0) thus residual should really be r=b(x)-A(x)x
      ! Instead we use r=b(x0)-A(x0)x0 which unfortunately is one step behind.
      !--------------------------------------------------------------------------
      IF(PRESENT(Matrix)) THEN
        A => Matrix
      ELSE
        A => Solver % Matrix
      END IF

      IF(PRESENT(RHS)) THEN
        b => RHS
      ELSE
        b => Solver % Matrix % rhs
      END IF
      
      ALLOCATE(r(n))
      r=0._dp

      IF (Parenv % Pes>1) THEN
        ALLOCATE( TmpRHSVec(n), TmpXVec(n) )

        nn = A % ParMatrix % SplittedMatrix % InsideMatrix % NumberOfRows

        TmpRhsVec = b
        CALL ParallelInitSolve( A, tmpXVec, TmpRhsVec, r)

        tmpXvec = x0(1:n)
        CALL ParallelVector(a,TmpXvec)
        CALL ParallelVector(A,tmpRhsvec)

        CALL ParallelMatrixVector(A, TmpXvec, r)
        DO i=1,nn
          r(i) = r(i) - tmprhsvec(i)
        END DO

        Change = ParallelNorm(nn,r)
        bNorm =  ParallelNorm(nn,tmpRhsVec)
      ELSE
        CALL MatrixVectorMultiply( A, x0, r)
        DO i=1,n
          r(i) = r(i) - b(i)
        END DO
        Change = ComputeNorm(Solver, n, r)
        bNorm  = ComputeNorm(Solver, n, b)
      END IF


      IF(.NOT. ConvergenceAbsolute) THEN
        IF(bNorm > 0.0) THEN
          Change = Change / bNorm
        END IF
      END IF
      DEALLOCATE(r)
      
    CASE('linear system residual')
      !--------------------------------------------------------------------------
      ! Here the true linear system residual r=b(x0)-A(x0)x is computed.
      ! This option is useful for certain special solvers.  
      !--------------------------------------------------------------------------
      A => Solver % Matrix
      b => Solver % Matrix % rhs
      
      IF (ParEnv % Pes > 1) THEN

        ALLOCATE( TmpRHSVec(n), TmpXVec(n), TmpRVec(n) )
        TmpRHSVec(1:n) = b(1:n)
        TmpXVec(1:n) = x(1:n)
        TmpRVec(1:n) = 0.0d0

        CALL ParallelVector(A, TmpRHSVec)
        CALL ParallelVector(A, TmpXVec)       
        CALL SParMatrixVector( TmpXVec, TmpRVec, ipar )
 
        nn = A % ParMatrix % SplittedMatrix % InsideMatrix % NumberOfRows

        DO i=1, nn
          TmpRVec(i) = TmpRHSVec(i) - TmpRVec(i)
        END DO

        Change = ParallelNorm( nn, TmpRVec )

        IF(.NOT. ConvergenceAbsolute) THEN
          bNorm = ParallelNorm( nn, TmpRHSVec )
          IF(bNorm > 0.0) THEN
            Change = Change / bNorm
          END IF
        END IF
        DEALLOCATE( TmpRHSVec, TmpXVec, TmpRVec )
      ELSE	
        ALLOCATE(r(n)) 
        CALL MatrixVectorMultiply( A, x, r)
        DO i=1,n
          r(i) = r(i) - b(i)
        END DO
        Change = SQRT( DOT_PRODUCT( r(1:n), r(1:n) ) )
        IF(.NOT. ConvergenceAbsolute) THEN
          bNorm = SQRT( DOT_PRODUCT( b(1:n), b(1:n) ) )
          IF(bNorm > 0.0) THEN
            Change = Change / bNorm
          END IF
        END IF
        DEALLOCATE(r)	
      END IF
      
    CASE('solution')      

      ALLOCATE(r(n))
      r = x(1:n)-x0(1:n)
      Change = ComputeNorm(Solver, n, r)
      IF( .NOT. ConvergenceAbsolute ) THEN
        IF( Norm + PrevNorm > 0.0) THEN
          Change = Change * 2.0_dp/ (Norm+PrevNorm)
        END IF
      END IF
      DEALLOCATE(r)      

    CASE('norm')

      Change = ABS( Norm-PrevNorm )
      IF( .NOT. ConvergenceAbsolute .AND. Norm + PrevNorm > 0.0) THEN
        Change = Change * 2.0_dp/ (Norm+PrevNorm)
      END IF
      
    CASE DEFAULT
      CALL Warn(Caller,'Unknown convergence measure: '//TRIM(ConvergenceType))    
      
    END SELECT
    
    !--------------------------------------------------------------------------
    ! Check for convergence: 0/1
    !--------------------------------------------------------------------------
    IF(SteadyState) THEN
      PrevChange = Solver % Variable % SteadyChange
      Solver % Variable % SteadyChange = Change
      Tolerance = ListGetCReal( SolverParams,'Steady State Convergence Tolerance',Stat)
      IF( Stat ) THEN
        IF( Change <= Tolerance ) THEN
          Solver % Variable % SteadyConverged = 1
        ELSE
          Solver % Variable % SteadyConverged = 0
        END IF          
      END IF
      
      Tolerance = ListGetCReal( SolverParams,'Steady State Divergence Limit',Stat)
      IF( Stat .AND. Change > Tolerance ) THEN
        IF( IterNo > 1 .AND. Change > PrevChange ) THEN
          CALL Info(Caller,'Steady state iteration diverged over tolerance')
          Solver % Variable % SteadyConverged = 2
        END IF
      END IF
      
      Tolerance = ListGetCReal( SolverParams,'Steady State Exit Condition',Stat)
      IF( Stat .AND. Tolerance > 0.0 ) THEN
        CALL Info(Caller,'Nonlinear iteration condition enforced by exit condition',Level=6)
        Solver % Variable % SteadyConverged = 3
      END IF

    ELSE
      PrevChange = Solver % Variable % NonlinChange 
      Solver % Variable % NonlinChange = Change
      Solver % Variable % NonlinConverged = 0

      MaxIter = ListGetInteger( SolverParams,'Nonlinear System Max Iterations',Stat)            
      
      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Convergence Tolerance',Stat)
      IF( Stat ) THEN
        IF( Change <= Tolerance ) THEN
          Solver % Variable % NonlinConverged = 1
        ELSE IF( IterNo >= MaxIter ) THEN
          IF( ListGetLogical( SolverParams,'Nonlinear System Abort Not Converged',Stat ) ) THEN
            CALL Fatal(Caller,'Nonlinear iteration did not converge to tolerance')
          ELSE
            CALL Info(Caller,'Nonlinear iteration did not converge to tolerance',Level=6)
            ! Solver % Variable % NonlinConverged = 2            
          END IF
        END IF
      END IF

      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Divergence Limit',Stat)
      IF( Stat .AND. Change > Tolerance ) THEN        
        IF( ( IterNo > 1 .AND. Change > PrevChange ) .OR. ( IterNo >= MaxIter ) ) THEN
          IF( ListGetLogical( SolverParams,'Nonlinear System Abort Diverged',Stat ) ) THEN
            CALL Fatal(Caller,'Nonlinear iteration diverged over limit')
          ELSE
            CALL Info(Caller,'Nonlinear iteration diverged over limit',Level=6)
            Solver % Variable % NonlinConverged = 2
          END IF
        END IF
      END IF

      Tolerance = ListGetCReal( SolverParams,'Nonlinear System Exit Condition',Stat)
      IF( Stat .AND. Tolerance > 0.0 ) THEN
        CALL Info(Caller,'Nonlinear iteration condition enforced by exit condition',Level=6)
        Solver % Variable % NonlinConverged = 3
      END IF
      
      IF( Solver % Variable % NonlinConverged > 1 ) THEN
        MinIter = ListGetInteger( SolverParams,'Nonlinear System Min Iterations',Stat)
        IF( Stat .AND. IterNo < MinIter ) THEN
          CALL Info(Caller,'Enforcing continuation of iteration')
          Solver % Variable % NonlinConverged = 0
        END IF
      END IF
      
      IF( .NOT. Solver % NewtonActive ) THEN
        Tolerance = ListGetCReal( SolverParams, 'Nonlinear System Newton After Tolerance',Stat )
        IF( Stat .AND. Change < Tolerance ) Solver % NewtonActive = .TRUE.
      END IF     
    END IF

    
    IF(Relax .AND. .NOT. RelaxBefore) THEN
      x(1:n) = (1-Relaxation)*x0(1:n) + Relaxation*x(1:n)
      IF( RelativeP ) x(dofs:n:dofs) = x(dofs:n:dofs) + Poffset
      Solver % Variable % Norm = ComputeNorm(Solver,n,x)
    END IF
    
    ! Steady state output is done in MainUtils
    SolverName = ListGetString( SolverParams, 'Equation',Stat)
    IF(.NOT. Stat) SolverName = Solver % Variable % Name
 
    IF(SteadyState) THEN        
      WRITE( Message, '(a,g15.8,g15.8,a)') &
         'SS (ITER='//TRIM(i2s(IterNo))//') (NRM,RELC): (',Norm, Change,&
          ' ) :: '// TRIM(SolverName)
    ELSE
      WRITE( Message, '(a,g15.8,g15.8,a)') &
         'NS (ITER='//TRIM(i2s(IterNo))//') (NRM,RELC): (',Norm, Change,&
          ' ) :: '// TRIM(SolverName)
    END IF
    CALL Info( Caller, Message, Level=3 )

    ! This provides a way to directly save the convergence data into an external
    ! file making it easier to follow the progress of Elmer simulation in other software.
    !------------------------------------------------------------------------------------    
    IF( ListGetLogical( CurrentModel % Simulation,'Convergence Monitor',Stat ) ) THEN
      CALL WriteConvergenceInfo()  
    END IF
    
    ! Optional a posteriori scaling for the computed fields
    ! May be useful for some floating systems where one want to impose some intergral 
    ! constraints without actually using them. Then first use just one Dirichlet point
    ! and then fix the level a posteriori using this condition. 
    !----------------------------------------------------------------------------------
    DoIt = .FALSE.
    IF( SteadyState ) THEN 
      DoIt = ListGetLogical( SolverParams,&
          'Nonlinear System Set Average Solution',Stat)
    ELSE 
      DoIt = ListGetLogical( SolverParams,&
          'Linear System Set Average Solution',Stat)
    END IF
    IF( DoIt ) THEN
      IF( ParEnv % PEs > 1 ) THEN
        CALL Fatal(Caller,'Setting average value not implemented in parallel!')
      END IF
      Ctarget = ListGetCReal( SolverParams,'Average Solution Value',Stat)      
      str = ListGetString( SolverParams,'Average Solution Weight Variable',Stat)
      IF( Stat ) THEN
        WeightVar => VariableGet( Solver % Mesh % Variables, str )
        IF( .NOT. ASSOCIATED( WeightVar ) ) THEN
          CALL Fatal(Caller,'> Average Solution Weight < missing: '//TRIM(str))
        END IF
        IF( SIZE(x) /= SIZE(WeightVar % Values ) ) THEN
          CALL Fatal(Caller,'Field and weight size mismatch: '//TRIM(str))          
        END IF
        Ctarget = Ctarget - SUM( WeightVar % Values * x ) / SUM( WeightVar % Values )
      ELSE
        Ctarget = Ctarget - SUM(x) / SIZE(x)
      END IF
      x = x + Ctarget
    END IF


    ! Calculate derivative a.k.a. sensitivity
    DoIt = .FALSE.
    IF( SteadyState ) THEN
      DoIt = ListGetLogical( SolverParams,'Calculate Derivative',Stat )
    ELSE
      DoIt = ListGetLogical( SolverParams,'Nonlinear Calculate Derivative',Stat )
    END IF

    IF( DoIt ) THEN
      IF( SteadyState ) THEN
        iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
        IterNo = NINT( iterVar % Values(1) )
      ELSE
        iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
        IterNo = NINT( iterVar % Values(1) )
      END IF
      
      Eps = 1.0_dp
      IF( IterNo > 1 ) THEN
        dtVar => VariableGet( Solver % Mesh % Variables, 'derivative eps' )
        IF( ASSOCIATED( dtVar ) ) THEN
          eps = dtVar % Values(1)
          Stat = .TRUE.
        ELSE
          eps = ListGetCReal( SolverParams,'derivative eps',Stat)
        END IF
        IF(.NOT. Stat) THEN
          Eps = 1.0_dp
          CALL Info(Caller,'Derivative Eps not given, using one',Level=7)
        END IF
      END IF
      
      str = GetVarname(Solver % Variable) // ' Derivative'
      VeloVar => VariableGet( Solver % Mesh % Variables, str )
      IF( ASSOCIATED( VeloVar ) ) THEN
        CALL Info(Caller,'Computing variable:'//TRIM(str))
        VeloVar % Values = (x - x0) / eps
      ELSE
        CALL Warn(Caller,'Derivative variable not present')
      END IF

    END IF
    
    IF(.NOT. SteadyState ) THEN    
      CALL UpdateDependentObjects( Solver, .FALSE. )        
    END IF


  CONTAINS

    SUBROUTINE WriteConvergenceInfo()

      INTEGER :: ConvInds(5),ConvUnit
      CHARACTER(LEN=MAX_NAME_LEN) :: ConvFile
      LOGICAL, SAVE :: ConvVisited = .FALSE.

      IF( ParEnv % MyPe /= 0 ) RETURN

      ConvFile = ListGetString(CurrentModel % Simulation,&
          'Convergence Monitor File',Stat)
      IF(.NOT. Stat) ConvFile = 'convergence.dat'

      IF( ConvVisited ) THEN
        OPEN(NEWUNIT=ConvUnit, FILE=ConvFile,STATUS='old',POSITION='append')
      ELSE
        OPEN(NEWUNIT=ConvUnit, File=ConvFile)
        WRITE(ConvUnit,'(A)') '! solver  ss/ns  timestep  coupled  nonlin  norm  change'
        ConvVisited = .TRUE.
      END IF

      ConvInds = 0
      ConvInds(1) = Solver % SolverId

      IF( SteadyState ) ConvInds(2) = 1 

      iterVar => VariableGet( Solver % Mesh % Variables, 'timestep' )
      ConvInds(3) = NINT( iterVar % Values(1) )

      iterVar => VariableGet( Solver % Mesh % Variables, 'coupled iter' )
      ConvInds(4) = NINT( iterVar % Values(1) )

      iterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      ConvInds(5) = NINT( iterVar % Values(1) )

      WRITE(ConvUnit,'(5I8,2G16.8)') ConvInds,Norm,Change
      CLOSE(ConvUnit)
      
    END SUBROUTINE WriteConvergenceInfo
        
!------------------------------------------------------------------------------
  END SUBROUTINE ComputeChange
!------------------------------------------------------------------------------
    



!------------------------------------------------------------------------------
!> Adaptive version for getting gaussian integration points.
!> Also saves some time in initializations.
!> Note: the routine uses the pointer to Solver to check whether definitions
!> need to be remade. 
!----------------------------------------------------------------------------------------------

  FUNCTION GaussPointsAdapt( Element, Solver, PReferenceElement ) RESULT(IntegStuff)

    IMPLICIT NONE
    TYPE(Element_t) :: Element
    TYPE(Solver_t), OPTIONAL, TARGET :: Solver
    LOGICAL, OPTIONAL :: PReferenceElement           !< For switching to the p-version reference element
    TYPE( GaussIntegrationPoints_t ) :: IntegStuff   !< Structure holding the integration points

    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, GaussDef
    TYPE(Solver_t), POINTER :: pSolver, prevSolver => NULL()
    TYPE(Variable_t), POINTER :: IntegVar
    INTEGER :: AdaptOrder, AdaptNp, Np, RelOrder
    REAL(KIND=dp) :: MinLim, MaxLim, MinV, MaxV, V
    LOGICAL :: UseAdapt, Found,ElementalRule
    INTEGER :: i,n,ElementalNp(8),prevVisited = -1
    LOGICAL :: Debug, InitDone 
    
    SAVE prevSolver, UseAdapt, MinLim, MaxLim, IntegVar, AdaptOrder, AdaptNp, RelOrder, Np, &
        ElementalRule, ElementalNp, prevVisited

    IF( PRESENT( Solver ) ) THEN
      pSolver => Solver
    ELSE
      pSolver => CurrentModel % Solver
    END IF
       
    !Debug = ( Element % ElementIndex == 1)

    InitDone = ASSOCIATED( pSolver, prevSolver ) .AND. ( prevVisited == pSolver % TimesVisited )
        
    IF( .NOT. InitDone ) THEN
      RelOrder = ListGetInteger( pSolver % Values,'Relative Integration Order',Found )
      AdaptNp = 0
      Np = ListGetInteger( pSolver % Values,'Number of Integration Points',Found )
      
      GaussDef = ListGetString( pSolver % Values,'Element Integration Points',ElementalRule )
      IF( ElementalRule ) THEN
        CALL ElementalGaussRules( GaussDef )
      END IF
      
      VarName = ListGetString( pSolver % Values,'Adaptive Integration Variable',UseAdapt )
      IF( UseAdapt ) THEN
        CALL Info('GaussPointsAdapt','Using adaptive gaussian integration rules')
        IntegVar => VariableGet( pSolver % Mesh % Variables, VarName )
        IF( .NOT. ASSOCIATED( IntegVar ) ) THEN
          CALL Fatal('GaussPointsAdapt','> Adaptive Integration Variable < does not exist')
        END IF
        IF( IntegVar % TYPE /= Variable_on_nodes ) THEN
          CALL Fatal('GaussPointsAdapt','Wrong type of integration variable!')
        END IF
        MinLim = ListGetCReal( pSolver % Values,'Adaptive Integration Lower Limit' )
        MaxLim = ListGetCReal( pSolver % Values,'Adaptive Integration Upper Limit' )
        AdaptNp = ListGetInteger( pSolver % Values,'Adaptive Integration Points',Found )
        IF(.NOT. Found ) THEN
          AdaptOrder = ListGetInteger( pSolver % Values,'Adaptive Integration Order',Found )        
        END IF
        IF(.NOT. Found ) AdaptOrder = 1
        !PRINT *,'Adaptive Integration Strategy:',MinV,MaxV,AdaptOrder,AdaptNp
      END IF

      prevSolver => pSolver
      prevVisited = pSolver % TimesVisited
    END IF

    IF( UseAdapt ) THEN
      RelOrder = 0
      Np = 0

      n = Element % TYPE % NumberOfNodes        
      MinV = MINVAL( IntegVar % Values( IntegVar % Perm( Element % NodeIndexes(1:n) ) ) )
      MaxV = MAXVAL( IntegVar % Values( IntegVar % Perm( Element % NodeIndexes(1:n) ) ) )
      
      IF( .NOT. ( MaxV < MinLim .OR. MinV > MaxLim ) ) THEN
        RelOrder = AdaptOrder
        Np = AdaptNp
      END IF
    END IF
      
    !IF( Debug ) PRINT *,'Adapt',UseAdapt,Element % ElementIndex, n,MaxV,MinV,MaxLim,MinLim,Np,RelOrder

    IF( ElementalRule ) THEN
      Np = ElementalNp( Element % TYPE % ElementCode / 100 ) 
    END IF
    
    IF( Np > 0 ) THEN
      IntegStuff = GaussPoints( Element, Np = Np, PReferenceElement = PReferenceElement ) 
    ELSE IF( RelOrder /= 0 ) THEN
      IntegStuff = GaussPoints( Element, RelOrder = RelOrder, PReferenceElement = PReferenceElement ) 
    ELSE      
      IntegStuff = GaussPoints( Element, PReferenceElement = PReferenceElement ) 
    END IF

    !IF( Debug ) PRINT *,'Adapt real nodes',IntegStuff % n
    

  CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE ElementalGaussRules(GaussDef)
!------------------------------------------------------------------------------
      CHARACTER(LEN=*) :: GaussDef
!------------------------------------------------------------------------------
      INTEGER  :: i,j,k,n,m
      

      n = LEN_TRIM(GaussDef)
      ElementalNp = 0

      ! PRINT *,'gauss def:',GaussDef(1:n)
      
      DO i=2,8
        j = 0
        SELECT CASE( i )
        CASE( 2 )
          j =  INDEX( GaussDef(1:n), '-line' ) ! position of string "-line"
          m = 5 ! length of string "-line"
        CASE( 3 )
          j =  INDEX( GaussDef(1:n), '-tri' ) 
          m = 4
        CASE( 4 )
          j =  INDEX( GaussDef(1:n), '-quad' )
          m = 5
        CASE( 5 )
          j =  INDEX( GaussDef(1:n), '-tetra' )
          m = 6
        CASE( 6 )
          j =  INDEX( GaussDef(1:n), '-pyramid' )
          m = 8
        CASE( 7 )
          j =  INDEX( GaussDef(1:n), '-prism' )
          m = 6
        CASE( 8 )
          j =  INDEX( GaussDef(1:n), '-brick' )
          m = 6
        END SELECT

        IF( j > 0 ) THEN
          READ( GaussDef(j+m:n), * ) k
          ElementalNp(i) = k
        END IF
      END DO

      ! PRINT *,'Elemental Gauss Rules:',ElementalNp
      
!------------------------------------------------------------------------------
    END SUBROUTINE ElementalGaussRules
!------------------------------------------------------------------------------
    
    
  END FUNCTION GaussPointsAdapt
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
!> Checks stepsize of a linear system so that the error has decreased.
!> Various indicatators and search algorithms have been implemented,
!------------------------------------------------------------------------------
  FUNCTION CheckStepSize(Solver,FirstIter,&
      nsize,values,values0) RESULT ( ReduceStep ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t) :: Solver
    LOGICAL :: FirstIter
    INTEGER, OPTIONAL :: nsize
    REAL(KIND=dp), OPTIONAL, TARGET :: values(:), values0(:)
    LOGICAL :: ReduceStep
!------------------------------------------------------------------------------
    INTEGER :: MaxTests=0,tests,MaxNonlinIter,NonlinIter, Dofs
    REAL(KIND=dp) :: Residual0, Residual1, Residual
    INTEGER :: i,n,m,ForceDof, SearchMode, CostMode, iter = 0
    TYPE(Matrix_t), POINTER :: A, MP
    TYPE(Variable_t), POINTER :: IterVar, Var
    REAL(KIND=dp), POINTER :: b(:), x(:), x0(:), r(:), x1(:), x2(:), mr(:), mx(:), mb(:)
    REAL(KIND=dp) :: Norm, PrevNorm, rNorm, bNorm, Relaxation, Alpha, Myy, &
        NonlinTol, LineTol, Cost0(4), Cost1(4), Cost(4), OrthoCoeff, x0norm, x1norm, Change, &
        LinTol
    REAL(KIND=dp), ALLOCATABLE :: TempRHS(:)
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: Stat, Init, Newton, Ortho, Debug, SaveToFile
    CHARACTER(LEN=MAX_NAME_LEN) :: SolverName, ConvergenceType, FileName
    TYPE(ValueList_t), POINTER :: SolverParams


    SAVE SolverParams, Alpha, Myy, Relaxation, MaxTests, tests, &
        Residual, NonlinTol, LinTol, x1, x0, LineTol, CostMode, SearchMode, &
        Cost0, Residual0, Cost1, n, Dofs, ForceDof, Ortho, Newton, &
        ConvergenceType, Norm, PrevNorm, iter, FileName, SaveToFile

    Debug = .FALSE.
    
    SolverParams => Solver % Values
    Var => Solver % Variable
    Dofs = Var % Dofs
 
    IF(PRESENT(values)) THEN
      x => values
    ELSE 
      x => Var % Values      
    END IF


    ! Assembly the vectors, if needed, and 
    ! also at first time get the line search parameters.
    !----------------------------------------------------
    IF( FirstIter ) THEN
      CALL Info('CheckStepSize','Initializing step-size search',Level=6)

      IF(PRESENT(nsize)) THEN
        n = nsize
      ELSE
        n = SIZE(x)
      END IF
      
      IF( ASSOCIATED( x0 ) ) THEN
        IF( SIZE(x0) /= n ) DEALLOCATE( x0 )
      END IF
      
      IF( PRESENT( values0 ) ) THEN
        x0 => values0 
      ELSE
        IF( .NOT. ASSOCIATED( x0 ) ) THEN
          ALLOCATE( x0(n) )
        END IF
      END IF
      
      IF( ASSOCIATED( x1 ) ) THEN
        IF( SIZE(x1) /= n ) DEALLOCATE( x1 )
      END IF
      IF( .NOT. ASSOCIATED( x1 ) ) THEN
        ALLOCATE( x1(n) )
      END IF

      Norm = 0.0_dp
      Var % NonlinConverged = 0
      Var % NonlinChange = 1.0_dp
      
      ! 1 - Residual norm : |Ax-b| 
      ! 2 - Quadratic functional : x^T(Ax-2b)/2
      ! 3 - Weighted residual : x^T(Ax-b)
      ! 4 - Lumped force : SUM(r_i)
      !------------------------------------------------------------
      CostMode = ListGetInteger( SolverParams,'Nonlinear System Linesearch Cost Mode',Stat)
      IF(.NOT. Stat) CostMode = 1
      
      ! 1 - Armijo-Goldstein criterion & successive relaxation 
      ! 2 - Minimize cost by bisection
      ! 3 - Find the zero cost by bisection
      !------------------------------------------------------------
      SearchMode = ListGetInteger( SolverParams,'Nonlinear System Linesearch Search Mode',Stat)
      IF(.NOT. Stat) SearchMode = 1

      ! Should the search direction be orthogonalized 
      !-----------------------------------------------------------
      Ortho = ListGetLogical( SolverParams,'Nonlinear System Linesearch Orthogonal',Stat)

      ! Is the outer ieration performed by Newton i.e. the search 
      ! should always be differential. 
      !-----------------------------------------------------------
      Newton = ListGetLogical( SolverParams,'Nonlinear System Linesearch Newton',Stat)

      NonlinTol = ListGetConstReal( SolverParams, &
          'Nonlinear System Convergence Tolerance', Stat )
      LinTol = ListGetConstReal( SolverParams, &
          'Linear System Convergence Tolerance', Stat )

      MaxNonlinIter = ListGetInteger( SolverParams,&
            'Nonlinear System Max Iterations',Stat)
      IF( MaxNonlinIter <= 2 ) THEN
        CALL Warn('CheckStepSize','For linesearch to work the nonlin iterations should be larger: '&
            //I2S(MaxNonlinIter))
      END IF

      ConvergenceType = ListGetString(SolverParams,&
          'Nonlinear System Convergence Measure',Stat)
      IF(.NOT. Stat) ConvergenceType = 'norm'

      ! Parameters related to line search algorithms
      !------------------------------------------------
      MaxTests = ListGetInteger( SolverParams,&
          'Nonlinear System Linesearch Iterations',Stat)
      IF( .NOT. Stat ) MaxTests = 10

      Myy = ListGetConstReal( SolverParams, &
          'Nonlinear System Linesearch Limit', Stat )
      IF(.NOT. Stat) Myy = 0.5_dp

      Relaxation = ListGetConstReal( SolverParams, &
          'Nonlinear System Linesearch Factor', Stat )
      IF(.NOT. Stat) Relaxation = 0.5_dp

      LineTol = ListGetConstReal( SolverParams, &
          'Nonlinear System Linesearch Tolerance', Stat )

      ForceDof = ListGetInteger( SolverParams, &
          'Nonlinear System Linesearch Force Index', Stat )
      IF(.NOT. Stat) ForceDof = Dofs

      FileName = ListGetString( SolverParams, &
          'Nonlinear System Linesearch Filename', SaveToFile )
      
      ! Computation of nonlinear change is now done with this routine
      ! so skip computing the change in the standard slot.
      !---------------------------------------------------------------
      CALL ListAddLogical(SolverParams,&
          'Skip Compute Nonlinear Change',.TRUE.)
    END IF

    !--------------------------------------------------------------------------
    ! This is the real residual: r=b-Ax
    ! We hope to roughly minimize L2 norm of r, or some related quantity
    !--------------------------------------------------------------------------
    A => Solver % Matrix
    b => Solver % Matrix % rhs

    ALLOCATE(r(n))
    IF (Parenv % Pes>1) THEN
      ALLOCATE(TempRHS(n))
      r = 0._dp
      TempRHS(1:n) = b(1:n)
      CALL ParallelInitSolve( A, x, TempRHS, r )

      MP => ParallelMatrix(A,mx,mb,mr)
      m = MP % NumberOfRows

      CALL ParallelMatrixVector( A, mx, r)
      r(1:m) = r(1:m) - TempRHS(1:m)
      Residual= ParallelNorm(n,r)
    ELSE
      CALL MatrixVectorMultiply( A, x, r)
      r(1:n) = r(1:n) - b(1:n)
      Residual = ComputeNorm(Solver, n, r)
    END IF

    ! Currently we compute all the costs to make it easier to study the 
    ! behavior of different measures when doing linesearch.
    IF( .TRUE. ) THEN
      Cost(1) = Residual
      Cost(2) = SUM( 0.5_dp * x(1:n) * ( r(1:n) - b(1:n) ) )
      Cost(3) = SUM( x(1:n) * r(1:n) )
      Cost(4) = SUM( r(ForceDof::Dofs) )
    ELSE      
      IF( CostMode == 1 ) THEN
        Cost(1) = Residual
      ELSE IF( CostMode == 2 ) THEN
        Cost(2) = SUM( 0.5_dp * x(1:n) * ( r(1:n) - b(1:n) ) )
      ELSE IF( CostMode == 3 ) THEN
        Cost(3) = SUM( x(1:n) * r(1:n) )
      ELSE IF( CostMode == 4 ) THEN
        Cost(4) = SUM( r(ForceDof::Dofs) )
      ELSE
        CALL Fatal('CheckStepSize','Unknown CostMode: '//TRIM(I2S(SearchMode)))
      END IF
      DEALLOCATE(r)
    END IF

    WRITE( Message,'(A,4ES15.7)') 'Cost: ',Cost
    CALL Info('CheckStepSize',Message,Level=8)

    ! At first iteration we cannot really do anything but go further 
    ! and save the reference residual for comparison.
    !-----------------------------------------------------------------------------
    IF( FirstIter ) THEN
      Tests = 0
      ReduceStep = .FALSE.
      x0(1:n) = x(1:n)
      Cost0 = Cost
      Residual0 = Residual

      IF( Debug ) THEN
        PRINT *,'x0 range: ',MINVAL(x0),MAXVAL(x0)
        PRINT *,'b0 range: ',MINVAL(b),MAXVAL(b)
        PRINT *,'Cost0: ',Cost0
      END IF

      IF( SaveToFile ) THEN
        CALL Info('CheckStepSize','Saving step information into file: '&
            //TRIM(FileName),Level=10)
        OPEN( 10, FILE = FileName, STATUS='UNKNOWN' )
        i = 0
        WRITE (10,'(2I6,5ES15.7)') Tests,i,Alpha,Cost
        CLOSE( 10 )
      END IF


      RETURN
    END IF

    Tests = Tests + 1

    IF( Tests == 1 ) THEN
      ! Started with no relaxation
      !---------------------------
      x1 = x
      Alpha = 1.0_dp
      Cost1 = Cost

      ! This is just debugging code waiting to be reused
      IF( .FALSE. ) THEN
        iter = iter + 1

        PRINT *,'Iter: ',iter
        NULLIFY( x2 ) 
        ALLOCATE( x2(n/2) ) 
        x2 = x(1::2)
        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
            'xiter '//TRIM(I2S(iter)),1,x2,Solver % Variable % Perm ) 
        PRINT *,'Xiter range:',MINVAL(x2),MAXVAL(x2)
        NULLIFY(x2)
        
!        NULLIFY( x2 ) 
!        ALLOCATE( x2(n/2) ) 
!        x2 = x(2::2)
!        CALL VariableAdd( Solver % Mesh % Variables, Solver % Mesh, Solver, &
!            'yiter '//TRIM(I2S(iter)),1,x2,Solver % Variable % Perm ) 
!        NULLIFY(x2)
      END IF
      
      IF( Debug ) THEN
        PRINT *,'b1 range: ',MINVAL(b),MAXVAL(b)
        PRINT *,'x1 range: ',MINVAL(x1),MAXVAL(x1)
        PRINT *,'Cost1: ',Cost1
      END IF

      ! Orthonormalization:
      ! The direction 'x0' has already been exhausted so remove that from 'x1'
      !-----------------------------------------------------------------------
      x0norm = ComputeNorm( Solver, n, x0 )
      IF( Ortho ) THEN
        IF( x0norm > EPSILON( x0norm ) ) THEN
          OrthoCoeff = SUM(x1*x0) / ( x0norm**2 )
          x1 = x1 - OrthoCoeff * x0
        END IF
      ELSE
        ! This basically checks whether the new and old solution is so 
        ! close that there is no point of finding better solution.
        x1 = x1 - x0 
        x1norm = ComputeNorm(Solver, n, x1)
        IF( x1norm < LinTol * x0norm ) THEN
          ReduceStep = .FALSE.
          GOTO 100
        END IF
      END IF

      IF( Debug ) THEN
        PRINT *,'x1 range orto: ',MINVAL(x1),MAXVAL(x1)
      END IF
    END IF

    ! Armijo-GoldStein Criterion for accepting stepsize
    !-----------------------------------------------------------------
    IF( SearchMode == 1 ) THEN
      ReduceStep = ArmijoGoldsteinSearch(Tests, Alpha )
    ELSE IF( SearchMode == 2 ) THEN
      ReduceStep = BisectMinimumSearch(Tests, Alpha) 
    ELSE IF( SearchMode == 3 ) THEN
      ReduceStep = BisectZeroSearch(Tests, Alpha)       
    ELSE
      CALL Fatal('CheckStepSize','Unknown SearchMode: '//TRIM(I2S(SearchMode)))
    END IF


    IF( SaveToFile ) THEN
      CALL Info('CheckStepSize','Saving step information into file: '&
          //TRIM(FileName),Level=10)
      OPEN( 10, FILE = FileName, POSITION='APPEND',STATUS='OLD' )
      IF( ReduceStep ) THEN
        i = 0
      ELSE
        i = 1
      END IF

      WRITE (10,'(2I6,5ES13.6)') Tests,i,Alpha,Cost
      CLOSE( 10 )
    END IF



100 IF( ReduceStep ) THEN
      IF( Tests >= MaxTests .AND. ReduceStep ) THEN
        CALL Fatal('CheckStepSize','Maximum number of linesearch steps taken without success!')
        ReduceStep = .FALSE.
      END IF
      
      ! New candidate 
      x(1:n) = x0(1:n) + Alpha * x1(1:n)

      WRITE(Message,'(A,I0,A,g15.6)') 'Step ',Tests,' rejected, trying new extent: ',Alpha
      CALL Info( 'CheckStepSize',Message,Level=6 )
    ELSE ! accept step
      WRITE(Message,'(A,I0,A,g15.6)') 'Step ',Tests,' accepted with extent: ',Alpha
      CALL Info( 'CheckStepSize',Message,Level=6 )
      
      ! Chosen candidate
      x(1:n) = x0(1:n) + Alpha * x1(1:n)

      PrevNorm = Norm
      Norm = ComputeNorm(Solver, n, x)

      IF( ConvergenceType == 'residual') THEN
        bNorm = ComputeNorm(Solver, n, b)
        IF( bNorm > 0.0_dp ) Change = Residual / bNorm
      ELSE
        Change = ABS( Norm-PrevNorm )
        IF( Norm + PrevNorm > 0.0) THEN
          Change = Change * 2.0_dp / ( Norm + PrevNorm )
        END IF
      END IF

      Solver % Variable % NonlinChange = Change
      Solver % Variable % Norm = Norm
     
      IF( Solver % Variable % NonlinChange <  NonlinTol ) THEN
        Solver % Variable % NonlinConverged = 1
      END IF
      
      SolverName = ListGetString( SolverParams, 'Equation',Stat)
      IF(.NOT. Stat) SolverName = Solver % Variable % Name
            
      IterVar => VariableGet( Solver % Mesh % Variables, 'nonlin iter')
      m = NINT(IterVar % Values(1))
      
      ! This replaces the standard error output usually written by the ComputeChange
      WRITE( Message, '(a,g15.8,g15.8,a)') &
          'NS (ITER='//TRIM(i2s(m))//') (NRM,RELC): (',Norm, Change, &
          ' ) :: '// TRIM(SolverName)
      CALL Info( 'CheckStepSize', Message, Level=3 )       
      
      WRITE(Message,'(A,I0,A,g15.6)') 'Step accepted after ',tests,' trials: ',Alpha
      CALL Info( 'CheckStepSize',Message,Level=5 )
      WRITE(Message,'(A,g15.6)') 'Previous cost:',Cost0(CostMode)
      CALL Info( 'CheckStepSize',Message,Level=6 )
      WRITE(Message,'(A,g15.6)') 'Initial cost: ',Cost1(CostMode)
      CALL Info( 'CheckStepSize',Message,Level=6 )
      WRITE(Message,'(A,g15.6)') 'Final cost:   ',Cost(CostMode)
      CALL Info( 'CheckStepSize',Message,Level=6 )
      
      Tests = 0
      x0 = x
      
      IF( Debug ) THEN
        PRINT *,'x0 range: ',MINVAL(x0),MAXVAL(x0)
        PRINT *,'Cost0: ',Cost0
        PRINT *,'Residual0: ',Residual0
      END IF

      IF( Newton ) FirstIter = .TRUE.

    END IF



  CONTAINS

!-----------------------------------------------------------------
!> Armijo-GoldStein Criterion for accepting stepsize
!-----------------------------------------------------------------
    FUNCTION ArmijoGoldsteinSearch(Tests,Alpha) RESULT ( ReduceStep )

      INTEGER :: Tests 
      REAL(KIND=dp) :: Alpha
      LOGICAL :: ReduceStep

      ReduceStep = ( Cost(CostMode) > ( 1.0_dp - Myy * Alpha ) * Cost0(CostMode) )
      IF( ReduceStep ) THEN
        Alpha = Alpha * Relaxation
      ELSE
        Cost0 = Cost
        Residual0 = Residual
      END IF

    END FUNCTION ArmijoGoldsteinSearch


!-------------------------------------------------------------------------------
!> Choose next parameter set from 1D bisection search
!-------------------------------------------------------------------------------

    FUNCTION BisectMinimumSearch(Tests, Alpha) RESULT ( ReduceStep ) 

      INTEGER :: Tests 
      REAL(KIND=dp) :: Alpha
      LOGICAL :: ReduceStep
      
      INTEGER :: i,j,k
      REAL(KIND=dp) :: step, p(3),c(3),r(3),raid,beta 
      
      SAVE step, p, c, r

      ReduceStep = .TRUE.
      
      IF(Tests == 1) THEN
        p(1) = 0.0_dp
        c(1) = Cost0(CostMode)
        r(1) = Residual0

        p(2) = 1.0_dp
        c(2) = Cost(CostMode)
        r(2) = Residual
        
        step = 0.25_dp
        Alpha = 0.5_dp
        RETURN
      ELSE 
        p(3) = Alpha
        c(3) = Cost(CostMode) 
        r(3) = Residual
      END IF

      
     ! Order the previous points so that p1 < p2 < p3
      DO k=1,2 
        DO i=k+1,3
          IF(p(i) < p(k)) THEN
            raid = p(k)
            p(k) = p(i)
            p(i) = raid

            raid = c(k)
            c(k) = c(i)
            c(i) = raid

            raid = r(k)
            r(k) = r(i)
            r(i) = raid
          END IF
        END DO
      END DO
      
      IF( Debug ) THEN
        PRINT *,'Bisect p:',p
        PRINT *,'Bisect c:',c
        PRINT *,'Bisect r:',r
      END IF

      ! The value of alpha already known accurately
      IF( MAXVAL(p)-MINVAL(p) < LineTol ) THEN
        ! PRINT *,'cond1'
        ReduceStep = .FALSE.
      END IF

      ! The value of cost function small compared to absolute value of it
      IF( MAXVAL(c)-MINVAL(c) < LineTol * MINVAL( ABS(c) ) ) THEN
        ! PRINT *,'cond2'
        ReduceStep = .FALSE.
      END IF

      ! We can also use the residual as criterion for stopping
      IF( Residual < LineTol * Residual0 ) THEN
        ! PRINT *,'cond3'
        ReduceStep = .FALSE.
      END IF

      ! Of these choose the one with smallest cost
      IF( .NOT. ReduceStep ) THEN
        i = 1
        DO k=2,3
          IF( c(k) < c(i) ) i = k
        END DO

        Alpha = p(i)
        Residual0 = r(i)
        Cost0(CostMode) = c(i)
        ! PRINT *,'Choosing i',i,Alpha,Residual0,Cost0

        RETURN
      END IF


      ! Monotonic line segment
      IF( (c(2)-c(1))*(c(3)-c(2)) > 0.0) THEN
        IF(c(3) < c(1)) THEN
          Alpha = p(3) + SIGN(step,p(3)-p(1))
          c(1) = c(3)
          p(1) = p(3)
          r(1) = r(3)
        ELSE
          Alpha = p(1) + SIGN(step,p(1)-p(3))
        END IF
      ELSE IF(c(2) < c(1) .OR. c(2) < c(3)) THEN 
        IF(c(3) < c(1)) THEN
          c(1) = c(3)
          p(1) = p(3)
          r(1) = r(3)
        END IF
        step = (p(2)-p(1))/2.0d0
        Alpha = p(1) + SIGN(step,p(2)-p(1))
      ELSE  
        IF( Debug ) THEN
          PRINT *,'p:',p
          PRINT *,'c:',c,Cost0
          PRINT *,'r:',r,Residual0
          PRINT *,'dc',c(2)-c(1),c(3)-c(2)
        END IF

        IF( MINVAL ( c ) < Cost0(CostMode) ) THEN
          i = 1
          DO k=2,3
            IF( c(k) < c(i) ) i = k
          END DO
          Alpha = p(i)
          Cost0(CostMode) = c(i)
          Residual0 = r(i)
         
          CALL Warn('BisectSearch','Bisection method improved but faced local maximium')
          ReduceStep = .FALSE.
        ELSE 
          IF( MINVAL ( r ) < Residual0 ) THEN
            CALL Warn('BisectSearch','Bisection method improved but faced local maximium')
          ELSE 
            CALL Warn('BisectSearch','Bisection method cannot handle local maxima')
          END IF

          i = 1
          DO k=2,3
            IF( r(k) < r(i) ) i = k
          END DO
          Alpha = p(i)
          Cost0(CostMode) = c(i)
          Residual0 = r(i)         
        END IF

        ReduceStep = .FALSE.
      END IF

      ! Because alpha should be in limit [0,1] make the corrections
      ! If the orthogonalization is used then we don't have the luxury
      ! of setting the extent as nicely.
      !------------------------------------------------------------
      IF( .NOT. Ortho ) THEN
        beta = alpha
        k = 0
        IF( Alpha < -EPSILON( Alpha ) ) THEN
          IF( p(1) < EPSILON( Alpha ) ) THEN
            step = (p(2)-p(1))/2.0_dp
            Alpha = p(1) + step
            k = 1
          ELSE
            Alpha = 0.0_dp
            k = 1
          END IF
        ELSE IF( Alpha > 1.0_dp + EPSILON( Alpha ) ) THEN
          IF( p(3) > 1.0_dp - EPSILON( Alpha ) ) THEN
            step = (p(3)-p(2))/2.0_dp
            Alpha = p(2) + step
            k = 2
          ELSE
            Alpha = 1.0_dp
            k = 3 
          END IF
        END IF
        
!        IF( ABS( beta-alpha) > TINY(alpha)) PRINT *,'Extent change',Beta,Alpha
      END IF

    END FUNCTION BisectMinimumSearch


!-------------------------------------------------------------------------------
!> Choose next parameter set from 1D bisection search
!-------------------------------------------------------------------------------
    FUNCTION BisectZeroSearch(Tests, Alpha) RESULT ( ReduceStep ) 

      INTEGER :: Tests 
      REAL(KIND=dp) :: Alpha
      LOGICAL :: ReduceStep
      
      INTEGER :: i,j,k
      REAL(KIND=dp) :: step, p(3),c(3),paid,caid,beta 
      
      SAVE step, p, c

      ReduceStep = .TRUE.
      
      IF(Tests == 1) THEN
        p(1) = 0.0_dp
        c(1) = Cost0(CostMode)
        
        p(2) = 1.0_dp
        c(2) = Cost1(CostMode)
        
        IF( Cost0(CostMode) * Cost1(CostMode) > 0.0_dp ) THEN
          CALL Warn('CostSearch','Lumped forces should have different sign!')
        END IF

        Alpha = 0.5_dp
        RETURN
      ELSE 
        p(3) = Alpha
        c(3) = Cost(CostMode) 
      END IF
      
     ! Order the previous points so that p1 < p2 < p3
      DO k=1,2 
        DO i=k+1,3
          IF(p(i) < p(k)) THEN
            paid = p(k)
            p(k) = p(i)
            p(i) = paid
            caid = c(k)
            c(k) = c(i)
            c(i) = caid
          END IF
        END DO
      END DO

      IF( Debug ) THEN
        PRINT *,'Cost p:',p
        PRINT *,'Cost c:',c
      END IF

      IF( p(3)-p(1) < LineTol ) THEN
        ReduceStep = .FALSE.
        RETURN
      END IF

      ! Zero value is between 1st interval
      IF( c(1)*c(2) < 0.0_dp ) THEN
        Alpha = (p(1)+p(2))/2.0_dp
      ELSE IF ( c(2)*c(3) < 0.0_dp ) THEN
        Alpha = (p(2)+p(3))/2.0_dp

        ! We don't need 1st values, but we do need 3rd 
        p(1) = p(3)
        c(1) = c(3)
      ELSE
        CALL Fatal('ForceSearch','Lumped forces should have different sign!')
      END IF
      
    END FUNCTION BisectZeroSearch

!------------------------------------------------------------------------------
  END FUNCTION CheckStepSize
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Apply Anderson acceleration to the solution of nonlinear system.
!> Also may apply acceleration to the linear system. 
!------------------------------------------------------------------------------
  SUBROUTINE NonlinearAcceleration(A,x,b,Solver,PreSolve,NoSolve)    
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) CONTIG :: b(:),x(:)
    TYPE(Solver_t) :: Solver
    LOGICAL :: PreSolve
    LOGICAL, OPTIONAL :: NoSolve
    !------------------------------------------------------------------------------
    ! We have a special stucture for the iterates and residuals so that we can
    ! cycle over the pointers instead of the values. 
    TYPE AndersonVect_t
      LOGICAL :: Additive
      REAL(KIND=dp), POINTER :: Iterate(:), Residual(:), Ax(:)
      INTEGER :: tag
    END TYPE AndersonVect_t
    TYPE(AndersonVect_t), ALLOCATABLE :: AndersonBasis(:), AndersonTmp
    INTEGER :: AndersonInterval, ItersCnt, AndersonVecs, VecsCnt, iter, n,i,j,k
    TYPE(Variable_t), POINTER :: iterV, Svar
    REAL(KIND=dp), ALLOCATABLE :: Alphas(:),AxTable(:,:),TmpVec(:) 
    REAL(KIND=dp) :: Nrm, AndersonRelax
    LOGICAL :: Found, DoRelax, KeepBasis, Visited = .FALSE., Parallel    
    INTEGER :: LinInterval
    INTEGER :: PrevSolverId = -1
    
    SAVE AndersonBasis, TmpVec, Alphas, ItersCnt, AndersonInterval, VecsCnt, AndersonVecs, &
        PrevSolverId, AxTable, AndersonRelax, DoRelax, Visited, KeepBasis, LinInterval
        
    IF( PreSolve ) THEN
      CALL Info('NonlinearAcceleration','Performing pre-solution steps',Level=8)
    ELSE
      CALL Info('NonlinearAcceleration','Performing post-solution steps',Level=8)
    END IF

    Parallel = ( ParEnv % PEs > 1 ) 
        
    iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
    iter = NINT(iterV % Values(1))

    IF(PRESENT(NoSolve)) NoSolve = .FALSE.
    
    n = A % NumberOfRows
          
    IF(.NOT. Visited ) THEN
      PrevSolverId = Solver % SolverId
      CALL Info('NonlinearAcceleration','Allocating structures for solution history',Level=6)

      AndersonInterval = ListGetInteger( Solver % Values,&
          'Nonlinear System Acceleration Interval',Found)      
      LinInterval = ListGetInteger( Solver % Values,&
          'Linear System Acceleration Interval',Found)      

      AndersonVecs = MAX( AndersonInterval, LinInterval )
      IF( AndersonVecs == 0 ) THEN
        CALL Fatal('NonlinearAcceleration','Both acceleration intervals are zero!')
      END IF
            
      AndersonRelax = ListGetCReal( Solver % Values,&
          'Nonlinear System Acceleration Relaxation',DoRelax)
      KeepBasis = ListGetLogical( Solver % Values,&
          'Nonlinear System Acceleration Keep Vectors',Found)            

      ItersCnt = 0    ! relates to "AndersonInterval"
      VecsCnt = 0     ! relates to "AndersonVecs"
      
      IF(.NOT. ALLOCATED( AndersonBasis ) ) THEN
        ALLOCATE( AndersonBasis(AndersonVecs) )
        DO i=1,AndersonVecs
          ALLOCATE( AndersonBasis(i) % Residual(n), &
              AndersonBasis(i) % Iterate(n) )
          AndersonBasis(i) % Residual = 0.0_dp
          AndersonBasis(i) % Iterate = 0.0_dp          
        END DO
        ALLOCATE( TmpVec(n), Alphas(AndersonVecs) )
      END IF
      Visited = .TRUE.
    END IF
    
    IF( PrevSolverId /= Solver % SolverId ) THEN
      CALL Fatal('NonlinearAcceleration','Current implementation only supports one solver!')
    END IF
      
    
    IF( PreSolve ) THEN           
      IF( iter == 1 ) THEN
        ItersCnt = 0
        IF( .NOT. KeepBasis ) VecsCnt = 0
      END IF

      ItersCnt = ItersCnt + 1
      VecsCnt = VecsCnt + 1

      ! Calculate the residual of the matrix equation
      ! Here 'x' comes before being modified hence A(x) is consistent. 
      CALL MatrixVectorMultiply( A, x, TmpVec )
      TmpVec = TmpVec - b

      ! Add the iterate and residual to the basis vectors.
      ! This is fast as we operate with pointers mainly.
      AndersonTmp = AndersonBasis(AndersonVecs)        
      DO i=AndersonVecs,2,-1
        AndersonBasis(i) = AndersonBasis(i-1)
      END DO
      AndersonBasis(1) = AndersonTmp
      AndersonBasis(1) % Residual = TmpVec
      AndersonBasis(1) % Iterate = x 

      ! Pure Anderson sweep is done every AndersonInterval iterations if we have full basis.
      IF(.NOT. DoRelax .AND. AndersonInterval > 0 ) THEN
        IF( VecsCnt >= AndersonVecs .AND. ItersCnt >= AndersonInterval ) THEN
          CALL AndersonMinimize( )
          ItersCnt = 0
          IF(PRESENT(NoSolve)) NoSolve = .TRUE.
          RETURN
        END IF
      END IF
      
      IF( LinInterval > 0 ) THEN
        CALL AndersonGuess()
      END IF
    ELSE
      ! Relaxation strategy is done after each linear solve.
      IF( DoRelax ) THEN
        CALL Info('NonlinearAcceleration','Minimizing residual using history data',Level=6)
        CALL AndersonMinimize( )
      END IF
    END IF

  CONTAINS 


    !------------------------------------------------------------------------------
    FUNCTION Mydot( n, x, y ) RESULT(s)
      !------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp)  :: s
      REAL(KIND=dp) CONTIG :: x(:)
      REAL(KIND=dp) CONTIG, OPTIONAL :: y(:)
      !------------------------------------------------------------------------------
      IF ( .NOT. Parallel ) THEN
        IF( PRESENT( y ) ) THEN
          s = DOT_PRODUCT( x(1:n), y(1:n) )
        ELSE
          s = DOT_PRODUCT( x(1:n), x(1:n) )
        END IF
      ELSE
        IF( PRESENT( y ) ) THEN
          s = ParallelDot( n, x, y )
        ELSE
          s = ParallelDot( n, x, x )
        END IF
      END IF
      !------------------------------------------------------------------------------
    END FUNCTION Mydot
    !------------------------------------------------------------------------------


    !------------------------------------------------------------------------------
    SUBROUTINE Mymv( A, x, b, Update )
      !------------------------------------------------------------------------------
      REAL(KIND=dp) CONTIG :: x(:), b(:)
      TYPE(Matrix_t), POINTER :: A
      LOGICAL, OPTIONAL :: Update
      !------------------------------------------------------------------------------
      IF ( .NOT. Parallel ) THEN
        CALL CRS_MatrixVectorMultiply( A, x, b )
      ELSE
        IF ( PRESENT( Update ) ) THEN
          CALL ParallelMatrixVector( A,x,b,Update,ZeroNotOwned=.TRUE. )
        ELSE
          CALL ParallelMatrixVector( A,x,b,ZeroNotOwned=.TRUE. )
        END IF
      END IF
      !------------------------------------------------------------------------------
    END SUBROUTINE Mymv
    !------------------------------------------------------------------------------

    
    ! Given set of basis vectors and residuals find a new suggestion for solution.
    ! Either use as such or combine it to solution when relaxation is used.
    ! This is applied to boost nonlinear iteration. 
    !------------------------------------------------------------------------------
    SUBROUTINE AndersonMinimize()
      INTEGER ::m, n, AndersonMinn
      REAL(KIND=dp) :: rr, rb
      
      m = MIN( ItersCnt, AndersonInterval )      
      
      AndersonMinN = ListGetInteger( Solver % Values,&
          'Nonlinear System Acceleration First Iteration',Found )
      IF(.NOT. (Found .OR. DoRelax)) AndersonMinN = AndersonInterval
            
      ! Nothing to do 
      IF( m < AndersonMinN ) RETURN
      
      ! If size of our basis is just one, there is not much to do...
      ! We can only perform classical relaxation. 
      IF( m == 1 ) THEN
        x = AndersonRelax * x + (1-AndersonRelax) * AndersonBasis(1) % Iterate
        RETURN
      END IF
      
      ! If we are converged then the solution should already be the 1st component.
      ! Hence use that as the basis. 
      Alphas(1) = 1.0_dp     
      TmpVec = AndersonBasis(1) % Residual
      
      ! Minimize the residual
      n = SIZE( AndersonBasis(1) % Residual ) 
      DO k=2,m
        rr = MyDot( n, AndersonBasis(k) % Residual ) 
        rb = MyDot( n, AndersonBasis(k) % Residual, TmpVec )         
        Alphas(k) = -rb / rr 
        TmpVec = TmpVec + Alphas(k) * AndersonBasis(k) % Residual
      END DO

      ! Normalize the coefficients such that the sum equals unity
      ! This way for example, Dirichlet BCs will be honored.
      Alphas = Alphas / SUM( Alphas(1:m) )

      IF( InfoActive(10) ) THEN
        DO i=1,m
          WRITE(Message,'(A,I0,A,ES12.3)') 'Alpha(',i,') = ',Alphas(i)
          CALL Info('NonlinearAcceleration',Message)
        END DO
      END IF
              
      ! Create the new suggestion for the solution vector
      ! We take part of the suggested new solution vector 'x' and
      ! part of minimized residual that was used in anderson acceleration.
      IF( DoRelax ) THEN
        Alphas = Alphas * (1-AndersonRelax)
        x = AndersonRelax * x
        DO k=1,m
          x = x + Alphas(k) * AndersonBasis(k) % Iterate
        END DO
      ELSE
        x = Alphas(1) * AndersonBasis(1) % Iterate
        DO k=2,m
          x = x + Alphas(k) * AndersonBasis(k) % Iterate
        END DO
      END IF
        
    END SUBROUTINE AndersonMinimize

    
    ! Given set of basis vectors and a linear system
    ! find a combincation of the vectors that minimizes the norm of the linear
    ! system. This may be used to provide a better initial guess for a linear system.
    !--------------------------------------------------------------------------------
    SUBROUTINE AndersonGuess()
      INTEGER :: AndersonMinN

      REAL(KIND=dp), POINTER, SAVE ::Betas(:), Ymat(:,:)
      LOGICAL, SAVE :: AllocationsDone = .FALSE.
      INTEGER :: i,j,m
      
      IF(.NOT. AllocationsDone ) THEN
        m = LinInterval
        DO i=1,LinInterval
          ALLOCATE( AndersonBasis(i) % Ax(n) )
          AndersonBasis(i) % Ax = 0.0_dp
        END DO
        ALLOCATE(Betas(m),Ymat(m,m))
        AllocationsDone = .TRUE.
      END IF
      
      m = MIN( VecsCnt, LinInterval )      

      ! Calculate the residual of the matrix equation
      DO i=1,m
        CALL Mymv( A, AndersonBasis(i) % Iterate, AndersonBasis(i) % Ax )
      END DO

      DO i=1,m
        DO j=i,m
          Ymat(i,j) = SUM( AxTable(:,i) * AxTable(:,j) )
          Ymat(j,i) = Ymat(i,j)
        END DO
        Betas(i) = SUM( AxTable(:,i) * b )
      END DO
      
      CALL LUSolve(m, YMat(1:m,1:m), Betas(1:m) )

      IF( InfoActive(10) ) THEN
        DO i=1,m
          WRITE(Message,'(A,I0,A,ES12.3)') 'Beta(',i,') = ',Betas(i)
          CALL Info('LinearAcceleration',Message)
        END DO
      END IF
                                
      x = Betas(m) * AndersonBasis(m) % Iterate
      DO i=1,m-1
        x = x + Betas(i) * AndersonBasis(i) % Iterate
      END DO

    END SUBROUTINE AndersonGuess
    
  END SUBROUTINE NonlinearAcceleration
!------------------------------------------------------------------------------

  

  !> Calculate the number of separature pieces in a serial mesh.
  !> This could be used to detect problems in mesh when suspecting
  !> floating parts not fixed by any BC, for example.
  !---------------------------------------------------------------------------------
  SUBROUTINE CalculateMeshPieces( Mesh )

    TYPE(Mesh_t), POINTER :: Mesh

    LOGICAL :: Ready
    INTEGER :: i,j,k,n,t,MinIndex,MaxIndex,Loop,NoPieces
    INTEGER, ALLOCATABLE :: MeshPiece(:),PiecePerm(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: Indexes(:)
    TYPE(Variable_t), POINTER :: Var

    IF( ParEnv % PEs > 1 ) THEN
      CALL Warn('CalculateMeshPieces','Implemented only for serial meshes, doing nothing!')
      RETURN
    END IF

    n = Mesh % NumberOfNodes
    ALLOCATE( MeshPiece( n ) ) 
    MeshPiece = 0

    ! Only set the piece for the nodes that are used by some element
    ! For others the marker will remain zero. 
    DO t = 1, Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)        
      Indexes => Element % NodeIndexes
      MeshPiece( Indexes ) = 1
    END DO    
    j = 0
    DO i = 1, n
      IF( MeshPiece(i) > 0 ) THEN
        j = j + 1
        MeshPiece(i) = j
      END IF
    END DO

    CALL Info('CalculateMeshPieces','Number of non-body nodes in mesh is '//TRIM(I2S(n-j)))
    
    ! We go through the elements and set all the piece indexes to minimimum index
    ! until the mesh is unchanged. Thereafter the whole piece will have the minimum index
    ! of the piece.
    Ready = .FALSE.
    Loop = 0
    DO WHILE(.NOT. Ready) 
      Ready = .TRUE.
      DO t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)        
        Indexes => Element % NodeIndexes

        MinIndex = MINVAL( MeshPiece( Indexes ) )
        MaxIndex = MAXVAL( MeshPiece( Indexes ) )
        IF( MaxIndex > MinIndex ) THEN
          MeshPiece( Indexes ) = MinIndex
          Ready = .FALSE.
        END IF
      END DO
      Loop = Loop + 1
    END DO
    CALL Info('CalculateMeshPieces','Mesh coloring loops: '//TRIM(I2S(Loop)))

    ! If the maximum index is one then for sure there is only one body
    IF( MaxIndex == 1 ) THEN
      CALL Info('CalculateMeshPieces','Mesh consists of single body!')
      RETURN
    END IF

    ! Compute the true number of different pieces
    ALLOCATE( PiecePerm( MaxIndex ) ) 
    PiecePerm = 0
    NoPieces = 0
    DO i = 1, n
      j = MeshPiece(i) 
      IF( j == 0 ) CYCLE
      IF( PiecePerm(j) == 0 ) THEN
        NoPieces = NoPieces + 1
        PiecePerm(j) = NoPieces 
      END IF
    END DO
    CALL Info('CalculateMeshPieces','Number of separate pieces in mesh is '//TRIM(I2S(NoPieces)))


    ! Save the mesh piece field to > mesh piece < 
    Var => VariableGet( Mesh % Variables,'Mesh Piece' )
    IF(.NOT. ASSOCIATED( Var ) ) THEN      
      CALL VariableAddVector ( Mesh % Variables,Mesh, CurrentModel % Solver,'Mesh Piece' )
      Var => VariableGet( Mesh % Variables,'Mesh Piece' )
    END IF

    IF( .NOT. ASSOCIATED( Var ) ) THEN
      CALL Fatal('CalculateMeshPieces','Could not get handle to variable > Mesh Piece <')
    END IF

    DO i = 1, n
      j = Var % Perm( i ) 
      IF( j == 0 ) CYCLE
      IF( MeshPiece(i) > 0 ) THEN
        Var % Values( j ) = 1.0_dp * PiecePerm( MeshPiece( i ) )
      ELSE
        Var % Values( j ) = 0.0_dp
      END IF
    END DO
    CALL Info('CalculateMeshPieces','Saving mesh piece field to: mesh piece')
  
  END SUBROUTINE CalculateMeshPieces



  ! Create a boundary matrix and at calculate step compute the boundary loads
  ! for one given body. This is not called by default but the user needs to
  ! include it in the code, both at assembly and after solution.
  !-----------------------------------------------------------------------------
  SUBROUTINE BCLoadsAssembly( Solver, Element, LocalMatrix, LocalRhs )

    TYPE(Solver_t) :: Solver
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: LocalMatrix(:,:)
    REAL(KIND=dp) :: LocalRhs(:)

    LOGICAL :: FirstStep
    INTEGER :: i,j,k,l,n,Row,Col,Dofs,ElemNo,TargetBody=-1
    TYPE(Matrix_t), POINTER :: BCMat
    REAL(KIND=dp) :: Val
    LOGICAL :: Found
    INTEGER, POINTER :: Perm(:), BCPerm(:)
    CHARACTER(MAX_NAME_LEN) :: Name   
    TYPE(Variable_t), POINTER :: BCVar


    SAVE :: BCMat, TargetBody, BCPerm, Perm, Dofs


    FirstStep = ( Solver % ActiveElements(1) == Element % ElementIndex )

    IF( FirstStep ) THEN
      CALL Info('BCLoadsAssembly','Visiting first element')
 
      BCMat => Solver % Matrix % EMatrix
      IF(.NOT. ASSOCIATED( BCMat ) ) THEN
        TargetBody = ListGetInteger( Solver % Values,'Boundary Loads Target Body',Found )
        IF( Found ) THEN
          CALL Info('BCLoadsAssembly','Target body set to: '//TRIM(I2S(TargetBody)))       
        ELSE
          TargetBody = -1
          RETURN
        END IF

        CALL Info('BCLoadsAssembly','Allocating structures for load computation')
        IF ( ParEnv % PEs > 1 ) THEN
          CALL Warn('BCLoadsAssembly','Not implemented in parallel')
        END IF

        ! Mark the active nodes
        ALLOCATE( BCPerm( Solver % Mesh % NumberOfNodes ) )
        BCPerm = 0

        ElemNo = 0
        k = Solver % Mesh % NumberOfBulkElements
        DO i = k+1,k + Solver % Mesh % NumberOfBoundaryElements
          Element => Solver % Mesh % Elements(i)
          Found = .FALSE.             
          IF( ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
            Found = ( Element % BoundaryInfo % Left % BodyId == TargetBody )
          END IF
          IF(.NOT. Found ) THEN
            IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
              Found = ( Element % BoundaryInfo % Right % BodyId == TargetBody )
            END IF
          END IF
          IF( Found ) THEN
            ElemNo = ElemNo + 1
            BCPerm( Element % NodeIndexes ) = 1
          END IF
        END DO

        CALL Info('BCLoadsAssembly','Number of related boundary elements: '//TRIM(I2S(ElemNo)))

        n = 0
        DO i=1,Solver % Mesh % NumberOfNodes
          IF( BCPerm(i) > 0 ) THEN
            n = n + 1
            BCPerm(i) = n
          END IF
        END DO
        CALL Info('BCLoadsAssembly','Number of active nodes: '//TRIM(I2S(n)))

        ! Create the list matrix 
        BCMat => AllocateMatrix()
        BCMat % Format = MATRIX_LIST           
        CALL AddToMatrixElement( BCMat, n, n, 0.0_dp )
        Solver % Matrix % EMatrix => BCMat

        ALLOCATE( BCMat % Rhs(n) )
        BCMat % Rhs = 0.0_dp
      END IF

      ! When visiting the routine after the 1st iteration the matrix for is already CRS 
      IF( BCMat % Format == MATRIX_CRS ) THEN
        BCMat % Values = 0.0_dp
        BCMat % Rhs = 0.0_dp
      END IF

      Dofs = Solver % Variable % Dofs
      Perm => Solver % Variable % Perm

      Name = TRIM(Solver % Variable % Name)//' BCLoads'
      BCVar => VariableGet( Solver % Mesh % Variables, TRIM( Name ) )
      IF(.NOT. ASSOCIATED( BCVar ) ) THEN
        CALL Info('CalculateBCLoads','Creating variable: '//TRIM(Name))
        CALL VariableAddVector( Solver % Mesh % Variables,&
            Solver % Mesh, Solver, Name, DOFs, Perm = BCPerm )
      END IF
      
    END IF

    IF( Element % BodyId == TargetBody ) THEN       
      n = Element % TYPE % NumberOfNodes
      DO i=1,n
        IF ( BCPerm( Element % NodeIndexes(i) ) == 0 ) CYCLE
        DO k=0,Dofs-1
          Row = Dofs * BCPerm( Element % NodeIndexes(i) ) - k
          BCMat % rhs(Row) = BCMat % rhs(Row) + LocalRhs(Dofs*i-k)
          DO j=1,n
            DO l=0,Dofs-1
              Col = Dofs * Perm( Element % NodeIndexes(j) ) - l
              Val = LocalMatrix(Dofs*i-k,Dofs*j-l)
              CALL AddToMatrixElement(BCMat,Row,Col,Val)
            END DO
          END DO
        END DO
      END DO
    END IF


  END SUBROUTINE BCLoadsAssembly


  ! Calculate the boundary loads resulting from the action of boundary matrix.
  !-----------------------------------------------------------------------------
  SUBROUTINE BCLoadsComputation( Solver )

    TYPE(Solver_t) :: Solver

    TYPE(Matrix_t), POINTER :: BCMat
    CHARACTER(MAX_NAME_LEN) :: Name   
    TYPE(Variable_t), POINTER :: BCVar


    BCMat => Solver % Matrix % EMatrix
    IF(.NOT. ASSOCIATED( BCMat ) ) THEN
      CALL Fatal('BCLoadsComputation','We should have the boundary matrix!')
    END IF
        
    CALL Info('CalculateBCLoads','Computing boundary loads')
    IF( BCMat % FORMAT == MATRIX_LIST ) THEN
      CALL List_ToCRSMatrix( BCMat )
      CALL Info('CalculateBCLoads','Matrix format changed to CRS')
    END IF

    Name = TRIM(Solver % Variable % Name)//' BCLoads'
    BCVar => VariableGet( Solver % Mesh % Variables, TRIM( Name ) )
    IF(.NOT. ASSOCIATED( BCVar ) ) THEN
      CALL Fatal('CalculateBCLoads','Variable not present: '//TRIM(Name))
    END IF
    
    CALL MatrixVectorMultiply( BCMat, Solver % Variable % Values, BCVar % Values )
    BCVar % Values = BCVar % Values - BCMat % rhs

    CALL Info('CalculateBCLoads','All done')

  END SUBROUTINE BCLoadsComputation


    
!------------------------------------------------------------------------------
!> Prints the values of the CRS matrix to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintMatrix( A, Parallel, CNumbering,SaveMass, SaveDamp, SaveStiff)
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A            !< Structure holding matrix
    LOGICAL :: Parallel    !< are we in parallel mode?
    LOGICAL :: CNumbering  !< Continuous numbering ?
    LOGICAL, OPTIONAL :: SaveMass  !< Should we save the mass matrix
    LOGICAL, OPTIONAL :: SaveDamp  !< Should we save the damping matrix
    LOGICAL, OPTIONAL :: SaveStiff !< Should we save the stiffness matrix
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,IndMass,IndDamp,IndStiff,IndMax,row,col
    LOGICAL :: DoMass, DoDamp, DoStiff, Found
    REAL(KIND=dp) :: Vals(3)
    INTEGER, ALLOCATABLE :: Owner(:)

    DoMass = .FALSE.
    IF( PRESENT( SaveMass ) ) DoMass = SaveMass
    IF( DoMass .AND. .NOT. ASSOCIATED( A % MassValues ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting mass matrix')
      DoMass = .FALSE. 
    END IF

    DoDamp = .FALSE.
    IF( PRESENT( SaveDamp ) ) DoDamp = SaveDamp
    IF( DoDamp .AND. .NOT. ASSOCIATED( A % DampValues ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting damp matrix')
      DoDamp = .FALSE. 
    END IF

    DoStiff = .TRUE.
    IF( PRESENT( SaveStiff ) ) DoStiff = SaveStiff
    IF( DoStiff .AND. .NOT. ASSOCIATED( A % Values ) ) THEN
      CALL Warn('CRS_PrintMatrix','Cannot save nonexisting stiff matrix')
      DoStiff = .FALSE. 
    END IF

    IF(.NOT. (DoStiff .OR. DoDamp .OR. DoMass ) ) THEN
      CALL Warn('CRS_PrintMatrix','Saving just the topology!')
    END IF
    
    IndStiff = 0
    IndDamp = 0
    IndMass = 0

    IF( DoStiff ) IndStiff = 1
    IF( DoDamp ) IndDamp = IndStiff + 1
    IF( DoMass ) IndMass = MAX( IndStiff, IndDamp ) + 1
    IndMax = MAX( IndStiff, IndDamp, IndMass )

    IF (Parallel.AND.Cnumbering) THEN
      n = SIZE(A % ParallelInfo % GlobalDOFs)
  
      ALLOCATE( A % Gorder(n), Owner(n) )
      CALL ContinuousNumbering( A % ParallelInfo, &
          A % Perm, A % Gorder, Owner )
    END IF

    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF
      DO j = A % Rows(i),A % Rows(i+1)-1

        col = A % Cols(j)
        IF(Parallel) THEN
          IF(Cnumbering) THEN
            col = A % Gorder(col)
          ELSE 
            col = A % ParallelInfo % GlobalDOFs(col)
          END IF
        END IF

        WRITE(1,'(I0,A,I0,A)',ADVANCE='NO') row,' ',col,' '

        IF( DoStiff ) THEN
          Vals(IndStiff) = A % Values(j)
        END IF
        IF( DoDamp ) THEN
          Vals(IndDamp) = A % DampValues(j)
        END IF
        IF( DoMass ) THEN
          Vals(IndMass) = A % MassValues(j)
        END IF

        IF( IndMax > 0 ) THEN
          WRITE(1,*) Vals(1:IndMax)          
        ELSE
          WRITE(1,'(A)') ' '
        END IF
      END DO
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE  PrintMatrix
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Prints the values of the right-hand-side vector to standard output.
!------------------------------------------------------------------------------
  SUBROUTINE PrintRHS( A, Parallel, CNumbering )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A  !< Structure holding matrix
    LOGICAL :: Parallel, CNumbering
!------------------------------------------------------------------------------
    INTEGER :: i, row
    REAL(KIND=dp) :: Val

    DO i=1,A % NumberOfRows
      row = i
      IF(Parallel) THEN
        IF(Cnumbering) THEN
          row = A % Gorder(i)
        ELSE 
          row = A % ParallelInfo % GlobalDOFs(i)
        END IF
      END IF

      Val = A % Rhs(i)
      WRITE(1,'(I0,A)',ADVANCE='NO') row,' '
      IF( ABS( Val ) <= TINY( Val ) ) THEN
        WRITE(1,'(A)') '0.0'
      ELSE
        WRITE(1,*) Val
      END IF
    END DO

  END SUBROUTINE PrintRHS
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Solves a linear system and also calls the necessary preconditioning routines.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SolveLinearSystem( A, b, &
       x, Norm, DOFs, Solver, BulkMatrix )
!------------------------------------------------------------------------------
    USE EigenSolve

    REAL(KIND=dp) CONTIG :: b(:), x(:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: A
    INTEGER :: DOFs
    TYPE(Solver_t), TARGET :: Solver
    TYPE(Matrix_t), OPTIONAL, POINTER :: BulkMatrix
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, NodalLoads
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: Relax,GotIt,Stat,ScaleSystem, EigenAnalysis, HarmonicAnalysis,&
               BackRotation, ApplyRowEquilibration, ApplyLimiter, Parallel, &
               SkipZeroRhs, ComplexSystem, ComputeChangeScaled, ConstraintModesAnalysis, &
               RecursiveAnalysis, CalcLoads
    INTEGER :: n,i,j,k,l,ii,m,DOF,istat,this,mn
    CHARACTER(LEN=MAX_NAME_LEN) :: Method, Prec, ProcName, SaveSlot
    INTEGER(KIND=AddrInt) :: Proc
    REAL(KIND=dp), ALLOCATABLE, TARGET :: Px(:), &
                TempVector(:), TempRHS(:), NonlinVals(:)
    REAL(KIND=dp), POINTER :: Diag(:)
    REAL(KIND=dp) :: s,Relaxation,Beta,Gamma,bnorm,Energy,xn,bn
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Matrix_t), POINTER :: Aaid, Projector, MP
    REAL(KIND=dp), POINTER :: mx(:), mb(:), mr(:)
    TYPE(Variable_t), POINTER :: IterV
    LOGICAL :: NormalizeToUnity, AndersonAcc, AndersonScaled, NoSolve
    
    INTERFACE 
       SUBROUTINE VankaCreate(A,Solver)
          USE Types
          TYPE(Matrix_t) :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE VankaCreate

       SUBROUTINE CircuitPrecCreate(A,Solver)
          USE Types
          TYPE(Matrix_t), TARGET :: A
          TYPE(Solver_t) :: Solver
       END SUBROUTINE CircuitPrecCreate

       SUBROUTINE FetiSolver(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
       END SUBROUTINE FetiSolver
 
       SUBROUTINE BlockSolveExt(A,x,b,Solver)
          USE Types
          TYPE(Matrix_t), POINTER :: A
          TYPE(Solver_t) :: Solver
          REAL(KIND=dp) :: x(:), b(:)
       END SUBROUTINE BlockSolveExt
    END INTERFACE
!------------------------------------------------------------------------------

    Params => Solver % Values
 
    ComplexSystem = ListGetLogical( Params, 'Linear System Complex', GotIt )
    IF ( GotIt ) A % COMPLEX = ComplexSystem
    
    ScaleSystem = ListGetLogical( Params, 'Linear System Scaling', GotIt )
    IF ( .NOT. GotIt  ) ScaleSystem = .TRUE.
    
    IF( ListGetLogical( Params,'Linear System Skip Complex',GotIt ) ) THEN
      CALL Info('SolveLinearSystem','This time skipping complex treatment',Level=20)
      A % COMPLEX = .FALSE.
      ComplexSystem = .FALSE.
    END IF

    IF( ListGetLogical( Params,'Linear System Skip Scaling',GotIt ) ) THEN     
      CALL Info('SolveLinearSystem','This time skipping scaling',Level=20)
      ScaleSystem = .FALSE.
    END IF
   
    IF( A % COMPLEX ) THEN
      CALL Info('SolveLinearSystem','Assuming complex valued linear system',Level=6)
    ELSE
      CALL Info('SolveLinearSystem','Assuming real valued linear system',Level=8)
    END IF

!------------------------------------------------------------------------------
!   If parallel execution, check for parallel matrix initializations
!------------------------------------------------------------------------------
    IF ( ParEnv % Pes>1.AND..NOT. ASSOCIATED(A % ParMatrix) ) THEN
      CALL ParallelInitMatrix( Solver, A )
    END IF

    IF ( ListGetLogical( Solver % Values, 'Linear System Save',GotIt )) THEN
      saveslot = ListGetString( Solver % Values,'Linear System Save Slot', GotIt )
      IF(SaveSlot == 'linear solve') CALL SaveLinearSystem( Solver, A )
    END IF

!------------------------------------------------------------------------------

    n = A % NumberOfRows

    BackRotation = ListGetLogical(Params,'Back Rotate N-T Solution',GotIt)
    IF (.NOT.GotIt) BackRotation=.TRUE.
    BackRotation = BackRotation .AND. ASSOCIATED(Solver % Variable % Perm)

    IF ( Solver % Matrix % Lumped .AND. Solver % TimeOrder == 1 ) THEN
       Method = ListGetString( Params, 'Timestepping Method', GotIt)
       IF (  Method == 'runge-kutta' .OR. Method == 'explicit euler' ) THEN
         ALLOCATE(Diag(n), TempRHS(n))

         TempRHS= b(1:n)
         Diag = A % Values(A % Diag)

         IF(ParEnv % Pes>1) THEN
           CALL ParallelSumVector(A,Diag)
           CALL ParallelSumVector(A,TempRHS)
         END IF

         DO i=1,n
            IF ( ABS(Diag(i)) /= 0._dp ) x(i) = TempRHS(i) / Diag(i)
         END DO

         DEALLOCATE(Diag, TempRHS)

         IF (BackRotation) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
         Norm = ComputeNorm(Solver, n, x) 
         RETURN
       END IF
    END IF
    
!------------------------------------------------------------------------------
!  These definitions are needed if chanching the iterative solver on-the-fly

    Solver % MultiGridSolver = ( ListGetString( Params, &
        'Linear System Solver', GotIt ) == 'multigrid' )
    Solver % MultiGridTotal = MAX( Solver % MultiGridTotal, &
        ListGetInteger( Params,'MG Levels', GotIt, minv=1 ) )
    Solver % MultiGridTotal = MAX( Solver % MultiGridTotal, &
        ListGetInteger( Params,'Multigrid Levels', GotIt, minv=1 ) )
    Solver % MultiGridLevel = Solver % MultigridTotal
!------------------------------------------------------------------------------

    EigenAnalysis = Solver % NOFEigenValues > 0 .AND. &
        ListGetLogical( Params, 'Eigen Analysis',GotIt )
    
    ConstraintModesAnalysis = ListGetLogical( Params, &
        'Constraint Modes Analysis',GotIt )

    HarmonicAnalysis = ( Solver % NOFEigenValues > 0 ) .AND. &
        ListGetLogical( Params, 'Harmonic Analysis',GotIt )

    ! These analyses types may require recursive strategies and may also have zero rhs
    RecursiveAnalysis = HarmonicAnalysis .OR. EigenAnalysis .OR. ConstraintModesAnalysis


    ApplyLimiter = ListGetLogical( Params,'Apply Limiter',GotIt ) 
    SkipZeroRhs = ListGetLogical( Params,'Skip Zero Rhs Test',GotIt )
#ifdef HAVE_FETI4I
    IF ( C_ASSOCIATED(A % PermonMatrix) ) THEN
      ScaleSystem = .FALSE.
      SkipZeroRhs = .TRUE.
    END IF
#endif

    IF ( .NOT. ( RecursiveAnalysis .OR. ApplyLimiter .OR. SkipZeroRhs ) ) THEN
      bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))      
      IF ( bnorm <= TINY( bnorm) ) THEN
        CALL Info('SolveSystem','Solution trivially zero!')
        x = 0.0d0

        ! Increase the nonlinear counter since otherwise some stuff may stagnate
        ! Normally this is done within ComputeChange
        iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
        Solver % Variable % NonlinIter = iterV % Values(1)
        iterV % Values(1) = iterV % Values(1) + 1 
        Solver % Variable % Norm = 0.0_dp
        Solver % Variable % NonlinConverged = 1
     
        RETURN
      END IF
    END IF

    IF ( Solver % MultiGridLevel == -1  ) RETURN

    ! Set the flags to false to allow recursive strategies for these analysis types, little dirty...
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.FALSE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.FALSE.)
      IF( ConstraintModesAnalysis ) CALL ListAddLogical( Solver % Values,'Constraint Modes Analysis',.FALSE.)
    END IF


!------------------------------------------------------------------------------
!   If solving harmonic analysis go there:
!   --------------------------------------
    IF ( HarmonicAnalysis ) THEN
      CALL SolveHarmonicSystem( A, Solver )
    END IF


!   If solving eigensystem go there:
!   --------------------------------
    IF ( EigenAnalysis ) THEN
      IF ( ScaleSystem ) CALL ScaleLinearSystem(Solver, A )

      CALL SolveEigenSystem( &
          A, Solver %  NOFEigenValues, &
          Solver % Variable % EigenValues,       &
          Solver % Variable % EigenVectors, Solver )
      
      IF ( ScaleSystem ) CALL BackScaleLinearSystem( Solver, A, EigenScaling = .TRUE. ) 
      IF ( BackRotation ) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )

      Norm = ComputeNorm(Solver,n,x)
      Solver % Variable % Norm = Norm
      
      NormalizeToUnity = ListGetLogical( Solver % Values, &
          'Eigen System Normalize To Unity',Stat)         

      IF(NormalizeToUnity .OR. ListGetLogical( Solver % Values,  &
          'Eigen System Mass Normalize', Stat ) ) THEN

        CALL ScaleEigenVectors( A, Solver % Variable % EigenVectors, &
            SIZE(Solver % Variable % EigenValues), NormalizeToUnity ) 
      END IF
      
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
          Solver % Variable % Name )
    END IF


!   If solving constraint modes analysis go there:
!   ----------------------------------------------
    IF ( ConstraintModesAnalysis ) THEN      
      CALL SolveConstraintModesSystem( A, x, b , Solver )
     
      IF ( BackRotation ) CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
      
      Norm = ComputeNorm(Solver,n,x)
      Solver % Variable % Norm = Norm
      
      CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
          Solver % Variable % Name )
    END IF
   
    
    ! We have solved {harmonic,eigen,constraint} system and no need to continue further
    IF( RecursiveAnalysis ) THEN
      IF( HarmonicAnalysis ) CALL ListAddLogical( Solver % Values,'Harmonic Analysis',.TRUE.)
      IF( EigenAnalysis ) CALL ListAddLogical( Solver % Values,'Eigen Analysis',.TRUE.)
      IF( ConstraintModesAnalysis ) CALL ListAddLogical( Solver % Values,'Constraint Modes Analysis',.TRUE.)
      RETURN
    END IF


! Check whether b=0 since then equation Ax=b has only the trivial solution, x=0. 
! In case of a limiter one still may need to check the limiter for contact.
!-----------------------------------------------------------------------------
    bnorm = SQRT(ParallelReduction(SUM(b(1:n)**2)))
    IF ( bnorm <= TINY( bnorm) .AND..NOT.SkipZeroRhs) THEN
      CALL Info('SolveSystem','Solution trivially zero!')
      x = 0.0d0

      ! Increase the nonlinear counter since otherwise some stuff may stagnate
      ! Normally this is done within ComputeChange
      iterV => VariableGet( Solver % Mesh % Variables, 'nonlin iter' )
      Solver % Variable % NonlinIter = iterV % Values(1)
      iterV % Values(1) = iterV % Values(1) + 1 
      Solver % Variable % Norm = 0.0_dp
      Solver % Variable % NonlinConverged = 1

      RETURN
    END IF

    AndersonAcc = ListGetLogical( Params,'Nonlinear System Acceleration',GotIt ) 
    AndersonScaled = ListgetLogical( Params,'Nonlinear System Acceleration Scaled',GotIt ) 
    
    IF( AndersonAcc .AND. .NOT. AndersonScaled ) THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .TRUE., NoSolve )
      IF(NoSolve) GOTO 120
    END IF
    
!   Convert rhs & initial value to the scaled system:
!   -------------------------------------------------
    IF ( ScaleSystem ) THEN
      ApplyRowEquilibration = ListGetLogical(Params,'Linear System Row Equilibration',GotIt)
      IF ( ApplyRowEquilibration ) THEN
        Parallel = ParEnv % PEs > 1
        CALL RowEquilibration(A, b, Parallel)
      ELSE
        CALL ScaleLinearSystem(Solver, A, b, x, &
            RhsScaling = (bnorm/=0._dp), ConstraintScaling=.TRUE. )
      END IF
    END IF

    ComputeChangeScaled = ListGetLogical(Params,&
        'Nonlinear System Compute Change in Scaled System',GotIt)
    IF(.NOT.GotIt) ComputeChangeScaled = .FALSE.

    IF(ComputeChangeScaled) THEN
       ALLOCATE(NonlinVals(SIZE(x)))
       NonlinVals = x
       IF (ASSOCIATED(Solver % Variable % Perm)) & 
           CALL RotateNTSystemAll(NonlinVals, Solver % Variable % Perm, DOFs)
    END IF

    IF( AndersonAcc .AND. AndersonScaled ) THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .TRUE., NoSolve )
      IF( NoSolve ) GOTO 110
    END IF
    
    ! Sometimes the r.h.s. may abruptly diminish in value resulting to significant 
    ! convergence issues or it may be that the system scales linearly with the source. 
    ! This flag tries to improve on the initial guess of the linear solvers, and may 
    ! sometimes even result to the exact solution.
    IF( ListGetLogical( Params,'Linear System Normalize Guess',GotIt ) ) THEN
      ALLOCATE( TempVector(A % NumberOfRows) )

      IF ( ParEnv % PEs > 1 ) THEN
        IF( .NOT. ALLOCATED( TempRHS ) ) THEN
          ALLOCATE( TempRHS(A % NumberOfRows) ); TempRHS=0._dp
        END IF

        Tempvector = 0._dp
        TempRHS(1:n) = b(1:n)
        CALL ParallelInitSolve( A, x, TempRHS, Tempvector )

        MP => ParallelMatrix(A,mx,mb,mr)
        mn = MP % NumberOfRows

        TempVector = 0._dp
        CALL ParallelMatrixVector( A, mx, TempVector )

        bn = ParallelDot( mn, TempVector, mb )
        xn = ParallelDot( mn, TempVector, TempVector )
        DEALLOCATE( TempRHS )
      ELSE
        CALL MatrixVectorMultiply( A, x, TempVector )
        xn = SUM( TempVector(1:n)**2 )
        bn = SUM( TempVector(1:n) * b(1:n) )
      END IF

      IF( xn > TINY( xn ) ) THEN
        x(1:n) = x(1:n) * ( bn / xn )
        WRITE( Message,'(A,ES12.3)') 'Linear System Normalizing Factor: ',bn/xn
        CALL Info('SolveLinearSystem',Message,Level=8) 
      END IF
      DEALLOCATE( TempVector )
    END IF

    IF( ListGetLogical( Params,'Linear System Nullify Guess',GotIt ) ) THEN
      x(1:n) = 0.0_dp
    END IF

    Method = ListGetString(Params,'Linear System Solver',GotIt)
    CALL Info('SolveSystem','Linear System Solver: '//TRIM(Method),Level=8)

    IF (Method=='multigrid' .OR. Method=='iterative' ) THEN
      Prec = ListGetString(Params,'Linear System Preconditioning',GotIt)
      IF( GotIt ) THEN
        CALL Info('SolveSystem','Linear System Preconditioning: '//TRIM(Prec),Level=8)
        IF ( Prec=='vanka' ) CALL VankaCreate(A,Solver)
        IF ( Prec=='circuit' ) CALL CircuitPrecCreate(A,Solver)
      END IF
    END IF

    IF ( ParEnv % PEs <= 1 ) THEN
      CALL Info('SolveSystem','Serial linear System Solver: '//TRIM(Method),Level=8)
      
      SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL IterSolver( A, x, b, Solver )
      CASE('feti')
        CALL Fatal('SolveLinearSystem', &
            'Feti solver available only in parallel.')
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
      CASE DEFAULT
        CALL DirectSolver( A, x, b, Solver )        
      END SELECT
    ELSE
      CALL Info('SolveSystem','Parallel linear System Solver: '//TRIM(Method),Level=8)

    SELECT CASE(Method)
      CASE('multigrid')
        CALL MultiGridSolve( A, x, b, &
            DOFs, Solver, Solver % MultiGridLevel )
      CASE('iterative')
        CALL ParallelIter( A, A % ParallelInfo, DOFs, &
            x, b, Solver, A % ParMatrix )
      CASE('feti')
        CALL FetiSolver( A, x, b, Solver )
      CASE('block')
        CALL BlockSolveExt( A, x, b, Solver )
     CASE DEFAULT
        CALL DirectSolver( A, x, b, Solver )
      END SELECT
    END IF

110 IF( AndersonAcc .AND. AndersonScaled )  THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .FALSE.)
    END IF
    
    IF(ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, NonlinVals, Matrix=A, RHS=b )
      DEALLOCATE(NonlinVals)
    END IF

    IF ( ScaleSystem ) THEN
      IF ( ApplyRowEquilibration ) THEN
        CALL ReverseRowEquilibration( A, b )
      ELSE
        CALL BackScaleLinearSystem( Solver, A, b, x, ConstraintScaling=.TRUE. )
      END IF
    END IF

120 IF( AndersonAcc .AND. .NOT. AndersonScaled )  THEN
      CALL NonlinearAcceleration( A, x, b, Solver, .FALSE.)
    END IF
    
    Aaid => A
    IF (PRESENT(BulkMatrix)) THEN
      IF (ASSOCIATED(BulkMatrix) ) Aaid=>BulkMatrix
    END IF
    
    NodalLoads => VariableGet( Solver % Mesh % Variables, &
        GetVarName(Solver % Variable) // ' Loads' )
    IF( ASSOCIATED( NodalLoads ) ) THEN
      ! Nodal loads may be allocated but the user may have toggled
      ! the 'calculate loads' flag such that no load computation should be performed.
      CalcLoads = ListGetLogical( Solver % Values,'Calculate Loads',GotIt )
      IF( .NOT. GotIt ) CalcLoads = .TRUE.
      IF( CalcLoads ) THEN
        CALL Info('SolveLinearSystem','Calculating nodal loads',Level=6)
        CALL CalculateLoads( Solver, Aaid, x, Dofs, .TRUE., NodalLoads ) 
      END IF
    END IF

    IF (BackRotation) THEN
      CALL BackRotateNTSystem( x, Solver % Variable % Perm, DOFs )
      IF( ASSOCIATED( NodalLoads ) ) THEN
        CALL BackRotateNTSystem(NodalLoads % Values,NodalLoads % Perm,DOFs)
      END IF
    END IF

!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
! Compute the change of the solution with different methods 
!------------------------------------------------------------------------------
    IF(.NOT.ComputeChangeScaled) THEN
      CALL ComputeChange(Solver,.FALSE.,n, x, Matrix=A, RHS=b )
    END IF
    Norm = Solver % Variable % Norm

!------------------------------------------------------------------------------
 
   Solver % Variable % PrimaryMesh => Solver % Mesh
   CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
         GetVarName(Solver % Variable) )
   
   IF ( ASSOCIATED( NodalLoads ) ) THEN
     NodalLoads % PrimaryMesh => Solver % Mesh
     CALL InvalidateVariable( CurrentModel % Meshes, Solver % Mesh, &
                  GetVarName(NodalLoads) )
   END IF

!------------------------------------------------------------------------------
! In order to be able to change the preconditioners or solvers the old matrix structures
! must be deallocated on request.

    IF( ListGetLogical( Params, 'Linear System Preconditioning Deallocate', GotIt) ) THEN
       ! ILU preconditioning
       IF( ASSOCIATED(A % ILUValues) ) THEN
          IF(  SIZE( A % ILUValues) /= SIZE(A % Values) ) &
             DEALLOCATE(A % ILUCols, A % ILURows, A % ILUDiag)
          DEALLOCATE(A % ILUValues)
       END IF
          
       ! Multigrid solver / preconditioner
       IF( Solver % MultigridLevel > 0 ) THEN
          Aaid => A 
          IF(ASSOCIATED( Aaid % Parent) ) THEN
             DO WHILE( ASSOCIATED( Aaid % Parent ) )
                Aaid => Aaid % Parent
             END DO
             DO WHILE( ASSOCIATED( Aaid % Child) )
                Aaid => Aaid % Child
                IF(ASSOCIATED(Aaid % Parent)) DEALLOCATE(Aaid % Parent )
                IF(ASSOCIATED(Aaid % Ematrix)) DEALLOCATE(Aaid % Ematrix )
             END DO
          END IF
       END IF
    END IF

  END SUBROUTINE SolveLinearSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  SUBROUTINE LinearSystemResidual( A, b, x, r )

    REAL(KIND=dp) CONTIG :: b(:)   
    REAL(KIND=dp) CONTIG :: x(:)   
    TYPE(Matrix_t), POINTER :: A   
    REAL(KIND=dp), POINTER :: r(:)
    REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: TmpXVec, TmpRVec, TmpRHSVec

    INTEGER :: i,n,nn 

    n = A % NumberOfRows 

    IF (Parenv % Pes>1) THEN
      CALL ParallelInitSolve(A,x,b,r)
      CALL ParallelMatrixVector(A,x,r,.TRUE.)
    ELSE
      CALL MatrixVectorMultiply( A, x, r)
    END IF

    DO i=1,n
      r(i) = b(i) - r(i)
    END DO

  END SUBROUTINE LinearSystemResidual



!------------------------------------------------------------------------------
!> Given a linear system Ax=b make a change of variables such that we will 
!> be solving for the residual Adx=b-Ax0 where dx=x-x0.
!------------------------------------------------------------------------------
  FUNCTION LinearSystemMaskedResidualNorm( A, b, x, ActiveRow, ActiveCol ) RESULT ( Nrm )

    REAL(KIND=dp) CONTIG :: b(:)   
    REAL(KIND=dp) CONTIG :: x(:)   
    TYPE(Matrix_t), POINTER :: A
    LOGICAL, DIMENSION(:) :: ActiveRow(:), ActiveCol(:)
    REAL(KIND=dp) :: Nrm
    
    REAL(KIND=dp), ALLOCATABLE :: r(:)
    INTEGER :: i,n,totn
    REAL(KIND=dp) :: r2sum

    n = A % NumberOfRows 

    ALLOCATE(r(n))
   
    IF (Parenv % Pes>1) THEN
      CALL Fatal('LinearSystemMaskedResidualNorm','Not implemented in parallel yet!')
!      CALL ParallelMatrixVector(A, x, r, .TRUE.)
    ELSE
      CALL MaskedMatrixVectorMultiply( A, x, r, ActiveRow, ActiveCol )
    END IF

    DO i=1,n
      IF( ActiveRow(i) ) THEN
        r(i) = b(i) - r(i)
      END IF
    END DO

    totn = NINT( ParallelReduction(1.0_dp * n) )

    r2sum = SUM( r**2 )
    Nrm = SQRT( ParallelReduction(r2sum) / totn )

    DEALLOCATE( r ) 
    
  END FUNCTION LinearSystemMaskedResidualNorm



  FUNCTION HaveConstraintMatrix( A ) RESULT( HaveConstraint ) 

    TYPE(Matrix_t), POINTER :: A
    LOGICAL :: HaveConstraint

    INTEGER :: n

    IF( .NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('EnquireConstraintMatrix','Matrix A not associated!')
    END IF

    n = 0
    IF ( ASSOCIATED(A % ConstraintMatrix) )  THEN
      IF ( A % ConstraintMatrix % NumberOFRows > 0 ) n = n + 1 
    END IF
         
    IF ( ASSOCIATED(A % AddMatrix) )  THEN
      IF ( A % AddMatrix % NumberOFRows > 0 ) n = n + 1
    END IF
    
    n = ParallelReduction(n*1._dp)

    HaveConstraint = ( n > 0 ) 
    
  END FUNCTION HaveConstraintMatrix

  
  
!------------------------------------------------------------------------------
!> Solve a system. Various additional utilities are included and 
!> naturally a call to the linear system solver.
!------------------------------------------------------------------------------
  RECURSIVE SUBROUTINE SolveSystem( A,ParA,b,x,Norm,DOFs,Solver )
!------------------------------------------------------------------------------
    REAL(KIND=dp) CONTIG, TARGET :: b(:)   !< The RHS vector
    REAL(KIND=dp) CONTIG :: x(:)   !< Previous solution on entry, new solution on exit (hopefully)
    REAL(KIND=dp) :: Norm          !< L2 Norm of solution
    TYPE(Matrix_t), POINTER :: A   !< The coefficient matrix
    INTEGER :: DOFs                !< Number of degrees of freedom per node for this equation
    TYPE(Solver_t), TARGET :: Solver                 !< Holds various solver options.
    TYPE(SParIterSolverGlobalD_t), POINTER :: ParA   !< holds info for parallel solver, 
                                                     !< if not executing in parallel this is just a dummy.
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var, NodalLoads
    TYPE(Mesh_t), POINTER :: Mesh, SaveMEsh
    LOGICAL :: Relax, Found, NeedPrevSol, Timing, ResidualMode,ConstraintMode, BlockMode, GloNum
    INTEGER :: n,i,j,k,l,m,istat,nrows,ncols,colsj,rowoffset
    CHARACTER(LEN=MAX_NAME_LEN) :: Method, ProcName, VariableName
    INTEGER(KIND=AddrInt) :: Proc
    REAL(KIND=dp) :: Relaxation,Beta,Gamma
    REAL(KIND=dp), ALLOCATABLE :: Diag(:), TempVector(:)
    REAL(KIND=dp), POINTER :: bb(:),Res(:)
    REAL(KIND=dp) :: t0,rt0,rst,st,ct
    TYPE(ValueList_t), POINTER :: Params

    INTERFACE
      SUBROUTINE BlockSolveExt(A,x,b,Solver)
        USE Types
        TYPE(Matrix_t), POINTER :: A
        TYPE(Solver_t) :: Solver
        REAL(KIND=dp) :: x(:), b(:)
      END SUBROUTINE BlockSolveExt
    END INTERFACE   


!------------------------------------------------------------------------------
    Params => Solver % Values

    CALL Info('SolveSystem','Solving linear system',Level=10)

    Timing = ListCheckPrefix(Params,'Linear System Timing')
    IF( Timing ) THEN
      t0 = CPUTime(); rt0 = RealTime()
    END IF

    n = A % NumberOfRows

    ResidualMode = ListGetLogical( Params,'Linear System Residual Mode',Found )
    
    BlockMode = ListGetLogical( Params,'Linear System Block Mode',Found ) 
      
!------------------------------------------------------------------------------
! The allocation of previous values has to be here in order to 
! work properly with the Dirichlet elimination.
!------------------------------------------------------------------------------
    NeedPrevSol = ResidualMode

    IF(.NOT. NeedPrevSol ) THEN
      Relaxation = ListGetConstReal( Params, &
          'Nonlinear System Relaxation Factor', Found )
      IF( Found ) NeedPrevSol = (Relaxation /= 1.0_dp)
    END IF

    IF(.NOT. NeedPrevSol ) THEN
      Method = ListGetString( Params, &
        'Nonlinear System Convergence Measure', Found ) 
      NeedPrevSol = ( Method == 'residual' .OR. Method == 'solution' )
    END IF

    IF( NeedPrevSol ) THEN
      CALL Info('SolveSystem','Previous solution must be stored before system is solved',Level=10)
      Found = ASSOCIATED(Solver % Variable % NonlinValues)
      IF( Found ) THEN
        IF ( SIZE(Solver % Variable % NonlinValues) /= n) THEN
          DEALLOCATE(Solver % Variable % NonlinValues)
          Found = .FALSE.
        END IF
      END IF
      IF(.NOT. Found) THEN
        ALLOCATE( Solver % Variable % NonlinValues(n), STAT=istat ) 
        IF ( istat /= 0 ) CALL Fatal( 'SolveSystem', 'Memory allocation error.' )
      END IF
      Solver % Variable % NonlinValues = x(1:n)
    END IF

    IF ( Solver % LinBeforeProc /= 0 ) THEN
      CALL Info('SolveSystem','Calling procedure before solving system',Level=7)
      istat = ExecLinSolveProcs( Solver % LinBeforeProc,CurrentModel,Solver, &
                       A, b, x, n, DOFs, Norm )
       IF ( istat /= 0 ) GOTO 10
    END IF

    ! If residual mode is requested make change of variables:
    ! Ax=b -> Adx = b-Ax0 = r
    IF( ResidualMode ) THEN
      CALL Info('SolveSystem','Changing the equation to residual based mode',Level=10)
      ALLOCATE( Res(n) ) 

      ! If needed move the current solution to N-T coordinate system
      ! before computing the residual.
      IF (ASSOCIATED(Solver % Variable % Perm)) &
          CALL RotateNTSystemAll(x, Solver % Variable % Perm, DOFs)

      CALL LinearSystemResidual( A, b, x, res )
      bb => res
      ! Set the initial guess for the redidual system to zero
      x = 0.0_dp
    ELSE
      bb => b
    END IF

    ConstraintMode = HaveConstraintMatrix( A ) 

    ! Here activate constraint solve only if constraints are not treated as blocks
    IF( BlockMode .AND. ConstraintMode ) THEN
      CALL Warn('SolveSystem','Matrix is constraint and block matrix, giving precedence to block nature!')
    END IF
      
    IF( BlockMode ) THEN
      !ScaleSystem = ListGetLogical( Params,'Linear System Scaling', Found )
      !IF(.NOT. Found ) ScaleSystem = .TRUE.
      !IF ( ScaleSystem ) CALL ScaleLinearSystem(Solver, A )
      CALL BlockSolveExt( A, x, bb, Solver )
      !IF ( ScaleSystem ) CALL BackScaleLinearSystem( Solver, A )
  
    ELSE IF ( ConstraintMode ) THEN
      CALL Info('SolveSystem','Solving linear system with constraint matrix',Level=10)
      IF( ListGetLogical( Params,'Save Constraint Matrix',Found ) ) THEN
        GloNum = ListGetLogical( Params,'Save Constaint Matrix Global Numbering',Found )
        CALL SaveProjector(A % ConstraintMatrix,.TRUE.,'cm',Parallel=GloNum)
      END IF
      CALL SolveWithLinearRestriction( A,bb,x,Norm,DOFs,Solver )
    ELSE ! standard mode
      CALL Info('SolveSystem','Solving linear system without constraint matrix',Level=12)
      CALL SolveLinearSystem( A,bb,x,Norm,DOFs,Solver )
    END IF
    CALL Info('SolveSystem','System solved',Level=12)

    ! Even in the residual mode the system is reverted back to complete vectors 
    ! and we may forget about the residual.
    IF( ResidualMode ) DEALLOCATE( Res ) 

!------------------------------------------------------------------------------

10  CONTINUE

    IF ( Solver % LinAfterProc /= 0 ) THEN
      CALL Info('SolveSystem','Calling procedure after solving system',Level=7)
      istat = ExecLinSolveProcs( Solver % LinAfterProc, CurrentModel, Solver, &
              A, b, x, n, DOFs, Norm )
    END IF

    IF ( Solver % TimeOrder == 2 ) THEN
      CALL Info('SolveSystem','Setting up PrevValues for 2nd order transient equations',Level=12)

      IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        Gamma =  0.5d0 - Solver % Alpha
        Beta  = (1.0d0 - Solver % Alpha)**2 / 4.0d0
        DO i=1,n
          Solver % Variable % PrevValues(i,2) = &
             (1.0d0/(Beta*Solver % dt**2))* &
               (x(i)-Solver % Variable % PrevValues(i,3)) -  &
                  (1.0d0/(Beta*Solver % dt))*Solver % Variable % PrevValues(i,4)+ &
                        (1.0d0-1.0d0/(2*Beta))*Solver % Variable % PrevValues(i,5)

          Solver % Variable % PrevValues(i,1) = &
            Solver % Variable % PrevValues(i,4) + &
               Solver % dt*((1.0d0-Gamma)*Solver % Variable % PrevValues(i,5)+&
                  Gamma*Solver % Variable % PrevValues(i,2))
        END DO
      END IF
    END IF

    IF( Timing ) THEN
      st  = CPUTime() - t0;
      rst = RealTime() - rt0

      WRITE(Message,'(a,f8.2,f8.2,a)') 'Linear system time (CPU,REAL) for '&
          //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
      CALL Info('SolveSystem',Message)    
      
      IF( ListGetLogical(Params,'Linear System Timing',Found)) THEN
        CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF
      
      IF( ListGetLogical(Params,'Linear System Timing Cumulative',Found)) THEN
        ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),Found)
        st = st + ct
        ct = ListGetConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),Found)
        rst = rst + ct
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys cpu time '&
            //GetVarName(Solver % Variable),st)
        CALL ListAddConstReal(CurrentModel % Simulation,'res: cum linsys real time '&
            //GetVarName(Solver % Variable),rst)
      END IF

    END IF

    CALL Info('SolveSystem','Finished solving the system',Level=12)

!------------------------------------------------------------------------------
END SUBROUTINE SolveSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Solve a linear eigen system.
!------------------------------------------------------------------------------
SUBROUTINE SolveEigenSystem( StiffMatrix, NOFEigen, &
        EigenValues, EigenVectors,Solver )
!------------------------------------------------------------------------------
    USE EigenSolve
!------------------------------------------------------------------------------
    COMPLEX(KIND=dp) :: EigenValues(:),EigenVectors(:,:)
    REAL(KIND=dp) :: Norm
    TYPE(Matrix_t), POINTER :: StiffMatrix
    INTEGER :: NOFEigen
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER :: n
    !------------------------------------------------------------------------------
    n = StiffMatrix % NumberOfRows

    IF ( .NOT. Solver % Matrix % COMPLEX ) THEN
      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolve( Solver, StiffMatrix, n, NOFEigen, &
                EigenValues, EigenVectors )
      ELSE
        CALL ParallelArpackEigenSolve( Solver, StiffMatrix, n, NOFEigen, &
                EigenValues, EigenVectors )
      END IF
    ELSE
      IF ( ParEnv % PEs <= 1 ) THEN
        CALL ArpackEigenSolveComplex( Solver, StiffMatrix, n/2, &
              NOFEigen, EigenValues, EigenVectors )
      ELSE
        CALL ParallelArpackEigenSolveComplex( Solver, StiffMatrix, n/2, NOFEigen, &
                EigenValues, EigenVectors )
      END IF
    END IF
   
    
!------------------------------------------------------------------------------
END SUBROUTINE SolveEigenSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Solve a linear system with permutated constraints.
!------------------------------------------------------------------------------
SUBROUTINE SolveConstraintModesSystem( A, x, b, Solver )
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: A
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) CONTIG :: x(:),b(:)
!------------------------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var
    INTEGER :: i,j,k,n,m
    LOGICAL :: PrecRecompute, Stat, Found, ComputeFluxes, Symmetric
    REAL(KIND=dp), POINTER CONTIG :: PValues(:)
    REAL(KIND=dp), ALLOCATABLE :: Fluxes(:), FluxesMatrix(:,:)
    CHARACTER(LEN=MAX_NAME_LEN) :: MatrixFile
    !------------------------------------------------------------------------------
    n = A % NumberOfRows
    
    Var => Solver % Variable
    IF( SIZE(x) /= n ) THEN
      CALL Fatal('SolveConstraintModesSystem','Conflicting sizes for matrix and variable!')
    END IF

    m = Var % NumberOfConstraintModes
    IF( m == 0 ) THEN
      CALL Fatal('SolveConstraintModesSystem','No constraint modes?!')
    END IF

    ComputeFluxes = ListGetLogical( Solver % Values,'Constraint Modes Fluxes',Found) 
    IF( ComputeFluxes ) THEN
      CALL Info('SolveConstraintModesSystem','Allocating for lumped fluxes',Level=10)
      ALLOCATE( Fluxes( n ) )
      ALLOCATE( FluxesMatrix( m, m ) )
      FluxesMatrix = 0.0_dp
    END IF
    

    DO i=1,m
      CALL Info('SolveConstraintModesSystem','Solving for mode: '//TRIM(I2S(i)))

      IF( i == 2 ) THEN
        CALL ListAddLogical( Solver % Values,'No Precondition Recompute',.TRUE.)
      END IF

      ! The matrix has been manipulated already before. This ensures
      ! that the system has values 1 at the constraint mode i.
      WHERE( Var % ConstraintModesIndeces == i ) b = A % Values(A % Diag)
        
      CALL SolveSystem( A,ParMatrix,b,x,Var % Norm,Var % DOFs,Solver )

      WHERE( Var % ConstraintModesIndeces == i ) b = 0._dp

      Var % ConstraintModes(i,:) = x

      IF( ComputeFluxes ) THEN
        CALL Info('SolveConstraintModesSystem','Computing lumped fluxes',Level=8)
        PValues => A % Values
        A % Values => A % BulkValues
        Fluxes = 0.0_dp
        CALL MatrixVectorMultiply( A, x, Fluxes ) 
        A % Values => PValues

        DO j=1,n
          k = Var % ConstraintModesIndeces(j)
          IF( k > 0 ) THEN
            IF( i /= k ) THEN
              FluxesMatrix(i,k) = FluxesMatrix(i,k) - Fluxes(j)
            END IF
            FluxesMatrix(i,i) = FluxesMatrix(i,i) + Fluxes(j)
          END IF
        END DO
      END IF
    END DO

    
    IF( ComputeFluxes ) THEN
      Symmetric = ListGetLogical( Solver % Values,&
          'Constraint Modes Fluxes Symmetric', Found ) 
      IF( Symmetric ) THEN
        FluxesMatrix = 0.5_dp * ( FluxesMatrix + TRANSPOSE( FluxesMatrix ) )
      END IF
      
      CALL Info( 'SolveConstraintModesSystem','Constraint Modes Fluxes', Level=5 )
      DO i=1,m
        DO j=1,m
          IF( Symmetric .AND. j < i ) CYCLE
          WRITE( Message, '(I3,I3,ES15.5)' ) i,j,FluxesMatrix(i,j)
          CALL Info( 'SolveConstraintModesSystem', Message, Level=5 )
        END DO
      END DO
      
      MatrixFile = ListGetString(Solver % Values,'Constraint Modes Fluxes Filename',Found )
      IF( Found ) THEN
        OPEN (10, FILE=MatrixFile)
        DO i=1,m
          DO j=1,m
            WRITE (10,'(ES17.9)',advance='no') FluxesMatrix(i,j)
          END DO
          WRITE(10,'(A)') ' '
        END DO
        CLOSE(10)     
        CALL Info( 'SolveConstraintModesSystem',&
            'Constraint modes fluxes was saved to file '//TRIM(MatrixFile),Level=5)
      END IF
      
      DEALLOCATE( Fluxes )
    END IF

    CALL ListAddLogical( Solver % Values,'No Precondition Recompute',.FALSE.)
    
!------------------------------------------------------------------------------
  END SUBROUTINE SolveConstraintModesSystem
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> A parser of the variable name that returns the true variablename
!> where the inline options have been interpreted.
!------------------------------------------------------------------------------
SUBROUTINE VariableNameParser(var_name, NoOutput, Global, Dofs, IpVariable, ElemVariable, DgVariable )

  CHARACTER(LEN=*)  :: var_name
  LOGICAL, OPTIONAL :: NoOutput, Global
  INTEGER, OPTIONAL :: Dofs
  LOGICAL, OPTIONAL :: IpVariable
  LOGICAL, OPTIONAL :: ElemVariable
  LOGICAL, OPTIONAL :: DgVariable
  
  INTEGER :: i,j,k,m

  IF(PRESENT(NoOutput)) NoOutput = .FALSE.
  IF(PRESENT(Global)) Global = .FALSE.
  IF(PRESENT(Dofs)) Dofs = 0
  IF(PRESENT(IpVariable)) IpVariable = .FALSE.
  
  DO WHILE( var_name(1:1) == '-' )

    m = 0
    IF ( SEQL(var_name, '-nooutput ') ) THEN
      IF(PRESENT(NoOutput)) NoOutput = .TRUE.
      m = 10

    ELSE IF ( SEQL(var_name, '-global ') ) THEN
      IF(PRESENT(Global)) Global = .TRUE.
      m = 8

    ELSE IF ( SEQL(var_name, '-ip ') ) THEN
      IF(PRESENT(IpVariable)) IpVariable = .TRUE.      
      m = 4

    ELSE IF ( SEQL(var_name, '-dg ') ) THEN
      IF(PRESENT(DgVariable)) DgVariable = .TRUE.      
      m = 4

    ELSE IF ( SEQL(var_name, '-elem ') ) THEN
      IF(PRESENT(ElemVariable)) ElemVariable = .TRUE.      
      m = 6
    END IF
    
    IF( m > 0 ) THEN
      var_name(1:LEN(var_name)-m) = var_name(m+1:)
    END IF
   
    IF ( SEQL(var_name, '-dofs ') ) THEN
      IF(PRESENT(DOFs)) READ( var_name(7:), * ) DOFs     
      j = LEN_TRIM( var_name )
      k = 7
      DO WHILE( var_name(k:k) /= ' '  )
        k = k + 1
        IF ( k > j ) EXIT
      END DO
      var_name(1:LEN(var_name)-(k+2)) = var_name(k+1:)
    END IF
  END DO

END SUBROUTINE VariableNameParser


   !> Create permutation for fields on integration points, optionally with mask.
   !> The non-masked version is saved to Solver structure for reuse while the
   !> masked version may be unique to every variable. 
   !-----------------------------------------------------------------------------------
   SUBROUTINE CreateIpPerm( Solver, MaskPerm, MaskName, SecName, UpdateOnly )

     TYPE(Solver_t), POINTER :: Solver
     INTEGER, POINTER, OPTIONAL :: MaskPerm(:)
     CHARACTER(LEN=MAX_NAME_LEN), OPTIONAL :: MaskName, SecName      
     LOGICAL, OPTIONAL :: UpdateOnly
     
     TYPE(Mesh_t), POINTER :: Mesh
     TYPE(GaussIntegrationPoints_t) :: IP
     TYPE(Element_t), POINTER :: Element
     INTEGER :: t, n, IpCount , RelOrder, nIp
     CHARACTER(LEN=MAX_NAME_LEN) :: EquationName
     LOGICAL :: Found, ActiveElem, ActiveElem2
     INTEGER, POINTER :: IpOffset(:) 
     TYPE(ValueList_t), POINTER :: BF
     LOGICAL :: UpdatePerm
     
     n = 0
     IF( PRESENT( MaskPerm ) ) n = n + 1
     IF( PRESENT( MaskName ) ) n = n + 1
     IF( PRESENT( SecName ) ) n = n + 1
     IF( PRESENT( UpdateOnly ) ) n = n + 1
     
     ! Currently a lazy check
     IF( n /= 0 .AND. n /= 3 .AND. n /= 2) THEN
       CALL Fatal('CreateIpPerm','Only some optional parameter combinations are possible')
     END IF

     UpdatePerm = .FALSE.
     IF( PRESENT( UpdateOnly ) ) UpdatePerm = UpdateOnly

     IF( UpdatePerm ) THEN
       CALL Info('CreateIpPerm','Updating IP permutation table',Level=8)       
     ELSE IF( PRESENT( MaskPerm ) ) THEN
       CALL Info('CreateIpPerm','Creating masked permutation for integration points',Level=8)
     ELSE       
       IF( ASSOCIATED( Solver % IpTable ) ) THEN
         CALL Info('CreateIpPerm','IpTable already allocated, returning')
       END IF
       CALL Info('CreateIpPerm','Creating permutation for integration points',Level=8)
     END IF

     EquationName = ListGetString( Solver % Values, 'Equation', Found)
     IF( .NOT. Found ) THEN
       CALL Fatal('CreateIpPerm','Equation not present!')
     END IF     
     
     Mesh => Solver % Mesh
     NULLIFY( IpOffset ) 

     n = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements

     IF( UpdatePerm ) THEN
       IpOffset => MaskPerm
       ActiveElem = (IpOffset(2)-IpOffset(1) > 0 )
       IF( n >= 2 ) ActiveElem2 = (IpOffset(3)-IpOffset(2) > 0 )
     ELSE
       ALLOCATE( IpOffset( n + 1) )     
       IpOffset = 0
       IF( PRESENT( MaskPerm ) ) MaskPerm => IpOffset
     END IF
     IpCount = 0

     nIp = ListGetInteger( Solver % Values,'Gauss Points on Ip Variables', Found ) 
     
     DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
       Element => Mesh % Elements(t)
            
       IF( .NOT. UpdatePerm ) THEN
         ActiveElem = .FALSE.
         IF( Element % PartIndex == ParEnv % myPE ) THEN
           IF ( CheckElementEquation( CurrentModel, Element, EquationName ) ) THEN             
             IF( PRESENT( MaskName ) ) THEN
               BF => ListGetSection( Element, SecName )
               ActiveElem = ListGetLogicalGen( BF, MaskName )
             ELSE
               ActiveElem = .TRUE.
             END IF
           END IF
         END IF
       END IF
         
       IF( ActiveElem ) THEN
         IF( nIp > 0 ) THEN
           IpCount = IpCount + nIp
         ELSE
           IP = GaussPointsAdapt( Element )
           IpCount = IpCount + Ip % n
         END IF
       END IF

       ! We are reusing the permutation table hence we must be one step ahead 
       IF( UpdatePerm .AND. n >= t+1) THEN
         ActiveElem = ActiveElem2
         ActiveElem2 = (IpOffset(t+2)-IpOffset(t+1) > 0 )
       END IF
         
       IpOffset(t+1) = IpCount
     END DO

     IF( .NOT. PRESENT( MaskPerm ) ) THEN
       ALLOCATE( Solver % IpTable ) 
       Solver % IpTable % IpOffset => IpOffset
       Solver % IpTable % IpCount = IpCount
     END IF

     IF( UpdatePerm ) THEN
       CALL Info('CreateIpPerm','Updated permutation for IP points: '//TRIM(I2S(IpCount)),Level=8)  
     ELSE       
       CALL Info('CreateIpPerm','Created permutation for IP points: '//TRIM(I2S(IpCount)),Level=8)  
     END IF
       
   END SUBROUTINE CreateIpPerm

   
   SUBROUTINE UpdateIpPerm( Solver, Perm )

     TYPE(Solver_t), POINTER :: Solver
     INTEGER, POINTER :: Perm(:)

     CALL CreateIpPerm( Solver, Perm, UpdateOnly = .TRUE.)

   END SUBROUTINE UpdateIpPerm



!------------------------------------------------------------------------------
!> Updates values for exported variables which are typically auxiliary variables derived
!> from the solution.
!------------------------------------------------------------------------------
  SUBROUTINE UpdateExportedVariables( Solver )  
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver
  
  INTEGER :: i,j,k,l,n,m,t,bf_id,dofs,nsize,i1,i2,NoGauss
  CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name,tmpname,condname
  REAL(KIND=dp), POINTER :: Values(:), Solution(:), LocalSol(:), LocalCond(:)
  INTEGER, POINTER :: Indexes(:), VarIndexes(:), Perm(:)
  LOGICAL :: Found, Conditional, GotIt, Stat, StateVariable, AllocationsDone = .FALSE.
  LOGICAL, POINTER :: ActivePart(:),ActiveCond(:)
  TYPE(Variable_t), POINTER :: ExpVariable
  TYPE(ValueList_t), POINTER :: ValueList
  TYPE(Element_t),POINTER :: Element  
  TYPE(GaussIntegrationPoints_t) :: IP
  TYPE(Nodes_t) :: Nodes
  REAL(KIND=dp), ALLOCATABLE :: Basis(:)
  REAL(KIND=dp) :: detJ
  TYPE(ValueHandle_t) :: LocalSol_h
  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(Solver_t), POINTER :: pSolver

  
  SAVE LocalSol_h

  CALL Info('UpdateExportedVariables','Updating variables, if any!',Level=20)

  AllocationsDone = .FALSE.
  Mesh => Solver % Mesh
  
  l = 0
  DO WHILE( .TRUE. )
    l = l + 1

    str = ComponentName( 'exported variable', l )    

    var_name = ListGetString( Solver % Values, str, GotIt )    
    IF(.NOT. GotIt) EXIT

    CALL Info('UpdateExportedVariables','Trying to set values for variable: '//TRIM(Var_name),Level=20)
    
    CALL VariableNameParser( var_name ) 
    
    ExpVariable => VariableGet( Mesh % Variables, Var_name )
    IF( .NOT. ASSOCIATED(ExpVariable)) CYCLE
      
    CALL Info('UpdateExportedVariables','Setting values for variable: '//TRIM(Var_name),Level=20)
      
    IF(.NOT. AllocationsDone ) THEN
      m = CurrentModel % NumberOFBodyForces
      ALLOCATE( ActivePart(m), ActiveCond(m) )
      
      m = Mesh % MaxElementDOFs
      ALLOCATE( LocalSol(m), LocalCond(m))
      
      m =  CurrentModel % MaxElementNodes
      ALLOCATE( Basis(m), Nodes % x(m), Nodes % y(m), Nodes % z(m) )

      AllocationsDone = .TRUE.
    END IF

    Dofs = ExpVariable % DOFs
    Values => ExpVariable % Values
    Perm => ExpVariable % Perm
    n = LEN_TRIM( var_name )
    
    StateVariable = ( SIZE( Values ) == DOFs ) .OR. ( ExpVariable % Type == Variable_Global ) 
    IF( StateVariable ) THEN
      CALL Info('UpdateExportedVariables','Updating state variable',Level=20)
      IF( Dofs > 1 ) THEN
        tmpname = ComponentName( var_name(1:n), j )
        Solution => Values( j:j )
      ELSE
        tmpname = var_name(1:n)
        Solution => Values
      END IF
 
      DO bf_id=1,CurrentModel % NumberOFBodyForces
        IF( ListCheckPresent( &
            CurrentModel % BodyForces(bf_id) % Values,TmpName ) ) THEN
          CALL Info('UpdateExportedVariables',&
              'Found a proper definition for state variable',Level=6)
          Solution = ListGetCReal( CurrentModel % BodyForces(bf_id) % Values,TmpName)
          EXIT
        END IF
      END DO
      CYCLE
    END IF

    CALL Info('UpdateExportedVariables','Updating field variable with dofs: '//TRIM(I2S(DOFs)),Level=12)

        
    DO j=1,DOFs
      
100   Values => ExpVariable % Values
      IF( Dofs > 1 ) THEN
        tmpname = ComponentName( var_name(1:n), j )
        Solution => Values( j:: DOFs ) 
      ELSE
        tmpname = var_name(1:n)
        Solution => Values
      END IF
      condname = TRIM(tmpname) //' Condition' 
      
      !------------------------------------------------------------------------------
      ! Go through the Dirichlet conditions in the body force lists
      !------------------------------------------------------------------------------      
      ActivePart = .FALSE.
      ActiveCond = .FALSE.

      DO bf_id=1,CurrentModel % NumberOFBodyForces
        ActivePart(bf_id) = ListCheckPresent( &
            CurrentModel % BodyForces(bf_id) % Values,TmpName ) 
        ActiveCond(bf_id) = ListCheckPresent( &
            CurrentModel % BodyForces(bf_id) % Values,CondName )      
      END DO
      
      IF ( .NOT. ANY( ActivePart ) ) CYCLE

      CALL Info('UpdateExportedVariables','Found a proper definition in body forces',Level=6)

      
      IF( ExpVariable % TYPE == Variable_on_gauss_points ) THEN 
        ! Initialize handle when doing values on Gauss points!
        CALL ListInitElementKeyword( LocalSol_h,'Body Force',TmpName )
      END IF

      DO t = 1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

        Element => Mesh % Elements(t) 
        IF( Element % BodyId <= 0 ) CYCLE
        bf_id = ListGetInteger( CurrentModel % Bodies(Element % BodyId) % Values,&
            'Body Force',GotIt)

        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id)) CYCLE
        Conditional = ActiveCond(bf_id)
        
        CurrentModel % CurrentElement => Element
        m = Element % TYPE % NumberOfNodes
        Indexes => Element % NodeIndexes
        ValueList => CurrentModel % BodyForces(bf_id) % Values

        IF( ExpVariable % TYPE == Variable_on_gauss_points ) THEN 
    
          i1 = Perm( Element % ElementIndex )
          i2 = Perm( Element % ElementIndex + 1 )
          NoGauss = i2 - i1
          
          ! This is not active here
          IF( NoGauss == 0 ) CYCLE
          
          IP = GaussPointsAdapt( Element, Solver )

          IF( NoGauss /= IP % n ) THEN
            
            CALL Info('UpdateExportedVariables','Number of Gauss points has changed, redoing permutations!',Level=8)

            pSolver => Solver
            CALL UpdateIpPerm( pSolver, Perm )
            nsize = MAXVAL( Perm )

            CALL Info('UpdateExportedVariables','Total number of new IP dofs: '//TRIM(I2S(nsize)))

            IF( SIZE( ExpVariable % Values ) /= ExpVariable % Dofs * nsize ) THEN
              DEALLOCATE( ExpVariable % Values )
              ALLOCATE( ExpVariable % Values( nsize * ExpVariable % Dofs ) )
            END IF
            ExpVariable % Values = 0.0_dp
            GOTO 100 
          END IF
          
          Nodes % x(1:m) = Mesh % Nodes % x(Indexes)
          Nodes % y(1:m) = Mesh % Nodes % y(Indexes)
          Nodes % z(1:m) = Mesh % Nodes % z(Indexes)

          IF( Conditional ) THEN
            CALL Warn('UpdateExportedVariable','Elemental variable cannot be conditional!')
          END IF

          DO k=1,IP % n
            stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
                IP % W(k), detJ, Basis )
            Solution(i1+k) = ListGetElementReal( LocalSol_h,Basis,Element,Found,GaussPoint=k) 
          END DO
          
        ELSE IF( ExpVariable % TYPE == Variable_on_elements ) THEN
          IF( Conditional ) THEN
            CALL Warn('UpdateExportedVariable','Elemental variables not conditional!')
          END IF
          LocalSol(1:m) = ListGetReal(ValueList, TmpName, m, Indexes(1:m) )
          i = Perm( Element % ElementIndex ) 
          IF( i > 0 ) Solution(i) = SUM( LocalSol(1:m) ) / m
          
        ELSE
          IF( ExpVariable % TYPE == Variable_on_nodes_on_elements ) THEN
            VarIndexes => Element % DGIndexes
          ELSE
            VarIndexes => Indexes
          END IF
          
          LocalSol(1:m) = ListGetReal(ValueList, TmpName, m, Indexes(1:m) )
          
          IF( Conditional ) THEN
            LocalCond(1:m) = ListGetReal(ValueList, CondName, m, Indexes(1:m) )
            DO i=1,m
              IF( LocalCond(i) > 0.0_dp ) THEN
                IF( Perm(VarIndexes(i)) > 0 ) THEN
                  Solution( Perm(VarIndexes(i)) ) = LocalSol(i)
                END IF
              END IF
            END DO
          ELSE
            IF( ALL( Perm(VarIndexes(1:m)) > 0 ) ) THEN
              Solution( Perm(VarIndexes(1:m)) ) = LocalSol(1:m)
            END IF
          END IF
          
        END IF
      END DO
      
    END DO
  END DO

  IF( AllocationsDone ) THEN
    DEALLOCATE(ActivePart, ActiveCond, LocalSol, LocalCond, Basis, &
        Nodes % x, Nodes % y, Nodes % z )
  END IF
    
END SUBROUTINE UpdateExportedVariables


!------------------------------------------------------------------------------
!> Derivates values for exported variables to come up with velocity and
!> acceleration fields.
!------------------------------------------------------------------------------
  SUBROUTINE DerivateExportedVariables( Solver )  
!------------------------------------------------------------------------------
  TYPE(Solver_t), TARGET :: Solver

  TYPE(Mesh_t), POINTER :: Mesh
  TYPE(ValueList_t), POINTER :: Params
  TYPE(Variable_t), POINTER :: Var, DerVar, dtVar
  CHARACTER(LEN=MAX_NAME_LEN) :: str, var_name
  INTEGER :: VarNo
  LOGICAL :: Found, DoIt
  REAL(KIND=dp) :: dt
  
  
  CALL Info('DerivateExportedVariables','Derivating variables, if any!',Level=20)

  Mesh => Solver % Mesh
  Params => Solver % Values

  VarNo = 0
  DO WHILE( .TRUE. )
    VarNo = VarNo + 1

    str = ComponentName( 'exported variable', VarNo )    
    
    var_name = ListGetString( Solver % Values, str, Found )    
    IF(.NOT. Found) EXIT
    
    CALL VariableNameParser( var_name ) 

    Var => VariableGet( Mesh % Variables, Var_name )
    IF( .NOT. ASSOCIATED(Var)) CYCLE
    IF( .NOT. ASSOCIATED(Var % PrevValues) ) CYCLE
    
    str = TRIM( ComponentName(Var_name) )//' Calculate Velocity'
    DoIt = ListGetLogical( Params, str, Found )        
    IF( DoIt ) THEN
      str = TRIM( ComponentName(var_name) ) // ' Velocity'
      DerVar => VariableGet( Solver % Mesh % Variables, str )        
      IF(.NOT. ASSOCIATED(DerVar)) THEN
        CALL Warn('DerivatingExportedVariables','Variable does not exist:'//TRIM(str))
        CYCLE
      END IF

      dtVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
      dt = dtVar % Values(1) 
      
      CALL Info('DerivatingExportedVariables','Computing numerical derivative for:'//TRIM(str),Level=8)     
      DerVar % Values = (Var % Values(:) - Var % PrevValues(:,1)) / dt
    END IF

    str = TRIM( ComponentName(Var_name) )//' Calculate Acceleration'
    DoIt = ListGetLogical( Params, str, Found )        
    IF( DoIt ) THEN
      str = TRIM( ComponentName(var_name) ) // ' Acceleration'
      DerVar => VariableGet( Solver % Mesh % Variables, str )        
      IF(.NOT. ASSOCIATED(DerVar)) THEN
        CALL Warn('DerivatingExportedVariables','Variable does not exist:'//TRIM(str))
        CYCLE
      END IF

      dtVar => VariableGet( Solver % Mesh % Variables, 'timestep size' )
      dt = dtVar % Values(1) 

      CALL Info('DerivatingExportedVariables','Computing numerical derivative for:'//TRIM(str),Level=8)     
      DerVar % Values = (Var % Values(:) - 2*Var % PrevValues(:,1) - Var % PrevValues(:,2)) / dt**2
    END IF

  END DO
    
END SUBROUTINE DerivateExportedVariables




!------------------------------------------------------------------------------
!> Solves a harmonic system.
!------------------------------------------------------------------------------
SUBROUTINE SolveHarmonicSystem( G, Solver )
!------------------------------------------------------------------------------
    USE DirichletUtils
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), TARGET :: G
!------------------------------------------------------------------------------
    TYPE(Matrix_t), POINTER :: BMatrix, A => NULL()
    INTEGER :: i,j,k,n, kr, ki, DOFs, ne, niter
    LOGICAL :: stat, Found, OptimizeBW, Real_given,Imag_given
    CHARACTER(LEN=MAX_NAME_LEN) :: Name
    REAL(KIND=dp) :: Omega, norm, s
    REAL(KIND=dp), POINTER :: Freqv(:,:)
    REAL(KIND=dp), ALLOCATABLE :: x(:)
    REAL(KIND=dp), POINTER :: b(:)
    REAL(KIND=dp) :: frequency
    INTEGER :: Nfrequency
    TYPE(ValueList_t), POINTER :: BC

    CALL Info( 'HarmonicSolve', 'Solving initially transient style system as harmonic one', Level=5)
    
    n = Solver % Matrix % NumberofRows
    DOFs = Solver % Variable % DOFs * 2

    A => G
    DO WHILE( ASSOCIATED(A) )
      BMatrix => A
      A => A % EMatrix
      IF ( ASSOCIATED(A) ) THEN
        IF ( A % COMPLEX ) THEN
          CALL Info('SolveHarmonicSystem','Reusing existing harmonic system',Level=10)
          EXIT
        END IF
      END IF
    END DO

    IF ( .NOT. ASSOCIATED(A) ) THEN      
      CALL Info('SolveHarmonicSystem','Creating new matrix for harmonic system',Level=10)      

      OptimizeBW = ListGetLogical(Solver % Values, 'Optimize Bandwidth', Found)
      IF ( .NOT. Found ) OptimizeBW = .TRUE.
      
      A => CreateMatrix( CurrentModel, Solver, Solver % Mesh,   &
              Solver % Variable % Perm, DOFs, MATRIX_CRS, OptimizeBW, &
              ListGetString( Solver % Values, 'Equation') )
      A % COMPLEX = .TRUE.
      BMatrix % EMatrix => A
      ALLOCATE( A % rhs(2*n) )
      
      DO j=1,Solver % Variable % DOFs
        Name = ComponentName( Solver % Variable % Name, j ) 
        DO i=1,CurrentModel % NumberOFBCs
          BC => CurrentModel % BCs(i) % Values
          real_given = ListCheckPresent( BC, Name )
          imag_given = ListCheckPresent( BC, TRIM(Name) // ' im' )
          
          IF ( real_given .AND. .NOT. imag_given ) THEN
            CALL ListAddConstReal( BC, TRIM(Name) // ' im', 0._dp)
          ELSE IF ( imag_given .AND. .NOT. real_given ) THEN
            CALL ListAddConstReal( BC, Name, 0._dp )
          END IF
        END DO
      END DO
    END IF

    b => A % rhs
    ALLOCATE( x(2*n) )
    x = 0
    
    b(1:2*n:2) = G % RHS(1:n)
    b(2:2*n:2) = G % RHS_im(1:n)

    
    Nfrequency = ListGetInteger( Solver % Values,'Harmonic System Values',Found )
    IF( Nfrequency > 1 ) THEN
      freqv => ListGetConstRealArray( Solver % Values, 'Frequency' )
    ELSE
      Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
      IF( .NOT. Found ) THEN
        CALL Fatal( 'AddEquation', '> Frequency < must be given for harmonic analysis.' )
      END IF
      
      Nfrequency = 1
      ! Add the number of frequencies even for case of one for some postprocessing stuff to work 
      CALL ListAddInteger( Solver % Values,'Harmonic System Values',Nfrequency )
    END IF
    
    niter = MIN(Nfrequency,Solver % NOFEigenValues)
    ne=Solver % NofEigenValues
    Solver % NofEigenValues=0

    DO i=1,niter
      IF( Nfrequency > 1 ) THEN
        Frequency = freqv(i,1)
        WRITE( Message, '(a,i5,e12.3)' ) 'Frequency sweep: ', i, frequency
      ELSE
        WRITE( Message, '(a,e12.3)' ) 'Frequency value: ', frequency
      END IF
      CALL Info( 'HarmonicSolve', Message )

      omega = 2 * PI * Frequency
      DO k=1,n
        kr = A % Rows(2*(k-1)+1)
        ki = A % Rows(2*(k-1)+2)
        DO j=G % Rows(k),G % Rows(k+1)-1
          A % Values(kr)   =  G % Values(j)
          IF (ASSOCIATED(G % MassValues)) A % Values(kr) = &
              A % Values(kr) - omega**2*G % MassValues(j)
          IF (ASSOCIATED(G % DampValues)) THEN
            A % Values(kr+1) = -G % Dampvalues(j) * omega
            A % Values(ki)   =  G % Dampvalues(j) * omega
          END IF
          A % Values(ki+1) =  G % Values(j)
          IF (ASSOCIATED(G % MassValues)) A % Values(ki+1) = &
            A % Values(ki+1) - omega**2*G % MassValues(j)
          kr = kr + 2
          ki = ki + 2
        END DO
      END DO

      
      DO j=1,Solver % Variable % DOFs
        Name = ComponentName( Solver % Variable % Name, j ) 

        CALL SetDirichletBoundaries( CurrentModel, A, b, Name, &
                2*j-1, DOFs, Solver % Variable % Perm )

        CALL SetDirichletBoundaries( CurrentModel, A, b, TRIM(Name) // ' im', &
                2*j, DOFs, Solver % Variable % Perm )
      END DO

      CALL EnforceDirichletConditions( Solver, A, b )
 
      
      CALL SolveLinearSystem( A, b, x, Norm, DOFs, Solver )
      
      DO j=1,n
        Solver % Variable % EigenVectors(i,j) = &
                 CMPLX( x(2*(j-1)+1),x(2*(j-1)+2),KIND=dp )
      END DO
    END DO

    Solver % NOFEigenValues = ne

    DEALLOCATE( x )
!------------------------------------------------------------------------------
 END SUBROUTINE SolveHarmonicSystem
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Just toggles the initial system to harmonic one and back
!------------------------------------------------------------------------------
SUBROUTINE ChangeToHarmonicSystem( Solver, BackToReal )
!------------------------------------------------------------------------------
  USE DirichletUtils
  TYPE(Solver_t) :: Solver
  LOGICAL, OPTIONAL :: BackToReal
  !------------------------------------------------------------------------------
  TYPE(Matrix_t), POINTER :: Are => NULL(), Aharm => NULL(), SaveMatrix 
  INTEGER :: i,j,k,n, kr, ki, DOFs
  LOGICAL :: stat, Found, OptimizeBW, Real_given, Imag_given
  CHARACTER(LEN=MAX_NAME_LEN) :: Name
  REAL(KIND=dp) :: Omega, s, val
  REAL(KIND=dp), POINTER :: b(:), TmpVals(:)
  REAL(KIND=dp) :: frequency
  TYPE(ValueList_t), POINTER :: BC
  TYPE(Variable_t), POINTER :: TmpVar, ReVar, HarmVar, SaveVar
  LOGICAL :: ToReal, ParseName, AnyDirichlet, Diagonal, HarmonicReal
  
  
  IF( .NOT. ASSOCIATED( Solver % Variable ) ) THEN
    CALL Warn('CgangeToHarmonicSystem','Not applicable without a variable')
    RETURN    
  END IF

  IF( .NOT. ASSOCIATED( Solver % Matrix ) ) THEN
    CALL Warn('CgangeToHarmonicSystem','Not applicable without a matrix')
    RETURN    
  END IF
  
  ToReal = .FALSE.
  IF( PRESENT( BackToReal ) ) ToReal = BackToReal

  IF( ToReal ) THEN
    IF( ASSOCIATED( Solver % Variable % Evar ) ) THEN
      IF( Solver % Variable % Evar % Dofs < Solver % Variable % Dofs ) THEN
        CALL Info('ChangeToHarmonicSystem','Changing the harmonic results back to real system!',Level=5)

        SaveVar => Solver % Variable
        SaveMatrix => Solver % Matrix 

        Solver % Variable => Solver % Variable % Evar
        Solver % Variable % Evar => SaveVar

        Solver % Matrix => Solver % Matrix % EMatrix
        Solver % Matrix % Ematrix => SaveMatrix

        ! Eliminate cyclic dependence that is a bummer when deallocating stuff
        NULLIFY( Solver % Matrix % EMatrix % Ematrix )
      END IF
    END IF
    RETURN
  END IF


  CALL Info('ChangeToHarmonicSystem','Changing the real transient system to harmonic one!',Level=5)

  SaveMatrix => Solver % Matrix
  SaveVar => Solver % Variable     

  n = Solver % Matrix % NumberofRows
  DOFs = SaveVar % Dofs
  Are => Solver % Matrix

  CALL Info('ChangeToHarmonicSystem','Number of real system rows: '//TRIM(I2S(n)),Level=16)
  
  ! Obtain the frequency, it may depend on iteration step etc. 
  Frequency = ListGetAngularFrequency( Solver % Values, Found ) / (2*PI)
  IF( .NOT. Found ) THEN
    CALL Fatal( 'ChangeToHarmonicSystem', '> Frequency < must be given for harmonic analysis.' )
  END IF
  WRITE( Message, '(a,e12.3)' ) 'Frequency value: ', frequency
  CALL Info( 'ChangeToHarmonicSystem', Message )
  omega = 2 * PI * Frequency

  
  CALL ListAddConstReal( CurrentModel % Simulation, 'res: frequency', Frequency )

  
  HarmonicReal = ListGetLogical( Solver % Values,'Harmonic Mode Real',Found ) 
  IF( HarmonicReal ) THEN
    CALL Info('ChangeToHarmonicSystem','Enforcing harmonic system to be real valued',Level=8)
    IF (ASSOCIATED(Are % MassValues)) THEN
      ARe % Values = Are % Values - omega**2* Are % MassValues
    ELSE
      CALL Fatal('ChangeToHarmonicSystem','Harmonic system requires mass!')
    END IF
    ! This is set outside so that it can be called more flexibilly
    CALL EnforceDirichletConditions( Solver, Are, Are % rhs  )
    RETURN
  END IF

 
  Diagonal = ListGetLogical( Solver % Values,'Harmonic Mode Block Diagonal',Found )  
  IF(.NOT. Found ) Diagonal = .NOT. ASSOCIATED(Are % DampValues)
  IF( Diagonal ) THEN
    CALL Info('ChangeToHarmonicSystem','Undamped system is assumed to be block diagonal')
  END IF

  
  ! Find whether the matrix already exists
  Aharm => Are % EMatrix
  IF( ASSOCIATED( Aharm ) ) THEN
    CALL Info('ChangeToHarmonicSystem','Found existing harmonic system',Level=10)
    IF( ALLOCATED( Aharm % ConstrainedDOF ) ) Aharm % ConstrainedDOF = .FALSE.
  ELSE    
    ! Create the matrix if it does not
    
    Aharm => CreateChildMatrix( Are, Dofs, 2*Dofs, CreateRhs = .TRUE., Diagonal = Diagonal )

    IF( ParEnv % PEs > 1 ) THEN
      CALL Warn('ChangeToHarmonicSystem','ParallelInfo may not have been generated properly!')
    END IF

    Aharm % COMPLEX = ListGetLogical( Solver % Values,'Linear System Complex', Found ) 
    IF( .NOT. Found ) Aharm % COMPLEX = .NOT. Diagonal !TRUE. 
  END IF


  ! Set the harmonic system r.h.s
  b => Aharm % rhs
  
  IF( ASSOCIATED( Are % Rhs ) ) THEN
    b(1:2*n:2) = Are % RHS(1:n)
  ELSE
    b(1:2*n:2) = 0.0_dp
  END IF
  
  IF( ASSOCIATED( Are % Rhs_im ) ) THEN
    b(2:2*n:2) = Are % RHS_im(1:n)            
  ELSE
    b(2:2*n:2) = 0.0_dp
  END IF

  IF( ASSOCIATED(Are % MassValues) ) THEN
    CALL Info('ChangeToHarmonicSystem','We have mass matrix values',Level=12)
  ELSE
    CALL Warn('ChangeToHarmonicSystem','We do not have mass matrix values!')
  END IF

  IF( ASSOCIATED(Are % DampValues) ) THEN
    CALL Info('ChangeToHarmonicSystem','We have damp matrix values',Level=12)
    IF( Diagonal ) THEN
      CALL Fatal('ChangeToHarmonicSystem','Damping matrix cannot be block diagonal!')
    END IF
  ELSE
    CALL Info('ChangeToHarmonicSystem','We do not have damp matrix values',Level=12)
  END IF


  ! Set the harmonic system matrix
  IF( Diagonal ) THEN
    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        val = Are % Values(j)
        IF (ASSOCIATED(Are % MassValues)) val = val - omega**2* Are % MassValues(j)
        
        Aharm % Values(kr) = val 
        Aharm % Values(ki) = val 
        kr = kr + 1
        ki = ki + 1
      END DO
    END DO
  ELSE
    DO k=1,n
      kr = Aharm % Rows(2*(k-1)+1)
      ki = Aharm % Rows(2*(k-1)+2)
      DO j=Are % Rows(k),Are % Rows(k+1)-1
        val = Are % Values(j)
        IF (ASSOCIATED(Are % MassValues)) val = val - omega**2* Are % MassValues(j)

        Aharm % Values(kr) = val
        Aharm % Values(ki+1) = val     
        
        IF (ASSOCIATED(Are % DampValues)) THEN
          Aharm % Values(kr+1) = -Are % Dampvalues(j) * omega
          Aharm % Values(ki)   =  Are % Dampvalues(j) * omega
        END IF

        kr = kr + 2
        ki = ki + 2
      END DO
    END DO
  END IF
    
  AnyDirichlet = .FALSE.
  
  ! Finally set the Dirichlet conditions for the solver    
  DO j=1,DOFs
    Name = ComponentName( Solver % Variable % Name, j ) 
    DO i=1,CurrentModel % NumberOFBCs
      BC => CurrentModel % BCs(i) % Values
      real_given = ListCheckPresent( BC, Name )
      imag_given = ListCheckPresent( BC, TRIM(Name) // ' im' )

      IF( real_given .OR. imag_given ) AnyDirichlet = .TRUE.

      IF ( real_given .AND. .NOT. imag_given ) THEN
        CALL Info('ChangeToHarmonicSystem','Setting zero >'//TRIM(Name)//' im< on BC '//TRIM(I2S(i)),Level=12)
        CALL ListAddConstReal( BC, TRIM(Name) // ' im', 0._dp)
      ELSE IF ( imag_given .AND. .NOT. real_given ) THEN
        CALL Info('ChangeToHarmonicSystem','Setting zero >'//TRIM(Name)//'< on BC '//TRIM(I2S(i)),Level=12)
        CALL ListAddConstReal( BC, Name, 0._dp )
      END IF
    END DO
  END DO



  IF( AnyDirichlet ) THEN
    DO j=1,DOFs
      Name = ComponentName( SaveVar % Name, j ) 
      
      CALL SetDirichletBoundaries( CurrentModel, Aharm, b, Name, &
          2*j-1, 2*DOFs, SaveVar % Perm )

      CALL SetDirichletBoundaries( CurrentModel, Aharm, b, TRIM(Name) // ' im', &
          2*j, 2*DOFs, SaveVar % Perm )
    END DO

    CALL EnforceDirichletConditions( Solver, Aharm, b )
  END IF


  
  ! Create the new fields, the total one and the imaginary one
  !-------------------------------------------------------------
  k = INDEX( SaveVar % name, '[' )
  ParseName = ( k > 0 ) 

  ! Name of the full complex variable not used for postprocessing
  IF( ParseName ) THEN
    Name = TRIM(SaveVar % Name(1:k-1))//' complex'
  ELSE
    Name = TRIM( SaveVar % Name )//' complex'
  END IF

  CALL Info('ChangeToHarmonicSystem','Harmonic system full name: '//TRIM(Name),Level=12)


  HarmVar => VariableGet( Solver % Mesh % Variables, Name )
  IF( ASSOCIATED( HarmVar ) ) THEN
    CALL Info('ChangeToHarmonicSystem','Reusing full system harmonic dofs',Level=12)
  ELSE
    CALL Info('ChangeToHarmonicSystem','Creating full system harmonic dofs',Level=12)
    CALL VariableAddVector( Solver % Mesh % Variables,Solver % Mesh,Solver, &
        Name,2*DOFs,Perm=SaveVar % Perm,Output=.FALSE.)
    HarmVar => VariableGet( Solver % Mesh % Variables, Name )
    IF(.NOT. ASSOCIATED( HarmVar ) ) CALL Fatal('ChangeToHarmonicSystem','New created variable should exist!')

    ! Repoint the values of the original solution vector
    HarmVar % Values(1:2*n:2) = SaveVar % Values(1:n)

    ! It beats me why this cannot be deallocated without some NaNs later
    !DEALLOCATE( SaveVar % Values )
    SaveVar % Values => HarmVar % Values(1:2*n:2)
    SaveVar % Secondary = .TRUE.

    ! Repoint the components of the original solution
    IF( Dofs > 1 ) THEN
      DO i=1,Dofs
        TmpVar => VariableGet( Solver % Mesh % Variables, ComponentName( SaveVar % Name, i ) )
        IF( ASSOCIATED( TmpVar ) ) THEN
          TmpVar % Values => HarmVar % Values(2*i-1::HarmVar % Dofs)
        ELSE
          CALL Fatal('ChangeToHarmonicSystem','Could not find re component '//TRIM(I2S(i)))
        END IF
      END DO
    END IF

    IF( ParseName ) THEN
      Name = ListGetString( Solver % Values,'Imaginary Variable',Found )
      IF(.NOT. Found ) THEN
        CALL Fatal('ChangeToHarmonicSystem','We need > Imaginary Variable < to create harmonic system!')
      END IF
    ELSE
      Name = TRIM( SaveVar % Name )//' im'
      CALL Info('ChangeToHarmonicSystem','Using derived name for imaginary component: '//TRIM(Name),Level=12)
    END IF

    TmpVals => HarmVar % Values(2:2*n:2)
    CALL VariableAdd( Solver % Mesh % Variables,Solver % Mesh,Solver, &
        Name, DOFs,TmpVals, Perm=SaveVar % Perm,Output=.TRUE.,Secondary=.TRUE.)        

    IF( Dofs > 1 ) THEN
      DO i=1,Dofs
        TmpVals => HarmVar % Values(2*i:2*n:2*Dofs)
        CALL VariableAdd( Solver % Mesh % Variables,Solver % Mesh,Solver, &
            ComponentName(Name,i),1,TmpVals,Perm=SaveVar % Perm,Output=.TRUE.,Secondary=.TRUE.)        
      END DO
    END IF
    
  END IF

  ! Now change the pointers such that when we visit the linear solver
  ! the system will automatically be solved as complex
  Solver % Variable => HarmVar
  Solver % Matrix => Aharm
  
  ! Save the original matrix and variable in Ematrix and Evar
  Solver % Matrix % Ematrix => SaveMatrix
  Solver % Variable % Evar => SaveVar    

  ! Eliminate cyclic dependence that is a bummer when deallocating stuff
  ! We are toggling {Are,Aharm} in {Solver % Matrix, Solver % Matrix % Ematrix}
  NULLIFY( Solver % Matrix % EMatrix % Ematrix )

!------------------------------------------------------------------------------
END SUBROUTINE ChangeToHarmonicSystem
!------------------------------------------------------------------------------

 

!------------------------------------------------------------------------------
!>  This subroutine will solve the system with some linear restriction.
!>  The restriction matrix is assumed to be in the ConstraintMatrix-field of 
!>  the StiffMatrix. The restriction vector is the RHS-field of the
!>  ConstraintMatrix.
!>  NOTE: Only serial solver implemented so far ...
!------------------------------------------------------------------------------
RECURSIVE SUBROUTINE SolveWithLinearRestriction( StiffMatrix, ForceVector, Solution, &
        Norm, DOFs, Solver )
!------------------------------------------------------------------------------  
  IMPLICIT NONE
  TYPE(Matrix_t), POINTER :: StiffMatrix !< Linear equation matrix information. 
                                         !< The restriction matrix is assumed to be in the EMatrix-field
  REAL(KIND=dp),TARGET :: ForceVector(:)        !< The right hand side of the linear equation
  REAL(KIND=dp),TARGET :: Solution(:)           !< Previous solution as input, new solution as output.
  REAL(KIND=dp) :: Norm                  !< The L2 norm of the solution.
  INTEGER :: DOFs                        !< Number of degrees of freedom of the equation.
  TYPE(Solver_t), TARGET :: Solver       !< Linear equation solver options.
!------------------------------------------------------------------------------
  TYPE(Solver_t), POINTER :: SolverPointer
  TYPE(Matrix_t), POINTER :: CollectionMatrix, RestMatrix, AddMatrix, &
       RestMatrixTranspose, TMat, XMat
  REAL(KIND=dp), POINTER CONTIG :: CollectionVector(:), RestVector(:),&
     AddVector(:), Tvals(:), Vals(:)
  REAL(KIND=dp), POINTER  :: MultiplierValues(:)
  REAL(KIND=dp), ALLOCATABLE, TARGET :: CollectionSolution(:), TotValues(:)
  INTEGER :: NumberOfRows, NumberOfValues, MultiplierDOFs, istat, NoEmptyRows 
  INTEGER :: i, j, k, l, m, n, p,q, ix, Loop
  TYPE(Variable_t), POINTER :: MultVar
  REAL(KIND=dp) :: scl, rowsum
  LOGICAL :: Found, ExportMultiplier, NotExplicit, Refactorize, EnforceDirichlet, EliminateDiscont, &
              NonEmptyRow, ComplexSystem, ConstraintScaling, UseTranspose, EliminateConstraints, &
              SkipConstraints
  SAVE MultiplierValues, SolverPointer

  CHARACTER(LEN=MAX_NAME_LEN) :: MultiplierName, str
  TYPE(ListMatrix_t), POINTER :: cList
  TYPE(ListMatrixEntry_t), POINTER :: cPtr, cPrev, cTmp

  INTEGER, ALLOCATABLE, TARGET :: SlavePerm(:), SlaveIPerm(:), MasterPerm(:), MasterIPerm(:)
  INTEGER, POINTER :: UsePerm(:), UseIPerm(:)
  REAL(KIND=dp), POINTER :: UseDiag(:)
  TYPE(ListMatrix_t), POINTER :: Lmat(:)
  LOGICAL  :: EliminateFromMaster, EliminateSlave, Parallel, UseTreeGauge
  REAL(KIND=dp), ALLOCATABLE, TARGET :: SlaveDiag(:), MasterDiag(:), DiagDiag(:)
  LOGICAL, ALLOCATABLE :: TrueDof(:)
  INTEGER, ALLOCATABLE :: Iperm(:)
  CHARACTER(*), PARAMETER :: Caller = 'SolveWithLinearRestriction'

  
!------------------------------------------------------------------------------
  CALL Info( Caller, ' ', Level=5 )
  SolverPointer => Solver
 
  
  Parallel = (ParEnv % PEs > 1 )

  NotExplicit = ListGetLogical(Solver % Values,'No Explicit Constrained Matrix',Found)
  IF(.NOT. Found) NotExplicit=.FALSE.

  RestMatrix => NULL()
  IF(.NOT.NotExplicit) &
        RestMatrix => StiffMatrix % ConstraintMatrix
  RestVector => Null()
  IF(ASSOCIATED(RestMatrix)) RestVector => RestMatrix % RHS

  AddMatrix => StiffMatrix % AddMatrix
  AddVector => NULL()
  IF(ASSOCIATED(AddMatrix)) &
    AddVector => AddMatrix % RHS

  NumberOfRows = StiffMatrix % NumberOfRows
  
  CollectionMatrix => StiffMatrix % CollectionMatrix
  Refactorize = ListGetLogical(Solver % Values,'Linear System Refactorize',Found)
  IF(.NOT.Found) Refactorize = .TRUE.

  IF(ASSOCIATED(CollectionMatrix)) THEN
    IF(Refactorize.AND..NOT.NotExplicit) THEN
      CALL Info( Caller,'Freeing previous collection matrix structures',Level=10)
      CALL FreeMatrix(CollectionMatrix)
      CollectionMatrix => NULL()
    ELSE
      CALL Info( Caller,'Keeping previous collection matrix structures',Level=10)
    END IF
  END IF

  IF(.NOT.ASSOCIATED(CollectionMatrix)) THEN
    CollectionMatrix => AllocateMatrix()
    CollectionMatrix % FORMAT = MATRIX_LIST
  ELSE
    DEALLOCATE(CollectionMatrix % RHS)
    CollectionMatrix % Values = 0.0_dp
  END IF
  IF(NotExplicit) CollectionMatrix % ConstraintMatrix => StiffMatrix % ConstraintMatrix  
  
  NumberOfRows = StiffMatrix % NumberOfRows
  IF(ASSOCIATED(AddMatrix)) NumberOfRows = MAX(NumberOfRows,AddMatrix % NumberOfRows)
  EliminateConstraints = ListGetLogical( Solver % Values, 'Eliminate Linear Constraints', Found)
  IF(ASSOCIATED(RestMatrix)) THEN
    IF(.NOT.EliminateConstraints) &
      NumberOfRows = NumberOFRows + RestMatrix % NumberOfRows
  END IF

  ALLOCATE( CollectionMatrix % RHS( NumberOfRows ), &
       CollectionSolution( NumberOfRows ), STAT = istat )
  IF ( istat /= 0 ) CALL Fatal( Caller, 'Memory allocation error.' )

  CollectionVector => CollectionMatrix % RHS
  CollectionVector = 0.0_dp
  CollectionSolution = 0.0_dp

!------------------------------------------------------------------------------
! If multiplier should be exported,  allocate memory and export the variable.
!------------------------------------------------------------------------------

  ExportMultiplier = ListGetLogical( Solver % Values, 'Export Lagrange Multiplier', Found )
  IF ( .NOT. Found ) ExportMultiplier = .FALSE.


  IF ( ExportMultiplier ) THEN
     MultiplierName = ListGetString( Solver % Values, 'Lagrange Multiplier Name', Found )
     IF ( .NOT. Found ) THEN
        CALL Info( Caller, &
              'Lagrange Multiplier Name set to LagrangeMultiplier', Level=5 )
        MultiplierName = "LagrangeMultiplier"
     END IF

     MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)
     j = 0
     IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberofRows
     IF(ASSOCIATED(AddMatrix))  j = j+MAX(0,AddMatrix % NumberofRows-StiffMatrix % NumberOfRows)

     IF ( .NOT. ASSOCIATED(MultVar) ) THEN
       ALLOCATE( MultiplierValues(j), STAT=istat )
       IF ( istat /= 0 ) CALL Fatal(Caller,'Memory allocation error.')

       MultiplierValues = 0.0_dp
       CALL VariableAdd(Solver % Mesh % Variables, Solver % Mesh, SolverPointer, &
                  MultiplierName, 1, MultiplierValues)
     END IF
     MultVar => VariableGet(Solver % Mesh % Variables, MultiplierName)

     MultiplierValues => MultVar % Values

     IF (j>SIZE(MultiplierValues)) THEN
       ALLOCATE(MultiplierValues(j)); MultiplierValues=0._dp
       MultiplierValues(1:SIZE(MultVar % Values)) = MultVar % Values
       DEALLOCATE(MultVar % Values)
       MultVar % Values => MultiplierValues
     END IF
  ELSE
     MultiplierValues => NULL()
  END IF

  UseTreeGauge = ListGetlogical( Solver % Values, 'Use Tree Gauge', Found )

!------------------------------------------------------------------------------
! Put the RestMatrix to lower part of CollectionMatrix
!------------------------------------------------------------------------------

  EnforceDirichlet = ListGetLogical( Solver % Values, 'Enforce Exact Dirichlet BCs',Found)
  IF(.NOT.Found) EnforceDirichlet = .TRUE.
  EnforceDirichlet = EnforceDirichlet .AND. ALLOCATED(StiffMatrix % ConstrainedDOF)

  ComplexSystem = StiffMatrix % COMPLEX
  ComplexSystem = ComplexSystem .OR. ListGetLogical( Solver % Values, &
           'Linear System Complex', Found )

  UseTranspose = ListGetLogical( Solver % Values, 'Use Transpose values', Found)

  IF(ASSOCIATED(RestMatrix).AND..NOT.EliminateConstraints) THEN

    CALL Info(Caller,'Adding ConstraintMatrix into CollectionMatrix',Level=8)
    CALL Info(Caller,'Number of Rows in constraint matrix: '&
        //TRIM(I2S(RestMatrix % NumberOfRows)),Level=12)
    CALL Info(Caller,'Number of Nofs in constraint matrix: '&
        //TRIM(I2S(SIZE(RestMatrix % Values))),Level=12)

    NoEmptyRows = 0
    ConstraintScaling = ListGetLogical(Solver % Values, 'Constraint Scaling',Found)
    IF(ConstraintScaling) THEN
      rowsum = ListGetConstReal( Solver % Values, 'Constraint Scale', Found)
      IF(Found) RestMatrix % Values = RestMatrix % Values * rowsum
    END IF

    ALLOCATE( iperm(SIZE(Solver % Variable % Perm)) )
    iperm = 0
    DO i=1,SIZE(Solver % Variable % Perm)
      IF ( Solver % Variable % Perm(i)>0) Iperm(Solver % Variable % Perm(i))=i
    END DO

    DO i=RestMatrix % NumberOfRows,1,-1

      k=StiffMatrix % NumberOfRows
      IF(ASSOCIATED(AddMatrix)) k=MAX(k,AddMatrix % NumberOfRows)
      k=k+i

      CALL AddToMatrixElement( CollectionMatrix,k,k,0._dp )
      IF(ComplexSystem) THEN
        IF(MOD(k,2)==0) THEN
          CALL AddToMatrixElement( CollectionMatrix,k,k-1,0._dp )
        ELSE
          CALL AddToMatrixElement( CollectionMatrix,k,k+1,0._dp )
        END IF
      END IF
      NonEmptyRow = .FALSE.

      rowsum = 0._dp
      l = -1
      DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
        IF(RestMatrix % Cols(j)==k) l=j
        rowsum = rowsum + ABS(RestMatrix % Values(j))
      END DO

      IF(rowsum>EPSILON(1._dp)) THEN
        IF(ConstraintScaling) THEN
          IF(l<=0.OR.l>0.AND.RestMatrix % Values(l)==0) THEN
            DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
              RestMatrix % Values(j) = RestMatrix % values(j)/rowsum
            END DO
            RestMatrix % RHS(i) = RestMatrix % RHS(i) / rowsum
          END IF
        END IF

        DO j=RestMatrix % Rows(i+1)-1,RestMatrix % Rows(i),-1
          Found = .TRUE.

          ! Skip non-positive column indexes
          IF( RestMatrix % Cols(j) <= 0 ) CYCLE
          IF ( .NOT. ComplexSystem ) THEN
            IF( ABS(RestMatrix % Values(j)) < EPSILON(1._dp)*rowsum ) CYCLE
          END IF

          IF (EnforceDirichlet .AND. RestMatrix % Cols(j) <= StiffMatrix % NumberOfRows) &
                  Found = .NOT.StiffMatrix % ConstrainedDOF(RestMatrix % Cols(j))

          IF(Found) THEN
            IF (ASSOCIATED(RestMatrix % TValues)) THEN
              CALL AddToMatrixElement( CollectionMatrix, &
                 RestMatrix % Cols(j), k, RestMatrix % TValues(j))
            ELSE
              CALL AddToMatrixElement( CollectionMatrix, &
                 RestMatrix % Cols(j), k, RestMatrix % Values(j))
            END IF

            IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
              CALL AddToMatrixElement( CollectionMatrix, &
                       k, RestMatrix % Cols(j), RestMatrix % TValues(j))
              NonEmptyRow = NonEmptyRow .OR. RestMatrix % TValues(j) /= 0
            ELSE
              CALL AddToMatrixElement( CollectionMatrix, &
                      k, RestMatrix % Cols(j), RestMatrix % Values(j))
              NonEmptyRow = NonEmptyRow .OR. RestMatrix % Values(j) /= 0
            END IF
          ELSE
            IF (UseTranspose .AND. ASSOCIATED(RestMatrix % TValues)) THEN
              CollectionVector(k) = CollectionVector(k) - &
                        RestMatrix % TValues(j) * ForceVector(RestMatrix % Cols(j)) / &
                           StiffMatrix % Values(StiffMatrix % Diag(RestMatrix % Cols(j)))
!            CALL AddToMatrixElement( CollectionMatrix, &
!                 k, RestMatrix % Cols(j), RestMatrix % TValues(j) )
!            NonEmptyRow = NonEmptyRow .OR. RestMatrix % TValues(j) /= 0
            ELSE
              CollectionVector(k) = CollectionVector(k) - &
                        RestMatrix % Values(j) * ForceVector(RestMatrix % Cols(j)) / &
                           StiffMatrix % Values(StiffMatrix % Diag(RestMatrix % Cols(j)))
!             CALL AddToMatrixElement( CollectionMatrix, &
!                 k, RestMatrix % Cols(j), RestMatrix % Values(j) )
!             NonEmptyRow = NonEmptyRow .OR. RestMatrix % Values(j) /= 0
            END IF
          END IF

        END DO
      END IF
 
      Found = .TRUE.
      IF (EnforceDirichlet) THEN
        IF(ASSOCIATED(RestMatrix % InvPerm)) THEN
          l = RestMatrix % InvPerm(i)
          IF(l>0) THEN
            l = MOD(l-1,StiffMatrix % NumberOfRows)+1
            IF(StiffMatrix % ConstrainedDOF(l)) THEN
              l = iperm((l-1)/Solver % Variable % DOFs+1) 
              IF (l<=Solver % Mesh % NumberOfNodes) THEN
                Found = .FALSE.
                CALL ZeroRow(CollectionMatrix,k)
                CollectionVector(k) = 0
                CALL SetMatrixElement(CollectionMatrix,k,k,1._dp)
              END IF
            END IF
          END IF
        END IF
      END IF

      ! If there is no matrix entry, there can be no non-zero r.h.s.
      IF ( Found ) THEN
        IF( .NOT.NonEmptyRow ) THEN
          NoEmptyRows = NoEmptyRows + 1
          CollectionVector(k) = 0._dp
!          might not be the right thing to do in parallel!!
          IF(UseTreeGauge) THEN
            CALL SetMatrixElement( CollectionMatrix,k,k,1._dp )
          END IF 
        ELSE
          IF( ASSOCIATED( RestVector ) ) CollectionVector(k) = CollectionVector(k) + RestVector(i)
        END IF
      END IF
    END DO

    IF( NoEmptyRows > 0 ) THEN
      CALL Info(Caller,&
          'Constraint Matrix in partition '//TRIM(I2S(ParEnv % MyPe))// &
          ' has '//TRIM(I2S(NoEmptyRows))// &
          ' empty rows out of '//TRIM(I2S(RestMatrix % NumberOfRows)), &
	  Level=6 )
    END IF

    CALL Info(Caller,'Finished Adding ConstraintMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the AddMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  IF(ASSOCIATED(AddMatrix)) THEN

    CALL Info(Caller,'Adding AddMatrix into CollectionMatrix',Level=10)

    DO i=AddMatrix % NumberOfRows,1,-1

      Found = .TRUE.
      IF (EnforceDirichlet .AND. i<=StiffMatrix % NumberOFRows) &
         Found = .NOT.StiffMatrix % ConstrainedDOF(i)

      IF(Found) THEN
        Found = .FALSE.
        DO j=AddMatrix % Rows(i+1)-1,AddMatrix % Rows(i),-1
            CALL AddToMatrixElement( CollectionMatrix, &
               i, AddMatrix % Cols(j), AddMatrix % Values(j))
            IF (i == AddMatrix % Cols(j)) Found = .TRUE.
        END DO

        CollectionVector(i) = CollectionVector(i) + AddVector(i)
        IF (.NOT.Found) THEN
          CALL AddToMatrixElement( CollectionMatrix, i, i, 0._dp )
          IF(ComplexSystem) THEN
            IF(MOD(i,2)==0) THEN
              CALL AddToMatrixElement( CollectionMatrix,i,i-1,0._dp )
            ELSE
              CALL AddToMatrixElement( CollectionMatrix,i,i+1,0._dp )
            END IF
          END IF
        END IF
      END IF
    END DO
    CALL Info(Caller,'Finished Adding AddMatrix',Level=12)
  END IF

!------------------------------------------------------------------------------
! Put the StiffMatrix to upper part of CollectionMatrix
!------------------------------------------------------------------------------
  CALL Info(Caller,'Adding Stiffness Matrix into CollectionMatrix',Level=10)

  DO i=StiffMatrix % NumberOfRows,1,-1
    DO j=StiffMatrix % Rows(i+1)-1,StiffMatrix % Rows(i),-1
      CALL AddToMatrixElement( CollectionMatrix, &
        i, StiffMatrix % Cols(j), StiffMatrix % Values(j) )
    END DO
    CollectionVector(i) = CollectionVector(i) + ForceVector(i)
  END DO

!------------------------------------------------------------------------------
! Eliminate constraints instead of adding the Lagrange coefficient equations.
! Assumes biorthogonal basis for Lagrange coefficient interpolation, but not
! necessarily biorthogonal constraint equation test functions.
!------------------------------------------------------------------------------
  IF (ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
    CALL Info(Caller,'Eliminating Constraints from CollectionMatrix',Level=10)

    n = StiffMatrix % NumberOfRows
    m = RestMatrix % NumberOfRows

    ALLOCATE(SlaveDiag(m),MasterDiag(m),SlavePerm(n),MasterPerm(n),&
        SlaveIPerm(m),MasterIPerm(m),DiagDiag(m))
    SlavePerm  = 0; SlaveIPerm  = 0; 
    MasterPerm = 0; MasterIPerm = 0

    Tvals => RestMatrix % TValues
    IF (.NOT.ASSOCIATED(Tvals)) Tvals => RestMatrix % Values 

    ! Extract diagonal entries for constraints:
    !------------------------------------------
    CALL Info(Caller,'Extracting diagonal entries for constraints',Level=15)
    DO i=1, RestMatrix % NumberOfRows
      m = RestMatrix % InvPerm(i)

      IF( m == 0 ) THEN
        PRINT *,'InvPerm is zero:',ParEnv % MyPe, i
        CYCLE
      END IF

      m = MOD(m-1,n) + 1
      SlavePerm(m)  = i
      SlaveIperm(i) = m

      DO j=RestMatrix % Rows(i), RestMatrix % Rows(i+1)-1
        k = RestMatrix % Cols(j)
        IF(k>n) THEN
           DiagDiag(i) = Tvals(j)
           CYCLE
        END IF

        IF( ABS( TVals(j) ) < TINY( 1.0_dp ) ) THEN
          PRINT *,'Tvals too small',ParEnv % MyPe,j,i,k,RestMatrix % InvPerm(i),Tvals(j)
        END IF

        IF(k == RestMatrix % InvPerm(i)) THEN
           SlaveDiag(i) = Tvals(j)
        ELSE
           MasterDiag(i) = Tvals(j)
           MasterPerm(k)  = i
           MasterIperm(i) = k
        END IF
      END DO
    END DO
  END IF

  IF (ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
    EliminateSlave = ListGetLogical( Solver % values, 'Eliminate Slave',Found )
    EliminateFromMaster = ListGetLogical( Solver % values, 'Eliminate From Master',Found )

    IF(EliminateFromMaster) THEN
      CALL Info(Caller,'Eliminating from master',Level=15)      
      UsePerm  => MasterPerm 
      UseDiag  => MasterDiag
      UseIPerm => MasterIPerm 
    ELSE
      CALL Info(Caller,'Eliminating from slave',Level=15)            
      UsePerm  => SlavePerm
      UseDiag  => SlaveDiag
      UseIPerm => SlaveIPerm
    END IF

    IF(UseTranspose) THEN
      Vals => Tvals
    ELSE
      Vals => RestMatrix % Values
    END IF
  END IF

  IF ( ParEnv % Pes>1 ) THEN
    EliminateDiscont =  ListGetLogical( Solver % values, 'Eliminate Discont',Found )
    IF( EliminateDiscont ) THEN
      CALL totv( StiffMatrix, SlaveDiag, SlaveIPerm )
      CALL totv( StiffMatrix, DiagDiag, SlaveIPerm )
      CALL totv( StiffMatrix, MasterDiag, MasterIPerm )
      CALL tota( StiffMatrix, TotValues, SlavePerm )
    END IF
  ELSE
    EliminateDiscont = .FALSE.
  END IF

  IF (ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
    ! Replace elimination equations by the constraints (could done be as a postprocessing
    ! step, if eq's totally eliminated from linsys.)
    ! ----------------------------------------------------------------------------------
    CALL Info(Caller,'Deleting rows from equation to be eliminated',Level=15)

    Lmat => CollectionMatrix % ListMatrix
    DO m=1,RestMatrix % NumberOfRows
      i = UseIPerm(m)
      CALL List_DeleteRow(Lmat, i, Keep=.TRUE.)
    END DO

    CALL Info(Caller,'Copying rows from constraint matrix to eliminate dofs',Level=15)
    DO m=1,RestMatrix % NumberOfRows
      i = UseIPerm(m)
      DO l=RestMatrix % Rows(m+1)-1, RestMatrix % Rows(m), -1
        j = RestMatrix % Cols(l)

        ! skip l-coeffient entries, handled separately afterwards:
        ! --------------------------------------------------------
        IF(j > n) CYCLE

        CALL List_AddToMatrixElement( Lmat, i, j, Vals(l) )
      END DO
      CollectionVector(i) = RestVector(m)
    END DO

    ! Eliminate slave dof cycles:
    ! ---------------------------
    Xmat => RestMatrix
    Found = .TRUE.
    Loop = 0
    DO WHILE(Found)
      DO i=Xmat % NumberofRows,1,-1
        q = 0
        DO j = Xmat % Rows(i+1)-1, Xmat % Rows(i),-1
          k = Xmat % Cols(j)
          IF(k>n) CYCLE
          IF(UsePerm(k)>0 .AND. ABS(TVals(j))>AEPS) q=q+1
        END DO
        IF(q>1) EXIT
      END DO
      Found = q>1

      Tmat => Xmat
      IF(Found) THEN
        Loop = Loop + 1
        CALL Info(Caller,'Recursive elimination round: '//TRIM(I2S(Loop)),Level=15)

        Tmat => AllocateMatrix()
        Tmat % Format = MATRIX_LIST

        DO i=Xmat % NumberofRows,1,-1
          DO j = Xmat % Rows(i+1)-1, Xmat % Rows(i),-1
            k = Xmat % Cols(j)
            IF ( ABS(Tvals(j))>AEPS ) &
              CALL List_AddToMatrixElement(Tmat % ListMatrix, i, k, TVals(j))
          END DO
        END DO


        DO m=1,Xmat % NumberOfRows
          i = UseIPerm(m)
          DO j=Xmat % Rows(m), Xmat % Rows(m+1)-1
            k = Xmat % Cols(j)

            ! The size of SlavePerm is often exceeded but I don't really undersrtand the operation...
            ! so this is just a dirty fix.
            IF( k > SIZE( SlavePerm ) ) CYCLE

            l = SlavePerm(k)

            IF(l>0 .AND. k/=i) THEN
              IF(ABS(Tvals(j))<AEPS) CYCLE
              scl = -TVals(j) / SlaveDiag(l)

              CALL List_DeleteMatrixElement( Tmat % ListMatrix, m, k )

              DO q=Xmat % Rows(l+1)-1, Xmat % Rows(l),-1
                IF(ABS(Tvals(q))<AEPS) CYCLE
                ix = Xmat % Cols(q)
                IF ( ix/=k ) &
                  CALL List_AddToMatrixElement( Tmat % ListMatrix, m, ix, scl * TVals(q) )
              END DO
            END IF
          END DO
        END DO

        CALL List_ToCRSMatrix(Tmat)
        Tvals => Tmat % Values
        IF(.NOT.ASSOCIATED(Xmat,RestMatrix)) CALL FreeMatrix(Xmat)
      END IF
      Xmat => TMat
    END DO

    ! Eliminate Lagrange Coefficients:
    ! --------------------------------

    CALL Info(Caller,'Eliminating Largrange Coefficients',Level=15)

    DO m=1,Tmat % NumberOfRows
      i = UseIPerm(m)
      IF( ABS( UseDiag(m) ) < TINY( 1.0_dp ) ) THEN
        PRINT *,'UseDiag too small:',m,ParEnv % MyPe,UseDiag(m)
        CYCLE
      END IF

      DO j=TMat % Rows(m), TMat % Rows(m+1)-1
        k = TMat % Cols(j)
        IF(k<=n) THEN
          IF(UsePerm(k)/=0) CYCLE

          IF ( EliminateDiscont ) THEN
            IF (EliminateFromMaster) THEN
              scl = -SlaveDiag(SlavePerm(k)) / UseDiag(m)
            ELSE
              scl = -MasterDiag(MasterPerm(k)) / UseDiag(m)
            END IF
          ELSE
            scl = -Tvals(j) / UseDiag(m)
          END IF
        ELSE
          k = UseIPerm(k-n)
          ! multiplied by 1/2 in GenerateConstraintMatrix()
          IF (EliminateDiscont) THEN
            scl = -2*DiagDiag(m) / UseDiag(m)
          ELSE
            scl = -2*Tvals(j) / UseDiag(m)
          END IF
        END IF

        DO l=StiffMatrix % Rows(i+1)-1, StiffMatrix % Rows(i),-1
          CALL List_AddToMatrixElement( Lmat, k, &
              StiffMatrix % Cols(l), scl * StiffMatrix % Values(l) )
        END DO
        CollectionVector(k) = CollectionVector(k) + scl * ForceVector(i)
      END DO
    END DO

    IF ( .NOT.ASSOCIATED(Tmat, RestMatrix ) ) CALL FreeMatrix(Tmat)

    ! Eliminate slave dofs, using the constraint equations:
    ! -----------------------------------------------------
    IF ( EliminateSlave ) THEN
      CALL Info(Caller,'Eliminate slave dofs using constraint equations',Level=15)

      IF(EliminateDiscont) THEN
        DO i=1,StiffMatrix % NumberOfRows
          IF ( UsePerm(i)/=0 ) CYCLE

          DO m=StiffMatrix % Rows(i), StiffMatrix % Rows(i+1)-1
             j = SlavePerm(StiffMatrix % Cols(m))
             IF ( j==0 ) CYCLE
             scl = -TotValues(m) / SlaveDiag(j)

             ! Delete elimination entry:
             ! -------------------------
             CALL List_DeleteMatrixElement(Lmat,i,StiffMatrix % Cols(m))

             k = UseIPerm(j)
             cTmp => Lmat(k) % Head
             DO WHILE(ASSOCIATED(cTmp))
                l = cTmp % Index
                IF ( l /= SlaveIPerm(j) ) &
                   CALL List_AddToMatrixElement( Lmat, i, l, scl*cTmp % Value )
              cTmp => cTmp % Next
            END DO
            CollectionVector(i) = CollectionVector(i) + scl * CollectionVector(k)
          END DO
        END DO
      ELSE

        CALL List_ToCRSMatrix(CollectionMatrix)
        Tmat => AllocateMatrix()
        Tmat % Format = MATRIX_LIST

        DO i=1,StiffMatrix % NumberOfRows
          IF(UsePerm(i)/=0) CYCLE

          DO m = CollectionMatrix % Rows(i), CollectionMatrix % Rows(i+1)-1
            j = SlavePerm(CollectionMatrix % Cols(m))

            IF(j==0) THEN
              CYCLE
            END IF
            IF( ABS( SlaveDiag(j) ) < TINY( 1.0_dp ) ) THEN
              PRINT *,'SlaveDiag too small:',j,ParEnv % MyPe,SlaveDiag(j)
              CYCLE
            END IF

            scl = -CollectionMatrix % Values(m) / SlaveDiag(j)
            CollectionMatrix % Values(m) = 0._dp

            ! ... and add replacement values:
            ! -------------------------------
            k = UseIPerm(j)
            DO p=CollectionMatrix % Rows(k+1)-1, CollectionMatrix % Rows(k), -1
               l = CollectionMatrix % Cols(p)
               IF ( l /= SlaveIPerm(j) ) &
                 CALL List_AddToMatrixElement( Tmat % listmatrix, i, l, scl*CollectionMatrix % Values(p) )
            END DO
            CollectionVector(i) = CollectionVector(i) + scl * CollectionVector(k)
          END DO
        END DO

        CALL List_ToListMatrix(CollectionMatrix)
        Lmat => CollectionMatrix % ListMatrix

        CALL List_ToCRSMatrix(Tmat)
        DO i=TMat % NumberOfRows,1,-1
          DO j=TMat % Rows(i+1)-1,TMat % Rows(i),-1
            CALL List_AddToMatrixElement( Lmat, i, TMat % cols(j), TMat % Values(j) )
          END DO
        END DO
        CALL FreeMatrix(Tmat)
      END IF
    END IF

    ! Optimize bandwidth, if needed:
    ! ------------------------------
    IF(EliminateFromMaster) THEN
      CALL Info(Caller,&
          'Optimizing bandwidth after elimination',Level=15)
      DO i=1,RestMatrix % NumberOfRows
        j = SlaveIPerm(i)
        k = MasterIPerm(i)

        Ctmp => Lmat(j) % Head
        Lmat(j) % Head => Lmat(k) % Head
        Lmat(k) % Head => Ctmp

        l = Lmat(j) % Degree
        Lmat(j) % Degree = Lmat(k) % Degree
        Lmat(k) % Degree = l

        scl = CollectionVector(j)
        CollectionVector(j) = CollectionVector(k)
        CollectionVector(k) = scl
      END DO
    END IF

    CALL Info(Caller,'Finished Adding ConstraintMatrix',Level=12)
  END IF

  CALL Info(Caller,'Reverting CollectionMatrix back to CRS matrix',Level=10)
  IF(CollectionMatrix % FORMAT==MATRIX_LIST) &
      CALL List_toCRSMatrix(CollectionMatrix)

  CALL Info( Caller, 'CollectionMatrix done', Level=5 )

!------------------------------------------------------------------------------
! Assign values to CollectionVector
!------------------------------------------------------------------------------

  j = StiffMatrix % NumberOfRows  
  CollectionSolution(1:j) = Solution(1:j)
  
  i = StiffMatrix % NumberOfRows+1
  j = SIZE(CollectionSolution)
  CollectionSolution(i:j) = 0._dp
  IF(ExportMultiplier) CollectionSolution(i:j) = MultiplierValues(1:j-i+1)

  CollectionMatrix % ExtraDOFs = CollectionMatrix % NumberOfRows - &
                  StiffMatrix % NumberOfRows

  CollectionMatrix % ParallelDOFs = 0
  IF(ASSOCIATED(AddMatrix)) &
    CollectionMatrix % ParallelDOFs = MAX(AddMatrix % NumberOfRows - &
                  StiffMatrix % NumberOfRows,0)

  CALL Info( Caller, 'CollectionVector done', Level=5 )

!------------------------------------------------------------------------------
! Solve the Collection-system 
!------------------------------------------------------------------------------

! Collectionmatrix % Complex = StiffMatrix % Complex

  ! We may want to skip the constraints for norm if we use certain other options
  SkipConstraints = ListGetLogical( Solver % values, &
      'Nonlinear System Convergence Without Constraints',Found )
  IF(.NOT. Found ) THEN
    SkipConstraints = ListGetLogical( Solver % values, 'Linear System Residual Mode',Found ) 
    IF( SkipConstraints ) THEN
      CALL Info(Caller,'Linear system residual mode must skip constraints',Level=5)
    ELSE
      SkipConstraints = ListGetLogical( Solver % values, 'NonLinear System Consistent Norm',Found ) 
      IF( SkipConstraints ) THEN
        CALL Info(Caller,'Nonlinear system consistent norm must skip constraints',Level=5)
      END IF
    END IF
    str = ListGetString( Solver % values, 'NonLinear System Convergence Measure',Found )
    IF( str == 'solution' ) THEN
      SkipConstraints = .TRUE.
      CALL Info(Caller,&
          'Nonlinear system convergence measure == "solution" must skip constraints',Level=5)
    END IF
    IF( SkipConstraints ) THEN
      CALL Info(Caller,'Enforcing convergence without constraints to True',Level=5)
      CALL ListAddLogical( Solver % Values, &
           'Nonlinear System Convergence Without Constraints',.TRUE.)
    END IF
  END IF

  !------------------------------------------------------------------------------
  ! Look at the nonlinear system previous values again, not taking the constrained
  ! system into account...
  !------------------------------------------------------------------------------
  Found = ASSOCIATED(Solver % Variable % NonlinValues)
  IF( Found .AND. .NOT. SkipConstraints ) THEN
    k = CollectionMatrix % NumberOfRows
    IF ( SIZE(Solver % Variable % NonlinValues) /= k) THEN
      DEALLOCATE(Solver % Variable % NonlinValues)
      ALLOCATE(Solver % Variable % NonlinValues(k))
    END IF
    Solver % Variable % NonlinValues(1:k) = CollectionSolution(1:k)
  END IF

  CollectionMatrix % Comm = StiffMatrix % Comm

  CALL Info(Caller,'Now going for the coupled linear system',Level=10)

  CALL SolveLinearSystem( CollectionMatrix, CollectionVector, &
      CollectionSolution, Norm, DOFs, Solver, StiffMatrix )

  
  !-------------------------------------------------------------------------------
  ! For restricted systems study the norm without some block components.
  ! For example, excluding gauge constraints may give valuable information
  ! of the real accuracy of the unconstrained system. Currently just for info.
  !-------------------------------------------------------------------------------
  IF( ListGetLogical( Solver % Values,'Restricted System Norm',Found ) ) THEN
    ALLOCATE( TrueDof( CollectionMatrix % NumberOfRows ) )
    TrueDof = .TRUE.
    
    Norm = LinearSystemMaskedResidualNorm( CollectionMatrix, CollectionVector, &
        CollectionSolution, TrueDof, TrueDof )
    
    WRITE( Message,'(A,ES13.6)') 'Residual norm of the original system:',Norm
    CALL Info(Caller,Message, Level = 5 )
    
    IF( ListGetLogical( Solver % Values,'Restricted System Norm Skip Nodes',Found ) ) THEN
      i = 1
      j = MAXVAL( Solver % Variable % Perm(1:Solver % Mesh % NumberOfNodes) )
      CALL Info(Caller,'Skipping nodal dof range: '&
          //TRIM(I2S(i))//'-'//TRIM(I2S(j)),Level=8)
      TrueDof(i:j) = .FALSE.
    END IF

    IF( ListGetLogical( Solver % Values,'Restricted System Norm Skip Constraints',Found ) ) THEN
      i = StiffMatrix % NumberOfRows + 1
      j = CollectionMatrix % NumberOfRows      
      CALL Info(Caller,'Skipping constraints dof range: '&
          //TRIM(I2S(i))//'-'//TRIM(I2S(j)),Level=8)
      TrueDof(i:j) = .FALSE.
    END IF
    
    Norm = LinearSystemMaskedResidualNorm( CollectionMatrix, CollectionVector, &
        CollectionSolution, TrueDof, TrueDof )
    
    WRITE( Message,'(A,ES13.6)') 'Residual norm of the masked system:',Norm
    CALL Info(Caller,Message, Level = 5 )
    
    DEALLOCATE( TrueDof )
  END IF
    

  
!------------------------------------------------------------------------------
! Separate the solution from CollectionSolution
!------------------------------------------------------------------------------
    CALL Info(Caller,'Picking solution from collection solution',Level=10)

    Solution = 0.0_dp
    i = 1
    j = StiffMatrix % NumberOfRows
    Solution(i:j) = CollectionSolution(i:j)

    IF ( ExportMultiplier ) THEN
      i = StiffMatrix % NumberOfRows
      j=0
      IF(ASSOCIATED(RestMatrix)) j = RestMatrix % NumberOfRows
      IF(ASSOCIATED(AddMatrix)) &
        j=j+MAX(0,AddMatrix % NumberOfRows - StiffMatrix % NumberOFRows)

      MultiplierValues = 0.0_dp
      IF(ASSOCIATED(RestMatrix).AND.EliminateConstraints) THEN
        ! Compute eliminated l-coefficient values:
        ! ---------------------------------------
        DO i=1,RestMatrix % NumberOfRows
          scl = 1._dp / UseDiag(i)
          m = UseIPerm(i)
          MultiplierValues(i) = scl * ForceVector(m)
          DO j=StiffMatrix % Rows(m), StiffMatrix % Rows(m+1)-1
            MultiplierValues(i) = MultiplierValues(i) - &
              scl * StiffMatrix % Values(j) * Solution(StiffMatrix % Cols(j))
          END DO
        END DO
      ELSE
        MultiplierValues(1:j) = CollectionSolution(i+1:i+j)
      END IF

      IF(EliminateConstraints.AND.EliminateDiscont) THEN
        IF (EliminateFromMaster) THEN
          CALL totv(StiffMatrix,MultiplierValues,MasterIPerm)
        ELSE
          CALL totv(StiffMatrix,MultiplierValues,SlaveIPerm)
        END IF
      END IF
    END IF

!------------------------------------------------------------------------------

    StiffMatrix % CollectionMatrix => CollectionMatrix
    DEALLOCATE(CollectionSolution)
    CollectionMatrix % ConstraintMatrix => NULL()

    CALL Info( Caller, 'All done', Level=5 )

CONTAINS

  SUBROUTINE totv( A, totvalues, perm )
    type(matrix_t), pointer :: A
    real(kind=dp) :: totvalues(:)
    integer, allocatable :: perm(:)

    real(kind=dp), ALLOCATABLE :: x(:),r(:)
    INTEGER :: i,j,ng

    ng = A % NumberOfRows
!   ng = ParallelReduction(1._dp*MAXVAL(A % ParallelInfo % GLobalDOfs))
    ALLOCATE(x(ng),r(ng))

    x = 0._dp
    IF(ALLOCATED(perm)) THEN
      DO i=1,SIZE(perm)
        j = Perm(i)
        !j = a % parallelinfo % globaldofs(j)
        x(j) = totvalues(i)
      END DO
    END IF

    CALL ParallelSumVector(A, x)
!   CALL MPI_ALLREDUCE( x,r, ng, MPI_DOUBLE_PRECISION, MPI_SUM, ELMER_COMM_WORLD, i ); x=r

    IF(ALLOCATED(perm)) THEN
      DO i=1,SIZE(perm)
        j = Perm(i)
        !j = A % parallelinfo % globaldofs(j)
        totvalues(i) = x(j)
      END DO
    END IF
  END SUBROUTINE Totv
    

  SUBROUTINE Tota( A, TotValues, cperm )
     type(matrix_t), pointer :: A
     integer, allocatable :: cperm(:)
     real(kind=dp), ALLOCATABLE :: totvalues(:)

     INTEGER, POINTER :: Diag(:), Rows(:), Cols(:)
     LOGICAL ::  found
     INTEGER :: status(MPI_STATUS_SIZE)
     REAL(KIND=dp), ALLOCATABLE, TARGET :: rval(:)
     INTEGER, ALLOCATABLE :: cnt(:), rrow(:),rcol(:), perm(:)
     INTEGER :: i,j,k,l,m,ii,jj,proc,rcnt,nn, dof, dofs, Active, n, nm,ierr

     TYPE Buf_t
        REAL(KIND=dp), ALLOCATABLE :: gval(:)
        INTEGER, ALLOCATABLE :: grow(:),gcol(:)
     END TYPE Buf_t
     TYPE(Buf_t), POINTER :: buf(:)

     Diag => A % Diag
     Rows => A % Rows
     Cols => A % Cols

     n = A % NumberOfRows

     ALLOCATE(TotValues(SIZE(A % Values))); TotValues=A % Values

     IF (ParEnv  % PEs>1 ) THEN
       ALLOCATE(cnt(0:ParEnv % PEs-1))
       cnt = 0
       DO i=1,n
         DO j=Rows(i),Rows(i+1)-1
!          IF(Cols(j)<=nm .OR. Cols(j)>nm+n) CYCLE
           iF ( ALLOCATED(CPerm)) THEN
             IF(cperm(Cols(j))==0) CYCLE
           END IF
           IF(TotValues(j)==0) CYCLE

           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             DO k=1,SIZE(A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours)
               m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(k)
               IF ( m==ParEnv % myPE ) CYCLE
               cnt(m) = cnt(m)+1
             END DO 
           END IF
         END DO
       END DO

       ALLOCATE( buf(0:ParEnv % PEs-1) )
       DO i=0,ParEnv % PEs-1
         IF ( cnt(i) > 0 ) &
           ALLOCATE( Buf(i) % gval(cnt(i)), Buf(i) % grow(cnt(i)), Buf(i) % gcol(cnt(i)) )
       END DO

       cnt = 0
       DO i=1,n
         DO j=Rows(i),Rows(i+1)-1
!          IF(Cols(j)<=nm .OR. Cols(j)>nm+n) CYCLE
           iF ( ALLOCATED(CPerm)) THEN
             IF(cperm(Cols(j))==0) CYCLE
           END IF
           IF(TotValues(j)==0) CYCLE

           IF ( A % ParallelInfo % Interface(Cols(j)) ) THEN
             DO k=1,SIZE(A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours)
               m = A % ParallelInfo % NeighbourList(Cols(j)) % Neighbours(k)
               IF ( m==ParEnv % myPE ) CYCLE
               cnt(m) = cnt(m)+1
               Buf(m) % gcol(cnt(m)) = A % ParallelInfo % GlobalDOFs(Cols(j))
               Buf(m) % gval(cnt(m)) = TotValues(j)
               Buf(m) % grow(cnt(m)) = A % ParallelInfo % GlobalDOFs(i)
             END DO
           END IF
         END DO
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( ParEnv % IsNeighbour(i+1) ) THEN
           CALL MPI_BSEND( cnt(i), 1, MPI_INTEGER, i, 7001, ELMER_COMM_WORLD, status, ierr )
           IF ( cnt(i)>0 ) THEN
             CALL MPI_BSEND( Buf(i) % grow, cnt(i), MPI_INTEGER, &
                 i, 7002, ELMER_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gcol, cnt(i), MPI_INTEGER, &
                 i, 7003, ELMER_COMM_WORLD, status, ierr )

             CALL MPI_BSEND( Buf(i) % gval, cnt(i), MPI_DOUBLE_PRECISION, &
                 i, 7004, ELMER_COMM_WORLD, status, ierr )
           END IF
         END IF
       END DO

       DO i=0,ParEnv % PEs-1
         IF ( cnt(i)>0 ) &
           DEALLOCATE( Buf(i) % gval, Buf(i) % grow, Buf(i) % gcol )
       END DO
       DEALLOCATE( cnt,Buf )

       DO i=1,ParEnv % NumOfNeighbours
         CALL MPI_RECV( rcnt, 1, MPI_INTEGER, &
           MPI_ANY_SOURCE, 7001, ELMER_COMM_WORLD, status, ierr )

         IF ( rcnt>0 ) THEN
           IF(.NOT.ALLOCATED(rrow)) THEN
             ALLOCATE( rrow(rcnt), rcol(rcnt), rval(rcnt) )
           ELSE IF(SIZE(rrow)<rcnt) THEN
             DEALLOCATE(rrow,rcol,rval)
             ALLOCATE( rrow(rcnt), rcol(rcnt), rval(rcnt) )
           ENDIF

           proc = status(MPI_SOURCE)
           CALL MPI_RECV( rrow, rcnt, MPI_INTEGER, &
              proc, 7002, ELMER_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rcol, rcnt, MPI_INTEGER, &
              proc, 7003, ELMER_COMM_WORLD, status, ierr )

           CALL MPI_RECV( rval, rcnt, MPI_DOUBLE_PRECISION, &
              proc, 7004, ELMER_COMM_WORLD, status, ierr )

           DO j=1,rcnt
             l = SearchNode(A % ParallelInfo,rcol(j),Order=A % ParallelInfo % Gorder )
             IF ( l>0 ) THEN
               k = SearchNode(A % ParallelInfo,rrow(j),Order=A % ParallelInfo % Gorder )
               IF ( k>0 ) THEN
                 IF ( l>=k ) THEN
                   DO m=Diag(k),Rows(k+1)-1
                     IF ( Cols(m) == l ) THEN
                       TotValues(m) = TotValues(m) + rval(j)
                       EXIT
                     ELSE IF( Cols(m)>l) THEN
                       EXIT
                     END IF
                   END DO
                 ELSE
                   DO m=Rows(k),Diag(k)-1
                     IF ( Cols(m)==l ) THEN
                       TotValues(m) = TotValues(m) + rval(j)
                       EXIT
                     ELSE IF( Cols(m)>l) THEN
                       EXIT
                     END IF
                   END DO
                 END IF
               END IF
             END IF
           END DO
         END IF
       END DO
     END IF
   END SUBROUTINE tota

!------------------------------------------------------------------------------
  END SUBROUTINE SolveWithLinearRestriction
!------------------------------------------------------------------------------
      

!------------------------------------------------------------------------------
  SUBROUTINE SaveLinearSystem( Solver, Ain )
!------------------------------------------------------------------------------
    TYPE( Solver_t ) :: Solver
    TYPE(Matrix_t), POINTER, OPTIONAL :: Ain
!------------------------------------------------------------------------------    
    TYPE(Matrix_t), POINTER :: A
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: dumpfile, dumpprefix
    INTEGER, POINTER :: Perm(:)
    REAL(KIND=dp), POINTER :: Sol(:)
    INTEGER :: i
    LOGICAL :: SaveMass, SaveDamp, SavePerm, SaveSol, Found , Parallel, CNumbering
    CHARACTER(*), PARAMETER :: Caller = 'SaveLinearSystem'
!------------------------------------------------------------------------------

    CALL Info(Caller,'Saving linear system',Level=4)

    Parallel = ParEnv % PEs > 1

    Params => Solver % Values
    IF(.NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'Parameter list not associated!')
    END IF

    CNumbering = ListGetLogical(Params, 'Linear System Save Continuous Numbering',Found)

    IF( PRESENT(Ain)) THEN
      A => Ain
    ELSE
      A => Solver % Matrix
    END IF

    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal(Caller,'Matrix not assciated!')
    END IF

    SaveMass = ListGetLogical( Params,'Linear System Save Mass',Found)

    SaveDamp = ListGetLogical( Params,'Linear System Save Damp',Found)   

    dumpprefix = ListGetString( Params, 'Linear System Save Prefix', Found)
    IF(.NOT. Found ) dumpprefix = 'linsys'

    dumpfile = TRIM(dumpprefix)//'_a.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
    CALL Info(Caller,'Saving matrix to: '//TRIM(dumpfile))
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    CALL PrintMatrix(A,Parallel,Cnumbering,SaveMass=SaveMass,SaveDamp=SaveDamp)
    CLOSE(1)

    dumpfile = TRIM(dumpprefix)//'_b.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
    CALL Info(Caller,'Saving matrix rhs to: '//TRIM(dumpfile))
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    CALL PrintRHS(A, Parallel, CNumbering)
    CLOSE(1)
    
    SavePerm = ListGetLogical( Params,'Linear System Save Perm',Found)
    IF( SavePerm ) THEN
      Perm => Solver % Variable % Perm
      IF( .NOT. ASSOCIATED( Perm ) ) THEN
        CALL Warn(Caller,'Permuation not associated!')
        SavePerm = .FALSE.
      ELSE
        dumpfile = TRIM(dumpprefix)//'_perm.dat'
        IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
        CALL Info(Caller,'Saving permutation to: '//TRIM(dumpfile))
        OPEN(1,FILE=dumpfile, STATUS='Unknown')
        DO i=1,SIZE(Perm)
          WRITE(1,'(I0,A,I0)') i,' ',Perm(i)
        END DO
        CLOSE( 1 ) 
      END IF
    END IF


    SaveSol = ListGetLogical( Params,'Linear System Save Solution',Found)
    IF( SaveSol ) THEN
      Sol => Solver % Variable % Values
      IF( .NOT. ASSOCIATED( Sol ) ) THEN
        CALL Warn(Caller,'Solution not associated!')
        SaveSol = .FALSE.
      ELSE
        dumpfile = TRIM(dumpprefix)//'_sol.dat'
        IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
        CALL Info(Caller,'Saving solution to: '//TRIM(dumpfile))
        OPEN(1,FILE=dumpfile, STATUS='Unknown')
        DO i=1,SIZE(Sol)
          WRITE(1,'(I0,ES15.6)') i,Sol(i)
        END DO
        CLOSE( 1 ) 
      END IF
    END IF
    
    
    dumpfile = TRIM(dumpprefix)//'_sizes.dat'
    IF(Parallel) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))
    CALL Info(Caller,'Saving matrix sizes to: '//TRIM(dumpfile))
    OPEN(1,FILE=dumpfile, STATUS='Unknown')
    WRITE(1,*) A % NumberOfRows
    WRITE(1,*) SIZE(A % Values)
    IF( SavePerm ) WRITE(1,*) SIZE( Perm )    
    CLOSE(1)

    IF(Parallel) THEN
      dumpfile = TRIM(dumpprefix)//'_sizes.dat'
      CALL Info(Caller,'Saving matrix sizes to: '//TRIM(dumpfile))
      OPEN(1,FILE=dumpfile, STATUS='Unknown')
      WRITE(1,*) NINT(ParallelReduction(1._dP*A % ParMatrix % &
                           SplittedMatrix % InsideMatrix % NumberOfRows))
      WRITE(1,*) NINT(ParallelReduction(1._dp*SIZE(A % Values)))
      IF( SavePerm ) WRITE(1,*) NINT(ParallelReduction(1._dp*SIZE( Perm )))
      CLOSE(1)
    END IF
    
    IF( ListGetLogical( Params,'Linear System Save and Stop',Found ) ) THEN
      CALL Info(Caller,'Just saved matrix and stopped!',Level=4)
      STOP EXIT_OK
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SaveLinearSystem
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Assemble Laplace matrix related to a solver and permutation vector. 
!------------------------------------------------------------------------------
  SUBROUTINE LaplaceMatrixAssembly( Solver, Perm, A )
    
    TYPE(Solver_t) :: Solver
    INTEGER, POINTER :: Perm(:)
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    !------------------------------------------------------------------------------

    INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
    INTEGER :: i,j,k,n,t,istat,BoundaryNodes
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    CHARACTER(LEN=MAX_NAME_LEN) :: BoundaryName
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
    REAL(KIND=dp) :: detJ, val
    LOGICAL :: Stat
    
    
    Mesh => Solver % Mesh
        
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), dBasisdx(n, 3), FORCE(N), STIFF(N,N), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)
    
    IF(.FALSE.) THEN
      N = Mesh % NumberOfNodes
      ALLOCATE( BoundaryPerm(n) )
      BoundaryPerm = 0
      BoundaryNodes = 0
      BoundaryName = 'Laplace Boundary'
      CALL MakePermUsingMask( CurrentModel,Solver,Mesh,BoundaryName, &
          .FALSE., BoundaryPerm, BoundaryNodes )
    END IF


    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes
      IF( ANY( Perm(Indexes) == 0 ) ) CYCLE

      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)

      STIFF = 0.0d0
      FORCE = 0.0d0
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )
      DO k=1,IP % n
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis, dBasisdx )
        
        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        DO i=1,n
          val = IP % s(k) * DetJ 
          
          ! This condition removes the natural boundary condition that would 
          ! try to fix the normal gradient of the field to zero.
          !--------------------------------------------------------------------
          IF(.FALSE.) THEN
            IF( BoundaryPerm( Indexes(i) ) > 0 ) CYCLE
          END IF

          DO j=1,n
            STIFF(i,j) = STIFF(i,j) + val * &
                SUM( dBasisdx(i,:) * dBasisdx(j,:) ) 
          END DO
        END DO
      END DO
      
      CALL UpdateGlobalEquations( A,STIFF,A % rhs,FORCE,n,1,Perm(Indexes(1:n)) )
    END DO
    
    DEALLOCATE( Basis, dBasisdx, FORCE, STIFF, & 
        Nodes % x, Nodes % y, Nodes % z)

  END SUBROUTINE LaplaceMatrixAssembly

  
!------------------------------------------------------------------------------
!> Assemble mass matrix related to a solver and permutation vector. 
!------------------------------------------------------------------------------
  SUBROUTINE MassMatrixAssembly( Solver, Perm, A )
    
    TYPE(Solver_t) :: Solver
    INTEGER, POINTER :: Perm(:)
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    !------------------------------------------------------------------------------

    INTEGER, POINTER :: Indexes(:)
    INTEGER :: i,j,k,n,t,istat
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: STIFF(:,:), FORCE(:)
    REAL(KIND=dp), POINTER :: Basis(:),rhs(:)
    REAL(KIND=dp) :: detJ, val
    LOGICAL :: Stat
    
    
    Mesh => Solver % Mesh
        
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), FORCE(N), STIFF(N,N), &
        Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)

    ALLOCATE( rhs(A % NumberOfRows) )
    rhs = 0.0_dp

    DO t=1,Mesh % NumberOfBulkElements
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes
      IF( ANY( Perm(Indexes) == 0 ) ) CYCLE

      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)

      STIFF = 0.0d0
      FORCE = 0.0d0
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )

      DO k=1,IP % n

        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis )
        
        ! Finally, the elemental matrix & vector:
        !----------------------------------------
        DO i=1,n
          val = IP % s(k) * DetJ 
          DO j=1,n
            STIFF(i,j) = STIFF(i,j) + val * Basis(i) * Basis(j)
          END DO
        END DO
      END DO

      CALL UpdateGlobalEquations( A,STIFF,rhs,FORCE,n,1,Perm(Indexes(1:n)) )
    END DO

    DEALLOCATE( Basis, FORCE, STIFF, & 
        Nodes % x, Nodes % y, Nodes % z)
    DEALLOCATE( rhs )

  END SUBROUTINE MassMatrixAssembly



!------------------------------------------------------------------------------
!> Create diagonal matrix from P (not square) by summing the entries together
!> and multiplying with a constant.
!------------------------------------------------------------------------------
  SUBROUTINE DiagonalMatrixSumming( Solver, P, A )
    
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER :: P, A
    !------------------------------------------------------------------------------
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: val, rowsum, minsum, maxsum, sumsum

    CALL Info('DiagonalMatrixSumming','Creating diagonal matrix from absolute rowsums')

    IF(.NOT. ASSOCIATED(P) ) CALL Fatal('DiagonalMatrixSumming','Matrix P not associated!')
    IF(.NOT. ASSOCIATED(A) ) CALL Fatal('DiagonalMatrixSumming','Matrix A not associated!')

    
    n = P % NumberOfRows 
    CALL Info('DiagonalMatrixSumming','Number of rows in matrix: '//TRIM(I2S(n)),Level=10)

    A % FORMAT = MATRIX_CRS
    
    A % NumberOfRows = n
    ALLOCATE( A % Cols(n), A % Rows(n+1), A % Values(n) )

    A % Cols = 0
    A % Rows = 0
    A % Values = 0.0_dp
    
    minsum = HUGE(minsum)
    maxsum = 0.0_dp
    sumsum = 0.0_dp
    
    DO i = 1, n
      rowsum = 0.0_dp
      DO j=P % Rows(i),P % Rows(i+1)-1
        k = P % Cols(j)
        val = P % Values(j) 
        rowsum = rowsum + ABS( val )
      END DO

      A % Values(i) = rowsum
      A % Cols(i) = i
      A % Rows(i) = i

      minsum = MIN(minsum, rowsum)
      maxsum = MAX(maxsum, rowsum)
      sumsum = sumsum + rowsum
    END DO
    A % Rows(n+1) = n+1

    PRINT *,'diagonal sums:',minsum,maxsum,sumsum/n
    
  END SUBROUTINE DiagonalMatrixSumming



!------------------------------------------------------------------------------
!> Assemble coupling matrix related to fluid-structure interaction
!------------------------------------------------------------------------------
  SUBROUTINE FsiCouplingAssembly( Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
      IsPlate, IsShell, IsNS )
    
    TYPE(Solver_t) :: Solver          ! leading solver
    TYPE(Variable_t), POINTER :: FVar ! fluid variable
    TYPE(Variable_t), POINTER :: SVar ! structure variable
    TYPE(Matrix_t), POINTER :: A_fs, A_sf, A_f, A_s
    LOGICAL :: IsPlate, IsShell, IsNS
   !------------------------------------------------------------------------------
    LOGICAL, POINTER :: ConstrainedF(:), ConstrainedS(:)
    INTEGER, POINTER :: FPerm(:), SPerm(:)
    INTEGER :: FDofs, SDofs
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Indexes(:), pIndexes(:)
    INTEGER :: i,j,ii,jj,k,n,t,istat,pn,ifluid,jstruct,pcomp
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Solver_t), POINTER :: PSolver
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp), ALLOCATABLE :: MASS(:,:)
    REAL(KIND=dp), POINTER :: Basis(:)
    REAL(KIND=dp) :: detJ, val, c(3), pc(3), Normal(3), coeff, Omega, Rho, area, fdiag
    LOGICAL :: Stat, IsHarmonic
    INTEGER :: dim,mat_id,tcount
    LOGICAL :: FreeF, FreeS, FreeFim, FreeSim, UseDensity, Found
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    REAL(KIND=dp) :: MultSF, MultFS
    CHARACTER(*), PARAMETER :: Caller = 'FsiCouplingAssembly'
   
    
    CALL Info(Caller,'Creating coupling matrix for harmonic FSI',Level=6)

    
    IF( A_fs % FORMAT /= MATRIX_LIST ) THEN
      A_fs % Values = 0.0_dp
      A_sf % Values = 0.0_dp      
    END IF
      
    
    Mesh => Solver % Mesh
    FPerm => FVar % Perm
    SPerm => SVar % Perm
    
    fdofs = FVar % Dofs
    sdofs = SVar % Dofs

    IF( IsPlate ) CALL Info(Caller,'Assuming structure to be plate',Level=8)

    IF( IsShell ) CALL Info(Caller,'Assuming structure to be shell',Level=8)

    IF( IsNS ) CALL Info(Caller,'Assuming fluid to have velocities',Level=8)


    UseDensity = .FALSE.
    DO i=1,CurrentModel % NumberOfSolvers
      PSolver => CurrentModel % Solvers(i)      
      IF( ASSOCIATED( PSolver % Variable, FVar ) ) THEN
        UseDensity = ListGetLogical( PSolver % Values,'Use Density',Found ) 
        EXIT
      END IF
    END DO
    IF( UseDensity ) THEN
      CALL Info(Caller,'The Helmholtz equation is multiplied by density',Level=10)
    END IF
    
    
    ConstrainedF => A_f % ConstrainedDof
    ConstrainedS => A_s % ConstrainedDof
    
    
    ! Here we assume harmonic coupling if there are more then 3 structure dofs
    dim = 3
    IsHarmonic = .FALSE.
    IF( IsPlate ) THEN
      IF( sdofs == 6 ) THEN
        IsHarmonic = .TRUE.
      ELSE IF( sdofs /= 3 ) THEN
        CALL Fatal(Caller,'Invalid number of dofs in plate solver: '//TRIM(I2S(sdofs)))
      END IF
    ELSE IF( IsShell ) THEN
      IF( sdofs == 12 ) THEN
        IsHarmonic = .TRUE.
      ELSE IF( sdofs /= 6 ) THEN
        CALL Fatal(Caller,'Invalid number of dofs in shell solver: '//TRIM(I2S(sdofs)))
      END IF
    ELSE
      IF( sdofs == 4 .OR. sdofs == 6 ) THEN
        IsHarmonic = .TRUE.
      ELSE IF( sdofs /= 2 .AND. sdofs /= 3 ) THEN
        CALL Fatal(Caller,'Invalid number of dofs in elasticity solver: '//TRIM(I2S(sdofs)))
      END IF
      IF( sdofs == 4 .OR. sdofs == 2 ) dim = 2
    END IF

    ! The elasticity solver defines whether the system is real or harmonic
    IF( IsHarmonic ) THEN
      CALL Info(Caller,'Assuming harmonic coupling matrix',Level=10)
    ELSE
      CALL Info(Caller,'Assuming real valued coupling matrix',Level=10)
    END IF

    
    ! The fluid system must be consistent with elasticity system
    IF( IsNS ) THEN
      IF( IsHarmonic ) THEN
        IF( fdofs /= 2*(dim+2) .AND. fdofs /= 2*(dim+1) ) THEN
          CALL Fatal(Caller,&
              'Inconsistent number of harmonic dofs in NS solver: '//TRIM(I2S(fdofs)))
        END IF
        ! pressure component
        pcomp = fdofs / 2
      ELSE
        IF( fdofs /= (dim+2) .AND. fdofs /= (dim+1) ) THEN
          CALL Fatal(Caller,&
              'Inconsistent number of real dofs in NS solver: '//TRIM(I2S(fdofs)))
        END IF
        pcomp = fdofs
      END IF
      ALLOCATE( NodeDone(MAXVAL(FPerm)) )
      NodeDone = .FALSE.
    ELSE
      IF( IsHarmonic ) THEN
        IF( fdofs /= 2 ) CALL Fatal(Caller,&
            'Inconsistent number of harmonic dofs in pressure solver: '//TRIM(I2S(fdofs)))
      ELSE
        IF( fdofs /= 1 ) CALL Fatal(Caller,&
            'Inconsistent number of real dofs in pressure solver: '//TRIM(I2S(fdofs)))
      END IF
      pcomp = 1
    END IF

    
    IF( IsHarmonic ) THEN
      Omega = 2 * PI * ListGetCReal( CurrentModel % Simulation,'Frequency',Stat ) 
      IF( .NOT. Stat) THEN
        CALL Fatal(Caller,'Frequency in Simulation list not found!')
      END IF
    ELSE
      Omega = 0.0_dp
    END IF
    
    i = SIZE( FVar % Values ) 
    j = SIZE( SVar % Values ) 
    
    CALL Info(Caller,'Fluid dofs '//TRIM(I2S(i))//&
        ' with '//TRIM(I2S(fdofs))//' components',Level=10)
    CALL Info(Caller,'Structure dofs '//TRIM(I2S(j))//&
        ' with '//TRIM(I2S(sdofs))//' components',Level=10)   
    CALL Info(Caller,'Assuming '//TRIM(I2S(dim))//&
        ' active dimensions',Level=10)   

    ! Add the lasrgest entry that allocates the whole list matrix structure
    CALL AddToMatrixElement(A_fs,i,j,0.0_dp)
    CALL AddToMatrixElement(A_sf,j,i,0.0_dp)
    
    N = Mesh % MaxElementNodes 
    ALLOCATE( Basis(n), MASS(N,N), Nodes % x(n), Nodes % y(n), Nodes % z(n), &
        STAT=istat)

    tcount = 0
    area = 0.0_dp
    
    MultFS = ListGetCReal( Solver % Values,'FS multiplier',Found)
    IF( .NOT. Found ) MultFS = 1.0_dp

    MultSF = ListGetCReal( Solver % Values,'SF multiplier',Found)
    IF( .NOT. Found ) MultSF = 1.0_dp
    
    FreeS = .TRUE.; FreeSim = .TRUE.
    FreeF = .TRUE.; FreeFim = .TRUE.    
    
    
    DO t=Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      
      Element => Mesh % Elements(t)
      n = Element % TYPE % NumberOfNodes
      Indexes => Element % NodeIndexes
      
      IF( ANY( FPerm(Indexes) == 0 ) ) CYCLE
      IF( ANY( SPerm(Indexes) == 0 ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      
      Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
      Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
      Nodes % z(1:n) = Mesh % Nodes % z(Indexes)
      
      Normal = NormalVector( Element, Nodes )

    
      ! The following is done in order to check that the normal points to the fluid      
      Parent => Element % BoundaryInfo % Left
      IF( ASSOCIATED( Parent ) ) THEN
        IF( ANY( FPerm(Parent % NodeIndexes) == 0 ) ) Parent => NULL()
      END IF

      IF(.NOT. ASSOCIATED( Parent ) ) THEN
        Parent => Element % BoundaryInfo % Right
        IF( ASSOCIATED( Parent ) ) THEN
          IF( ANY( FPerm(Parent % NodeIndexes) == 0 ) ) Parent => NULL()
        END IF
      END IF
                
      ! Could not find a proper fluid element to define the normal 
      IF(.NOT. ASSOCIATED( Parent ) ) CYCLE

      tcount = tcount + 1

      
      pn = Parent % TYPE % NumberOfNodes
      pIndexes => Parent % NodeIndexes
      
      c(1) =  SUM( Nodes % x(1:n) ) / n
      c(2) =  SUM( Nodes % y(1:n) ) / n
      c(3) =  SUM( Nodes % z(1:n) ) / n
      
      pc(1) =  SUM( Mesh % Nodes % x(pIndexes) ) / pn
      pc(2) =  SUM( Mesh % Nodes % y(pIndexes) ) / pn
      pc(3) =  SUM( Mesh % Nodes % z(pIndexes) ) / pn
      
      ! The normal vector has negative projection to the vector drawn from center of
      ! boundary element to the center of bulk element. 
      IF( SUM( (pc-c) * Normal ) < 0.0_dp ) THEN
        Normal = -Normal
      END IF
      
      MASS(1:n,1:n) = 0.0_dp
      
      mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,'Material' )
      rho = ListGetConstReal( CurrentModel % Materials(mat_id) % Values,'Density',Stat)
      IF(.NOT. Stat) rho = ListGetConstReal( CurrentModel % Materials(mat_id) % Values, &
          'Equilibrium Density',Stat)

      IF( .NOT. Stat) THEN
        CALL Fatal(Caller,'Fluid density not found in material :'//TRIM(I2S(mat_id)))
      END IF
      
      ! The sign depends on the convection of the normal direction
      ! If density is divided out already in the Helmholtz equation the multiplier will
      ! be different. 
      IF( UseDensity ) THEN
        coeff = omega**2
      ELSE
        coeff = rho * omega**2
      END IF
      
      ! Numerical integration:
      !----------------------
      IP = GaussPoints( Element )
      
      DO k=1,IP % n        
        ! Basis function values & derivatives at the integration point:
        !--------------------------------------------------------------
        stat = ElementInfo( Element, Nodes, IP % U(k), IP % V(k), &
            IP % W(k),  detJ, Basis )
        
        ! The mass matrix of the boundary element
        !----------------------------------------
        val = IP % s(k) * detJ
        DO i=1,n
          DO j=1,n
            MASS(i,j) = MASS(i,j) + val * Basis(i) * Basis(j)
          END DO
        END DO
        area = area + val
      END DO

      ! A: fs
      ! Effect of structure on fluid           
      IF( IsNs ) THEN
        ! For the N-S equation the condition applies directly on the velocity components
        
        DO i=1,n
          ii = Indexes(i)
          j = i
          jj = Indexes(j) ! one-to-one mapping


          IF( NodeDone( Fperm(ii) ) ) CYCLE
          NodeDone( FPerm(ii) ) = .TRUE.
          
          
          DO k=1,dim
            
            ! The velocity component of the fluid
            IF( IsHarmonic ) THEN
              ifluid = fdofs*(FPerm(ii)-1)+2*(k-1)+1
              !IF( ASSOCIATED( ConstrainedF ) ) THEN
              !  FreeF = .NOT. ConstrainedF(ifluid)
              !  FreeFim = .NOT. ConstrainedF(ifluid+1)            
              !END IF
            ELSE
              ifluid = fdofs*(FPerm(ii)-1)+k
              !IF( ASSOCIATED( ConstrainedF ) ) THEN
              !  FreeF = .NOT. ConstrainedF(ifluid)
              !END IF
            END IF
                          
            ! Shell and 3D elasticity are both treated with the same routine
            IF( .NOT. IsPlate ) THEN

              IF( IsHarmonic ) THEN
                val = omega
                jstruct = sdofs*(SPerm(jj)-1)+2*(k-1)+1  
              ELSE
                CALL Fatal(Caller,'NS coupling only done for harmonic system!')               
              END IF

                
            ELSE ! If IsPlate
              IF( IsHarmonic ) THEN              
                val = omega * Normal(k)
                
                ! By default the plate should be oriented so that normal points to z
                ! If there is a plate then fluid is always 3D
                IF( Normal(3) < 0 ) val = -val

                jstruct = sdofs*(SPerm(jj)-1)+1
              ELSE
                CALL Fatal(Caller,'NS coupling only done for harmonic system!')               
              END IF
            END IF

            IF( IsHarmonic ) THEN                                             
              ! Structure load on the fluid: v = i*omega*u
              fdiag = A_f % Values( A_f % diag(ifluid) )
              IF( FreeF ) THEN
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,MultFS*val*fdiag)     ! Re 
              ELSE
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,0.0_dp)
              END IF
              
              fdiag = A_f % Values( A_f % diag(ifluid+1) )
              IF( FreeFim ) THEN
                CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,-MultFS*val*fdiag)      ! Im
              ELSE                
                CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp )
              END IF

              ! These must be created for completeness because the matrix topology of complex
              ! matrices must be the same for all components.
              CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
              CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp)
            ELSE
              CALL Fatal(Caller,'NS coupling only done for harmonic system!')
            END IF
          END DO
        END DO

      ELSE ! .NOT. IsNS
        ! For pressure equations (Helmholtz) the structure applies a Neumann condition
        
        DO i=1,n
          ii = Indexes(i)

          ! The pressure component of the fluid
          IF( IsHarmonic ) THEN
            ifluid = fdofs*(FPerm(ii)-1)+2*(pcomp-1)+1
            IF( ASSOCIATED( ConstrainedF ) ) THEN
              FreeF = .NOT. ConstrainedF(ifluid)
              FreeFim = .NOT. ConstrainedF(ifluid+1)            
            END IF
          ELSE
            ifluid = fdofs*(FPerm(ii)-1)+pcomp
            IF( ASSOCIATED( ConstrainedF ) ) THEN
              FreeF = .NOT. ConstrainedF(ifluid)
            END IF
          END IF


          DO j=1,n
            jj = Indexes(j)

            ! Shell and 3D elasticity are both treated with the same routine
            IF( .NOT. IsPlate ) THEN

              DO k=1,dim

                val = MASS(i,j) * Normal(k)

                IF( IsHarmonic ) THEN
                  jstruct = sdofs*(SPerm(jj)-1)+2*(k-1)+1  

                  ! Structure load on the fluid: This assembles
                  !
                  !    -1/rho <dp/dn,v> = -omega^2 <u.n,v> = omega^2 <u.m,v> 
                  !
                  ! with the normal vectors satisfying m = -n. Note that the density (rho) 
                  ! must be defined for Helmholtz solver to make it assemble a system
                  ! consistent with the boundary integral -1/rho <dp/dn,v>.
                  IF( FreeF ) THEN
                    CALL AddToMatrixElement(A_fs,ifluid,jstruct,MultFS*val*coeff)     ! Re 
                  ELSE
                    CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
                  END IF

                  IF( FreeFim ) THEN
                    CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,MultFS*val*coeff) ! Im
                  ELSE                
                    CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp )
                  END IF

                  ! These must be created for completeness because the matrix topology of complex
                  ! matrices must be the same for all components.
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,0.0_dp)     
                  CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,0.0_dp)
                ELSE
                  jstruct = sdofs*(SPerm(jj)-1)+k

                  ! Structure load on the fluid: dp/dn = -u. (This seems strange???)
                  IF( FreeF ) THEN
                    CALL AddToMatrixElement(A_fs,ifluid,jstruct,-MultFS*val)           
                  END IF
                END IF
              END DO

            ELSE ! If IsPlate

              val = MASS(i,j) 

              ! By default the plate should be oriented so that normal points to z
              ! If there is a plate then fluid is always 3D
              IF( Normal(3) < 0 ) val = -val

              IF( IsHarmonic ) THEN
                jstruct = sdofs*(SPerm(jj)-1)+1

                ! Structure load on the fluid: -1/rho dp/dn = -omega^2 u.n = omega^2 u.m
                IF( FreeF ) THEN
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct,MultFS*val*coeff)     ! Re 
                ELSE
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct,0.0_dp)
                END IF

                IF( FreeFim ) THEN
                  CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,MultFS*val*coeff) ! Im
                ELSE                
                  CALL AddToMatrixElement(A_fs,ifluid+1,jstruct+1,0.0_dp )
                END IF

                ! These must be created for completeness because the matrix topology of complex
                ! matrices must be the same for all components.
                CALL AddToMatrixElement(A_fs,ifluid,jstruct+1,0.0_dp)
                CALL AddToMatrixElement(A_fs,ifluid+1,jstruct,0.0_dp)
              ELSE
                jstruct = sdofs*(SPerm(jj)-1)+1

                ! Structure load on the fluid: dp/dn = -u. (This seems strange???)
                IF( FreeF ) THEN
                  CALL AddToMatrixElement(A_fs,ifluid,jstruct,-MultFS*val)           
                END IF
              END IF

            END IF

          END DO
        END DO
      END IF
        

      ! A_sf:
      ! Effect of fluid (pressure) on structure.
      ! Each component get the normal component of the pressure as a r.h.s. term.
      ! The plate equation just gets the full load and is treated separately. 
      !----------------------------------------------------------------------------
      DO i=1,n
        ii = Indexes(i)

        ! The pressure component of the fluid
        IF( IsHarmonic ) THEN
          ifluid = fdofs*(FPerm(ii)-1)+2*(pcomp-1)+1
        ELSE
          ifluid = fdofs*(FPerm(ii)-1)+pcomp
        END IF

        DO j=1,n
          jj = Indexes(j)
          
          ! Shell and 3D elasticity are both treated with the same routine
          IF( .NOT. IsPlate ) THEN

            DO k=1,dim
              
              val = MASS(i,j) * Normal(k)
              
              IF( IsHarmonic ) THEN
                jstruct = sdofs*(SPerm(jj)-1)+2*(k-1)+1  
                
                IF( ASSOCIATED( ConstrainedS ) ) THEN
                  FreeS = .NOT. ConstrainedS(jstruct)
                  FreeSim = .NOT. ConstrainedS(jstruct+1)
                END IF

                ! Fluid load on the structure: tau \cdot n = p * n
                IF( FreeS ) THEN
                  CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           ! Re terms coupling
                ELSE
                  CALL AddToMatrixElement(A_sf,jstruct,ifluid,0.0_dp)
                END IF
                
                IF( FreeSim ) THEN
                  CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,MultSF*val)       ! Im
                ELSE
                  CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,0.0_dp)
                END IF

                ! These must be created for completeness because the matrix topology of complex
                ! matrices must be the same for all components.
                CALL AddToMatrixElement(A_sf,jstruct,ifluid+1,0.0_dp)
                CALL AddToMatrixElement(A_sf,jstruct+1,ifluid,0.0_dp)
              ELSE
                jstruct = sdofs*(SPerm(jj)-1)+k
                
                IF( ASSOCIATED( ConstrainedS ) ) THEN
                  FreeS = .NOT. ConstrainedS(jstruct)
                END IF

                ! Fluid load on the structure: tau \cdot n = p * n
                IF( FreeS ) THEN
                  CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           
                END IF
                
              END IF
            END DO
            
          ELSE ! If IsPlate
            
            val = MASS(i,j) 
            
            ! By default the plate should be oriented so that normal points to z
            ! If there is a plate then fluid is always 3D
            IF( Normal(3) < 0 ) val = -val
            
            IF( IsHarmonic ) THEN
              jstruct = sdofs*(SPerm(jj)-1)+1
              
              IF( ASSOCIATED( ConstrainedS ) ) THEN
                FreeS = .NOT. ConstrainedS(jstruct)
                FreeSim = .NOT. ConstrainedS(jstruct+1)
              END IF
              
              ! Fluid load on the structure: tau \cdot n = p * n
              IF( FreeS ) THEN
                CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           ! Re terms coupling
              ELSE
                CALL AddToMatrixElement(A_sf,jstruct,ifluid,0.0_dp)
              END IF

              IF( FreeSim ) THEN
                CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,MultSF*val)       ! Im
              ELSE
                CALL AddToMatrixElement(A_sf,jstruct+1,ifluid+1,0.0_dp)
              END IF

              ! These must be created for completeness because the matrix topology of complex
              ! matrices must be the same for all components.
              CALL AddToMatrixElement(A_sf,jstruct,ifluid+1,0.0_dp)
              CALL AddToMatrixElement(A_sf,jstruct+1,ifluid,0.0_dp)
            ELSE
              jstruct = sdofs*(SPerm(jj)-1)+1

              IF( ASSOCIATED( ConstrainedS ) ) THEN
                FreeS = .NOT. ConstrainedS(jstruct)
              END IF
              
              ! Fluid load on the structure: tau \cdot n = p * n
              IF( FreeS ) THEN
                CALL AddToMatrixElement(A_sf,jstruct,ifluid,MultSF*val)           
              END IF
            END IF
            
          END IF

        END DO
      END DO

    END DO ! Loop over boundary elements
    
      
    DEALLOCATE( Basis, MASS, Nodes % x, Nodes % y, Nodes % z)

    IF( A_fs % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(A_fs)
      CALL List_toCRSMatrix(A_sf)
    END IF
      
    !PRINT *,'interface area:',area
    !PRINT *,'interface fs sum:',SUM(A_fs % Values), SUM( ABS( A_fs % Values ) )
    !PRINT *,'interface sf sum:',SUM(A_sf % Values), SUM( ABS( A_sf % Values ) )

    CALL Info(Caller,'Number of elements on interface: '&
        //TRIM(I2S(tcount)),Level=10)    
    CALL Info(Caller,'Number of entries in fluid-structure matrix: '&
        //TRIM(I2S(SIZE(A_fs % Values))),Level=10)
    CALL Info(Caller,'Number of entries in structure-fluid matrix: '&
        //TRIM(I2S(SIZE(A_sf % Values))),Level=10)
    
    CALL Info(Caller,'All done',Level=20)

    
  END SUBROUTINE FsiCouplingAssembly





  
  ! The following function is a copy from ShellSolver.F90.
  ! The suffix Int is added for unique naming.
  !-------------------------------------------------------------------------------
  FUNCTION GetElementalDirectorInt(Mesh, Element, &
      ElementNodes, node) RESULT(DirectorValues) 
  !-------------------------------------------------------------------------------    
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER, INTENT(IN) :: Element
    TYPE(Nodes_t), OPTIONAL, INTENT(IN) :: ElementNodes
    INTEGER, OPTIONAL :: node
    REAL(KIND=dp), POINTER :: DirectorValues(:)
    !-------------------------------------------------------------------------------
    TYPE(Nodes_t) :: Nodes
    LOGICAL :: Visited = .FALSE., UseElementProperty = .FALSE., UseNormalSolver = .FALSE.
    REAL(KIND=dp) :: Normal(3)
    REAL(KIND=dp), POINTER :: NodalNormals(:)
    REAL(KIND=dp), POINTER :: DirectorAtNode(:)
    REAL(KIND=dp), POINTER :: PropertyValues(:)
    INTEGER :: n    

    SAVE Visited, UseElementProperty, NodalNormals, DirectorAtNode, Nodes
    !-------------------------------------------------------------------------------
        
    IF (.NOT. Visited) THEN
      DirectorValues => GetElementPropertyInt('director', Element)
      UseElementProperty = ASSOCIATED( DirectorValues ) 

      IF (.NOT. UseElementProperty) THEN
        n = CurrentModel % MaxElementNodes
        ALLOCATE( NodalNormals(3*n), Nodes % x(n), Nodes % y(n), Nodes % z(n) ) 
      END IF
      ALLOCATE( DirectorAtNode(3) )
      Visited = .TRUE.
    END IF

    IF ( UseElementProperty ) THEN    
      PropertyValues => GetElementPropertyInt('director', Element)
      IF( PRESENT( node ) ) THEN        
        DirectorAtNode(1:3) = PropertyValues(3*(node-1)+1:3*node)
        DirectorValues => DirectorAtNode
      ELSE
        DirectorValues => PropertyValues
      END IF
      
    ELSE
      IF( PRESENT( ElementNodes ) ) THEN
        Normal = NormalVector( Element, ElementNodes, Check = .TRUE. ) 
      ELSE
        Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
        Normal = NormalVector( Element, Nodes, Check = .TRUE. ) 
      END IF

      IF( PRESENT( node ) ) THEN
        !PRINT *,'Normal:',Normal
        DirectorAtNode(1:3) = Normal(1:3)
        DirectorValues => DirectorAtNode
      ELSE              
        n = Element % TYPE % NumberOfNodes
        NodalNormals(1:3*n:3) = Normal(1)
        NodalNormals(2:3*n:3) = Normal(2)
        NodalNormals(3:3*n:3) = Normal(3)      
        DirectorValues => NodalNormals
      END IF
    END IF

  CONTAINS
        
    FUNCTION GetElementPropertyInt( Name, Element ) RESULT(Values)
      CHARACTER(LEN=*) :: Name
      TYPE(Element_t), POINTER :: Element
      REAL(KIND=dp), POINTER :: Values(:)

      TYPE(ElementData_t), POINTER :: p

      Values => NULL()
      p=> Element % PropertyData

      DO WHILE( ASSOCIATED(p) )
        IF ( Name==p % Name ) THEN
          Values => p % Values
          RETURN
        END IF
        p => p % Next
      END DO
    END FUNCTION GetElementPropertyInt
    
  !-------------------------------------------------------------------------------    
  END FUNCTION GetElementalDirectorInt
  !-------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Assemble coupling matrices related to structure-structure interaction.
!> A possible scenario is that the diagonal blocks are the matrices of the 
!> solvers listed using the keyword "Block Solvers". The (1,1)-block is then
!> tied up with the value of the first entry in the "Block Solvers" array. 
!------------------------------------------------------------------------------
  SUBROUTINE StructureCouplingAssembly(Solver, FVar, SVar, A_f, A_s, A_fs, A_sf, &
      IsSolid, IsPlate, IsShell, IsBeam )
!------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver          !< The leading solver defining block structure 
    TYPE(Variable_t), POINTER :: FVar !< "Slave" structure variable
    TYPE(Variable_t), POINTER :: SVar !< "Master" structure variable
    TYPE(Matrix_t), POINTER :: A_f    !< (2,2)-block for the "slave" variable
    TYPE(Matrix_t), POINTER :: A_s    !< (1,1)-block for the "master" variable
    TYPE(Matrix_t), POINTER :: A_fs   !< (2,1)-block for interaction
    TYPE(Matrix_t), POINTER :: A_sf   !< (1,2)-block for interaction
    LOGICAL :: IsSolid, IsPlate, IsShell, IsBeam !< The type of the slave variable
   !------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, POINTER :: ConstrainedF(:), ConstrainedS(:)
    LOGICAL :: DoDamp, DoMass
    INTEGER, POINTER :: FPerm(:), SPerm(:)
    INTEGER :: FDofs, SDofs
    INTEGER :: i,j,k,jf,js,kf,ks,nf,ns,dim,ncount
    REAL(KIND=dp) :: vdiag
    CHARACTER(*), PARAMETER :: Caller = 'StructureCouplingAssembly'
   !------------------------------------------------------------------------------

    CALL Info(Caller,'Creating coupling matrix for structures',Level=6)
    
    Mesh => Solver % Mesh
    dim = Mesh % MeshDim

    ! S refers to the first and F to the second block (was fluid):
    FPerm => FVar % Perm
    SPerm => SVar % Perm
    
    fdofs = FVar % Dofs
    sdofs = SVar % Dofs

    IF( IsSolid ) CALL Info(Caller,'Assuming coupling with solid solver',Level=8)
    IF( IsBeam )  CALL Info(Caller,'Assuming coupling with beam solver',Level=8)
    IF( IsPlate ) CALL Info(Caller,'Assuming coupling with plate solver',Level=8)
    IF( IsShell ) CALL Info(Caller,'Assuming coupling with shell solver',Level=8)
    
    ConstrainedF => A_f % ConstrainedDof
    ConstrainedS => A_s % ConstrainedDof
                  
    nf = SIZE( FVar % Values ) 
    ns = SIZE( SVar % Values ) 
    
    CALL Info(Caller,'Slave structure dofs '//TRIM(I2S(nf))//&
        ' with '//TRIM(I2S(fdofs))//' components',Level=10)
    CALL Info(Caller,'Master structure dofs '//TRIM(I2S(ns))//&
        ' with '//TRIM(I2S(sdofs))//' components',Level=10)   
    CALL Info(Caller,'Assuming '//TRIM(I2S(dim))//&
        ' active dimensions',Level=10)   

    IF( A_fs % FORMAT == MATRIX_LIST ) THEN
      ! Add the largest entry that allocates the whole list matrix structure
      CALL AddToMatrixElement(A_fs,nf,ns,0.0_dp)
      CALL AddToMatrixElement(A_sf,ns,nf,0.0_dp)
    ELSE
      ! If we are revisiting then initialize the CRS matrices to zero
      A_fs % Values = 0.0_dp
      A_sf % Values = 0.0_dp      
    END IF

    DoMass = .FALSE.
    IF( ASSOCIATED( A_f % MassValues ) ) THEN
      IF( ASSOCIATED( A_s % MassValues ) ) THEN
        DoMass = .TRUE.        
      ELSE
        CALL Warn(Caller,'Both solid and shell should have MassValues!')
      END IF
    END IF

    DoDamp = ASSOCIATED( A_f % DampValues )
    IF( DoDamp ) THEN
      CALL Warn(Caller,'Damping matrix values at shell interface will be dropped!')
    END IF

    ! This is still under development and not used for anything
    ! Probably this will not be needed at all but rather we need the director.
    !IF( IsShell ) CALL DetermineCouplingNormals()

    ! For the shell equation enforce the directional derivative of the displacement
    ! field in implicit manner from solid displacements. The interaction conditions
    ! for the corresponding forces are also created.
    IF (IsShell) THEN
      BLOCK
        INTEGER, POINTER :: Indexes(:)
        INTEGER, ALLOCATABLE :: NodeHits(:), InterfacePerm(:), InterfaceElems(:,:)
        INTEGER :: InterfaceN, hits
        INTEGER :: p,lf,ls,ii,jj,n,m,t
        REAL(KIND=dp), POINTER :: Director(:)
        REAL(KIND=dp), POINTER :: Basis(:), dBasisdx(:,:)
        REAL(KIND=dp), ALLOCATABLE :: A_f0(:)
        REAL(KIND=dp) :: u,v,w,weight,detJ,val
        REAL(KIND=dp) :: x, y, z 

        TYPE(Element_t), POINTER :: Element, ShellElement
        TYPE(Nodes_t) :: Nodes
        LOGICAL :: Stat

        n = Mesh % MaxElementNodes 
        ALLOCATE( Basis(n), dBasisdx(n,3), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

        ! Memorize the original values
        ALLOCATE( A_f0( SIZE( A_f % Values ) ) )
        A_f0 = A_f % Values

        ALLOCATE( NodeHits( Mesh % NumberOfNodes ), InterfacePerm( Mesh % NumberOfNodes ) )
        NodeHits = 0
        InterfacePerm = 0

        ! First, zero the rows related to directional derivative dofs, 
        ! i.e. the components 4,5,6. "s" refers to solid and "f" to shell.
        InterfaceN = 0
        DO i=1,Mesh % NumberOfNodes
          jf = FPerm(i)      
          js = SPerm(i)
          IF( jf == 0 .OR. js == 0 ) CYCLE

          ! also number the interface
          InterfaceN = InterfaceN + 1
          InterfacePerm(i) = InterfaceN

          DO lf = 4, 6
            kf = fdofs*(jf-1)+lf

            IF( ConstrainedF(kf) ) CYCLE

            DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
              A_f % Values(k) = 0.0_dp
              IF (DoMass) THEN
                A_f % MassValues(k) = 0.0_dp
              END IF
              IF( DoDamp) THEN
                A_f % DampValues(k) = 0.0_dp
              END IF
            END DO
            A_f % rhs(kf) = 0.0_dp
          END DO
        END DO

        CALL Info(Caller,'Number of nodes at interface: '//TRIM(I2S(InterfaceN)),Level=12)

        ALLOCATE( InterfaceElems(InterfaceN,2) )
        InterfaceElems = 0

        ! Then go through shell elements and associate each interface node with a list of
        ! shell elements defined in terms of the interface node
        DO t=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(t)
          Indexes => Element % NodeIndexes 

          n = Element % TYPE % ElementCode
          IF( n > 500 .OR. n < 300 ) CYCLE

          ! We must have shell equation present everywhere and solid equation at least in two nodes
          IF(ANY( FPerm(Indexes) == 0 ) ) CYCLE 
          k = COUNT( SPerm(Indexes) > 0 )
          IF( k < 2 ) CYCLE

          n = Element % Type % NumberOfNodes            
          DO i=1,n
            j = Indexes(i)
            k = InterfacePerm(j)
            IF( k == 0) CYCLE

            ! Assuming just two shell parents currently
            IF( InterfaceElems(k,1) == 0 ) THEN
              InterfaceElems(k,1) = t
            ELSE IF( InterfaceElems(k,2) == 0 ) THEN
              InterfaceElems(k,2) = t
            ELSE
              CALL Fatal(Caller,'Tree interface elems?')
            END IF
          END DO
        END DO


        ! Then go through solid elements associated with the interface and count
        ! how many solid elements share each interface node.
        NodeHits = 0
        DO t=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(t)
          Indexes => Element % NodeIndexes 

          ! We must have solid equation present everywhere and shell at least at one node.
          IF(ANY( SPerm(Indexes) == 0 ) ) CYCLE
          IF(ALL( FPerm(Indexes) == 0 ) ) CYCLE

          n = COUNT( FPerm(Indexes) > 0 )
          IF( n /= 2 ) THEN
            CALL Fatal(Caller,'Currently we can only deal with two hits!')
          END IF

          DO i=1,SIZE(Indexes)
            j = FPerm(Indexes(i))
            IF( j == 0 ) CYCLE
            NodeHits(j) = NodeHits(j) + 1              
          END DO
        END DO

        !PRINT *,'Maximum node hits:',MAXVAL(NodeHits)

        ! Then go through each solid element associated with the interface and
        ! create matrix entries defining the interaction conditions for the
        ! directional derivatives and corresponding forces. 
        DO t=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(t)
          Indexes => Element % NodeIndexes 

          ! We must have solid equation present everywhere and shell at least at one node.
          IF(ANY( SPerm(Indexes) == 0 ) ) CYCLE
          IF(ALL( FPerm(Indexes) == 0 ) ) CYCLE

          n = Element % TYPE % NumberOfNodes
          Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
          Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
          Nodes % z(1:n) = Mesh % Nodes % z(Indexes)

          x = 0.0_dp; y = 0.0_dp; z = 0.0_dp
          DO i=1,n
            IF( FPerm(Indexes(i)) == 0 ) CYCLE
            x = x + Nodes % x(i)
            y = y + Nodes % y(i)
            z = z + Nodes % z(i)
          END DO
          x = x/2; y = y/2; z = z/2;

          ! TO DO: The following call may not work for p-elements!
          CALL GlobalToLocal( u, v, w, x, y, z, Element, Nodes )

          ! Integration at the center of the edge
          stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis, dBasisdx )          

          DO ii = 1, n
            i = Indexes(ii)           
            jf = FPerm(i)      
            IF( jf == 0 ) CYCLE

            Weight = 1.0_dp / NodeHits(jf)

            !PRINT *,'Weight:',Weight

            DO j=1,2
              ShellElement => Mesh % Elements(InterfaceElems(InterfacePerm(i),j))
              hits = 0
              DO k=1,ShellElement % TYPE % NumberOfNodes
                IF( ANY( Indexes == ShellElement % NodeIndexes(k) ) ) hits = hits + 1
              END DO
              IF( hits >= 2 ) EXIT
            END DO

            ! Retrieve the director of the shell element: 
            m = ShellElement % TYPE % NumberOfNodes
            DO jj=1,m
              IF(Element % NodeIndexes(ii) == ShellElement % NodeIndexes(jj) ) EXIT
            END DO
            Director => GetElementalDirectorInt(Mesh,ShellElement,node=jj)

            !PRINT *,'Director:',ShellElement % ElementIndex,jj,Director            

            DO lf = 4, 6                            
              kf = fdofs*(jf-1)+lf

              IF( ConstrainedF(kf) ) CYCLE

              DO ls = 1, dim
                !
                ! Directional derivative dofs of the shell equations: 
                ! We try to enforce the condition d_{i+3}=-<(grad u)n,e_i> 
                ! where i=1,2,3; i+3=lf, n is director, e_i is unit vector, and 
                ! u is the displacement field of the solid. 
                DO p=1,n
                  js = SPerm(Indexes(p))
                  ks = sdofs*(js-1)+lf-3
                  val = Director(ls) * dBasisdx(p,ls)

                  CALL AddToMatrixElement(A_fs,kf,ks,weight*val)

                  ! Here the idea is to distribute the implicit moments of the shell solver
                  ! to forces for the solid solver. So even though the stiffness matrix related to the
                  ! directional derivatives is nullified, the forces are not forgotten.
                  ! This part may be thought of as being based on two (Råback's) conjectures: 
                  ! in the first place the Lagrange variable formulation should bring us to a symmetric 
                  ! coefficient matrix and the values of Lagrange variables can be estimated as nodal 
                  ! reactions obtained by performing a matrix-vector product.
                  !
                  ! Note that no attempt is currently made to transfer external moment
                  ! loads of the shell model to loads of the coupled model. Likewise
                  ! rotational inertia terms of the shell model are not transformed
                  ! to inertia terms of the coupled model. Neglecting the rotational
                  ! inertia might be acceptable in many cases.
                  !
                  ! Note that the minus sign of the entries is correct here:
                  DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
                    CALL AddToMatrixElement(A_sf,ks,A_f % Cols(k),-weight*val*A_f0(k)) 
                  END DO
                END DO
              END DO

              ! This should sum up to unity!
              CALL AddToMatrixElement(A_f,kf,kf,weight)
            END DO
          END DO
        END DO
        DEALLOCATE( Basis, dBasisdx, Nodes % x, Nodes % y, Nodes % z )
        DEALLOCATE(A_f0, NodeHits, InterfacePerm)
      END BLOCK
    END IF
    
    ! Note: we may have to recheck this coupling if visiting for 2nd time!
    !
    ! Three DOFs for both shells and solids are the real Cartesian components of
    ! the displacement. Hence we can deal with the common parts of solid-solid and 
    ! solid-shell coupling in same subroutine.
    !
    IF( IsSolid .OR. IsShell ) THEN  
      ncount = 0
      DO i=1,Mesh % NumberOfNodes
        jf = FPerm(i)      
        js = SPerm(i)

        ! We set coupling at nodes that belong to both equations.
        IF( jf == 0 .OR. js == 0 ) CYCLE
        ncount = ncount + 1

        ! For the given node go through all displacement components. 
        DO j = 1, dim
          ! Indices for matrix rows
          kf = fdofs*(jf-1)+j
          ks = sdofs*(js-1)+j

          ! This is the original diagonal entry of the stiffness matrix.
          ! Let's keep it so that Dirichlet conditions are ideally set. 
          vdiag = A_f % Values( A_f % Diag(kf) ) 

          ! Copy the force from rhs from "F" to "S" and zero it
          A_s % rhs(ks) = A_s % rhs(ks) + A_f % rhs(kf)
          A_f % rhs(kf) = 0.0_dp

          ! Copy the force in implicit form from "F" to "SF" coupling matrix, and zero it.
          ! Now the solid equation includes forces of both equations. 
          DO k=A_f % Rows(kf),A_f % Rows(kf+1)-1
            IF( .NOT. ConstrainedS(ks) ) THEN        
              CALL AddToMatrixElement(A_sf,ks,A_f % Cols(k), A_f % Values(k) )
            END IF
            A_f % Values(k) = 0.0_dp

            ! We zero the mass associated to the Dirichlet conditions since
            ! otherwise the inertia will affect the condition.
            ! We use mass lumping since not all dofs of shell are present in the solid. 
            IF( DoMass ) THEN
              A_s % MassValues(A_s % Diag(ks)) = A_s % MassValues(A_s % Diag(ks)) + A_f % MassValues(k)
              A_f % MassValues(k) = 0.0_dp
            END IF
            IF( DoDamp) THEN
              A_f % DampValues(k) = 0.0_dp
            END IF            
          END DO
          
          ! Set Dirichlet Condition to "F" such that it is equal to "S".
          ! Basically we could eliminate displacement condition and do this afterwards
          ! but this is much more economical. 
          A_f % Values( A_f % Diag(kf)) = vdiag          
          CALL AddToMatrixElement(A_fs,kf,ks, -vdiag )
          
        END DO
      END DO
    ELSE
      CALL Fatal(Caller,'Coupling type not implemented yet!')
    END IF

      
    IF( A_fs % FORMAT == MATRIX_LIST ) THEN
      CALL List_toCRSMatrix(A_fs)
      CALL List_toCRSMatrix(A_sf)
    END IF
      
    !PRINT *,'interface fs sum:',SUM(A_fs % Values), SUM( ABS( A_fs % Values ) )
    !PRINT *,'interface sf sum:',SUM(A_sf % Values), SUM( ABS( A_sf % Values ) )

    CALL Info(Caller,'Number of nodes on interface: '&
        //TRIM(I2S(ncount)),Level=10)    
    CALL Info(Caller,'Number of entries in slave-master coupling matrix: '&
        //TRIM(I2S(SIZE(A_fs % Values))),Level=10)
    CALL Info(Caller,'Number of entries in master-slave coupling matrix: '&
        //TRIM(I2S(SIZE(A_sf % Values))),Level=10)
    
    CALL Info(Caller,'All done',Level=20)
    
  CONTAINS


    ! This routine determines normals of the solid at the common nodes with shell solver.
    ! The normals are determined by summing up potential outer normals and thereafter
    ! subtracting projections to the shell normals.
    !------------------------------------------------------------------------------------
    SUBROUTINE DetermineCouplingNormals()
      INTEGER, ALLOCATABLE :: CouplingPerm(:)
      REAL(KIND=dp), ALLOCATABLE, TARGET :: CouplingNormals(:,:)
      REAL(KIND=dp), POINTER :: WallNormal(:)
      REAL(KIND=dp) :: Normal(3), sNormal
      INTEGER :: CouplingNodes, n, t, nbulk, nbound
      TYPE(Element_t), POINTER :: Element, Parent1, Parent2
      TYPE(Nodes_t), SAVE :: Nodes
      LOGICAL :: Solid1,Solid2
      

      ! allocate elemental stuff
      n = Mesh % MaxElementNodes
      IF ( .NOT. ASSOCIATED( Nodes % x ) ) THEN
        ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
      END IF
     
      ! Generate the permutation for the common nodes
      n = Mesh % NumberOfNodes
      ALLOCATE(CouplingPerm(n))
      WHERE( FVar % Perm(1:n) > 0 .AND. SVar % Perm(1:n) > 0 )
        CouplingPerm = 1
      ELSE WHERE
        CouplingPerm = 0
      END WHERE
      j = 0 
      DO i=1,n
        IF( CouplingPerm(i) > 0 ) THEN
          j = j + 1
          CouplingPerm(i) = j
        END IF
      END DO
      CouplingNodes = j      
      PRINT *,'number of common nodes:',j

      ALLOCATE( CouplingNormals(j,3) )
      CouplingNormals = 0.0_dp
      
      nbulk = Mesh % NumberOfBulkElements
      nbound = Mesh % NumberOfBoundaryElements
      
      ! Sum up all the wall normals associated to coupling nodes together
      DO t=nbulk+1, nbulk+nbound
        Element => Mesh % Elements(t)

        ! If there a node for which we need normal? 
        IF( COUNT( CouplingPerm( Element % NodeIndexes ) > 0 ) < 2 ) CYCLE

        IF( ANY( SVar % Perm( Element % NodeIndexes ) == 0 ) ) CYCLE
        
        ! This needs to be an element where normal can be defined
        !IF( GetElementDim(Element) /= 2 ) CYCLE
        IF( Element % TYPE % ElementCode > 500 ) CYCLE
        IF( Element % TYPE % ElementCode < 300 ) CYCLE
   
        IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
        
        n = Element % TYPE % NumberOfNodes

        !CALL GetElementNodes(Nodes,Element)
        Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
        
        ! Determine whether parents also are active on the solid
        Solid1 = .FALSE.
        Parent1 => Element % BoundaryInfo % Left
        IF( ASSOCIATED( Parent1 ) ) THEN
          Solid1 = ALL(  SVar % Perm( Parent1 % NodeIndexes ) > 0 )
        END IF
        
        Solid2 = .FALSE.
        Parent2 => Element % BoundaryInfo % Right
        IF( ASSOCIATED( Parent2 ) ) THEN
          Solid2 = ALL(  SVar % Perm( Parent2 % NodeIndexes ) > 0 )
        END IF        

        ! Only consider external walls with just either parent in solid
        IF( .NOT. XOR( Solid1, Solid2 ) ) CYCLE
        
        ! Check that the normal points outward of the solid
        IF( Solid1 ) THEN
          Normal = NormalVector(Element,Nodes,Parent=Parent1)
        ELSE
          Normal = NormalVector(Element,Nodes,Parent=Parent2)
        END IF
        
        n = Element % TYPE % NumberOfNodes
        DO i=1,n          
          j = CouplingPerm( Element % NodeIndexes(i) )
          IF( j == 0 ) CYCLE
          
          ! Note that we assume that normals are consistent in a way that they can be summed up
          ! and do not cancel each other
          WallNormal => CouplingNormals(j,:)
          WallNormal = WallNormal + Normal 
        END DO                  
      END DO

      ! Remove the shell normal from the wall normal
      DO t=1, nbulk+nbound
        Element => Mesh % Elements(t)

        ! If there a node for which we need normal? 
        IF( COUNT( CouplingPerm( Element % NodeIndexes ) > 0 ) < 2 ) CYCLE

        ! Shell must be active for all nodes
        IF( ANY( FVar % Perm( Element % NodeIndexes ) == 0 ) ) CYCLE

        ! This needs to be an element where shell can be solved
        !IF( GetElementDim(Element) /= 2 ) CYCLE
        IF( Element % TYPE % ElementCode > 500 ) CYCLE
        IF( Element % TYPE % ElementCode < 300 ) CYCLE
        
        n = Element % TYPE % NumberOfNodes

        !CALL GetElementNodes(Nodes,Element)
        Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
        Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
        Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)

        ! Normal vector for shell, no need check the sign
        Normal = NormalVector(Element,Nodes)
 
        DO i=1,n
          j = CouplingPerm( Element % NodeIndexes(i) )
          IF( j == 0 ) CYCLE
          WallNormal => CouplingNormals(j,:)
          WallNormal = WallNormal - SUM( WallNormal * Normal ) * Normal
        END DO
      END DO

      ! Finally normalize the normals such that their length is one
      j = 0
      DO i=1,CouplingNodes
        WallNormal => CouplingNormals(i,:)
        sNormal = SQRT( SUM( WallNormal**2) )
        IF( sNormal > 1.0d-3 ) THEN
          WallNormal = WallNormal / sNormal 
          PRINT *,'WallNormal:',WallNormal
        ELSE
          j = j + 1
        END IF
      END DO
      
      IF( j > 0 ) THEN
        CALL Fatal('DetermineCouplingNormals','Could not define normals count: '//TRIM(I2S(j)))
      END IF
       
      
    END SUBROUTINE DetermineCouplingNormals


  END SUBROUTINE StructureCouplingAssembly

  
!---------------------------------------------------------------------------------
!> Multiply a linear system by a constant or a given scalar field.
!
!> There are three multiplication modes:
!> 1) Multiply matrix or rhs with a constant factor
!> 2) Multiply matrix or rhs with a constant factor but only blockwise
!> 3) Multiply matrix or rhs with a vector retrieved by a field variable
!
!> And also three things to multiply:
!> a) The right-hand-side of the linear system
!> b) The matrix part of the linear system
!> c) The diagonal entries of the matrix
!
!> Possible uses of the routine include cases where the user wants to introduce diagonal
!> implicit relaxation to the linear system, or to eliminate some coupling terms in 
!> monolithic systems that make the solution of the linear problems more difficult.
!----------------------------------------------------------------------------------
  SUBROUTINE LinearSystemMultiply( Solver )
!----------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER, POINTER :: Perm(:),Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:),Rhs(:)
    TYPE(Variable_t), POINTER :: ThisVar,CoeffVar
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Coeff,Coeff2
    INTEGER :: i,j,j2,k,l,jk,n,Mode,Dofs
    LOGICAL :: Found, UpdateRhs, Symmetric
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: str, VarName

    Params => Solver % Values
    Mesh => Solver % Mesh
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a Mesh!')
    END IF
    A => Solver % Matrix
    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a matrix equation!')
    END IF
    ThisVar => Solver % Variable
    IF(.NOT. ASSOCIATED( ThisVar ) ) THEN
      CALL Fatal('LinearSystemMultiply','Subroutine requires a default variable to exist!')
    END IF

    Perm => ThisVar % Perm
    Dofs = ThisVar % Dofs
    n = A % NumberOfRows
    Cols => A % Cols
    Rows => A % Rows
    Rhs => A % Rhs
    Values => A % Values
        
    UpdateRhs = ListGetLogical( Params,'Linear System Multiply Consistent',Found)
    Symmetric = ListGetLogical( Params,'Linear System Multiply Symmetric',Found)

    ! First, Multiply the k:th piece of the r.h.s. vector if requested
    !-----------------------------------------------------------
    DO k=1,Dofs
      Mode = 0
      
      WRITE( str,'(A)') 'Linear System Rhs Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the rhs with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        WRITE( str,'(A,I0)') TRIM(str)//' ',k
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Multiplying component ',k,' of the rhs with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        END IF
      END IF
      IF( Mode == 0 ) THEN
        str = 'Linear System Rhs Variable'
        VarName = ListGetString( Params,str,Found ) 
        NULLIFY( CoeffVar ) 
        IF( Found ) THEN
          CoeffVar => VariableGet( Mesh % Variables, VarName )
        ELSE
          WRITE( str,'(A,I0)') TRIM(str)//' ',k
          VarName = ListGetString( Params,str,Found )
          IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
        END IF
        IF( ASSOCIATED( CoeffVar ) ) THEN
          IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
            CALL Fatal('LinearSystemMultiply','Permutations should be the same')
          END IF
          Mode = 3
          WRITE( Message,'(A,I0,A)') 'Multiplying component ',k,' of the rhs with > '//TRIM(VarName)//' <'
          CALL Info('LinearSystemMultiply',Message, Level=6 )

          !PRINT *,'Range:',Mode,MINVAL(CoeffVar % Values),MAXVAL(CoeffVar % Values)
        END IF
      END IF
      IF( Mode == 0 ) CYCLE
 
      IF( Mode == 1 ) THEN
        IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
          Rhs = Coeff * Rhs
        END IF
        EXIT
      ELSE 
        DO j=1,SIZE( Perm ) 
          jk = Dofs*(j-1)+k
          IF( Mode == 3 ) Coeff = CoeffVar % Values(j)
          Rhs( jk ) = Coeff * Rhs( jk )
        END DO
      END IF
    END DO
    ! End of r.h.s. multiplication

    ! Secondly, Multiply the kl block of the matrix
    !------------------------------------------------
    DO k=1,Dofs
      DO l=1,Dofs
        Mode = 0
        str = 'Linear System Matrix Factor'
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 1
          WRITE( Message,'(A,ES12.3)') 'Multiplying the matrix with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        ELSE
          WRITE( str,'(A,I0,I0)') TRIM(str)//' ',k,l
          Coeff = ListGetCReal( Params, str, Found )
          IF( Found ) THEN
            Mode = 2
            WRITE( Message,'(A,I0,I0,A,ES12.3)') 'Multiplying block (',k,l,') of the matrix with ',Coeff
            CALL Info('LinearSystemMultiply',Message, Level=6 )
          END IF
        END IF
        IF( Mode == 0 ) THEN
          str = 'Linear System Matrix Variable'
          VarName = ListGetString( Params,str,Found )
          NULLIFY( CoeffVar ) 
          IF( Found ) THEN
            CoeffVar => VariableGet( Mesh % Variables, str )                                    
          ELSE
            WRITE( str,'(A,I0,I0)') TRIM(str)//' ',k,l
            VarName = ListGetString( Params,str,Found )
            IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
          END IF
          IF( ASSOCIATED( CoeffVar ) ) THEN
            IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
              CALL Fatal('LinearSystemMultiply','Permutations should be the same')
            END IF
            Mode = 3
            WRITE( Message,'(A,I0,I0,A)') 'Multiplying block (',k,l,') of the matrix with > '//TRIM(VarName)//' <'
            CALL Info('LinearSystemMultiply',Message, Level=6 )
          END IF
        END IF

        IF( Mode == 0 ) CYCLE

        IF( Mode == 1 ) THEN
          IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
            Values = Coeff * Values
          END IF
        ELSE 
          DO j=1,SIZE( Perm ) 
            jk = Dofs*(j-1)+k
            IF( Mode == 3 ) Coeff = CoeffVar % Values(j)
            DO i=Rows(jk),Rows(jk+1)-1 
              IF( MODULO( Cols(i), Dofs ) == MODULO( l, Dofs ) ) THEN
                IF( Mode == 3 .AND. Symmetric ) THEN          
                  j2 = (Cols(i)-1)/Dofs + 1 
                  Coeff2 = CoeffVar % Values(j2)
                  Values( i ) = SQRT( Coeff * Coeff2 ) * Values( i )
                ELSE
                  Values( i ) = Coeff * Values( i )
                END IF
              END IF
            END DO
          END DO
        END IF
      END DO
      IF( Mode == 1 ) EXIT
    END DO
    ! end of matrix multiplication


    ! Finally, Multiply the diagonal of the matrix
    !------------------------------------------------
    DO k=1,Dofs
      Mode = 0

      str = 'Linear System Diagonal Factor'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Multiplying the diagonal with ',Coeff
        CALL Info('LinearSystemMultiply',Message, Level=6 )
      ELSE
        WRITE( str,'(A,I0)') TRIM(str)//' ',k
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Multiplying component ',k,' of the matrix diagonal with ',Coeff
          CALL Info('LinearSystemMultiply',Message, Level=6 )          
        END IF
      END IF

      IF( Mode == 0 ) THEN
        str = 'Linear System Diagonal Variable'
        VarName = ListGetString( Params,str,Found )
        NULLIFY( CoeffVar )
        IF( Found ) THEN
          CoeffVar => VariableGet( Mesh % Variables, VarName )
        ELSE
          WRITE( str,'(A,I0)') TRIM(str)//' ',k
          VarName = ListGetString( Params,str,Found )
          IF( Found ) CoeffVar => VariableGet( Mesh % Variables, VarName )
        END IF
        IF( ASSOCIATED( CoeffVar ) ) THEN
          IF( ANY( CoeffVar % Perm /= Perm ) ) THEN
            CALL Fatal('LinearSystemMultiply','Permutations should be the same')
          END IF
          Mode = 3
          WRITE( Message,'(A,I0,A)') 'Multiplying component ',k,' of the diagonal with > '//TRIM(VarName)//' <'
          CALL Info('LinearSystemMultiply',Message, Level=6 )
        END IF
      END IF
      
      IF( Mode == 0 ) CYCLE

      IF( Mode == 1 ) THEN
        IF( ABS( Coeff - 1.0_dp ) > EPSILON( Coeff ) ) THEN
          IF( UpdateRhs ) Rhs = Rhs + ( Coeff - 1 ) * Values( A % Diag ) * ThisVar % Values
          Values( A % Diag ) = Coeff * Values( A % Diag )
        END IF
        EXIT
      ELSE 
        DO j=1,SIZE( Perm ) 
          jk = Dofs*(j-1)+k
          IF( Mode == 3 ) Coeff = CoeffVar % Values(j)

          IF( UpdateRhs ) Rhs( jk ) = Rhs( jk ) + (Coeff - 1) * Values(A % Diag(jk)) * ThisVar % Values(jk)
          Values( A % Diag(jk) ) = Coeff * Values( A % Diag(jk) )
        END DO
      END IF
    END DO
    ! end of diagonal multiplication


  END SUBROUTINE LinearSystemMultiply






!---------------------------------------------------------------------------------
!> Set the diagonal entry to given minimum.
!----------------------------------------------------------------------------------
  SUBROUTINE LinearSystemMinDiagonal( Solver )
!----------------------------------------------------------------------------------    
    TYPE(Solver_t) :: Solver
    !------------------------------------------------------------------------------
    INTEGER, POINTER :: Perm(:),Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:),Rhs(:)
    TYPE(Variable_t), POINTER :: ThisVar
    TYPE(Matrix_t), POINTER :: A
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Coeff
    INTEGER :: i,j,j2,k,l,jk,n,Mode,Dofs
    LOGICAL :: Found, UpdateRhs, Symmetric
    TYPE(ValueList_t), POINTER :: Params
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    INTEGER :: NoSet
    REAL(KIND=dp) :: DiagSum, val, DiagMax

    Params => Solver % Values
    Mesh => Solver % Mesh
    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a Mesh!')
    END IF
    A => Solver % Matrix
    IF(.NOT. ASSOCIATED( A ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a matrix equation!')
    END IF
    ThisVar => Solver % Variable
    IF(.NOT. ASSOCIATED( ThisVar ) ) THEN
      CALL Fatal('LinearSystemMinDiagonal','Subroutine requires a default variable to exist!')
    END IF

    Perm => ThisVar % Perm
    Dofs = ThisVar % Dofs
    n = A % NumberOfRows
    Cols => A % Cols
    Rows => A % Rows
    Rhs => A % Rhs
    Values => A % Values

    ! Set the minimum value for each component, only nodel dofs considered
    !---------------------------------------------------------------------
    NoSet = 0
    DiagMax = 0.0_dp
    DiagSum = 0.0_dp
    n = MAXVAL( Perm ( 1:Mesh % NumberOfNodes ) )

    DO k=1,Dofs
      Mode = 0

      str = 'Linear System Diagonal Min'
      Coeff = ListGetCReal( Params, str, Found )
      IF( Found ) THEN
        Mode = 1
        WRITE( Message,'(A,ES12.3)') 'Setting minimum of the diagonal to ',Coeff
        CALL Info('LinearSystemMinDiagonal',Message, Level=6 )
      ELSE
        WRITE( str,'(A,I0)') TRIM(str)//' ',k
        Coeff = ListGetCReal( Params, str, Found )
        IF( Found ) THEN
          Mode = 2 
          WRITE( Message,'(A,I0,A,ES12.3)') 'Setting minimum of diagonal component ',k,' to ',Coeff
          CALL Info('LinearSystemMinDiagonal',Message, Level=6 )          
        END IF
      END IF
      
      IF( Mode == 0 ) CYCLE
      
      DO j=1,n
        jk = Dofs*(j-1)+k
        l = A % Diag(jk) 
        IF( l == 0 ) CYCLE

        ! Only add the diagonal to the owned dof
        IF( ParEnv % PEs > 1 ) THEN
          IF( A % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPE ) CYCLE
        END IF

        val = ABS( Values( l ) )
        DiagSum = DiagSum + val
        DiagMax = MAX( DiagMax, val )
        IF( val < Coeff ) THEN
          Values( A % Diag(jk) ) = Coeff
          NoSet = NoSet + 1
        END IF
      END DO
    END DO

    CALL Info('LinearSystemMinDiagonal','Number of diagonal values set to minimum: '//TRIM(I2S(NoSet)),Level=5)
    WRITE( Message,'(A,ES12.3)') 'Average abs(diagonal) entry: ',DiagSum / n
    CALL Info('LinearSystemMinDiagonal',Message, Level=6 )

    WRITE( Message,'(A,ES12.3)') 'Maximum abs(diagonal) entry: ',DiagMax
    CALL Info('LinearSystemMinDiagonal',Message, Level=6 )


  END SUBROUTINE LinearSystemMinDiagonal





  !----------------------------------------------------------------------
  !> Make the high-order flux corrected transport (FCT) correction after 
  !> the low order approximation has been solved. 
  !
  !> For more information see, for example, 
  !> Dmitri Kuzmin (2008): "Explicit and implicit FEM-FCT algorithms with flux linearization"
  !----------------------------------------------------------------------
  SUBROUTINE FCT_Correction( Solver  )

    TYPE(Solver_t), POINTER :: Solver

    TYPE(Valuelist_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: n,i,j,k,k2
    INTEGER, POINTER :: Rows(:),Cols(:),Perm(:)
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp), POINTER :: BulkValues(:),u(:),M_L(:),udot(:), &
        pp(:),pm(:),qp(:),qm(:),corr(:),ku(:),ulow(:)
    REAL(KIND=dp), ALLOCATABLE :: mc_udot(:), fct_u(:)
    REAL(KIND=dp), POINTER CONTIG :: M_C(:), SaveValues(:)
    REAL(KIND=dp) :: rsum, Norm,m_ij,k_ij,du,d_ij,f_ij,c_ij,Ceps,CorrCoeff,&
        rmi,rmj,rpi,rpj,dt
    TYPE(Variable_t), POINTER :: Var, Variables
    LOGICAL :: Found, Symmetric, SaveFields, SkipCorrection
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, TmpVarName

    REAL(KIND=dp), POINTER :: mmc(:), mmc_h(:), fct_d(:)
    TYPE(Element_t), POINTER :: Element
    LOGICAL, ALLOCATABLE :: ActiveNodes(:)

    Params => Solver % Values

    SkipCorrection = ListGetLogical( Params,'FCT Correction Skip',Found )
    Symmetric = ListGetLogical( Params,'FCT Correction Symmetric',Found )
    SaveFields = ListGetLogical( Params,'FCT Correction Save',Found )

    IF( SkipCorrection .AND. .NOT. SaveFields ) THEN
      CALL Info('FCT_Correction','Skipping the computation of FCT correction',Level=5)
    END IF

    CALL Info('FCT_Correction','Computing corrector for the low order solution',Level=5)

    ! PRINT *,'FCT Norm Before Correction:',SQRT( SUM( Solver % Variable % Values**2) )
 
    Mesh => Solver % Mesh
    Variables => Solver % Mesh % Variables
 
    ! Set pointers 
    A => Solver % Matrix
    n = A % NumberOfRows
    Rows => A % Rows
    Cols => A % Cols

    BulkValues => A % BulkValues
    M_C => A % MassValues
    Perm => Solver % Variable % Perm

    M_L => A % MassValuesLumped 
    IF (ParEnv % PEs>1) CALL ParallelSumVector(A,M_L)

    Var => VariableGet( Variables,'timestep size')
    dt = Var % Values(1) 

    ! low order solution at the start, high order in the end
    u => Solver % Variable % Values
    VarName = GetVarName(Solver % Variable)

    ! Here a bunch of vectors are stored for visualization and debugging purposes
    !----------------------------------------------------------------------------

    ! low order solution is the solution without corrections
    ! This is created and saved only if requested
    !---------------------------------------------------------------------------
    IF( SaveFields ) THEN
      TmpVarName = TRIM( VarName )//' fctlow'    
      Var => VariableGet( Variables, TmpVarName )
      IF( .NOT. ASSOCIATED(Var) ) THEN
        CALL VariableAddVector( Variables, Mesh, Solver,&
            TmpVarName, Perm = Perm, Output = SaveFields )
        Var => VariableGet( Variables, TmpVarName )
      END IF
      ulow => Var % Values
      ulow = u
    END IF

    ! Create auxiliary vectors for the fct algorithm
    !---------------------------------------------------------------------------

    ! r.h.s. term
    TmpVarName = TRIM( VarName )//' fctku'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    ku => Var % Values

    ! time derivative from lower order analysis
    TmpVarName = TRIM( VarName )//' fctudot'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    udot => Var % Values

    ! Fields related to flux limiters
    TmpVarName = TRIM( VarName )//' fctpp'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    pp => Var % Values
    
    TmpVarName = TRIM( VarName )//' fctpm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    pm => Var % Values
    
    TmpVarName = TRIM( VarName )//' fctqp'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    qp => Var % Values

    TmpVarName = TRIM( VarName )//' fctqm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    qm => Var % Values

    TmpVarName = TRIM( VarName )//' fctmm'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    Var % Values = M_L

    ! higher order correction 
    TmpVarName = TRIM( VarName )//' fctcorr'    
    Var => VariableGet( Variables, TmpVarName )
    IF( .NOT. ASSOCIATED(Var) ) THEN
      CALL VariableAddVector( Variables, Mesh, Solver,&
          TmpVarName, Perm = Perm, Output = SaveFields )
      Var => VariableGet( Variables, TmpVarName )
    END IF
    corr => Var % Values


    ! 1) Compute the nodal time derivatives
    ! M_C*udot=K*ulow  (M_C is the consistent mass matrix)
    !----------------------------------------------------------------------
    CALL Info('FCT_Correction','Compute nodal time derivatives',Level=10)
    ! Compute: ku = K*ulow
#if 0
    DO i=1,n
      rsum = 0.0_dp
      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        K_ij = BulkValues(k)
        rsum = rsum + K_ij * u(j) 
      END DO
      ku(i) = rsum
    END DO
    ! Solve the linear system for udot
    ! The stiffness matrix is momentarily replaced by the consistent mass matrix M_C
    ! Also the namespace is replaced to 'fct:' so that different strategies may 
    ! be applied to the mass matrix solution.
    CALL ListPushNameSpace('fct:')
    SaveValues => A % Values
    A % Values => M_C
    CALL SolveLinearSystem( A, ku, udot, Norm, 1, Solver )
    A % Values => SaveValues
    CALL ListPopNamespace()
#else

  BLOCK
    REAL(KIND=dp), ALLOCATABLE :: TmpRhsVec(:), TmpXVec(:)

    SaveValues => A % Values

    IF (Parenv % PEs>1) THEN
      ALLOCATE(TmpRHSVec(n), TmpXVec(n))
      TmpxVec = 0._dp; tmpRHSVec = 0._dp

      A % Values => BulkValues
      CALL ParallelInitSolve(A,TmpXVec,TmpRhsVec,u)
      CALL ParallelVector(A,TmpRhsvec,u)

      CALL ParallelMatrixVector(A,TmpRhsvec,TmpXVec)

      CALL PartitionVector(A,Ku,TmpXVec)
      DEALLOCATE(TmpRhsVec, TmpXVec)
    ELSE
      DO i=1,n
        rsum = 0._dp
        DO k=Rows(i),Rows(i+1)-1
          j = Cols(k)
          K_ij = BulkValues(k)
          rsum = rsum + K_ij * u(j) 
        END DO
        ku(i) = rsum
      END DO
    END IF

    CALL ListPushNameSpace('fct:')
    A % Values => M_C
    udot = 0._dp
    CALL SolveLinearSystem(A,Ku,Udot,Norm,1,Solver)
    CALL ListPopNamespace()

    A % Values => SaveValues
  END BLOCK
#endif

    ! Computation of correction factors (Zalesak's limiter)
    ! Code derived initially from Kuzmin's subroutine   
    !---------------------------------------------------------
    CALL Info('FCT_Correction','Compute correction factors',Level=10)
    pp = 0 
    pm = 0
    qp = 0 
    qm = 0

    IF(ParEnv % PEs>1) THEN
      fct_d => A % FCT_D
      mmc    => A % MassValues
      mmc_h  => A % HaloMassValues

      ALLOCATE(ActiveNodes(n)); activeNodes=.FALSE.
      DO i=1,Solver % NumberOfActiveElements
        Element => Solver % Mesh % Elements(Solver % ActiveElements(i))
        IF ( Element % PartIndex /= ParEnv % MyPE ) CYCLE
        Activenodes(Solver % Variable % Perm(Element % NodeIndexes)) = .TRUE.
      END DO
    ELSE
      fct_d => A % FCT_D
      mmc => A % MassValues
    END IF
    DO i=1,n
      IF (ParEnv % PEs > 1 ) THEN
        IF ( .NOT. ActiveNodes(i) ) CYCLE
      end if

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)

        IF (ParEnv % PEs>1) THEN
          IF ( .NOT.ActiveNodes(j) ) CYCLE
        END IF

        ! Compute the raw antidiffusive fluxes
        ! f_ij=m_ij*[udot(i)-udot(j)]+d_ij*[ulow(i)-ulow(j)]
        !-----------------------------------------------------
        ! d_ij and m_ij are both symmetric
        ! Hence F_ji = -F_ij
           
        f_ij = mmc(k)*(udot(i)-udot(j)) + fct_d(k)*(u(i)-u(j))
        IF ( ParEnv % PEs>1 ) f_ij=f_ij+mmc_h(k)*(udot(i)-udot(j))
        ! Compared to Kuzmin's paper F_ij=-F_ij since d_ij and 
        ! udot have different signs. 
        f_ij = -f_ij

        ! Antidiffusive fluxes to be limited
        du = u(j)-u(i)

        ! Prelimiting of antidiffusive fluxes i.e. du and the flux have different signs
        IF (f_ij*du >= TINY(du)) THEN
          f_ij = 0._dp
        ELSE        
          ! Positive/negative edge contributions
          pp(i) = pp(i) + MAX(0._dp,f_ij)
          pm(i) = pm(i) + MIN(0._dp,f_ij)
        END IF

        ! Maximum/minimum solution increments
        qp(i) = MAX(qp(i),du)
        qm(i) = MIN(qm(i),du)
      END DO
    END DO

    ! Computation of nodal correction factors
    ! These are eliminated as vectors to save some space
    ! and replaced by rpi, rpj, rmi, rmj
    ! DO i=1,n
    !  IF( pp(i) > Ceps ) THEN
    !    rp(i) = MIN( 1.0_dp, M_L(i)*qp(i)/pp(i) )
    !  END IF
    !  IF( pm(i) < -Ceps ) THEN
    !    rm(i) = MIN( 1.0_dp, M_L(i)*qm(i)/pm(i) )
    !  END IF
    ! END DO

    ! Correct the low-order solution
    ! (M_L*ufct)_i=(M_L*ulow)_i+dt*sum(alpha_ij*f_ij)
    !-------------------------------------------------
    ! Symmetric flux limiting
    ! Correction of the right-hand side


!   IF (ParEnv % PEs>1) THEN
!     CALL ParallelSumVector(A,pm)
!     CALL ParallelSumVector(A,pp)
!     CALL ParallelSumVector(A,qm)
!     CALL ParallelSumVector(A,qp)
!   END IF

    CorrCoeff = ListGetCReal( Params,'FCT Correction Coefficient',Found )
    IF( .NOT. Found ) CorrCoeff = 1._dp

    Ceps = TINY( Ceps )
    corr = 0._dp
    DO i=1,n
      IF (ParEnv % PEs>1) THEN
        IF( .NOT. ActiveNodes(i)) CYCLE
      END IF

      IF( pp(i) > Ceps ) THEN
        rpi = MIN( 1._dp, M_L(i)*qp(i)/pp(i) )
      ELSE
        rpi = 0._dp
      END IF

      IF( pm(i) < -Ceps ) THEN
        rmi = MIN( 1._dp, M_L(i)*qm(i)/pm(i) )
      ELSE
        rmi = 0._dp
      END IF

      DO k=Rows(i),Rows(i+1)-1
        j = Cols(k)
        IF(ParEnv % PEs>1) THEN
          IF(.NOT.ActiveNodes(j)) CYCLE
        END IF

        f_ij = mmc(k)*(udot(i)-udot(j))  + fct_d(k)*(u(i)-u(j))
        IF (ParEnv % PEs>1) f_ij = f_ij + mmc_h(k)*(udot(i)-udot(j))
        f_ij = -f_ij

        IF (f_ij > 0) THEN 
          IF( pm(j) < -Ceps ) THEN
            rmj = MIN( 1.0_dp, M_L(j)*qm(j)/pm(j) )
          ELSE
            rmj = 0._dp
          END IF
          c_ij = MIN(rpi,rmj)
        ELSE 
          IF( pp(j) > Ceps ) THEN
            rpj = MIN( 1._dp, M_L(j)*qp(j)/pp(j) )
          ELSE
            rpj = 0._dp
          END IF
          c_ij = MIN(rmi,rpj)
        END IF
        corr(i) = corr(i) + c_ij * f_ij
      END DO
      corr(i) = CorrCoeff * corr(i) / M_L(i)
    END DO

    IF (ParEnv % PEs>1) THEN
!     CALL ParallelSumVector(A,corr)
      DEALLOCATE(A % HaloValues, A % HaloMassValues)
      A % HaloValues => Null(); A % HaloMassValues => Null()
    END IF

    ! Optionally skip applying the correction, just for debugging purposes
    IF( SkipCorrection ) THEN
      CALL Info('FCT_Correction','Skipping Applying corrector',Level=4)
    ELSE
      CALL Info('FCT_Correction','Applying corrector for the low order solution',Level=10)

      u = u + corr

      ! PRINT *,'FCT Norm After Correction:',SQRT( SUM( Solver % Variable % Values**2) )
    END IF

  END SUBROUTINE FCT_Correction



  ! Create Linear constraints from mortar BCs:
  ! -------------------------------------------   
  SUBROUTINE GenerateProjectors(Model,Solver,Nonlinear,SteadyState) 
    
     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver
     LOGICAL, OPTIONAL :: Nonlinear, SteadyState

     LOGICAL :: IsNonlinear,IsSteadyState,Timing, RequireNonlinear, ContactBC
     LOGICAL :: ApplyMortar, ApplyContact, Found
     INTEGER :: i,j,k,l,n,dsize,size0,col,row,dim
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: CM, CMP, CM0, CM1
     TYPE(Variable_t), POINTER :: DispVar
     REAL(KIND=dp) :: t0,rt0,rst,st,ct
     CHARACTER(*), PARAMETER :: Caller = 'GenerateProjectors'
     
     ApplyMortar = ListGetLogical(Solver % Values,'Apply Mortar BCs',Found) 
     ApplyContact = ListGetLogical(Solver % Values,'Apply Contact BCs',Found) 

     IF( .NOT. ( ApplyMortar .OR. ApplyContact) ) RETURN
     
     i = ListGetInteger( Solver % Values,'Mortar BC Master Solver',Found ) 
     IF( Found ) THEN
       Solver % MortarBCs => CurrentModel % Solvers(i) % MortarBCs
       IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
         CALL Fatal(Caller,'Could not reuse projectors from solver: '//TRIM(I2S(i)))
       END IF
       CALL Info(Caller,'Reusing projectors from solver: '//TRIM(I2S(i)),Level=8)
       RETURN
     END IF

     CALL Info(Caller,'Generating mortar projectors',Level=8)

     Timing = ListCheckPrefix(Solver % Values,'Projector Timing')
     IF( Timing ) THEN
       t0 = CPUTime(); rt0 = RealTime()      
     END IF

     IsNonlinear = .FALSE.
     IF( PRESENT( Nonlinear ) ) IsNonlinear = Nonlinear
     IsSteadyState = .NOT. IsNonlinear

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) THEN
       ALLOCATE( Solver % MortarBCs( Model % NumberOfBCs ) )
       DO i=1, Model % NumberOfBCs
         Solver % MortarBCs(i) % Projector => NULL()
       END DO
     END IF
     
     dim = CoordinateSystemDimension()

     DO i=1,Model % NumberOFBCs
       BC => Model % BCs(i) % Values
       
       ContactBC = .FALSE.
       j = ListGetInteger( BC,'Mortar BC',Found)       
       IF( .NOT. Found ) THEN
         j = ListGetInteger( BC,'Contact BC',Found)       
         ContactBC = Found
       END IF
       IF( .NOT. Found ) CYCLE

       RequireNonlinear = ListGetLogical( BC,'Mortar BC Nonlinear',Found)
       IF( .NOT. Found ) THEN
         RequireNonlinear = ContactBC .AND. .NOT. ListGetLogical( BC,'Tie Contact',Found )
       END IF

       IF( IsNonlinear ) THEN
         IF( .NOT. RequireNonlinear ) CYCLE
       ELSE
         IF( RequireNonlinear ) CYCLE
       END IF             

       IF( ASSOCIATED( Solver % MortarBCs(i) % Projector ) ) THEN
         IF( ListGetLogical( BC,'Mortar BC Static',Found) ) CYCLE         
         
         IF( ASSOCIATED( Solver % MortarBCs(i) % Projector % Ematrix ) ) THEN
           CALL FreeMatrix( Solver % MortarBCs(i) % Projector % Ematrix )
         END IF
         CALL FreeMatrix( Solver % MortarBCs(i) % Projector )
       END IF
       
       Solver % MortarBCs(i) % Projector => &
           PeriodicProjector(Model,Solver % Mesh,i,j,dim,.TRUE.)
       
       IF( ASSOCIATED( Solver % MortarBCs(i) % Projector ) ) THEN
         Solver % MortarBCsChanged = .TRUE.
       END IF

     END DO


     IF( Timing ) THEN
       st  = CPUTime() - t0;
       rst = RealTime() - rt0
       
       WRITE(Message,'(a,f8.2,f8.2,a)') 'Projector creation time (CPU,REAL) for '&
           //GetVarName(Solver % Variable)//': ',st,rst,' (s)'
       CALL Info(Caller,Message)    
       
       IF( ListGetLogical(Solver % Values,'Projector Timing',Found)) THEN
         CALL ListAddConstReal(CurrentModel % Simulation,'res: projector cpu time '&
             //GetVarName(Solver % Variable),st)
         CALL ListAddConstReal(CurrentModel % Simulation,'res: projector real time '&
             //GetVarName(Solver % Variable),rst)
       END IF

       IF( ListGetLogical(Solver % Values,'Projector Timing Cumulative',Found)) THEN
         ct = ListGetConstReal(CurrentModel % Simulation,'res: cum projector cpu time '&
             //GetVarName(Solver % Variable),Found)
         st = st + ct
         ct = ListGetConstReal(CurrentModel % Simulation,'res: cum projector real time '&
             //GetVarName(Solver % Variable),Found)
         rst = rst + ct
         CALL ListAddConstReal(CurrentModel % Simulation,'res: cum projector cpu time '&
             //GetVarName(Solver % Variable),st)
         CALL ListAddConstReal(CurrentModel % Simulation,'res: cum projector real time '&
             //GetVarName(Solver % Variable),rst)
       END IF
     END IF
     
   END SUBROUTINE GenerateProjectors



   ! Generate constraint matrix from mortar projectors. 
   ! This routine takes each boundary projector and applies it 
   ! to the current field variable (scalar or vector) merging 
   ! all into one single projector. 
   !---------------------------------------------------------
   SUBROUTINE GenerateConstraintMatrix( Model, Solver )

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     INTEGER, POINTER :: Perm(:)
     INTEGER :: i,j,j2,k,k2,l,l2,dofs,maxperm,permsize,bc_ind,constraint_ind,row,col,col2,mcount,bcount,kk
     TYPE(Matrix_t), POINTER :: Atmp,Btmp, Ctmp
     LOGICAL :: AllocationsDone, CreateSelf, ComplexMatrix, TransposePresent, Found, &
         SetDof, SomeSet, SomeSkip, SumProjectors, NewRow, SumThis
     INTEGER, ALLOCATABLE :: SumPerm(:),SumCount(:)
     LOGICAL, ALLOCATABLE :: ActiveComponents(:), SetDefined(:)
     TYPE(ValueList_t), POINTER :: BC
     TYPE(MortarBC_t), POINTER :: MortarBC
     REAL(KIND=dp) :: wsum, Scale
     INTEGER :: rowoffset, arows, sumrow, EliminatedRows, NeglectedRows, sumrow0, k20
     CHARACTER(LEN=MAX_NAME_LEN) :: Str
     LOGICAL :: ThisIsMortar, Reorder
     LOGICAL :: AnyPriority
     INTEGER :: Priority, PrevPriority
     INTEGER, ALLOCATABLE :: BCOrdering(:), BCPriority(:)
     LOGICAL :: NeedToGenerate, ComplexSumRow 

     LOGICAL :: HaveMortarDiag, LumpedDiag, PerFlipActive
     REAL(KIND=dp) :: MortarDiag, val, valsum, EpsVal
     LOGICAL, POINTER :: PerFlip(:)
     CHARACTER(*), PARAMETER :: Caller = 'GenerateConstraintMatrix'

     
     ! Should we genarete the matrix
     NeedToGenerate = Solver % MortarBCsChanged

     PerFlipActive = Solver % PeriodicFlipActive
     IF( PerFlipActive ) THEN
       CALL Info(Caller,'Periodic flip is active',Level=8)
       PerFlip => Solver % Mesh % PeriodicFlip           
     END IF
     
     ! Set pointers to save the initial constraint matrix
     ! We assume that the given ConstraintMatrix is constant but we have consider it the 1st time
     IF(.NOT. Solver % ConstraintMatrixVisited ) THEN       
       IF( ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
         CALL Info(Caller,'Saving initial constraint matrix to Solver',Level=12)
         Solver % ConstraintMatrix => Solver % Matrix % ConstraintMatrix
         Solver % Matrix % ConstraintMatrix => NULL()
         NeedToGenerate = .TRUE. 
       END IF
       Solver % ConstraintMatrixVisited = .TRUE.
     END IF
     
     IF( NeedToGenerate ) THEN
       CALL Info(Caller,'Building constraint matrix',Level=12)
     ELSE     
       CALL Info(Caller,'Nothing to do for now',Level=12)
       RETURN
     END IF
       
     
     ! Compute the number and size of initial constraint matrices
     !-----------------------------------------------------------
     row    = 0
     mcount = 0
     bcount = 0
     Ctmp => Solver % ConstraintMatrix
     IF( ASSOCIATED( Ctmp ) ) THEN
       DO WHILE(ASSOCIATED(Ctmp))
         mcount = mcount + 1
         row = row + Ctmp % NumberOfRows
         Ctmp => Ctmp % ConstraintMatrix
       END DO
       CALL Info(Caller,'Number of initial constraint matrices: '//TRIM(I2S(mcount)),Level=12)       
     END IF
       
     
     ! Compute the number and size of mortar matrices
     !-----------------------------------------------
     IF( ASSOCIATED( Solver % MortarBCs ) ) THEN
       DO bc_ind=1,Model % NumberOFBCs
         Atmp => Solver % MortarBCs(bc_ind) % Projector
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         bcount = bcount + 1
         row = row + Atmp % NumberOfRows
       END DO
       CALL Info(Caller,'Number of mortar matrices: '//TRIM(I2S(bcount)),Level=12)       
     END IF
     
     IF( row==0 ) THEN
       CALL Info(Caller,'Nothing to do since there are no constrained dofs!',Level=12)       
       RETURN
     END IF

     MortarDiag = ListGetCReal( Solver % Values,'Mortar Diag',HaveMortarDiag )
     LumpedDiag = ListGetLogical( Solver % Values,'Lumped Diag',Found )

     IF( HaveMortarDiag ) THEN
       CALL Info(Caller,&
           'Adding diagonal entry to mortar constraint!',Level=12)              
     END IF
     
     IF( mcount == 1 .AND. bcount == 0 .AND. .NOT. HaveMortarDiag ) THEN
       CALL Info(Caller,'Using initial constraint matrix',Level=12)       
       Solver % Matrix % ConstraintMatrix => Solver % ConstraintMatrix
       RETURN
     END IF

     ! Now we are generating something more complex and different than last time
     IF( ASSOCIATED( Solver % Matrix % ConstraintMatrix ) ) THEN
       CALL Info(Caller,'Releasing previous constraint matrix',Level=12)
       CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
       Solver % Matrix % ConstraintMatrix => NULL()
     END IF
       
     EpsVal = ListGetConstReal( Solver % Values,&
         'Minimum Projector Value', Found )
     IF(.NOT. Found ) EpsVal = 1.0d-8
     
     
     SumProjectors = ListGetLogical( Solver % Values,&
         'Mortar BCs Additive', Found )
     IF( .NOT. Found ) THEN
       IF( bcount > 1 .AND. ListGetLogical( Solver % Values, &
           'Eliminate Linear Constraints',Found ) ) THEN
         CALL Info(Caller,&
             'Enforcing > Mortar BCs Additive < to True to enable elimination',Level=8)
         SumProjectors = .TRUE.
       END IF       
       IF( .NOT. SumProjectors .AND. ListGetLogical( Solver % Values, &
           'Apply Conforming BCs',Found ) ) THEN
         CALL Info(Caller,&
             'Enforcing > Mortar BCs Additive < to True because of conforming BCs',Level=8)
         SumProjectors = .TRUE.
       END IF
     END IF
     EliminatedRows = 0

     CALL Info(Caller,'There are '&
         //TRIM(I2S(row))//' initial rows in constraint matrices',Level=10)
     
     dofs = Solver % Variable % DOFs
     Perm => Solver % Variable % Perm
     permsize = SIZE( Perm )
     maxperm  = MAXVAL( Perm )
     AllocationsDone = .FALSE.
     arows = Solver % Matrix % NumberOfRows
     
     ALLOCATE( ActiveComponents(dofs), SetDefined(dofs) ) 
     
     IF( SumProjectors ) THEN
       ALLOCATE( SumPerm( dofs * permsize ) )
       SumPerm = 0
       ALLOCATE( SumCount( arows ) )
       SumCount = 0
     END IF
     
     ComplexMatrix = Solver % Matrix % Complex
     ComplexSumRow = .FALSE.
     
     IF( ComplexMatrix ) THEN
       IF( MODULO( Dofs,2 ) /= 0 ) CALL Fatal(Caller,&
           'Complex matrix should have even number of components!')
     ELSE
       ! Currently complex matrix is enforced if there is an even number of 
       ! entries since it seems that we cannot rely on the flag to be set.
       ComplexMatrix = ListGetLogical( Solver % Values,'Linear System Complex',Found )
       IF( .NOT. Found ) ComplexMatrix = ( MODULO( Dofs,2 ) == 0 )
     END IF

     
     AnyPriority = ListCheckPresentAnyBC( Model,'Projector Priority') 
     IF( AnyPriority ) THEN
       IF(.NOT. SumProjectors ) THEN
         CALL Warn(Caller,'Priority has effect only in additive mode!')
         AnyPriority = .FALSE.
       ELSE
         CALL Info(Caller,'Using priority for projector entries',Level=7)
         ALLOCATE( BCPriority(Model % NumberOfBCs), BCOrdering( Model % NumberOfBCs) )
         BCPriority = 0; BCOrdering = 0
         DO bc_ind=1, Model % NumberOFBCs
           Priority = ListGetInteger( Model % BCs(bc_ind) % Values,'Projector Priority',Found)
           BCPriority(bc_ind) = -bc_ind + Priority * Model % NumberOfBCs 
           BCOrdering(bc_ind) = bc_ind
         END DO
         CALL SortI( Model % NumberOfBCs, BCPriority, BCOrdering )
       END IF
     END IF
     NeglectedRows = 0


100  sumrow = 0
     k2 = 0
     rowoffset = 0
     Priority = -1
     PrevPriority = -1
     sumrow0 = 0
     k20 = 0
     
     TransposePresent = .FALSE.
     Ctmp => Solver % ConstraintMatrix

     DO constraint_ind = Model % NumberOFBCs+mcount,1,-1
       
       ! This is the default i.e. all components are applied mortar BCs
       ActiveComponents = .TRUE.
       
       IF(constraint_ind > Model % NumberOfBCs) THEN
         ThisIsMortar = .FALSE.
         SumThis = .FALSE.
         Atmp => Ctmp
         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE
         Ctmp => Ctmp % ConstraintMatrix
         IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
           IF(.NOT. AllocationsDone ) THEN
             CALL Warn(Caller,'InvPerm is expected, using identity!')
           END IF
         END IF
         CALL Info(Caller,'Adding initial constraint matrix: '&
             //TRIM(I2S(constraint_ind - Model % NumberOfBCs)),Level=8)         
       ELSE
         ThisIsMortar = .TRUE.
         SumThis = SumProjectors
         IF( AnyPriority ) THEN
           bc_ind = BCOrdering(constraint_ind)
         ELSE
           bc_ind = constraint_ind 
         END IF

         MortarBC => Solver % MortarBCs(bc_ind) 
         Atmp => MortarBC % Projector

         IF( .NOT. ASSOCIATED( Atmp ) ) CYCLE

         IF(.NOT. AllocationsDone ) THEN
           CALL Info(Caller,'Adding projector for BC: '//TRIM(I2S(bc_ind)),Level=8)
         END IF
           
         IF( .NOT. ASSOCIATED( Atmp % InvPerm ) ) THEN
           CALL Fatal(Caller,'InvPerm is required!')
         END IF

         IF( AnyPriority ) THEN
           Priority = ListGetInteger( Model % BCs(bc_ind) % Values,'Projector Priority',Found)
         END IF
           
         ! Enable that the user can for vector valued cases either set some 
         ! or skip some field components. 
         SomeSet = .FALSE.
         SomeSkip = .FALSE.
         DO i=1,Dofs
           IF( Dofs > 1 ) THEN
             str = ComponentName( Solver % Variable, i )
           ELSE
             str = Solver % Variable % Name 
           END IF

           SetDof = ListGetLogical( Model % BCs(bc_ind) % Values,'Mortar BC '//TRIM(str),Found )

           SetDefined(i) = Found
           IF(Found) THEN
             ActiveComponents(i) = SetDof
             IF( SetDof ) THEN
               SomeSet = .TRUE.
             ELSE
               SomeSkip = .TRUE.
             END IF
           END IF
         END DO
         
         ! By default all components are applied mortar BC and some are turned off.
         ! If the user does the opposite then the default for other components is True.
         IF( SomeSet .AND. .NOT. ALL(SetDefined) ) THEN
           IF( SomeSkip ) THEN
             CALL Fatal(Caller,'Do not know what to do with all components')
           ELSE
             CALL Info(Caller,'Unspecified components will not be set for BC '//TRIM(I2S(bc_ind)),Level=10)
             DO i=1,Dofs
               IF( .NOT. SetDefined(i) ) ActiveComponents(i) = .FALSE.
             END DO
           END IF
         END IF
       END IF

       TransposePresent = TransposePresent .OR. ASSOCIATED(Atmp % Child)
       IF( TransposePresent ) THEN
         CALL Info(Caller,'Transpose matrix is present')
       END IF

       ! If the projector is of type x_s=P*x_m then generate a constraint matrix
       ! of type [D-P]x=0 where D is diagonal unit matrix. 
       CreateSelf = ( Atmp % ProjectorType == PROJECTOR_TYPE_NODAL ) 
       
       IF( SumThis .AND. CreateSelf ) THEN
         CALL Fatal(Caller,'It is impossible to sum up nodal projectors!')
       END IF

       ! Assume the mortar matrices refer to unordered mesh dofs
       ! and existing ConstraintMatrix to already ordered entities. 
       Reorder = ThisIsMortar
       
       ComplexSumRow = ListGetLogical( Solver % Values,'Complex Sum Row ', Found )
       IF(.NOT. Found ) THEN       
         ComplexSumRow = ( dofs == 2 .AND. ComplexMatrix .AND. .NOT. CreateSelf .AND. &
             SumThis .AND. .NOT. (ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) )
       END IF
         
       IF( Dofs == 1 ) THEN         

         IF( .NOT. ActiveComponents(1) ) THEN
           CALL Info(Caller,'Skipping component: '//TRIM(I2S(1)),Level=12)
           CYCLE
         END IF

         ! Number the rows. 
         IF( SumThis ) THEN
           DO i=1,Atmp % NumberOfRows                               
             ! Skip empty row
             IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE 

             ! If the mortar boundary is not active at this round don't apply it
             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Active ) ) THEN
                 IF( .NOT. MortarBC % Active(i) ) CYCLE
               END IF
             END IF
             
             ! Use InvPerm if it is present
             IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
               k = Atmp % InvPerm(i)
               ! Node does not have an active dof to be constrained
               IF( k == 0 ) CYCLE
             ELSE
               k = i
             END IF

             kk = k             
             IF( Reorder ) THEN
               kk = Perm(k)
               IF( kk == 0 ) CYCLE
             END IF
             
             NewRow = ( SumPerm(kk) == 0 )
             IF( NewRow ) THEN
               sumrow = sumrow + 1                
               SumPerm(kk) = sumrow 
             ELSE IF(.NOT. AllocationsDone ) THEN
               IF( Priority /= PrevPriority .AND. SumPerm(kk) < 0 ) THEN
                 NeglectedRows = NeglectedRows + 1
               ELSE
                 EliminatedRows = EliminatedRows + 1
               END IF
             END IF
           END DO
         END IF
         
         IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
           IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN                   
             k = MAXVAL( Atmp % Cols )
             ALLOCATE( MortarBC % Perm(k) )
             MortarBC % Perm = 0
             DO k=1,SIZE(Atmp % InvPerm )
               j = Atmp % InvPerm(k)
               MortarBC % Perm( j ) = k
             END DO
           END IF
         END IF
         
         
         DO i=1,Atmp % NumberOfRows           

           IF( Atmp % Rows(i) >= Atmp % Rows(i+1) ) CYCLE ! skip empty rows

           ! If the mortar boundary is not active at this round don't apply it
           IF( ThisIsMortar ) THEN
             IF( ASSOCIATED( MortarBC % Active ) ) THEN
               IF( .NOT. MortarBC % Active(i) ) CYCLE
             END IF
           END IF
             
           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF
            
           kk = k
           IF( Reorder ) THEN
             kk = Perm(k) 
             IF( kk == 0 ) CYCLE
           END IF
             
           IF( SumThis ) THEN             
             row = SumPerm(kk)
               
             ! Mark this for future contributions so we know this is already set
             ! and can skip this above.
             IF( AnyPriority ) THEN
               IF( row < 0 ) CYCLE
               IF( Priority /= PrevPriority ) SumPerm(kk) = -SumPerm(kk)
             END IF
             
             IF( row <= 0 ) THEN
               CALL Fatal(Caller,'Invalid row index: '//TRIM(I2S(row)))
             END IF
           ELSE
             sumrow = sumrow + 1
             row = sumrow
           END IF

           IF( AllocationsDone ) THEN
             Btmp % InvPerm(row) = rowoffset + kk
           END IF

           
           wsum = 0.0_dp
           

           valsum = 0.0_dp
           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1             
             valsum = valsum + ABS( Atmp % Values(l) ) 
           END DO
             

           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1
             
             col = Atmp % Cols(l) 
             val = Atmp % Values(l)

             IF( ABS( val ) < EpsVal * valsum ) CYCLE

             
             IF( Reorder ) THEN
               IF( col <= permsize ) THEN
                 col2 = Perm(col)
                 IF( col2 == 0 ) CYCLE
               ELSE
                 CALL Fatal(Caller,'col index too large: '//TRIM(I2S(col)))
               END IF
             ELSE
               col2 = col
             END IF
               
             IF( AllocationsDone ) THEN
               ! By Default there is no scaling
               Scale = 1.0_dp
               IF( ThisIsMortar ) THEN
                 IF( CreateSelf ) THEN
                   ! We want to create [D-P] hence the negative sign
                   Scale = MortarBC % MasterScale
                   wsum = wsum + val
                 ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                   ! Look if the component refers to the slave
                   IF( MortarBC % Perm( col ) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                     wsum = wsum + val
                   ELSE
                     Scale = MortarBC % MasterScale
                   END IF
                 ELSE
                   wsum = wsum + val
                 END IF

                 ! If we sum up to anti-periodic dof then use different sign
                 ! - except if the target is also antiperiodic.
                 IF( PerFlipActive ) THEN
                   IF( XOR( PerFlip(col),PerFlip(k) ) ) Scale = -Scale
                 END IF
                 
               END IF

               ! Add a new column index to the summed up row               
               ! At the first sweep we need to find the first unset position
               IF( SumThis ) THEN
                 k2 = Btmp % Rows(row)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF
               
               Btmp % Cols(k2) = col2
               Btmp % Values(k2) = Scale * val
               IF(ASSOCIATED(Btmp % TValues)) THEN
                 IF(ASSOCIATED(Atmp % Child)) THEN
                   Btmp % TValues(k2) = Scale * Atmp % Child % Values(l)
                 ELSE
                   Btmp % TValues(k2) = Scale * val
                 END IF
               END IF
             ELSE
               k2 = k2 + 1
               IF( SumThis ) THEN
                 SumCount(row) = SumCount(row) + 1
               END IF
             END IF
           END DO
           
           ! Add the self entry as in 'D'
           IF( CreateSelf ) THEN
             k2 = k2 + 1
             IF( AllocationsDone ) THEN
               Btmp % Cols(k2) = Perm( Atmp % InvPerm(i) )
               Btmp % Values(k2) = MortarBC % SlaveScale * wsum
             ELSE               
               IF( SumThis) SumCount(row) = SumCount(row) + 1
             END IF
           END IF
           
           ! Add a diagonal entry if requested. When this is done at the final stage
           ! all the hassle with the right column index is easier.
           IF( ThisIsMortar ) THEN
             diag: IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
               IF( .NOT. HaveMortarDiag ) THEN
                 MortarDiag = MortarBC % Diag(i)
                 LumpedDiag = MortarBC % LumpedDiag
               END IF
              
               IF( LumpedDiag ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = row + arows 
                   ! The factor 0.5 comes from the fact that the 
                   ! contribution is summed twice, 2nd time as transpose
                   ! For Nodal projector the entry is 1/(weight*coeff)
                   ! For Galerkin projector the is weight/coeff 
                   Btmp % Values(k2) = Btmp % Values(k2) - 0.5_dp * MortarDiag * wsum
                 ELSE
                   IF( SumThis) SumCount(row) = SumCount(row) + 1
                 END IF
               ELSE
                 IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN                   
                   CALL Fatal(Caller,'MortarBC % Perm required, try lumped')
                 END IF
                 
                 DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                   col = Atmp % Cols(l) 

                   IF( col > permsize ) THEN
                     PRINT *,'col too large',col,permsize
                     CYCLE
                   END IF
                   col2 = Perm(col)
                   IF( col2 == 0 ) CYCLE
                     
                   IF( CreateSelf ) THEN
                     Scale = -MortarBC % MasterScale
                   ELSE
                     IF( MortarBC % Perm( col ) > 0 ) THEN
                       Scale = MortarBC % SlaveScale 
                     ELSE
                       CYCLE                     
                     END IF
                   END IF
                   
                   k2 = k2 + 1
                   IF( AllocationsDone ) THEN                                        
                     IF( SumThis ) THEN
                       l2 = ABS( SumPerm( col2) )
                     ELSE
                       l2 = MortarBC % Perm(col)
                     END IF
                     
                     Btmp % Cols(k2) = l2 + arows + rowoffset
                     Btmp % Values(k2) = Btmp % Values(k2) - 0.5_dp * val * MortarDiag
                   ELSE
                     IF( SumThis) SumCount(row) = SumCount(row) + 1
                   END IF
                 END DO
               END IF
             END IF diag
           END IF

           IF( AllocationsDone ) THEN
             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                 Btmp % Rhs(row) = Btmp % Rhs(row) + wsum * MortarBC % rhs(i)
               END IF
             END IF

             ! If every component is uniquely summed we can compute the row indexes simply
             IF( .NOT. SumThis ) THEN
               Btmp % Rows(row+1) = k2 + 1
             END IF
           END IF
         END DO
         
       ELSE IF( ComplexSumRow ) THEN

         CALL Info(Caller,'Using simplified complex summing!',Level=6)
         ComplexSumRow = .TRUE.
         
         ! In case of a vector valued problem create a projector that acts on all 
         ! components of the vector. Otherwise follow the same logic.
         IF( SumThis ) THEN
           DO i=1,Atmp % NumberOfRows                        
             
             IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
               k = Atmp % InvPerm(i)
               IF( k == 0 ) CYCLE
             ELSE
               k = i
             END IF
             
             kk = k
             IF( Reorder ) THEN
               kk = Perm(k)
               IF( kk == 0 ) CYCLE
             END IF
             
             NewRow = ( SumPerm(kk) == 0 )
             IF( NewRow ) THEN
               sumrow = sumrow + 1                
               SumPerm(kk) = sumrow 
             ELSE IF(.NOT. AllocationsDone ) THEN
               EliminatedRows = EliminatedRows + 1
             END IF
           END DO
         END IF
           
         
         DO i=1,Atmp % NumberOfRows           
           
           IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
             k = Atmp % InvPerm(i)
             IF( k == 0 ) CYCLE
           ELSE
             k = i
           END IF
            
           kk = k
           IF( Reorder ) THEN
             kk = Perm(k) 
             IF( kk == 0 ) CYCLE
           END IF
             
           IF( SumThis ) THEN             
             row = SumPerm(kk)
           ELSE
             sumrow = sumrow + 1
             row = sumrow
           END IF

           ! For complex matrices 
           IF( AllocationsDone ) THEN
             Btmp % InvPerm(2*row-1) = rowoffset + 2*(kk-1)+1
             Btmp % InvPerm(2*row) = rowoffset + 2*kk
           END IF

           wsum = 0.0_dp
                        

           DO l=Atmp % Rows(i),Atmp % Rows(i+1)-1
             
             col = Atmp % Cols(l) 
             val = Atmp % Values(l)
             
             IF( Reorder ) THEN
               col2 = Perm(col)
               IF( col2 == 0 ) CYCLE
             ELSE
               col2 = col
             END IF
               
             IF( AllocationsDone ) THEN
               ! By Default there is no scaling
               Scale = 1.0_dp
               IF( ThisIsMortar ) THEN
                 IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                   ! Look if the component refers to the slave
                   IF( MortarBC % Perm( col ) > 0 ) THEN
                     Scale = MortarBC % SlaveScale 
                     wsum = wsum + val
                   ELSE
                     Scale = MortarBC % MasterScale
                   END IF
                 ELSE
                   wsum = wsum + val
                 END IF
                 
                 ! If we sum up to anti-periodic dof then use different sign
                 ! - except if the target is also antiperiodic.
                 IF( PerFlipActive ) THEN
                   IF( XOR( PerFlip(col),PerFlip(k) ) ) Scale = -Scale
                 END IF
                 
               END IF

               ! Add a new column index to the summed up row               
               ! At the first sweep we need to find the first unset position
               ! Real part
               IF( SumThis ) THEN
                 k2 = Btmp % Rows(2*row-1)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF
                                            
               Btmp % Cols(k2) = 2 * col2 - 1
               Btmp % Values(k2) = Scale * val

               k2 = k2 + 1
               Btmp % Cols(k2) = 2 * col2
               Btmp % Values(k2) = 0.0

               ! Complex part
               IF( SumThis ) THEN
                 k2 = Btmp % Rows(2*row)
                 DO WHILE( Btmp % Cols(k2) > 0 )
                   k2 = k2 + 1
                 END DO
               ELSE
                 k2 = k2 + 1
               END IF

               Btmp % Cols(k2) = 2 * col2 - 1 
               Btmp % Values(k2) = 0.0
             
               k2 = k2 + 1
               Btmp % Cols(k2) = 2 * col2 
               Btmp % Values(k2) = Scale * val
             ELSE
               k2 = k2 + 4
               IF( SumThis ) THEN
                 SumCount(2*row-1) = SumCount(2*row-1) + 2
                 SumCount(2*row) = SumCount(2*row) + 2
               END IF
             END IF
           END DO
           
           IF( AllocationsDone ) THEN
             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                 Btmp % Rhs(2*row-1) = Btmp % Rhs(2*row-1) + wsum * MortarBC % rhs(i)
               END IF
             END IF
           END IF
         END DO
         
       ELSE
         
         ! dofs > 1
         ! In case of a vector valued problem create a projector that acts on all 
         ! components of the vector. Otherwise follow the same logic.
         DO i=1,Atmp % NumberOfRows           
           DO j=1,Dofs
             
             IF( .NOT. ActiveComponents(j) ) THEN
               CALL Info(Caller,'Skipping component: '//TRIM(I2S(j)),Level=12)
               CYCLE
             END IF
             
             ! For complex matrices both entries mist be created
             ! since preconditioning benefits from 
             IF( ComplexMatrix ) THEN
               IF( MODULO( j, 2 ) == 0 ) THEN
                 j2 = j-1
               ELSE 
                 j2 = j+1
               END IF
             ELSE
               j2 = 0
             END IF

             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Active ) ) THEN
                 IF( .NOT. MortarBC % Active(Dofs*(i-1)+j) ) CYCLE
               END IF
             END IF

             IF( ASSOCIATED( Atmp % InvPerm ) ) THEN
               k = Atmp % InvPerm(i)
               IF( k == 0 ) CYCLE
             ELSE
               k = i
             END IF

             kk = k
             IF( Reorder ) THEN
               kk = Perm(k)
               IF( kk == 0 ) CYCLE
             END IF

             IF( SumThis ) THEN
               IF( Dofs*(k-1)+j > SIZE(SumPerm) ) THEN
                 CALL Fatal(Caller,'Index out of range')
               END IF
               NewRow = ( SumPerm(Dofs*(kk-1)+j) == 0 )
               IF( NewRow ) THEN
                 sumrow = sumrow + 1                
                 IF( Priority /= 0 ) THEN
                   ! Use negative sign to show that this has already been set by priority
                   SumPerm(Dofs*(kk-1)+j) = -sumrow 
                 ELSE
                   SumPerm(Dofs*(kk-1)+j) = sumrow 
                 END IF
               ELSE IF( Priority /= PrevPriority .AND. SumPerm(Dofs*(kk-1)+j) < 0 ) THEN
                 IF(.NOT. AllocationsDone ) THEN
                   NeglectedRows = NeglectedRows + 1
                 END IF                 
                 CYCLE
               ELSE
                 IF(.NOT. AllocationsDone ) THEN
                   EliminatedRows = EliminatedRows + 1
                 END IF
               END IF
               row = ABS( SumPerm(Dofs*(kk-1)+j) )
             ELSE
               sumrow = sumrow + 1
               row = sumrow
             END IF

             IF( AllocationsDone ) THEN
               Btmp % InvPerm(row) = rowoffset + Dofs * ( kk - 1 ) + j
             END IF

             
             wsum = 0.0_dp

             DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1             

               col = Atmp % Cols(k)                

               IF( Reorder ) THEN                 
                 IF( col <= permsize ) THEN
                   col2 = Perm(col)
                   IF( col2 == 0 ) CYCLE
                 ELSE 
                   PRINT *,'col too large',col,permsize
                   CYCLE
                 END IF
               ELSE
                 col2 = col
               END IF

                 
               k2 = k2 + 1
               
               IF( AllocationsDone ) THEN
                 Scale = 1.0_dp
                 IF( ThisIsMortar ) THEN
                   IF( CreateSelf ) THEN
                     Scale = MortarBC % MasterScale
                     wsum = wsum + Atmp % Values(k)
                   ELSE IF( ASSOCIATED( MortarBC % Perm ) ) THEN
                     IF( MortarBC % Perm(col) > 0 ) THEN
                       Scale = MortarBC % SlaveScale 
                       wsum = wsum + Atmp % Values(k) 
                     ELSE
                       Scale = MortarBC % MasterScale
                     END IF
                   END IF

                   ! If we sum up to anti-periodic dof then use different sign
                   ! - except if the target is also antiperiodic.
                   IF( PerFlipActive ) THEN
                     IF( XOR( PerFlip(col),PerFlip(k) ) ) Scale = -Scale
                   END IF

                 END IF
                 
                 Btmp % Cols(k2) = Dofs * ( col2 - 1) + j
                 Btmp % Values(k2) = Scale * Atmp % Values(k)
                 IF(ASSOCIATED(Btmp % Tvalues)) THEN
                   IF(ASSOCIATED(Atmp % Child)) THEN
                     Btmp % TValues(k2) = Scale * Atmp % Child % Values(k)
                   ELSE
                     Btmp % TValues(k2) = Scale * Atmp % Values(k)
                   END IF
                 END IF
               ELSE
                 IF( SumThis ) THEN
                   SumCount(row) = SumCount(row) + 1
                 END IF                 
               END IF
             END DO
             
             ! Add the self entry as in 'D'
             IF( CreateSelf ) THEN
               k2 = k2 + 1
               IF( AllocationsDone ) THEN
                 Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j
                 Btmp % Values(k2) = MortarBC % SlaveScale * wsum
               END IF
             END IF
             
             ! Create the imaginary part (real part) corresponding to the 
             ! real part (imaginary part) of the projector. 
             IF( j2 /= 0 ) THEN
               DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1             

                 col = Atmp % Cols(k)                

                 IF( Reorder ) THEN
                   IF( col <= permsize ) THEN
                     col2 = Perm(col)
                     IF( col2 == 0 ) CYCLE
                   END IF
                 ELSE
                   col2 = col
                 END IF

                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = Dofs * ( col2 - 1) + j2
                 ELSE
                   IF( SumThis ) THEN
                     SumCount(row) = SumCount(row) + 1
                   END IF
                 END IF
               END DO

               IF( CreateSelf ) THEN
                 k2 = k2 + 1
                 IF( AllocationsDone ) THEN
                   Btmp % Cols(k2) = Dofs * ( Perm( Atmp % InvPerm(i) ) -1 ) + j2
                 END IF
               END IF
             END IF


             IF( ThisIsMortar ) THEN
               IF( ASSOCIATED( MortarBC % Diag ) .OR. HaveMortarDiag ) THEN
                 IF( .NOT. HaveMortarDiag ) THEN
                   MortarDiag = MortarBC % Diag(Dofs*(i-1)+j)
                   LumpedDiag = MortarBC % LumpedDiag
                 END IF

                 IF( LumpedDiag ) THEN
                   k2 = k2 + 1
                   IF( AllocationsDone ) THEN
                     Btmp % Cols(k2) = row + arows
                     Btmp % Values(k2) = -0.5_dp * wsum * MortarDiag
                   END IF
                 ELSE
                   DO k=Atmp % Rows(i),Atmp % Rows(i+1)-1                 
                     col = Atmp % Cols(k) 

                     IF( col > permsize ) CYCLE
                     col2 = Perm(col)

                     IF( CreateSelf ) THEN
                       Scale = -MortarBC % MasterScale
                     ELSE 
                       IF( MortarBC % Perm( col ) > 0 ) THEN
                         Scale = MortarBC % SlaveScale 
                       ELSE
                         CYCLE                     
                       END IF
                     END IF

                     k2 = k2 + 1
                     IF( AllocationsDone ) THEN                   
                       Btmp % Cols(k2) = Dofs*(MortarBC % Perm( col )-1)+j + arows + rowoffset
                       Btmp % Values(k2) = -0.5_dp * Atmp % Values(k) * MortarDiag
                     END IF
                   END DO
                 END IF
               END IF
             END IF

               
             IF( AllocationsDone ) THEN
               IF( ThisIsMortar ) THEN
                 IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
                   Btmp % Rhs(row) = wsum * MortarBC % rhs(Dofs*(i-1)+j)
                 END IF
               END IF
               IF(.NOT. SumThis ) THEN
                 Btmp % Rows(row+1) = k2 + 1
               END IF
             END IF

           END DO
         END DO
       END IF ! dofs > 1
       
       IF( .NOT. SumThis ) THEN
         rowoffset = rowoffset + Arows
         IF( SumProjectors ) THEN
           CALL Info(Caller,'Not summed up size is ' &
           //TRIM(I2S(sumrow))//' rows and '//TRIM(I2S(k2))//' nonzeros',Level=8)
           sumrow0 = sumrow
           k20 = k2
         END IF
       END IF
         
       PrevPriority = Priority 
     END DO ! constrain_ind

     IF( k2 == 0 ) THEN
       CALL Info(Caller,'No entries in constraint matrix!',Level=6)
!      Solver % Matrix % ConstraintMatrix => NULL()
       RETURN
     END IF

     ! Allocate the united matrix of all the boundary matrices
     !-------------------------------------------------------
     IF( .NOT. AllocationsDone ) THEN
       CALL Info(Caller,'Allocating '//&
           TRIM(I2S(sumrow))//' rows and '//TRIM(I2S(k2))//' nonzeros',&
           Level=6)

       IF( ComplexSumRow ) THEN
         sumrow = 2 * sumrow
       END IF
       
       Btmp => AllocateMatrix()
       ALLOCATE( Btmp % RHS(sumrow), Btmp % Rows(sumrow+1), &
           Btmp % Cols(k2), Btmp % Values(k2), &
           Btmp % InvPerm(sumrow) )

       Btmp % Rhs = 0.0_dp
       Btmp % Rows = 0
       Btmp % Cols = 0
       Btmp % Values = 0.0_dp
       Btmp % NumberOFRows = sumrow 
       Btmp % InvPerm = 0
       Btmp % Rows(1) = 1

       IF(TransposePresent) THEN
         ALLOCATE(Btmp % TValues(k2))
         Btmp % Tvalues = 0._dp
       END IF

       IF( SumProjectors ) THEN
         Btmp % Rows(sumrow0+1) = k20+1 
         DO i=sumrow0+2,sumrow+1
           Btmp % Rows(i) = Btmp % Rows(i-1) + SumCount(i-1)
         END DO
         SumPerm = 0
         DEALLOCATE( SumCount ) 
       END IF

       AllocationsDone = .TRUE.

       GOTO 100
     END IF
     
     CALL Info(Caller,'Used '//TRIM(I2S(sumrow))//&
         ' rows and '//TRIM(I2S(k2))//' nonzeros',Level=6)
          
     ! Eliminate entries
     IF( SumProjectors ) THEN
       CALL Info(Caller,'Number of eliminated rows: '//TRIM(I2S(EliminatedRows)))
       IF( EliminatedRows > 0 ) CALL CRS_PackMatrix( Btmp ) 
     END IF

     IF( NeglectedRows > 0 ) THEN
       CALL Info(Caller,'Number of neglected rows: '//TRIM(I2S(NeglectedRows)))
     END IF
        
     Solver % Matrix % ConstraintMatrix => Btmp     
     Solver % MortarBCsChanged = .FALSE.
     
     CALL Info(Caller,'Finished creating constraint matrix',Level=12)

   END SUBROUTINE GenerateConstraintMatrix
     

   SUBROUTINE ReleaseConstraintMatrix(Solver) 
     TYPE(Solver_t) :: Solver

     CALL FreeMatrix(Solver % Matrix % ConstraintMatrix)
     Solver % Matrix % ConstraintMatrix => NULL()

   END SUBROUTINE ReleaseConstraintMatrix


   SUBROUTINE ReleaseProjectors(Model, Solver) 

     TYPE(Model_t) :: Model
     TYPE(Solver_t) :: Solver

     TYPE(ValueList_t), POINTER :: BC
     TYPE(Matrix_t), POINTER :: Projector
     INTEGER :: i
     

     IF( .NOT. ASSOCIATED( Solver % MortarBCs ) ) RETURN

     DO i=1,Model % NumberOFBCs
       BC => Model % BCs(i) % Values
       Projector => Solver % MortarBCs(i) % Projector 
       IF( ASSOCIATED( Projector ) ) THEN
         IF( ASSOCIATED( Projector % EMatrix ) ) THEN
           CALL FreeMatrix( Projector % Ematrix ) 
         END IF
         CALL FreeMatrix( Projector )
         Solver % MortarBCs(i) % Projector => NULL()
       END IF
     END DO

   END SUBROUTINE ReleaseProjectors


   !> Defines and potentially creates output directory.
   !> The output directory may given in different ways, and even be part of the
   !> filename, or be relative to home directory. We try to parse every possible
   !> scenario here that user might have in mind.
   !-----------------------------------------------------------------------------
   SUBROUTINE SolverOutputDirectory( Solver, Filename, OutputDirectory, &
       MakeDir, UseMeshDir  )

     TYPE(Solver_t) :: Solver
     CHARACTER(LEN=MAX_NAME_LEN) :: Filename, OutputDirectory
     LOGICAL, OPTIONAL :: MakeDir, UseMeshDir

     LOGICAL :: Found, AbsPathInName, DoDir, PartitioningSubDir
     INTEGER :: nd, nf, n
     CHARACTER(LEN=MAX_NAME_LEN) :: Str

     IF( PRESENT( MakeDir ) ) THEN
       DoDir = MakeDir
     ELSE
       DoDir = ( Solver % TimesVisited == 0 ) .AND. ( ParEnv % MyPe == 0 )
     END IF

     ! Output directory is obtained in order
     ! 1) solver section
     ! 2) simulation section
     ! 3) header section
     OutputDirectory = ListGetString( Solver % Values,'Output Directory',Found) 
     IF(.NOT. Found) OutputDirectory = ListGetString( CurrentModel % Simulation,&
         'Output Directory',Found) 
     IF(.NOT. Found) OutputDirectory = TRIM(OutputPath)          
     nd = LEN_TRIM(OutputDirectory)

     ! If the path is just working directory then that is not an excude
     ! to not use the mesh name, or directory that comes with the filename 
     IF(.NOT. Found .AND. nd == 1 .AND. OutputDirectory(1:1)=='.') nd = 0

     ! If requested by the optional parameter use the mesh directory when
     ! no results directory given. This is an old convection used in some solvers. 
     IF( nd == 0 .AND. PRESENT( UseMeshDir ) ) THEN
       IF( UseMeshDir ) THEN
         OutputDirectory = TRIM(CurrentModel % Mesh % Name)
         nd = LEN_TRIM(OutputDirectory)       
       END IF
     END IF
     
     ! Use may have given part or all of the path in the filename.
     ! This is not preferred, but we cannot trust the user.
     nf = LEN_TRIM(Filename)        
     n = INDEX(Filename(1:nf),'/')
     AbsPathInName = INDEX(FileName,':')>0 .OR. (Filename(1:1)=='/') &
         .OR. (Filename(1:1)==Backslash)

     IF( nd > 0 .AND. .NOT. AbsPathInName ) THEN
       ! Check that we have not given the path relative to home directory
       ! because the code does not understand the meaning of tilde.
       IF( OutputDirectory(1:2) == '~/') THEN
         CALL GETENV('HOME',Str)
         OutputDirectory = TRIM(Str)//'/'//OutputDirectory(3:nd)
         nd = LEN_TRIM(OutputDirectory)
       END IF
       ! To be on the safe side create the directory. If it already exists no harm done.
       ! Note that only one direcory may be created. Hence if there is a path with many subdirectories
       ! that will be a problem. Fortran does not have a standard ENQUIRE for directories hence
       ! we just try to make it. 
       IF( DoDir ) THEN
         CALL Info('SolverOutputDirectory','Creating directory: '//TRIM(OutputDirectory(1:nd)),Level=8)
         CALL MakeDirectory( OutputDirectory(1:nd) // CHAR(0) )      
       END IF
     END IF

     ! In this case the filename includes also path and we remove it from there and
     ! add it to the directory. 
     IF( n > 2 ) THEN    
       CALL Info('SolverOutputDirectory','Parcing path from filename: '//TRIM(Filename(1:n)),Level=10)
       IF( AbsPathInName .OR. nd == 0) THEN
         ! If the path is absolute then it overruns the given path!
         OutputDirectory = Filename(1:n-1)
         nd = n-1
       ELSE
         ! If path is relative we add it to the OutputDirectory and take it away from Filename
         OutputDirectory = OutputDirectory(1:nd)//'/'//Filename(1:n-1)        
         nd = nd + n 
       END IF
       Filename = Filename(n+1:nf)      

       IF( DoDir ) THEN
         CALL Info('SolverOutputDirectory','Creating directory: '//TRIM(OutputDirectory(1:nd)),Level=8)
         CALL MakeDirectory( OutputDirectory(1:nd) // CHAR(0) )
       END IF
     END IF

     ! Finally, on request save each partitioning to different directory.
     PartitioningSubDir = ListGetLogical( Solver % Values,'Output Partitioning Directory',Found)
     IF(.NOT. Found ) THEN
       PartitioningSubDir = ListGetLogical( CurrentModel % Simulation,'Output Partitioning Directory',Found)
     END IF
     IF( PartitioningSubDir ) THEN
       OutputDirectory = TRIM(OutputDirectory)//'/np'//TRIM(I2S(ParEnv % PEs))
       nd = LEN_TRIM(OutputDirectory)             
       IF( DoDir ) THEN
         CALL Info('SolverOutputDirectory','Creating directory: '//TRIM(OutputDirectory(1:nd)),Level=8)
         CALL MakeDirectory( OutputDirectory(1:nd) // CHAR(0) )
       END IF
      END IF
       

     
     !PRINT *,'Filename:',TRIM(Filename)
     !PRINT *,'OutputDirectory:',TRIM(OutputDirectory) 

   END SUBROUTINE SolverOutputDirectory
   !-----------------------------------------------------------------------------
 

END MODULE SolverUtils

!> \}
