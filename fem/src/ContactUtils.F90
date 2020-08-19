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
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/

!> Utilities for nonequality conditions (limiteres & contact) used by
!> individual solvers. 
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{


MODULE ContactUtils

   USE Types
!   USE LoadMod
!   USE DirectSolve
!   USE Multigrid
!   USE IterSolve
!   USE ElementUtils
!   USE ComponentUtils
!   USE TimeIntegrate
!   USE ModelDescription
   USE MeshUtils
   USE AssemblyUtils
   USE ScalingUtils
!   USE SolverUtils
!   USE MortarProjector
!   USE ParallelUtils
!   USE ParallelEigenSolve
   USE ListMatrix
   USE CRSMatrix
   
CONTAINS

   
!------------------------------------------------------------------------------
!> Determine soft limiters set. This is called after the solution.
!> and can therefore be active only on the 2nd nonlinear iteration round.
!------------------------------------------------------------------------------
   SUBROUTINE DetermineSoftLimiter( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(variable_t), POINTER :: Var, LoadVar, IterV, LimitVar
     TYPE(Element_t), POINTER :: Element
     INTEGER :: i,j,k,n,t,ind,dofs, dof, bf, bc, Upper, Removed, Added, &
         ElemFirst, ElemLast, totsize, i2, j2, ind2
     REAL(KIND=dp), POINTER :: FieldValues(:), LoadValues(:), &
         ElemLimit(:),ElemInit(:), ElemActive(:)
     REAL(KIND=dp) :: LimitSign, EqSign, ValEps, LoadEps, val
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF, GotInit, GotActive
     LOGICAL, ALLOCATABLE :: LimitDone(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params, Entity
     CHARACTER(LEN=MAX_NAME_LEN) :: Name, LimitName, InitName, ActiveName
     LOGICAL, ALLOCATABLE :: InterfaceDof(:)
     INTEGER :: ConservativeAfterIters, NonlinIter, CoupledIter, DownStreamDirection
     LOGICAL :: Conservative, ConservativeAdd, ConservativeRemove, &
         DoAdd, DoRemove, DirectionActive, FirstTime, DownStreamRemove
     TYPE(Mesh_t), POINTER :: Mesh
     CHARACTER(*), PARAMETER :: Caller = 'DetermineSoftLimiter'
     
     Model => CurrentModel
     Var => Solver % Variable
     Mesh => Solver % Mesh
     

     ! Check the iterations counts and determine whether this is the first 
     ! time with this solver. 
     !------------------------------------------------------------------------
     FirstTime = .TRUE.
     iterV => VariableGet( Mesh % Variables,'nonlin iter')
     IF( ASSOCIATED( iterV ) ) THEN
       NonlinIter =  NINT( iterV % Values(1) ) 
       IF( NonlinIter > 1 ) FirstTime = .FALSE.
     END IF

     iterV => VariableGet( Mesh % Variables,'coupled iter')
     IF( ASSOCIATED( iterV ) ) THEN
       CoupledIter = NINT( iterV % Values(1) )
       IF( CoupledIter > 1 ) FirstTime = .FALSE.
     END IF
          
     ! Determine variable for computing the contact load used to determine the 
     ! soft limit set.
     !------------------------------------------------------------------------
     CALL Info(Caller,'Determining soft limiter problems',Level=8)
     LoadVar => VariableGet( Model % Variables, &
         GetVarName(Var) // ' Contact Load',ThisOnly = .TRUE. )
     CALL CalculateLoads( Solver, Solver % Matrix, Var % Values, Var % DOFs, .FALSE., LoadVar ) 

     IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
       CALL Fatal(Caller, &
           'No Loads associated with variable '//GetVarName(Var) )
       RETURN
     END IF
     LoadValues => LoadVar % Values


     ! The variable to be constrained by the soft limiters
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     Params => Solver % Values

     ConservativeAdd = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Add After Iterations',Conservative ) 
     IF( Conservative ) THEN
       ConservativeAdd = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeAdd ) THEN
         CALL Info(Caller,'Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeRemove = .FALSE.
     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',Found )      
     IF( Found ) THEN
       Conservative = .TRUE.  
       ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info(Caller,'Removing dofs in conservative fashion',Level=8)
       END IF
     END IF

     DownStreamRemove = ListGetLogical( Params,'Apply Limiter Remove Downstream',Found)
     IF( DownStreamRemove ) THEN
       CALL Info(Caller,'Removing contact dofs only in downstream',Level=8)      
       ConservativeRemove = .TRUE.
       Conservative = .TRUE.
       DownStreamDirection = ListGetInteger( Params,'Apply Limiter Downstream Direction',Found)
       IF(.NOT. Found ) DownStreamDirection = 1
     END IF
       
     LoadEps = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps = EPSILON( LoadEps )
         
     ValEps = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps = EPSILON( ValEps )

     ! The user may want to toggle the sign for various kinds of equations
     ! The default sign that come from standard formulation of Laplace equation.
     !---------------------------------------------------------------------------       
     IF( ListGetLogical( Params,'Limiter Load Sign Negative',Found) ) THEN
       EqSign = -1.0_dp
     ELSE
       EqSign = 1.0_dp
     END IF

     ! Loop through upper and lower limits     
     !------------------------------------------------------------------------
     DO Upper=0,1

       DirectionActive = .FALSE.

       ! If we have both upper and lower limiter then these logical vectors need to be 
       ! reinitialized for the 2nd sweep.
       IF( ALLOCATED( LimitDone) ) LimitDone = .FALSE.
       IF( ALLOCATED( InterfaceDof ) ) InterfaceDof = .FALSE. 

       ! Upper and lower limits have different sign for testing
       !----------------------------------------------------------------------       
       IF( Upper == 0 ) THEN
         LimitSign = -EqSign
       ELSE
         LimitSign = EqSign
       END IF       
       
       ! Go through the components of the field, if many
       !-------------------------------------------------
       DO DOF = 1,dofs

         name = Var % name
         IF ( Var % DOFs > 1 ) name = ComponentName(name,DOF)

         ! The keywords for the correct lower or upper limit of the variable
         !------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           LimitName = TRIM(name)//' Lower Limit'           
           InitName = TRIM(name)//' Lower Initial'
           ActiveName = TRIM(name)//' Lower Active'
         ELSE
           LimitName = TRIM(name)//' Upper Limit' 
           InitName = TRIM(name)//' Upper Initial' 
           ActiveName = TRIM(name)//' Upper Active' 
         END IF

         AnyLimitBC = ListCheckPresentAnyBC( Model, LimitName )
         AnyLimitBF = ListCheckPresentAnyBodyForce( Model, LimitName )

         ! If there is no active keyword then there really is nothing to do
         !----------------------------------------------------------------
         IF( .NOT. ( AnyLimitBC .OR. AnyLimitBF ) ) CYCLE
         DirectionActive = .TRUE.
         
         CALL Info(Caller,'Applying limit: '//TRIM(LimitName),Level=8)

         ! OK: Do contact for a particular dof and only upper or lower limit
         !------------------------------------------------------------------------

         ! Define the range of elements for which the limiters are active
         !---------------------------------------------------------------
         ElemFirst = Model % NumberOfBulkElements + 1           
         ElemLast = Model % NumberOfBulkElements 
        
         IF( AnyLimitBF ) ElemFirst = 1
         IF( AnyLimitBC ) ElemLast = Model % NumberOfBulkElements + &
             Model % NumberOfBoundaryElements 
         
         IF(.NOT. ALLOCATED( LimitDone) ) THEN
           n = Model % MaxElementNodes
           ALLOCATE( LimitDone( totsize ), ElemLimit(n), ElemInit(n), ElemActive(n) )
           LimitDone = .FALSE.
         END IF

         ! Check that active set vectors for limiters exist, otherwise allocate
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           IF( .NOT. ASSOCIATED(Var % LowerLimitActive ) ) THEN
             ALLOCATE( Var % LowerLimitActive( totsize ) )
             Var % LowerLimitActive = .FALSE.
           END IF
           LimitActive => Var % LowerLimitActive
         ELSE
           IF( .NOT. ASSOCIATED( Var % UpperLimitActive ) ) THEN
             ALLOCATE( Var % UpperLimitActive( totsize ) )
             Var % UpperLimitActive = .FALSE.
           END IF
           LimitActive => Var % UpperLimitActive
         END IF
 
         Removed = 0
         Added = 0        
         IF(.NOT. ALLOCATED( LimitDone) ) THEN
           n = Model % MaxElementNodes
           ALLOCATE( LimitDone( totsize ), ElemLimit(n), ElemInit(n), ElemActive(n) )
           LimitDone = .FALSE.
         END IF


         IF( FirstTime ) THEN
           ! In the first time set the initial set 
           !----------------------------------------------------------------------
           DO t = ElemFirst, ElemLast
             
             Element => Model % Elements(t)
             Model % CurrentElement => Element

             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             Found = .FALSE.
             IF( t > Model % NumberOfBulkElements ) THEN
               DO bc = 1,Model % NumberOfBCs
                 IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                   Found = .TRUE.
                   Entity => Model % BCs(bc) % Values
                   EXIT
                 END IF
               END DO
               IF(.NOT. Found ) CYCLE
             ELSE             
               bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                   'Body Force', Found)
               IF(.NOT. Found ) CYCLE
               Entity => Model % BodyForces(bf) % Values
             END IF
             
             ElemLimit(1:n) = ListGetReal( Entity, &
                 LimitName, n, NodeIndexes, Found)             
             IF(.NOT. Found) CYCLE

             ElemInit(1:n) = ListGetReal( Entity, &
                 InitName, n, NodeIndexes, GotInit)
             ElemActive(1:n) = ListGetReal( Entity, &
                 ActiveName, n, NodeIndexes, GotActive)
             IF(.NOT. ( GotInit .OR. GotActive ) ) CYCLE


             DO i=1,n
               j = FieldPerm( NodeIndexes(i) )
               IF( j == 0 ) CYCLE
               ind = Dofs * ( j - 1) + Dof

               IF( LimitDone(ind) ) CYCLE
             
               ! Go through the active set and free nodes with wrong sign in contact force
               !--------------------------------------------------------------------------       
               IF( GotInit .AND. ElemInit(i) > 0.0_dp ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE.
               ELSE IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE.
               ELSE
                 LimitActive(ind) = .FALSE.
               END IF

               ! Enforce the values to limits because nonlinear material models
               ! may otherwise lead to divergence of the iteration
               !--------------------------------------------------------------
               IF( LimitActive(ind) ) THEN
                 IF( Upper == 0 ) THEN
                   Var % Values(ind) = MAX( val, ElemLimit(i) )
                 ELSE
                   Var % Values(ind) = MIN( val, ElemLimit(i) )
                 END IF
               END IF
               
               LimitDone(ind) = .TRUE.             
             END DO
           END DO

           CYCLE
         END IF


         IF( Conservative ) THEN
           IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
             ALLOCATE( InterfaceDof( totsize ) )
             InterfaceDof = .FALSE. 
           END IF

           
           ! Mark limited and unlimited neighbours and thereby make a 
           ! list of interface dofs. 
           !----------------------------------------------------------------------
           DO t = ElemFirst, ElemLast
             
             Element => Model % Elements(t)
             Model % CurrentElement => Element
             n = Element % TYPE % NumberOfNodes
             NodeIndexes => Element % NodeIndexes
             
             Found = .FALSE.
             IF( t > Model % NumberOfBulkElements ) THEN
               DO bc = 1,Model % NumberOfBCs
                 IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                   Found = .TRUE.
                   Entity => Model % BCs(bc) % Values
                   EXIT
                 END IF
               END DO
               IF(.NOT. Found ) CYCLE
             ELSE             
               bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                   'Body Force', Found)
               IF(.NOT. Found ) CYCLE
               Entity => Model % BodyForces(bf) % Values
             END IF

             ElemLimit(1:n) = ListGetReal( Entity, &
                 LimitName, n, NodeIndexes, Found)
             IF(.NOT. Found) CYCLE


             IF( DownStreamRemove ) THEN
               ! This includes only interface dofs donwstream from
               ! non-contact zone.
               BLOCK
                 REAL(kind=DP) :: r1(3),r2(3),dr(3),reps=1.0d-6
                 
                 DO i=1,n
                   j = FieldPerm( NodeIndexes(i) )
                   IF( j == 0 ) CYCLE
                   ind = Dofs * ( j - 1) + Dof
                   
                   ! Downstream of non-contact zone
                   IF(LimitActive(ind)) CYCLE
                                      
                   DO i2 = i,n
                     IF( i2 == i ) CYCLE                   
                     j2 = FieldPerm( NodeIndexes(i2) )
                     IF( j2 == 0 ) CYCLE
                     ind2 = Dofs * ( j2 - 1) + Dof
                     
                     IF( LimitActive(ind2) ) THEN
                       r2(1) =  Mesh % Nodes % x(NodeIndexes(i2))
                       r2(2) =  Mesh % Nodes % y(NodeIndexes(i2))
                       r2(3) =  Mesh % Nodes % z(NodeIndexes(i2))
                       
                       r1(1) = Mesh % Nodes % x(NodeIndexes(i))
                       r1(2) = Mesh % Nodes % y(NodeIndexes(i))
                       r1(3) = Mesh % Nodes % z(NodeIndexes(i))

                       k = DownStreamDirection 
                       IF( k > 0 ) THEN
                         dr = r2 - r1
                       ELSE
                         dr = r1 - r2
                         k = -k
                       END IF
                       
                       IF( dr(k) < reps ) CYCLE
                       
                       IF( dr(k) > 0.5*SQRT(SUM(dr*dr)) ) THEN
                         InterfaceDof(ind2) = .TRUE.
                         !PRINT *,'downstream coord:',dr
                       END IF
                     END IF
                   END DO
                 END DO
               END BLOCK
             ELSE
               ! This includes all interface dofs
               DO i=1,n
                 j = FieldPerm( NodeIndexes(i) )
                 IF( j == 0 ) CYCLE
                 ind = Dofs * ( j - 1) + Dof
                 
                 DO i2 = i+1,n
                   j2 = FieldPerm( NodeIndexes(i2) )
                   IF( j2 == 0 ) CYCLE
                   ind2 = Dofs * ( j2 - 1) + Dof
                   
                   IF( LimitActive(ind) .NEQV. LimitActive(ind2) ) THEN
                     InterfaceDof(ind) = .TRUE.
                     InterfaceDof(ind2) = .TRUE.
                   END IF
                 END DO
               END DO
             END IF
           END DO

           CALL Info(Caller,&
               'Number of interface dofs: '//TRIM(I2S(COUNT(InterfaceDof))),Level=8)
         END IF

         IF( DownStreamRemove ) THEN
           t = COUNT(InterfaceDof)
           CALL Info(Caller,'Downstream contact set dofs:'//TRIM(I2S(t)),Level=8)
         END IF
         
       
         ! Add and release dofs from the contact set:
         ! If it is removed it cannot be added. 
         !----------------------------------------------------------------------
         DO t = ElemFirst, ElemLast

           Element => Model % Elements(t)
           Model % CurrentElement => Element
           n = Element % TYPE % NumberOfNodes
           NodeIndexes => Element % NodeIndexes
           
           Found = .FALSE.
           IF( t > Model % NumberOfBulkElements ) THEN
             DO bc = 1,Model % NumberOfBCs
               IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
                 Found = .TRUE.
                 Entity => Model % BCs(bc) % Values
                 EXIT
               END IF
             END DO
             IF(.NOT. Found ) CYCLE
           ELSE             
             bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                 'Body Force', Found)
             IF(.NOT. Found ) CYCLE
             Entity => Model % BodyForces(bf) % Values
           END IF
           
           ElemLimit(1:n) = ListGetReal( Entity, &
               LimitName, n, NodeIndexes, Found)             
           IF(.NOT. Found) CYCLE
           
           ElemActive(1:n) = ListGetReal( Entity, &
               ActiveName, n, NodeIndexes, GotActive)

           DO i=1,n
             j = FieldPerm( NodeIndexes(i) )
             IF( j == 0 ) CYCLE
             ind = Dofs * ( j - 1) + Dof

             IF( LimitDone(ind) ) CYCLE
             
             ! Go through the active set and free nodes with wrong sign in contact force
             !--------------------------------------------------------------------------       
             IF( GotActive .AND. ElemActive(i) > 0.0_dp ) THEN
               IF(.NOT. LimitActive( ind ) ) THEN
                 added = added + 1
                 LimitActive(ind) = .TRUE. 
               END IF
             ELSE IF( LimitActive( ind ) ) THEN
               DoRemove = ( LimitSign * LoadValues(ind) > LimitSign * LoadEps ) 
               IF( DoRemove ) THEN
                 ! In the conservative mode only release nodes from contact set 
                 ! when they are adjacent to dofs that previously was not in the set.
                 ! This means that set is released only at the boundaries. 
                 IF( ConservativeRemove ) DoRemove = InterfaceDof( ind ) 
                 IF( DoRemove ) THEN
                   removed = removed + 1
                   LimitActive(ind) = .FALSE.
                   CYCLE
                 END IF
               END IF
             ELSE
               ! Go through the dofs that are beyond the contact surface.
               !-----------------------------------------------------------
               val = Var % Values(ind) 
               IF( Upper == 0 ) THEN
                 DoAdd = ( val < ElemLimit(i) - ValEps )
               ELSE
                 DoAdd = ( val > ElemLimit(i) + ValEps )
               END IF
               
               IF( DoAdd ) THEN
                 IF( ConservativeAdd ) DoAdd = InterfaceDof( ind ) 
                 IF( DoAdd ) THEN
                   IF( .NOT. LimitActive(ind) ) THEN
                     added = added + 1
                     LimitActive(ind) = .TRUE.
                   END IF
                 END IF
               END IF
             END IF

             ! Enforce the values to limits because nonlinear material models
             ! may otherwise lead to divergence of the iteration
             !--------------------------------------------------------------
             IF( LimitActive(ind) ) THEN
               IF( Upper == 0 ) THEN
                 Var % Values(ind) = MAX( val, ElemLimit(i) )
               ELSE
                 Var % Values(ind) = MIN( val, ElemLimit(i) )
               END IF
             END IF

             LimitDone(ind) = .TRUE.             
           END DO
         END DO
       END DO

       IF( DirectionActive ) THEN      
         ! Output some information before exiting
         !---------------------------------------------------------------------
         IF( Upper == 0 ) THEN
           CALL Info(Caller,'Determined lower soft limit set',Level=5)
         ELSE
           CALL Info(Caller,'Determined upper soft limit set',Level=5)
         END IF

         WRITE(Message,'(A,I0)') 'Number of limited dofs for '&
             //TRIM(GetVarName(Var))//': ',COUNT( LimitActive )
         CALL Info(Caller,Message,Level=5)
         
         IF(added >= 0) THEN
           WRITE(Message,'(A,I0,A)') 'Added ',added,' dofs to the set'
           CALL Info(Caller,Message,Level=5)
         END IF
         
         IF(removed >= 0) THEN
           WRITE(Message,'(A,I0,A)') 'Removed ',removed,' dofs from the set'
           CALL Info(Caller,Message,Level=5)
         END IF
       END IF
     END DO

     ! Optionally save the limiters as a field variable so that 
     ! lower limit is given value -1.0 and upper limit value +1.0.
     IF( ListGetLogical( Params,'Save Limiter',Found ) ) THEN
       
       LimitVar => VariableGet( Model % Variables, &
           GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       IF(.NOT. ASSOCIATED( LimitVar ) ) THEN
         CALL Info(Caller,'Creating field for contact: '//TRIM(GetVarName(Var)),Level=7)
         CALL VariableAddVector( Model % Variables, Solver % Mesh, Solver,&
             GetVarName(Var) //' Contact Active', Perm = FieldPerm )
         LimitVar => VariableGet( Model % Variables, &
             GetVarName(Var) // ' Contact Active',ThisOnly = .TRUE. )
       END IF

       ! Currently the visulized limit is always scalar even though the limited field could be a vector!
       DO i = 1, SIZE( LimitVar % Values ) 
         LimitVar % Values(i) = 0.0_dp
         DO j=1,Var % Dofs
           IF( ASSOCIATED( Var % LowerLimitActive ) ) THEN
             IF( Var % LowerLimitActive(Var%Dofs*(i-1)+j) ) LimitVar % Values(i) = -1.0_dp
           END IF
           IF( ASSOCIATED( Var % UpperLimitActive ) ) THEN
             IF( Var % UpperLimitActive(Var%Dofs*(i-1)+j) ) LimitVar % Values(i) = 1.0_dp
           END IF
         END DO
       END DO
     END IF

     IF( ALLOCATED( LimitDone ) ) THEN
       DEALLOCATE( LimitDone, ElemLimit, ElemInit, ElemActive ) 
     END IF
     
     IF( ALLOCATED( InterfaceDof ) ) THEN
       DEALLOCATE( InterfaceDof )
     END IF

     CALL Info(Caller,'All done',Level=12)

  END SUBROUTINE DetermineSoftLimiter
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  FUNCTION GetSolutionRotation(A,n) RESULT(rotated)
!------------------------------------------------------------------------------
    INTEGER :: n
    LOGICAL :: rotated
    REAL(KIND=dp) :: A(3,3)
!------------------------------------------------------------------------------
    INTEGER :: k,dim
!------------------------------------------------------------------------------
    dim = CoordinateSystemDimension()

    Rotated=.FALSE.

    A=0._dp
    A(1,1)=1._dp
    A(2,2)=1._dp
    A(3,3)=1._dp
    IF (NormalTangentialNOFNodes<=0) RETURN

    k = BoundaryReorder(n)
    IF (k>0) THEN
      Rotated = .TRUE.
      IF (dim==2) THEN
        A(1,1)= BoundaryNormals(k,1)
        A(1,2)=-BoundaryNormals(k,2)
        A(2,1)= BoundaryNormals(k,2)
        A(2,2)= BoundaryNormals(k,1)
      ELSE
        A(:,1)=BoundaryNormals(k,:)
        A(:,2)=BoundaryTangent1(k,:)
        A(:,3)=BoundaryTangent2(k,:)
      END IF
    END IF
!------------------------------------------------------------------------------
  END FUNCTION GetSolutionRotation
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Subroutine for determine the contact set and create the necessary data
!> for setting up the contact conditions. As input the mortar projectors,
!> the current solution, and the stiffness matrix are used.  
!------------------------------------------------------------------------------
   SUBROUTINE DetermineContact( Solver )
!------------------------------------------------------------------------------
     TYPE(Solver_t) :: Solver
!-----------------------------------------------------------------------------
     TYPE(Model_t), POINTER :: Model
     TYPE(variable_t), POINTER :: Var, LoadVar, IterVar
     TYPE(Variable_t), POINTER :: DistVar, NormalLoadVar, SlipLoadVar, VeloVar, &
         WeightVar, NormalActiveVar, StickActiveVar, GapVar
     TYPE(Element_t), POINTER :: Element
     TYPE(Mesh_t), POINTER :: Mesh
     INTEGER :: i,j,k,l,n,m,t,ind,dofs, bf, Upper, &
         ElemFirst, ElemLast, totsize, i2, j2, ind2, bc_ind, master_ind, &
         DistSign, LimitSign, DofN, DofT1, DofT2, Limited, LimitedMin, TimeStep
     REAL(KIND=dp), POINTER :: FieldValues(:), LoadValues(:), ElemLimit(:),pNormal(:,:),&
         RotatedField(:)
     REAL(KIND=dp) :: ValEps, LoadEps, val, ContactNormal(3), &
         ContactT1(3), ContactT2(3), LocalT1(3), LocalT2(3), &
         LocalNormal(3), NodalForce(3), wsum, coeff, &
         Dist, DistN, DistT1, DistT2, NTT(3,3), RotVec(3), dt
     INTEGER, POINTER :: FieldPerm(:), NodeIndexes(:)
     LOGICAL :: Found,AnyLimitBC, AnyLimitBF
     LOGICAL, ALLOCATABLE :: LimitDone(:),InterfaceDof(:)
     LOGICAL, POINTER :: LimitActive(:)
     TYPE(ValueList_t), POINTER :: Params
     CHARACTER(LEN=MAX_NAME_LEN) :: Str, LimitName, VarName, ContactType
     INTEGER :: ConservativeAfterIters, ActiveDirection, NonlinIter, CoupledIter
     LOGICAL :: ConservativeAdd, ConservativeRemove, &
         DoAdd, DoRemove, DirectionActive, Rotated, FlatProjector, PlaneProjector, &
         RotationalProjector, NormalProjector, FirstTime = .TRUE., &
         AnyRotatedContact, ThisRotatedContact, StickContact, TieContact, FrictionContact, SlipContact, &
         CalculateVelocity, NodalNormal, ResidualMode, AddDiag, SkipFriction, DoIt
     TYPE(MortarBC_t), POINTER :: MortarBC
     TYPE(Matrix_t), POINTER :: Projector, DualProjector
     TYPE(ValueList_t), POINTER :: BC, MasterBC
     REAL(KIND=dp), POINTER :: nWrk(:,:)
     LOGICAL :: CreateDual
    CHARACTER(*), PARAMETER :: Caller = 'DetermineContact'

     
     SAVE FirstTime

     CALL Info(Caller,'Setting up contact conditions',Level=8)
     
     Model => CurrentModel
     Var => Solver % Variable
     VarName = GetVarName( Var ) 
     Mesh => Solver % Mesh

     ! Is any boundary rotated or not
     AnyRotatedContact = ( NormalTangentialNOFNodes > 0 ) 

     ! The variable to be constrained by the contact algorithm
     ! Here it is assumed to be some "displacement" i.e. a vector quantity
     FieldValues => Var % Values
     FieldPerm => Var % Perm
     totsize = SIZE( FieldValues )
     dofs = Var % Dofs
     Params => Solver % Values

     IterVar => VariableGet( Model % Variables,'coupled iter')
     CoupledIter = NINT( IterVar % Values(1) )

     IterVar => VariableGet( Model % Variables,'nonlin iter')
     NonlinIter = NINT( IterVar % Values(1) )
     
     IterVar => VariableGet( Model % Variables,'timestep')
     Timestep = NINT( IterVar % Values(1) )

     IterVar => VariableGet( Mesh % Variables,'timestep size')
     IF( ASSOCIATED( IterVar ) ) THEN
       dt = IterVar % Values(1)
     ELSE
       dt = 1.0_dp
     END IF

     !FirstTime = ( NonlinIter == 1 .AND. CoupledIter == 1 ) 

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Add After Iterations',ConservativeAdd ) 
     IF( ConservativeAdd ) THEN
       IF( CoupledIter == 1 ) ConservativeAdd = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeAdd ) THEN
         CALL Info(Caller,'Adding dofs in conservative fashion',Level=8)
       END IF
     END IF

     ConservativeAfterIters = ListGetInteger(Params,&
         'Apply Limiter Conservative Remove After Iterations',ConservativeRemove ) 
     IF( ConservativeRemove ) THEN
       IF( CoupledIter == 1 ) ConservativeRemove = ( ConservativeAfterIters < NonlinIter )
       IF( ConservativeRemove ) THEN
         CALL Info(Caller,'Removing dofs in conservative fashion',Level=8)
       END IF
     END IF
         
     ResidualMode = ListGetLogical(Params,&
         'Linear System Residual Mode',Found )

     CalculateVelocity = ListGetLogical(Params,&
         'Apply Contact Velocity',Found )
     IF(.NOT. Found ) THEN
       Str = ListGetString( CurrentModel % Simulation, 'Simulation Type' )
       CalculateVelocity =  ( Str == 'transient' ) 
     END IF

     NodalNormal = ListGetLogical(Params,&
         'Use Nodal Normal',Found )

     LoadEps = ListGetConstReal(Params,'Limiter Load Tolerance',Found ) 
     IF(.NOT. Found ) LoadEps = EPSILON( LoadEps )
         
     ValEps = ListGetConstReal(Params,'Limiter Value Tolerance',Found ) 
     IF(.NOT. Found ) ValEps = EPSILON( ValEps )

     IF( .NOT. ASSOCIATED( Model % Solver % MortarBCs ) ) THEN
       CALL Fatal(Caller,'Cannot apply contact without projectors!')
     END IF

     ! a) Create rotateted contact if needed
     CALL RotatedDisplacementField() 

     ! b) Create and/or obtain pointers to boundary variables 
     CALL GetContactFields( FirstTime )

     ! c) Calculate the contact loads to the normal direction
     LoadVar => CalculateContactLoad() 
     LoadValues => LoadVar % Values

     
     ! Loop over each contact pair
     !--------------------------------------------------------------
     DO bc_ind = 1, Model % NumberOfBCs
       
       MortarBC => Model % Solver % MortarBCs(bc_ind)  
       IF( .NOT. ASSOCIATED( MortarBC ) ) CYCLE

       Projector => MortarBC % Projector
       IF(.NOT. ASSOCIATED(Projector) ) CYCLE

       BC => Model % BCs(bc_ind) % Values

       CALL Info(Caller,'Set contact for boundary: '&
           //TRIM(I2S(bc_ind)),Level=8)
       Model % Solver % MortarBCsChanged = .TRUE.
       
       FlatProjector = ListGetLogical( BC, 'Flat Projector',Found ) 
       PlaneProjector = ListGetLogical( BC, 'Plane Projector',Found )
       RotationalProjector = ListGetLogical( BC, 'Rotational Projector',Found ) .OR. &
           ListGetLogical( BC, 'Cylindrical Projector',Found )
       NormalProjector = ListGetLogical( BC, 'Normal Projector',Found )
       
       ! Is the current boundary rotated or not
       ThisRotatedContact = ListGetLogical( BC,'Normal-Tangential '//TRIM(VarName),Found)

       IF( FlatProjector ) THEN
         ActiveDirection = ListGetInteger( BC, 'Flat Projector Coordinate',Found )
         IF( .NOT. Found ) ActiveDirection = dofs       
       ELSE IF( PlaneProjector ) THEN
         pNormal => ListGetConstRealArray( BC,'Plane Projector Normal',Found)
         IF( ThisRotatedContact ) THEN
           ActiveDirection = 1
         ELSE
           ActiveDirection = 1
           DO i=2,3
             IF( ABS( pnormal(i,1) ) > ABS( pnormal(ActiveDirection,1) ) ) THEN
               ActiveDirection = i
             END IF
           END DO
           CALL Info(Caller,'Active direction set to: '//TRIM(I2S(ActiveDirection)))
         END IF
       ELSE IF( RotationalProjector .OR. NormalProjector ) THEN
         ActiveDirection = 1
         IF( .NOT. ThisRotatedContact ) THEN
           CALL Warn(Caller,'Rotational and normal projectors should only work with N-T coordinates!')
         END IF
       ELSE
         CALL Fatal(Caller,'Projector must be current either flat, plane, cylinder or rotational!')
       END IF
      

       ! Get the pointer to the other side i.e. master boundary  
       master_ind = ListGetInteger( BC,'Mortar BC',Found )
       IF( .NOT. Found ) master_ind = ListGetInteger( BC,'Contact BC',Found )
       MasterBC => Model % BCs(master_ind) % Values
       
       ! If we have dual projector we may use it to map certain quantities directly to master nodes
       DualProjector => Projector % Ematrix
       CreateDual = ASSOCIATED( DualProjector )
       IF( CreateDual ) THEN
         CALL Info(Caller,'Using also the dual projector',Level=8)
       END IF
      
       ! If we have N-T system then the mortar condition for the master side
       ! should have reverse sign as both normal displacement diminish the gap.
       IF( ThisRotatedContact ) THEN
         IF( master_ind > 0 ) THEN
           IF( .NOT. ListGetLogical( MasterBC, &
               'Normal-Tangential '//TRIM(VarName),Found) ) THEN
             CALL Fatal(Caller,'Master boundary '//TRIM(I2S(master_ind))//&
                 ' should also have N-T coordinates!')
           END IF
         END IF

         CALL Info(Caller,'We have a normal-tangential system',Level=6)
         MortarBC % MasterScale = -1.0_dp
         DofN = 1
       ELSE                 
         DofN = ActiveDirection 
       END IF

       ! Get the degrees of freedom related to the normal and tangential directions
       DofT1 = 0; DofT2 = 0
       DO i=1,dofs
         IF( i == DofN ) CYCLE
         IF( DofT1 == 0 ) THEN
           DofT1 = i 
           CYCLE
         END IF
         IF( DofT2 == 0 ) THEN
           DofT2 = i
           CYCLE
         END IF
       END DO

       ! This is the normal that is used to detect the signed distance
       ! and tangent vectors used to detect surface velocity
       IF( PlaneProjector ) THEN
         ContactNormal = pNormal(1:3,1)
       ELSE
         ContactNormal = 0.0_dp
         ContactNormal(ActiveDirection) = 1.0_dp
       END IF
       ContactT1 = 0.0_dp
       ContactT1(DofT1) = 1.0_dp
       ContactT2 = 0.0_dp
       IF(DofT2>0) ContactT2(DofT2) = 1.0_dp

       ! Get the contact type. There are four possibilities currently. 
       ! Only one is active at a time while others are false. 
       StickContact = .FALSE.; TieContact = .FALSE.
       FrictionContact = .FALSE.; SlipContact = .FALSE.

       ContactType = ListGetString( BC,'Contact Type',Found ) 
       IF( Found ) THEN
         SELECT CASE ( ContactType )
         CASE('stick')
           StickContact = .TRUE.
         CASE('tie')
           TieContact = .TRUE.
         CASE('friction')
           FrictionContact = .TRUE.
         CASE('slide')
           SlipContact = .TRUE.
         CASE Default
           CALL Fatal(Caller,'Unknown contact type: '//TRIM(ContactType))
         END SELECT
       ELSE
         StickContact = ListGetLogical( BC,'Stick Contact',Found )
         IF(.NOT. Found ) TieContact = ListGetLogical( BC,'Tie Contact',Found )
         IF(.NOT. Found ) FrictionContact = ListGetLogical( BC,'Friction Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slip Contact',Found )
         IF(.NOT. Found ) SlipContact = ListGetLogical( BC,'Slide Contact',Found )
         IF(.NOT. Found ) THEN 
           CALL Warn(Caller,'No contact type given, assuming > Slip Contact <')
           SlipContact = .TRUE.
         END IF
       END IF

       IF( StickContact ) CALL Info(Caller,'Using stick contact for displacement',Level=10)
       IF( TieContact ) CALL Info(Caller,'Using tie contact for displacement',Level=10)
       IF( FrictionContact ) CALL Info(Caller,'Using friction contact for displacement',Level=10)
       IF( SlipContact ) CALL Info(Caller,'Using slip contact for displacement',Level=10)
       

       ! At the start it may be beneficial to assume initial tie contact
       IF( (FrictionContact .OR. StickContact .OR. SlipContact ) .AND. &
           (TimeStep == 1 .AND. NonlinIter == 1 ) ) THEN
         DoIt = ListGetLogical(BC,'Initial Tie Contact',Found )
         IF( DoIt ) THEN
           FrictionContact = .FALSE.; StickContact = .FALSE.; SlipContact = .FALSE.
           TieContact = .TRUE.
           CALL Info(Caller,'Assuming initial tie contact',Level=10)
         END IF
       END IF
         
       ! At the first time it may be beneficial to assume frictionless initial contact.
       SkipFriction = .FALSE.
       IF( (FrictionContact .OR. StickContact .OR. SlipContact ) .AND. TimeStep == 1 ) THEN
         DoIt = .NOT. ListGetLogical(BC,'Initial Contact Friction',Found )
         IF( DoIt ) THEN
           FrictionContact = .FALSE.; StickContact = .FALSE.
           SlipContact = .TRUE.
           SkipFriction = .TRUE.
           CALL Info(Caller,'Assuming frictionless initial contact',Level=10)
         END IF
       ELSE IF( ( FrictionContact .OR. SlipContact) .AND. NonlinIter == 1 ) THEN
         DoIt = ListGetLogical(BC,'Nonlinear System Initial Stick',Found )
         IF(.NOT. Found ) THEN
           ! If contact velocity is not given then it is difficult to determine the direction at 
           ! start of nonlinear iteration when the initial guess still reflects the old displacements. 
           DoIt = .NOT. ListCheckPresent( BC,'Contact Velocity') .AND. &
               ListCheckPresent( BC,'Dynamic Friction Coefficient')
         END IF
         IF( DoIt ) THEN
           FrictionContact = .FALSE.
           SlipContact = .FALSE.
           StickContact = .TRUE.
           CALL Info(Caller,'Assuming sticking in first iteration initial contact',Level=10)
         END IF        
       END IF

       ! If we have stick contact then create a diagonal entry to the projection matrix.
       IF( StickContact .OR. FrictionContact ) THEN
         AddDiag = ListCheckPresent( BC,'Stick Contact Coefficient')      
       ELSE
         AddDiag = .FALSE.
       END IF

       ! d) allocate and initialize all necessary vectors for the contact 
       !------------------------------------------------------------------
       CALL InitializeMortarVectors()
     
       ! e) If the contact set is set up in a conservative fashion we need to mark interface nodes
       !------------------------------------------------------------------
       IF( ConservativeAdd .OR. ConservativeRemove ) THEN
         CALL MarkInterfaceDofs()
       END IF

       ! f) Compute the normal load used to determine whether contact should be released.
       !    Also check the direction to which the signed distance should be computed
       !------------------------------------------------------------------
       CALL CalculateContactPressure()

       ! g) Calculate the distance used to determine whether contact should be added
       !------------------------------------------------------------------
       CALL CalculateMortarDistance()
        
       ! h) Determine the contact set in normal direction
       !------------------------------------------------------------------
       CALL NormalContactSet()

       ! i) If requested ensure a minimum number of contact nodes
       !-------------------------------------------------------------------
       LimitedMin = ListGetInteger( BC,'Contact Active Set Minimum',Found)
       IF( Found ) CALL IncreaseContactSet( LimitedMin )

       ! j) Determine the stick set in tangent direction
       !------------------------------------------------------------------
       CALL TangentContactSet()
       
       ! k) Add the stick coefficient if present
       !------------------------------------------------------------------
       IF( AddDiag ) THEN
         CALL StickCoefficientSet()
       END IF

       ! l) We can map information from slave to master either by creating a dual projector
       !    or using the transpose of the original projector to map field from slave to master.
       !-----------------------------------------------------------------------------------
       IF(.NOT. CreateDual ) THEN
         CALL ProjectFromSlaveToMaster()
       END IF

       ! m) If we have dynamic friction then add it 
       IF( .NOT. SkipFriction .AND. ( SlipContact .OR. FrictionContact ) ) THEN
         CALL SetSlideFriction()
       END IF

       IF( ConservativeAdd .OR. ConservativeRemove ) THEN
         DEALLOCATE( InterfaceDof )
       END IF
     END DO
     
     ! Use N-T coordinate system for the initial guess
     ! This is mandatory if using the residual mode linear solvers 
     IF( AnyRotatedContact ) THEN
       DEALLOCATE( RotatedField ) 
     END IF
     

     FirstTime = .FALSE.
     CALL Info(Caller,'All done',Level=10)

   CONTAINS


     ! Given the cartesian solution compute the rotated solution.
     !-------------------------------------------------------------------------
     SUBROUTINE RotatedDisplacementField( ) 

       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,m

       IF( .NOT. AnyRotatedContact ) RETURN

       CALL Info(Caller,'Rotating displacement field',Level=8)
       ALLOCATE( RotatedField(Solver % Matrix % NumberOfRows ) )
       RotatedField = Var % Values
       
       DO i=1,Solver % Mesh % NumberOfNodes
         j = Solver % Variable % Perm(i)
         IF( j == 0 ) CYCLE
         m = BoundaryReorder(i)
         IF( m == 0 ) CYCLE
         
         RotVec = 0._dp
         DO k=1,Var % DOFs
           RotVec(k) = RotatedField(Var % DOfs*(j-1)+k)
         END DO
         CALL RotateNTSystem( RotVec, i )
         DO k=1,Var % DOFs
           RotatedField(Var % Dofs*(j-1)+k) = RotVec( k )
         END DO
       END DO

     END SUBROUTINE RotatedDisplacementField


     ! Given the previous solution and the current stiffness matrix 
     ! computes the load normal to the surface i.e. the contact load.
     ! If we have normal-tangential coordinate system then also the load is in 
     ! the same coordinate system. 
     !-------------------------------------------------------------------------
     FUNCTION CalculateContactLoad( ) RESULT ( LoadVar )

       TYPE(Variable_t), POINTER :: LoadVar
       REAL(KIND=dp), POINTER :: TempX(:)
       REAL(KIND=dp) :: RotVec(3)
       INTEGER :: i,j,k,m


       CALL Info(Caller,'Determining contact load for contact problems',Level=10)

       LoadVar => VariableGet( Model % Variables, &
           TRIM(VarName) // ' Contact Load',ThisOnly = .TRUE. )
       IF( .NOT. ASSOCIATED( LoadVar ) ) THEN
         CALL Fatal(Caller, &
             'No Loads associated with variable: '//GetVarName(Var) )
       END IF

       IF( AnyRotatedContact ) THEN
         TempX => RotatedField 
       ELSE
         TempX => FieldValues
       END IF

       CALL CalculateLoads( Solver, Solver % Matrix, TempX, Var % DOFs, .FALSE., LoadVar ) 

     END FUNCTION CalculateContactLoad


     ! Create fields where the contact information will be saved.
     ! Create the fields both for slave and master nodes at each 
     ! contact pair.
     !--------------------------------------------------------------
     SUBROUTINE GetContactFields( DoAllocate )

       LOGICAL :: DoAllocate
       INTEGER, POINTER :: BoundaryPerm(:), Indexes(:)
       INTEGER :: i,j,k,t,n
       TYPE(Element_t), POINTER :: Element
       LOGICAL, ALLOCATABLE :: ActiveBCs(:)


       IF( DoAllocate ) THEN
         CALL Info(Caller,'Creating contact fields',Level=8)
        
         ALLOCATE( BoundaryPerm(Mesh % NumberOfNodes) )
         BoundaryPerm = 0

         ALLOCATE( ActiveBCs(Model % NumberOfBcs ) )
         ActiveBCs = .FALSE.

         DO i=1,Model % NumberOfBCs 
           j = ListGetInteger( Model % BCs(i) % Values,'Mortar BC',Found ) 
           IF(.NOT. Found ) THEN
             j = ListGetInteger( Model % BCs(i) % Values,'Contact BC',Found ) 
           END IF
           IF( j > 0 ) THEN
             ActiveBCs(i) = .TRUE.
             ActiveBCs(j) = .TRUE. 
           END IF
         END DO

         DO t=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

           Element => Mesh % Elements( t )                 
           DO i = 1, Model % NumberOfBCs
             IF ( Element % BoundaryInfo % Constraint == Model % BCs(i) % Tag ) THEN
               IF( ActiveBCs(i) ) THEN
                 BoundaryPerm( Element % NodeIndexes ) = 1
               END IF
             END IF
           END DO
         END DO

         DEALLOCATE( ActiveBCs )

         j = 0
         DO i=1,Mesh % NumberOfNodes
           IF( BoundaryPerm(i) > 0 ) THEN
             j = j + 1
             BoundaryPerm(i) = j
           END IF
         END DO

         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Distance',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Gap',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Normalload',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Slipload',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Weight',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Active',1,Perm = BoundaryPerm )
         CALL VariableAddVector( Model % Variables,Mesh,Solver,&
             TRIM(VarName)//' Contact Stick',1,Perm = BoundaryPerm )
         IF( CalculateVelocity ) THEN
           CALL VariableAddVector( Model % Variables,Mesh,Solver,&
               TRIM(VarName)//' Contact Velocity',Dofs,Perm = BoundaryPerm )
         END IF
       END IF

       DistVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Distance')
       GapVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Gap')
       NormalLoadVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Normalload')
       SlipLoadVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Slipload')
       WeightVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Weight')
       NormalActiveVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Active')
       StickActiveVar => VariableGet( Model % Variables,&
           TRIM(VarName)//' Contact Stick') 
       IF( CalculateVelocity ) THEN
         VeloVar => VariableGet( Model % Variables,&
             TRIM(VarName)//' Contact Velocity')
       END IF

       NormalActiveVar % Values = -1.0_dp
       StickActiveVar % Values = -1.0_dp

     END SUBROUTINE GetContactFields



     ! Allocates the vectors related to the mortar contact surface, if needed.
     ! Initialize the mortar vectors and mortar permutation future use. 
     ! As the geometry changes the size of the projectors may also change. 
     !----------------------------------------------------------------------------
     SUBROUTINE InitializeMortarVectors()

       INTEGER :: onesize, totsize
       INTEGER, POINTER :: Perm(:)
       LOGICAL, POINTER :: Active(:)
       REAL(KIND=dp), POINTER :: Diag(:)
       LOGICAL :: SamePerm, SameSize

       
       onesize = Projector % NumberOfRows
       totsize = Dofs * onesize

       IF( .NOT. AddDiag .AND. ASSOCIATED(MortarBC % Diag) ) THEN
         DEALLOCATE( MortarBC % Diag ) 
       END IF


       ! Create the permutation that is later need in putting the diag and rhs to correct position
       ALLOCATE( Perm( SIZE( FieldPerm ) ) )
       Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j == 0 ) CYCLE
         Perm( j ) = i
       END DO

       ! First time nothing is allocated
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info(Caller,'Allocating projector mortar vectors',Level=10)
         ALLOCATE( MortarBC % Active( totsize ), MortarBC % Rhs( totsize) )
         MortarBC % Active = .FALSE.
         MortarBC % Rhs = 0.0_dp
         MortarBC % Perm => Perm 

         IF( AddDiag ) THEN
           ALLOCATE( MortarBC % Diag( totsize ) )
           MortarBC % Diag = 0.0_dp
         END IF

         RETURN
       END IF

       
       ! If permutation has changed we need to change the vectors also
       SamePerm = ANY( Perm /= MortarBC % Perm )
       SameSize = ( SIZE(MortarBC % Rhs) == totsize )

       ! Permutation unchanged, just return
       IF( SamePerm ) THEN
         DEALLOCATE( Perm ) 
         RETURN
       END IF
       
       ! Permutation changes, and also sizes changed
       IF(.NOT. SameSize ) THEN
         DEALLOCATE( MortarBC % Rhs )
         ALLOCATE( MortarBC % Rhs( totsize ) )
         MortarBC % Rhs = 0.0_dp
       END IF

       ! .NOT. SamePerm
       ALLOCATE(Active(totsize))
       Active = .FALSE.

       IF( AddDiag ) THEN
         ALLOCATE( Diag(totsize) )
         Diag = 0.0_dp
       END IF


       DO i=1,SIZE( Perm ) 
         j = Perm(i)
         IF( j == 0 ) CYCLE

         k = MortarBC % Perm(i)
         IF( k == 0 ) CYCLE

         DO l=1,Dofs
           Active(Dofs*(j-1)+l) = MortarBC % Active(Dofs*(k-1)+l)
         END DO
       END DO

       DEALLOCATE( MortarBC % Active ) 
       DEALLOCATE( MortarBC % Perm ) 
       MortarBC % Active => Active 
       MortarBC % Perm => Perm 

       IF( AddDiag ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) THEN
           DEALLOCATE( MortarBC % Diag ) 
         END IF
         MortarBC % Diag => Diag 
       END IF

       CALL Info(Caller,'Copied > Active < flag to changed projector',Level=8)

     END SUBROUTINE InitializeMortarVectors

     

     ! Make a list of interface dofs to allow conservative algorithms. 
     ! There only nodes that are at the interface are added or removed from the set.
     !------------------------------------------------------------------------------
     SUBROUTINE MarkInterfaceDofs()
       
       INTEGER :: i,j,i2,j2,k,k2,l,n,ind,ind2,elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(Element_t), POINTER :: Element
       
       CALL Info(Caller,'Marking interface dofs for conservative adding/removal',Level=8)

       IF(.NOT. ALLOCATED( InterfaceDof ) ) THEN
         ALLOCATE( InterfaceDof( SIZE(MortarBC % Active) ) )
       END IF
       InterfaceDof = .FALSE. 


       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         
         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE
         
         n = Element % TYPE % NumberOfNodes
         Indexes => Element % NodeIndexes
         
         DO i=1,n
           j = FieldPerm( Indexes(i) )
           IF( j == 0 ) CYCLE
           k = MortarBC % Perm( Indexes(i) )
           
           DO i2 = i+1,n
             j2 = FieldPerm( Indexes(i2) )
             IF( j2 == 0 ) CYCLE
             k2 = MortarBC % perm( Indexes(i2) )
             
             DO l=1,Dofs             
               ind = Dofs * ( k - 1 ) + l
               ind2 = Dofs * ( k2 - 1) + l
               
               IF( MortarBC % Active(ind) .NEQV. MortarBC % Active(ind2) ) THEN
                 InterfaceDof(ind) = .TRUE.
                 InterfaceDof(ind2) = .TRUE.
               END IF
             END DO
           END DO
         END DO
       END DO

       n = COUNT(InterfaceDof)
       CALL Info(Caller,&
           'Number of interface dofs: '//TRIM(I2S(n)),Level=8)
     END SUBROUTINE MarkInterfaceDofs
     

     ! Calculates the signed distance that is used to define whether we have contact or not.
     ! If distance is negative then we can later add the corresponding node to the contact set
     ! Also computes the right-hand-side of the mortar equality constrained which is the 
     ! desired distance in the active direction. Works also for residual mode which greatly 
     ! improves the convergence for large displacements.  
     !----------------------------------------------------------------------------------------
     SUBROUTINE CalculateMortarDistance()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactVec(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3), CartVec(3), ContactDist
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, wsum, wsumM, mult
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL :: IsSlave, IsMaster, DistanceSet
       LOGICAL, ALLOCATABLE :: SlaveNode(:), MasterNode(:), NodeDone(:)
       INTEGER, POINTER :: Indexes(:)
       INTEGER :: elemcode, CoeffSign
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:)
       INTEGER :: l2,elem,i1,i2,j1,j2
       LOGICAL :: LinearContactGap, DebugNormals
       
       
       CALL Info(Caller,'Computing distance between mortar boundaries',Level=14)

       DispVals => Solver % Variable % Values
       IF( .NOT. ASSOCIATED( DispVals ) ) THEN
         CALL Fatal('CalculateMortarDistance','Displacement variable not associated!')
       END IF

       IF( CalculateVelocity ) THEN
         IF( .NOT. ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
           CALL Fatal('CalculateMortarDistance','Displacement PrevValues not associated!')         
         END IF
         IF( Solver % TimeOrder == 1 ) THEN
           PrevDispVals => Solver % Variable % PrevValues(:,1)
         ELSE
           PrevDispVals => Solver % Variable % PrevValues(:,3)
         END IF
         IF(.NOT. ASSOCIATED( PrevDispVals ) ) CALL Fatal(Caller,&
             'Previous displacement field required!')
       END IF

       LinearContactGap = ListGetLogical( Model % Simulation,&
           'Contact BCs linear gap', Found )      

       ALLOCATE( SlaveNode( Mesh % NumberOfNodes ) ) 
       SlaveNode = .FALSE.

       IF( CreateDual ) THEN
         ALLOCATE( MasterNode( Mesh % NumberOfNodes ) ) 
         MasterNode = .FALSE.
       END IF

       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) THEN
           SlaveNode( Element % NodeIndexes ) = .TRUE.
         END IF
         IF( CreateDual ) THEN
           IF ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) THEN
             MasterNode( Element % NodeIndexes ) = .TRUE.
           END IF
         END IF
       END DO

       ! First create the master, then the slave if needed
       IsSlave = .TRUE.
       IsMaster = .NOT. IsSlave
       ActiveProjector => Projector
       MinDist = HUGE(MinDist)
       MaxDist = -HUGE(MaxDist)

       IF( .NOT. ASSOCIATED( ActiveProjector ) ) THEN
         CALL Fatal('CalculateMortarDistance','Projector not associated!')
       END IF

       DebugNormals = ListGetLogical( Params,'Debug Normals',Found ) 

       IF( DebugNormals ) THEN
         PRINT *,'Flags:',TieContact,ResidualMode,ThisRotatedContact,NodalNormal,StickContact,RotationalProjector
       END IF

       
100    CONTINUE

       DO i = 1,ActiveProjector % NumberOfRows

         j = ActiveProjector % InvPerm(i)

         IF( j == 0 ) CYCLE

         wsum = 0.0_dp
         wsumM = 0.0_dp
         Dist = 0.0_dp
         DistN = 0.0_dp
         DistT1 = 0.0_dp
         DistT2 = 0.0_dp
         ContactVelo = 0.0_dp
         ContactVec = 0.0_dp
         DistanceSet = .FALSE.
         ContactDist = 0.0_dp
         CartVec = 0.0_dp
         
         ! This is the most simple contact condition. We just want no slip on the contact.
         IF( TieContact .AND. .NOT. ResidualMode ) GOTO 200

         ! Get the normal of the slave surface.
         IF( ThisRotatedContact ) THEN
           Rotated = GetSolutionRotation(NTT, j )
           LocalNormal = NTT(:,1)
           LocalNormal0 = LocalNormal
           LocalT1 = NTT(:,2)
           IF( Dofs == 3 ) LocalT2 = NTT(:,3)
         ELSE
           LocalNormal = ContactNormal
           LocalT1 = ContactT1
           IF( Dofs == 3 ) LocalT2 = ContactT2 
         END IF

         ! Compute normal of the master surface from the average sum of normals
         IF( NodalNormal ) THEN
           LocalNormal = 0.0_dp
           LocalT1 = 0.0_dp
           LocalT2 = 0.0_dp

           DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
             k = ActiveProjector % Cols(j)

             l = FieldPerm( k ) 
             IF( l == 0 ) CYCLE

             coeff = ActiveProjector % Values(j)             
             Rotated = GetSolutionRotation(NTT, k )

             ! Weighted direction for the unit vectors
             LocalNormal = LocalNormal + coeff * NTT(:,1)
             LocalT1 = LocalT1 + coeff * NTT(:,2)
             IF( Dofs == 3 ) LocalT2 = LocalT2 + coeff * NTT(:,3)
           END DO

           ! Normalize the unit vector length to one
           LocalNormal = LocalNormal / SQRT( SUM( LocalNormal**2 ) )
           LocalT1 = LocalT1 / SQRT( SUM( LocalT1**2 ) )
           IF( Dofs == 3 ) LocalT2 = LocalT2 / SQRT( SUM( LocalT1**2 ) )

           !PRINT *,'NodalNormal:',i,j,LocalNormal0,LocalNormal
         END IF

         ! For debugging reason, check that normals are roughly opposite
         IF( DebugNormals ) THEN
           DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
             k = ActiveProjector % Cols(j)
             
             l = FieldPerm( k ) 
             IF( l == 0 ) CYCLE

             Rotated = GetSolutionRotation(NTT, k )
             coeff = SUM( LocalNormal * NTT(:,1) ) + SUM( LocalT1*NTT(:,2)) + SUM(LocalT2*NTT(:,3))
             IF( SlaveNode(k) .AND. coeff < 2.5_dp ) THEN
               Found = .TRUE.
               !PRINT *,'Slave Normal:',i,j,k,Rotated,coeff
             ELSE IF( .NOT. SlaveNode(k) .AND. coeff > -2.5_dp ) THEN
               Found = .TRUE.
               !PRINT *,'Master Normal:',i,j,k,Rotated,coeff
             ELSE
               Found = .FALSE.
             END IF
             IF( Found ) THEN
               !PRINT *,'Prod:',SUM( LocalNormal * NTT(:,1) ), SUM( LocalT1*NTT(:,2)), SUM(LocalT2*NTT(:,3))
               !PRINT *,'N:',LocalNormal,NTT(:,1)
               !PRINT *,'T1:',LocalT1,NTT(:,2)
               !PRINT *,'T2:',LocalT2,NTT(:,3)
             END IF
           END DO
         END IF

         ! Compute the weigted distance in the normal direction.
         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           coeff = ActiveProjector % Values(j)
                  
           ! Only compute the sum related to the active projector
           IF( SlaveNode(k) ) THEN
             wsum = wsum + coeff
           ELSE
             wsumM = wsumM + coeff
           END IF           
         END DO

         IF( ABS( wsum ) <= TINY( wsum ) ) THEN
           CALL Fatal('CalculateMortarDistance','wsum seems to be almost zero!')
         END IF
         IF( ABS( wsumM ) <= TINY( wsumM ) ) THEN
           CALL Fatal('CalculateMortarDistance','wsumM seems to be almost zero!')
         END IF

         ! Slave and master multipliers should sum up to same value
         mult = ABS( wsum / wsumM ) 

         ! Compute the weigted distance in the normal direction.
         DO j = ActiveProjector % Rows(i),ActiveProjector % Rows(i+1)-1
           k = ActiveProjector % Cols(j)

           l = FieldPerm( k ) 
           IF( l == 0 ) CYCLE

           ! This includes only the coordinate since the displacement
           ! is added to the coordinate!
           coeff = ActiveProjector % Values(j)

           CoeffSign = 1
           
           ! Only compute the sum related to the active projector
           IF( .NOT. SlaveNode(k) ) THEN
             coeff = mult * coeff
             IF( ThisRotatedContact ) CoeffSign = -1
           END IF
             
           IF( dofs == 2 ) THEN
             disp(1) = DispVals( 2 * l - 1)
             disp(2) = DispVals( 2 * l )
             disp(3) = 0.0_dp
           ELSE
             disp(1) = DispVals( 3 * l - 2)
             disp(2) = DispVals( 3 * l - 1 )
             disp(3) = DispVals( 3 * l )
           END IF

           ! If nonlinear analysis is used we may need to cancel the introduced gap due to numerical errors 
           IF( TieContact .AND. ResidualMode ) THEN
             IF( ThisRotatedContact ) THEN
               ContactVec(1) = ContactVec(1) + coeff * SUM( LocalNormal * Disp )
               ContactVec(2) = ContactVec(2) + coeff * SUM( LocalT1 * Disp )
               IF( Dofs == 3) ContactVec(3) = ContactVec(3) + coeff * SUM( LocalT2 * Disp )
             ELSE
               ContactVec(1) = ContactVec(1) + coeff * SUM( ContactNormal * Disp )
               ContactVec(2) = ContactVec(2) + coeff * SUM( ContactT1 * Disp )
               IF( Dofs == 3 ) ContactVec(3) = ContactVec(3) + coeff * SUM( ContactT2 * Disp ) 
             END IF
             CYCLE
           END IF

           coord(1) = Mesh % Nodes % x( k ) 
           coord(2) = Mesh % Nodes % y( k ) 
           coord(3) = Mesh % Nodes % z( k ) 

           IF( CalculateVelocity ) THEN
             IF( dofs == 2 ) THEN
               PrevDisp(1) = PrevDispVals( 2 * l - 1)
               PrevDisp(2) = PrevDispVals( 2 * l )
               PrevDisp(3) = 0.0_dp
             ELSE
               PrevDisp(1) = PrevDispVals( 3 * l - 2)
               PrevDisp(2) = PrevDispVals( 3 * l - 1 )
               PrevDisp(3) = PrevDispVals( 3 * l )
             END IF
           END IF

           ! If the linear system is in residual mode also set the current coordinate in residual mode too!
           ! Note that displacement field is given always in cartesian coordinates!
           IF( ResidualMode ) THEN
             Coord = Coord + Disp
           END IF

           ! DistN is used to give the distance that we need to move the original coordinates
           ! in the wanted direction in order to have contact.
           IF( ThisRotatedContact ) THEN
             ContactVec(1) = ContactVec(1) + coeff * SUM( LocalNormal * Coord )
           ELSE
             ContactVec(1) = ContactVec(1) + coeff * SUM( ContactNormal * Coord )
           END IF

           ! Tangential distances needed to move the original coordinates to the contact position
           ! If stick is required then we want to keep the tangential slip zero. 
           IF( StickContact ) THEN             
             SlipCoord = -PrevDisp 
             IF( ResidualMode ) SlipCoord = SlipCoord + Disp 

             IF( ThisRotatedContact ) THEN
               ContactVec(2) = ContactVec(2) + coeff * SUM( LocalT1 * SlipCoord )
               IF( Dofs == 3) ContactVec(3) = ContactVec(3) + coeff * SUM( LocalT2 * SlipCoord )
             ELSE
               ContactVec(2) = ContactVec(2) + coeff * SUM( ContactT1 * SlipCoord )
               IF( Dofs == 3 ) ContactVec(3) = ContactVec(3) + coeff * SUM( ContactT2 * SlipCoord )
             END IF
           END IF

           ! If not in the residual mode still take into account the displacement for the condition
           IF( .NOT. ResidualMode ) Coord = Coord + Disp

           ! Dist is used to compute the current signed distance that is used to determine
           ! whether we have contact or not. 
           IF( RotationalProjector ) THEN
             Dist = Dist + coeff * SQRT( SUM( Coord**2 ) )
           ELSE IF( NormalProjector ) THEN
             Dist = Dist + coeff * SUM( LocalNormal * Coord )
           ELSE             
             Dist = Dist + coeff * SUM( ContactNormal * Coord )
           END IF

           CartVec = CartVec + coeff * Coord
           
           IF( CalculateVelocity ) THEN
             Velo = ( Disp - PrevDisp ) !/ dt
             ContactVelo(1) = ContactVelo(1) + coeff * SUM( Velo * LocalNormal ) 
             ContactVelo(2) = ContactVelo(2) + coeff * SUM( Velo * LocalT1 )
             ContactVelo(3) = ContactVelo(3) + coeff * SUM( Velo * LocalT2 ) 
           END IF
           DistanceSet = .TRUE.
         END DO

         ! Divide by weight to get back to real distance in the direction of the normal
         ContactVec = ContactVec / wsum 
         Dist = DistSign * Dist / wsum
         IF( CalculateVelocity ) THEN
           ContactVelo = ContactVelo / wsum
         END IF
         CartVec = CartVec / wsum
         
200      IF( IsSlave ) THEN
           MortarBC % Rhs(Dofs*(i-1)+DofN) = -ContactVec(1)
           IF( StickContact .OR. TieContact ) THEN
             MortarBC % Rhs(Dofs*(i-1)+DofT1) = -ContactVec(2) 
             IF( Dofs == 3 ) THEN
               MortarBC % Rhs(Dofs*(i-1)+DofT2) = -ContactVec(3)
             END IF
           END IF

           MinDist = MIN( Dist, MinDist ) 
           MaxDist = MAX( Dist, MaxDist )
         END IF

         IF( IsMaster ) THEN
           Dist = -Dist
           ContactVelo = -ContactVelo
         END IF

         ! We use the same permutation for all boundary variables
         IF(ActiveProjector % InvPerm(i) <= 0 ) CYCLE
         j = DistVar % Perm( ActiveProjector % InvPerm(i) )

         DistVar % Values( j ) = Dist

         GapVar % Values( j ) = ContactVec(1)

         IF( CalculateVelocity ) THEN
           DO k=1,Dofs             
             VeloVar % Values( Dofs*(j-1)+k ) = ContactVelo(k) 
           END DO
         END IF
       END DO

       IF( IsSlave ) THEN
         IF( CreateDual ) THEN
           IsSlave = .FALSE.
           IsMaster = .NOT. IsSlave
           ActiveProjector => DualProjector
           GOTO 100
         END IF
       END IF

       
       IF( LinearContactGap ) THEN       
         DO elem=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
           
           Element => Mesh % Elements( elem )         
           
           IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
           IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 
           IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
           
           Indexes => Element % NodeIndexes         
           ElemCode = Element % TYPE % ElementCode           

           SELECT CASE ( ElemCode )

           CASE( 408 )
             DO i=5,8
               i1=i-4
               i2=i1+1
               IF(i2==5) i2=1
               j = DistVar % Perm(Indexes(i))
               IF( j == 0 ) CYCLE
               j1 = DistVar % Perm(Indexes(i1))
               j2 = DistVar % Perm(Indexes(i2))

               DistVar % Values(j) = 0.5_dp * &
                   ( DistVar % Values(j1) + DistVar % Values(j2))
               GapVar % Values(j) = 0.5_dp * &
                   ( GapVar % Values(j1) + GapVar % Values(j2))
               
               IF( CalculateVelocity ) THEN
                 DO k=1,Dofs
                   VeloVar % Values( Dofs*(j-1)+k ) = 0.5_dp * &
                       ( VeloVar % Values(Dofs*(j1-1)+k) + VeloVar % Values(Dofs*(j2-1)+k))  
                 END DO
               END IF
             END DO
             
           CASE DEFAULT
             CALL Fatal(Caller,'Implement linear gaps for: '//TRIM(I2S(ElemCode)))
           END SELECT
         END DO
       END IF
     
       DEALLOCATE( SlaveNode )
       IF( CreateDual ) DEALLOCATE( MasterNode )

       IF( InfoActive(20 ) ) THEN
         PRINT *,'Distance Range:',MinDist, MaxDist
         PRINT *,'Distance Offset:',MINVAL( MortarBC % Rhs ), MAXVAL( MortarBC % Rhs )
       END IF
         
     END SUBROUTINE CalculateMortarDistance



     ! Calculates the contact pressure in the normal direction from the nodal loads.
     ! The nodal loads may be given either in cartesian or n-t coordinate system. 
     !-------------------------------------------------------------------------------
     SUBROUTINE CalculateContactPressure()
       
       INTEGER :: elem
       INTEGER, POINTER :: Indexes(:)
       TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
       TYPE(Nodes_t) :: Nodes
       INTEGER :: i,j,k,t,CoordSys, NormalSign0, NormalSign, NormalCount
       REAL(KIND=dp) :: s, x, DetJ, u, v, w, Normal(3),NodalForce(3),DotProd, &
           NormalForce, SlipForce
       REAL(KIND=dp), ALLOCATABLE :: Basis(:)
       LOGICAL :: Stat, IsSlave, IsMaster
       TYPE(Matrix_t), POINTER :: ActiveProjector
       LOGICAL, ALLOCATABLE :: NodeDone(:)
       LOGICAL :: LinearContactLoads
       INTEGER :: i1,i2,j1,j2,ElemCode
       
       n = Mesh % MaxElementNodes
       ALLOCATE(Basis(n), Nodes % x(n), Nodes % y(n), Nodes % z(n) )

       CoordSys = CurrentCoordinateSystem()
       NodalForce = 0.0_dp

       NormalSign0 = 0
       NormalCount = 0
       
       ALLOCATE( NodeDone( Mesh % NumberOfNodes ) )
       NodeDone = .FALSE.

       LinearContactLoads = ListGetLogical( Model % Simulation,&
           'Contact BCs linear loads', Found )

       
100    DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         
         Element => Mesh % Elements( elem )         

         IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
         IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 

         IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
         
         Indexes => Element % NodeIndexes
         
         n = Element % TYPE % NumberOfNodes
         Nodes % x(1:n) = Mesh % Nodes % x(Indexes)
         Nodes % y(1:n) = Mesh % Nodes % y(Indexes)
         Nodes % z(1:n) = Mesh % Nodes % z(Indexes)
         
         IntegStuff = GaussPoints( Element, n )

         DO t=1,IntegStuff % n        
           U = IntegStuff % u(t)
           V = IntegStuff % v(t)
           W = IntegStuff % w(t)
           
           stat = ElementInfo( Element, Nodes, U, V, W, detJ, Basis )
           S = DetJ * IntegStuff % s(t)
           
           IF ( CoordSys /= Cartesian ) THEN
             X = SUM( Nodes % X(1:n) * Basis(1:n) )
             s = s * x
           END IF
           
           Normal = NormalVector( Element,Nodes,u,v,.TRUE. )

           ! Check the consistency of sign in the projector
           IF( IsSlave .AND. ( FlatProjector .OR. PlaneProjector .OR. NormalProjector ) ) THEN
             DotProd = SUM( Normal * ContactNormal ) 
             IF( DotProd < 0.0 ) THEN
               NormalSign = 1
             ELSE
               NormalSign = -1 
             END IF
             IF( NormalSign0 == 0 ) THEN
               NormalSign0 = NormalSign
             ELSE
               IF( NormalSign0 /= NormalSign ) NormalCount = NormalCount + 1
             END IF
           END IF

           DO i=1,n
             j = NormalLoadVar % Perm( Indexes(i) ) 

             IF( .NOT. NodeDone( Indexes(i) ) ) THEN             
               NodeDone( Indexes(i) ) = .TRUE.
               WeightVar % Values(j) = 0.0_dp
               NormalLoadVar % Values(j) = 0.0_dp
               SlipLoadVar % Values(j) = 0.0_dp
             END IF

             k = FieldPerm( Indexes(i) )
             IF( k == 0 .OR. j == 0 ) CYCLE
             DO l=1,dofs
               NodalForce(l) = LoadValues(dofs*(k-1)+l)
             END DO

             IF( ThisRotatedContact ) THEN
               NormalForce = NodalForce(1)
             ELSE
               NormalForce = SUM( NodalForce * Normal ) 
             END IF
             SlipForce = SQRT( SUM( NodalForce**2 ) - NormalForce**2 )

             NormalLoadVar % Values(j) = NormalLoadVar % Values(j) - &
                 s * Basis(i) * NormalForce
             SlipLoadVar % Values(j) = SlipLoadVar % Values(j) + &
                 s * Basis(i) * SlipForce
             
             WeightVar % Values(j) = WeightVar % Values(j) + s * Basis(i)
           END DO
           
         END DO
       END DO

       ! Normalize the computed normal loads such that the unit will be that of pressure
       DO i=1,Mesh % NumberOfNodes
         IF( NodeDone( i ) ) THEN             
           j = WeightVar % Perm(i)
           SlipLoadVar % Values(j) = SlipLoadVar % Values(j) / WeightVar % Values(j)**2
           NormalLoadVar % Values(j) = NormalLoadVar % Values(j) / WeightVar % Values(j)**2
         END IF
       END DO

       IF( LinearContactLoads ) THEN       
         DO elem=Mesh % NumberOfBulkElements + 1, &
             Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
           
           Element => Mesh % Elements( elem )         

           IsSlave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) 
           IsMaster = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag ) 
           IF( .NOT. ( IsSlave .OR. ( CreateDual .AND. IsMaster ) ) ) CYCLE
           
           Indexes => Element % NodeIndexes         
           ElemCode = Element % TYPE % ElementCode           

           SELECT CASE ( ElemCode )

           CASE( 408 )
             DO i=5,8
               i1=i-4
               i2=i1+1
               IF(i2==5) i2=1
               j = SlipLoadVar % Perm(Indexes(i))
               j1 = SlipLoadVar % Perm(Indexes(i1))
               j2 = SlipLoadVar % Perm(Indexes(i2))
               SlipLoadVar % Values(j) = 0.5_dp * &
                   ( SlipLoadVar % Values(j1) + SlipLoadVar % Values(j2))
               NormalLoadVar % Values(j) = 0.5_dp * &
                   ( NormalLoadVar % Values(j1) + NormalLoadVar % Values(j2))
             END DO
               
           CASE DEFAULT
             CALL Fatal(Caller,'Implement linear loads for: '//TRIM(I2S(ElemCode)))
           END SELECT
         END DO
       END IF
       
       IF( FlatProjector .OR. PlaneProjector .OR. NormalProjector ) THEN
         IF( NormalCount == 0 ) THEN
           CALL Info(Caller,'All normals are consistently signed',Level=10)
         ELSE
           CALL Warn(Caller,'There are normals with conflicting signs: '&
               //TRIM(I2S(NormalCount) ) )
           NormalSign = 1
         END IF
         CALL Info(Caller,'Normal direction for distance measure: '&
             //TRIM(I2S(NormalSign)),Level=8)
         DistSign = NormalSign 
       END IF

       ! Check whether the normal sign has been enforced
       IF( ListGetLogical( BC,'Normal Sign Negative',Found ) ) DistSign = -1
       IF( ListGetLogical( BC,'Normal Sign Positive',Found ) ) DistSign = 1


       DEALLOCATE( Basis, Nodes % x, Nodes % y, Nodes % z, NodeDone )

     END SUBROUTINE CalculateContactPressure

  

     ! Sets the contact in the normal direction by looking at the signed distance and 
     ! contact force. The initial contact set may be enlarged to eliminate null-space
     ! related to rigid-body motion.
     !----------------------------------------------------------------------------------
     SUBROUTINE NormalContactSet() 
       
       INTEGER :: LimitSign, Removed, Added
       REAL(KIND=dp) :: DistOffSet, MinLoad, MaxLoad, NodeLoad, MinDist, MaxDist, NodeDist
       INTEGER :: i,j,k,ind
       LOGICAL :: Found

       ! This is related to the formulation of the PDE and is probably fixed for all elasticity solvers
       LimitSign = -1

       Removed = 0
       Added = 0        
       MinLoad = HUGE(MinLoad)
       MaxLoad = -HUGE(MaxLoad)
       MinDist = HUGE(MinDist)
       MaxDist = -HUGE(MaxDist)

       Found = .FALSE.
       IF( FirstTime ) THEN
         DistOffset = ListGetCReal( BC,&
             'Mortar BC Initial Contact Depth',Found)
         IF(.NOT. Found ) DistOffset = ListGetCReal( BC,&
             'Contact Depth Offset Initial',Found)
       END IF
       IF( .NOT. Found ) DistOffset = ListGetCReal( BC,&
           'Contact Depth Offset',Found)

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = NormalLoadVar % Perm(j)

         ind = Dofs * (i-1) + DofN

         ! Tie contact should always be in contact - if we have found a counterpart
         IF( TieContact ) THEN
           MortarBC % Active(ind) = .TRUE.
           CYCLE
         END IF

         ! Enforce contact 
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Contact Active Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           MortarBC % Active(ind) = .TRUE.
           CYCLE
         END IF

         ! Enforce no contact
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Contact Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           MortarBC % Active(ind) = .FALSE.
           CYCLE
         END IF

         ! Free nodes with wrong sign in contact force
         !--------------------------------------------------------------------------       
         IF( MortarBC % Active( ind ) ) THEN
           NodeLoad = NormalLoadVar % Values(k)
           MaxLoad = MAX( MaxLoad, NodeLoad )
           MinLoad = MIN( MinLoad, NodeLoad )
           DoRemove = ( LimitSign * NodeLoad > LimitSign * LoadEps ) 
           IF( DoRemove .AND. ConservativeRemove ) THEN
             DoRemove = InterfaceDof(ind) 
           END IF
           IF( DoRemove ) THEN
             removed = removed + 1
             MortarBC % Active(ind) = .FALSE.
           END IF
         ELSE 
           NodeDist = DistVar % Values(k)
           MaxDist = MAX( MaxDist, NodeDist ) 
           MinDist = MIN( MinDist, NodeDist )

           DoAdd = ( NodeDist < -ValEps + DistOffset )            
           IF( DoAdd .AND. ConservativeAdd ) THEN
             DoAdd = InterfaceDof(ind)
           END IF
           IF( DoAdd ) THEN
             added = added + 1
             MortarBC % Active(ind) = .TRUE.
           END IF
         END IF
       END DO

       IF( InfoActive(20) ) THEN
         IF ( -HUGE(MaxDist) /= MaxDist ) THEN
           IF( MaxDist - MinDist >= 0.0_dp ) THEN
             PRINT *,'NormalContactSet Dist:',MinDist,MaxDist
           END IF
         END IF
         IF ( -HUGE(MaxLoad) /= MaxLoad) THEN
           IF( MaxLoad - MinLoad >= 0.0_dp ) THEN
             PRINT *,'NormalContactSet Load:',MinLoad,MaxLoad
           END IF
         END IF
       END IF

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the set'
         CALL Info(Caller,Message,Level=5)
       END IF
       
       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' nodes from the set'
         CALL Info(Caller,Message,Level=5)
       END IF

     END SUBROUTINE NormalContactSet



     ! If requested add new nodes to the contact set
     ! This would be typically done in order to make the elastic problem well defined
     ! Without any contact the bodies may float around.
     !---------------------------------------------------------------------------------
     SUBROUTINE IncreaseContactSet( LimitedMin )
       INTEGER :: LimitedMin

       REAL(KIND=dp), ALLOCATABLE :: DistArray(:)
       INTEGER, ALLOCATABLE :: IndArray(:)
       REAL(KIND=dp) :: Dist
       INTEGER :: i,j,ind,LimitedNow,NewNodes

       ! Nothing to do 
       IF( LimitedMin <= 0 ) RETURN

       LimitedNow = COUNT( MortarBC % active(DofN::Dofs) )      
       NewNodes = LimitedMin - LimitedNow
       IF( NewNodes <= 0 ) RETURN

       WRITE(Message,'(A,I0)') 'Initial number of contact nodes for '&
           //TRIM(VarName)//': ',LimitedNow 
       CALL Info(Caller,Message,Level=5)

       CALL Info(Caller,&
           'Setting '//TRIM(I2S(NewNodes))//' additional contact nodes',Level=5)

       ALLOCATE( DistArray( NewNodes ), IndArray( NewNodes ) ) 
       DistArray = HUGE( DistArray ) 
       IndArray = 0

       ! Find additional contact nodes from the closest non-contact nodes
       DO i = 1,Projector % NumberOfRows
         ind = Dofs * (i-1) + DofN
         IF( MortarBC % Active(ind)  ) CYCLE

         IF( Projector % InvPerm(i) == 0 ) CYCLE
         j = DistVar % Perm(Projector % InvPerm(i))
         Dist = DistVar % Values(j)
 
         IF( Dist < DistArray(NewNodes) ) THEN
           DistArray(NewNodes) = Dist
           IndArray(NewNodes) = i

           ! Order the new nodes such that the last node always has the largest distance
           ! This way we only need to compare to the one distance when adding new nodes.
           DO j=1,NewNodes-1
             IF( DistArray(j) > DistArray(NewNodes) ) THEN
               Dist = DistArray(NewNodes)
               DistArray(NewNodes) = DistArray(j)
               DistArray(j) = Dist                
               ind = IndArray(NewNodes)
               IndArray(NewNodes) = IndArray(j)
               IndArray(j) = ind                
             END IF
           END DO
         END IF
       END DO

       IF( ANY( IndArray == 0 ) ) THEN
         CALL Fatal(Caller,'Could not define sufficient number of new nodes!')
       END IF

       WRITE(Message,'(A,ES12.4)') 'Maximum distance needed for new nodes:',DistArray(NewNodes)
       CALL Info(Caller,Message,Level=8)

       MortarBC % Active( Dofs*(IndArray-1)+DofN ) = .TRUE.

       DEALLOCATE( DistArray, IndArray ) 

     END SUBROUTINE IncreaseContactSet



     ! Sets the contact in the tangent direction(s) i.e. the stick condition.
     ! Stick condition in 1st and 2nd tangent condition are always the same. 
     !----------------------------------------------------------------------------------
     SUBROUTINE TangentContactSet() 
       
       INTEGER :: Removed0, Removed, Added
       REAL(KIND=dp) :: NodeLoad, TangentLoad, mustatic, mudynamic, stickcoeff, &
           Fstatic, Fdynamic, Ftangent, du(3), Slip
       INTEGER :: i,j,k,l,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       IF( FrictionContact .AND. &
           ListGetLogical( BC,'Stick Contact Global',Found ) ) THEN
        
         ! Sum up global normal and slide forces
         DO i = 1,Projector % NumberOfRows
           j = Projector % InvPerm( i ) 
           IF( j == 0 ) CYCLE
           k = FieldPerm( j ) 
           IF( k == 0 ) CYCLE
           k = NormalLoadVar % Perm(j)
                      
           ! If there is no contact there can be no stick either
           indN = Dofs * (i-1) + DofN
           IF( .NOT. MortarBC % Active(indN) ) CYCLE

           NodeLoad = NormalLoadVar % Values(k)
           TangentLoad = SlipLoadVar % Values(k)
         
           mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
           mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )
           IF( mustatic <= mudynamic ) THEN
             CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
           END IF
           
           Fstatic = Fstatic + mustatic * ABS( NodeLoad ) 
           Fdynamic = Fdynamic + mudynamic * ABS( NodeLoad )
           Ftangent = Ftangent + ABS( TangentLoad ) 
           IF( Ftangent > Fstatic ) THEN
             SlipContact = .TRUE.
             FrictionContact = .FALSE.
           ELSE 
             GOTO 100
           END IF
         END DO
       END IF

       
       ! For stick and tie contact inherit the active flag from the normal component
       IF( SlipContact ) THEN
         MortarBC % Active( DofT1 :: Dofs ) = .FALSE.
         IF( Dofs == 3 ) THEN
            MortarBC % Active( DofT2 :: Dofs ) = .FALSE.
          END IF
          GOTO 100 
       ELSE IF( StickContact .OR. TieContact ) THEN
         MortarBC % Active( DofT1 :: Dofs ) = MortarBC % Active( DofN :: Dofs )
         IF( Dofs == 3 ) THEN
           MortarBC % Active( DofT2 :: Dofs ) = MortarBC % Active( DofN :: Dofs ) 
         END IF
         GOTO 100
       END IF

       CALL Info('TangentContactSet','Setting the stick set tangent components',Level=10)

       Removed0 = 0
       Removed = 0
       Added = 0        

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = NormalLoadVar % Perm(j)

         indN = Dofs * (i-1) + DofN
         indT1 = ind - DofN + DofT1
         IF(Dofs == 3 ) indT2 = ind - DofN + DofT2

         ! If there is no contact there can be no stick either
         IF( .NOT. MortarBC % Active(indN) ) THEN
           IF( MortarBC % Active(indT1) ) THEN
             removed0 = removed0 + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( Dofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
           CYCLE
         END IF

         ! Ok, we have normal contact what about stick
         ! Enforce stick condition
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Active Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( .NOT. MortarBC % Active(indT1) ) added = added + 1
           MortarBC % Active(indT1) = .TRUE.
           IF( Dofs == 3 ) MortarBC % Active(indT2) = .TRUE.
           CYCLE
         END IF

         ! Enforce no-stick condition (=slip)
         !------------------------------------------------------
         coeff = ListGetRealAtNode( BC,'Stick Passive Condition', j, Found )
         IF( Found .AND. coeff > 0.0_dp ) THEN
           IF( MortarBC % Active(IndT1) ) removed = removed + 1
           MortarBC % Active(indT1) = .FALSE.
           IF( Dofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           CYCLE
         END IF

         ! Remove nodes with too large tangent force
         !--------------------------------------------------------------------------       

         NodeLoad = NormalLoadVar % Values(k)
         TangentLoad = SlipLoadVar % Values(k)
         
         mustatic = ListGetRealAtNode( BC,'Static Friction Coefficient', j )
         mudynamic = ListGetRealAtNode( BC,'Dynamic Friction Coefficient', j )

         IF( mustatic <= mudynamic ) THEN
           CALL Warn('TangentContactSet','Static friction coefficient should be larger than dynamic!')
         END IF

         IF( MortarBC % Active(IndT1) ) THEN
           IF( TangentLoad > mustatic * ABS( NodeLoad ) ) THEN
             removed = removed + 1
             MortarBC % Active(indT1) = .FALSE.
             IF( Dofs == 3 ) MortarBC % Active(indT2) = .FALSE.
           END IF
         ELSE              
           stickcoeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j, Found )
           IF( Found ) THEN
             DO l=1,Dofs             
               du(l) = VeloVar % Values( Dofs*(k-1)+l ) 
             END DO
             IF( Dofs == 3 ) THEN
               Slip = SQRT(du(dofT1)**2 + du(DofT2)**2)
             ELSE
               Slip = ABS( du(dofT1) )
             END IF
             IF( stickcoeff * slip  < mudynamic * ABS( NodeLoad ) ) THEN
               added = added + 1
               MortarBC % Active(indT1) = .TRUE.
               IF( Dofs == 3 ) MortarBC % Active(indT2) = .TRUE.
             END IF
           END IF
         END IF
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' nodes to the stick set'
         CALL Info(Caller,Message,Level=5)
       END IF
       
       IF(removed0 > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed0,' non-contact nodes from the stick set'
         CALL Info(Caller,Message,Level=5)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' sliding nodes from the stick set'
         CALL Info(Caller,Message,Level=5)
       END IF


100    CALL Info(Caller,'Creating fields out of normal and stick contact sets',Level=10)

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         k = NormalActiveVar % Perm(j)

         IF( MortarBC % Active(Dofs*(i-1)+DofN) ) THEN
           NormalActiveVar % Values(k) = 1.0_dp
         ELSE
           NormalActiveVar % Values(k) = -1.0_dp
         END IF

         IF( MortarBC % Active(Dofs*(i-1)+DofT1) ) THEN
           StickActiveVar % Values(k) = 1.0_dp
         ELSE
           StickActiveVar % Values(k) = -1.0_dp
         END IF
       END DO

     END SUBROUTINE TangentContactSet



     ! Sets the diagonal entry for slip in the tangent direction(s).
     ! This coefficient may be used to relax the stick condition, and also to
     ! revert back nodes from slip to stick set. 
     !----------------------------------------------------------------------------------
     SUBROUTINE StickCoefficientSet() 
       
       REAL(KIND=dp) :: NodeLoad, TangentLoad
       INTEGER :: i,j,k,ind,IndN, IndT1, IndT2
       LOGICAL :: Found

       CALL Info('StickCoefficienttSet','Setting the stick coefficient entry for tangent components at stick',Level=10)

       ! Determine now whether we have contact or not
       DO i = 1,Projector % NumberOfRows
         j = Projector % InvPerm( i ) 
         IF( j == 0 ) CYCLE
         k = FieldPerm( j ) 
         IF( k == 0 ) CYCLE
         k = NormalLoadVar % Perm(j)

         indN = Dofs * (i-1) + DofN
         indT1 = Dofs * (i-1) + DofT1
         IF(Dofs == 3 ) indT2 = Dofs * (i-1) + DofT2

         IF( .NOT. MortarBC % Active(indN) ) THEN
           ! If there is no contact there can be no stick either
           coeff = 0.0_dp            
         ELSE IF( .NOT. MortarBC % Active(indT1) ) THEN
           ! If there is no stick there can be no stick coefficient either
           coeff = 0.0_dp
         ELSE
           ! Get the stick coefficient
           coeff = ListGetRealAtNode( BC,'Stick Contact Coefficient', j )
         END IF

         MortarBC % Diag(indT1) = coeff
         IF( Dofs == 3 ) MortarBC % Diag(indT2) = coeff
       END DO
       
     END SUBROUTINE StickCoefficientSet



     ! Here we eliminate the middle nodes from the higher order elements if they 
     ! are different than both nodes of which they are associated with.
     ! There is no way geometric information could be accurate enough to allow
     ! such contacts to exist.
     !---------------------------------------------------------------------------
     SUBROUTINE QuadraticContactSet()

       LOGICAL :: ElemActive(8)
       INTEGER :: i,j,k,n,added, removed, elem, elemcode, ElemInds(8)
       INTEGER, POINTER :: Indexes(:)

       added = 0
       removed = 0

       DO elem=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( elem )         

         IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         Indexes => Element % NodeIndexes         
         n = Element % TYPE % NumberOfNodes
         elemcode = Element % Type % ElementCode

         DO i=1,n
           ElemActive(i) = MortarBC % Active( ElemInds(i) ) 
           IF(j>0) THEN
             ElemInds(i) = Dofs * ( j - 1) + DofN
             ElemActive(i) = MortarBC % Active( ElemInds(i) ) 
           ELSE
             ElemActive(i) = .FALSE.
           END IF
         END DO

         SELECT CASE ( elemcode ) 

         CASE( 202, 303, 404 ) 
           CONTINUE

         CASE( 203 )
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(3) ) ) THEN
             MortarBC % Active( ElemInds(3) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE( 306 ) 
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(4) ) ) THEN
             MortarBC % Active( ElemInds(4) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(2) .EQV. ElemActive(3) ) &
               .AND. ( ElemActive(2) .NEQV. ElemActive(5) ) ) THEN
             MortarBC % Active( ElemInds(5) ) = ElemActive(2) 
             IF( ElemActive(2) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(3) .EQV. ElemActive(1) ) &
               .AND. ( ElemActive(3) .NEQV. ElemActive(6) ) ) THEN
             MortarBC % Active( ElemInds(6) ) = ElemActive(3) 
             IF( ElemActive(3) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE( 408 ) 
           IF( ( ElemActive(1) .EQV. ElemActive(2) ) &
               .AND. ( ElemActive(1) .NEQV. ElemActive(5) ) ) THEN
             MortarBC % Active( ElemInds(5) ) = ElemActive(1) 
             IF( ElemActive(1) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(2) .EQV. ElemActive(3) ) &
               .AND. ( ElemActive(2) .NEQV. ElemActive(6) ) ) THEN
             MortarBC % Active( ElemInds(6) ) = ElemActive(2) 
             IF( ElemActive(2) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(3) .EQV. ElemActive(4) ) &
               .AND. ( ElemActive(3) .NEQV. ElemActive(7) ) ) THEN
             MortarBC % Active( ElemInds(7) ) = ElemActive(3) 
             IF( ElemActive(3) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

           IF( ( ElemActive(4) .EQV. ElemActive(1) ) &
               .AND. ( ElemActive(4) .NEQV. ElemActive(8) ) ) THEN
             MortarBC % Active( ElemInds(8) ) = ElemActive(4) 
             IF( ElemActive(4) ) THEN
               added = added + 1
             ELSE
               removed = removed + 1
             END IF
           END IF

         CASE DEFAULT
           CALL Fatal('NormalContactSet','Cannot deal with element: '//TRIM(I2S(elemcode)))

         END SELECT
       END DO

       IF(added > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Added ',added,' quadratic nodes to contact set'
         CALL Info(Caller,Message,Level=5)
       END IF

       IF(removed > 0) THEN
         WRITE(Message,'(A,I0,A)') 'Removed ',removed,' quadratic nodes from contact set'
         CALL Info(Caller,Message,Level=5)
       END IF
         
     END SUBROUTINE QuadraticContactSet


     ! Project contact fields from slave to master
     !----------------------------------------------------------------------------------------
     SUBROUTINE ProjectFromSlaveToMaster()

       REAL(KIND=dp) :: Disp(3), Coord(3), PrevDisp(3), Velo(3), ContactVelo(3), &
           LocalNormal0(3), SlipCoord(3)
       REAL(KIND=dp), POINTER :: DispVals(:), PrevDispVals(:) 
       REAL(KIND=dp) :: MinDist, MaxDist, CoeffEps
       LOGICAL, ALLOCATABLE :: SlaveNode(:), NodeDone(:)
       REAL(KIND=dp), ALLOCATABLE :: CoeffTable(:), RealActive(:)
       INTEGER :: i,j,k,l,l2

       CALL Info(Caller,'Mapping entities from slave to master',Level=10)

       ALLOCATE( SlaveNode( Mesh % NumberOfNodes ) ) 
       SlaveNode = .FALSE.

       DO i=Mesh % NumberOfBulkElements + 1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

         Element => Mesh % Elements( i )                  
         IF( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag ) THEN
           SlaveNode( Element % NodeIndexes ) = .TRUE.
         END IF
       END DO
 
       n = SIZE( DistVar % Values )
       ALLOCATE( CoeffTable( n ), NodeDone( n ) )
           

       CoeffTable = 0.0_dp
       NodeDone = .FALSE.
       

       DO i = 1,Projector % NumberOfRows             
         
         IF( Projector % InvPerm(i) == 0 ) CYCLE          
         l = DistVar % Perm( Projector % InvPerm(i) )
         
         DO j = Projector % Rows(i),Projector % Rows(i+1)-1
           k = Projector % Cols(j)
           
           IF( FieldPerm( k ) == 0 ) CYCLE               
           IF( SlaveNode( k ) ) CYCLE
           
           coeff = Projector % Values(j)
           
           l2 = DistVar % Perm( k )                
           
           IF(.NOT. NodeDone( l2 ) ) THEN
             DistVar % Values( l2 ) = 0.0_dp
             GapVar % Values( l2 ) = 0.0_dp
             NormalActiveVar % Values( l2 ) = 0.0_dp
             StickActiveVar % Values( l2 ) = 0.0_dp
             NormalLoadVar % Values( l2 ) = 0.0_dp
             SlipLoadVar % Values( l2 ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,Dofs             
                 VeloVar % Values( Dofs*(l2-1)+k ) = 0.0_dp
               END DO
             END IF
             NodeDone( l2 ) = .TRUE.
           END IF

           CoeffTable( l2 ) = CoeffTable( l2 ) + coeff           
           DistVar % Values( l2 ) = DistVar % Values( l2 ) + coeff * DistVar % Values( l ) 
           GapVar % Values( l2 ) = GapVar % Values( l2 ) + coeff * GapVar % Values( l )
           NormalActiveVar % Values( l2 ) = NormalActiveVar % Values( l2 ) + coeff * NormalActiveVar % Values( l ) 
           StickActiveVar % Values( l2 ) = StickActiveVar % Values( l2 ) + coeff * StickActiveVar % Values( l ) 
           NormalLoadVar % Values( l2 ) = NormalLoadVar % Values( l2 ) + coeff * NormalLoadVar % Values( l ) 
           SlipLoadVar % Values( l2 ) = SlipLoadVar % Values( l2 ) + coeff * SlipLoadVar % Values( l ) 
           IF( CalculateVelocity ) THEN
             DO k=1,Dofs             
               VeloVar % Values( Dofs*(l2-1)+k ) = VeloVar % Values( Dofs*(l2-1)+k ) + &
                   coeff * VeloVar % Values( Dofs*(l-1)+k)
             END DO
           END IF
         END DO
       END DO
       
       CoeffEps = 1.0d-8 * MAXVAL( ABS( CoeffTable ) )
       DO i=1,SIZE( CoeffTable )            
         IF( NodeDone( i ) .AND. ( ABS( CoeffTable(i) ) > CoeffEps ) ) THEN
           DistVar % Values( i ) = DistVar % Values( i ) / CoeffTable( i ) 
           GapVar % Values( i ) = GapVar % Values( i ) / CoeffTable( i ) 
           NormalActiveVar % Values( i ) = NormalActiveVar % Values( i ) / CoeffTable( i ) 
           StickActiveVar % Values( i ) = StickActiveVar % Values( i ) / CoeffTable( i ) 

           IF( NormalActiveVar % Values( i ) >= 0.0_dp ) THEN
             NormalLoadVar % Values( i ) = NormalLoadVar % Values( i ) / CoeffTable( i ) 
             SlipLoadVar % Values( i ) = SlipLoadVar % Values( i ) / CoeffTable( i ) 
             IF( CalculateVelocity ) THEN
               DO k=1,Dofs
                 VeloVar % Values( Dofs*(i-1)+k ) = VeloVar % Values( Dofs*(i-1)+k ) / CoeffTable( i ) 
               END DO
             END IF
           ELSE
             NormalLoadVar % Values( i ) = 0.0_dp
             SlipLoadVar % Values( i ) = 0.0_dp
             IF( CalculateVelocity ) THEN
               DO k=1,Dofs
                 VeloVar % Values( Dofs*(i-1)+k ) = 0.0_dp
               END DO
             END IF             
           END IF

         END IF
       END DO

       DO i = 1, Projector % NumberOfRows
         j = Projector % InvPerm(i)
         IF( j == 0 ) CYCLE
         k = NormalActiveVar % Perm(j)
         
         IF( NormalActiveVar % Values( k ) < 0.0_dp ) THEN
           IF( CalculateVelocity ) THEN
             DO l=1,Dofs
               VeloVar % Values( Dofs*(k-1)+l ) = 0.0_dp
             END DO
           END IF
         END IF
       END DO


     END SUBROUTINE ProjectFromSlaveToMaster
   


     ! Set the friction in an implicit manner by copying matrix rows of the normal component
     ! to matrix rows of the tangential component multiplied by friction coefficient and 
     ! direction vector. 
     !---------------------------------------------------------------------------------------
     SUBROUTINE SetSlideFriction()

       REAL(KIND=dp), POINTER :: Values(:)
       LOGICAL, ALLOCATABLE :: NodeDone(:)
       REAL(KIND=dp) :: Coeff, ActiveLimit
       TYPE(Element_t), POINTER :: Element
       INTEGER, POINTER :: NodeIndexes(:)
       INTEGER :: i,j,k,k2,k3,l,l2,l3,n,t
       TYPE(Matrix_t), POINTER :: A
       LOGICAL :: Slave, Master, GivenDirection
       REAL(KIND=dp), POINTER :: VeloDir(:,:)
       REAL(KIND=dp) :: VeloCoeff(3),AbsVeloCoeff
       INTEGER :: VeloSign = 1


       IF(.NOT. ListCheckPresent( BC, 'Dynamic Friction Coefficient') ) RETURN
      
       CALL Info(Caller,'Setting contact friction for boundary',Level=10)

       GivenDirection = ListCheckPresent( BC,'Contact Velocity')
       IF(.NOT. GivenDirection ) THEN
         IF(.NOT. ASSOCIATED( VeloVar ) ) THEN
           CALL Fatal(Caller,'Contact velocity must be given in some way')
         END IF
       END IF

       ActiveLimit = 0.0_dp
      
       Values => Solver % Matrix % values              
       ALLOCATE( NodeDone( SIZE( FieldPerm ) ) )
       A => Solver % Matrix        

       NodeDone = .FALSE.
       Coeff = 0.0_dp

       DO t = Mesh % NumberOfBulkElements+1, &
           Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
         Element => Mesh % Elements(t)
         
         Model % CurrentElement => Element
         
         Slave = ( Element % BoundaryInfo % Constraint == Model % BCs(bc_ind) % Tag )
         Master = ( Element % BoundaryInfo % Constraint == Model % BCs(master_ind) % Tag )
         
         IF( .NOT. ( Slave .OR. Master ) ) CYCLE

         NodeIndexes => Element % NodeIndexes
         n = Element % TYPE % NumberOfNodes
         
         DO i = 1, n
           j = Nodeindexes(i) 

           IF( NodeDone( j ) ) CYCLE
           IF( FieldPerm( j ) == 0 ) CYCLE

           ! Skipping the nodes not in the boundary
           k = NormalActiveVar % Perm( j )
           IF( k == 0 ) CYCLE
           
           ! Skipping the nodes not in contact. 
           IF( NormalActiveVar % Values( k ) <= -ActiveLimit ) CYCLE
           
           ! skipping the nodes in tangent stick
           IF( StickActiveVar % Values( k ) >= ActiveLimit ) CYCLE

           NodeDone( j ) = .TRUE.

           IF( Slave ) THEN
             Coeff = ListGetRealAtNode( BC,& 
                 'Dynamic Friction Coefficient', j, Found )
           ELSE
             Coeff = ListGetRealAtNode( MasterBC,& 
                 'Dynamic Friction Coefficient', j, Found )
             ! If friction not found in master then use the friction coefficient of the slave 
             ! Ideally they should be the same. 
             IF(.NOT. Found ) THEN
               Coeff = ListGetRealAtNode( BC,& 
                   'Dynamic Friction Coefficient', j, Found )               
             END IF
           END IF

           ! There is no point of setting too small friction coefficient
           IF(ABS(Coeff) < 1.0d-10) CYCLE

           IF( ThisRotatedContact ) THEN
             Rotated = GetSolutionRotation(NTT, j )
             LocalNormal = NTT(:,1)
             LocalT1 = NTT(:,2)
             IF( Dofs == 3 ) LocalT2 = NTT(:,3)
           ELSE
             Rotated = .FALSE.
             LocalNormal = ContactNormal
             LocalT1 = ContactT1
             IF( Dofs == 3 ) LocalT2 = ContactT2 
           END IF
           
           VeloCoeff = 0.0_dp           
           VeloSign = 1

           IF( GivenDirection ) THEN
             IF( Slave ) THEN
               VeloDir => ListGetConstRealArray( BC, &
                   'Contact Velocity', Found)
             ELSE
               VeloDir => ListGetConstRealArray( MasterBC, &
                   'Contact Velocity', Found)
               IF(.NOT. Found ) THEN
                 ! If velocity direction not found in master then use the opposite velocity of the slave
                 VeloDir => ListGetConstRealArray( BC, &
                     'Contact Velocity', Found)
                 VeloSign = -1
               END IF
             END IF
             VeloCoeff(DofT1) = SUM( VeloDir(1:3,1) * LocalT1 )
             IF( Dofs == 3 ) THEN
               VeloCoeff(DofT2) = SUM( VeloDir(1:3,1) * LocalT2 )
             END IF
           ELSE
             VeloCoeff(DofT1) = VeloVar % Values(Dofs*(k-1)+DofT1) 
             IF(Dofs==3) VeloCoeff(DofT2) = VeloVar % Values(Dofs*(k-1)+DofT2) 
             IF( .NOT. Slave .AND. .NOT. Rotated ) THEN
               VeloSign = -1
             END IF
           END IF

           ! Normalize coefficient to unity so that it only represents the direction of the force

           AbsVeloCoeff = SQRT( SUM( VeloCoeff**2 ) )
           IF( AbsVeloCoeff > TINY(AbsVeloCoeff) ) THEN
             VeloCoeff = VeloSign * VeloCoeff / AbsVeloCoeff
           ELSE
             CYCLE
           END IF

           ! Add the friction coefficient 
           VeloCoeff = Coeff * VeloCoeff 

           j = FieldPerm( j ) 
           k = DOFs * (j-1) + DofN 

           k2 = DOFs * (j-1) + DofT1 
           A % Rhs(k2) = A % Rhs(k2) - VeloCoeff(DofT1) * A % Rhs(k)

           IF( Dofs == 3 ) THEN
             k3 = DOFs * (j-1) + DofT2
             A % Rhs(k3) = A % Rhs(k3) - VeloCoeff(DofT2) * A % Rhs(k)             
           END IF

           DO l = A % Rows(k),A % Rows(k+1)-1
             DO l2 = A % Rows(k2), A % Rows(k2+1)-1
               IF( A % Cols(l2) == A % Cols(l) ) EXIT
             END DO

             A % Values(l2) = A % Values(l2) - VeloCoeff(DofT1) * A % Values(l)
             
             IF( Dofs == 3 ) THEN
               DO l3 = A % Rows(k3), A % Rows(k3+1)-1
                 IF( A % Cols(l3) == A % Cols(l) ) EXIT
               END DO
               A % Values(l3) = A % Values(l3) - VeloCoeff(DofT2) * A % Values(l)
             END IF
           END DO
         END DO
       END DO
       
       n = COUNT( NodeDone ) 
       CALL Info('SetSlideFriction','Number of friction nodes: '//TRIM(I2S(n)),Level=10)
       
       DEALLOCATE( NodeDone )

     END SUBROUTINE SetSlideFriction
     

   END SUBROUTINE DetermineContact


 END MODULE ContactUtils
