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
! *  Utilities for setting Dirichlet conditions. 
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
! *  Original Date: 28 Sep 1998
! *
! *****************************************************************************/

!> Utilities for setting Dirichlet conditions.
!------------------------------------------------------------------------------

!> \ingroup ElmerLib
!> \{

MODULE DirichletUtils

#include "../config.h"
  
  USE Types
  USE ElementUtils
  USE MeshUtils
  USE AssemblyUtils
  USE ListMatrix
  USE CRSMatrix
  USE ScalingUtils
  USE ParallelUtils
   
  IMPLICIT NONE


CONTAINS 


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

  

!> Sets one Dirichlet condition to the desired value
!------------------------------------------------------------------------------
   SUBROUTINE UpdateDirichletDof( A, dof, dval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    REAL(KIND=dp) :: dval

    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT. ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % Dvalues( dof ) = dval
    A % ConstrainedDOF( dof ) = .TRUE.
    
  END SUBROUTINE UpdateDirichletDof
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
   SUBROUTINE UpdateDirichletDofC( A, dof, cval )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    COMPLEX(KIND=dp) :: cval

    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT. ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % Dvalues( 2*dof-1 ) = REAL( cval )
    A % ConstrainedDOF( 2*dof-1 ) = .TRUE.

    A % Dvalues( 2*dof ) = AIMAG( cval )
    A % ConstrainedDOF( 2*dof ) = .TRUE.
    
  END SUBROUTINE UpdateDirichletDofC
!------------------------------------------------------------------------------



  
!> Releases one Dirichlet condition 
!------------------------------------------------------------------------------
   SUBROUTINE ReleaseDirichletDof( A, dof )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: A
    INTEGER :: dof
    REAL(KIND=dp) :: dval
      
    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
      A % ConstrainedDOF = .FALSE.
    END IF
    
    IF(.NOT.ALLOCATED(A % Dvalues)) THEN
      ALLOCATE(A % Dvalues(A % NumberOfRows))
      A % Dvalues = 0._dp
    END IF
    
    A % ConstrainedDOF( dof ) = .FALSE.
    
  END SUBROUTINE ReleaseDirichletDof
!------------------------------------------------------------------------------


  
!> Release the range or min/max values of Dirichlet values.
!------------------------------------------------------------------------------
  FUNCTION DirichletDofsRange( Solver, Oper ) RESULT ( val ) 
!------------------------------------------------------------------------------
    TYPE(Solver_t), OPTIONAL :: Solver
    CHARACTER(LEN=*), OPTIONAL :: Oper 
    REAL(KIND=dp) :: val
    
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: minv,maxv
    LOGICAL :: FindMin, FindMax
    INTEGER :: i,OperNo
    
    IF( PRESENT( Solver ) ) THEN
      A => Solver % Matrix
    ELSE
      A => CurrentModel % Solver % Matrix
    END IF
    
    val = 0.0_dp
    
    ! Defaulting to range
    OperNo = 0

    IF( PRESENT( Oper ) ) THEN
      IF( Oper == 'range' ) THEN
        OperNo = 0
      ELSE IF( Oper == 'min' ) THEN
        OperNo = 1 
      ELSE IF( Oper == 'max' ) THEN
        OperNo = 2
      ELSE
        CALL Fatal('DirichletDofRange','Unknown operator: '//TRIM(Oper))
      END IF
    END IF
          
    IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
      RETURN
    END IF
  
    IF( OperNo == 0 .OR. OperNo == 1 ) THEN
      minv = HUGE( minv ) 
      DO i=1,SIZE( A % ConstrainedDOF )
        IF( A % ConstrainedDOF(i) ) minv = MIN( A % DValues(i), minv ) 
      END DO
      minv = ParallelReduction( minv, 1 ) 
    END IF

    IF( OperNo == 0 .OR. OperNo == 2 ) THEN
      maxv = -HUGE( maxv ) 
      DO i=1,SIZE( A % ConstrainedDOF )
        IF( A % ConstrainedDOF(i) ) maxv = MAX( A % DValues(i), maxv ) 
      END DO
      maxv = ParallelReduction( maxv, 2 ) 
    END IF
    
    IF( OperNo == 0 ) THEN    
      val = maxv - minv
    ELSE IF( OperNo == 1 ) THEN
      val = minv
    ELSE
      val = maxv
    END IF
      
  END FUNCTION DirichletDofsRange
!------------------------------------------------------------------------------
  

!------------------------------------------------------------------------------
!> Set dirichlet boundary condition for given dof. The conditions are
!> set based on the given name and applied directly to the matrix structure
!> so that a row is zeroed except for the diagonal which is set to one. 
!> Then the r.h.s. value determines the value of the field variable 
!> in the solution of the linear system.
!------------------------------------------------------------------------------
   SUBROUTINE SetDirichletBoundaries( Model, A, b, Name, DOF, NDOFs, Perm, &
       PermOffSet, OffDiagonalMatrix )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model          !< The current model structure
    TYPE(Matrix_t), POINTER :: A    !< The global matrix
    REAL(KIND=dp) :: b(:)           !< The global RHS vector
    CHARACTER(LEN=*) :: Name        !< Name of the dof to be set
    INTEGER :: DOF                  !< The order number of the dof
    INTEGER :: NDOFs                !< The total number of DOFs for this equation
    INTEGER :: Perm(:)              !< The node reordering info, this has been generated at the beginning of the 
                                    !< simulation for bandwidth optimization
    INTEGER, OPTIONAL :: PermOffSet  !< If the matrix and permutation vectors are not in sync the offset may used as a remedy. 
                                     !< Needed in fully coupled systems.
    LOGICAL, OPTIONAL :: OffDiagonalMatrix  !< For block systems the only the diagonal matrix should be given non-zero 
                                            !< entries for matrix and r.h.s., for off-diagonal matrices just set the row to zero.
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:), IndNodes(:), BCOrder(:)
    INTEGER, ALLOCATABLE :: Indexes(:), PassPerm(:)
    INTEGER :: BC,i,j,j2,k,l,m,n,nd,p,t,k1,k2,OffSet
    LOGICAL :: GotIt, periodic, OrderByBCNumbering, ReorderBCs
    REAL(KIND=dp), POINTER :: MinDist(:)
    REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
    REAL(KIND=dp) ::  s

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver

    LOGICAL :: Conditional
    LOGICAL, ALLOCATABLE :: DonePeriodic(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: CondName, DirName, PassName, PassCondName

    INTEGER :: NoNodes,NoDims,bf_id,nlen, NOFNodesFound, dim, &
        bndry_start, bndry_end, Upper
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), Condition(:), Work(:)!,DiagScaling(:)
    REAL(KIND=dp) :: GlobalMinDist,Dist, Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActiveCond(:), ActivePartAll(:)
    TYPE(ValueList_t), POINTER :: ValueList, Params
    LOGICAL :: NodesFound, Passive, OffDiagonal, ApplyLimiter
    LOGICAL, POINTER :: LimitActive(:)
    TYPE(Variable_t), POINTER :: Var

    TYPE(Element_t), POINTER :: Parent

    INTEGER :: ind, ElemFirst, ElemLast, bf, BCstr, BCend, BCinc
    REAL(KIND=dp) :: SingleVal
    LOGICAL :: AnySingleBC, AnySingleBF
    LOGICAL, ALLOCATABLE :: LumpedNodeSet(:)
    LOGICAL :: NeedListMatrix
    INTEGER, ALLOCATABLE :: Rows0(:), Cols0(:)
    REAL(KIND=dp), POINTER :: BulkValues0(:)
    INTEGER :: DirCount
    CHARACTER(*), PARAMETER :: Caller = 'SetDirichletBoundaries'
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER, POINTER :: PlaneInds(:)
    
!------------------------------------------------------------------------------
! These logical vectors are used to minimize extra effort in setting up different BCs
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    n = MAX( Model % NumberOfBodyForces,Model % NumberOfBCs)
    IF( n == 0 ) THEN
      CALL Info(Caller,'No BCs or Body Forces present, exiting early...',Level=12)
    END IF

    ALLOCATE( ActivePart(n), ActivePartAll(n), ActiveCond(n))
    CondName = Name(1:nlen) // ' Condition'
    PassName = Name(1:nlen) // ' Passive'
    PassCondName = Name(1:nlen) // ' Condition' // ' Passive'

    OffSet = 0
    OffDiagonal = .FALSE.
    IF( PRESENT( PermOffSet) ) OffSet = PermOffSet
    IF( PRESENT( OffDiagonalMatrix ) ) OffDiagonal = OffDiagonalMatrix

    Mesh => Model % Mesh
    ALLOCATE( Indexes(Mesh % MaxElementDOFs) )
!------------------------------------------------------------------------------
! Go through the periodic BCs and set the linear dependence
!------------------------------------------------------------------------------

   ActivePart = .FALSE.
   DO BC=1,Model % NumberOfBCs
     IF ( ListGetLogical( Model % BCs(BC) % Values, &
         'Periodic BC ' // Name(1:nlen), GotIt ) ) ActivePart(BC) = .TRUE.
     IF ( ListGetLogical( Model % BCs(BC) % Values, &
         'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) ActivePart(BC) = .TRUE.
     IF ( ListCheckPresent( Model % BCs(BC) % Values, &
         'Periodic BC Scale ' // Name(1:nlen) ) ) ActivePart(BC) = .TRUE.
   END DO
   
   IF( ANY(ActivePart) ) THEN    
     IF( Offset > 0 ) THEN
       CALL Fatal(Caller,'Periodicity not considered with offset')
     END IF

     ALLOCATE( DonePeriodic( Mesh % NumberOFNodes ) )
     DonePeriodic = .FALSE.
     DO BC=1,Model % NumberOfBCs
       IF( ActivePart(BC) ) THEN
         CALL SetPeriodicBoundariesPass1( Model, A, b, Name, DOF, &
             NDOFs, Perm, BC, DonePeriodic )
       END IF
     END DO
     
     DonePeriodic = .FALSE.
     DO BC=1,Model % NumberOfBCs
       IF(ActivePart(BC)) THEN       
         CALL SetPeriodicBoundariesPass2( Model, A, b, Name, DOF, &
             NDOFs, Perm, BC, DonePeriodic )
       END IF
     END DO

     IF( InfoActive(12) ) THEN
       CALL Info(Caller,'Number of periodic points set: '&
           //TRIM(I2S(COUNT(DonePeriodic))),Level=12)
     END IF

     DEALLOCATE( DonePeriodic ) 

   END IF
   

! Add the possible friction coefficient
!----------------------------------------------------------
   IF ( ListCheckPresentAnyBC( Model,'Friction BC ' // Name(1:nlen) ) ) THEN
     CALL SetFrictionBoundaries( Model, A, b, Name, NDOFs, Perm )
   END IF


! Add the possible nodal jump in case of mortar projectors
!---------------------------------------------------------------
   IF( ListGetLogical( Model % Solver % Values,'Apply Mortar BCs',GotIt ) ) THEN
     CALL SetWeightedProjectorJump( Model, A, b, &
                      Name, DOF, NDOFs, Perm )
   END IF


!------------------------------------------------------------------------------
! Go through the normal Dirichlet BCs applied on the boundaries
!------------------------------------------------------------------------------

    ActivePart = .FALSE.
    ActiveCond = .FALSE.
    ActivePartAll = .FALSE.
    DO BC=1,Model % NumberOfBCs
      ActivePartAll(BC) = ListCheckPresent( &
            Model % BCs(bc) % Values, Name(1:nlen) // ' DOFs' )
      ActivePart(BC) = ListCheckPresent( Model % BCs(bc) % Values, Name ) 
      ActiveCond(BC) = ListCheckPresent( Model % BCs(bc) % Values, CondName )      
    END DO

    OrderByBCNumbering = ListGetLogical( Model % Simulation, &
       'Set Dirichlet BCs by BC Numbering', gotIt)

    BCOrder => ListGetIntegerArray( Model % Solver % Values, &
         'Dirichlet BC Order', ReorderBCs)
    IF(ReorderBCs) THEN
       IF(.NOT. OrderByBCNumbering) THEN
          CALL Warn(Caller,"Requested 'Dirichlet BC Order' but &
               &not 'Set Dirichlet BCs by BC Numbering', ignoring...")
       ELSE IF(SIZE(BCOrder) /= Model % NumberOfBCs) THEN
          CALL Fatal(Caller,"'Dirichlet BC Order' is the wrong length!")
       END IF
    END IF

    bndry_start = Model % NumberOfBulkElements+1
    bndry_end   = bndry_start+Model % NumberOfBoundaryElements-1
    DirCount = 0

    ! check and set some flags for nodes belonging to n-t boundaries
    ! getting set by other bcs:
    ! --------------------------------------------------------------
    IF ( NormalTangentialNOFNodes>0 ) THEN
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                   Model % BCs(BC) % Tag ) CYCLE

            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element

            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
          
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                 Model % BCs(BC) % Tag ) CYCLE
          
            ValueList => Model % BCs(BC) % Values
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            CALL CheckNTelement(n,t)
          END DO
        END DO
      END IF

      IF ( DOF<= 0 ) THEN
        DO t=bndry_start,bndry_end
          Element => Model % Elements(t)
          n = Element % TYPE % NumberOfNodes
          DO j=1,n
            k = BoundaryReorder(Element % NodeIndexes(j))
            IF (k>0) THEN
              NTelement(k,:)=0
              NTzeroing_done(k,:) = .FALSE.
            END IF
          END DO
        END DO
      END IF
    END IF

    
    ! Set the Dirichlet BCs from active boundary elements, if any...:
    !----------------------------------------------------------------
    IF( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN    
      IF ( OrderByBCNumbering ) THEN
        DO i=1,Model % NumberOfBCs
          BC = i
          IF(ReorderBCs) BC = BCOrder(BC)
          IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
          Conditional = ActiveCond(BC)

          DO t = bndry_start, bndry_end
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                Model % BCs(BC) % Tag ) CYCLE
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p) == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n,t)
          END DO
        END DO
      ELSE
        DO t = bndry_start, bndry_end
          DO BC=1,Model % NumberOfBCs
            IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE
            Conditional = ActiveCond(BC)
            
            Element => Model % Elements(t)
            IF ( Element % BoundaryInfo % Constraint /= &
                Model % BCs(BC) % Tag ) CYCLE
            
            Model % CurrentElement => Element
            IF ( ActivePart(BC) ) THEN
              n = Element % TYPE % NumberOfNodes
              IF ( Model % Solver % DG ) THEN
                 Parent => Element % BoundaryInfo % Left
                 DO p=1,Parent % Type % NumberOfNodes
                   DO j=1,n
                      IF (Parent % NodeIndexes(p)  == Element % NodeIndexes(j) ) THEN
                        Indexes(j) = Parent % DGIndexes(p); EXIT
                      END IF
                   END DO
                 END DO
              ELSE
                Indexes(1:n) = Element % NodeIndexes
              END IF
            ELSE
              n = SgetElementDOFs( Indexes )
            END IF
            ValueList => Model % BCs(BC) % Values
            CALL SetElementValues(n,t)
          END DO
        END DO
      END IF
    END IF


!------------------------------------------------------------------------------
! Go through the Dirichlet conditions in the body force lists
!------------------------------------------------------------------------------
    
    ActivePart = .FALSE.
    ActiveCond = .FALSE.
    ActivePartAll = .FALSE.
    Passive = .FALSE.
    DO bf_id=1,Model % NumberOFBodyForces
      ValueList => Model % BodyForces(bf_id) % Values

      ActivePartAll(bf_id) = ListCheckPresent(ValueList, Name(1:nlen) // ' DOFs' ) 
      ActiveCond(bf_id) = ListCheckPresent( ValueList,CondName )      
      ActivePart(bf_id) = ListCheckPresent(ValueList, Name(1:nlen) ) 

      Passive = Passive .OR. ListCheckPresent(ValueList, PassName)
    END DO
    
    IF ( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN
      Solver => Model % Solver
      Mesh   => Solver % Mesh

      ALLOCATE(PassPerm(Mesh % NumberOfNodes),NodeIndexes(1));PassPerm=0
      DO i=0,Mesh % PassBCCnt-1
        j=Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements-i
        PassPerm(Mesh % Elements(j) % NodeIndexes)=1
      END DO

      DO t = 1, Mesh % NumberOfBulkElements 
        Element => Mesh % Elements(t)
        IF( Element % BodyId <= 0 .OR. Element % BodyId > Model % NumberOfBodies ) THEN
          CALL Warn(Caller,'Element body id beyond body table!')
          CYCLE
        END IF
                    
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id)) CYCLE
        Conditional = ActiveCond(bf_id)

        Model % CurrentElement => Element
        
        IF ( ActivePart(bf_id) ) THEN
          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes
        ELSE
          n = SgetElementDOFs( Indexes )
        END IF

        ValueList => Model % BodyForces(bf_id) % Values
        IF(.NOT. ASSOCIATED( ValueList ) ) CYCLE
        
        IF (ListGetLogical(ValueList,PassCondName,GotIt)) THEN
          IF (.NOT.CheckPassiveElement(Element)) CYCLE
          DO j=1,n
            NodeIndexes(1) = Indexes(j)
            IF(PassPerm(NodeIndexes(1))==0) CALL SetPointValues(1)
          END DO
        ELSE
          CALL SetElementValues( n,t )
        END IF
        
      END DO
      
      DEALLOCATE(NodeIndexes,PassPerm)
    END IF
    
    DEALLOCATE(ActivePart, ActiveCond)

    
!------------------------------------------------------------------------------
! Go through the pointwise Dirichlet BCs that are created on-the-fly
! Note that it is best that the coordinates are transformed to nodes using 
! the right variable. Otherwise it could point to nodes that are not active.
!------------------------------------------------------------------------------
     
    DO BC=1,Model % NumberOfBCs
      
      ValueList => Model % BCs(BC) % Values
      IF( .NOT. ListCheckPresent( ValueList,Name )) CYCLE
      NodesFound = ListCheckPresent( ValueList,'Target Nodes' )

      ! The coordinates are only requested for a body that has no list of nodes.
      ! At the first calling the list of coordinates is transformed to list of nodes.
      IF(.NOT. NodesFound) THEN
        CoordNodes => ListGetConstRealArray(ValueList,'Target Coordinates',GotIt)
        IF(GotIt) THEN
          Eps = ListGetConstReal( ValueList, 'Target Coordinates Eps', Gotit )
          IF ( .NOT. GotIt ) THEN
            Eps = HUGE(Eps)
          ELSE
            ! We are looking at square of distance
            Eps = Eps**2
          END IF

          NoNodes = SIZE(CoordNodes,1)
          NoDims = SIZE(CoordNodes,2)
          
          IF(NoNodes > 0) THEN               
            ALLOCATE( IndNodes(NoNodes), MinDist(NoNodes) )
            IndNodes = -1
            MinDist = HUGE( Dist )
            DO j=1,NoNodes
              DO i=1,Model % NumberOfNodes
                IF( Perm(i) == 0) CYCLE
                
                Dist = (Mesh % Nodes % x(i) - CoordNodes(j,1))**2 
                IF(NoDims >= 2) Dist = Dist + (Mesh % Nodes % y(i) - CoordNodes(j,2))**2
                IF(NoDims == 3) Dist = Dist + (Mesh % Nodes % z(i) - CoordNodes(j,3))**2
                Dist = SQRT(Dist)
                
                IF(Dist < MinDist(j) .AND. Dist <= Eps ) THEN
                  MinDist(j) = Dist
                  IndNodes(j) = i
                END IF
              END DO
            END DO

            ! In parallel case eliminate all except the nearest node. 
            ! This relies on the fact that for each node partition the 
            ! distance to nearest node is computed accurately. 
            DO j=1,NoNodes
              GlobalMinDist = ParallelReduction( MinDist(j), 1 )
              IF( ABS( GlobalMinDist - MinDist(j) ) > TINY(Dist) ) THEN
                IndNodes(j) = 0
              END IF
            END DO

            NOFNodesFound = 0
            DO j=1,NoNodes
               IF ( IndNodes(j)>0 ) THEN
                 NOFNodesFound = NOFNodesFound+1
                 IndNodes(NOFNodesFound) = IndNodes(j)
               END IF
            END DO
            
            ! In the first time add the found nodes to the list structure
            IF ( NOFNodesFound > 0 ) THEN
              CALL ListAddIntegerArray( ValueList,'Target Nodes', &
                  NOFNodesFound, IndNodes) 
              NodesFound = .TRUE.            
            ELSE
              ! If no nodes found, add still an empty list and make sure the 
              ! zero is not treated later on. Otherwise this search would be 
              ! retreated each time. 
              CALL ListAddIntegerArray( ValueList,'Target Nodes', &
                  1, IndNodes) 
            END IF

            ! Finally deallocate the temporal vectors
            DEALLOCATE( IndNodes, MinDist ) 
          END IF
        END IF
      END IF
      
      ! If the target coordinates has already been assigned to an empty list 
      ! cycle over it by testing the 1st node. 
      IF( NodesFound ) THEN
        NodeIndexes => ListGetIntegerArray( ValueList,'Target Nodes')
        IF( NodeIndexes(1) == 0 ) NodesFound = .FALSE. 
      END IF

      IF(NodesFound) THEN           
        Conditional = ListCheckPresent( ValueList, CondName )      
        n = SIZE(NodeIndexes)
        CALL SetPointValues(n)
      END IF
    END DO

    
!------------------------------------------------------------------------------
!   Go through soft upper and lower limits
!------------------------------------------------------------------------------
    Params => Model % Solver % Values
    ApplyLimiter = ListGetLogical( Params,'Apply Limiter',GotIt) 

    IF( Dof/=0 .AND. ApplyLimiter ) THEN
      CALL Info(Caller,'Applying limiters',Level=10)

      DO Upper=0,1
        
        ! The limiters have been implemented only componentwise
        !-------------------------------------------------------
        
        NULLIFY( LimitActive ) 
        Var => Model % Solver % Variable
        IF( Upper == 0 ) THEN
          IF( ASSOCIATED( Var % LowerLimitActive ) ) &
              LimitActive => Var % LowerLimitActive
        ELSE
          IF( ASSOCIATED( Var % UpperLimitActive ) ) &
              LimitActive => Var % UpperLimitActive
        END IF
        
        IF( .NOT. ASSOCIATED( LimitActive ) ) CYCLE
        
        IF( Upper == 0 ) THEN
          CondName = TRIM(name)//' Lower Limit' 
        ELSE
          CondName = TRIM(name)//' Upper Limit' 
        END IF
        
        ! check and set some flags for nodes belonging to n-t boundaries
        ! getting set by other bcs:
        ! --------------------------------------------------------------
        DO t = 1, Model % NumberOfBulkElements+Model % NumberOfBoundaryElements
          Element => Model % Elements(t)
          Model % CurrentElement => Element
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          
          IF( t > Model % NumberOfBulkElements ) THEN
            DO bc = 1,Model % NumberOfBCs
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
              ValueList => Model % BCs(BC) % Values
              CALL SetLimiterValues(n)
            END DO
          ELSE             
            bf_id = ListGetInteger( Model % Bodies(Element % bodyid) % Values, &
                'Body Force', GotIt)
            IF(.NOT. GotIt ) CYCLE
            ValueList => Model % Bodyforces(bf_id) % Values
            CALL SetLimiterValues(n)
          END IF
        END DO
      END DO
    END IF
    
    
    ! Check the boundaries and body forces for possible single nodes BCs that are used to fixed 
    ! the domain for undetermined equations. The loop is slower than optimal in the case that there is 
    ! a large amount of different boundaries that have a node to set. 
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Single Node'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

    IF( AnySingleBC .OR. AnySingleBF ) THEN
      Solver => Model % Solver
      Mesh   => Solver % Mesh

      DO bc = 1,Model % NumberOfBCs  + Model % NumberOfBodyForces    

        ! Make a distinction between BCs and BFs. 
        ! These are treated in the same loop because most of the logic is still the same. 
        IF( bc <= Model % NumberOfBCs ) THEN
          IF(.NOT. AnySingleBC ) CYCLE
          ValueList => Model % BCs(BC) % Values
          ElemFirst =  Mesh % NumberOfBulkElements + 1 
          ElemLast =  Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        ELSE
          IF( .NOT. AnySingleBF ) CYCLE
          ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
          ElemFirst =  1
          ElemLast =  Mesh % NumberOfBulkElements
        END IF

        SingleVal = ListGetCReal( ValueList,DirName, GotIt) 
        IF( .NOT. GotIt ) CYCLE
        ind = ListGetInteger( ValueList,TRIM(Name)//' Single Node Index',GotIt )     
        
        ! On the first time find a one single uniquely defined node for setting 
        ! the value. In parallel it will be an unshared node with the highest possible 
        ! node number 
        IF( .NOT. GotIt ) THEN        
          
          ind = 0
          DO t = ElemFirst, ElemLast
            Element => Mesh % Elements(t)
            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              j = Element % BodyId
              IF( j < 0 .OR. j > Model % NumberOfBodies ) CYCLE
              bf = ListGetInteger( Model % Bodies(j) % Values,'Body Force',GotIt)
              IF(.NOT. GotIt) CYCLE
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF
            
            DO i=1,n
              j = NodeIndexes(i)
              IF( Perm(j) == 0) CYCLE
              IF( ParEnv % PEs > 1 ) THEN
                IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
              END IF
              ind = j 
              EXIT
            END DO
            IF( ind > 0 ) EXIT
          END DO

          k = NINT( ParallelReduction( 1.0_dp * ind, 2 ) ) 
           
          ! Find the maximum partition that owns a suitable node. 
          ! It could be minimum also, just some convection is needed. 
          IF( ParEnv % PEs > 1 ) THEN
            k = -1
            IF( ind > 0 ) k = ParEnv % MyPe          
            k = NINT( ParallelReduction( 1.0_dp * k, 2 ) ) 
            IF( k == -1 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            ELSE
              IF( k /= ParEnv % MyPe ) ind = 0                         
              IF( InfoActive(8) ) THEN
                ind = NINT( ParallelReduction( 1.0_dp * ind, 2 ) )               
                CALL Info(Caller,'Fixing single node '&
                    //TRIM(I2S(ind))//' at partition '//TRIM(I2S(k)),Level=8)
                IF( k /= ParEnv % MyPe ) ind = 0
              END IF
            END IF
          ELSE
            IF( ind == 0 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            ELSE              
              CALL Info(Caller,'Fixing single node '//TRIM(I2S(ind)),Level=8)
            END IF
          END IF
            
          CALL ListAddInteger( ValueList,TRIM(Name)//' Single Node Index', ind )          
        END IF

        ! Ok, if this is the partition where the single node to eliminate the floating should 
        ! be eliminated then set it here. Index equal to zero tells that we are in a wrong partition.        
        IF( ind > 0 ) THEN
          CALL SetSinglePoint(ind,DOF,SingleVal,.TRUE.)
        END IF
      END DO
    END IF

    
!------------------------------------------------------------------------------
!   Take care of the matrix entries of passive elements
!------------------------------------------------------------------------------

    IF ( Passive ) THEN
      Solver => Model % Solver
      Mesh => Solver % Mesh
      DO i=1,Solver % NumberOfActiveElements
        Element => Mesh % Elements(Solver % ActiveElements(i))
        IF (CheckPassiveElement(Element)) THEN
          n = sGetElementDOFs(Indexes,UElement=Element)
          DO j=1,n
            k=Indexes(j)
            IF (k<=0) CYCLE

            k=Perm(k)
            IF (k<=0) CYCLE

            s=0._dp
            DO l=1,NDOFs
              m=NDOFs*(k-1)+l
              s=s+ABS(A % Values(A % Diag(m)))
            END DO
            IF (s>EPSILON(s)) CYCLE
 
            DO l=1,NDOFs
              m = NDOFs*(k-1)+l
              IF(A % ConstrainedDOF(m)) CYCLE
              CALL SetSinglePoint(k,l,Solver % Variable % Values(m),.FALSE.)
            END DO
          END DO
        END IF
      END DO
    END IF


    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Constant'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    AnySingleBF = ListCheckPresentAnyBodyForce( Model, DirName )

    IF( AnySingleBC .OR. AnySingleBF ) THEN
      ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

      IF( AnySingleBC ) CALL Info(Caller,'Found BC constraint: '//TRIM(DirName))
      IF( AnySingleBF ) CALL Info(Caller,'Found BodyForce constraint: '//TRIM(DirName))

      ! Improve the logic in future
      ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
      IF( AnySingleBC ) THEN
        NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Constant Node Index')
      ELSE 
        NeedListMatrix = .NOT. ListCheckPresentAnyBodyForce( Model, TRIM(Name)//' Constant Node Index')
      END IF
      
      ! Move the list matrix because of its flexibility
      IF( NeedListMatrix ) THEN
        CALL Info(Caller,'Using List maxtrix to set constant constraints',Level=8)
        CALL Info('SetDircihletBoundaries','Original matrix non-zeros: '&
            //TRIM(I2S(SIZE( A % Cols ))),Level=8)
        IF( ASSOCIATED( A % BulkValues ) ) THEN
          ALLOCATE( Cols0( SIZE( A % Cols ) ), Rows0( SIZE( A % Rows ) ) )
          Cols0 = A % Cols
          Rows0 = A % Rows
        END IF
        CALL List_toListMatrix(A)
      END IF

      DO bc = 1,Model % NumberOfBCs + Model % NumberOfBodyForces

        ! Make a distinction between BCs and BFs. 
        ! These are treated in the same loop because most of the logic is still the same. 
        IF( bc <= Model % NumberOfBCs ) THEN
          IF(.NOT. AnySingleBC ) CYCLE
          ValueList => Model % BCs(BC) % Values
          ElemFirst =  Model % NumberOfBulkElements + 1 
          ElemLast =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements
        ELSE
          IF(.NOT. AnySingleBF ) CYCLE
          ValueList => Model % BodyForces(bc - Model % NumberOfBCs) % Values
          ElemFirst =  1
          ElemLast =  Model % NumberOfBulkElements
        END IF

        IF( .NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE

        ind = ListGetInteger( ValueList,TRIM(Name)//' Constant Node Index',GotIt )     
        

        ! On the first time find a one single uniquely defined node for setting 
        ! the value. In parallel it will be an unshared node with the highest possible 
        ! node number 
        IF( .NOT. GotIt ) THEN        

          ind = 0
          DO t = ElemFirst, ElemLast
            Element => Model % Elements(t)

            IF( bc <= Model % NumberOfBCs ) THEN
              IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
            ELSE
              bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
              IF( bc - Model % NumberOfBCs /= bf ) CYCLE
            END IF            

            n = Element % TYPE % NumberOfNodes
            NodeIndexes => Element % NodeIndexes

            DO i=1,n
              j = NodeIndexes(i)
              IF( Perm(j) == 0) CYCLE
              IF( ParEnv % PEs > 1 ) THEN
                IF( SIZE( Mesh % ParallelInfo % NeighbourList(j) % Neighbours) > 1 ) CYCLE               
                IF( Mesh % ParallelInfo % NeighbourList(j) % Neighbours(1) /= ParEnv % MyPe ) CYCLE               
               END IF
              ind = j 
              EXIT
            END DO
            IF( ind > 0 ) EXIT
          END DO

          ! Find the maximum partition that owns the node. 
          ! It could be minimum also, just some convection is needed. 
          IF( ParEnv % PEs > 1 ) THEN
            k = -1
            IF( ind > 0 ) k = ParEnv % MyPe          
            k = NINT( ParallelReduction( 1.0_dp * k, 2 ) ) 
            IF( k == -1 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            END IF
            IF( k /= ParEnv % MyPe ) ind = 0
          ELSE
            IF( ind == 0 ) THEN
              CALL Warn(Caller,'Could not find node to set: '//TRIM(DirName))
            END IF
          END IF

          CALL ListAddInteger( ValueList,TRIM(Name)//' Constant Node Index', ind )
          NeedListMatrix = .TRUE.
        END IF

        IF( ParEnv % PEs > 1 ) CALL Warn(Caller,'Node index not set properly in parallel')
        IF( ind == 0 ) CYCLE

        ! Ok, now sum up the rows to the corresponding nodal index
        LumpedNodeSet = .FALSE.

        ! Don't lump the "supernode" and therefore mark it set already
        LumpedNodeSet(ind) = .TRUE.

        DO t = ElemFirst, ElemLast
          Element => Model % Elements(t)

          IF( bc <= Model % NumberOfBCs ) THEN
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
          ELSE
            bf = ListGetInteger( Model % Bodies(Element % bodyid) % Values,'Body Force',GotIt)
            IF( bc - Model % NumberOfBCs /= bf ) CYCLE
          END IF

          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes

          CALL SetLumpedRows(ind,n)
        END DO

        n = COUNT( LumpedNodeSet ) 
        CALL Info(Caller,'Number of lumped nodes set: '//TRIM(I2S(n)),Level=10)
      END DO

      IF( NeedListMatrix ) THEN
        DEALLOCATE( LumpedNodeSet )

        ! Revert back to CRS matrix
        CALL List_ToCRSMatrix(A)

        ! This is needed in order to copy the old BulkValues to a vector that 
        ! has the same size as the new matrix. Otherwise the matrix vector multiplication
        ! with the new Rows and Cols will fail. 
        IF( ASSOCIATED( A % BulkValues ) ) THEN
          BulkValues0 => A % BulkValues
          NULLIFY( A % BulkValues ) 
          ALLOCATE( A % BulkValues( SIZE( A % Values ) ) )
          A % BulkValues = 0.0_dp

          DO i=1,A % NumberOfRows
            DO j = Rows0(i), Rows0(i+1)-1
              k = Cols0(j) 
              DO j2 = A % Rows(i), A % Rows(i+1)-1
                k2 = A % Cols(j2)
                IF( k == k2 ) THEN
                  A % BulkValues(j2) = BulkValues0(j)
                  EXIT
                END IF
              END DO
            END DO
          END DO

          DEALLOCATE( Cols0, Rows0, BulkValues0 ) 
        END IF
        
        CALL Info('SetDircihletBoundaries','Modified matrix non-zeros: '&
            //TRIM(I2S(SIZE( A % Cols ))),Level=8)
      END IF
    END IF



    ! Check the boundaries and body forces for possible single nodes BCs that must have a constant
    ! value on that boundary / body force.
    !--------------------------------------------------------------------------------------------
    DirName = TRIM(Name)//' Plane'
    AnySingleBC = ListCheckPresentAnyBC( Model, DirName )
    
    IF( AnySingleBC ) THEN
      dim = CoordinateSystemDimension()
      
      ALLOCATE( LumpedNodeSet( SIZE( Perm ) ) )

      CALL Info(Caller,'Found BC constraint: '//TRIM(DirName))

      ! Improve the logic in future
      ! Now we assume that if the "supernode" has been found then also the matrix has the correct topology. 
      NeedListMatrix = .NOT. ListCheckPresentAnyBC( Model, TRIM(Name)//' Plane Node Indices')
      
      ! Move the list matrix because of its flexibility
      IF( NeedListMatrix ) THEN
        CALL Info(Caller,'Using List maxtrix to set constant constraints',Level=8)
        CALL Info('SetDircihletBoundaries','Original matrix non-zeros: '&
            //TRIM(I2S(SIZE( A % Cols ))),Level=8)
        IF( ASSOCIATED( A % BulkValues ) ) THEN
          ALLOCATE( Cols0( SIZE( A % Cols ) ), Rows0( SIZE( A % Rows ) ) )
          Cols0 = A % Cols
          Rows0 = A % Rows
        END IF
        CALL List_toListMatrix(A)
      END IF

      ElemFirst =  Model % NumberOfBulkElements + 1 
      ElemLast =  Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

      DO bc = 1,Model % NumberOfBCs 

        ValueList => Model % BCs(BC) % Values
        IF( .NOT. ListGetLogical( ValueList,DirName, GotIt) ) CYCLE

        PlaneInds => ListGetIntegerArray( ValueList,TRIM(Name)//' Plane Node Indices',GotIt )     

        IF(.NOT. GotIt ) THEN
          IF(.NOT. ALLOCATED(CandNodes) ) THEN
            ALLOCATE( CandNodes( Mesh % NumberOfNodes ) )        
          END IF
          CandNodes = .FALSE.

          ! Add nodes to the set that are associated with this BC only.
          DO t = ElemFirst, ElemLast
            Element => Model % Elements(t)            
            IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
              NodeIndexes => Element % NodeIndexes
              CandNodes(NodeIndexes) = .TRUE.
            END IF
          END DO

          ! Remove nodes from the set that may be set by other BCs also. 
          DO t = ElemFirst, ElemLast
            Element => Model % Elements(t)            
            IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) THEN
              NodeIndexes => Element % NodeIndexes
              CandNodes(NodeIndexes) = .FALSE.
            END IF
          END DO

          ALLOCATE(PlaneInds(3))
          CALL FindExtremumNodes(Mesh,CandNodes,dim,PlaneInds) 
          
          CALL ListAddIntegerArray( ValueList,TRIM(Name)//' Plane Node Indices',dim, PlaneInds )
          NeedListMatrix = .TRUE.
        END IF

        IF( ParEnv % PEs > 1 ) CALL Warn(Caller,'Node index perhaps not set properly in parallel')
        ! IF( ind == 0 ) CYCLE

        ! Ok, now sum up the rows to the corresponding nodal index
        LumpedNodeSet = .FALSE.

        ! Don't lump the "supernodes" and therefore mark it set already
        LumpedNodeSet(PlaneInds) = .TRUE.

        DO t = ElemFirst, ElemLast
          Element => Model % Elements(t)

          IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
            CALL SetRigidRows(PlaneInds,bc,n)
          END IF
        END DO

        n = COUNT( LumpedNodeSet ) 
        CALL Info(Caller,'Number of lumped nodes set: '//TRIM(I2S(n)),Level=10)
      END DO

      IF( NeedListMatrix ) THEN
        DEALLOCATE( LumpedNodeSet )

        ! Revert back to CRS matrix
        CALL List_ToCRSMatrix(A)

        ! This is needed in order to copy the old BulkValues to a vector that 
        ! has the same size as the new matrix. Otherwise the matrix vector multiplication
        ! with the new Rows and Cols will fail. 
        IF( ASSOCIATED( A % BulkValues ) ) THEN
          BulkValues0 => A % BulkValues
          NULLIFY( A % BulkValues ) 
          ALLOCATE( A % BulkValues( SIZE( A % Values ) ) )
          A % BulkValues = 0.0_dp

          DO i=1,A % NumberOfRows
            DO j = Rows0(i), Rows0(i+1)-1
              k = Cols0(j) 
              DO j2 = A % Rows(i), A % Rows(i+1)-1
                k2 = A % Cols(j2)
                IF( k == k2 ) THEN
                  A % BulkValues(j2) = BulkValues0(j)
                  EXIT
                END IF
              END DO
            END DO
          END DO

          DEALLOCATE( Cols0, Rows0, BulkValues0 ) 
        END IF
        
        CALL Info('SetDircihletBoundaries','Modified matrix non-zeros: '&
            //TRIM(I2S(SIZE( A % Cols ))),Level=8)
      END IF
    END IF

    
    IF( InfoActive(12) )  THEN
      DirCount = NINT( ParallelReduction( 1.0_dp * DirCount ) )
      CALL Info(Caller,'Number of dofs set for '//TRIM(Name)//': '&
          //TRIM(I2S(DirCount)),Level=12)
    END IF
      
    
!------------------------------------------------------------------------------

  CONTAINS

     ! Check n-t node setting element
     !-------------------------------
    SUBROUTINE CheckNTElement(n,elno)
      INTEGER :: n,elno
      INTEGER :: i,j,k,l,m,dim,kmax
      LOGICAL :: found
      REAL(KIND=dp) :: Condition(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF <= 0 ) RETURN
      IF ( ALL(BoundaryReorder(Indexes(1:n))<1) ) RETURN
      IF ( .NOT. ListCheckPresent(ValueList, Name) ) RETURN
      IF ( ListGetLogical(ValueList,NormalTangentialName,Found) ) RETURN

      IF ( Conditional ) THEN
        Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
        Conditional = Conditional .AND. GotIt
      END IF

      !
      ! Check for nodes belonging to n-t boundary getting set by other bcs.
      ! -------------------------------------------------------------------
      DO j=1,n
        IF ( Conditional .AND. Condition(j)<0.0_dp ) CYCLE
        k = Perm(Indexes(j))
        IF ( k > 0 ) THEN          
          k = k + OffSet
          m = BoundaryReorder(Indexes(j))
          IF ( m>0 ) THEN
            RotVec = 0._dp
            RotVec(DOF) = 1._dp
            CALL RotateNTSystem( RotVec, Indexes(j) )
            kmax = 1
            DO k=1,dim
              IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) kmax = k
            END DO
            NTelement(m,kmax)=elno
          END IF
        END IF
      END DO
!------------------------------------------------------------------------------
    END SUBROUTINE CheckNTElement
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to a specific boundary or bulk element.
!------------------------------------------------------------------------------
    SUBROUTINE SetElementValues(n,elno)
!------------------------------------------------------------------------------
      INTEGER :: n,elno
      INTEGER :: i,j,k,l,m,dim,kmax,lmax
      LOGICAL :: CheckNT,found
      REAL(KIND=dp) :: Condition(n), Work(n), RotVec(3)
      
      dim = CoordinateSystemDimension()

      IF ( DOF > 0 ) THEN
        IF (Model % Solver % DG) THEN
          Work(1:n)  = ListGetReal( ValueList, Name, n, Element % NodeIndexes, gotIt )
        ELSE
          Work(1:n)  = ListGetReal( ValueList, Name, n, Indexes, gotIt )
        END IF
        IF ( .NOT. GotIt ) THEN
          Work(1:n)  = ListGetReal( ValueList, Name(1:nlen) // ' DOFs', n, Indexes, gotIt )
        END IF
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, Indexes, gotIt )
      END IF
      
      IF ( gotIt ) THEN
        IF ( Conditional ) THEN
          IF (Model % Solver % DG) THEN
            Condition(1:n) = ListGetReal( ValueList, CondName, n, Element % NodeIndexes, gotIt )
          ELSE
            Condition(1:n) = ListGetReal( ValueList, CondName, n, Indexes, gotIt )
          END IF
          Conditional = Conditional .AND. GotIt
        END IF

       !
       ! Check for nodes belonging to n-t boundary getting set by other bcs.
       ! -------------------------------------------------------------------
        CheckNT = .FALSE.
        IF ( NormalTangentialNOFNodes>0 .AND. DOF>0 ) THEN
          CheckNT = .TRUE.
          IF ( ALL(BoundaryReorder(Indexes(1:n))<1) ) CheckNT = .FALSE.
          IF ( ListGetLogical(ValueList,NormalTangentialName,Found)) CheckNT=.FALSE.
        END IF
        
        DO j=1,n
          IF ( Conditional .AND. Condition(j) < 0.0d0 ) CYCLE

          k = Perm(Indexes(j))
          IF ( k > 0 ) THEN
            
            IF ( DOF>0 ) THEN
              m = 0
              IF ( NormalTangentialNOFNodes>0 ) m=BoundaryReorder(Indexes(j))
              IF ( m>0 .AND. CheckNT ) THEN
                RotVec = 0._dp
                RotVec(DOF) = 1._dp
                CALL RotateNTSystem( RotVec, Indexes(j) )

                ! When cartesian component "DOF" is defined set the N-T component
                ! closest to its direction. 
                kmax = 1 
                DO k=2,dim
                  IF ( ABS(RotVec(k)) > ABS(RotVec(kmax)) ) THEN
                    kmax = k
                  END IF
                END DO

                lmax = NDOFs * (Perm(Indexes(j))-1) + kmax
                IF ( .NOT. NTZeroing_done(m,kmax) ) THEN
                  NTZeroing_done(m,kmax) = .TRUE.
                  b(lmax) = 0._dp

                  IF( .NOT. OffDiagonal ) THEN
                    b(lmax) = b(lmax) + Work(j) !/DiagScaling(lmax)
                  END IF

                  ! Consider all components of the cartesian vector mapped to the 
                  ! N-T coordinate system. Should this perhaps have scaling included?
                  DirCount = DirCount + 1
                  CALL ZeroRow( A,lmax )
                  IF( .NOT. OffDiagonal) THEN
                    DO k=1,dim
                      l = NDOFs * (Perm(Indexes(j))-1) + k
                      CALL SetMatrixElement( A,lmax,l,RotVec(k) )
                    END DO
                  END IF
                  NTZeroing_done(m,kmax)   = .TRUE.
                  A % ConstrainedDOF(lmax) = .FALSE.
                END IF
              ELSE
                DirCount = DirCount + 1
                CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
              END IF
            ELSE
              DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                DirCount = DirCount + 1
                CALL SetSinglePoint(k,l,WorkA(l,1,j),.FALSE.)
              END DO
            END IF
          END IF
        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetElementValues
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Set values related to a specific boundary or bulk element.
!> If scaling has been applied the rows need to be scaled when
!> they are moved.
!------------------------------------------------------------------------------
    SUBROUTINE SetLumpedRows(ind0,n)
!------------------------------------------------------------------------------
      INTEGER :: ind0,n
      INTEGER :: ind,i,j,k,k0
      REAL(KIND=dp) :: Coeff
      ! -------------------------------------------------------------------        

      
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.

        IF ( DOF > 0 ) THEN
          k0 = Offset + NDOFs * (Perm(ind0)-1) + DOF
          k = OffSet + NDOFs * (Perm(ind)-1) + DOF

          Coeff = 1.0_dp
          
          CALL CRS_MoveRow( A, k, k0, Coeff )
          b(k0) = b(k0) + Coeff * b(k)

          CALL AddToMatrixElement( A, k, k, 1.0_dp )
          CALL AddToMatrixElement( A, k, k0, -Coeff )
          b(k) = 0.0_dp
        ELSE
          DO l = 1, NDOFs
            k0 = Offset + NDOFs + (Perm(ind0)-1) * DOF
            k = OffSet + NDOFs * (Perm(ind)-1) + l

            Coeff = 1.0_dp
            
            CALL CRS_MoveRow( A, k, k0, Coeff )
            b(k0) = b(k0) + Coeff * b(k)
          
            CALL AddToMatrixElement( A, k, k, 1.0_dp )
            CALL AddToMatrixElement( A, k, k0, -Coeff )
          END DO
        END IF
      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE SetLumpedRows
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to a rigid plane boundary such that any node on the boundary
!> is expressed as linear combinination of the selected three (or two if on line)
!> nodes.
!------------------------------------------------------------------------------
    SUBROUTINE SetRigidRows(inds0,bcind,n)
!------------------------------------------------------------------------------
      INTEGER :: inds0(:)
      INTEGER :: bcind
      INTEGER :: n

      INTEGER :: bcind0 =  0
      INTEGER :: ind,i,j,k,k0
      REAL(KIND=dp) :: Coeff, Weights(3)
      REAL(KIND=dp) :: BaseCoord(3,3),r1(3),r2(3),Coord(3),dCoord(3),Amat(2,2),A0mat(2,2),bvec(2)
      
      SAVE bcind0, BaseCoord, A0mat, r1, r2
!-------------------------------------------------------------------        

      IF(bcind /= bcind0 ) THEN
        BaseCoord = 0.0_dp
        DO i=1,dim
          j = inds0(i)
          BaseCoord(i,1) = Mesh % Nodes % x(j)
          BaseCoord(i,2) = Mesh % Nodes % y(j)
          BaseCoord(i,3) = Mesh % Nodes % z(j)
        END DO
        bcind0 = bcind

        r1 = BaseCoord(2,:) - BaseCoord(1,:)
        Amat(1,1) = SUM(r1*r1)
        
        IF( dim == 3 ) THEN
          r2 = BaseCoord(3,:) - BaseCoord(1,:)
          Amat(1,2) = SUM(r1*r2)
          Amat(2,1) = Amat(1,2)
          Amat(2,2) = SUM(r2*r2)
        END IF

        A0mat = Amat
        bcind0 = bcind
      END IF
                   
      DO j=1,n
        ind = Indexes(j)

        IF( LumpedNodeSet(ind) ) CYCLE
        LumpedNodeSet(ind) = .TRUE.
        
        Coord(1) = Mesh % Nodes % x(ind)
        Coord(2) = Mesh % Nodes % y(ind)
        Coord(3) = Mesh % Nodes % z(ind)

        dCoord = Coord - BaseCoord(1,:)
        
        bvec(1) = SUM( dCoord * r1 )
        IF( dim == 3 ) THEN
          bvec(2) = SUM( dCoord * r2 )
        END IF

        IF( dim == 2 ) THEN
          bvec(1) = bvec(1) / A0mat(1,1)
          Weights(2) = bvec(1)
          Weights(1) = 1.0_dp - Weights(2)
        ELSE
          Amat = A0mat          
          CALL LUSolve(2,Amat,bvec)          
          Weights(2:3) = bvec(1:2)
          Weights(1) = 1.0_dp - SUM(bvec(1:2))
        END IF

        DO l = 1, dim
          k = OffSet + NDOFs * (Perm(ind)-1) + l    

          ! Distribute row in accordance with the weights
          DO m = 1, dim
            k0 = Offset + NDOFs * (Perm(inds0(m))-1) + l          
            Coeff = Weights(m)

            b(k0) = b(k0) + Coeff * b(k)
            IF( m < dim ) THEN
              ! This does not nullify the row
              CALL CRS_MoveRow( A, k, k0, Coeff, 1.0_dp )
            ELSE
              ! Now also nullify the row
              CALL CRS_MoveRow( A, k, k0, Coeff )
              b(k) = 0.0_dp
            END IF
          END DO

          ! Express the node as linear combination of the base nodes
          DO m = 1,dim
            k0 = Offset + NDOFs * (Perm(inds0(m))-1) + l          
            Coeff = Weights(m)            
            CALL AddToMatrixElement( A, k, k0, -Coeff )
          END DO
          CALL AddToMatrixElement( A, k, k, 1.0_dp )
        END DO

      END DO

!------------------------------------------------------------------------------
    END SUBROUTINE SetRigidRows
!------------------------------------------------------------------------------

    
!------------------------------------------------------------------------------
!> Set values related to individual points.
!------------------------------------------------------------------------------
    SUBROUTINE SetPointValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: Work(n), Condition(n)        

      INTEGER :: i,j,k,k1,l

      IF ( DOF > 0 ) THEN
        Work(1:n) = ListGetReal( ValueList, Name, n, NodeIndexes, gotIt )
      ELSE
        CALL ListGetRealArray( ValueList, Name, WorkA, n, NodeIndexes, gotIt )
      END IF

      IF ( gotIt ) THEN

        Condition(1:n) = 1.0d0
        IF ( Conditional ) THEN
          Condition(1:n) = ListGetReal( ValueList, CondName, n, NodeIndexes, gotIt )
          Conditional = Conditional .AND. GotIt
        END IF

        DO j=1,n
          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF ( Conditional .AND. Condition(j) < 0.0d0 ) CYCLE
          IF ( NodeIndexes(j) > SIZE(Perm) .OR. NodeIndexes(j) < 1 ) THEN
            CALL Warn(Caller,'Invalid Node Number')
            CYCLE
          END IF

          IF ( DOF>0 ) THEN
            CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
          ELSE
            DO l=1,MIN( NDOFs, SIZE(Worka,1) )
              CALL SetSinglePoint(k,l,WorkA(l,1,j),.FALSE.)
            END DO
          END IF

        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetPointValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to one single point
!------------------------------------------------------------------------------
    SUBROUTINE SetSinglePoint(ind,DOF,val,ApplyPerm)
!------------------------------------------------------------------------------
      LOGICAL :: ApplyPerm
      INTEGER :: ind, DOF
      REAL(KIND=dp) :: val

      REAL(KIND=dp) :: s
      INTEGER :: i,j,k,k1,l


      IF(.NOT. ALLOCATED(A % ConstrainedDOF)) THEN
        ALLOCATE(A % ConstrainedDOF(A % NumberOfRows))
        A % ConstrainedDOF = .FALSE.
      END IF
      
      IF(.NOT. ALLOCATED(A % Dvalues)) THEN
        ALLOCATE(A % Dvalues(A % NumberOfRows))
        A % Dvalues = 0._dp
      END IF
      
      k = ind
      IF (ApplyPerm) k = Perm(ind)
      IF( k == 0 ) RETURN
      
      k = OffSet + NDOFs * (k-1) + DOF

      ! Do not add non-zero entries to pure halo nodes which are not associated with the partition.
      ! These are nodes could be created by the -halobc flag in ElmerGrid.
      IF( ParEnv % PEs > 1 ) THEN
        IF( .NOT. ANY( A % ParallelInfo % NeighbourList(k) % Neighbours == ParEnv % MyPe ) ) THEN
           RETURN
        END IF
      END IF

      DirCount = DirCount + 1
      
      IF( .NOT. OffDiagonal ) THEN
        A % Dvalues(k) =  val
      END IF
      A % ConstrainedDOF(k) = .TRUE.

!------------------------------------------------------------------------------
    END SUBROUTINE SetSinglePoint
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Set values related to upper and lower limiters.
!------------------------------------------------------------------------------
    SUBROUTINE SetLimiterValues(n)
!------------------------------------------------------------------------------
      INTEGER :: n
      REAL(KIND=dp) :: Work(n)

      Work(1:n)  = ListGetReal( ValueList, CondName, n, NodeIndexes, gotIt )

      IF ( gotIt ) THEN
        DO j=1,n
          k = Perm(NodeIndexes(j))
          IF( k == 0 ) CYCLE

          IF( .NOT. LimitActive(nDofs*(k-1)+dof)) CYCLE
          CALL SetSinglePoint(k,DOF,Work(j),.FALSE.)
        END DO
      END IF
!------------------------------------------------------------------------------
    END SUBROUTINE SetLimiterValues
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> At first pass sum together the rows related to the periodic dofs.
!------------------------------------------------------------------------------
   SUBROUTINE SetPeriodicBoundariesPass1( Model, A, b, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    LOGICAL :: Done(:)            !< Has the node already been done. 
    INTEGER :: This               !< Number of the current boundary.
    INTEGER :: DOF                !< The order number of the dof
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: p,q,i,j,k,l,m,n,nn,ii,nlen,jmp,size0
    INTEGER, POINTER :: PPerm(:)
    LOGICAL :: GotIt, Found, Jump
    REAL(KIND=dp) :: Scale, weight, coeff
    TYPE(Matrix_t), POINTER :: F, G, Projector, Projector1
    TYPE(Variable_t), POINTER :: Var, WeightVar
    TYPE(ValueList_t), POINTER :: BC
    TYPE(MortarBC_t), POINTER :: MortarBC
!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    BC => Model % BCs(This) % Values

    IF ( ListGetLogical( BC,& 
        'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      IF( ListGetLogical( BC,'Antisymmetric BC',GotIt ) ) THEN
        Scale = 1.0_dp
      ELSE
        Scale = -1.0_dp
      END IF
    ELSE IF ( ListGetLogical( BC, &
        'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = 1.0d0
    ELSE 
      Scale = ListGetConstReal( BC, &
          'Periodic BC Scale ' // Name(1:nlen), GotIt) 
      IF(.NOT. GotIt ) RETURN      
    END IF
    
    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN
    
!   For explicit conditions just create the dependency almost like a normal Dirichlet BC, 
!   For implicit one (otherwise) do the assembly of the projector:
!   ---------------------------------
    IF ( ListGetLogical( BC, &
        'Periodic BC Explicit', Found ) ) THEN
      
      Var => VariableGet( Model % Variables,Name(1:nlen) ) 
      
      DO i=1,Projector % NumberOfRows
        ii = Projector % InvPerm(i)
        IF( ii == 0 ) CYCLE
        k = Perm(ii)
        IF ( .NOT. Done(ii) .AND. k>0 ) THEN
          k = NDOFs * (k-1) + DOF
          !IF( .NOT.A % NoDirichlet ) THEN
          !  CALL ZeroRow( A,k )
          !  CALL AddToMatrixElement( A, k, k, 1.0_dp )
          !ELSE
          !END IF
          !IF(ALLOCATED(A % Dvalues))
          A % Dvalues(k) = 0._dp
          !IF(ALLOCATED(A % ConstrainedDOF))
          A % ConstrainedDOF(k) = .TRUE.
          
          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
            IF ( Projector % Cols(l) <= 0 ) CYCLE
            m = Perm( Projector % Cols(l) )
            IF ( m > 0 ) THEN
              m = NDOFs * (m-1) + DOF
              !IF(ALLOCATED(A % Dvalues)) THEN
                A % Dvalues(k) = A % Dvalues(k) - Scale * Projector % Values(l) * &
                    Var % Values(m) !/DiagScaling(k)
              !ELSE
              !  b(k) = b(k) - Scale * Projector % Values(l) * &
              !      Var % Values(m)/DiagScaling(k)
              !END IF
            END IF
          END DO
        END IF
      END DO
      
    ELSE IF ( ListGetLogical( BC, &
        'Periodic BC Use Lagrange Coefficient', Found ) ) THEN

      Jump = ListCheckPresent( BC, &
          'Periodic BC Coefficient '//Name(1:nlen))
      
      IF( .NOT. ASSOCIATED( Model % Solver % MortarBCs ) ) THEN
        CALL Info('SetPeriodicBoundariesPass1',&
            'Allocating mortar BCs for solver',Level=5)
        ALLOCATE( Model % Solver % MortarBCs( Model % NumberOfBCs ) )
        DO i=1, Model % NumberOfBCs
          Model % Solver % MortarBCs(i) % Projector => NULL()
        END DO
      END IF
      
      IF( ASSOCIATED( Projector, &
          Model % Solver % MortarBCs(This) % Projector) ) THEN
        CALL Info('SetPeridociBoundariesPass1','Using existing projector: '&
            //TRIM(I2S(This)),Level=8)
        RETURN
      END IF
      
      Model % Solver % MortarBCs(This) % Projector => Projector
      CALL Info('SetPeridociBoundariesPass1','Using projector as mortar constraint: '&
          //TRIM(I2S(This)),Level=5)

      MortarBC => Model % Solver % MortarBCs(This)      
      IF( Jump ) THEN
        IF( ASSOCIATED( MortarBC % Diag ) ) THEN
          IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Diag ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
          CALL Info('SetWeightedPeridocBCsPass1','Allocating projector mortar diag',Level=10)
          ALLOCATE( MortarBC % Diag( NDofs * Projector % NumberOfRows ) )
          MortarBC % Diag = 0.0_dp
        ELSE
          MortarBC % Diag( DOF::NDOFs ) = 0.0_dp
        END IF

        IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
          IF( SIZE( MortarBC % Rhs ) < NDofs * Projector % NumberOfRows ) THEN
            DEALLOCATE( MortarBC % Rhs ) 
          END IF
        END IF
        IF( .NOT. ASSOCIATED( MortarBC % Rhs ) ) THEN
          CALL Info('SetWeightedProjectorPass1','Allocating projector mortar rhs',Level=10)
          ALLOCATE( MortarBC % Rhs( NDofs * Projector % NumberOfRows ) )
          MortarBC % Rhs = 0.0_dp
        ELSE
          MortarBC % Rhs( DOF::NDOFs ) = 0.0_dp
        END IF
      END IF

      ! Create the permutation that is later need in putting the diag and rhs to correct position
      IF( ASSOCIATED( MortarBC % Perm ) ) THEN
        IF( SIZE( MortarBC % Perm ) < SIZE( Perm ) ) THEN
          DEALLOCATE( MortarBC % Perm ) 
        END IF
      END IF
      IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
        CALL Info('SetWeightedProjectorPass1','Allocating projector mortar perm',Level=10)
        ALLOCATE( MortarBC % Perm( SIZE( Perm ) ) )
      END IF
      
      MortarBC % Perm = 0
      DO i=1,SIZE( Projector % InvPerm )
        j = Projector % InvPerm(i) 
        IF( j > 0 .AND. j <= SIZE( Perm ) ) THEN
          MortarBC % Perm( j ) = i
        END IF
      END DO
      
      ! We can use directly the nodal projector
      MortarBC % Projector => Projector
      MortarBC % SlaveScale = -Scale
      MortarBC % MasterScale = -1.0_dp
 
      IF( Jump ) THEN
        PPerm => Perm
        CALL CalculateNodalWeights(Model % Solver,.TRUE.,&
            PPerm,'Periodic BC Coefficient '//Name(1:nlen),WeightVar )
        IF(.NOT. ASSOCIATED( WeightVar ) ) THEN
          CALL Fatal('SetPeriodicBoundariesPass1','Nodal weights needed for setting jumps!')
        END IF
        
        DO i=1,Projector % NumberOfRows
          k = Projector % InvPerm(i)
          IF ( k<=0 ) CYCLE
          
          ! Add the diagonal unity projector (scaled)
          weight = WeightVar % Values( PPerm( k ) )
          coeff = ListGetRealAtNode( BC,'Periodic BC Coefficient '&
              //Name(1:nlen), k, Found )

          ! For Nodal projector the entry is 1/(weight*coeff)
          ! For Galerkin projector the is weight/coeff 
          IF( Found ) THEN
            MortarBC % Diag( NDOFS* (i-1) + DOF ) = 1.0_dp / ( weight * coeff ) 
          END IF
        END DO
      END IF

      Model % Solver % MortarBCsChanged = .TRUE.
      
    ELSE

      ALLOCATE(F)
      F = A
      F % RHS => F % BulkRHS
      F % Values => F % BulkValues

      DO i=1,Projector % NumberOfRows
         ii = Projector % InvPerm(i)
         IF( ii == 0 ) CYCLE
         k = Perm(ii)
         IF ( .NOT. Done(ii) .AND. k>0 ) THEN
            k = NDOFs * (k-1) + DOF
            DO l=Projector % Rows(i),Projector % Rows(i+1)-1
              IF ( Projector % Cols(l) <= 0 .OR. Projector % Values(l)==0.0d0 ) CYCLE

              m = Perm(Projector % Cols(l))
              IF ( m > 0 ) THEN
                m = NDOFs*(m-1) + DOF
                DO nn=A % Rows(k),A % Rows(k+1)-1
                   CALL AddToMatrixElement( A, m, A % Cols(nn), &
                          -scale*Projector % Values(l) * A % Values(nn) ) !/ &
                          !DiagScaling(k) * DiagScaling(m))
                   IF (ASSOCIATED(F % Values)) THEN
                     CALL AddToMatrixElement( F, m, F % Cols(nn), &
                          -scale*Projector % Values(l) * F % Values(nn) )
                   END IF
                END DO
                b(m)=b(m) - scale*Projector % Values(l)*b(k) !*&
                        !DiagScaling(m) / DiagScaling(k)
                IF (ASSOCIATED(F % RHS)) THEN
                  F % RHS(m) = F % RHS(m) - scale*Projector % Values(l)*F % RHS(k)
                END IF
              END IF
            END DO
         END IF
         Done(ii) = .TRUE.
      END DO
      DEALLOCATE(F)
    END IF

!------------------------------------------------------------------------------
   END SUBROUTINE SetPeriodicBoundariesPass1
!------------------------------------------------------------------------------


!> At second pass add the [...1 .. -1 ...] row structure that results to the 
!> equality of the periodic dofs. 
!------------------------------------------------------------------------------
   SUBROUTINE SetPeriodicBoundariesPass2( Model, A, b, &
                      Name, DOF, NDOFs, Perm, This, Done )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    LOGICAL :: Done(:)            !< Has the node already been done. 
    INTEGER :: This               !< Number of the current boundary.
    INTEGER :: DOF                !< The order number of the dof
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,nn,ii,nlen
    LOGICAL :: GotIt, Jump, Found
    REAL(KIND=dp) :: Scale,s,ValueOffset,val,coeff,weight
    TYPE(Matrix_t), POINTER :: Projector
    INTEGER, POINTER :: PPerm(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Variable_t), POINTER :: WeightVar
!------------------------------------------------------------------------------


    BC =>  Model % BCs(This) % Values
    IF ( ListGetLogical( BC, &
         'Periodic BC Use Lagrange Coefficient', GotIt ) ) RETURN

    IF ( ListGetLogical( BC, &
         'Periodic BC Explicit', GotIt ) ) RETURN

    nlen = LEN_TRIM(Name)

    IF ( ListGetLogical( BC, &
       'Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      IF( ListGetLogical( BC,'Antisymmetric BC',GotIt ) ) THEN
        Scale = 1.0_dp
      ELSE
        Scale = -1.0_dp
      END IF
    ELSE IF ( ListGetLogical( BC, &
        'Anti Periodic BC ' // Name(1:nlen), GotIt ) ) THEN
      Scale = 1.0d0
    ELSE 
      Scale = ListGetCReal( BC, &
          'Periodic BC Scale ' // Name(1:nlen), GotIt) 
      IF(.NOT. GotIt ) RETURN      
    END IF

    ValueOffset = ListGetCReal( BC, &
          'Periodic BC Offset ' // Name(1:nlen), GotIt) 

    Jump = ListCheckPresent( BC, &
        'Periodic BC Coefficient '//Name(1:nlen))
    IF( Jump ) THEN
      PPerm => Perm
      CALL CalculateNodalWeights(Model % Solver,.TRUE.,&
          PPerm,'Periodic BC Coefficient '//Name(1:nlen),WeightVar )
      IF(.NOT. ASSOCIATED( WeightVar ) ) THEN
        CALL Fatal('SetPeriodicBoundariesPass1','Nodal weights needed for setting jumps!')
      END IF
    END IF

    Projector => Model % BCs(This) % PMatrix
    IF ( .NOT. ASSOCIATED(Projector) ) RETURN


!   Do the assembly of the projector:
!   ---------------------------------
    DO i=1,Projector % NumberOfRows
       ii = Projector % InvPerm(i)
       IF( ii == 0 ) CYCLE

       k = Perm( ii )
       IF ( .NOT. Done(ii) .AND. k > 0 ) THEN

         IF( Jump ) THEN
           weight = WeightVar % Values( k )
           coeff = ListGetRealAtNode( BC,'Periodic BC Coefficient '&
               //Name(1:nlen),ii, Found )
           val = -weight * coeff !* DiagScaling(k)**2
           scale = -1.0
         ELSE         
           val = 1.0_dp
         END IF

          k = NDOFs * (k-1) + DOF
          IF(.NOT. Jump) THEN
            CALL ZeroRow( A,k )
            b(k) = 0.0_dp
          END IF

          DO l = Projector % Rows(i), Projector % Rows(i+1)-1
             IF ( Projector % Cols(l) <= 0 ) CYCLE
             m = Perm( Projector % Cols(l) )
             IF ( m > 0 ) THEN
               m = NDOFs * (m-1) + DOF
               CALL AddToMatrixElement( A,k,m,val * Projector % Values(l) ) !* &
                   !( DiagScaling(m) / DiagScaling(k) ) )
             END IF
          END DO

          !IF(.NOT. Jump) THEN
          !  A % ConstrainedDof(k) = .TRUE.
          !  A % DValues(k) = -ValueOffset
          !ELSE          
            b(k) = b(k) - ValueOffset !/ DiagScaling(k)
          !END IF
          CALL AddToMatrixElement( A,k,k,scale*val )
          
        END IF
       Done(ii) = .TRUE.
    END DO
!------------------------------------------------------------------------------
   END SUBROUTINE SetPeriodicBoundariesPass2
!------------------------------------------------------------------------------




!> At first pass sum together the rows related to the periodic dofs.
!------------------------------------------------------------------------------
   SUBROUTINE SetFrictionBoundaries( Model, A, b, &
                      Name, NDOFs, Perm )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: t,u,j,k,k2,l,l2,n,bc_id,nlen,NormalInd
    LOGICAL :: Found
    REAL(KIND=dp),ALLOCATABLE :: Coeff(:)
    LOGICAL, ALLOCATABLE :: NodeDone(:)
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
!------------------------------------------------------------------------------


    CALL Info('SetFrictionBoundaries','Setting friction boundaries for variable: '//TRIM(Name),&
        Level=8)

    IF( NDOFs /= 2 ) THEN
      CALL Warn('SetFrictionBoundaries','Assumes friction only in 2D system')
    END IF

    nlen = LEN_TRIM(Name)
    Mesh => Model % Mesh

    ALLOCATE( NodeDone( SIZE( Perm ) ) )
    ALLOCATE( Coeff( Mesh % MaxElementNodes ) )

    NodeDone = .FALSE.
    Coeff = 0.0_dp
    
    DO t = Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      
      Model % CurrentElement => Element
            
      DO bc_id = 1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
      END DO
      IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE     
      BC => Model % BCs(bc_id) % Values

      IF( .NOT. ListGetLogical( BC,& 
          'Friction BC ' // Name(1:nlen), Found ) ) CYCLE

      NodeIndexes => Element % NodeIndexes
      n = Element % TYPE % NumberOfNodes
      
      Coeff(1:n) = ListGetReal( BC,& 
          'Friction Coefficient ' // Name(1:nlen), n, NodeIndexes )
      IF( ListGetLogical( BC,& 
          'Normal-Tangential ' // Name(1:nlen), Found ) ) THEN
        NormalInd = 1 
      ELSE
        NormalInd = ListGetInteger( BC,& 
            'Friction Normal Component ' // Name(1:nlen) )
      END IF

      DO i = 1, n
        j = Perm( Nodeindexes(i) )
        IF( NodeDone( j ) ) CYCLE

        k = NDOFs * (j-1) + NormalInd
        k2 = NDOFs * (j-1) + ( 3 - NormalInd ) 

        DO l = A % Rows(k),A % Rows(k+1)-1
          DO l2 = A % Rows(k2), A % Rows(k2+1)-1
            IF( A % Cols(l2) == A % Cols(l) ) EXIT
          END DO
          A % Values(l2) = A % Values(l2) - Coeff(i) * A % Values(l)
        END DO
        A % Rhs(k2) = A % Rhs(k2) - Coeff(i) * A % Rhs(k)
        NodeDone( j ) = .TRUE.
      END DO
    END DO

    n = COUNT( NodeDone ) 
    CALL Info('SetFrictionBoundaries','Number of friction nodes: '//TRIM(I2S(n)),Level=10)

    DEALLOCATE( NodeDone, Coeff )

!------------------------------------------------------------------------------
  END SUBROUTINE SetFrictionBoundaries
!------------------------------------------------------------------------------


!> Set the diagonal entry related to mortar BCs.
!> This implements the implicit jump condition. 
!------------------------------------------------------------------------------
   SUBROUTINE SetWeightedProjectorJump( Model, A, b, &
       Name, DOF, NDOFs, Perm )
     !------------------------------------------------------------------------------
     TYPE(Model_t) :: Model        !< The current model structure
     TYPE(Matrix_t), POINTER :: A  !< The global matrix
     REAL(KIND=dp) :: b(:)         !< The global RHS vector
     CHARACTER(LEN=*) :: Name      !< name of the dof to be set
     INTEGER :: DOF                !< The order number of the dof
     INTEGER :: NDOFs              !< the total number of DOFs for this equation
     INTEGER, TARGET :: Perm(:)    !< The node reordering info
     !------------------------------------------------------------------------------
     INTEGER :: i,j,k,i2,j2,k2,n,u,v,node,totsize,nodesize,bc_ind,&
         nnodes,nlen,TargetBC
     INTEGER, POINTER :: PPerm(:)
     INTEGER, ALLOCATABLE :: InvInvPerm(:)
     LOGICAL :: Found, AddRhs, AddCoeff, AddRes
     LOGICAL, ALLOCATABLE :: NodeDone(:)
     REAL(KIND=dp) :: coeff, weight, voff, res
     TYPE(Matrix_t), POINTER :: Projector
     TYPE(ValueList_t), POINTER :: BC
     TYPE(Element_t), POINTER :: Element, Left, Right
     LOGICAL :: SomethingDone
     TYPE(MortarBC_t), POINTER :: MortarBC
     !------------------------------------------------------------------------------

     ! If there is no mortar projector then nothing to do
     SomethingDone = .FALSE.

     ! Go through the projectors and check for jumps
     ! If there is a jump add an entry to the diagonal-to-be
     DO bc_ind=1,Model % NumberOFBCs

       MortarBC => Model % Solver % MortarBCs(bc_ind) 

       Projector => MortarBC % Projector
       IF( .NOT. ASSOCIATED( Projector ) ) CYCLE

       ! For this boundary there should also be a coefficient 
       ! otherwise nothing needs to be done. 
       nlen = LEN_TRIM(Name)
       BC => Model % BCs(bc_ind) % Values

       AddCoeff = ListCheckPresent( BC,'Mortar BC Coefficient '//Name(1:nlen))
       AddRes = ListCheckPresent( BC,'Mortar BC Resistivity '//Name(1:nlen))
       AddRhs = ListCheckPresent( BC,'Mortar BC Offset '//Name(1:nlen))

       IF( .NOT. (AddCoeff .OR. AddRes .OR. AddRhs) ) CYCLE

       Model % Solver % MortarBCsChanged = .TRUE.
       
       IF( .NOT. ASSOCIATED( Projector % InvPerm ) ) THEN
         CALL Fatal('SetWeightedProjectorJump','The > Projector % InvPerm < is really needed here!')
       END IF

       totsize = MAXVAL( Perm )
       nodesize = MAXVAL( Perm(1:Model % Mesh % NumberOfNodes ) )

       IF( AddCoeff .OR. AddRes ) THEN
         IF( ASSOCIATED( MortarBC % Diag ) ) THEN
           IF( SIZE( MortarBC % Diag ) < NDofs * Projector % NumberOfRows ) THEN
             DEALLOCATE( MortarBC % Diag ) 
           END IF
         END IF
         IF( .NOT. ASSOCIATED( MortarBC % Diag ) ) THEN
           CALL Info('SetWeightedProjectorJump','Allocating projector mortar diag',Level=10)
           ALLOCATE( MortarBC % Diag( NDofs * Projector % NumberOfRows ) )
           MortarBC % Diag = 0.0_dp
         ELSE
           MortarBC % Diag(DOF::NDOFs) = 0.0_dp
         END IF
       END IF

       IF( AddRhs ) THEN
         IF( ASSOCIATED( MortarBC % Rhs ) ) THEN
           IF( SIZE( MortarBC % Rhs ) < NDofs * Projector % NumberOfRows ) THEN
             DEALLOCATE( MortarBC % Rhs ) 
           END IF
         END IF
         IF( .NOT. ASSOCIATED( MortarBC % Rhs ) ) THEN
           CALL Info('SetWeightedProjectorJump','Allocating projector mortar rhs',Level=10)
           ALLOCATE( MortarBC % Rhs( NDofs * Projector % NumberOfRows ) )
           MortarBC % Rhs = 0.0_dp
         ELSE
           MortarBC % Rhs(DOF::NDOFs) = 0.0_dp
         END IF
       END IF

       ! Create the permutation that is later need in putting the diag and rhs to correct position
       IF( ASSOCIATED( MortarBC % Perm ) ) THEN
         IF( SIZE( MortarBC % Perm ) < SIZE( Perm ) ) THEN
           DEALLOCATE( MortarBC % Perm ) 
         END IF
       END IF
       IF( .NOT. ASSOCIATED( MortarBC % Perm ) ) THEN
         CALL Info('SetWeightedProjectorJump','Allocating projector mortar perm',Level=10)
         ALLOCATE( MortarBC % Perm( SIZE( Perm ) ) )
       END IF

       MortarBC % Perm = 0
       DO i=1,SIZE( Projector % InvPerm )
         j = Projector % InvPerm(i) 
         IF( j > 0 .AND. j <= nodesize ) THEN
           MortarBC % Perm( j ) = i
         END IF
       END DO


       TargetBC = ListGetInteger( BC,'Mortar BC',Found ) 

       CALL Info('SetWeightedProjectorJump','Setting jump to mortar projector in BC '&
           //TRIM(I2S(bc_ind)),Level=7)
    
       ! Create a table that shows how the additional degrees of freedom map
       ! to their corresponding regular dof. This is needed when creating the jump.
       ALLOCATE( NodeDone( Projector % NumberOfRows ) )
       NodeDone = .FALSE.
       
       ! Looping through elements rather than looping through projector rows directly
       ! is done in order to be able to refer to boundary properties associated 
       ! with the element. 
       DO t=1,Model % Mesh % NumberOfBoundaryElements
         Element => Model % Mesh % Elements( t + Model % Mesh % NumberOfBulkElements )

         IF( Element % BoundaryInfo % Constraint /= Model % BCs(bc_ind) % Tag ) CYCLE

         ! Outside code this tells the active element
         Model % CurrentElement => Element

         Left => Element % BoundaryInfo % Left
         Right => Element % BoundaryInfo % Right 
        
         IF( TargetBC > 0 ) THEN
           IF( ASSOCIATED( Left ) ) THEN
             IF( Left % PartIndex /= ParEnv % myPE ) CYCLE
           ELSE IF ( ASSOCIATED( Right ) ) THEN
             IF( Left % PartIndex /= ParEnv % myPE ) CYCLE
           ELSE
             CYCLE
           END IF
         ELSE
           ! This case is for the case when TargetBC = 0 i.e. for Discontinuous BC
           ! These are conditions that resulted to creation of zero 
           ! constraint matrix entries in this partition so no need to do them.
           IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
             CYCLE
           END IF
           
           ! For this we have a zero mass matrix entry so don't bother to add zero
!          IF( Left % PartIndex /= ParEnv % myPE .AND. &
!              Right % PartIndex /= ParEnv % myPe ) THEN
!            CYCLE
!          END IF
         END IF

         nnodes = Element % TYPE % NumberOfNodes
         DO u=1, nnodes
           node = Element % NodeIndexes(u)

           IF( Perm( node ) == 0 ) CYCLE

           i = MortarBC % Perm( node ) 
           IF( i == 0 ) CYCLE

           IF( NodeDone( i ) ) CYCLE
           NodeDone( i ) = .TRUE. 

           Found = .FALSE.

           IF( AddCoeff ) THEN
             coeff = ListGetRealAtNode( BC,'Mortar BC Coefficient '&
                 //Name(1:nlen),node, Found )        
             res = 1.0_dp / coeff
           END IF

           IF( AddRes ) THEN
             res = ListGetRealAtNode( BC,'Mortar BC Resistivity '&
                 //Name(1:nlen),node, Found )
           END IF

           ! For Nodal projector the entry is 1/(weight*coeff)
           ! For Galerkin projector the is weight/coeff 
           IF( Found ) THEN 
             IF( AddCoeff .OR. Addres ) THEN
               MortarBC % Diag(NDOFs*(i-1)+DOF) = res
             END IF
           END IF

           IF( AddRhs ) THEN
             voff = ListGetRealAtNode( BC,'Mortar BC Offset '&
                 //Name(1:nlen),node, Found )        
             IF( Found ) THEN
               MortarBC % Rhs(NDofs*(i-1)+DOF) = voff
             END IF
           END IF

         END DO
       END DO
       
       SomethingDone = .TRUE.

       DEALLOCATE( NodeDone )
     END DO

     IF( SomethingDone ) THEN
       CALL Info('setWeightedProjectorJump','Created a jump for weighted projector',Level=7)
     END IF
 
!------------------------------------------------------------------------------
   END SUBROUTINE SetWeightedProjectorJump
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletBoundaries
!------------------------------------------------------------------------------

  
!> Prepare to set Dirichlet conditions for attachment DOFs in the case of
!> component mode synthesis
!------------------------------------------------------------------------------
  SUBROUTINE SetConstraintModesBoundaries( Model, A, b, &
      Name, NDOFs, Perm )
    !------------------------------------------------------------------------------
    TYPE(Model_t) :: Model        !< The current model structure
    TYPE(Matrix_t), POINTER :: A  !< The global matrix
    REAL(KIND=dp) :: b(:)         !< The global RHS vector
    CHARACTER(LEN=*) :: Name      !< name of the dof to be set
    INTEGER :: NDOFs              !< the total number of DOFs for this equation
    INTEGER, TARGET :: Perm(:)    !< The node reordering info, this has been generated at the
                                  !< beginning of the simulation for bandwidth optimization
!------------------------------------------------------------------------------
    INTEGER :: i,t,u,j,k,k2,l,l2,n,bc_id,nlen,NormalInd
    LOGICAL :: Found
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Element_t), POINTER :: Element
    TYPE(Variable_t), POINTER :: Var
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER, ALLOCATABLE :: BCPerm(:)

!------------------------------------------------------------------------------

    nlen = LEN_TRIM(Name)
    Mesh => Model % Mesh
    Solver => Model % Solver
    Var => Solver % Variable 
    
    ! This needs to be allocated only once, hence return if already set
    IF( Var % NumberOfConstraintModes > 0 ) RETURN

    CALL Info('SetConstraintModesBoundaries','Setting constraint modes boundaries for variable: '&
        //TRIM(Name),Level=7)

    ! Allocate the indeces for the constraint modes
    ALLOCATE( Var % ConstraintModesIndeces( A % NumberOfRows ) )
    Var % ConstraintModesIndeces = 0
    
    ALLOCATE( BCPerm( Model % NumberOfBCs ) ) 
    BCPerm = 0
    
    j = 0 
    DO bc_id = 1,Model % NumberOfBCs
      BC => Model % BCs(bc_id) % Values        
      k = ListGetInteger( BC,&
          'Constraint Mode '// Name(1:nlen), Found )
      IF( Found ) THEN
        IF( k == 0 ) k = -1  ! Ground gets negative value
        BCPerm(bc_id) = k        
      ELSE IF( ListGetLogical( BC,& 
          'Constraint Modes ' // Name(1:nlen), Found ) ) THEN
        j = j + 1
        BCPerm(bc_id) = j
      END IF
    END DO
    
    j = MAXVAL( BCPerm ) 
    CALL Info('SetConstraintModesBoundaries','Number of active constraint modes boundaries: '&
        //TRIM(I2S(j)),Level=7)
    IF( j == 0 ) THEN
      CALL Fatal('SetConstraintModesBoundaries',&
          'Constraint Modes Analysis requested but no constrained BCs given!')
    END IF

    Var % NumberOfConstraintModes = NDOFS * j 
    
    
    DO t = Mesh % NumberOfBulkElements+1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(t)
      
      DO bc_id = 1,Model % NumberOfBCs
        IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc_id) % Tag ) EXIT
      END DO
      IF( bc_id > Model % NumberOfBCs ) CYCLE      
      IF( BCPerm(bc_id) == 0 ) CYCLE
      
      NodeIndexes => Element % NodeIndexes

      ! For vector valued problems treat each component as separate dof
      DO k=1,NDOFs
        Var % ConstraintModesIndeces( NDOFS*(Perm(NodeIndexes)-1)+k) = NDOFS*(BCPerm(bc_id)-1)+k
      END DO
    END DO
    
    ! The constraint modes can be either lumped or not.
    ! If they are not lumped then mark each individually
    IF( .NOT. ListGetLogical(Solver % Values,'Constraint Modes Lumped',Found ) ) THEN
      j = 0
      DO i=1,A % NumberOfRows
        IF( Var % ConstraintModesIndeces(i) > 0 ) THEN
          j = j + 1
          Var % ConstraintModesIndeces(i) = j
        END IF
      END DO
      CALL Info('SetConstraintModesBoundaries','Number of active constraint modes: '&
          //TRIM(I2S(j)),Level=7)
      Var % NumberOfConstraintModes = j 
    END IF
    

    ! Manipulate the boundaries such that we need to modify only the r.h.s. in the actual linear solver
    DO k=1,A % NumberOfRows       
      IF( Var % ConstraintModesIndeces(k) == 0 ) CYCLE
      A % ConstrainedDOF(k) = .TRUE.
      A % DValues(k) = 0.0_dp
    END DO
    
    ALLOCATE( Var % ConstraintModes( Var % NumberOfConstraintModes, A % NumberOfRows ) )
    Var % ConstraintModes = 0.0_dp

    DEALLOCATE( BCPerm ) 
    
    CALL Info('SetConstraintModesBoundarues','All done',Level=10)

!------------------------------------------------------------------------------
  END SUBROUTINE SetConstraintModesBoundaries
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Sets just one Dirichlet point in contrast to setting the whole field.
!> This is a lower order routine that the previous one. 
!------------------------------------------------------------------------------
  SUBROUTINE SetDirichletPoint( A, b,DOF, NDOFs, Perm, NodeIndex, NodeValue) 
!------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: b(:)
    REAL(KIND=dp) :: NodeValue
    INTEGER :: DOF, NDOFs, Perm(:), NodeIndex
!------------------------------------------------------------------------------

    REAL(KIND=dp) :: s
    INTEGER :: PermIndex

!------------------------------------------------------------------------------
    PermIndex = Perm(NodeIndex)
    IF ( PermIndex > 0 ) THEN
      PermIndex = NDOFs * (PermIndex-1) + DOF
      A % ConstrainedDOF(PermIndex) = .TRUE.
      A % DValues(PermIndex) = NodeValue      
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE SetDirichletPoint
!------------------------------------------------------------------------------



!------------------------------------------------------------------------------
!> Sets nodal loads directly to the matrix structure. 
!> The intended use for this, is for example, in multiphysics coupling where
!> the nodal loads may have been computed by another solver. 
!------------------------------------------------------------------------------
   SUBROUTINE SetNodalLoads( Model, A, b, Name, DOF, NDOFs, Perm )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model         !< The current model structure
    TYPE(Matrix_t), POINTER :: A   !< The global matrix
    REAL(KIND=dp) :: b(:)          !< The global RHS vector
    CHARACTER(LEN=*) :: Name       !< Name of the dof to be set
    INTEGER :: DOF                 !< The order number of the dof
    INTEGER :: NDOFs               !< The total number of DOFs for this equation
    INTEGER :: Perm(:)             !< The node reordering info, this has been generated at the
                                   !< beginning of the simulation for bandwidth optimization.
!------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: Element
    INTEGER, ALLOCATABLE :: Indexes(:)
    INTEGER, POINTER :: NodeIndexes(:), Neigh(:)
    INTEGER :: BC,i,j,k,l,n,t,k1,k2
    LOGICAL :: GotIt
    REAL(KIND=dp), POINTER :: WorkA(:,:,:) => NULL()
    REAL(KIND=dp) ::  s

    LOGICAL :: Conditional
    CHARACTER(LEN=MAX_NAME_LEN) :: LoadName

    INTEGER, POINTER :: IndNodes(:)
    INTEGER :: NoNodes,NoDims,bf_id,nlen,NOFNodesFound
    REAL(KIND=dp), POINTER :: CoordNodes(:,:), DiagScaling(:),MinDist(:)
    REAL(KIND=dp) :: GlobalMinDist,Dist,Eps
    LOGICAL, ALLOCATABLE :: ActivePart(:), ActivePartAll(:), DoneLoad(:)
    LOGICAL :: NodesFound
    TYPE(ValueList_t), POINTER :: ValueList

    LoadName = TRIM(Name) // ' Load'
    nlen = LEN_TRIM(LoadName)
    
    CALL Info('SetNodalLoads','Checking for nodal loads for variable: '//TRIM(Name),Level=12)

    n = MAX(Model % NumberOfBCs, Model % NumberOFBodyForces) 
    ALLOCATE( ActivePart(n), ActivePartAll(n) )

    ALLOCATE( Indexes(Model % Solver % Mesh % MaxElementDOFs) )
!------------------------------------------------------------------------------
! Go through the boundaries
!------------------------------------------------------------------------------

    !DiagScaling => A % DiagScaling
    !IF (.NOT.ASSOCIATED(DiagScaling)) THEN
    !  ALLOCATE(DiagScaling(A % NumberOFRows))
    !  DiagScaling=1._dp
    !END IF

    ActivePart = .FALSE.
    ActivePartAll = .FALSE.
    DO BC=1,Model % NumberOfBCs
      IF(.NOT. ListCheckPresent( Model % BCs(BC) % Values,'Target Boundaries')) CYCLE
      ActivePart(BC) = ListCheckPresent( Model % BCs(BC) % Values, LoadName )
      ActivePartAll(BC) = ListCheckPresent( &
          Model % BCs(BC) % Values, LoadName(1:nlen) // ' DOFs' )
    END DO

    IF ( ANY(ActivePart) .OR. ANY(ActivePartAll) ) THEN
      CALL Info('SetNodalLoads','Setting nodal loads on boundaries: '//TRIM(LoadName),Level=9)
      ALLOCATE(DoneLoad( SIZE(b)/NDOFs) )
      DoneLoad = .FALSE.

      DO BC=1,Model % NumberOfBCs
        IF(.NOT. ActivePart(BC) .AND. .NOT. ActivePartAll(BC) ) CYCLE

        DO t = Model % NumberOfBulkElements + 1, &
          Model % NumberOfBulkElements + Model % NumberOfBoundaryElements

          Element => Model % Elements(t)
          IF ( Element % BoundaryInfo % Constraint /= Model % BCs(BC) % Tag ) CYCLE
          
          Model % CurrentElement => Element
          IF ( ActivePart(BC) ) THEN
            n = Element % TYPE % NumberOfNodes
            Indexes(1:n) = Element % NodeIndexes
          ELSE
            n = SgetElementDOFs( Indexes )
          END IF
          ValueList => Model % BCs(BC) % Values

          CALL SetElementLoads( n )
        END DO
      END DO
    END IF

!------------------------------------------------------------------------------
! Go through the nodal load conditions for the body force list
!------------------------------------------------------------------------------

    ActivePart = .FALSE.
    ActivePartAll = .FALSE.
    DO bf_id=1,Model % NumberOFBodyForces
      ActivePart(bf_id) = ListCheckPresent( Model % BodyForces(bf_id) % Values, LoadName ) 
      ActivePartAll(bf_id) = ListCheckPresent( &
            Model % BodyForces(bf_id) % Values, LoadName(1:nlen) // ' DOFs' ) 
    END DO

    IF ( ANY( ActivePart ) .OR. ANY(ActivePartAll) ) THEN
      CALL Info('SetNodalLoads','Setting nodal loads on body force: '//TRIM(LoadName),Level=9)
      IF(.NOT. ALLOCATED(DoneLoad)) ALLOCATE(DoneLoad( SIZE(b)/NDOFs) )      
      DoneLoad = .FALSE.

      DO t = 1, Model % NumberOfBulkElements 
        Element => Model % Elements(t)
        bf_id = ListGetInteger( Model % Bodies(Element % BodyId) % Values,'Body Force', GotIt)
        
        IF(.NOT. GotIt) CYCLE
        IF(.NOT. ActivePart(bf_id) .AND. .NOT. ActivePartAll(bf_id) ) CYCLE

        Model % CurrentElement => Element
        IF ( ActivePart(bf_id) ) THEN
          n = Element % TYPE % NumberOfNodes
          Indexes(1:n) = Element % NodeIndexes
        ELSE
          n = SgetElementDOFs( Indexes )
        END IF
        ValueList => Model % BodyForces(bf_id) % Values

        CALL SetElementLoads( n )
      END DO
    END IF
   
    DEALLOCATE(ActivePart)
    IF(ALLOCATED(DoneLoad)) DEALLOCATE(DoneLoad)


!------------------------------------------------------------------------------
! Go through the point loads that are created on-the-fly
!------------------------------------------------------------------------------

    DO BC=1,Model % NumberOfBCs
      ValueList => Model % BCs(BC) % Values
      IF( .NOT. ListCheckPresent( ValueList,LoadName )) CYCLE
      NodesFound = ListCheckPresent(ValueList,'Target Nodes')

      ! At the first calling the list of coordinates is transformed to list of nodes.
      IF(.NOT. NodesFound) THEN
        CoordNodes => ListGetConstRealArray(ValueList, 'Target Coordinates',GotIt)
        IF(GotIt) THEN
          Eps = ListGetConstReal( ValueList, 'Target Coordinates Eps', Gotit )
          IF ( .NOT. GotIt ) THEN
            Eps = HUGE(Eps)
          ELSE
            ! We are looking at square of distance
            Eps = Eps**2
          END IF

          NoNodes = SIZE(CoordNodes,1)
          NoDims = SIZE(CoordNodes,2)
          
          IF(NoNodes > 0) THEN               
            ALLOCATE( IndNodes(NoNodes), MinDist(NoNodes) )
            IndNodes = -1
            MinDist = HUGE( Dist )
            DO j=1,NoNodes
              DO i=1,Model % NumberOfNodes
                IF( Perm(i) == 0) CYCLE
                
                Dist = (Model % Mesh % Nodes % x(i) - CoordNodes(j,1))**2 
                IF(NoDims >= 2) Dist = Dist + (Model % Mesh % Nodes % y(i) - CoordNodes(j,2))**2
                IF(NoDims == 3) Dist = Dist + (Model % Mesh % Nodes % z(i) - CoordNodes(j,3))**2
                Dist = SQRT(Dist)
                
                IF(Dist < MinDist(j) .AND. Dist <= Eps ) THEN
                  MinDist(j) = Dist
                  IndNodes(j) = i
                END IF
              END DO
            END DO

            ! In parallel case eliminate all except the nearest node. 
            ! This relies on the fact that for each node partition the 
            ! distance to nearest node is computed accurately. 
            DO j=1,NoNodes
              GlobalMinDist = ParallelReduction( MinDist(j), 1 )
              IF(ABS(GlobalMinDist - MinDist(j) )>TINY(Dist)) THEN
                IndNodes(j) = 0
              ELSE IF(ParEnv % PEs>1) THEN
                ! In parallel apply load only on the owner partition:
                ! ---------------------------------------------------
                neigh=>Model % Mesh % ParallelInfo % NeighbourList(IndNodes(j)) % Neighbours
                DO i=1,SIZE(Neigh)
                  IF(ParEnv % Active(neigh(i))) EXIT
                END DO
                IF(neigh(i)/=ParEnv % MyPE) IndNodes(j) = 0
              END IF
            END DO

            NOFNodesFound = 0
            DO j=1,NoNodes
               IF ( IndNodes(j)>0 ) THEN
                 NOFNodesFound = NOFNodesFound+1
                 IndNodes(NOFNodesFound) = IndNodes(j)
               END IF
            END DO
            
            ! In the first time add the found nodes to the list structure
            IF ( NOFNodesFound > 0 ) THEN
              CALL ListAddIntegerArray( ValueList,'Target Nodes', &
                  NOFNodesFound, IndNodes) 
              NodesFound = .TRUE.            
            ELSE
              ! If no nodes found, add still an empty list and make sure the 
              ! zero is not treated later on. Otherwise this search would be 
              ! retreated each time. 
              CALL ListAddIntegerArray( ValueList,'Target Nodes', 0, IndNodes) 
            END IF

            ! Finally deallocate the temporal vectors
            DEALLOCATE( IndNodes, MinDist ) 
          END IF
        END IF
      END IF
      
      IF(NodesFound) THEN           
        CALL Info('SetNodalLoads','Setting nodal loads on target nodes: '//TRIM(Name),Level=9)
        NodeIndexes => ListGetIntegerArray( ValueList,'Target Nodes')
        n = SIZE(NodeIndexes)
        CALL SetPointLoads(n)
      END IF

    END DO

    DEALLOCATE( Indexes )
    !IF(.NOT.ASSOCIATED(A % DiagScaling,DiagScaling)) DEALLOCATE(DiagScaling)

    CALL Info('SetNodalLoads','Finished checking for nodal loads',Level=12)


CONTAINS

     SUBROUTINE SetElementLoads(n)
       INTEGER :: n
       REAL(KIND=dp) :: Work(n)
       
       NodeIndexes => Element % NodeIndexes(1:n)
       
       IF ( DOF > 0 ) THEN
         Work(1:n) = ListGetReal( ValueList, LoadName, n, Indexes, gotIt )
         IF ( .NOT. Gotit ) THEN
           Work(1:n) = ListGetReal( ValueList, LoadName(1:nlen) // ' DOFs', n, Indexes, gotIt )
         END IF
       ELSE
         CALL ListGetRealArray( ValueList, LoadName, WorkA, n, Indexes, gotIt )
       END IF

       IF ( gotIt ) THEN

         DO j=1,n
           k = Perm(Indexes(j))
           
           IF ( k > 0 ) THEN
             IF ( DoneLoad(k) ) CYCLE
             DoneLoad(k) = .TRUE.

             IF ( DOF>0 ) THEN
               k = NDOFs * (k-1) + DOF
               IF( ParEnv % Pes > 1 ) THEN
                  IF(  A % ParallelInfo % NeighbourList(k) % Neighbours(1) /= ParEnv % MyPe ) CYCLE
               END IF
               b(k) = b(k) + Work(j) !* DiagScaling(k)
             ELSE
               DO l=1,MIN( NDOFs, SIZE(Worka,1) )
                 k1 = NDOFs * (k-1) + l
                 b(k1) = b(k1) + WorkA(l,1,j) !* DiagScaling(k1)
               END DO
             END IF
           END IF
         END DO
       END IF
       
     END SUBROUTINE SetElementLoads
     
     
     SUBROUTINE SetPointLoads(n)
       INTEGER :: n
       REAL(KIND=dp) :: Work(n)

       IF(n<=0) RETURN
       
       IF ( DOF > 0 ) THEN
         Work(1:n) = ListGetReal( ValueList, LoadName, n, NodeIndexes, gotIt )
       ELSE
         CALL ListGetRealArray( ValueList, LoadName, WorkA, n, NodeIndexes, gotIt )
       END IF
       
       IF ( GotIt ) THEN
         DO j=1,n
           IF ( NodeIndexes(j) > SIZE(Perm) .OR. NodeIndexes(j) < 1 ) THEN
             CALL Warn('SetNodalLoads','Invalid Node Number')
             CYCLE
           END IF
         
           k = Perm(NodeIndexes(j))
           IF ( k > 0 ) THEN
             IF ( DOF>0 ) THEN
               k = NDOFs * (k-1) + DOF
               b(k) = b(k) + Work(j) !* DiagScaling(k)
             ELSE
               DO l=1,MIN( NDOFs, SIZE(WorkA,1) )
                 k1 = NDOFs * (k-1) + l
                 b(k1) = b(k1) + WorkA(l,1,j) !* DiagScaling(k1)
               END DO
             END IF
           END IF
         END DO
       END IF

     END SUBROUTINE SetPointLoads
     
!------------------------------------------------------------------------------
   END SUBROUTINE SetNodalLoads
!------------------------------------------------------------------------------

  !-------------------------------------------------------------------------------
  SUBROUTINE CommunicateDirichletBCs(A)
  !-------------------------------------------------------------------------------
     TYPE(Matrix_t) :: A

     REAL(KIND=dp), ALLOCATABLE :: d_e(:,:), g_e(:)
     INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
     INTEGER :: i,j,k,l,n,nn,ii(ParEnv % PEs), ierr, status(MPI_STATUS_SIZE)

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

     n = COUNT(A % ConstrainedDOF .AND. A % ParallelInfo % Interface)
     ALLOCATE( s_e(n, nn ), r_e(n) )
     ALLOCATE( d_e(n, nn ), g_e(n) )

     CALL CheckBuffer( nn*3*n )

     ii = 0
     DO i=1, A % NumberOfRows
       IF(A % ConstrainedDOF(i) .AND. A % ParallelInfo % Interface(i) ) THEN
          DO j=1,SIZE(A % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = A % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              d_e(ii(k),k) = A % DValues(i)
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
         CALL MPI_BSEND( d_e(1:ii(i),i),ii(i),MPI_DOUBLE_PRECISION,j-1,112,ELMER_COMM_WORLD,ierr )
       END IF
     END DO

     DO i=1, nn
       j = fneigh(i)
       CALL MPI_RECV( n,1,MPI_INTEGER,j-1,110,ELMER_COMM_WORLD, status,ierr )
       IF ( n>0 ) THEN
         IF( n>SIZE(r_e)) THEN
           DEALLOCATE(r_e,g_e)
           ALLOCATE(r_e(n),g_e(n))
         END IF

         CALL MPI_RECV( r_e,n,MPI_INTEGER,j-1,111,ELMER_COMM_WORLD,status,ierr )
         CALL MPI_RECV( g_e,n,MPI_DOUBLE_PRECISION,j-1,112,ELMER_COMM_WORLD, status,ierr )
         DO j=1,n
           k = SearchNode( A % ParallelInfo, r_e(j), Order=A % ParallelInfo % Gorder )
           IF ( k>0 ) THEN
             IF(.NOT. A % ConstrainedDOF(k)) THEN
               CALL ZeroRow(A, k )
               A % Values(A % Diag(k)) = 1._dp
               A % Dvalues(k) = g_e(j)
               A % ConstrainedDOF(k) = .TRUE.
             END IF
           END IF
         END DO
       END IF
     END DO
     DEALLOCATE(s_e, r_e, d_e, g_e)
  !-------------------------------------------------------------------------------
  END SUBROUTINE CommunicateDirichletBCs
  !-------------------------------------------------------------------------------

  
  !-------------------------------------------------------------------------------
  SUBROUTINE EnforceDirichletConditions( Solver, A, b, OffDiagonal ) 
  !------------------------------------------------------------------------------
    IMPLICIT NONE
    TYPE(Solver_t) :: Solver
    TYPE(Matrix_t), POINTER :: A
    REAL(KIND=dp) :: b(:)
    LOGICAL, OPTIONAL :: OffDiagonal
    
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL :: ScaleSystem, DirichletComm, Found, NoDiag
    REAL(KIND=dp), POINTER :: DiagScaling(:)
    REAL(KIND=dp) :: dval, s
    INTEGER :: i,j,k,n
    CHARACTER(*), PARAMETER :: Caller = 'EnforceDirichletConditions'

    
    Params => Solver % Values

    IF(.NOT. ALLOCATED( A % ConstrainedDOF ) ) THEN
      CALL Info(Caller,&
          'ConstrainedDOF not associated, returning...',Level=8)
      RETURN
    END IF

    
    n = COUNT( A % ConstrainedDOF )
    n = NINT( ParallelReduction(1.0_dp * n ) )

    IF( n == 0 ) THEN
      CALL Info(Caller,'No Dirichlet conditions to enforce, exiting!',Level=10)
      RETURN
    END IF
    
    
    IF( PRESENT( OffDiagonal ) ) THEN
      NoDiag = OffDiagonal
    ELSE
      NoDiag = .FALSE.
    END IF

    IF( NoDiag ) THEN
      ScaleSystem = .FALSE.
    ELSE
      ScaleSystem = ListGetLogical(Params,'Linear System Dirichlet Scaling',Found)
      IF(.NOT.Found) THEN
        ScaleSystem = ListGetLogical(Params,'Linear System Scaling',Found)
        IF(.NOT.Found) ScaleSystem=.TRUE.
      END IF
    END IF
          
    IF( ScaleSystem ) THEN
      CALL Info(Caller,'Applying Dirichlet conditions using scaled diagonal',Level=8)
      CALL ScaleLinearSystem(Solver,A,b,ApplyScaling=.FALSE.)
      DiagScaling => A % DiagScaling
    END IF
    
    ! Communicate the Dirichlet conditions for parallel cases since there may be orphans      
    IF ( ParEnv % PEs > 1 ) THEN
      DirichletComm = ListGetLogical( CurrentModel % Simulation, 'Dirichlet Comm', Found)
      IF(.NOT. Found) DirichletComm = .TRUE.
      IF( DirichletComm) CALL CommunicateDirichletBCs(A)
    END IF
    
    ! Eliminate all entries in matrix that may be eliminated in one sweep
    ! If this is an offdiagonal entry this cannot be done.  
    IF ( A % Symmetric .AND. .NOT. NoDiag ) THEN
      CALL CRS_ElimSymmDirichlet(A,b)
    END IF
 
    
    DO k=1,A % NumberOfRows

      IF ( A % ConstrainedDOF(k) ) THEN
        
        dval = A % Dvalues(k) 
        
        IF( ScaleSystem ) THEN
          s = DiagScaling(k)            
          IF( ABS(s) <= TINY(s) ) s = 1.0_dp
        ELSE
          s = 1.0_dp
        END IF
        s = 1._dp / s**2
          
        CALL ZeroRow(A, k)

        ! Off-diagonal entries for a block matrix are neglected since the code will
        ! also go through the diagonal entries where the r.h.s. target value will be set.
        IF(.NOT. NoDiag ) THEN
          CALL SetMatrixElement(A,k,k,s)
          b(k) = s * dval
        END IF

      END IF
    END DO

    ! Deallocate scaling since otherwise it could be misused out of context
    IF (ScaleSystem) DEALLOCATE( A % DiagScaling ) 
        
    CALL Info(Caller,'Dirichlet boundary conditions enforced', Level=12)
    
  END SUBROUTINE EnforceDirichletConditions
!-------------------------------------------------------------------------------

  
END MODULE DirichletUtils
