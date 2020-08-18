!*****************************************************************************/
! *
! *  Elmer, A Finite Element Software for Multiphysical Problems
! *
! *  Copyright 1st April 1995 - , CSC - IT Center for Science Ltd., Finland
! * 
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation; either
! * version 2.1 of the License, or (at your option) any later version.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! * 
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library (in file ../LGPL-2.1); if not, write 
! * to the Free Software Foundation, Inc., 51 Franklin Street, 
! * Fifth Floor, Boston, MA  02110-1301  USA
! *
! *****************************************************************************/
!
!/******************************************************************************
! *
! *  Authors: Juha Ruokolainen, Peter Råback
! *  Email:   Juha.Ruokolainen@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *  Original Date: 02 Apr 2001
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{

!------------------------------------------------------------------------------
!>  Mesh manipulation utilities for *Solver - routines
!------------------------------------------------------------------------------

MODULE MeshUtils

    USE LoadMod
    USE ElementUtils
    USE ElementDescription
    USE Interpolation
    USE ParallelUtils
    USE Types
    IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
!> Allocated one single element. 
!------------------------------------------------------------------------------
   FUNCTION AllocateElement() RESULT( Element )
!------------------------------------------------------------------------------
     TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------
    INTEGER :: istat
!------------------------------------------------------------------------------

     ALLOCATE( Element, STAT=istat )
     IF ( istat /= 0 ) &
        CALL Fatal( 'AllocateElement', 'Unable to allocate a few bytes of memory?' )
     Element % BDOFs    =  0
     Element % NDOFs    =  0
     Element % BodyId   = -1
     Element % Splitted =  0
     Element % hK = 0
     Element % ElementIndex = 0
     Element % StabilizationMk = 0
     NULLIFY( Element % TYPE )
     NULLIFY( Element % PDefs )
     NULLIFY( Element % BubbleIndexes )
     NULLIFY( Element % DGIndexes )
     NULLIFY( Element % NodeIndexes )
     NULLIFY( Element % EdgeIndexes )
     NULLIFY( Element % FaceIndexes )
     NULLIFY( Element % BoundaryInfo )
!------------------------------------------------------------------------------
   END FUNCTION AllocateElement
!------------------------------------------------------------------------------
 
!------------------------------------------------------------------------------
   SUBROUTINE AllocatePDefinitions(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ! Sanity check to avoid memory leaks
     IF (.NOT. ASSOCIATED(Element % PDefs)) THEN
        ALLOCATE(Element % PDefs, STAT=istat)
        IF ( istat /= 0) CALL Fatal('AllocatePDefinitions','Unable to allocate memory')
     ELSE
       CALL Info('AllocatePDefinitions','P element definitions already allocated',Level=10)
     END IF

     ! Initialize fields
     Element % PDefs % P = 0 
     Element % PDefs % TetraType = 0
     Element % PDefs % isEdge = .FALSE.
     Element % PDefs % pyramidQuadEdge = .FALSE.
     Element % PDefs % localNumber = 0
     Element % PDefs % GaussPoints = 0
!------------------------------------------------------------------------------
   END SUBROUTINE AllocatePDefinitions
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
   SUBROUTINE AllocateBoundaryInfo(Element)
!------------------------------------------------------------------------------
     IMPLICIT NONE
     INTEGER :: istat,n

     TYPE(Element_t) :: Element

     ALLOCATE(Element % BoundaryInfo, STAT=istat)
     IF ( istat /= 0) CALL Fatal('AllocateBoundaryInfo','Unable to allocate memory')

     Element % BoundaryInfo % Left => NULL()
     Element % BoundaryInfo % Right => NULL()
     Element % BoundaryInfo % GebhardtFactors => NULL()
     Element % BoundaryInfo % Constraint =  0

!------------------------------------------------------------------------------
   END SUBROUTINE AllocateBoundaryInfo
!------------------------------------------------------------------------------

!> Allocate mesh structure and return handle to it.
!------------------------------------------------------------------------------
   FUNCTION AllocateMesh(NumberOfBulkElements, NumberOfBoundaryElements, &
       NumberOfNodes, InitParallel ) RESULT(Mesh)
!------------------------------------------------------------------------------
     INTEGER, OPTIONAL :: NumberOfBulkElements, NumberOfBoundaryElements, NumberOfNodes
     LOGICAL, OPTIONAL :: InitParallel
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     INTEGER :: istat, i, n
     CHARACTER(*), PARAMETER :: Caller = 'AllocateMesh'
     
     ALLOCATE( Mesh, STAT=istat )
     IF ( istat /= 0 ) CALL Fatal( Caller, 'Unable to allocate a few bytes of memory?' )

!    Nothing computed on this mesh yet!
!    ----------------------------------
     Mesh % SavesDone    = 0
     Mesh % OutputActive = .FALSE.

     Mesh % AdaptiveDepth = 0
     Mesh % Changed   = .FALSE. !  TODO: Change this sometime
     Mesh % Stabilize = .FALSE.
     Mesh % MeshTag = 1

     Mesh % Variables => NULL()
     Mesh % Parent => NULL()
     Mesh % Child => NULL()
     Mesh % Next => NULL()
     Mesh % RootQuadrant => NULL()
     Mesh % Edges => NULL()
     Mesh % Faces => NULL()
     Mesh % Projector => NULL()
     Mesh % NumberOfEdges = 0
     Mesh % NumberOfFaces = 0

     Mesh % NumberOfBulkElements = 0
     Mesh % NumberOfBoundaryElements = 0
     Mesh % Elements => NULL()
     
     Mesh % DiscontMesh = .FALSE.
     Mesh % SingleMesh  = .FALSE.
     Mesh % InvPerm => NULL()

     Mesh % MinFaceDOFs = 1000
     Mesh % MinEdgeDOFs = 1000
     Mesh % MaxFaceDOFs = 0
     Mesh % MaxEdgeDOFs = 0
     Mesh % MaxBDOFs = 0
     Mesh % MaxElementDOFs  = 0
     Mesh % MaxElementNodes = 0

     Mesh % ViewFactors => NULL()

     ALLOCATE( Mesh % Nodes, STAT=istat )
     IF ( istat /= 0 ) CALL Fatal( Caller, 'Unable to allocate a few bytes of memory?' )
     
     NULLIFY( Mesh % Nodes % x )
     NULLIFY( Mesh % Nodes % y )
     NULLIFY( Mesh % Nodes % z )
     Mesh % Nodes % NumberOfNodes = 0
     Mesh % NumberOfNodes = 0
       
     Mesh % NodesOrig => Mesh % Nodes
     NULLIFY( Mesh % NodesMapped )

     Mesh % EntityWeightsComputed = .FALSE.
     Mesh % BCWeight => NULL()
     Mesh % BodyForceWeight => NULL()
     Mesh % BodyWeight => NULL()
     Mesh % MaterialWeight => NULL()
    
     Mesh % ParallelInfo % NumberOfIfDOFs =  0        
     NULLIFY( Mesh % ParallelInfo % GlobalDOFs )
     NULLIFY( Mesh % ParallelInfo % INTERFACE )
     NULLIFY( Mesh % ParallelInfo % NeighbourList )     

     i = 0
     IF( PRESENT( NumberOfBulkElements ) ) THEN       
       Mesh % NumberOfBulkElements = NumberOfBulkElements
       i = i + 1
     END IF
     
     IF( PRESENT( NumberOfBoundaryElements ) ) THEN
       Mesh % NumberOfBoundaryElements = NumberOfBoundaryElements
       i = i + 1
     END IF

     IF( PRESENT( NumberOfNodes ) ) THEN
       Mesh % NumberOfNodes = NumberOfNodes
       i = i + 1
     END IF
     
     IF( i > 0 ) THEN
       IF( i < 3 ) CALL Fatal(Caller,'Either give all or no optional parameters!')
       CALL InitializeMesh( Mesh, InitParallel )         
     END IF       
     
!------------------------------------------------------------------------------
   END FUNCTION AllocateMesh
!------------------------------------------------------------------------------


   ! Initialize mesh structures after the size information has been 
   ! retrieved.
   !----------------------------------------------------------------
   SUBROUTINE InitializeMesh(Mesh, InitParallel)     
     TYPE(Mesh_t), POINTER :: Mesh
     LOGICAL, OPTIONAL :: InitParallel
     
     INTEGER :: i,j,k,NoElems,istat
     TYPE(Element_t), POINTER :: Element
     CHARACTER(*), PARAMETER :: Caller = 'InitializeMesh'
     LOGICAL :: DoParallel
     
     IF( Mesh % NumberOfNodes == 0 ) THEN
       CALL Warn(Caller,'Mesh has zero nodes!')
       RETURN
     ELSE
       CALL Info(Caller,'Number of nodes in mesh: '&
           //TRIM(I2S(Mesh % NumberOfNodes)),Level=8)
     END IF

     CALL Info(Caller,'Number of bulk elements in mesh: '&
         //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=8)        

     CALL Info(Caller,'Number of boundary elements in mesh: '&
         //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=8)        

     Mesh % Nodes % NumberOfNodes = Mesh % NumberOfNodes          

     NoElems = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

     IF( NoElems == 0 ) THEN
       CALL Fatal('InitializeMesh','Mesh has zero elements!')
     END IF

     Mesh % MaxElementDOFs  = 0
     Mesh % MinEdgeDOFs     = 1000
     Mesh % MinFaceDOFs     = 1000
     Mesh % MaxEdgeDOFs     = 0
     Mesh % MaxFaceDOFs     = 0
     Mesh % MaxBDOFs        = 0

     Mesh % DisContMesh = .FALSE.
     Mesh % DisContPerm => NULL()
     Mesh % DisContNodes = 0

     CALL Info(Caller,'Initial number of max element nodes: '&
         //TRIM(I2S(Mesh % MaxElementNodes)),Level=10) 

     ! Allocate the elements
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Elements, NoElems, Caller )

     DO j=1,NoElems        
       Element => Mesh % Elements(j)        

       Element % DGDOFs = 0
       Element % BodyId = 0
       Element % TYPE => NULL()
       Element % BoundaryInfo => NULL()
       Element % PDefs => NULL()
       Element % DGIndexes => NULL()
       Element % EdgeIndexes => NULL()
       Element % FaceIndexes => NULL()
       Element % BubbleIndexes => NULL()
     END DO

     ! Allocate the nodes
     !-------------------------------------------------------------------------
     CALL AllocateVector( Mesh % Nodes % x, Mesh % NumberOfNodes, Caller )
     CALL AllocateVector( Mesh % Nodes % y, Mesh % NumberOfNodes, Caller )
     CALL AllocateVector( Mesh % Nodes % z, Mesh % NumberOfNodes, Caller )
     
     IF( .NOT. PRESENT( InitParallel ) ) RETURN
     IF( .NOT. InitParallel ) RETURN
     
     CALL Info( Caller,'Allocating parallel info',Level=12)
     
     ALLOCATE(Mesh % ParallelInfo % GlobalDOFs(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     ALLOCATE(Mesh % ParallelInfo % INTERFACE(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     ALLOCATE(Mesh % ParallelInfo % NeighbourList(Mesh % NumberOfNodes), STAT=istat )
     IF ( istat /= 0 ) &
         CALL Fatal( Caller, 'Unable to allocate Mesh % ParallelInfo % NeighbourList' )
     DO i=1,Mesh % NumberOfNodes
       NULLIFY(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
     END DO
     
   END SUBROUTINE InitializeMesh


   
!------------------------------------------------------------------------------
   SUBROUTINE GetMaxDefs(Model, Mesh, Element, ElementDef, SolverId, BodyId, Def_Dofs)
!------------------------------------------------------------------------------
     CHARACTER(*) :: ElementDef
     TYPE(Model_t) :: Model
     TYPE(MEsh_t) :: Mesh
     TYPE(Element_t) :: Element
     INTEGER :: SolverId, BodyId, Def_Dofs(:,:)

     TYPE(ValueList_t), POINTER :: Params
     INTEGER :: i, j,k,l, n, slen, Family
     INTEGER, POINTER :: Body_Dofs(:,:)
     LOGICAL  :: stat, Found
     REAL(KIND=dp) :: x,y,z
     TYPE(Solver_t), POINTER  :: Solver
     CHARACTER(MAX_NAME_LEN) :: str, RESULT

     TYPE(ValueList_t), POINTER :: BodyParams
     CHARACTER(MAX_NAME_LEN) :: ElementDefBody
     
     BodyParams => Model % Bodies(BodyId) % Values

     ElementDefBody=ListGetString(BodyParams,'Solver '//TRIM(i2s(SolverId))//': Element',Found )
     IF (Found) THEN
       CALL Info('GetMaxDefs','Element found for body '//TRIM(i2s(BodyId))//' with solver '//TRIM(i2s(SolverId)), Level=5) 
       CALL Info('GetMaxDefs','Default element type is: '//ElementDef, Level=5)
       CALL Info('GetMaxDefs','New element type for this body is now: '//ElementDefBody, Level=5)
       ElementDef=ElementDefBody
     END IF

     Solver => Model % Solvers(SolverId)
     Params => Solver % Values

     IF ( .NOT. ALLOCATED(Solver % Def_Dofs) ) THEN
       ALLOCATE(Solver % Def_Dofs(10,Model % NumberOfBodies,6))
       Solver % Def_Dofs=-1
       Solver % Def_Dofs(:,:,1)=1
     END IF
     Body_Dofs => Solver % Def_Dofs(1:8,BodyId,:)

     j = INDEX(ElementDef, '-') ! FIX this to include elementtypewise defs...
     IF ( j>0 ) RETURN

     j = INDEX( ElementDef, 'n:' )
     IF ( j>0 ) THEN
       READ( ElementDef(j+2:), * ) l
       Body_Dofs(:,1) = l
       Def_Dofs(:,1) = MAX(Def_Dofs(:,1), l)
     END IF
          
      j = INDEX( ElementDef, 'e:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,2) = l
        Def_Dofs(1:8,2) = MAX(Def_Dofs(1:8,2), l )
      END IF
          
      j = INDEX( ElementDef, 'f:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,3) = l
        Def_Dofs(1:8,3) = MAX(Def_Dofs(1:8,3), l )
      END IF
          
      j = INDEX( ElementDef, 'd:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(:,4) = l
        Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4), l )
      ELSE 
        IF ( ListGetLogical( Solver % Values, &
            'Discontinuous Galerkin', stat ) ) THEN
          Body_Dofs(:,4) = 0
          Def_Dofs(1:8,4) = MAX(Def_Dofs(1:8,4),0 )
        END IF
      END IF
          
      j = INDEX( ElementDef, 'b:' )
      IF ( j>0 ) THEN
        READ( ElementDef(j+2:), * ) l
        Body_Dofs(1:8,5) = l
        Def_Dofs(1:8,5) = MAX(Def_Dofs(1:8,5), l )
      END IF
          
      j = INDEX( ElementDef, 'p:' )
      IF ( j>0 ) THEN
        IF ( ElementDef(j+2:j+2) == '%' ) THEN
          n = Element % TYPE % NumberOfNodes
          x = SUM(Mesh % Nodes % x(Element % NodeIndexes))/n
          y = SUM(Mesh % Nodes % y(Element % NodeIndexes))/n
          z = SUM(Mesh % Nodes % z(Element % NodeIndexes))/n
!          WRITE( str, * ) 'cx= ',TRIM(i2s(Element % ElementIndex)),x,y,z
          WRITE( str, * ) 'cx= ',TRIM(i2s(Element % BodyId)),x,y,z
          str = TRIM(str) // '; ' // TRIM(ElementDef(j+3:))//'(cx)'
          slen = LEN_TRIM(str)
          CALL matc(str,RESULT,slen)
          READ(RESULT(1:slen),*) x
          Body_Dofs(:,6) = 0
          Def_Dofs(1:8,6)  = MAX(Def_Dofs(1:8,6),NINT(x))
          Family = Element % TYPE % ElementCode / 100
          Solver % Def_Dofs(Family, BodyId, 6) = &
              MAX(Solver % Def_Dofs(Family, BodyId, 6), NINT(x))
        ELSE
          READ( ElementDef(j+2:), * ) l
          Body_Dofs(:,6) = l
          Def_Dofs(1:8,6) = MAX(Def_Dofs(1:8,6), l )
        END IF
      END IF

!------------------------------------------------------------------------------
  END SUBROUTINE GetMaxDefs
!------------------------------------------------------------------------------


  SUBROUTINE MarkHaloNodes( Mesh, HaloNode, FoundHaloNodes )

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, POINTER :: HaloNode(:)
    LOGICAL :: FoundHaloNodes

    INTEGER :: n,t
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: Indexes(:)
    LOGICAL :: AllocDone

    ! Check whether we need to skip some elements and nodes on the halo boundary 
    ! We don't want to create additional nodes on the nodes that are on the halo only 
    ! since they just would create further need for new halo...
    FoundHaloNodes = .FALSE.
    IF( ParEnv % PEs > 1 ) THEN
      DO t = 1, Mesh % NumberOfBulkElements
        Element => Mesh % Elements(t)
        IF( ParEnv % MyPe /= Element % PartIndex ) THEN
          FoundHaloNodes = .TRUE.
          EXIT
        END IF
      END DO
    END IF


    ! If we have halo check the truly active nodes
    IF( FoundHaloNodes ) THEN
      CALL Info('MarkHaloNodes',&
          'Checking for nodes that are not really needed in bulk assembly',Level=12)

      IF( .NOT. ASSOCIATED( HaloNode ) ) THEN
        ALLOCATE( HaloNode( Mesh % NumberOfNodes ) )
        AllocDone = .TRUE.
      ELSE
        AllocDone = .FALSE.
      END IF

      ! Node is a halo node if it is not needed by any proper element
      HaloNode = .TRUE.
      DO t = 1, Mesh % NumberOfBulkElements     
        Element => Mesh % Elements(t)
        IF( ParEnv % MyPe == Element % PartIndex ) THEN
          Indexes => Element % NodeIndexes
          HaloNode( Indexes ) = .FALSE.
        END IF
      END DO

      n = COUNT( HaloNode ) 
      FoundHaloNodes = ( n > 0 ) 
      CALL Info('MarkHaloNodes','Number of passive nodes in the halo: '&
          //TRIM(I2S(n)),Level=10)

      ! If there are no halo nodes and the allocation was done within this subroutine
      ! then deallocate also. 
      IF( .NOT. FoundHaloNodes .AND. AllocDone ) THEN
        DEALLOCATE( HaloNode ) 
      END IF
    END IF

  END SUBROUTINE MarkHaloNodes



  ! Mark nodes that are associated with at least some boundary element.
  !------------------------------------------------------------------------------
  SUBROUTINE MarkBCNodes(Mesh,BCNode,NoBCNodes)
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, ALLOCATABLE :: BCNode(:)
    INTEGER :: NoBCNodes

    INTEGER :: elem
    TYPE(Element_t), POINTER :: Element

    CALL Info('MarkInterfaceNodes','Marking interface nodes',Level=8)

    IF(.NOT. ALLOCATED( BCNode ) ) THEN
      ALLOCATE( BCNode( Mesh % NumberOfNodes ) )
    END IF
    BCNode = .FALSE. 

    DO elem=Mesh % NumberOfBulkElements + 1, &
        Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements

      Element => Mesh % Elements( elem )         
      !IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE

      BCNode(Element % NodeIndexes) = .TRUE.
    END DO

    NoBCNodes = COUNT( BCNode )

    CALL Info('MarkBCNodes','Number of BC nodes: '//TRIM(I2S(NoBCNodes)),Level=8)

  END SUBROUTINE MarkBCNodes


  

!> Create a discontinuous mesh over requested boundaries.
!> The nodes are duplicated in order to facilitate the discontinuity.
!> The duplicate nodes are not created by default if the connectivity 
!> of the nodes is needed by other bulk elements than those directly 
!> associated with the discontinuous boundaries. 
!------------------------------------------------------------------------------
 SUBROUTINE CreateDiscontMesh( Model, Mesh, DoAlways )

   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL, OPTIONAL :: DoAlways

   INTEGER, POINTER :: DisContPerm(:)
   LOGICAL, ALLOCATABLE :: DisContNode(:), DisContElem(:), ParentUsed(:), &
       MovingNode(:), StayingNode(:)
   LOGICAL :: Found, DisCont, GreedyBulk, GreedyBC, Debug, DoubleBC, UseTargetBodies, &
       UseConsistantBody, LeftHit, RightHit, Moving, Moving2, Set, Parallel
   INTEGER :: i,j,k,l,n,m,t,bc
   INTEGER :: NoNodes, NoDisContElems, NoDisContNodes, &
       NoBulkElems, NoBoundElems, NoParentElems, NoMissingElems, &
       DisContTarget, NoMoving, NoStaying, NoStayingElems, NoMovingElems, &
       NoUndecided, PrevUndecided, NoEdges, Iter, ElemFamily, DecideLimit, &
       ActiveBCs, CandA, CandB, RightBody, LeftBody, ConflictElems
   INTEGER, TARGET :: TargetBody(1)
   INTEGER, POINTER :: Indexes(:),ParentIndexes(:),TargetBodies(:)
   TYPE(Element_t), POINTER :: Element, LeftElem, RightElem, ParentElem, OtherElem
   CHARACTER(MAX_NAME_LEN) :: DiscontFlag
   LOGICAL :: CheckForHalo
   LOGICAL, POINTER :: HaloNode(:)
   TYPE(ValueList_t), POINTER :: BCList

   LOGICAL :: DoneThisAlready = .FALSE.

   IF(.NOT.PRESENT(DoAlways)) THEN
     IF (DoneThisAlready) RETURN
   ELSE 
     IF(.NOT.DoAlways) THEN
       IF (DoneThisAlready) RETURN
     END IF
   END IF
   DoneThisAlready = .TRUE.

   Discont = .FALSE.
   DoubleBC = .FALSE.
   ActiveBCs = 0
   DO bc = 1,Model % NumberOfBCs
     DisCont = ListGetLogical( Model % BCs(bc) % Values,'Discontinuous Boundary',Found )
     ! If the target boundary / periodic bc / mortar bc is zero
     ! it refers to itself. Otherwise the boundary will be doubled.
     IF( DisCont ) THEN
       i = ListGetInteger( Model % BCs(bc) % Values,'Discontinuous BC',Found )
       j = ListGetInteger( Model % BCs(bc) % Values,'Periodic BC',Found )
       k = ListGetInteger( Model % BCs(bc) % Values,'Mortar BC',Found )
       l = ListGetInteger( Model % BCs(bc) % Values,'Contact BC',Found )
       DoubleBC = ( i + j + k + l > 0 )
       ActiveBCs = ActiveBCs + 1
       BCList => Model % BCs(bc) % Values
     END IF
   END DO
   IF(ActiveBCs == 0 ) RETURN
   
   CALL Info('CreateDiscontMesh','Creating discontinuous boundaries')

   IF( ActiveBCs > 1 ) THEN
     CALL Warn('CreateDiscontMesh','Be careful when using more than one > Discontinuous Boundary < !')
   END IF

   Parallel = ( ParEnv % PEs > 1 )

   NoNodes = Mesh % NumberOfNodes
   NoBulkElems = Mesh % NumberOfBulkElements
   NoBoundElems = Mesh % NumberOfBoundaryElements
   
   ALLOCATE( DisContNode(NoNodes))
   ALLOCATE( DisContElem(NoBoundElems))
   ALLOCATE( ParentUsed(NoBulkElems))
   DisContNode = .FALSE.
   DisContElem = .FALSE.
   ParentUsed = .FALSE.
   NoDisContElems = 0
   NoMissingElems = 0


   ! Check whether we need to skip some elements and nodes on the halo boundary 
   ! We might not want to create additional nodes on the nodes that are on the halo only 
   ! since they just would create further need for new halo...
   CheckForHalo = ListGetLogical( Model % Simulation,'No Discontinuous Halo',Found ) 
   IF(.NOT. Found ) CheckForHalo = .TRUE.
   IF( CheckForHalo ) THEN
     HaloNode => NULL()
     CALL MarkHaloNodes( Mesh, HaloNode, CheckForHalo ) 
   END IF

   ! Go over all boundary elements and mark nodes that should be 
   ! discontinuous and nodes that should be continuous 
   DO t = 1, NoBoundElems
     
     Element => Mesh % Elements(NoBulkElems + t)
     Indexes => Element % NodeIndexes
     n = Element % Type % NumberOfNodes

     DisCont = .FALSE.
     DO bc = 1,Model % NumberOfBCs
       IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
         DisCont = ListGetLogical( Model % BCs(bc) % Values,'Discontinuous Boundary',Found )
         IF( DisCont ) EXIT
       END IF
     END DO     
     IF(.NOT. DisCont ) CYCLE
     
     DO i=1,n
       j = Indexes(i) 
       IF( CheckForHalo ) THEN
         IF( HaloNode(j) ) CYCLE
       END IF
       DisContNode(j) = .TRUE.
     END DO
     DisContElem( t ) = .TRUE.
     
     LeftElem => Element % BoundaryInfo % Left
     IF( ASSOCIATED( LeftElem ) ) THEN
       ParentUsed( LeftElem % ElementIndex ) = .TRUE.
     ELSE
       NoMissingElems = NoMissingElems + 1 
     END IF
     
     RightElem => Element % BoundaryInfo % Right
     IF( ASSOCIATED( RightElem ) ) THEN
       ParentUsed( RightElem % ElementIndex ) = .TRUE.
     ELSE
       NoMissingElems = NoMissingElems + 1
     END IF
   END DO
   
   IF( NoMissingElems > 0 ) THEN
     CALL Warn('CreateDiscontMesh','Missing '//TRIM(I2S(NoMissingElems))// &
     ' parent elements in partition '//TRIM(I2S(ParEnv % MyPe))) 
   END IF

   ! Calculate the number of discontinuous nodes and the number of bulk elements 
   ! associated to them. 
   NoDisContElems = COUNT( DiscontElem )
   NoDisContNodes = COUNT( DisContNode ) 
   CALL Info('CreateDiscontMesh','Number of discontinuous boundary elements: '&
       //TRIM(I2S(NoDisContElems)),Level=7)
   CALL Info('CreateDiscontMesh','Number of candicate nodes: '&
       //TRIM(I2S(NoDisContNodes)),Level=7)

   ! By default all nodes that are associated to elements immediately at the discontinuous 
   ! boundary are treated as discontinuous. However, the user may be not be greedy and release
   ! some nodes from the list that are associated also with other non-discontinuous elements.   
   ConflictElems = 0
   IF( NoDiscontNodes > 0 ) THEN
     n = NoDiscontNodes
     
     GreedyBulk = ListGetLogical( Model % Simulation,'Discontinuous Bulk Greedy',Found ) 
     IF(.NOT. Found ) GreedyBulk = .TRUE.     
     
     GreedyBC = ListGetLogical( Model % Simulation,'Discontinuous Boundary Greedy',Found ) 
     IF(.NOT. Found ) GreedyBC = .TRUE.     
     
     IF( .NOT. ( GreedyBC .AND. GreedyBulk ) ) THEN
       CALL Info('CreateDiscontMesh','Applying non-greedy strategies for Discontinuous mesh',Level=12)

       DO t = 1,NoBulkElems+NoBoundElems
         Element => Mesh % Elements(t)

         IF( t <= NoBulkElems ) THEN
           IF( GreedyBulk ) CYCLE
           IF( ParentUsed(t) ) CYCLE
         ELSE
           IF( GreedyBC ) CYCLE
           IF( DiscontElem(t-NoBulkElems) ) CYCLE
           !IF( Element % BoundaryInfo % Constraint == 0 ) CYCLE
           ! Check that this is not an internal BC
           IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) CYCLE
           IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right) ) CYCLE
         END IF
         Indexes => Element % NodeIndexes

         IF( ANY( DisContNode( Indexes ) ) ) THEN
           !PRINT *,'t',Element % BoundaryInfo % Constraint, t,DisContElem(t), &
           !    Indexes, DisContNode( Indexes ) 
           DisContNode( Indexes ) = .FALSE.
           ConflictElems = ConflictElems + 1
         END IF
       END DO
       NoDisContNodes = COUNT( DisContNode ) 
     END IF

     IF( ConflictElems > 0 ) THEN
       CALL Info('CreateDiscontMesh','Conflicting discontinuity in elements: '&
           //TRIM(I2S(ConflictElems)))
     END IF

     IF( NoDiscontNodes < n ) THEN
       CALL Info('CreateDiscontMesh','Number of local discontinuous nodes: '&
           //TRIM(I2S(NoDisContNodes)), Level=12)
     ELSE
       CALL Info('CreateDiscontMesh','All candidate nodes used',Level=12)
     END IF
     
     IF( NoDiscontNodes == 0 ) THEN
       IF( n > 0 .AND. .NOT. GreedyBulk ) THEN
         CALL Info('CreateDiscontMesh','You might want to try the Greedy bulk strategy',Level=3)
       END IF
     END IF
   END IF
   
   i = NINT( ParallelReduction( 1.0_dp * NoDiscontNodes ) )
   CALL Info('CreateDiscontMesh','Number of discontinuous nodes: '&
       //TRIM(I2S(i)),Level=7)

   IF( i == 0 ) THEN
     CALL Warn('CreateDiscontMesh','Nothing to create, exiting...')
     IF( CheckForHalo ) DEALLOCATE( HaloNode ) 
     DEALLOCATE( DiscontNode, DiscontElem, ParentUsed )
     RETURN
   END IF

   ! Ok, we have marked discontinuous nodes, now give them an index. 
   ! This should also create the indexes in parallel.
   DisContPerm => NULL()
   ALLOCATE( DisContPerm(NoNodes) )
   DisContPerm = 0    

   ! We could end up here on an parallel case only
   ! Then we must make the parallel numbering, so jump to the end where this is done. 
   IF( NoDisContNodes == 0 ) THEN
     IF( DoubleBC ) THEN       
       Mesh % DiscontMesh = .FALSE.
       DEALLOCATE( DisContPerm ) 
     ELSE
       Mesh % DisContMesh = .TRUE.
       Mesh % DisContPerm => DisContPerm
       Mesh % DisContNodes = 0
     END IF
     GOTO 200
   END IF
   
   ! Create a table showing nodes that are related to the moving nodes by
   ! the moving elements. 
   ALLOCATE( MovingNode( NoNodes ), StayingNode( NoNodes ) ) 
   MovingNode = .FALSE.
   StayingNode = .FALSE.

   ! For historical reasons there is both single 'body' and multiple 'bodies'
   ! that define on which side of the discontinuity the new nodes will be. 
   DiscontFlag = 'Discontinuous Target Bodies'
   TargetBodies => ListGetIntegerArray( BCList, DiscontFlag, UseTargetBodies ) 
   IF(.NOT. UseTargetBodies ) THEN
     DiscontFlag = 'Discontinuous Target Body'
     TargetBodies => ListGetIntegerArray( BCList, DiscontFlag, UseTargetBodies ) 
   END IF

   ! If either parent is consistently one of the bodies then we can create a discontinuous 
   ! boundary. Note that this currently only works currently in serial!
   IF(.NOT. UseTargetBodies ) THEN
     IF( ParEnv % PEs > 1 ) THEN
       CALL Fatal('CreateDiscontMesh','Please give > Discontinuous Target Bodies < on the BC!')
     END IF
     
     CALL Info('CreateDiscontMesh','Trying to find a dominating parent body',Level=12)

     CandA = -1
     CandB = -1
     DO t=1, NoBoundElems
       IF(.NOT. DisContElem(t) ) CYCLE
       Element => Mesh % Elements(NoBulkElems + t)

       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) THEN
         CALL Fatal('CreateDiscontMesh','Alternative strategy requires all parent elements!')
       END IF
       IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
         CALL Fatal('CreateDiscontMesh','Alternative strategy requires all parent elements!')
       END IF

       LeftBody = Element % BoundaryInfo % Left % BodyId         
       RightBody = Element % BoundaryInfo % Right % BodyId

       IF( CandA == -1 ) THEN
         CandA = LeftBody 
       ELSE IF( CandA == 0 ) THEN
         CYCLE
       ELSE IF( CandA /= LeftBody .AND. CandA /= RightBody ) THEN
         CandA = 0
       END IF

       IF( CandB == -1 ) THEN
         CandB = RightBody
       ELSE IF( CandB == 0 ) THEN
         CYCLE
       ELSE IF( CandB /= LeftBody .AND. CandB /= RightBody ) THEN
         CandB = 0
       END IF
     END DO

     ! Choose the bigger one to honor the old convention
     ! This eliminates at the same time the unsuccessful case of zero.
     TargetBody(1) = MAX( CandA, CandB )

     IF( TargetBody(1) > 0 ) THEN
       CALL Info('CreateDiscontMesh',&
           'There seems to be a consistent discontinuous body: '&
           //TRIM(I2S(TargetBody(1))),Level=8)
       UseConsistantBody = .TRUE.
       TargetBodies => TargetBody
     ELSE
       CALL Fatal('CreateDiscontMesh',&
           'No simple rules available for determining discontinuous body')
     END IF
   END IF


   ! Assume we have only one active BC and we know the list of discontinuous 
   ! target bodies there. Hence we have all the info needed to set the 
   ! discontinuous elements also for other bulk elements. 
   ! This could be made more generic...
   NoUndecided = 0
   NoMovingElems = 0 
   NoStayingElems = 0

   DO t=1, NoBulkElems
     Element => Mesh % Elements(t)

     ! No need to treat halo elements
     !IF( CheckForHalo .AND. Element % PartIndex /= ParEnv % MyPe ) CYCLE

     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE
     Moving = ANY( TargetBodies == Element % BodyId )

     IF( Moving ) THEN
       NoMovingElems = NoMovingElems + 1 
       MovingNode(Indexes) = .TRUE.
     ELSE
       StayingNode(Indexes) = .TRUE.
       NoStayingElems = NoStayingElems + 1
     END IF
   END DO

   CALL Info('CreateDiscontMesh','Number of bulk elements moving: '&
       //TRIM(I2S(NoMovingElems)), Level=8)
   CALL Info('CreateDiscontMesh','Number of bulk elements staying: '&
       //TRIM(I2S(NoStayingElems)), Level=8)

   ! Set discontinuous nodes only if there is a real moving node associted with it
   ! Otherwise we would create a zero to the permutation vector. 
   ! If there is just a staying node then no need to create discontinuity at this node.
   DiscontNode = DiscontNode .AND. MovingNode 

   ! Create permutation numbering for the discontinuous nodes   
   ! Doubling will be done only for nodes that have both parents
   j = 0
   DO i=1,NoNodes
     IF( DisContNode(i) ) THEN
       j = j + 1
       DisContPerm(i) = j
     END IF
   END DO
   IF( j < NoDiscontNodes ) THEN
     PRINT *,'Some discontinuous nodes only needed on the other side:',&
         ParEnv % MyPe, NoDiscontNodes-j
     NoDiscontNodes = j 
   END IF


   ! Now set the new indexes for bulk elements
   ! In parallel skip the halo elements
   DO t=1, NoBulkElems
     Element => Mesh % Elements(t)

     ! No need to treat halo elements
     !IF( CheckForHalo .AND. Element % PartIndex /= ParEnv % MyPe ) CYCLE
     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE
     Moving = ANY( TargetBodies == Element % BodyId )

     IF( Moving ) THEN
       DO i=1, SIZE(Indexes) 
         j = DisContPerm(Indexes(i))
         IF( j > 0 ) Indexes(i) = NoNodes + j
       END DO
     END IF
   END DO

    
   ! Now set also the unset boundary elements by following the ownership of the parent elements
   ! or the majority opinion if this is conflicting.
   DO t=1, NoBoundElems

     Element => Mesh % Elements(NoBulkElems + t)

     ! If the element has no constraint then there is no need to treat it
     IF( Element % BoundaryInfo % Constraint == 0 ) CYCLE

     IF( DisContElem(t) ) THEN
       LeftElem => Element % BoundaryInfo % Left
       RightElem => Element % BoundaryInfo % Right

       IF( ASSOCIATED( LeftElem ) ) THEN
         Moving = ANY( TargetBodies == LeftElem % BodyId ) 
       ELSE
         Moving = .NOT. ANY( TargetBodies == RightElem % BodyId )
       END IF
       IF( Moving ) THEN
         Element % BoundaryInfo % Left => RightElem
         Element % BoundaryInfo % Right => LeftElem 
       END IF
       CYCLE
     END IF


     Indexes => Element % NodeIndexes

     IF( .NOT. ANY( DisContNode( Indexes ) ) ) CYCLE

     ElemFamily = Element % TYPE % ElementCode / 100 
     LeftElem => Element % BoundaryInfo % Left
     RightElem => Element % BoundaryInfo % Right

     ! The boundary element follows the parent element if it is clear what to do
     Set = .TRUE.
     IF( ASSOCIATED( LeftElem ) .AND. ASSOCIATED( RightElem ) ) THEN
       Moving = ANY( TargetBodies == LeftElem % BodyId )
       Moving2 = ANY( TargetBodies == RightElem % BodyId ) 
       IF( Moving .NEQV. Moving2) THEN
         CALL Warn('CreateDiscontMesh','Conflicting moving information')
         !PRINT *,'Moving:',t,Element % BoundaryInfo % Constraint, &
         !    Moving,Moving2,LeftElem % BodyId, RightElem % BodyId
         Set = .FALSE.
       ELSE
         IF( Moving ) THEN
           Element % BoundaryInfo % Left => RightElem
           Element % BoundaryInfo % Right => LeftElem 
         END IF
       END IF
     ELSE IF( ASSOCIATED( LeftElem ) ) THEN
       Moving = ANY( LeftElem % NodeIndexes > NoNodes ) 
     ELSE IF( ASSOCIATED( RightElem ) ) THEN
       Moving = ANY( RightElem % NodeIndexes > NoNodes )
     ELSE
       CALL Fatal('CreateDiscontMesh','Boundary BC has no parants!')
     END IF

     ! Otherwise we follow the majority rule
     IF( .NOT. Set ) THEN
       NoMoving = COUNT( MovingNode(Indexes) ) 
       NoStaying = COUNT( StayingNode(Indexes) ) 

       IF( NoStaying /= NoMoving ) THEN
         Moving = ( NoMoving > NoStaying )
         Set = .TRUE.
       END IF
     END IF

     ! Ok, finally set whether boundary element is moving or staying
     IF( Set ) THEN
       IF( Moving ) THEN
         NoMovingElems = NoMovingElems + 1 
         DO i=1, SIZE(Indexes) 
           j = DisContPerm(Indexes(i))
           IF( j > 0 ) Indexes(i) = NoNodes + j
         END DO
       ELSE
         NoStayingElems = NoStayingElems + 1
       END IF
     ELSE
       NoUndecided = NoUndecided + 1
     END IF
   END DO

   CALL Info('CreateDiscontMesh','Number of related elements moving: '&
       //TRIM(I2S(NoMovingElems)), Level=8 )
   CALL Info('CreateDiscontMesh','Number of related elements staying: '&
       //TRIM(I2S(NoStayingElems)), Level=8 )
   IF( NoUndecided == 0 ) THEN
     CALL Info('CreateDiscontMesh','All elements marked either moving or staying')
   ELSE
     CALL Info('CreateDiscontMesh','Number of related undecided elements: '//TRIM(I2S(NoUndecided)) )
     CALL Warn('CreateDiscontMesh','Could not decide what to do with some boundary elements!')
   END IF


   m = COUNT( DiscontNode .AND. .NOT. MovingNode )
   IF( m > 0 ) THEN
     PRINT *,'Number of discont nodes not moving: ',ParEnv % MyPe, m
   END IF

   m = COUNT( DiscontNode .AND. .NOT. StayingNode )
   IF( m > 0 ) THEN
     PRINT *,'Number of discont nodes not staying: ',ParEnv % MyPe, m
     DO i=1,SIZE(DisContNode)
       IF( DiscontNode(i) .AND. .NOT. StayingNode(i) ) THEN
         IF( ParEnv % PEs == 1 ) THEN
           PRINT *,'Node:',ParEnv % MyPe,i
         ELSE
           PRINT *,'Node:',ParEnv % MyPe,i,Mesh % ParallelInfo % GlobalDofs(i), &
               Mesh % ParallelInfo % NeighbourList(i) % Neighbours
         END IF
         PRINT *,'Coord:',ParEnv % MyPe, Mesh % Nodes % x(i), Mesh % Nodes % y(i)
       END IF
     END DO
   END IF

   !DEALLOCATE( MovingNode, StayingNode )

   ! Now add the new nodes also to the nodes structure
   ! and give the new nodes the same coordinates as the ones
   ! that they were derived from. 
   Mesh % NumberOfNodes = NoNodes + NoDisContNodes   
   CALL EnlargeCoordinates( Mesh ) 

   CALL Info('CreateDiscontMesh','Setting new coordinate positions',Level=12)
   DO i=1, NoNodes
     j = DisContPerm(i)
     IF( j > 0 ) THEN
       k = NoNodes + j
       Mesh % Nodes % x(k) = Mesh % Nodes % x(i)
       Mesh % Nodes % y(k) = Mesh % Nodes % y(i)
       Mesh % Nodes % z(k) = Mesh % Nodes % z(i)
     END IF
   END DO


   ! If the discontinuous boundary is duplicated then no information of it 
   ! is saved. The periodic and mortar conditions now need to perform
   ! searches. On the other hand the meshes may now freely move.,
   IF( DoubleBC ) THEN
     CALL Info('CreateDiscontMesh','Creating secondary boundary for Discontinuous gap',Level=10)

     CALL EnlargeBoundaryElements( Mesh, NoDiscontElems ) 

     NoDisContElems = 0
     DO t=1, NoBoundElems

       ! Is this a boundary to be doubled?
       IF(.NOT. DisContElem(t) ) CYCLE

       Element => Mesh % Elements(NoBulkElems + t)
       IF(.NOT. ASSOCIATED(Element) ) THEN
         CALL Fatal('CreateDiscontMesh','Element '//TRIM(I2S(NoBulkElems+t))//' not associated!')
       END IF
       Indexes => Element % NodeIndexes

       DisContTarget = 0
       Found = .FALSE.
       DO bc = 1,Model % NumberOfBCs
         IF ( Element % BoundaryInfo % Constraint == Model % BCs(bc) % Tag ) THEN
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Discontinuous BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Mortar BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Periodic BC',Found )
           IF( Found ) EXIT
           DisContTarget = ListGetInteger( Model % BCs(bc) % Values,&
               'Contact BC',Found )
           IF( Found ) EXIT
         END IF
       END DO
       IF( .NOT. Found .OR. DisContTarget == 0 ) THEN
         CALL Fatal('CreateDiscontMesh','Nonzero target boundary must be given for all, if any bc!')
       END IF

       RightElem => Element % BoundaryInfo % Right
       LeftElem => Element % BoundaryInfo % Left 

       NoDisContElems = NoDisContElems + 1              
       j = NoBulkElems + NoBoundElems + NoDisContElems 

       OtherElem => Mesh % Elements( j )
       IF(.NOT. ASSOCIATED(OtherElem) ) THEN
         CALL Fatal('CreateDiscontMesh','Other elem '//TRIM(I2S(j))//' not associated!')
       END IF

       OtherElem = Element 
       OtherElem % TYPE => Element % TYPE

       NULLIFY( OtherElem % BoundaryInfo ) 
       ALLOCATE( OtherElem % BoundaryInfo ) 
       OtherElem % BoundaryInfo % Left => Element % BoundaryInfo % Right

       ! Now both boundary elements are just one sided. Remove the associated to the other side. 
       NULLIFY( Element % BoundaryInfo % Right ) 
       NULLIFY( OtherElem % BoundaryInfo % Right )

       NULLIFY( OtherElem % NodeIndexes )
       n = SIZE( Element % NodeIndexes ) 
       ALLOCATE( OtherElem % NodeIndexes( n ) ) 

       ! Ok, we found the element to manipulate the indexes. 
       ! The new index is numbered on top of the old indexes. 
       DO i=1,n
         j = Element % NodeIndexes(i) 
         IF( DisContPerm(j) > 0 ) THEN
           OtherElem % NodeIndexes(i) = NoNodes + DisContPerm(j)
         ELSE 
           OtherElem % NodeIndexes(i) = j
         END IF
       END DO

       OtherElem % BoundaryInfo % Constraint = DisContTarget
     END DO

     CALL Info('CreateDiscontMesh','Number of original bulk elements: '&
         //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=10)
     CALL Info('CreateDiscontMesh','Number of original boundary elements: '&
         //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=10)
     CALL Info('CreateDiscontMesh','Number of additional boundary elements: '&
         //TRIM(I2S(NoDisContElems)),Level=10)

     Mesh % DiscontMesh = .FALSE.
   ELSE
     Mesh % DisContMesh = .TRUE.
     Mesh % DisContPerm => DisContPerm
     Mesh % DisContNodes = NoDisContNodes 
   END IF

200 CONTINUE


   CALL EnlargeParallelInfo(Mesh, DiscontPerm )
   IF( ParEnv % PEs > 1 ) THEN
     m = COUNT( Mesh % ParallelInfo % GlobalDofs == 0) 
     IF( m > 0 ) CALL Warn('CreateDiscontMesh','There are nodes with zero global dof index: '//TRIM(I2S(m)))
   END IF

   IF( DoubleBC .AND. NoDiscontNodes > 0 ) DEALLOCATE( DisContPerm )


   DEALLOCATE( DisContNode, DiscontElem )   
  
 END SUBROUTINE CreateDiscontMesh


!> Reallocate coordinate arrays for iso-parametric p-elements,
!> or if the size of nodes has been increased due to discontinuity. 
!> This does not seem to be necessary for other types of 
!> elements (face, edge, etc.)
! -----------------------------------------------------------    
 SUBROUTINE EnlargeCoordinates(Mesh)

   TYPE(Mesh_t) :: Mesh
   INTEGER :: n0, n
   REAL(KIND=dp), POINTER :: TmpCoord(:)

   INTEGER :: i
   LOGICAL :: pelementsPresent

   n = Mesh % NumberOfNodes + &
       Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
       Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
       Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements
   n0 = SIZE( Mesh % Nodes % x )

   pelementsPresent = .FALSE.
   DO i=1,Mesh % NumberOfBulkElements
     IF(isPelement(Mesh % Elements(i))) THEN
       pelementsPresent = .TRUE.; EXIT
     END IF
   END DO

   IF ( Mesh % NumberOfNodes > n0 .OR. n > n0 .AND. pelementsPresent ) THEN
     CALL Info('EnlargeCoordinates','Increasing number of nodes from '&
         //TRIM(I2S(n0))//' to '//TRIM(I2S(n)),Level=8)

     TmpCoord => Mesh % Nodes % x
     ALLOCATE( Mesh % Nodes % x(n) )
     Mesh % Nodes % x(1:n0) = TmpCoord
     Mesh % Nodes % x(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )

     TmpCoord => Mesh % Nodes % y
     ALLOCATE( Mesh % Nodes % y(n) )
     Mesh % Nodes % y(1:n0) = TmpCoord
     Mesh % Nodes % y(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )

     TmpCoord => Mesh % Nodes % z
     ALLOCATE( Mesh % Nodes % z(n) )
     Mesh % Nodes % z(1:n0) = TmpCoord
     Mesh % Nodes % z(n0 + 1:n) = 0.0_dp
     DEALLOCATE( TmpCoord )
   END IF

 END SUBROUTINE EnlargeCoordinates


 
 SUBROUTINE EnlargeBoundaryElements(Mesh, DoubleElements )

   TYPE(Mesh_t) :: Mesh
   INTEGER :: DoubleElements
   INTEGER :: n,n0,i,j
   REAL(KIND=dp), POINTER :: TmpCoord(:)
   TYPE(Element_t), POINTER :: NewElements(:),OldElements(:), Element

   IF( DoubleElements == 0 ) RETURN

   n0 = Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
   n = n0 + DoubleElements

   CALL Info('EnlargeBoundaryElements','Increasing number of elements from '&
       //TRIM(I2S(n0))//' to '//TRIM(I2S(n)),Level=8)

   OldElements => Mesh % Elements
   CALL AllocateVector( Mesh % Elements, n, 'EnlargeBoundaryElements' )
   DO i=1,n0
     Mesh % Elements(i) = OldElements(i)
     IF(ASSOCIATED(OldElements(i) % BoundaryInfo)) THEN
       IF (ASSOCIATED(OldElements(i) % BoundaryInfo % Left)) &
           Mesh % Elements(i) % BoundaryInfo % Left => &
           Mesh % Elements(OldElements(i) % BoundaryInfo % Left % ElementIndex)
       
       IF (ASSOCIATED(OldElements(i) % BoundaryInfo % Right)) &
           Mesh % Elements(i) % BoundaryInfo % Right => &
           Mesh % Elements(OldElements(i) % BoundaryInfo % Right % ElementIndex)
     END IF
   END DO

   DO i=n0+1,n
     Element => Mesh % Elements(i)

     Element % DGDOFs = 0
     Element % BodyId = 0
     Element % TYPE => NULL()
     Element % BoundaryInfo => NULL()
     Element % PDefs => NULL()
     Element % DGIndexes => NULL()
     Element % EdgeIndexes => NULL()
     Element % FaceIndexes => NULL()
     Element % BubbleIndexes => NULL()
   END DO

   DEALLOCATE( OldElements ) 
   Mesh % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements + DoubleElements

 END SUBROUTINE EnlargeBoundaryElements


 SUBROUTINE EnlargeParallelInfo( Mesh, DiscontPerm )

   TYPE(Mesh_t) :: Mesh
   INTEGER, POINTER :: DiscontPerm(:)

   INTEGER :: nmax,n0,n1,i,j,istat, goffset
   INTEGER, POINTER :: TmpGlobalDofs(:) 
   INTEGER, ALLOCATABLE :: Perm(:)
   LOGICAL, POINTER :: Intf(:)
   TYPE(NeighbourList_t), POINTER :: Nlist(:)

   IF ( ParEnv % PEs <= 1 ) RETURN

   ! As index offset use the number of nodes in the whole mesh
   goffset = ParallelReduction( MAXVAL(Mesh % ParallelInfo % GlobalDofs)*1._dp,2 )

   n0 = SIZE( Mesh % ParallelInfo % GlobalDofs )
   n1 = Mesh % NumberOfNodes 
   IF( n0 >= n1 ) THEN
     CALL Info('EnlargeParallelInfo','No need to grow: '&
         //TRIM(I2S(n0))//' vs. '//TRIM(I2S(n1)),Level=10)
     RETURN
   END IF
   
   CALL Info('EnlargeParallelInfo','Increasing global numbering size from '&
         //TRIM(I2S(n0))//' to '//TRIM(I2S(n1)),Level=8)

   ! Create permutation table for the added nodes
   ALLOCATE(Perm(n1)); Perm  = 0
   DO i=1,n0
     IF ( DiscontPerm(i) > 0 ) THEN
       Perm(DiscontPerm(i)+n0) = i
     END IF
   END DO

   ! Create the enlarged set of global nodes indexes
   ALLOCATE( TmpGlobalDofs(n1), STAT=istat )
   IF (istat /= 0) CALL Fatal('LoadMesh', 'Unable to allocate TmpGlobalDofs array.')
   TmpGlobalDofs = 0
   DO i=1,n0
     TmpGlobalDofs(i) = Mesh % ParallelInfo % GlobalDofs(i)
   END DO
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0) THEN
       TmpGlobalDofs(i) = TmpGlobalDOfs(j) + goffset
     END IF
   END DO
   DEALLOCATE(Mesh % ParallelInfo % GlobalDofs)
   Mesh % ParallelInfo % GlobalDOfs => TmpGlobalDofs

   ! Create the enlarged list of neighbours
   ALLOCATE(Nlist(n1))
   DO i=1,n0
     IF( ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) THEN
       Nlist(i) % Neighbours => &
           Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       Mesh % ParallelInfo % NeighbourList(i) % Neighbours => NULL()
     ELSE 
       Nlist(i) % Neighbours => NULL()
     END IF
   END DO

   DO i=n0+1,n1
     j = Perm(i)
     IF ( j > 0 ) THEN
       IF( ASSOCIATED( Nlist(j) % Neighbours ) ) THEN
         ALLOCATE( Nlist(i) % Neighbours(SIZE(Nlist(j) % Neighbours) ) )
         Nlist(i) % Neighbours = Nlist(j) % Neighbours
       ELSE
         Nlist(i) % Neighbours => NULL()
       END IF
     END IF
   END DO
   DEALLOCATE(Mesh % ParallelInfo % NeighbourList)
   Mesh % ParallelInfo % NeighbourList => Nlist


   ! Create logical table showing the interface nodes
   ALLOCATE( Intf(n1) )
   Intf = .FALSE.
   Intf(1:n0) = Mesh % ParallelInfo % INTERFACE(1:n0)
   DO i=n0+1,n1
     j = Perm(i)
     IF(j > 0 ) THEN
       Intf(i) = Intf(j) 
     END IF
   END DO
   DEALLOCATE( Mesh % ParallelInfo % INTERFACE )
   Mesh % ParallelInfo % Interface => Intf


 END SUBROUTINE EnlargeParallelInfo




 !> Fortran reader for Elmer ascii mesh file format.
 !> This is a Fortran replacement for the old C++ eio library. 
 !------------------------------------------------------------------------
 SUBROUTINE ElmerAsciiMesh(Step, PMesh, MeshNamePar, ThisPe, NumPEs, IsParallel )

   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel

   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: PrevStep=0, iostat
   INTEGER, PARAMETER :: FileUnit = 10
   CHARACTER(MAX_NAME_LEN) :: BaseName, FileName
   INTEGER :: i,j,k,n,BaseNameLen, SharedNodes = 0, mype = 0, numprocs = 0
   INTEGER, POINTER :: NodeTags(:), ElementTags(:), LocalPerm(:)
   INTEGER :: MinNodeTag = 0, MaxNodeTag = 0, istat
   LOGICAL :: ElementPermutation=.FALSE., NodePermutation=.FALSE., Parallel



   SAVE PrevStep, BaseName, BaseNameLen, Mesh, mype, Parallel, &
       NodeTags, ElementTags, LocalPerm

   CALL Info('ElmerAsciiMesh','Performing step: '//TRIM(I2S(Step)),Level=8)

   IF( Step - PrevStep /= 1 ) THEN
     CALL Fatal('ElmerAsciiMesh','The routine should be called in sequence: '// &
         TRIM(I2S(PrevStep))//' : '//TRIM(I2S(Step)) )
   END IF
   PrevStep = Step
   IF( PrevStep == 6 ) PrevStep = 0 

   IF( Step == 1 ) THEN
     IF(.NOT. PRESENT( MeshNamePar ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give MeshNamePar!')
     END IF
     BaseName = TRIM( MeshNamePar ) 
     IF(.NOT. PRESENT( PMesh ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give PMesh!')
     END IF
     Mesh => PMesh
     IF(.NOT. PRESENT( ThisPe ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give ThisPe!')
     END IF
     mype = ThisPe 
     IF(.NOT. PRESENT( NumPEs) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give NumPEs!')
     END IF
     numprocs = NumPEs
     IF(.NOT. PRESENT( IsParallel ) ) THEN
       CALL Fatal('ElmerAsciiMesh','When calling in mode one give IsParallel!')
     END IF
     Parallel = IsParallel

     i = LEN_TRIM(MeshNamePar)
     DO WHILE(MeshNamePar(i:i) == CHAR(0))
       i=i-1
     END DO
     BaseNameLen = i
     CALL Info('LoadMesh','Base mesh name: '//TRIM(MeshNamePar(1:BaseNameLen)))
   END IF


   SELECT CASE( Step ) 

   CASE(1)       
     CALL ReadHeaderFile()

   CASE(2)
     CALL ReadNodesFile()

   CASE(3)
     CALL ReadElementsFile()

   CASE(4)
     CALL ReadBoundaryFile()
     CALL PermuteNodeNumbering()

   CASE(5)
     CALL InitParallelInfo()
     CALL ReadSharedFile()

   CASE(6)
     IF( ASSOCIATED( LocalPerm) ) DEALLOCATE( LocalPerm ) 
     IF( ASSOCIATED( ElementTags) ) DEALLOCATE( ElementTags )

   END SELECT


 CONTAINS


   FUNCTION read_ints(s,j,halo) RESULT(n)
     INTEGER :: j(:)
     CHARACTER(LEN=*) :: s
     LOGICAL :: halo
     
     INTEGER :: i,k,l,m,n,ic
     INTEGER, PARAMETER :: ic0 = ICHAR('0'), ic9 = ICHAR('9'), icm = ICHAR('-'), &
         icd = ICHAR('/'), ics = ICHAR(' ')
     
     k = LEN_TRIM(s)
     l = 1
     n = 0
     halo = .FALSE.
     DO WHILE(l<=k.AND.n<SIZE(j))
       DO WHILE(l<=k)
         ic = ICHAR(s(l:l))
         IF( ic == ics ) THEN
           CONTINUE
         ELSE IF( ic == icd ) THEN
           halo = .TRUE.
         ELSE
           EXIT
         END IF
         l=l+1
       END DO
       IF(l>k) EXIT
       IF(.NOT.(ic==icm .OR. ic>=ic0 .AND. ic<=ic9)) EXIT
       
       m = l+1
       DO WHILE(m<=k)
         ic = ICHAR(s(m:m))
         IF(ic<ic0 .OR. ic>ic9) EXIT
         m=m+1
       END DO
       
       n = n + 1
       j(n) = s2i(s(l:m-1),m-l)
       l = m
     END DO
   END FUNCTION read_ints
   

   !---------------------------------------------------
   ! Read header file and allocate some mesh structures
   !---------------------------------------------------
   SUBROUTINE ReadHeaderFile()

     INTEGER :: TypeCount
     INTEGER :: Types(64),CountByType(64)

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(numprocs))//&
           '/part.'//TRIM(I2S(mype+1))//'.header'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.header'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading header info from file: '//TRIM(FileName),Level=10)
     END IF

     READ(FileUnit,*,IOSTAT=iostat) Mesh % NumberOfNodes, &
         Mesh % NumberOfBulkElements,&
         Mesh % NumberOfBoundaryElements
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not read header 1st line in file: '//TRIM(FileName))
     END IF

     Types = 0
     CountByType = 0
     READ(FileUnit,*,IOSTAT=iostat) TypeCount
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not read the type count in file: '//TRIM(FileName))
     END IF
     DO i=1,TypeCount
       READ(FileUnit,*,IOSTAT=iostat) Types(i),CountByType(i)
       IF( iostat /= 0 ) THEN
         CALL Fatal('LoadMesh','Could not read type count '&
             //TRIM(I2S(i))//'in file: '//TRIM(FileName))
       END IF
     END DO

     IF( Parallel ) THEN
       READ(FileUnit,*,IOSTAT=iostat) SharedNodes
       IF( iostat /= 0 ) THEN
         CALL Fatal('LoadMesh','Could not read shared nodes in file: '//TRIM(FileName))
       END IF
     ELSE
       SharedNodes = 0
     END IF

     Mesh % MaxElementNodes = 0
     DO i=1,TypeCount
       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes, MODULO( Types(i), 100) )
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadHeaderFile


   !-----------------------------------------------------------------------
   ! Read nodes file and create nodal permutation if needed
   !-----------------------------------------------------------------------
   SUBROUTINE ReadNodesFile()

     REAL(KIND=dp) :: Coords(3)
     INTEGER :: NodeTag

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(numprocs))//&
           '/part.'//TRIM(I2S(mype+1))//'.nodes'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.nodes'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ALLOCATE( NodeTags(Mesh % NumberOfNodes ) ) 
     NodeTags = 0

     NodePermutation = .FALSE.
     DO j = 1, Mesh % NumberOfNodes
       READ(FileUnit,*,IOSTAT=iostat) NodeTag, k, Coords
       IF( iostat /= 0 ) THEN
         CALL Fatal('LoadMesh','Problem load node '//TRIM(I2S(j))//' in file: '//TRIM(Filename))
       END IF

       IF( NodeTags(j) /= j ) NodePermutation = .TRUE.
 
       NodeTags(j) = NodeTag
       Mesh % Nodes % x(j) = Coords(1)
       Mesh % Nodes % y(j) = Coords(2)
       Mesh % Nodes % z(j) = Coords(3)
     END DO

     CLOSE(FileUnit)

   END SUBROUTINE ReadNodesFile


   !------------------------------------------------------------------------------
   ! Read elements file and create elemental permutation if needed 
   !------------------------------------------------------------------------------
   SUBROUTINE ReadElementsFile()
     TYPE(Element_t), POINTER :: Element
     INTEGER :: ElemType, Tag, Body, ElemNo, Ivals(64),nread, ioffset, partn
     CHARACTER(256) :: str
     LOGICAL :: halo


     CALL AllocateVector( ElementTags, Mesh % NumberOfBulkElements+1, 'LoadMesh')   
     ElementTags = 0
     ElementPermutation = .FALSE.

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)// &
          '/partitioning.'//TRIM(I2S(numprocs))//&
             '/part.'//TRIM(I2S(mype+1))//'.elements'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.elements'
     END IF

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', iostat=IOSTAT )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadElementsFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading bulk elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=1,Mesh % NumberOfBulkElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read start of element entry: '//TRIM(I2S(j)))
       END IF

       nread = read_ints(str,ivals,halo)

       tag = ivals(1)

       IF( halo ) THEN
         ioffset = 1
         partn = ivals(2) 
       ELSE
         ioffset = 0
         partn = 0 
       END IF
       body = ivals(ioffset+2)
       ElemType = ivals(ioffset+3)

       ElementTags(j) = tag
       IF( j /= tag ) ElementPermutation = .TRUE.             
       Element % ElementIndex = j
       Element % BodyId = body

       IF( partn > 0 ) THEN
         Element % PartIndex = partn-1
       ELSE
         Element % PartIndex = mype
       END IF

       Element % TYPE => GetElementType(ElemType)

       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadElementsFile','Element of type '&
             //TRIM(I2S(ElemType))//' could not be associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       IF( nread < n + ioffset + 3 ) THEN
         CALL Fatal('ReadElementsFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF

       CALL AllocateVector( Element % NodeIndexes, n )

       Element % NodeIndexes(1:n) = IVals(4+ioffset:nread)
     END DO
     CLOSE( FileUnit ) 

   END SUBROUTINE ReadElementsFile
   !------------------------------------------------------------------------------


   !------------------------------------------------------------------------------
   ! Read boundary elements file and remap the parents if needed.  
   !------------------------------------------------------------------------------
   SUBROUTINE ReadBoundaryFile()
     INTEGER, POINTER :: LocalEPerm(:)
     INTEGER :: MinEIndex, MaxEIndex, ElemNodes, i
     INTEGER :: Left, Right, bndry, tag, ElemType, IVals(64), nread, ioffset, partn
     TYPE(Element_t), POINTER :: Element
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF( Parallel ) THEN
       FileName = BaseName(1:BaseNameLen)//&
          '/partitioning.'//TRIM(I2S(numprocs))//&
           '/part.'//TRIM(I2S(mype+1))//'.boundary'
     ELSE
       FileName = BaseName(1:BaseNameLen)//'/mesh.boundary'
     END IF

     ! Create permutation for the elements. This is needed when the element 
     ! parents are mapped to the new order. This is needed for mapping of the 
     ! parents. Otherwise the element numbering is arbitrary. 
     !------------------------------------------------------------------------------
     IF( ElementPermutation ) THEN
       MinEIndex = MINVAL( ElementTags(1:Mesh % NumberOfBulkElements) )
       MaxEIndex = MAXVAL( ElementTags(1:Mesh % NumberOfBulkElements) )

       LocalEPerm => NULL()
       CALL AllocateVector( LocalEPerm, MaxEIndex - MinEIndex + 1, 'LoadMesh' )
       LocalEPerm = 0
       DO i=1,Mesh % NumberOfBulkElements
         LocalEPerm( ElementTags(i) - MinEIndex + 1 ) = i
       END DO
     ELSE
       MinEIndex = 1 
       MaxEIndex = Mesh % NumberOfBulkElements
     END IF


     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', iostat=IOSTAT )
     IF( iostat /= 0 ) THEN
       CALL Fatal('ReadBoundaryFile','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading boundary elements from file: '//TRIM(FileName),Level=10)
     END IF


     DO j=Mesh % NumberOfBulkElements+1, &
         Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

       Element => Mesh % Elements(j)
       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('ReadElementsFile','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read boundary element entry: '//TRIM(I2S(j)))
       END IF
       nread = read_ints(str,ivals,halo)
       
       tag = ivals(1)

       IF( halo ) THEN
         partn = ivals(2)
         ioffset = 1
       ELSE
         partn = 0
         ioffset = 0
       END IF

       bndry = ivals(ioffset+2)
       left = ivals(ioffset+3)
       right = ivals(ioffset+4)
       ElemType = ivals(ioffset+5)
       
       Element % ElementIndex = j
       Element % TYPE => GetElementType(ElemType)
       IF ( .NOT. ASSOCIATED(Element % TYPE) ) THEN
         CALL Fatal('ReadBoundaryFile','Element of type '//TRIM(I2S(ElemType))//'could not be associated!')
       END IF

       ElemNodes = Element % TYPE % NumberOfNodes
       Mesh % MaxElementNodes = MAX( Mesh % MaxElementNodes, ElemNodes )

       IF( partn == 0 ) THEN
         Element % PartIndex = mype
       ELSE
         Element % PartIndex = partn-1
       END IF

       CALL AllocateBoundaryInfo( Element ) 

       Element % BoundaryInfo % Constraint = bndry
       Element % BoundaryInfo % Left => NULL()
       Element % BoundaryInfo % Right => NULL()

       IF ( Left >= MinEIndex .AND. Left <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Left  = LocalEPerm(Left - MinEIndex + 1)
         END IF
       ELSE IF ( Left > 0 ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag, Left
         CALL Error( 'ReadBoundaryFile', Message )
         Left = 0
       END IF

       IF ( Right >= MinEIndex .AND. Right <= MaxEIndex ) THEN
         IF( ElementPermutation ) THEN
           Right = LocalEPerm(Right - MinEIndex + 1)
         END IF
       ELSE IF ( Right > 0 ) THEN
         WRITE( Message, * ) mype,'BOUNDARY PARENT out of range: ', Tag,Right
         CALL Error( 'ReadBoundaryFile', Message )
         Right = 0
       END IF

       IF ( Left >= 1 ) THEN
         Element % BoundaryInfo % Left => Mesh % Elements(left)
       END IF

       IF ( Right >= 1 ) THEN
         Element % BoundaryInfo % Right => Mesh % Elements(right)
       END IF

       n = Element % TYPE % NumberOfNodes
       CALL AllocateVector( Element % NodeIndexes, n )

       IF( nread < 5 + n + ioffset ) THEN
         CALL Fatal('ReadBoundaryFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF
       Element % NodeIndexes(1:n) = Ivals(6+ioffset:nread)
     END DO
     CLOSE( FileUnit )


     IF( ElementPermutation ) THEN
       DEALLOCATE( LocalEPerm ) 
     END IF

   END SUBROUTINE ReadBoundaryFile
   !------------------------------------------------------------------------------



   ! Make a permutation for the bulk and boundary element topology if 
   ! the nodes are permuted. This is always the case in parallel.
   ! The initial numbering is needed only when the nodes are loaded and 
   ! hence this is a local subroutine. 
   !----------------------------------------------------------------------
   SUBROUTINE PermuteNodeNumbering()

     TYPE(Element_t), POINTER :: Element

     IF( NodePermutation ) THEN
       CALL Info('LoadMesh','Performing node mapping',Level=6)

       MinNodeTag = MINVAL( NodeTags )
       MaxNodeTag = MAXVAL( NodeTags )

       CALL AllocateVector( LocalPerm, MaxNodeTag-MinNodeTag+1, 'LoadMesh' )
       LocalPerm = 0
       DO i=1,Mesh % NumberOfNodes
         LocalPerm(NodeTags(i) - MinNodeTag + 1) = i
       END DO

       DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements       
         Element => Mesh % Elements(i)
         n = Element % TYPE % NumberOfNodes

         DO j=1,n
           k = Element % NodeIndexes(j) 
           Element % NodeIndexes(j) = LocalPerm(k - MinNodeTag + 1)
         END DO
       END DO
     ELSE
       CALL Info('LoadMesh','Node mapping is continuous',Level=8)
     END IF

     ! Set the for now, if the case is truly parallel we'll have to revisit these
     ! when reading the parallel information. 
     Mesh % ParallelInfo % NumberOfIfDOFs = 0
     Mesh % ParallelInfo % GlobalDOFs => NodeTags

   END SUBROUTINE PermuteNodeNumbering


   ! Initialize some parallel structures once the non-nodal 
   ! element types are known. 
   ! Currently this is here mainly because the 
   ! Elemental and Nodal tags are local
   !-------------------------------------------------------
   SUBROUTINE InitParallelInfo()

     INTEGER, POINTER :: TmpGlobalDofs(:)

     ! These two have already been set, and if the case is serial
     ! case they can be as is.
     !Mesh % ParallelInfo % NumberOfIfDOFs = 0
     !Mesh % ParallelInfo % GlobalDOFs => NodeTags


     ! This also for serial runs ...
     DO i=1,Mesh % NumberOfBulkElements
       Mesh % Elements(i) % GElementIndex = ElementTags(i)
     END DO

     IF(.NOT. Parallel ) RETURN

     n = Mesh % NumberOfNodes + &
         Mesh % MaxEdgeDOFs * Mesh % NumberOFEdges + &
         Mesh % MaxFaceDOFs * Mesh % NumberOFFaces + &
         Mesh % MaxBDOFs    * Mesh % NumberOFBulkElements

     ALLOCATE( TmpGlobalDOFs(n) )
     TmpGlobalDOFs = 0
     TmpGlobalDOFs(1:Mesh % NumberOfNodes) = &
         Mesh % ParallelInfo % GlobalDOFs(1:Mesh % NumberOfNodes)
     DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs ) 
     Mesh % ParallelInfo % GlobalDofs => TmpGlobalDofs

     ALLOCATE(Mesh % ParallelInfo % NeighbourList(n), STAT=istat)
     IF (istat /= 0) CALL Fatal('LoadMesh', 'Unable to allocate NeighbourList array.')

     DO i=1,n
       NULLIFY( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
     END DO

     CALL AllocateVector( Mesh % ParallelInfo % INTERFACE, n, 'LoadMesh')
     Mesh % ParallelInfo % INTERFACE = .FALSE.       

   END SUBROUTINE InitParallelInfo


   ! Read the file that shows the shared nodes.
   !------------------------------------------------------------------------
   SUBROUTINE ReadSharedFile()

     INTEGER :: Ivals(64)
     INTEGER :: npart, tag, nread
     CHARACTER(256) :: str
     LOGICAL :: halo

     IF(.NOT. Parallel) RETURN

     FileName = BaseName(1:BaseNameLen)//&
       '/partitioning.'//TRIM(I2S(numprocs))//&
         '/part.'//TRIM(I2S(mype+1))//'.shared'

     OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT = iostat )
     IF( iostat /= 0 ) THEN
       CALL Fatal('LoadMesh','Could not open file: '//TRIM(Filename))
     ELSE
       CALL Info('LoadMesh','Reading nodes from file: '//TRIM(FileName),Level=10)
     END IF

     ! This loop could be made more effective, for example
     ! by reading tags and nparts to a temporal vector
     ! The operation using the str takes much more time.
     !-----------------------------------------------------
     DO i=1,SharedNodes          
       READ(FileUnit, '(a)', IOSTAT=iostat) str
       IF( iostat /= 0 ) THEN
         CALL Fatal('ReadElementsFile','Could not read shared nodes entry: '//TRIM(I2S(i)))
       END IF
       nread = read_ints(str,ivals,halo)

       tag = ivals(1)
       npart = ivals(2)       

       k = LocalPerm( tag-MinNodeTag+1 )
       Mesh % ParallelInfo % INTERFACE(k) = .TRUE.
       CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(k) % Neighbours,npart)

       IF( nread < 2 + npart ) THEN
         CALL Fatal('ReadSharedFile','Line '//TRIM(I2S(j))//' does not contain enough entries')
       END IF
       
       Mesh % ParallelInfo % NeighbourList(k) % Neighbours = ivals(3:nread) - 1

       ! this partition does not own the node
       IF ( ivals(3)-1 /= mype ) THEN
         Mesh % ParallelInfo % NumberOfIfDOFs = &
             Mesh % ParallelInfo % NumberOfIfDOFs + 1
       END IF
     END DO

     CLOSE( FileUnit )

   END SUBROUTINE ReadSharedFile

 END SUBROUTINE ElmerAsciiMesh



 !> An interface over potential mesh loading strategies. 
 !----------------------------------------------------------------- 
 SUBROUTINE LoadMeshStep( Step, PMesh, MeshNamePar, ThisPe, NumPEs,IsParallel ) 
   
   IMPLICIT NONE

   INTEGER :: Step
   CHARACTER(LEN=*), OPTIONAL :: MeshNamePar
   TYPE(Mesh_t), POINTER, OPTIONAL :: PMesh
   INTEGER, OPTIONAL :: ThisPe, NumPEs
   LOGICAL, OPTIONAL :: IsParallel

   ! Currently only one strategy to get the mesh is implemented 
   ! but there could be others.
   !
   ! This has not yet been tested in parallel and for sure
   ! it does not work for halo elements. 
   !-----------------------------------------------------------------
   CALL ElmerAsciiMesh( Step, PMesh, MeshNamePar, ThisPe, NumPEs, IsParallel ) 

 END SUBROUTINE LoadMeshStep

 !------------------------------------------------------------------------------
 ! Set the mesh dimension by studying the coordinate values.
 ! This could be less conservative also...
 !------------------------------------------------------------------------------    
 SUBROUTINE SetMeshDimension( Mesh )
   TYPE(Mesh_t), POINTER :: Mesh
   
   REAL(KIND=dp) :: x, y, z
   LOGICAL :: C(3)
   INTEGER :: i
   
   IF( Mesh % NumberOfNodes == 0 ) RETURN

   ! Compare value to some node, why not the 1st one
   x = Mesh % Nodes % x(1)
   y = Mesh % Nodes % y(1)
   z = Mesh % Nodes % z(1)
   
   C(1) = ANY( Mesh % Nodes % x /= x ) 
   C(2) = ANY( Mesh % Nodes % y /= y )  
   C(3) = ANY( Mesh % Nodes % z /= z )  

   ! This version is perhaps too liberal 
   Mesh % MeshDim = COUNT( C )
   Mesh % MaxDim = 0
   DO i=1,3
     IF( C(i) ) Mesh % MaxDim = i
   END DO
      
   CALL Info('SetMeshDimension','Dimension of mesh is: '//TRIM(I2S(Mesh % MeshDim)),Level=8)
   CALL Info('SetMeshDimension','Max dimension of mesh is: '//TRIM(I2S(Mesh % MaxDim)),Level=8)

 END SUBROUTINE SetMeshDimension

 
 !------------------------------------------------------------------------------
 !> Function to load mesh from disk.
 !------------------------------------------------------------------------------
 FUNCTION LoadMesh2( Model, MeshDirPar, MeshNamePar,&
     BoundariesOnly, NumProcs, MyPE, Def_Dofs, mySolver, &
     LoadOnly ) RESULT( Mesh )
   !------------------------------------------------------------------------------
   USE PElementMaps, ONLY : GetRefPElementNodes

   IMPLICIT NONE

   CHARACTER(LEN=*) :: MeshDirPar,MeshNamePar
   LOGICAL :: BoundariesOnly    
   INTEGER, OPTIONAL :: numprocs,mype,Def_Dofs(:,:), mySolver
   TYPE(Mesh_t),  POINTER :: Mesh
   TYPE(Model_t) :: Model
   LOGICAL, OPTIONAL :: LoadOnly 
   !------------------------------------------------------------------------------    
   INTEGER :: i,j,k,n
   INTEGER :: BaseNameLen, Save_Dim
   LOGICAL :: GotIt, Found, ForcePrep=.FALSE.
   CHARACTER(MAX_NAME_LEN) :: FileName
   TYPE(Element_t), POINTER :: Element
   TYPE(Matrix_t), POINTER :: Projector
   LOGICAL :: parallel, LoadNewMesh


   Mesh => Null()

   n = LEN_TRIM(MeshNamePar)
   DO WHILE (MeshNamePar(n:n)==CHAR(0).OR.MeshNamePar(n:n)==' ')
     n=n-1
   END DO
   IF(NumProcs<=1) THEN
     INQUIRE( FILE=MeshNamePar(1:n)//'/mesh.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Fatal('LoadMesh','Requested mesh > '//MeshNamePar(1:n)//' < does not exist!')
     END IF
   ELSE
     INQUIRE( FILE=MeshNamePar(1:n)//'/partitioning.'// & 
         TRIM(i2s(Numprocs))//'/part.1.header', EXIST=Found)
     IF(.NOT. Found ) THEN
       CALL Warn('LoadMesh','Requested mesh > '//MeshNamePar(1:n)//' < in partition '&
           //TRIM(I2S(Numprocs))//' does not exist!')
       RETURN
     END IF
   END IF

   CALL Info('LoadMesh','Starting',Level=8)

   Parallel = .FALSE.
   IF ( PRESENT(numprocs) .AND. PRESENT(mype) ) THEN
     IF ( numprocs > 1 ) Parallel = .TRUE.
   END IF

   Mesh => AllocateMesh()

   ! Get sizes of mesh structures for allocation
   !--------------------------------------------------------------------
   CALL LoadMeshStep( 1, Mesh, MeshNamePar, mype, numprocs, Parallel )

   ! Initialize and allocate mesh structures
   !---------------------------------------------------------------------
   IF( BoundariesOnly ) Mesh % NumberOfBulkElements = 0
   CALL InitializeMesh( Mesh )

   ! Get the (x,y,z) coordinates
   !--------------------------------------------------------------------------
   CALL LoadMeshStep( 2 )
   ! Permute and scale the coordinates.
   ! This also finds the mesh dimension. It is needed prior to getting the 
   ! elementtypes since wrong permutation or dimension may spoil that. 
   !-------------------------------------------------------------------
   CALL MapCoordinates()
   
   ! Get the bulk elements: element types, body index, topology
   !--------------------------------------------------------------------------
   CALL LoadMeshStep( 3 )

   ! Get the boundary elements: boundary types, boundary index, parents, topology
   !------------------------------------------------------------------------------
   CALL LoadMeshStep( 4 )

   ! Read elemental data - this is rarely used, parallel implementation lacking?
   !--------------------------------------------------------------------------
   i = LEN_TRIM(MeshNamePar)
   DO WHILE(MeshNamePar(i:i) == CHAR(0))
     i=i-1
   END DO
   BaseNameLen = i
   
   FileName = MeshNamePar(1:BaseNameLen)//'/mesh.elements.data'
   CALL ReadElementPropertyFile( FileName, Mesh )

   ! Read mesh.names - this could be saved by some mesh formats
   !--------------------------------------------------------------------------
   IF( ListGetLogical( Model % Simulation,'Use Mesh Names',Found ) ) THEN
     FileName = MeshNamePar(1:BaseNameLen)//'/mesh.names'
     CALL ReadTargetNames( Model, FileName )
   END IF


   ! Map bodies using Target Bodies and boundaries using Target Boundaries.
   ! This must be done before the element definitions are studied since
   ! then the pointer should be to the correct body index. 
   !------------------------------------------------------------------------
   CALL MapBodiesAndBCs()

   ! Read parallel mesh information: shared nodes
   !------------------------------------------------------------------
   CALL LoadMeshStep( 5 )

   ! Create the discontinuous mesh that accounts for the jumps in BCs
   ! This must be created after the whole mesh has been read in and 
   ! bodies and bcs have been mapped to full operation.
   ! To consider non-nodal elements it must be done before them.
   !--------------------------------------------------------------------
   CALL CreateDiscontMesh(Model,Mesh)

   ! Deallocate some stuff no longer needed
   !------------------------------------------------------------------
   CALL LoadMeshStep( 6 )

   CALL Info('LoadMesh','Loading mesh done',Level=8)

   ForcePrep = ListGetLogical( Model % Simulation,'Finalize Meshes Before Extrusion',Found)
   
   IF( PRESENT( LoadOnly ) ) THEN
     IF( LoadOnly ) THEN
       RETURN
     ELSE
       ForcePrep = .TRUE.
     END IF
   END IF

   ! Prepare the mesh for next steps.
   ! For example, create non-nodal mesh structures, periodic projectors etc. 
   IF( (ListCheckPresent( Model % Simulation,'Extruded Mesh Levels') .OR. &
       ListCheckPresent( Model % Simulation,'Extruded Mesh Layers')) .AND. (.NOT. ForcePrep) ) THEN
     CALL Info('LoadMesh','This mesh will be extruded, skipping finalization',Level=12)
     RETURN
   END IF

   CALL PrepareMesh(Model,Mesh,Parallel,Def_Dofs,mySolver)      
   CALL Info('LoadMesh','Preparing mesh done',Level=8)

   
 CONTAINS


   !------------------------------------------------------------------------------
   ! Map bodies and boundaries as prescirbed by the 'Target Bodies' and 
   ! 'Target Boundaries' keywords.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapBodiesAndBCs()

     TYPE(Element_t), POINTER :: Element
     INTEGER, ALLOCATABLE :: IndexMap(:), TmpIndexMap(:)
     INTEGER, POINTER :: Blist(:)
     INTEGER :: id,minid,maxid,body,bndry,DefaultTargetBC


     ! If "target bodies" is used map the bodies accordingly
     !------------------------------------------------------
     Found = .FALSE. 
     DO id=1,Model % NumberOfBodies
       IF( ListCheckPresent( Model % Bodies(id) % Values,'Target Bodies') ) THEN
         Found = .TRUE.
         EXIT
       END IF
     END DO

     IF( Found ) THEN
       CALL Info('LoadMesh','Remapping bodies',Level=8)      
       minid = HUGE( minid ) 
       maxid = -HUGE( maxid ) 
       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId
         minid = MIN( id, minid ) 
         maxid = MAX( id, maxid )
       END DO
       IF( minid > maxid ) THEN
         CALL Fatal('LoadMesh','Body indexes are screwed!')
       END IF
       CALL Info('LoadMesh','Minimum initial body index: '//TRIM(I2S(minid)),Level=6 )
       CALL Info('LoadMesh','Maximum initial body index: '//TRIM(I2S(maxid)),Level=6 )

       minid = MIN( 1, minid ) 
       maxid = MAX( Model % NumberOfBodies, maxid ) 
       ALLOCATE( IndexMap(minid:maxid) )
       IndexMap = 0

       DO id=1,Model % NumberOfBodies
         BList => ListGetIntegerArray( Model % Bodies(id) % Values, &
             'Target Bodies', GotIt ) 
         IF ( Gotit ) THEN
           DO k=1,SIZE(BList)
             body = Blist(k)
             IF( body > maxid .OR. body < minid ) THEN
#if 0
               CALL Warn('LoadMesh','Unused body entry in > Target Bodies <  : '&
                   //TRIM(I2S(body)) )              
#endif
             ELSE IF( IndexMap( body ) /= 0 ) THEN
               CALL Warn('LoadMesh','Multiple bodies have same > Target Bodies < entry : '&
                   //TRIM(I2S(body)))
             ELSE
               IndexMap( body ) = id 
             END IF
           END DO
         ELSE
           IF( IndexMap( id ) /= 0 ) THEN
             CALL Warn('LoadMesh','Unset body already set by > Target Boundaries < : '&
                 //TRIM(I2S(id)) )
           ELSE 
             IndexMap( id ) = id
           END IF
         END IF

       END DO

       IF( .FALSE. ) THEN
         PRINT *,'Body mapping'
         DO id=minid,maxid
           IF( IndexMap( id ) /= 0 ) PRINT *,id,' : ',IndexMap(id)
         END DO
       END IF

       DO i=1,Mesh % NumberOfBulkElements
         Element => Mesh % Elements(i)
         id = Element % BodyId
!        IF( IndexMap( id ) == 0 ) THEN
!          PRINT *,'Unmapped body: ',id
!          IndexMap(id) = id
!        END IF
         Element % BodyId = IndexMap( id ) 
       END DO

       DEALLOCATE( IndexMap )
     ELSE
       CALL Info('LoadMesh','Skipping remapping of bodies',Level=10)      
     END IF


     IF( Mesh % NumberOfBoundaryElements == 0 ) RETURN

     ! Target boundaries are usually given so this is not conditional
     !---------------------------------------------------------------
     CALL Info('LoadMesh','Remapping boundaries',Level=8)      
     minid = HUGE( minid ) 
     maxid = -HUGE( maxid ) 
     DO i=Mesh % NumberOfBulkElements+1,&
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)
       id = Element % BoundaryInfo % Constraint
       minid = MIN( id, minid ) 
       maxid = MAX( id, maxid )
     END DO


     CALL Info('LoadMesh','Minimum initial boundary index: '//TRIM(I2S(minid)),Level=6 )
     CALL Info('LoadMesh','Maximum initial boundary index: '//TRIM(I2S(maxid)),Level=6 )
     IF( minid > maxid ) THEN
       CALL Fatal('LoadMesh','Boundary indexes are screwed')
     END IF

     minid = MIN( minid, 1 ) 
     maxid = MAX( maxid, Model % NumberOfBCs ) 
     ALLOCATE( IndexMap(minid:maxid) )
     IndexMap = 0


     DO j=1,Model % NumberOfBoundaries
       id = ListGetInteger( Model % Boundaries(j) % Values, &
           'Boundary Condition',GotIt, minv=1, maxv=Model % NumberOFBCs )
       IF( id == 0 ) CYCLE
       bndry = Model % BoundaryId(j)
       IF( bndry > maxid ) THEN
         CALL Warn('LoadMesh','BoundaryId exceeds range')
       ELSE IF( bndry == 0 ) THEN
         CALL Warn('LoadMesh','BoundaryId is zero')
       ELSE
         IndexMap( bndry ) = id
       END IF
     END DO

     DefaultTargetBC = 0
     DO id=1,Model % NumberOfBCs
       IF(ListGetLogical( Model % BCs(id) % Values, &
           'Default Target', GotIt)) DefaultTargetBC = id       
       BList => ListGetIntegerArray( Model % BCs(id) % Values, &
           'Target Boundaries', GotIt )
       IF ( Gotit ) THEN
         DO k=1,SIZE(BList)
           bndry = Blist(k)
           IF( bndry > maxid ) THEN
#if 0
  in my opinion, this is quite usual ... Juha
             CALL Warn('LoadMesh','Unused BC entry in > Target Boundaries <  : '&
                 //TRIM(I2S(bndry)) )              
#endif
           ELSE IF( IndexMap( bndry ) /= 0 ) THEN
             CALL Warn('LoadMesh','Multiple BCs have same > Target Boundaries < entry : '&
                 //TRIM(I2S(bndry)) )
           ELSE 
             IndexMap( bndry ) = id 
           END IF
         END DO
       ELSE
         IF (ListCheckPresent(Model % BCs(id) % Values, 'Target Nodes') .OR. &
             ListCheckPresent(Model % BCs(id) % Values, 'Target Coordinates')) &
             CYCLE
         IF (IndexMap( id ) /= 0 .AND. id == DefaultTargetBC ) THEN ! DefaultTarget has been given
           CALL Warn('LoadMesh','Default Target is a Target Boundaries entry in > Boundary Condition < : '&
               //TRIM(I2S(IndexMap(id))) )
         END IF
         !
         !IF( IndexMap( id ) /= 0 .AND. id /= DefaultTargetBC ) THEN
         !  CALL Warn('LoadMesh','Unset BC already set by > Target Boundaries < : '&
         !      //TRIM(I2S(id)) )
         !ELSE 
         !  ! IndexMap( id ) = id
         !END IF
       END IF
     END DO

     IF( .FALSE. ) THEN
       PRINT *,'Boundary mapping'
       DO id=minid,maxid
         IF( IndexMap( id ) /= 0 ) PRINT *,id,' : ',IndexMap(id)
       END DO
     END IF

     IF( DefaultTargetBC /= 0 ) THEN
       CALL Info('LoadMesh','Default Target BC: '&
           //TRIM(I2S(DefaultTargetBC)),Level=8)
     END IF


     DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       n = Element % TYPE % NumberOfNodes
       bndry = Element % BoundaryInfo % Constraint 

       IF( bndry > maxid .OR. bndry < minid ) THEN
         CALL Warn('LoadMesh','Boundary index '//TRIM(I2S(bndry))&
             //' not in range: '//TRIM(I2S(minid))//','//TRIM(I2S(maxid)) )
       END IF

       IF( IndexMap( bndry ) < 0 ) THEN
         Element % BoundaryInfo % Constraint = 0
         CYCLE

       ELSE IF( IndexMap( bndry ) == 0 ) THEN
         IF( DefaultTargetBC /= 0 ) THEN
!          PRINT *,'Default boundary map: ',bndry,DefaultTargetBC
           IndexMap( bndry ) = DefaultTargetBC
         ELSE 
!          IF( bndry <= Model % NumberOfBCs ) THEN            
!            PRINT *,'Unmapped boundary: ',bndry
!          ELSE
!            PRINT *,'Unused boundary: ',bndry
!          END IF
           IndexMap( bndry ) = -1 
           Element % BoundaryInfo % Constraint = 0           
           CYCLE
         END IF
       END IF

       bndry = IndexMap( bndry ) 
       Element % BoundaryInfo % Constraint = bndry 

       IF( bndry <= Model % NumberOfBCs ) THEN
         Element % BodyId  = ListGetInteger( &
             Model % BCs(bndry) % Values, 'Body Id', Gotit, 1, Model % NumberOfBodies )
         Element % BoundaryInfo % OutBody = &
             ListGetInteger( Model % BCs(bndry) % Values, &
             'Normal Target Body', GotIt, maxv=Model % NumberOFBodies ) 
       END IF
     END DO

     DEALLOCATE( IndexMap ) 

   END SUBROUTINE MapBodiesAndBCs

   

   !------------------------------------------------------------------------------
   ! Map and scale coordinates, and increase the size of the coordinate
   ! vectors, if requested.
   !------------------------------------------------------------------------------    
   SUBROUTINE MapCoordinates()

     REAL(KIND=dp), POINTER CONTIG :: NodesX(:), NodesY(:), NodesZ(:)
     REAL(KIND=dp), POINTER :: Wrk(:,:)
     INTEGER, POINTER :: CoordMap(:)
     REAL(KIND=dp) :: CoordScale(3)
     INTEGER :: mesh_dim, model_dim
     
     ! Perform coordinate mapping
     !------------------------------------------------------------
     CoordMap => ListGetIntegerArray( Model % Simulation, &
         'Coordinate Mapping',GotIt )
     IF ( GotIt ) THEN
       CALL Info('LoadMesh','Performing coordinate mapping',Level=8)

       IF ( SIZE( CoordMap ) /= 3 ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'LoadMesh', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'LoadMesh', Message )
       END IF

       IF ( ALL( CoordMap(1:3) /= 1 ) .OR. ALL( CoordMap(1:3) /= 2 ) .OR. ALL( CoordMap(1:3) /= 3 ) ) THEN
         WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
         CALL Error( 'LoadMesh', Message )
         WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
         CALL Fatal( 'LoadMesh', Message )
       END IF

       IF( CoordMap(1) == 1 ) THEN
         NodesX => Mesh % Nodes % x
       ELSE IF( CoordMap(1) == 2 ) THEN
         NodesX => Mesh % Nodes % y
       ELSE
         NodesX => Mesh % Nodes % z
       END IF

       IF( CoordMap(2) == 1 ) THEN
         NodesY => Mesh % Nodes % x
       ELSE IF( CoordMap(2) == 2 ) THEN
         NodesY => Mesh % Nodes % y
       ELSE
         NodesY => Mesh % Nodes % z
       END IF

       IF( CoordMap(3) == 1 ) THEN
         NodesZ => Mesh % Nodes % x
       ELSE IF( CoordMap(3) == 2 ) THEN
         NodesZ => Mesh % Nodes % y
       ELSE
         NodesZ => Mesh % Nodes % z
       END IF

       Mesh % Nodes % x => NodesX
       Mesh % Nodes % y => NodesY
       Mesh % Nodes % z => NodesZ
     END IF

     ! Determine the mesh dimension 
     !----------------------------------------------------------------------------
     CALL SetMeshDimension( Mesh )
     
     mesh_dim = Mesh % MaxDim

     ! Scaling of coordinates
     !-----------------------------------------------------------------------------
     Wrk => ListGetConstRealArray( Model % Simulation,'Coordinate Scaling',GotIt )    
     IF( GotIt ) THEN            
       CoordScale = 1.0_dp
       DO i=1,mesh_dim
         j = MIN( i, SIZE(Wrk,1) )
         CoordScale(i) = Wrk(j,1)
       END DO
       WRITE(Message,'(A,3ES10.3)') 'Scaling coordinates:',CoordScale(1:3)
       CALL Info('LoadMesh',Message) 
       Mesh % Nodes % x = CoordScale(1) * Mesh % Nodes % x
       IF( mesh_dim > 1 ) Mesh % Nodes % y = CoordScale(2) * Mesh % Nodes % y
       IF( mesh_dim > 2 ) Mesh % Nodes % z = CoordScale(3) * Mesh % Nodes % z
     END IF

   END SUBROUTINE MapCoordinates

 !------------------------------------------------------------------------------
 END FUNCTION LoadMesh2
 !------------------------------------------------------------------------------


 !> Prepare a clean nodal mesh as it comes after being loaded from disk.
 !> Study the non-nodal elements (face, edge, DG, and p-elements)
 !> Create parallel info for the non-nodal elements
 !> Enlarge the coordinate vectors for p-elements.
 !> Generate static projector for periodic BCS.
 !-------------------------------------------------------------------
 SUBROUTINE PrepareMesh( Model, Mesh, Parallel, Def_Dofs, mySolver )

   TYPE(Model_t) :: Model
   TYPE(Mesh_t), POINTER :: Mesh
   LOGICAL :: Parallel
   INTEGER, OPTIONAL :: Def_Dofs(:,:), mySolver
   LOGICAL :: Found

   
   
   IF( Mesh % MaxDim == 0) THEN
     CALL SetMeshDimension( Mesh )
   END IF
   Model % DIMENSION = MAX( Model % DIMENSION, Mesh % MaxDim ) 
   
   CALL NonNodalElements()

   IF( Parallel ) THEN
     CALL ParallelNonNodalElements()
   END IF
     
   CALL EnlargeCoordinates( Mesh ) 
   
   IF( ListGetLogical( Model % Simulation,'Inspect Quadratic Mesh', Found ) ) THEN
     CALL InspectQuadraticMesh( Mesh ) 
   END IF
   
   IF( ListGetLogical( Model % Simulation,'Inspect Mesh',Found ) ) THEN
     CALL InspectMesh( Mesh ) 
   END IF

   IF(ListGetLogical( Model % Simulation, 'Parallel Reduce Element Max Sizes', Found ) ) THEN
     Mesh % MaxElementDOFs  = NINT( ParallelReduction( 1.0_dp*Mesh % MaxElementDOFs,2  ) )
     Mesh % MaxElementNodes = NINT( ParallelReduction( 1.0_dp*Mesh % MaxElementNodes,2 ) )
   END IF
   
   
 CONTAINS
     

   ! Check for the non-nodal element basis
   !--------------------------------------------------------
   SUBROUTINE NonNodalElements()

     INTEGER, POINTER :: EdgeDofs(:), FaceDofs(:)
     INTEGER :: i, j, k, l, s, n, DGIndex, body_id, body_id0, eq_id, solver_id, el_id, &
         mat_id
     LOGICAL :: NeedEdges, Found, FoundDef0, FoundDef, FoundEq, GotIt, MeshDeps, &
         FoundEqDefs, FoundSolverDefs(Model % NumberOfSolvers), &
         FirstOrderElements, InheritDG, Hit, Stat
     TYPE(Element_t), POINTER :: Element, Parent, pParent
     TYPE(Element_t) :: DummyElement
     TYPE(ValueList_t), POINTER :: Vlist
     INTEGER :: inDOFs(10,6)
     CHARACTER(MAX_NAME_LEN) :: ElementDef0, ElementDef
     
     
     EdgeDOFs => NULL()
     CALL AllocateVector( EdgeDOFs, Mesh % NumberOfBulkElements, 'LoadMesh' )
     FaceDOFs => NULL()
     CALL AllocateVector( FaceDOFs, Mesh % NumberOfBulkElements, 'LoadMesh' )     
    
     DGIndex = 0
     NeedEdges = .FALSE.

     InDofs = 0
     InDofs(:,1) = 1
     IF ( PRESENT(Def_Dofs) ) THEN
       inDofs = Def_Dofs
     ELSE
       DO s=1,Model % NumberOfSolvers
         DO i=1,6
           DO j=1,8
             inDofs(j,i) = MAX(Indofs(j,i),MAXVAL(Model % Solvers(s) % Def_Dofs(j,:,i)))
           END DO
         END DO
       END DO
     END IF

     ! P-basis only over 1st order elements:
     ! -------------------------------------
     FirstOrderElements = .TRUE.
     DO i=1,Mesh % NumberOfBulkElements
       IF (Mesh % Elements(i) % Type % BasisFunctionDegree>1) THEN
         FirstOrderElements = .FALSE.; EXIT
       END IF
     END DO

    !
    ! Check whether the "Element" definitions can depend on mesh
    ! -----------------------------------------------------------
    MeshDeps = .FALSE.; FoundEqDefs = .FALSE.;  FoundSolverDefs = .FALSE.

    !
    ! As a preliminary step, check if an element definition is given 
    ! an equation section. The more common way is give the element
    ! definition in a solver section.
    !
    DO eq_id=1,Model % NumberOFEquations
      Vlist => Model % Equations(eq_id) % Values
      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)
      FoundEqDefs = FoundEqDefs .OR. FoundDef0

      IF (FoundDef0) THEN
        !
        ! Check if the order of p-basis is defined by calling a special
        ! MATC function:
        !
        j = INDEX(ElementDef0,'p:')
        IF (j>0.AND. ElementDef0(j+2:j+2)=='%') MeshDeps = .TRUE.
      ELSE
        !
        ! Check if element definitions are given for each solver separately
        ! by using a special keyword construct and tag the corresponding
        ! entries in the list of the solvers. This was thought to serve
        ! the definition of bodywise p-orders, but it seems this doesn't
        ! work really. TO DO: REPAIR OR REMOVE
        ! 
        DO Solver_id=1,Model % NumberOfSolvers
          IF (PRESENT(mySolver)) THEN
            IF ( Solver_id /= mySolver ) CYCLE
          ELSE
            IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
          END IF

          ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)
          FoundSolverDefs(Solver_id) = FoundSolverDefs(solver_id) .OR. FoundDef

          j = INDEX(ElementDef,'p:')
          IF (j>0.AND. ElementDef0(j+2:j+2)=='%') MeshDeps = .TRUE.
        END DO
      END IF
    END DO

    !
    ! Tag solvers for which the element definition has been given in
    ! a solver section:
    !
    DO solver_id=1,Model % NumberOFSolvers
      Vlist => Model % Solvers(solver_id) % Values

      ElementDef0 = ListGetString(Vlist,'Element',FoundDef0)
      FoundSolverDefs(Solver_id) = FoundSolverDefs(solver_id) .OR. FoundDef0

      j = INDEX(ElementDef0,'p:')
      IF (j>0.AND. ElementDef0(j+2:j+2)=='%') meshdeps = .TRUE.
    END DO

    ! The basic case without the order of p-basis being defined by a MATC function:
    !
    IF (.NOT.MeshDeps) THEN
      ElementDef = ' '
      FoundDef0 = .FALSE.
      DO body_id=1,Model % NumberOfBodies
        ElementDef0 = ' '
        Vlist => Model % Bodies(body_id) % Values
        eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
        IF( FoundEq ) THEN
          Vlist => Model % Equations(eq_id) % Values
          IF(FoundEqDefs) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

          DO solver_id=1,Model % NumberOfSolvers

            IF(PRESENT(mySolver)) THEN
              IF ( Solver_id /= mySolver ) CYCLE
            ELSE
              IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
            END IF

            FoundDef = .FALSE.
            IF(FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)

            IF ( FoundDef ) THEN
              CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef, solver_id, body_id, Indofs )
            ELSE
              IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) &
                 ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

              CALL GetMaxDefs( Model, Mesh, DummyElement, ElementDef0, solver_id, body_id, Indofs )

              IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
            END IF
          END DO
        END IF
      END DO
    END IF

     ! non-nodal elements in bulk elements
     !------------------------------------------------------------
     body_id0 = -1; FoundDef=.FALSE.; FoundEq=.FALSE.
     ElementDef = ' '

     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       body_id = Element % BodyId
       n = Element % TYPE % NumberOfNodes
       
       ! Check the Solver specific element types
       IF( Meshdeps ) THEN
         IF ( body_id/=body_id0 ) THEN
           Vlist => Model % Bodies(body_id) % Values
           eq_id = ListGetInteger(Vlist,'Equation',FoundEq)
         END IF

         ElementDef0 = ' '
         IF( FoundEq ) THEN
           Vlist => Model % Equations(eq_id) % Values
           FoundDef0 = .FALSE.
           IF( FoundEqDefs.AND.body_id/=body_id0 ) ElementDef0 = ListGetString(Vlist,'Element',FoundDef0 )

           DO solver_id=1,Model % NumberOfSolvers
             IF(PRESENT(mySolver)) THEN
               IF ( Solver_id /= mySolver ) CYCLE
             ELSE
               IF (ListCheckPresent(Model % Solvers(Solver_id) % Values, 'Mesh')) CYCLE
             END IF

             FoundDef = .FALSE.
             IF (FoundSolverDefs(solver_id)) &
                ElementDef = ListGetString(Vlist,'Element{'//TRIM(i2s(solver_id))//'}',FoundDef)

             IF ( FoundDef ) THEN
               CALL GetMaxDefs( Model, Mesh, Element, ElementDef, solver_id, body_id, Indofs )
             ELSE
               IF(.NOT. FoundDef0.AND.FoundSolverDefs(solver_id)) &
                  ElementDef0 = ListGetString(Model % Solvers(solver_id) % Values,'Element',GotIt)

               CALL GetMaxDefs( Model, Mesh, Element, ElementDef0, solver_id, body_id, Indofs )

               IF(.NOT. FoundDef0.AND.FoundSolverDefs(Solver_id)) ElementDef0 = ' '
             END IF
           END DO
         END IF
         body_id0 = body_id
       END IF


       el_id = Element % TYPE % ElementCode / 100

       ! Apply the elementtypes
       IF ( inDOFs(el_id,1) /= 0 ) THEN
         Element % NDOFs = n
       ELSE
         Element % NDOFs = 0
       END IF

       EdgeDOFs(i) = MAX(0,inDOFs(el_id,2))
       FaceDOFs(i) = MAX(0,inDOFs(el_id,3))

       IF ( inDofs(el_id,4) == 0 ) THEN
         inDOFs(el_id,4) = n
       END IF
         
       NULLIFY( Element % DGIndexes )
       IF ( inDOFs(el_id,4) > 0 ) THEN
         CALL AllocateVector( Element % DGIndexes, inDOFs(el_id,4))
         DO j=1,inDOFs(el_id,4)
           DGIndex = DGIndex + 1
           Element % DGIndexes(j) = DGIndex
         END DO
       ELSE
         NULLIFY( Element % DGIndexes )
       END IF
       Element % DGDOFs = MAX(0,inDOFs(el_id,4))
       NeedEdges = NeedEdges .OR. ANY( inDOFs(el_id,2:4)>0 )
       
       ! Check if given element is a p element
       IF (FirstOrderElements .AND. inDOFs(el_id,6) > 0) THEN
         CALL AllocatePDefinitions(Element)
         NeedEdges = .TRUE.

         ! Calculate element bubble dofs and set element p
         Element % PDefs % P = inDOFs(el_id,6)
         IF ( inDOFs(el_id,5) > 0 ) THEN
           Element % BDOFs = inDOFs(el_id,5)
         ELSE
           Element % BDOFs = getBubbleDOFs(Element, Element % PDefs % P)
         END IF

         ! All elements in actual mesh are not edges
         Element % PDefs % pyramidQuadEdge = .FALSE.
         Element % PDefs % isEdge = .FALSE.

         ! If element is of type tetrahedron and is a p element, 
         ! do the Ainsworth & Coyle trick
         IF (Element % TYPE % ElementCode == 504) CALL ConvertToACTetra(Element)
         CALL GetRefPElementNodes( Element % Type,  Element % Type % NodeU, &
             Element % Type % NodeV, Element % Type % NodeW )
       ELSE 
         ! Clear P element definitions and set manual bubbles
         Element % PDefs => NULL()
         Element % BDOFs = MAX(0,inDOFs(el_id,5))
         ! WRITE (*,*) Element % BDOFs
       END IF

       Mesh % MaxElementNodes = MAX( &
           Mesh % MaxElementNodes,Element % TYPE % NumberOfNodes )
     END DO

     InheritDG = .FALSE.
     IF( dgindex > 0 ) THEN
       InheritDG = ListCheckPresentAnyMaterial( CurrentModel,'DG Parent Material')
     END IF
     
     ! non-nodal elements in boundary elements
     !------------------------------------------------------------    
     DO i = Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements 

       Element => Mesh % Elements(i)

       IF(.NOT. ASSOCIATED( Element ) ) THEN
         CALL Fatal('NonNodalElements','Element '//TRIM(I2S(i))//' not associated!')
       END IF

       IF(.NOT. ASSOCIATED( Element % TYPE ) ) THEN
         CALL Fatal('NonNodalElements','Type in Element '//TRIM(I2S(i))//' not associated!')
       END IF

       n = Element % TYPE % NumberOfNodes
       Element % NDOFs  = n
       el_id = ELement % TYPE % ElementCode / 100

       IF ( ASSOCIATED(Element % BoundaryInfo % Left) ) THEN
         IF( Element % BoundaryInfo % Left % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           Element % BDOFs = &
               EdgeDOFs(Element % BoundaryInfo % Left % ElementIndex)
         ELSE
           Element % BDOFs = FaceDOFs(Element % BoundaryInfo % Left % ElementIndex)
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       IF ( ASSOCIATED(Element % BoundaryInfo % Right) ) THEN
         IF ( Element % BoundaryInfo % Right % NDOFs == 0 ) THEN
           Element % NDOFs = 0
         END IF

         IF ( Element % TYPE % DIMENSION == 1 ) THEN
           Element % BDOFs = &
               EdgeDOFs(Element % BoundaryInfo % Right % ElementIndex)
         ELSE
           Element % BDOFs = FaceDOFs(Element % BoundaryInfo % Right % ElementIndex)
           Element % BDOFs = MAX(Element % BDOFs, MAX(0,InDOFs(el_id+6,5)))
         END IF
       END IF

       ! Optionally also set DG indexes for BCs
       ! It is easy for outside boundaries, but for internal boundaries
       ! we need a flag "DG Parent Material".
       IF( InheritDG ) THEN
         IF(.NOT. ASSOCIATED( Element % DGIndexes ) ) THEN
           ALLOCATE( Element % DGIndexes(n) )
           Element % DGIndexes = 0
         END IF
         
         Hit = .TRUE.
         k = 0
         DO l=1,2        
           IF(l==1) THEN
             Parent => Element % BoundaryInfo % Left
           ELSE
             Parent => Element % BoundaryInfo % Right
           END IF
           IF(.NOT. ASSOCIATED( Parent ) ) CYCLE
           k = k + 1
           pParent => Parent
           
           mat_id = ListGetInteger( CurrentModel % Bodies(Parent % BodyId) % Values,&
               'Material',Found )
           IF(mat_id > 0 ) THEN           
             VList => CurrentModel % Materials(mat_id) % Values
           END IF
           IF( ASSOCIATED(Vlist) ) THEN
             Hit = ListGetLogical(Vlist,'DG Parent Material',Found )
           END IF
           IF( Hit ) EXIT
         END DO
         
         IF( k == 0 ) THEN
           CALL Fatal('NonnodalElements','Cannot define DG indexes for BC!')
         ELSE IF( k == 1 ) THEN
           Parent => pParent        
         ELSE IF(.NOT. Hit ) THEN
           CALL Fatal('NonnodalElements','Cannot define DG indexes for internal BC!')       
         END IF
         
         DO l=1,n
           DO j=1, Parent % TYPE % NumberOfNodes
             IF( Element % NodeIndexes(l) == Parent % NodeIndexes(j) ) THEN
               Element % DGIndexes(l) = Parent % DGIndexes(j)
               EXIT
             END IF
           END DO
         END DO
       END IF
       
     END DO

     IF ( Mesh % MaxElementDOFs <= 0 ) Mesh % MaxElementDOFs = Mesh % MaxElementNodes 

     ! Override automated "NeedEdges" if requested by the user.
     !------------------------------------------------------------------------------------
     IF(PRESENT(mySolver)) THEN
       Stat = ListGetLogical(Model % Solvers(mySolver) % Values, 'Need Edges', Found)
       IF(Found) NeedEdges = Stat

       IF( ListGetLogical(Model % Solvers(mySolver) % Values, 'NeedEdges', Found) ) THEN
         IF(.NOT. NeedEdges) CALL Fatal('NonNodalElements','Use "Need Edges" instead of "NeedEdges"') 
       END IF
     END IF

     IF( Mesh % MeshDim == 2 ) THEN
       Stat = ListGetLogical(Model % Simulation, 'Need Edges 2D', Found)
       IF(Found) NeedEdges = Stat
     END IF

     IF( Mesh % MeshDim == 3 ) THEN
       Stat = ListGetLogical(Model % Simulation, 'Need Edges 3D', Found)
       IF(Found) NeedEdges = Stat
     END IF
     
     IF ( NeedEdges ) THEN
       CALL Info('NonNodalElements','Requested elements require creation of edges',Level=8)
       CALL SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs)
     END IF

     CALL SetMeshMaxDOFs(Mesh)

     IF( ASSOCIATED(EdgeDOFs) ) DEALLOCATE(EdgeDOFs )
     IF( ASSOCIATED(FaceDOFs) ) DEALLOCATE(FaceDOFs)

   END SUBROUTINE NonNodalElements


   ! When the parallel nodal neighbours have been found 
   ! perform numbering for face and edge elements as well.
   !-------------------------------------------------------------------    
   SUBROUTINE ParallelNonNodalElements()

     INTEGER :: i,n,mype     
     TYPE(Element_t), POINTER :: Element

     !IF(.NOT. Parallel ) RETURN

     n = SIZE( Mesh % ParallelInfo % NeighbourList )
     mype = ParEnv % Mype

     
     ! For unset neighbours just set the this partition to be the only owner
     DO i=1,n
       IF (.NOT.ASSOCIATED(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)) THEN
         CALL AllocateVector(Mesh % ParallelInfo % NeighbourList(i) % Neighbours,1)
         Mesh % ParallelInfo % NeighbourList(i) % Neighbours(1) = mype
       END IF
     END DO

     ! Create parallel numbering of faces
     CALL SParFaceNumbering(Mesh, .TRUE. )

     DO i=1,Mesh % NumberOfFaces
       Mesh % MinFaceDOFs = MIN(Mesh % MinFaceDOFs,Mesh % Faces(i) % BDOFs)
       Mesh % MaxFaceDOFs = MAX(Mesh % MaxFaceDOFs,Mesh % Faces(i) % BDOFs)
     END DO
     IF(Mesh % MinFaceDOFs > Mesh % MaxFaceDOFs) Mesh % MinFaceDOFs = Mesh % MaxFaceDOFs

     ! Create parallel numbering for edges
     CALL SParEdgeNumbering(Mesh, .TRUE.)

     DO i=1,Mesh % NumberOfEdges
       Mesh % MinEdgeDOFs = MIN(Mesh % MinEdgeDOFs,Mesh % Edges(i) % BDOFs)
       Mesh % MaxEdgeDOFs = MAX(Mesh % MaxEdgeDOFs,Mesh % Edges(i) % BDOFs)
     END DO
     IF(Mesh % MinEdgeDOFs > Mesh % MaxEdgeDOFs) Mesh % MinEdgeDOFs = Mesh % MaxEdgeDOFs

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)        
       Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
           Element % TYPE % NumberOfNodes + &
           Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
           Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
           Element % BDOFs, &
           Element % DGDOFs )
     END DO


   END SUBROUTINE ParallelNonNodalElements

   
 END SUBROUTINE PrepareMesh

 

 SUBROUTINE InspectMesh(Mesh)
   
   TYPE(Mesh_t), POINTER :: Mesh
   INTEGER :: i,j,mini,maxi
   INTEGER, POINTER :: Indexes(:)
   INTEGER, ALLOCATABLE :: ActiveCount(:)

   PRINT *,'Inspecting mesh for ranges and correctness'

   PRINT *,'No bulk elements:',Mesh % NumberOfBulkElements
   PRINT *,'No boundary elements:',Mesh % NumberOfBoundaryElements
   PRINT *,'No nodes:',Mesh % NumberOfNodes

   PRINT *,'Range:'
   PRINT *,'X:',MINVAL( Mesh % Nodes % x ), MAXVAL( Mesh % Nodes % x )
   PRINT *,'Y:',MINVAL( Mesh % Nodes % y ), MAXVAL( Mesh % Nodes % y )
   PRINT *,'Z:',MINVAL( Mesh % Nodes % z ), MAXVAL( Mesh % Nodes % z )

   ALLOCATE( ActiveCount( Mesh % NumberOfNodes ) )

   mini = HUGE(mini)
   maxi = 0
   ActiveCount = 0
   DO i=1,Mesh % NumberOfBulkElements
     Indexes => Mesh % Elements(i) % NodeIndexes
     mini = MIN(mini, MINVAL( Indexes ) )
     maxi = MAX(maxi, MAXVAL( Indexes ) )
     ActiveCount(Indexes) = ActiveCount(Indexes) + 1
   END DO
   PRINT *,'Bulk index range: ',mini,maxi
   PRINT *,'Bulk nodes:',COUNT(ActiveCount > 0 )
   PRINT *,'Bulk index count: ',MINVAL(ActiveCount),MAXVAL(ActiveCount)

   mini = HUGE(mini)
   maxi = 0
   ActiveCount = 0
   DO i=Mesh % NumberOfBulkElements+1, &
       Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements
     Indexes => Mesh % Elements(i) % NodeIndexes
     mini = MIN(mini, MINVAL( Indexes ) )
     maxi = MAX(maxi, MAXVAL( Indexes ) )
     ActiveCount(Indexes) = ActiveCount(Indexes) + 1
   END DO
   PRINT *,'Boundary index range: ',mini,maxi
   PRINT *,'Boundary nodes: ',COUNT(ActiveCount > 0)
   PRINT *,'Boundary index count: ',MINVAL(ActiveCount),MAXVAL(ActiveCount)

   DEALLOCATE( ActiveCount )

   PRINT *,'Done inspecting mesh'

 END SUBROUTINE InspectMesh



!------------------------------------------------------------------------------
  SUBROUTINE SetMeshEdgeFaceDOFs(Mesh,EdgeDOFs,FaceDOFs,inDOFs,NeedEdges)
!------------------------------------------------------------------------------
    INTEGER, OPTIONAL :: EdgeDOFs(:), FaceDOFs(:)
    TYPE(Mesh_t) :: Mesh
    INTEGER, OPTIONAL :: indofs(:,:)
    LOGICAL, OPTIONAL :: NeedEdges
!------------------------------------------------------------------------------
    INTEGER :: i,j,el_id
    TYPE(Element_t), POINTER :: Element, Edge, Face
    LOGICAL :: AssignEdges
!------------------------------------------------------------------------------

    CALL FindMeshEdges(Mesh)

    AssignEdges = .FALSE.
    IF (PRESENT(NeedEdges)) AssignEdges = NeedEdges
    
    ! Set edge and face polynomial degree and degrees of freedom for
    ! all elements
    DO i=1,Mesh % NumberOFBulkElements
       Element => Mesh % Elements(i)

       ! Iterate each edge of element
       DO j = 1,Element % TYPE % NumberOfEdges
          Edge => Mesh % Edges( Element % EdgeIndexes(j) ) 
          
          ! Set attributes of p element edges
          IF ( ASSOCIATED(Element % PDefs) ) THEN   
             ! Set edge polynomial degree and dofs
             Edge % PDefs % P = MAX( Element % PDefs % P, Edge % PDefs % P)
             Edge % BDOFs = MAX(Edge % BDOFs, Edge % PDefs % P - 1)
             Edge % PDefs % isEdge = .TRUE.
             ! Get gauss points for edge. If no dofs 2 gauss points are 
             ! still needed for integration of linear equation!
             Edge % PDefs % GaussPoints = (Edge % BDOFs+2)**Edge % TYPE % DIMENSION  

             IF (ASSOCIATED(Edge % BoundaryInfo % Left) ) THEN
               CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Left, Mesh)
             ELSE
               CALL AssignLocalNumber(Edge, Edge % BoundaryInfo % Right, Mesh)
             END IF
             
          ! Other element types, which need edge dofs
          ELSE IF(PRESENT(EdgeDOFs)) THEN
            Edge % BDOFs = MAX(EdgeDOFs(i), Edge % BDOFs)
          ELSE
            Edge % BDOFs = Max(1, Edge % BDOFs)
          END IF

          ! Get maximum dof for edges
          Mesh % MinEdgeDOFs = MIN(Edge % BDOFs, Mesh % MinEdgeDOFs)
          Mesh % MaxEdgeDOFs = MAX(Edge % BDOFs, Mesh % MaxEdgeDOFs)
       END DO
       IF ( Mesh % MinEdgeDOFs > Mesh % MaxEdgeDOFs ) Mesh % MinEdgeDOFs = MEsh % MaxEdgeDOFs

       ! Iterate each face of element
       DO j=1,Element % TYPE % NumberOfFaces
          Face => Mesh % Faces( Element % FaceIndexes(j) )

          ! Set attributes of p element faces
          IF ( ASSOCIATED(Element % PDefs) ) THEN
             ! Set face polynomial degree and dofs
             Face % PDefs % P = MAX(Element % PDefs % P, Face % PDefs % P)
             ! Get number of face dofs
             Face % BDOFs = MAX( Face % BDOFs, getFaceDOFs(Element, Face % PDefs % P, j) )
             Face % PDefs % isEdge = .TRUE.
             Face % PDefs % GaussPoints = getNumberOfGaussPointsFace( Face, Mesh )
             IF (ASSOCIATED(Face % BoundaryInfo % Left) ) THEN
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Left, Mesh)
             ELSE
               CALL AssignLocalNumber(Face, Face % BoundaryInfo % Right, Mesh)
             END IF
          ELSE IF (PRESENT(FaceDOFs)) THEN
             el_id = face % TYPE % ElementCode / 100
             Face % BDOFs = MAX(FaceDOFs(i), Face % BDOFs)
             IF ( PRESENT(inDOFs) ) Face % BDOFs = MAX(Face % BDOFs, InDOFs(el_id+6,5))
          END IF
             
          ! Get maximum dof for faces
          Mesh % MinFaceDOFs = MIN(Face % BDOFs, Mesh % MinFaceDOFs)
          Mesh % MaxFaceDOFs = MAX(Face % BDOFs, Mesh % MaxFaceDOFs)
       END DO
    END DO
    IF ( Mesh % MinFaceDOFs > Mesh % MaxFaceDOFs ) Mesh % MinFaceDOFs = MEsh % MaxFaceDOFs

    ! Set local edges for boundary elements
    DO i=Mesh % NumberOfBulkElements + 1, &
         Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Element => Mesh % Elements(i)

       ! Here set local number and copy attributes to this boundary element for left parent.
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          ! Local edges are only assigned for p elements
          IF (ASSOCIATED(Element % BoundaryInfo % Left % PDefs)) THEN
            CALL AllocatePDefinitions(Element)
            Element % PDefs % isEdge = .TRUE.
            CALL AssignLocalNumber(Element, Element % BoundaryInfo % Left, Mesh)
            ! CYCLE
          END IF
       END IF

       ! Here set local number and copy attributes to this boundary element for right parent
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          ! Local edges are only assigned for p elements
          IF (ASSOCIATED(Element % BoundaryInfo % Right % PDefs)) THEN
             CALL AllocatePDefinitions(Element)
             Element % PDefs % isEdge = .TRUE.
             CALL AssignLocalNumber(Element, Element % BoundaryInfo % Right, Mesh)
          END IF
       END IF

       IF (AssignEdges) THEN
         IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
           CALL AssignLocalNumber(Element,Element % BoundaryInfo % Left, Mesh, NoPE=.TRUE.)
         END IF
         IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
           CALL AssignLocalNumber(Element,Element % BoundaryInfo % Right, Mesh, NoPE=.TRUE.)
         END IF
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE SetMeshEdgeFaceDofs
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
 SUBROUTINE SetMeshMaxDOFs(Mesh)
!------------------------------------------------------------------------------
   TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
   TYPE(Element_t), POINTER :: Element
   INTEGER :: i,j,n

   ! Set gauss points for each p element
   DO i=1,Mesh % NumberOfBulkElements
     Element => Mesh % Elements(i)

     IF ( ASSOCIATED(Element % PDefs) ) THEN
       Element % PDefs % GaussPoints = getNumberOfGaussPoints( Element, Mesh )
     END IF

     ! Set max element dofs here (because element size may have changed
     ! when edges and faces have been set). This is the absolute worst case.
     ! Element which has MaxElementDOFs may not even be present as a 
     ! real element
     Mesh % MaxElementDOFs = MAX( Mesh % MaxElementDOFs, &
          Element % TYPE % NumberOfNodes + &
          Element % TYPE % NumberOfEdges * Mesh % MaxEdgeDOFs + &
          Element % TYPE % NumberOfFaces * Mesh % MaxFaceDOFs + &
          Element % BDOFs, &
          Element % DGDOFs )

     Mesh % MaxBDOFs = MAX( Element % BDOFs, Mesh % MaxBDOFs )
   END DO

   DO i=1,Mesh % NumberOFBulkElements
     Element => Mesh % Elements(i)
     IF ( Element % BDOFs > 0 ) THEN
       ALLOCATE( Element % BubbleIndexes(Element % BDOFs) )
       DO j=1,Element % BDOFs
         Element % BubbleIndexes(j) = Mesh % MaxBDOFs*(i-1)+j
       END DO
     END IF
   END DO
!------------------------------------------------------------------------------
 END SUBROUTINE SetMeshMaxDOFs
!------------------------------------------------------------------------------
 
 SUBROUTINE ReadTargetNames(Model,Filename)
     CHARACTER(LEN=*) :: FileName
     TYPE(Model_t) :: Model
!------------------------------------------------------------------------------
   INTEGER, PARAMETER :: FileUnit = 10
   INTEGER, PARAMETER :: A=ICHAR('A'),Z=ICHAR('Z'),U2L=ICHAR('a')-ICHAR('A')
   INTEGER :: i,j,k,iostat,i1,i2,i3,n
   INTEGER :: ivals(256)
   CHARACTER(LEN=1024) :: str, name0, name1
   TYPE(ValueList_t), POINTER :: Vlist
   LOGICAL :: Found, AlreadySet

   OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', IOSTAT=iostat )
   IF( iostat /= 0 ) THEN
     CALL Fatal('ReadTargetNames','Requested the use of entity names but this file does not exits: '//TRIM(FileName))
   END IF
   
   CALL Info('ReadTargetNames','Reading names info from file: '//TRIM(FileName))

   DO WHILE( .TRUE. ) 
     READ(FileUnit,'(A)',IOSTAT=iostat) str
     IF( iostat /= 0 ) EXIT
     i = INDEX( str,'$')     
     j = INDEX( str,'=')
     IF( i == 0 .OR. j == 0 ) CYCLE

     i = i + 1
     DO WHILE(i<=LEN_TRIM(str) .AND. str(i:i)==' ')
       i = i + 1
     END DO     
     
     i1 = i
     i2 = j-1
     i3 = j+1

     ! Move to lowercase since the "name" in sif file is also
     ! always in lowercase. 
     DO i=i1,i2
       j = i+1-i1
       k = ICHAR(str(i:i))
       IF ( k >= A .AND. k<= Z ) THEN
         name0(j:j) = CHAR(k+U2L)
       ELSE
         name0(j:j) = str(i:i)
       END IF
     END DO

     n = str2ints( str(i3:),ivals )
     IF( n == 0 ) THEN
       CALL Fatal('ReadTargetNames','Could not find arguments for: '//str(i1:i2))
     END IF

     AlreadySet = .FALSE.

     DO i=1,Model % NumberOfBCs
       Vlist => Model % BCs(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches BC '//TRIM(I2S(i))
         IF( AlreadySet ) THEN
           CALL Fatal('ReadTargetNames','Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Boundaries') ) THEN
           CALL Info('ReadTargetNames','> Target Boundaries < already defined for BC '&
               //TRIM(I2S(i)))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Boundaries',n,ivals(1:n))
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO

     DO i=1,Model % NumberOfBodies
       Vlist => Model % Bodies(i) % Values
       name1 = ListGetString( Vlist,'Name',Found )
       IF(.NOT. Found ) CYCLE
       IF( name0(1:i2-i1+1) == TRIM(name1) ) THEN
!        PRINT *,'Name > '//TRIM(name1)//' < matches body '//TRIM(I2S(i))
         IF( AlreadySet ) THEN
           CALL Fatal('ReadTargetNames','Mapping of name is not unique: '//TRIM(name1) )
         ELSE IF( ListCheckPresent( Vlist,'Target Bodies') ) THEN
           CALL Info('ReadTargetNames','> Target Bodies < already defined for Body '&
               //TRIM(I2S(i)))
         ELSE
           CALL ListAddIntegerArray( Vlist,'Target Bodies',n,ivals(1:n))
           AlreadySet = .TRUE.
         END IF
       END IF
     END DO
     
     IF(.NOT. AlreadySet ) THEN
       CALL Warn('ReadTargetNames','Could not map name to Body nor BC: '//name0(1:i2-i1+1) )
     END IF

   END DO

   CLOSE(FileUnit)
   
 END SUBROUTINE ReadTargetNames


!------------------------------------------------------------------------------
!> This subroutine reads elementwise input data from the file mesh.elements.data 
!> and inserts the data into the structured data variable 
!> Mesh % Elements(element_id) % PropertyData. The contents of the file should
!> be arranged as
!> 
!> element: element_id_1
!> data_set_name_1: a_1 a_2 ... a_n
!> data_set_name_2: b_1 b_2 ... b_m
!> data_set_name_3: ...
!> end
!> element: ...
!> ...
!> end
!------------------------------------------------------------------------------
  SUBROUTINE ReadElementPropertyFile(FileName,Mesh)
!------------------------------------------------------------------------------
     CHARACTER(LEN=*) :: FileName
     TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    INTEGER, PARAMETER :: MAXLEN=1024
    CHARACTER(LEN=:), ALLOCATABLE :: str
    INTEGER :: i,j,n
    INTEGER, PARAMETER :: FileUnit = 10
    REAL(KIND=dp) :: x
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementData_t), POINTER :: PD,PD1
!------------------------------------------------------------------------------
    ALLOCATE(CHARACTER(MAX_STRING_LEN)::str)

    OPEN( Unit=FileUnit, File=FileName, STATUS='OLD', ERR=10 )

    DO WHILE( ReadAndTrim(FileUnit,str) )
      READ( str(9:),*) i
      IF ( i < 0 .OR. i > Mesh % NumberOFBulkElements ) THEN
        CALL Fatal( 'ReadElementProperties', 'Element id out of range.' )
      END IF

      IF ( SEQL( str, 'element:') ) THEN
        Element => Mesh % Elements(i)
        PD => Element % PropertyData

        DO WHILE(ReadAndTrim(FileUnit,str))
          IF ( str == 'end' ) EXIT

          i = INDEX(str, ':')
          IF ( i<=0 ) CYCLE

          IF ( .NOT.ASSOCIATED(PD)  ) THEN
            ALLOCATE( Element % PropertyData )
            PD => Element % PropertyData
            PD % Name = TRIM(str(1:i-1))
          ELSE
            DO WHILE(ASSOCIATED(PD))
              IF ( PD % Name==TRIM(str(1:i-1)) ) EXIT
              PD1 => PD
              PD => PD % Next
            END DO
            
            IF (.NOT. ASSOCIATED(PD) ) THEN
              ALLOCATE(PD1 % Next)
              PD => PD1 % Next
              PD % Name = TRIM(str(1:i-1))
            END IF
          END IF

          j = i+1
          n = 0
          DO WHILE(j<=LEN_TRIM(str))
            READ( str(j:), *, END=20,ERR=20 ) x
            n = n + 1
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
              j = j + 1
            END DO
            DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
              j = j + 1
            END DO
          END DO
20        CONTINUE
          IF ( n>0 ) THEN
            ALLOCATE(PD % Values(n))
            j = i+1
            n = 1
            DO WHILE(j<=LEN_TRIM(str))
              READ( str(j:), *, END=30,ERR=30 ) PD % Values(n)
              n = n + 1
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)==' ')
                j = j + 1
              END DO
              DO WHILE(j<=LEN_TRIM(str) .AND. str(j:j)/=' ')
                j = j + 1
              END DO
            END DO
30          CONTINUE
          END IF
        END DO
      END IF
    END DO

    CLOSE(FileUnit)

10  CONTINUE

!------------------------------------------------------------------------------
  END SUBROUTINE ReadElementPropertyFile
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE MeshStabParams( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    TYPE(Solver_t), POINTER :: Solver
    INTEGER :: i,n, istat
    LOGICAL :: stat, Stabilize, UseLongEdge
    TYPE(Nodes_t) :: Nodes
    TYPE(Element_t), POINTER :: Element
!------------------------------------------------------------------------------

    CALL Info('MeshStabParams','Computing stabilization parameters',Level=7)
    CALL ResetTimer('MeshStabParams')

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal('MeshStabParams','Mesh not associated')
    END IF
    
    IF ( Mesh % NumberOfNodes <= 0 ) RETURN

    Stabilize = .FALSE.
    
    DO i=1,CurrentModel % NumberOfSolvers
      Solver => CurrentModel % Solvers(i)
      IF ( ASSOCIATED( Mesh, Solver % Mesh ) ) THEN
        Stabilize = Stabilize .OR. &
            ListGetLogical( Solver % Values, 'Stabilize', Stat )
        Stabilize = Stabilize .OR. &
            ListGetString( Solver % Values,  &
            'Stabilization Method', Stat )=='vms'
        Stabilize = Stabilize .OR. &
            ListGetString( Solver % Values,  &
            'Stabilization Method', Stat )=='stabilized'
      END IF
    END DO

    Mesh % Stabilize = Stabilize 
    
    IF( ListGetLogical(CurrentModel % Simulation, &
        "Skip Mesh Stabilization",Stat) ) RETURN
    
    !IF( .NOT. Stabilize ) THEN
    !  CALL Info('MeshStabParams','No need to compute stabilization parameters',Level=10)      
    !  RETURN      
    !END IF
    
    CALL AllocateVector( Nodes % x, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % y, Mesh % MaxElementNodes )
    CALL AllocateVector( Nodes % z, Mesh % MaxElementNodes )

    UseLongEdge = ListGetLogical(CurrentModel % Simulation, &
         "Stabilization Use Longest Element Edge",Stat)

    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)
       n = Element % TYPE % NumberOfNodes
       Nodes % x(1:n) = Mesh % Nodes % x(Element % NodeIndexes)
       Nodes % y(1:n) = Mesh % Nodes % y(Element % NodeIndexes)
       Nodes % z(1:n) = Mesh % Nodes % z(Element % NodeIndexes)
       IF ( Mesh % Stabilize ) THEN
          CALL StabParam( Element, Nodes,n, &
              Element % StabilizationMK, Element % hK, UseLongEdge=UseLongEdge)
       ELSE
          Element % hK = ElementDiameter( Element, Nodes, UseLongEdge=UseLongEdge)
       END IF
    END DO
 
    DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )

    CALL CheckTimer('MeshStabParams',Level=7,Delete=.TRUE.)
!----------------------------------------------------------------------------
  END SUBROUTINE MeshStabParams
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> The quadratic mesh should be such that the center nodes lie roughly between
!> the corner nodes. This routine checks that this is actually the case.
!> The intended use for the routine is different kind of mesh related debugging.
!------------------------------------------------------------------------------
  SUBROUTINE InspectQuadraticMesh( Mesh, EnforceToCenter ) 
    
    TYPE(Mesh_t), TARGET :: Mesh
    LOGICAL, OPTIONAL :: EnforceToCenter

    LOGICAL :: Enforce
    INTEGER :: i,n,k,k1,k2,k3,ElemCode,ElemFamily,ElemDegree,ErrCount,TotCount
    REAL(KIND=dp) :: Center(3),Ref(3),Dist,Length
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: CenterMap(:,:)
    INTEGER, TARGET  :: TriangleCenterMap(3,3), QuadCenterMap(4,3), &
        TetraCenterMap(6,3), BrickCenterMap(12,3), WedgeCenterMap(9,3), PyramidCenterMap(8,3) 
    
    CALL Info('InspectQuadraticMesh','Inspecting quadratic mesh for outliers')
    CALL Info('InspectQuadraticMesh','Number of nodes: '//TRIM(I2S(Mesh % NumberOfNodes)),Level=8)
    CALL Info('InspectQuadraticMesh','Number of bulk elements: '&
        //TRIM(I2S(Mesh % NumberOfBulkElements)),Level=8)
    CALL Info('InspectQuadraticMesh','Number of boundary elements: '&
        //TRIM(I2S(Mesh % NumberOfBoundaryElements)),Level=8)


    IF( PRESENT( EnforceToCenter ) ) THEN
      Enforce = EnforceToCenter
    ELSE
      Enforce = .FALSE.
    END IF

    TriangleCenterMap(1,:) = [ 1, 2, 4]
    TriangleCenterMap(2,:) = [ 2, 3, 5]
    TriangleCenterMap(3,:) = [ 3, 1, 6]
    
    QuadCenterMap(1,:) = [ 1, 2, 5]
    QuadCenterMap(2,:) = [ 2, 3, 6]
    QuadCenterMap(3,:) = [ 3, 4, 7]
    QuadCenterMap(4,:) = [ 4, 1, 8]
    
    TetraCenterMap(1,:) = [ 1, 2, 5]
    TetraCenterMap(2,:) = [ 2, 3, 6]
    TetraCenterMap(3,:) = [ 3, 1, 7]
    TetraCenterMap(4,:) = [ 1, 4, 8]
    TetraCenterMap(5,:) = [ 2, 4, 9]
    TetraCenterMap(6,:) = [ 3, 4, 10]

    BrickCenterMap(1,:) = [ 1, 2,  9 ]
    BrickCenterMap(2,:) = [ 2, 3,  10 ]
    BrickCenterMap(3,:) = [ 3, 4,  11 ]
    BrickCenterMap(4,:) = [ 4, 1,  12 ]
    BrickCenterMap(5,:) = [ 1, 5,  13 ]
    BrickCenterMap(6,:) = [ 2, 6,  14 ]
    BrickCenterMap(7,:) = [ 3, 7,  15 ]
    BrickCenterMap(8,:) = [ 4, 8,  16 ]
    BrickCenterMap(9,:) = [ 5, 6,  17 ]
    BrickCenterMap(10,:) = [ 6, 7, 18 ]
    BrickCenterMap(11,:) = [ 7, 8, 19 ]
    BrickCenterMap(12,:) = [ 8, 5, 20 ]
    
    WedgeCenterMap(1,:) = [ 1, 2, 7 ]
    WedgeCenterMap(2,:) = [ 2, 3, 8 ]
    WedgeCenterMap(3,:) = [ 3, 1, 9 ]
    WedgeCenterMap(4,:) = [ 4, 5, 10 ]
    WedgeCenterMap(5,:) = [ 5, 6, 11 ]
    WedgeCenterMap(6,:) = [ 6, 4, 12 ]
    WedgeCenterMap(7,:) = [ 1, 4, 13 ]
    WedgeCenterMap(8,:) = [ 2, 5, 14 ]
    WedgeCenterMap(9,:) = [ 3, 6, 15 ]
    
    PyramidCenterMap(1,:) = [ 1,2,6 ]
    PyramidCenterMap(2,:) = [ 2,3,7 ]
    PyramidCenterMap(3,:) = [ 3,4,8 ]
    PyramidCenterMap(4,:) = [ 4,1,9 ]
    PyramidCenterMap(5,:) = [ 1,5,10 ]
    PyramidCenterMap(6,:) = [ 2,5,11 ]
    PyramidCenterMap(7,:) = [ 3,5,12 ]
    PyramidCenterMap(8,:) = [ 4,5,13 ]
    
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z
    
    !   Loop over elements:
    !   -------------------
    ErrCount = 0
    TotCount = 0

    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      Element => Mesh % Elements(i)

      ElemCode = Element % TYPE % ElementCode 
      ElemFamily = ElemCode / 100
      ElemDegree = Element % TYPE % BasisFunctionDegree
      
      ! Only check quadratic elements!
      IF( ElemDegree /= 2 ) CYCLE
      
      SELECT CASE( ElemFamily ) 

      CASE(3)
        n = 3
        CenterMap => TriangleCenterMap
        
      CASE(4)
        n = 4
        CenterMap => QuadCenterMap
        
      CASE(5)
        n = 6
        CenterMap => TetraCenterMap
        
      CASE(6)
        n = 8
        CenterMap => PyramidCenterMap
        
      CASE(7)
        n = 9
        CenterMap => WedgeCenterMap
        
      CASE(8)
        n = 12
        CenterMap => BrickCenterMap
        
      CASE DEFAULT
        CALL Fatal('FindMeshEdges','Element type '//TRIM(I2S(ElemCode))//' not implemented!')

      END SELECT
      
      !      Loop over every edge of every element:
      !      --------------------------------------
       DO k=1,n
         k1 = Element % NodeIndexes( CenterMap(k,1) )
         k2 = Element % NodeIndexes( CenterMap(k,2) )
         k3 = Element % NodeIndexes( CenterMap(k,3) )
         
         Center(1) = ( x(k1) + x(k2) ) / 2.0_dp
         Center(2) = ( y(k1) + y(k2) ) / 2.0_dp
         Center(3) = ( z(k1) + z(k2) ) / 2.0_dp

         Ref(1) = x(k3)
         Ref(2) = y(k3) 
         Ref(3) = z(k3)

         Length = SQRT( (x(k1) - x(k2))**2.0 + (y(k1) - y(k2))**2.0 + (z(k1) - z(k2))**2.0 )
         Dist = SQRT( SUM( (Center - Ref)**2.0 ) )

         TotCount = TotCount + 1
         IF( Dist > 0.01 * Length ) THEN
           ErrCount = ErrCount + 1
           PRINT *,'Center Displacement:',i,ElemCode,n,k,Dist/Length
         END IF

         IF( Enforce ) THEN
           x(k3) = Center(1)
           y(k3) = Center(2)
           z(k3) = Center(3)
         END IF

       END DO
     END DO
         
     IF( TotCount > 0 ) THEN
       CALL Info('InspectQuadraticMesh','Number of outlier nodes is '&
           //TRIM(I2S(ErrCount))//' out of '//TRIM(I2S(TotCount)),Level=6)
     ELSE
       CALL Info('InspectQuadraticMesh','No quadratic elements to inspect',Level=8)
     END IF

  END SUBROUTINE InspectQuadraticMesh




!------------------------------------------------------------------------------
  FUNCTION Find_Face(Mesh,Parent,Element) RESULT(ptr)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Ptr
    TYPE(Mesh_t) :: Mesh
    TYPE(Element_t) :: Parent, Element

    INTEGER :: i,j,k,n

    Ptr => NULL()
    DO i=1,Parent % TYPE % NumberOfFaces
      Ptr => Mesh % Faces(Parent % FaceIndexes(i))
      n=0
      DO j=1,Ptr % TYPE % NumberOfNodes
        DO k=1,Element % TYPE % NumberOfNodes
          IF (Ptr % NodeIndexes(j) == Element % NodeIndexes(k)) n=n+1
        END DO
      END DO
      IF (n==Ptr % TYPE % NumberOfNodes) EXIT
    END DO
!------------------------------------------------------------------------------
  END FUNCTION Find_Face
!------------------------------------------------------------------------------

  
  !------------------------------------------------------------------------------------------------
  !> Finds nodes for which CandNodes are True such that their mutual distance is somehow
  !> maximized. We first find lower left corner, then the node that is furtherst apart from it,
  !> and continue as long as there are nodes to find. Typically we would be content with two nodes
  !> on a line, three nodes on a plane, and four nodes on a volume.
  !-------------------------------------------------------------------------------------------------
  SUBROUTINE FindExtremumNodes(Mesh,CandNodes,NoExt,Inds) 
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, ALLOCATABLE :: CandNodes(:)
    INTEGER :: NoExt
    INTEGER, POINTER :: Inds(:)

    REAL(KIND=dp) :: Coord(3),dCoord(3),dist,MinDist,MaxDist
    REAL(KIND=dp), ALLOCATABLE :: SetCoord(:,:)
    INTEGER :: i,j,k
    
    ALLOCATE( SetCoord(NoExt,3) )
    SetCoord = 0.0_dp
    Inds = 0
    
    ! First find the lower left corner
    MinDist = HUGE(MinDist) 
    DO i=1, Mesh % NumberOfNodes
      IF(.NOT. CandNodes(i) ) CYCLE
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
      Dist = SUM( Coord )
      IF( Dist < MinDist ) THEN
        Inds(1) = i
        MinDist = Dist
        SetCoord(1,:) = Coord
      END IF
    END DO
    
    ! Find more points such that their minimum distance to the previous point(s)
    ! is maximized.
    DO j=2,NoExt
      ! The maximum minimum distance of any node from the previously defined nodes
      MaxDist = 0.0_dp
      DO i=1, Mesh % NumberOfNodes
        IF(.NOT. CandNodes(i) ) CYCLE
        Coord(1) = Mesh % Nodes % x(i)
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        
        ! Minimum distance from the previously defined nodes
        MinDist = HUGE(MinDist)
        DO k=1,j-1
          dCoord = SetCoord(k,:) - Coord
          Dist = SUM( dCoord**2 )          
          MinDist = MIN( Dist, MinDist )
        END DO
        
        ! If the minimum distance is greater than in any other node, choose this
        IF( MaxDist < MinDist ) THEN
          MaxDist = MinDist 
          Inds(j) = i
          SetCoord(j,:) = Coord
        END IF
      END DO
    END DO

    PRINT *,'Extremum Inds:',Inds
    DO i=1,NoExt
      PRINT *,'Node:',Inds(i),SetCoord(i,:)
    END DO
    
  END SUBROUTINE FindExtremumNodes

  
  ! Save projector, mainly a utility for debugging purposes
  !--------------------------------------------------------
  SUBROUTINE SaveProjector(Projector,SaveRowSum,Prefix,InvPerm,Parallel)
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: SaveRowSum 
    CHARACTER(LEN=*) :: Prefix
    INTEGER, POINTER, OPTIONAL :: InvPerm(:)
    LOGICAL, OPTIONAL :: Parallel

    CHARACTER(LEN=MAX_NAME_LEN) :: Filename
    INTEGER :: i,j,ii,jj
    REAL(KIND=dp) :: rowsum, dia, val
    INTEGER, POINTER :: IntInvPerm(:)
    LOGICAL :: GlobalInds
    INTEGER, POINTER :: GlobalDofs(:)
    
    IF(.NOT.ASSOCIATED(Projector)) RETURN
    
    IF( PRESENT( InvPerm ) ) THEN
      IntInvPerm => InvPerm 
    ELSE
      IntInvPerm => Projector % InvPerm
    END IF

    GlobalInds = .FALSE.
    IF(ParEnv % PEs == 1 ) THEN
      FileName = TRIM(Prefix)//'.dat'
    ELSE
      FileName = TRIM(Prefix)//'_part'//&
          TRIM(I2S(ParEnv % MyPe))//'.dat'
      IF( PRESENT( Parallel ) ) GlobalInds = Parallel
    END IF

    IF( GlobalInds ) THEN
      NULLIFY( GlobalDofs ) 
      IF( ASSOCIATED( CurrentModel % Solver % Matrix ) ) THEN
        GlobalDofs => CurrentModel % Solver % Matrix % ParallelInfo % GlobalDofs
      END IF
      IF(.NOT. ASSOCIATED( GlobalDofs ) ) THEN
        CALL Info('SaveProjector','Cannot find GlobalDofs for Solver matrix')
        GlobalDofs => CurrentModel % Mesh % ParallelInfo % GlobalDofs
      END IF
    END IF
          
    OPEN(1,FILE=FileName,STATUS='Unknown')    
    DO i=1,projector % numberofrows
      IF( ASSOCIATED( IntInvPerm ) ) THEN
        ii = intinvperm(i)        
        IF( ii == 0) THEN
          PRINT *,'Projector InvPerm is zero:',ParEnv % MyPe, i, ii
          CYCLE
        END IF
      ELSE
        ii = i
      END IF
      IF( GlobalInds ) THEN
        IF( ii > SIZE( GlobalDofs ) ) THEN
          PRINT *,'ParEnv % MyPe, Projecor invperm is larger than globaldofs',&
              ii, SIZE( GlobalDofs ), i, Projector % NumberOfRows
          CYCLE
        END IF
        ii = GlobalDofs(ii)
      END IF
      IF( ii == 0) THEN
        PRINT *,'Projector global InvPerm is zero:',ParEnv % MyPe, i, ii
        CYCLE
      END IF
      DO j=projector % rows(i), projector % rows(i+1)-1
        jj = projector % cols(j)
        IF( jj == 0) THEN
          PRINT *,'Projector col is zero:',ParEnv % MyPe, i, ii, j, jj
          CYCLE
        END IF       
        val = projector % values(j)
        IF( GlobalInds ) THEN
          IF( jj > SIZE( GlobalDofs ) ) THEN
            PRINT *,'Projecor invperm is larger than globaldofs',&
                jj, SIZE( GlobalDofs )
            CYCLE
          END IF
          jj = GlobalDofs(jj)
          IF( jj == 0) THEN
            PRINT *,'Projector global col is zero:',ParEnv % MyPe, i, ii, j, jj
            CYCLE
          END IF
          WRITE(1,*) ii,jj,ParEnv % MyPe, val
        ELSE
          WRITE(1,*) ii,jj,val
        END IF
      END DO
    END DO
    CLOSE(1)     

    IF( SaveRowSum ) THEN
      IF(ParEnv % PEs == 1 ) THEN
        FileName = TRIM(Prefix)//'_rsum.dat'
      ELSE
        FileName = TRIM(Prefix)//'_rsum_part'//&
            TRIM(I2S(ParEnv % MyPe))//'.dat'
      END IF
      
      OPEN(1,FILE=FileName,STATUS='Unknown')
      DO i=1,projector % numberofrows
        IF( ASSOCIATED( IntInvPerm ) ) THEN
          ii = intinvperm(i)
          IF( ii == 0 ) CYCLE
        ELSE
          ii = i
        END IF
        rowsum = 0.0_dp
        dia = 0.0_dp

        DO j=projector % rows(i), projector % rows(i+1)-1          
          jj = projector % cols(j)
          val = projector % values(j)
          IF( ii == jj ) THEN
            dia = val
          END IF
          rowsum = rowsum + val
        END DO

        IF( GlobalInds ) THEN
          ii = GlobalDofs(ii)
          WRITE(1,*) ii, i, &
              projector % rows(i+1)-projector % rows(i), ParEnv % MyPe, dia, rowsum
        ELSE
          WRITE(1,*) ii, i, &
              projector % rows(i+1)-projector % rows(i),dia, rowsum
        END IF

      END DO
      CLOSE(1)     
    END IF

  END SUBROUTINE SaveProjector



  ! Set projector abs(rowsum) to unity
  !--------------------------------------------------------
  SUBROUTINE SetProjectorRowsum( Projector )
    TYPE(Matrix_t), POINTER :: Projector

    INTEGER :: i,j
    REAL(KIND=dp) :: rowsum

    DO i=1,projector % numberofrows
      rowsum = 0.0_dp
      DO j=projector % rows(i), projector % rows(i+1)-1
        rowsum = rowsum + ABS( projector % values(j) )
      END DO
      DO j=projector % rows(i), projector % rows(i+1)-1
        projector % values(j) = projector % values(j) / rowsum
      END DO
    END DO

  END SUBROUTINE SetProjectorRowsum




!------------------------------------------------------------------------------
!> Writes the mesh to disk. Note that this does not include the information
!> of shared nodes needed in parallel computation. This may be used for 
!> debugging purposes and for adaptive solution, for example. 
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk( NewMesh, Path )
!------------------------------------------------------------------------------
    CHARACTER(LEN=*) :: Path
    TYPE(Mesh_t), POINTER :: NewMesh
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,MaxNodes,ElmCode,Parent1,Parent2
!------------------------------------------------------------------------------

    OPEN( 1,FILE=TRIM(Path) // '/mesh.header',STATUS='UNKNOWN' )
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, NewMesh % NumberOfBoundaryElements
    
    WRITE( 1,'(i0)' ) 2
    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(k) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(k) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(k) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBoundaryElements

    MaxNodes = 0
    ElmCode  = 0
    DO i=1,NewMesh % NumberOfBulkElements
       IF ( NewMesh % Elements(i) % TYPE % NumberOfNodes > MaxNodes ) THEN
          ElmCode  = NewMesh % Elements(i) % TYPE % ElementCode
          MaxNodes = NewMesh % Elements(i) % TYPE % NumberOfNodes
       END IF
    END DO
    WRITE( 1,'(i0,x,i0)' ) ElmCode,NewMesh % NumberOfBulkElements
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.nodes', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       WRITE(1,'(i0,a,3e23.15)',ADVANCE='NO') i,' -1 ', &
            NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.elements', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       WRITE(1,'(3(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % Elements(i) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    OPEN( 1,FILE=TRIM(Path) // '/mesh.boundary', STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex
       WRITE(1,'(5(i0,x))',ADVANCE='NO') i, &
            NewMesh % Elements(k) % BoundaryInfo % Constraint, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % Elements(k) % NodeIndexes(j)
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)
!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDisk2(Model, NewMesh, Path, Partition )
!------------------------------------------------------------------------------
    USE Types
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: NewMesh
    CHARACTER(LEN=*) :: Path
    INTEGER, OPTIONAL :: Partition
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeList(100),ElmCodeCounts(100),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    INTEGER, POINTER :: BList(:)
    INTEGER, ALLOCATABLE :: ElementCodes(:)
    LOGICAL :: Parallel, WarnNoTarget, Found
    CHARACTER(LEN=MAX_NAME_LEN) :: headerFN, elementFN, nodeFN,&
         boundFN, sharedFN
!------------------------------------------------------------------------------

    IF(PRESENT(Partition)) THEN
       Parallel = .TRUE.
       WRITE(headerFN, '(A,I0,A)') '/part.',Partition+1,'.header'
       WRITE(elementFN, '(A,I0,A)') '/part.',Partition+1,'.elements'
       WRITE(nodeFN, '(A,I0,A)') '/part.',Partition+1,'.nodes'
       WRITE(boundFN, '(A,I0,A)') '/part.',Partition+1,'.boundary'
       WRITE(sharedFN, '(A,I0,A)') '/part.',Partition+1,'.shared'
    ELSE
       Parallel = .FALSE.
       headerFN = '/mesh.header'
       elementFN = '/mesh.elements'
       nodeFN = '/mesh.nodes'
       boundFN = '/mesh.boundary'
    END IF

    !Info for header file

    ElmCodeList = 0 !init array
    NumElmCodes = 0
    NumElements = NewMesh % NumberOfBoundaryElements + &
         NewMesh % NumberOfBulkElements
    ALLOCATE(ElementCodes(NumElements))

    !cycle to bring element code list into array-inquirable form
    DO i=1,NumElements
       ElementCodes(i) = NewMesh % Elements(i) % TYPE % ElementCode
    END DO

    DO i=NumElements,1,-1 !this should give element codes increasing value, which appears to be
                          !'standard' though I doubt it matters
       IF(ANY(ElmCodeList == ElementCodes(i))) CYCLE
       NumElmCodes = NumElmCodes + 1
       ElmCodeList(NumElmCodes) = ElementCodes(i)
    END DO

    DO j=1,NumElmCodes
       ElmCodeCounts(j) = COUNT(ElementCodes == ElmCodeList(j))
    END DO

    !Write header file
    OPEN( 1,FILE=TRIM(Path) // headerFN,STATUS='UNKNOWN' )
    WRITE( 1,'(i0,x,i0,x,i0)' ) NewMesh % NumberOfNodes, &
         NewMesh % NumberOfBulkElements, &
         NewMesh % NumberOfBoundaryElements

    WRITE( 1,'(i0)' ) NumElmCodes
    DO j=1,NumElmCodes
       WRITE( 1,'(i0,x,i0,x)' ) ElmCodeList(j),ElmCodeCounts(j)
    END DO
    IF(Parallel) THEN !need number of shared nodes
       NoShared = 0
       DO i=1,NewMesh % NumberOfNodes
          IF(SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
               Neighbours) > 1) THEN
             NoShared = NoShared + 1
          END IF
       END DO
       WRITE( 1,'(i0,x,i0)') NoShared, 0
    END IF
    CLOSE(1)

    !Write nodes file
    OPEN( 1,FILE=TRIM(Path) // nodeFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       IF (Parallel) THEN
          WRITE(1,'(i0,x)', ADVANCE='NO') &
               NewMesh % ParallelInfo % GlobalDOFs(i)
       ELSE
          WRITE(1,'(i0,x)', ADVANCE='NO') i
       END IF
       WRITE(1,'(a,x,ES17.10,x,ES17.10,x,ES17.10)',ADVANCE='NO') &
            ' -1 ', NewMesh % Nodes % x(i), &
            NewMesh % Nodes % y(i), NewMesh % Nodes % z(i)
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)

    !Write elements file
    OPEN( 1,FILE=TRIM(Path) // elementFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBulkElements
       IF(Parallel) THEN
          ElemID = NewMesh % Elements(i) % GElementIndex
       ELSE
          ElemID = i
       END IF
       WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') ElemID, &
            NewMesh % Elements(i) % BodyId, &
            NewMesh % Elements(i) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(i) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(i) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(i) % NodeIndexes(j)
          END IF
          WRITE(1,'(i0,x)', ADVANCE='NO') m
       END DO
       WRITE(1,*) ''
    END DO
    CLOSE(1)

    !Write boundary file
    WarnNoTarget = .FALSE.
    OPEN( 1,FILE=TRIM(Path) // boundFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfBoundaryElements
       k = i + NewMesh % NumberOfBulkElements
       parent1 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Left ) ) &
          parent1 = NewMesh % Elements(k) % BoundaryInfo % Left % ElementIndex
       parent2 = 0
       IF ( ASSOCIATED( NewMesh % Elements(k) % BoundaryInfo % Right ) ) &
          parent2 = NewMesh % Elements(k) % BoundaryInfo % Right % ElementIndex

       IF(Parallel) THEN
          IF(parent1 /= 0) parent1 = NewMesh % Elements(parent1) % GElementIndex
          IF(parent2 /= 0) parent2 = NewMesh % Elements(parent2) % GElementIndex
       END IF

       Constraint = NewMesh % Elements(k) % BoundaryInfo % Constraint
       BList => ListGetIntegerArray( Model % BCs(Constraint) % Values, &
            'Target Boundaries', Found )
       IF(Found) THEN
          IF(SIZE(BList) > 1) THEN
             CALL WARN("WriteMeshToDisk2",&
                  "A BC has more than one Target Boundary, SaveMesh output will not match input!")
          END IF
          meshBC = BList(1)
       ELSE
          WarnNoTarget = .TRUE.
          meshBC = Constraint
       END IF

       !This meshBC stuff will *only* work if each BC has only 1 target boundary
       WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            meshBC, Parent1,Parent2,&
            NewMesh % Elements(k) % TYPE % ElementCode
       DO j=1,NewMesh % Elements(k) % TYPE % NumberOfNodes
          IF(Parallel) THEN
             m = NewMesh % ParallelInfo % GlobalDOFs(&
                  NewMesh % Elements(k) % NodeIndexes(j))
          ELSE
             m = NewMesh % Elements(k) % NodeIndexes(j)
          END IF
          WRITE(1,'(x,i0)', ADVANCE='NO') m
       END DO
       WRITE(1,*) !blank write statement to create new line without extra space.
    END DO
    CLOSE(1)

    IF(WarnNoTarget) THEN
       CALL WARN("WriteMeshToDisk2","Couldn't find a Target Boundary, assuming mapping to self")
    END IF

    IF(.NOT. Parallel) RETURN

    !Write .shared file
    !Need to create part.n.shared from Mesh % ParallelInfo %
    !NeighbourList % Neighbours.
    OPEN( 1,FILE=TRIM(Path) // sharedFN, STATUS='UNKNOWN' )
    DO i=1,NewMesh % NumberOfNodes
       nneigh = SIZE(NewMesh % ParallelInfo % NeighbourList(i) % &
            Neighbours)
       IF(nneigh < 2) CYCLE
       WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') &
            NewMesh % ParallelInfo % GlobalDOFs(i),nneigh
       DO j=1,nneigh
          WRITE(1,'(I0, x)',ADVANCE='NO') NewMesh % ParallelInfo %&
               NeighbourList(i) % Neighbours(j) + 1
       END DO
       WRITE( 1,* ) ''
    END DO
    CLOSE(1)


!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDisk2
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Writes the mesh to disk, including detection of elementcodes and shared node
!> info necessary for parallel meshes.
!------------------------------------------------------------------------------
  SUBROUTINE WriteMeshToDiskPartitioned(Model, Mesh, Path, &
      ElementPart, NeighbourList )
!------------------------------------------------------------------------------
    USE Types
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: Path
    INTEGER, POINTER :: ElementPart(:)
    TYPE(NeighbourList_t),POINTER  :: NeighbourList(:)
!------------------------------------------------------------------------------
    TYPE(Element_t), POINTER :: Element
    INTEGER :: NoBoundaryElements, NoBulkElements, NoNodes, NoPartitions, Partition
    INTEGER :: i,j,k,m,MaxNodes,ElmCode,NumElmCodes,ElmCodeCounts(827),&
         Parent1,Parent2, ElemID, nneigh, Constraint, meshBC, NumElements, NoShared
    LOGICAL :: Found, Hit
    CHARACTER(LEN=MAX_NAME_LEN) :: DirectoryName, PrefixName
!------------------------------------------------------------------------------

    NoPartitions = MAXVAL( ElementPart ) 
    NumElmCodes = 0
    NumElements = Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements
        
    WRITE(DirectoryName, '(A,A,I0)') TRIM(PATH),'/partitioning.',NoPartitions
    CALL MakeDirectory( TRIM(DirectoryName) // CHAR(0) )
    CALL Info('WriteMeshToDiskPartitioned','Writing parallel mesh to disk: '//TRIM(DirectoryName))
   

    DO Partition = 1, NoPartitions 
      
      CALL Info('WriteMeshToDiskPartitioned','Writing piece to file: '//TRIM(I2S(Partition)),Level=12)
      
      WRITE( PrefixName,'(A,A,I0)') TRIM(DirectoryName),'/part.',Partition  

      CALL Info('WriteMeshToDiskPartitioned','Write nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.nodes', STATUS='UNKNOWN' )
      NoNodes = 0
      DO i=1,Mesh % NumberOfNodes
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          WRITE(1,'(I0,x,I0,x,3ES17.10)') i,-1, &
              Mesh % Nodes % x(i), Mesh % Nodes % y(i), Mesh % Nodes % z(i)
          NoNodes = NoNodes + 1
        END IF
      END DO
      CLOSE(1)
      

      CALL Info('WriteMeshToDiskPartitioned','Write shared nodes file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.shared', STATUS='UNKNOWN' )
      NoShared = 0
      DO i=1,Mesh % NumberOfNodes
        nneigh = SIZE( NeighbourList(i) % Neighbours )
        IF( nneigh <= 1 ) CYCLE
        
        IF( ANY( NeighbourList(i) % Neighbours == Partition ) ) THEN
          NoShared = NoShared + 1
          WRITE(1,'(i0, x, i0, x)',ADVANCE='NO') i,nneigh
          DO j=1,nneigh
            WRITE(1,'(I0, x)',ADVANCE='NO') NeighbourList(i) % Neighbours(j) 
          END DO
          WRITE( 1,* ) ''
        END IF
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write elements file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.elements', STATUS='UNKNOWN' )
      NoBulkElements = 0
      ElmCodeCounts = 0      
      DO i=1,Mesh % NumberOfBulkElements
        IF( ElementPart(i) /= Partition ) CYCLE

        Element => Mesh % Elements(i)
        WRITE(1,'(i0,x,i0,x,i0,x)',ADVANCE='NO') i, &
            Element % BodyId, Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(i0,x)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) ''
        
        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBulkElements = NoBulkElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write boundary file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.boundary', STATUS='UNKNOWN' )
      NoBoundaryElements = 0
      DO i=Mesh % NumberOfBulkElements +1 ,&
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        Element => Mesh % Elements(i)
       
        parent1 = 0
        parent2 = 0
        Constraint = 0
        
        IF( ASSOCIATED( Element % BoundaryInfo ) ) THEN
          IF ( ASSOCIATED( Element % BoundaryInfo % Left ) ) &
              parent1 = Element % BoundaryInfo % Left % ElementIndex
          IF ( ASSOCIATED( Element % BoundaryInfo % Right ) ) &
              parent2 = Element % BoundaryInfo % Right % ElementIndex        
          Constraint = Element % BoundaryInfo % Constraint
        END IF

        Hit = .FALSE.
        IF( parent1 > 0 ) THEN
          IF( ElementPart( parent1 ) == Partition ) Hit = .TRUE.
        END IF
        IF( parent2 > 0 ) THEN
          IF( ElementPart( parent2 ) == Partition ) Hit = .TRUE.
        END IF

        IF( .NOT. Hit ) CYCLE

        WRITE(1,'(i0,x,i0,x,i0,x,i0,x,i0)',ADVANCE='NO') i, & 
            Constraint, Parent1, Parent2,&
            Element % TYPE % ElementCode
        DO j=1,Element % TYPE % NumberOfNodes
          WRITE(1,'(x,i0)', ADVANCE='NO') Element % NodeIndexes(j)
        END DO
        WRITE(1,*) 

        ElmCode = Element % TYPE % ElementCode
        ElmCodeCounts( ElmCode ) = ElmCodeCounts( ElmCode ) + 1
        NoBoundaryElements = NoBoundaryElements + 1
      END DO
      CLOSE(1)


      CALL Info('WriteMeshToDiskPartitioned','Write header file',Level=12)
      OPEN( 1,FILE=TRIM(PrefixName) // '.header',STATUS='UNKNOWN' )
      NumElmCodes = COUNT( ElmCodeCounts > 0 ) 
      WRITE( 1,'(i0,x,i0,x,i0)' ) NoNodes, &
          NoBulkElements, NoBoundaryElements      
      WRITE( 1,'(i0)' ) NumElmCodes
      DO i=SIZE(ElmCodeCounts),1,-1
        IF( ElmCodeCounts(i) == 0 ) CYCLE
        WRITE( 1,'(i0,x,i0,x)' ) i,ElmCodeCounts(i)
      END DO
      WRITE( 1,'(i0,x,i0)') NoShared, 0
      CLOSE(1)
      
      CALL Info('WriteMeshToDiskPartitioned','Done writing partition',Level=12)
    END DO

    CALL Info('WriteMeshToDiskPartitioned','Done writing parallel mesh',Level=8)

!------------------------------------------------------------------------------
  END SUBROUTINE WriteMeshToDiskPartitioned
!------------------------------------------------------------------------------




!------------------------------------------------------------------------------
!> Generate element edge (faces in 3D) tables for given mesh.
!> Currently only for triangles and tetras. If mesh already
!> has edges do nothing.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges( Mesh, FindEdges)
!------------------------------------------------------------------------------
     TYPE(Mesh_t) :: Mesh
     LOGICAL, OPTIONAL :: FindEdges

     LOGICAL :: FindEdges3D
     INTEGER :: MeshDim, SpaceDim, MaxElemDim 

     IF(PRESENT(FindEdges)) THEN
       FindEdges3D = FindEdges
     ELSE
       FindEdges3D = .TRUE.
     END IF

!------------------------------------------------------------------------------

     SpaceDim = CoordinateSystemDimension()
     MeshDim = Mesh % MeshDim

     IF( MeshDim == 0 ) THEN
       CALL Fatal('FindMeshEdges','Mesh dimension is zero!')
     END IF
     IF( SpaceDim > MeshDim ) THEN
       CALL Warn('FindMeshEdges','Mesh dimension and space dimension differ: '&
           // TRIM(I2S(MeshDim))//' vs. '//TRIM(I2S(SpaceDim)))
     END IF

     MaxElemDim = EnsureElemDim( MeshDim ) 
     IF( MaxElemDim < MeshDim ) THEN
       CALL Warn('FindMeshEdges','Element dimension smaller than mesh dimension: '//&
           TRIM(I2S(MaxElemDim))//' vs '//TRIM(I2S(MeshDim)))
     END IF


     SELECT CASE( MaxElemDim )

     CASE(2)
       IF ( .NOT.ASSOCIATED( Mesh % Edges ) ) THEN
         CALL Info('FindMeshEdges','Determining edges in 2D mesh',Level=8)
         CALL FindMeshEdges2D( Mesh )
       END IF

     CASE(3)
       IF ( .NOT.ASSOCIATED( Mesh % Faces) ) THEN
         CALL Info('FindMeshEdges','Determining faces in 3D mesh',Level=8)
         CALL FindMeshFaces3D( Mesh )
       END IF
       IF(FindEdges3D) THEN
         IF ( .NOT.ASSOCIATED( Mesh % Edges) ) THEN
           CALL Info('FindMeshEdges','Determining edges in 3D mesh',Level=8)
           CALL FindMeshEdges3D( Mesh )
         END IF
       END IF
     END SELECT

     CALL AssignConstraints()

CONTAINS

  ! Check that the element dimension really follows the mesh dimension
  ! The default is the MeshDim so we return immediately after that is 
  ! confirmed. 
  !--------------------------------------------------------------------
    FUNCTION EnsureElemDim(MeshDim) RESULT (MaxElemDim)

      INTEGER :: MeshDim, MaxElemDim 
      INTEGER :: i,ElemDim, ElemCode

      MaxElemDim = 0

      DO i=1,Mesh % NumberOfBulkElements
        ElemCode = Mesh % Elements(i) % Type % ElementCode
        IF( ElemCode > 500 ) THEN
          ElemDim = 3 
        ELSE IF( ElemCode > 300 ) THEN
          ElemDim = 2
        ELSE IF( ElemCode > 200 ) THEN
          ElemDim = 1
        END IF
        MaxElemDim = MAX( MaxElemDim, ElemDim ) 
        IF( MaxElemDim == MeshDim ) EXIT
      END DO
          
    END FUNCTION EnsureElemDim


    SUBROUTINE AssignConstraints()

      INTEGER, POINTER :: FaceInd(:)
      INTEGER :: i,j,k,l,n,nd,nfound
      TYPE(Element_t), POINTER :: Element, Boundary, Face, Faces(:)

      DO i=1,Mesh % NumberOfBoundaryElements
        Boundary => Mesh % Elements(Mesh % NumberOfBulkElements+i)

        Element  => Boundary % BoundaryInfo % Left
        IF (.NOT.ASSOCIATED(Element) ) &
          Element  => Boundary % BoundaryInfo % Right
        IF (.NOT.ASSOCIATED(Element) ) CYCLE

        SELECT CASE(Boundary % TYPE % DIMENSION)
        CASE(1)
          nd = Element % TYPE % NumberOfEdges
          Faces   => Mesh % Edges
          FaceInd => Element % EdgeIndexes
        CASE(2)
          nd = Element % TYPE % NumberOfFaces
          Faces   => Mesh % Faces
          FaceInd => Element % FaceIndexes
        CASE DEFAULT
          Faces => NULL()
          FaceInd => NULL()
        END SELECT

        IF ( .NOT. ASSOCIATED(Faces) .OR. .NOT. ASSOCIATED(FaceInd) ) CYCLE

        DO j=1,nd
          Face => Faces(FaceInd(j))
          IF ( .NOT.ASSOCIATED(Face % TYPE,Boundary % TYPE) ) CYCLE

          n = Boundary % TYPE % NumberOfNodes
          nfound = 0
          DO k=1,n
            DO l=1,n
              IF ( Boundary % NodeIndexes(k)==Face % NodeIndexes(l) ) &
                nfound = nfound+1
            END DO
          END DO
          IF ( nfound==n ) THEN
            Face % BoundaryInfo % Constraint = Boundary % BoundaryInfo % Constraint; EXIT
          END IF
        END DO
      END DO
    END SUBROUTINE AssignConstraints
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Find 2D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges2D( Mesh, BulkMask )
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: BulkMask(:)
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node,Edge
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
     
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    TYPE(Element_t), POINTER :: Element, Edges(:)

    LOGICAL :: Found,Masked
    INTEGER :: i,j,k,n,NofEdges,Edge,Swap,Node1,Node2,istat,Degree,allocstat
!------------------------------------------------------------------------------
!
!   Initialize:
!   -----------
    CALL Info('FindMeshEdges2D','Allocating edge table of size: '&
        //TRIM(I2S(4*Mesh % NumberOfBulkElements)),Level=12)

    Masked = PRESENT(BulkMask)
    
    CALL AllocateVector( Mesh % Edges, 4*Mesh % NumberOfBulkElements )
    Edges => Mesh % Edges

    DO i=1,Mesh % NumberOfBulkElements
      IF(Masked) THEN
        IF(.NOT. BulkMask(i)) CYCLE
      END IF
       Element => Mesh % Elements(i)

       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector( Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    CALL Info('FindMeshEdges2D','Creating hash table of size '&
        //TRIM(I2S(Mesh % NumberOfNodes))//' for noto-to-node connectivity',Level=12)
    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

!   Loop over elements:
!   -------------------
    NofEdges = 0
    DO i=1,Mesh % NumberOfBulkElements

       IF(Masked) THEN
         IF(.NOT. BulkMask(i)) CYCLE
       END IF

       Element => Mesh % Elements(i)

       SELECT CASE( Element % TYPE % ElementCode / 100 )
         CASE(3)
            n = 3
         CASE(4)
            n = 4
       END SELECT

!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n
!         We use MIN(Node1,Node2) as the hash table key:
!         ----------------------------------------------
          Node1 = Element % NodeIndexes(k)
          IF ( k<n ) THEN
             Node2 = Element % NodeIndexes(k+1)
          ELSE
             Node2 = Element % NodeIndexes(1)
          END IF

          IF ( Node2 < Node1 ) THEN
             Swap  = Node1
             Node1 = Node2
             Node2 = Swap
          END IF

!         Look the edge from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.         
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node == Node2 ) THEN
                Found = .TRUE.
                Edge = HashPtr % Edge
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO

!         Existing edge, update structures:
!         ----------------------------------
          IF ( Found ) THEN
             Element % EdgeIndexes(k) = Edge
             Edges(Edge) % BoundaryInfo % Right => Element
          ELSE

!            Edge not yet there, create:
!            ---------------------------
             NofEdges = NofEdges + 1
             Edge = NofEdges

             Degree = Element % TYPE % BasisFunctionDegree

             Edges(Edge) % ElementIndex = Edge
             CALL AllocateVector( Edges(Edge) % NodeIndexes, Degree+1)
             ALLOCATE( Edges(Edge) % BoundaryInfo, STAT=allocstat )
             IF( allocstat /= 0 ) THEN
               CALL Fatal('FindMeshEdges2D','Allocation error for BoyndaryInfo alloction')
             END IF

             Edges(Edge) % TYPE => GetElementType( 201+Degree, .FALSE. )

             Edges(Edge) % NodeIndexes(1) = Element % NodeIndexes(k)
             IF ( k < n ) THEN
                Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(k+1)
             ELSE
                Edges(Edge) % NodeIndexes(2) = Element % NodeIndexes(1)
             END IF

             DO j=2,Degree
                Edges(Edge) % NodeIndexes(j+1) = Element % NodeIndexes(k+n+j-2)
             END DO
             
             ! Create P element definitions if needed
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
               CALL AllocatePDefinitions(Edges(Edge))
               Edges(Edge) % PDefs % P = 0
             ELSE
               NULLIFY( Edges(Edge) % PDefs )
             END IF

             Edges(Edge) % NDofs = 0
             IF (Element % NDOFs /= 0 ) &
                Edges(Edge) % NDOFs  = Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             NULLIFY( Edges(Edge) % EdgeIndexes )
             NULLIFY( Edges(Edge) % FaceIndexes )

             Element % EdgeIndexes(k) = Edge

             Edges(Edge) % BoundaryInfo % Left => Element
             NULLIFY( Edges(Edge) % BoundaryInfo % Right )
              
!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr, STAT=allocstat )
             IF( allocstat /= 0 ) THEN
               CALL Fatal('FindMeshEdges2D','Allocation error for HashPtr alloction')
             END IF

             HashPtr % Edge = Edge
             HashPtr % Node = Node2
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO

    Mesh % NumberOfEdges = NofEdges
    CALL Info('FindMeshEdges2D','Number of edges found: '//TRIM(I2S(NofEdges)),Level=10)

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )

    CALL Info('FindMeshEdges2D','All done',Level=12)

!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges2D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh faces.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshFaces3D( Mesh, BulkMask)
    USE PElementMaps, ONLY : GetElementFaceMap
    USE PElementBase, ONLY : isPTetra

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: BulkMask(:)
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node1,Node2,Face
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
    
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    LOGICAL :: Found,Masked
    INTEGER :: n1,n2,n3,n4
    INTEGER :: i,j,k,n,NofFaces,Face,Swap,Node1,Node2,Node3,istat,Degree
     
    TYPE(Element_t), POINTER :: Element, Faces(:)

    INTEGER, POINTER :: FaceMap(:,:)
    INTEGER, TARGET  :: TetraFaceMap(4,6), BrickFaceMap(6,9), &
         WedgeFaceMap(5,8), PyramidFaceMap(5,8)
    
    INTEGER :: nf(4)
!------------------------------------------------------------------------------
    
    CALL Info('FindMeshFaces3D','Finding mesh faces in 3D mesh',Level=12)

    Masked = PRESENT(BulkMask)

    TetraFaceMap(1,:) = [ 1, 2, 3, 5, 6, 7 ]
    TetraFaceMap(2,:) = [ 1, 2, 4, 5, 9, 8 ]
    TetraFaceMap(3,:) = [ 2, 3, 4, 6, 10, 9 ]
    TetraFaceMap(4,:) = [ 3, 1, 4, 7, 8,10 ]

    WedgeFaceMap(1,:) = [ 1, 2, 3, 7, 8, 9, -1, -1 ]
    WedgeFaceMap(2,:) = [ 4, 5, 6, 10, 11, 12, -1, -1 ]
    WedgeFaceMap(3,:) = [ 1, 2, 5, 4, 7, 14, 10, 13 ]
    WedgeFaceMap(4,:) = [ 3, 2, 5, 6, 8, 14, 11, 15 ]
    WedgeFaceMap(5,:) = [ 3, 1, 4, 6, 9, 13, 12, 15 ]

    PyramidFaceMap(1,:) = [ 1, 2, 3, 4,  6,  7,  8,  9 ]
    PyramidFaceMap(2,:) = [ 1, 2, 5, 6, 11, 10, -1, -1 ]
    PyramidFaceMap(3,:) = [ 2, 3, 5, 7, 12, 11, -1, -1 ]
    PyramidFaceMap(4,:) = [ 3, 4, 5, 8, 13, 12, -1, -1 ]
    PyramidFaceMap(5,:) = [ 4, 1, 5, 9, 10, 13, -1, -1 ]

    BrickFaceMap(1,:) = [ 1, 2, 3, 4,  9, 10, 11, 12, 25 ]
    BrickFaceMap(2,:) = [ 5, 6, 7, 8, 17, 18, 19, 20, 26 ]
    BrickFaceMap(3,:) = [ 1, 2, 6, 5,  9, 14, 17, 13, 21 ]
    BrickFaceMap(4,:) = [ 2, 3, 7, 6, 10, 15, 18, 14, 22 ]
    BrickFaceMap(5,:) = [ 3, 4, 8, 7, 11, 16, 19, 15, 23 ]
    BrickFaceMap(6,:) = [ 4, 1, 5, 8, 12, 13, 20, 16, 24 ]

!
!   Initialize:
!   -----------
    IF(Masked) THEN
      CALL AllocateVector( Mesh % Faces, 6*COUNT(BulkMask), 'FindMeshFaces3D' )
    ELSE
      CALL AllocateVector( Mesh % Faces, 6*Mesh % NumberOfBulkElements, 'FindMeshFaces3D' )
    END IF
    Faces => Mesh % Faces

    DO i=1,Mesh % NumberOfBulkElements
       IF(Masked) THEN
         IF(.NOT. BulkMask(i)) CYCLE
       END IF
       Element => Mesh % Elements(i)
       IF ( .NOT. ASSOCIATED( Element % FaceIndexes ) ) &
          CALL AllocateVector(Element % FaceIndexes, Element % TYPE % NumberOfFaces )
       Element % FaceIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

!   Loop over elements:
!   -------------------
    NofFaces = 0
    DO i=1,Mesh % NumberOfBulkElements
       IF(Masked) THEN
         IF(.NOT. BulkMask(i)) CYCLE
       END IF

       Element => Mesh % Elements(i)

       ! For P elements mappings are different
       IF ( ASSOCIATED(Element % PDefs) ) THEN
          CALL GetElementFaceMap(Element, FaceMap)
          n = Element % TYPE % NumberOfFaces
       ELSE
          SELECT CASE( Element % TYPE % ElementCode / 100 )
          CASE(5)
             n = 4
             FaceMap => TetraFaceMap
          CASE(6)
             n = 5
             FaceMap => PyramidFaceMap
          CASE(7)
             n = 5 
             FaceMap => WedgeFaceMap
          CASE(8)
             n = 6
             FaceMap => BrickFaceMap
          CASE DEFAULT
             CYCLE
             ! WRITE(Message,*) 'Element type',Element % Type % ElementCode,'not implemented.' 
             ! CALL Fatal('FindMeshFaces',Message)
          END SELECT
       END IF
 
!      Loop over every face of every element:
!      --------------------------------------
       DO k=1,n
          
          
!         We use MIN(Node1,Node2,Node3) as the hash table key:
!         ---------------------------------------------------
          SELECT CASE( Element % TYPE % ElementCode / 100 )
             CASE(5)
!
!               Tetras:
!               =======
                nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                CALL sort( 3, nf )

             CASE(6)
!
!               Pyramids:
!               =========
                IF ( k == 1 ) THEN
                   nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                   CALL sort( 4, nf )
                ELSE
                   nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                   CALL sort( 3, nf )
                END IF

             CASE(7)
!
!               Wedges:
!               =======
                IF ( k <= 2 ) THEN
                   nf(1:3) = Element % NodeIndexes(FaceMap(k,1:3))
                   CALL sort( 3, nf )
                ELSE
                   nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                   CALL sort( 4, nf )
                END IF
                
             CASE(8)
!
!               Bricks:
!               =======
                nf(1:4) = Element % NodeIndexes(FaceMap(k,1:4))
                CALL sort( 4, nf )

             CASE DEFAULT
                WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
                CALL Fatal('FindMeshFaces',Message)
          END SELECT

          Node1 = nf(1)
          Node2 = nf(2)
          Node3 = nf(3)
          
!         Look the face from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node1 == Node2 .AND. HashPtr % Node2 == Node3) THEN
                Found = .TRUE.
                Face = HashPtr % Face
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO

!         Existing face, update structures:
!         ----------------------------------
          IF ( Found ) THEN
             Element % FaceIndexes(k) = Face
             Faces(Face) % BoundaryInfo % Right => Element
          ELSE

!            Face not yet there, create:
!            ---------------------------
             NofFaces = NofFaces + 1
             Face = NofFaces
             Faces(Face) % ElementIndex = Face

             Degree = Element % TYPE % BasisFunctionDegree


             SELECT CASE( Element % TYPE % ElementCode / 100 )
             CASE(5)
               !
               !               for tetras:
               !               -----------
               SELECT CASE( Degree ) 
               CASE(1)
                 n1 = 3
               CASE(2)
                 n1 = 6
               CASE(3)
                 n1 = 10
               END SELECT

               Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )

             CASE(6)

               !               Pyramids ( 605 and 613 supported )
               !               -------------------------------
               IF ( k == 1 ) THEN
                 n1 = Degree * 4
                 Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
               ELSE
                 n1 = Degree * 3
                 Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
               END IF

             CASE(7)

               !               for wedges, 706 and 715 supported:
               !               -------------------------------
               IF ( k <= 2 ) THEN
                 n1 = Degree * 3
                 Faces(Face) % TYPE => GetElementType( 300+n1, .FALSE. )
               ELSE
                 n1 = Degree * 4
                 Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE. )
               END IF


             CASE(8)
               !
               !               for bricks:
               !               -----------
               SELECT CASE( Element % TYPE % NumberOfNodes ) 
               CASE(8)
                 n1 = 4
               CASE(20)
                 n1 = 8
               CASE(27)
                 n1 = 9
               END SELECT

               Faces(Face) % TYPE => GetElementType( 400+n1, .FALSE.)

             CASE DEFAULT
               WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
               CALL Fatal('FindMeshFaces',Message)

             END SELECT

             ! Allocate p structures for p elements
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
                CALL AllocatePDefinitions(Faces(Face))
                Faces(Face) % PDefs % P = 0
             ELSE
               NULLIFY( Faces(Face) % PDefs )
             END IF
             
             Faces(Face) % NDOFs  = 0
             IF (Element % NDOFs /= 0 ) &
                Faces(Face) % NDOFs  = Faces(Face) % TYPE % NumberOfNodes
             Faces(Face) % BDOFs  = 0
             Faces(Face) % DGDOFs = 0
             Faces(Face) % EdgeIndexes => NULL()
             Faces(Face) % FaceIndexes => NULL()

             CALL AllocateVector( Faces(Face) % NodeIndexes,n1 )
             DO n2=1,n1
                Faces(Face) % NodeIndexes(n2) = &
                         Element % NodeIndexes(FaceMap(k,n2)) 
             END DO

             Element % FaceIndexes(k) = Face

             ALLOCATE( Faces(Face) % BoundaryInfo )
             Faces(Face) % BoundaryInfo % Left => Element
             NULLIFY( Faces(Face) % BoundaryInfo % Right )
              
!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Face = Face
             HashPtr % Node1 = Node2
             HashPtr % Node2 = Node3
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO

    Mesh % NumberOfFaces = NofFaces
    CALL Info('FindMeshFaces3D','Number of faces found: '//TRIM(I2S(NofFaces)),Level=10)

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )

    CALL Info('FindMeshFaces3D','All done',Level=12)
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshFaces3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Find 3D mesh edges.
!------------------------------------------------------------------------------
  SUBROUTINE FindMeshEdges3D( Mesh )
    USE PElementMaps, ONLY : GetElementEdgeMap, GetElementFaceEdgeMap
    USE PElementBase, ONLY : isPPyramid

    IMPLICIT NONE
!------------------------------------------------------------------------------
    TYPE(Mesh_t) :: Mesh
!------------------------------------------------------------------------------
    TYPE HashEntry_t
       INTEGER :: Node1,Edge
       TYPE(HashEntry_t), POINTER :: Next
    END TYPE HashEntry_t

    TYPE HashTable_t
       TYPE(HashEntry_t), POINTER :: Head
    END TYPE HashTable_t
    
    TYPE(HashTable_t), ALLOCATABLE :: HashTable(:)
    TYPE(HashEntry_t), POINTER :: HashPtr, HashPtr1

    LOGICAL :: Found
    INTEGER :: n1,n2
    INTEGER :: i,j,k,n,NofEdges,Edge,Node1,Node2,istat,Degree,ii,jj
     
    TYPE(Element_t), POINTER :: Element, Edges(:), Face

    INTEGER, POINTER :: EdgeMap(:,:), FaceEdgeMap(:,:)
    INTEGER, TARGET  :: TetraEdgeMap(6,3), BrickEdgeMap(12,3), TetraFaceMap(4,6), &
      WedgeEdgeMap(9,3), PyramidEdgeMap(8,3), TetraFaceEdgeMap(4,3), &
      BrickFaceEdgeMap(8,4), WedgeFaceEdgeMap(6,4), PyramidFaceEdgeMap(5,4)
!------------------------------------------------------------------------------

    CALL Info('FindMeshEdges3D','Finding mesh edges in 3D mesh',Level=12)

    TetraFaceMap(1,:) = [ 1, 2, 3, 5, 6, 7 ]
    TetraFaceMap(2,:) = [ 1, 2, 4, 5, 9, 8 ]
    TetraFaceMap(3,:) = [ 2, 3, 4, 6,10, 9 ]
    TetraFaceMap(4,:) = [ 3, 1, 4, 7, 8,10 ]

    TetraFaceEdgeMap(1,:) = [ 1,2,3 ]
    TetraFaceEdgeMap(2,:) = [ 1,5,4 ]
    TetraFaceEdgeMap(3,:) = [ 2,6,5 ]
    TetraFaceEdgeMap(4,:) = [ 3,4,6 ]

    TetraEdgeMap(1,:) = [ 1,2,5 ]
    TetraEdgeMap(2,:) = [ 2,3,6 ]
    TetraEdgeMap(3,:) = [ 3,1,7 ]
    TetraEdgeMap(4,:) = [ 1,4,8 ]
    TetraEdgeMap(5,:) = [ 2,4,9 ]
    TetraEdgeMap(6,:) = [ 3,4,10 ]

    PyramidEdgeMap(1,:) = [ 1,2,1 ]
    PyramidEdgeMap(2,:) = [ 2,3,1 ]
    PyramidEdgeMap(3,:) = [ 3,4,1 ]
    PyramidEdgeMap(4,:) = [ 4,1,1 ]
    PyramidEdgeMap(5,:) = [ 1,5,1 ]
    PyramidEdgeMap(6,:) = [ 2,5,1 ]
    PyramidEdgeMap(7,:) = [ 3,5,1 ]
    PyramidEdgeMap(8,:) = [ 4,5,1 ]

    PyramidFaceEdgeMap(1,:) = [ 1,2,3,4 ]
    PyramidFaceEdgeMap(2,:) = [ 1,6,5,0 ]
    PyramidFaceEdgeMap(3,:) = [ 2,7,6,0 ]
    PyramidFaceEdgeMap(4,:) = [ 3,8,7,0 ]
    PyramidFaceEdgeMap(5,:) = [ 4,5,8,0 ]

    WedgeEdgeMap(1,:) = [ 1, 2, 1 ]
    WedgeEdgeMap(2,:) = [ 2, 3, 1 ]
    WedgeEdgeMap(3,:) = [ 1, 3, 1 ]
    WedgeEdgeMap(4,:) = [ 4, 5, 1 ]
    WedgeEdgeMap(5,:) = [ 5, 6, 1 ]
    WedgeEdgeMap(6,:) = [ 6, 4, 1 ]
    WedgeEdgeMap(7,:) = [ 1, 4, 1 ]
    WedgeEdgeMap(8,:) = [ 2, 5, 1 ]
    WedgeEdgeMap(9,:) = [ 3, 6, 1 ]

    WedgeFaceEdgeMap(1,:) = [ 1,2,3,0 ]
    WedgeFaceEdgeMap(2,:) = [ 4,5,6,0 ]
    WedgeFaceEdgeMap(3,:) = [ 1,8,4,7 ]
    WedgeFaceEdgeMap(4,:) = [ 2,9,5,8 ]
    WedgeFaceEdgeMap(5,:) = [ 3,7,6,9 ]

    BrickEdgeMap(1,:) = [ 1, 2,  9 ]
    BrickEdgeMap(2,:) = [ 2, 3,  10 ]
    BrickEdgeMap(3,:) = [ 4, 3,  11 ]
    BrickEdgeMap(4,:) = [ 1, 4,  12 ]
    BrickEdgeMap(5,:) = [ 5, 6,  13 ]
    BrickEdgeMap(6,:) = [ 6, 7,  14 ]
    BrickEdgeMap(7,:) = [ 8, 7,  15 ]
    BrickEdgeMap(8,:) = [ 5, 8,  16 ]
    BrickEdgeMap(9,:) = [ 1, 5,  17 ]
    BrickEdgeMap(10,:) = [ 2, 6, 18 ]
    BrickEdgeMap(11,:) = [ 3, 7, 19 ]
    BrickEdgeMap(12,:) = [ 4, 8, 20 ]

    BrickFaceEdgeMap(1,:) = [ 1,2,3,4   ]
    BrickFaceEdgeMap(2,:) = [ 5,6,7,8   ]    
    BrickFaceEdgeMap(3,:) = [ 1,10,5,9  ]
    BrickFaceEdgeMap(4,:) = [ 2,11,6,10 ]
    BrickFaceEdgeMap(5,:) = [ 3,12,7,11 ]
    BrickFaceEdgeMap(6,:) = [ 4,9,8,12  ]

!
!   Initialize:
!   -----------
    CALL AllocateVector( Mesh % Edges, 12*Mesh % NumberOfBulkElements )
    Edges => Mesh % Edges

    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)
       IF ( .NOT. ASSOCIATED( Element % EdgeIndexes ) ) &
          CALL AllocateVector(Element % EdgeIndexes, Element % TYPE % NumberOfEdges )
       Element % EdgeIndexes = 0
    END DO

    ALLOCATE( HashTable( Mesh % NumberOfNodes ) )
    DO i=1,Mesh % NumberOfNodes
       NULLIFY( HashTable(i) % Head )
    END DO
!------------------------------------------------------------------------------

!   Loop over elements:
!   -------------------
    NofEdges = 0
    DO i=1,Mesh % NumberOfBulkElements
       Element => Mesh % Elements(i)

       ! For P elements mappings are different
       IF ( ASSOCIATED(Element % PDefs) ) THEN
          CALL GetElementEdgeMap( Element, EdgeMap )
          CALL GetElementFaceEdgeMap( Element, FaceEdgeMap ) 
          n = Element % TYPE % NumberOfEdges
       ELSE 
          SELECT CASE( Element % TYPE % ElementCode / 100 )
          CASE(5)
             n = 6
             EdgeMap => TetraEdgeMap
             FaceEdgeMap => TetraFaceEdgeMap
          CASE(6)
             n = 8
             EdgeMap => PyramidEdgeMap
             FaceEdgeMap => PyramidFaceEdgeMap
          CASE(7)
             n = 9
             EdgeMap => WedgeEdgeMap
             FaceEdgeMap => WedgeFaceEdgeMap
          CASE(8)
             n = 12
             EdgeMap => BrickEdgeMap
             FaceEdgeMap => BrickFaceEdgeMap
          CASE DEFAULT
             CYCLE
             WRITE(Message,*) 'Element type',Element % TYPE % ElementCode,'not implemented.' 
             CALL Fatal('FindMeshEdges',Message)
          END SELECT
       END IF

!      Loop over every edge of every element:
!      --------------------------------------
       DO k=1,n

!         Use MIN(Node1,Node2) as key to hash table:
!         ------------------------------------------
          n1 = Element % NodeIndexes(EdgeMap(k,1))
          n2 = Element % NodeIndexes(EdgeMap(k,2))
          IF ( n1 < n2 ) THEN
             Node1 = n1
             Node2 = n2
          ELSE
             Node1 = n2
             Node2 = n1
          END IF
!
!         Look the edge from the hash table:
!         ----------------------------------
          HashPtr => HashTable(Node1) % Head
          Found = .FALSE.
          DO WHILE( ASSOCIATED( HashPtr ) )
             IF ( HashPtr % Node1 == Node2 ) THEN
                Found = .TRUE.
                Edge = HashPtr % Edge
                EXIT
             END IF
             HashPtr => HashPtr % Next
          END DO
!
!         Existing edge, update structures:
!         ---------------------------------
          IF ( Found ) THEN
             Element % EdgeIndexes(k) = Edge

             ! Mark edge as an edge of pydamid square face 
             IF (isPPyramid(Element) .AND. k < 5) THEN
                Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
             END IF

             IF ( ASSOCIATED(Mesh % Faces) ) THEN
               DO ii=1,Element % TYPE % NumberOfFaces
                 Face => Mesh % Faces(Element % FaceIndexes(ii))
                 IF ( .NOT. ASSOCIATED(Face % EdgeIndexes) ) THEN
                   ALLOCATE(Face % EdgeIndexes(Face % TYPE % NumberOfEdges))
                   Face % EdgeIndexes = 0
                 END IF
                 DO jj=1,Face % TYPE % NumberOfEdges
                    IF (FaceEdgeMap(ii,jj) == k) THEN
                       Face % EdgeIndexes(jj) = Edge
                       IF ( .NOT. ASSOCIATED(Edges(Edge) % BoundaryInfo % Left)) THEN
                          Edges(Edge) % BoundaryInfo % Left => Face
                       ELSE
                          Edges(Edge) % BoundaryInfo % Right => Face
                       END IF
                       EXIT
                    END IF
                 END DO
               END DO
             END IF
          ELSE

!            Edge not yet there, create:
!            ---------------------------
             NofEdges = NofEdges + 1
             Edge = NofEdges
             Edges(Edge) % ElementIndex = Edge
             Degree = Element % TYPE % BasisFunctionDegree

!            Edge is always a line segment with deg+1 nodes:
!            -----------------------------------------------
             Edges(Edge) % TYPE => GetElementType( 201 + degree, .FALSE.)

             Edges(Edge) % NDOFs  = 0
             IF (Element % NDOFs /= 0 ) &
                Edges(Edge) % NDOFs  = Edges(Edge) % TYPE % NumberOfNodes
             Edges(Edge) % BDOFs  = 0
             Edges(Edge) % DGDOFs = 0
             Edges(Edge) % EdgeIndexes => NULL()
             Edges(Edge) % FaceIndexes => NULL()

             CALL AllocateVector( Edges(Edge) % NodeIndexes, degree + 1 )
             DO n2=1,degree+1
               Edges(Edge) % NodeIndexes(n2) = &
                    Element % NodeIndexes(EdgeMap(k,n2))
             END DO

             Element % EdgeIndexes(k) = Edge
             ALLOCATE( Edges(Edge) % BoundaryInfo )
             Edges(Edge) % BoundaryInfo % Left  => NULL()
             Edges(Edge) % BoundaryInfo % Right => NULL()

             ! Allocate P element definitions 
             IF ( ASSOCIATED( Element % PDefs ) ) THEN
                CALL AllocatePDefinitions(Edges(Edge))
             
                Edges(Edge) % PDefs % P = 0
                Edges(Edge) % PDefs % pyramidQuadEdge = .FALSE.
                ! Here mark edge as edge of pyramid if needed (or set as not)
                IF (isPPyramid(Element) .AND. k < 5) THEN
                   Edges(Edge) % PDefs % pyramidQuadEdge = .TRUE.
                END IF
             ELSE
                NULLIFY( Edges(Edge) % PDefs )
             END IF

             IF ( ASSOCIATED(Mesh % Faces) ) THEN
               DO ii=1,Element % TYPE % NumberOfFaces
                 Face => Mesh % Faces( Element % FaceIndexes(ii) )
                 IF ( .NOT. ASSOCIATED(Face % EdgeIndexes) ) THEN
                    ALLOCATE( Face % EdgeIndexes( Face % TYPE % NumberOfEdges ) )
                    Face % EdgeIndexes = 0
                 END IF
                 DO jj=1,Face % TYPE % NumberOfEdges
                    IF ( FaceEdgeMap(ii,jj) == k ) THEN
                       Face % EdgeIndexes(jj) = Edge
                       IF (.NOT.ASSOCIATED( Edges(Edge) % BoundaryInfo % Left)) THEN
                          Edges(Edge) % BoundaryInfo % Left => Face
                       ELSE
                          Edges(Edge) % BoundaryInfo % Right => Face
                       END IF
                    END IF
                 END DO
               END DO
             END IF

!            Update the hash table:
!            ----------------------
             ALLOCATE( HashPtr )
             HashPtr % Edge = Edge
             HashPtr % Node1 = Node2
             HashPtr % Next => HashTable(Node1) % Head
             HashTable(Node1) % Head => HashPtr
          END IF
       END DO
    END DO

    Mesh % NumberOfEdges = NofEdges
    CALL Info('FindMeshEdges3D','Number of edges found: '//TRIM(I2S(NofEdges)),Level=10)

!   Delete the hash table:
!   ----------------------
    DO i=1,Mesh % NumberOfNodes
       HashPtr => HashTable(i) % Head
       DO WHILE( ASSOCIATED(HashPtr) )
          HashPtr1 => HashPtr % Next
          DEALLOCATE( HashPtr )
          HashPtr  => HashPtr1
       END DO
    END DO
    DEALLOCATE( HashTable )

    IF (ASSOCIATED(Mesh % Faces)) CALL FixFaceEdges()

    CALL Info('FindMeshEdges3D','All done',Level=12)

CONTAINS 

    SUBROUTINE FixFaceEdges()

      INTEGER :: i,j,k,n,swap,edgeind(4),i1(2),i2(2)

      DO i=1,Mesh % NumberOfFaces
        Face => Mesh % Faces(i)
        n = Face % TYPE % NumberOfEdges
        Edgeind(1:n) = Face % EdgeIndexes(1:n)
        DO j=1,n
          i1 = Mesh % Edges(Edgeind(j)) % NodeIndexes(1:2)
          IF ( i1(1)>i1(2) ) THEN
            swap=i1(1)
            i1(1)=i1(2)
            i1(2)=swap
          END IF
          DO k=1,n
            i2(1) = k
            i2(2) = k+1
            IF ( i2(2)>n ) i2(2)=1
            i2 = Face % NodeIndexes(i2)
            IF ( i2(1)>i2(2) ) THEN
              swap=i2(1)
              i2(1)=i2(2)
              i2(2)=swap
            END IF
            IF ( ALL(i1 == i2) ) THEN
              Face % EdgeIndexes(k) = edgeind(j)
              EXIT
            END IF
          END DO
        END DO
      END DO
    END SUBROUTINE FixFaceEdges
!------------------------------------------------------------------------------
  END SUBROUTINE FindMeshEdges3D
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Finds neighbours of the nodes in given direction.
!> The algorithm finds the neighbour that within 45 degrees of the 
!> given direction has the smallest distance.
!------------------------------------------------------------------------------
  SUBROUTINE FindNeighbourNodes( Mesh,Direction,Neighbours,EndNeighbours)
!------------------------------------------------------------------------------

  TYPE(Mesh_t) , POINTER :: Mesh 
  REAL(KIND=dp) :: Direction(:)
  INTEGER :: Neighbours(:)
  INTEGER, OPTIONAL :: EndNeighbours(:)

  TYPE(Nodes_t) :: ElementNodes
  TYPE(Element_t),POINTER :: CurrentElement
  REAL(KIND=dp), POINTER :: Distances(:)
  REAL(KIND=dp) :: rn(3), rs(3), ss, sn
  INTEGER, POINTER :: NodeIndexes(:)
  INTEGER :: i,j,k,n,t,DIM,istat

  IF(SIZE(Neighbours) < Mesh % NumberOfNodes) THEN
    CALL Warn('FindNeigbourNodes','SIZE of Neighbours should equal Number of Nodes!')
    RETURN
  END IF


  IF(PRESENT(EndNeighbours)) THEN
    IF(SIZE(EndNeighbours) < Mesh % NumberOfNodes) THEN
      CALL Warn('FindNeigbourNodes','SIZE of EndNeigbours should equal Number of Nodes!')
      RETURN
    END IF
  END IF


  DIM = CoordinateSystemDimension()
  N = Mesh % MaxElementNodes

  CALL AllocateVector( ElementNodes % x, n )
  CALL AllocateVector( ElementNodes % y, n )
  CALL AllocateVector( ElementNodes % z, n )
  CALL AllocateVector( Distances, Mesh % NumberOfNodes )

  Neighbours = 0
  Distances = HUGE(Distances)
 
  rn(1:DIM) = Direction(1:DIM)
  ss = SQRT(SUM(rn(1:DIM)**2))
  rn = rn / ss

  DO t=1,Mesh % NumberOfBulkElements

    CurrentElement => Mesh % Elements(t)
    n = CurrentElement % TYPE % NumberOfNodes
    NodeIndexes => CurrentElement % NodeIndexes
  
    ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
    ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
    IF(DIM == 3) THEN
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))
    END IF


    DO i=1,n
      DO j=i+1,n
        rs(1) = ElementNodes % x(j) - ElementNodes % x(i)
        rs(2) = ElementNodes % y(j) - ElementNodes % y(i)
        IF (DIM == 3) THEN
          rs(3) = ElementNodes % z(j) - ElementNodes % z(i)
        END IF
        
        ss = SQRT(SUM(rs(1:DIM)**2))
        sn = SUM(rs(1:DIM)*rn(1:DIM))

        IF(ss < SQRT(2.0) * ABS(sn)) THEN
          IF(sn > 0) THEN
            IF(ss < Distances(NodeIndexes(i))) THEN
              Distances(NodeIndexes(i)) = ss
              Neighbours(NodeIndexes(i)) = NodeIndexes(j)
            END IF
          ELSE
            IF(ss < Distances(NodeIndexes(j))) THEN
              Distances(NodeIndexes(j)) = ss
              Neighbours(NodeIndexes(j)) = NodeIndexes(i)
            END IF
          END IF
        END IF
      END DO
    END DO
  END DO

  ! This loop finds the final neighbour in the end of the chain 
  IF(PRESENT(EndNeighbours)) THEN
    EndNeighbours = Neighbours

    DO t=1,Mesh%NumberOfNodes
      j = Neighbours(t)
      DO WHILE(j /= 0)
        EndNeighbours(t) = j
        j = Neighbours(j)
      END DO
    END DO
  END IF
  DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z, Distances)
!------------------------------------------------------------------------------
END SUBROUTINE FindNeighbourNodes
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE UpdateSolverMesh( Solver, Mesh )
!------------------------------------------------------------------------------
     TYPE( Mesh_t ), POINTER :: Mesh
     TYPE( Solver_t ), TARGET :: Solver
!------------------------------------------------------------------------------
     INTEGER :: i,j,k,n,n1,n2,DOFs
     LOGICAL :: Found, OptimizeBandwidth
     TYPE(Matrix_t), POINTER   :: Matrix
     REAL(KIND=dp), POINTER :: Work(:)
     INTEGER, POINTER :: Permutation(:)
     TYPE(Variable_t), POINTER :: TimeVar, SaveVar, Var
     CHARACTER(LEN=MAX_NAME_LEN) :: str
!------------------------------------------------------------------------------
     SaveVar => Solver % Variable
     DOFs = SaveVar % DOFs

     Solver % Mesh => Mesh
     CALL SetCurrentMesh( CurrentModel, Mesh )
!
!    Create matrix and variable structures for
!    current equation on the new mesh:
!    -----------------------------------------
     Solver % Variable => VariableGet( Mesh % Variables, &
        Solver % Variable % Name, ThisOnly = .FALSE. )

     CALL AllocateVector( Permutation, SIZE(Solver % Variable % Perm) )

     OptimizeBandwidth = ListGetLogical( Solver % Values, 'Optimize Bandwidth', Found )
     IF ( .NOT. Found ) OptimizeBandwidth = .TRUE.

     Matrix => CreateMatrix( CurrentModel, Solver, &
        Mesh, Permutation, DOFs, MATRIX_CRS, OptimizeBandwidth, &
        ListGetString( Solver % Values, 'Equation' ) )

     Matrix % Symmetric = ListGetLogical( Solver % Values, &
             'Linear System Symmetric', Found )

     Matrix % Lumped = ListGetLogical( Solver % Values, &
             'Lumped Mass Matrix', Found )

     ALLOCATE( Work(SIZE(Solver % Variable % Values)) )
     Work = Solver % Variable % Values
     DO k=0,DOFs-1
        DO i=1,SIZE(Permutation)
           IF ( Permutation(i) > 0 ) THEN
              Solver % Variable % Values( DOFs*Permutation(i)-k ) = &
                 Work( DOFs*Solver % Variable % Perm(i)-k )
           END IF
        END DO
     END DO

     IF ( ASSOCIATED( Solver % Variable % PrevValues ) ) THEN
        DO j=1,SIZE(Solver % Variable % PrevValues,2)
           Work = Solver % Variable % PrevValues(:,j)
           DO k=0,DOFs-1
              DO i=1,SIZE(Permutation)
                 IF ( Permutation(i) > 0 ) THEN
                    Solver % Variable % PrevValues( DOFs*Permutation(i) - k,j ) =  &
                        Work( DOFs * Solver % Variable % Perm(i) - k )
                  END IF
              END DO
           END DO
        END DO
     END IF
     DEALLOCATE( Work )

     Solver % Variable % Perm = Permutation
     Solver % Variable % Solver => Solver

     DEALLOCATE( Permutation )
     CALL AllocateVector( Matrix % RHS, Matrix % NumberOfRows )

     IF ( ASSOCIATED(SaveVar % EigenValues) ) THEN
        n = SIZE(SaveVar % EigenValues)

        IF ( n > 0 ) THEN
           Solver % NOFEigenValues = n
           CALL AllocateVector( Solver % Variable % EigenValues,n )
           CALL AllocateArray( Solver % Variable % EigenVectors, n, &
                    SIZE(Solver % Variable % Values) ) 

           IF( Solver % Variable % Dofs > 1 ) THEN
             DO k=1,Solver % Variable % DOFs
               str = ComponentName( Solver % Variable % Name, k )
               Var => VariableGet( Solver % Mesh % Variables, str, .TRUE. )
               IF ( ASSOCIATED( Var ) ) THEN
                 Var % EigenValues => Solver % Variable % EigenValues
                 Var % EigenVectors =>  & 
                     Solver % Variable % EigenVectors(:,k::Solver % Variable % DOFs )
               END IF
             END DO
           END IF
           
           Solver % Variable % EigenValues  = 0.0d0
           Solver % Variable % EigenVectors = 0.0d0

           CALL AllocateVector( Matrix % MassValues, SIZE(Matrix % Values) )
           Matrix % MassValues = 0.0d0
        END IF
     ELSE IF ( ASSOCIATED( Solver % Matrix ) ) THEN
        IF( ASSOCIATED( Solver % Matrix % Force) ) THEN
           n1 = Matrix % NumberOFRows
           n2 = SIZE(Solver % Matrix % Force,2)
           ALLOCATE(Matrix % Force(n1,n2))
           Matrix % Force = 0.0d0
        END IF
     END IF

     Solver % Matrix => Matrix
     Solver % Mesh % Changed = .TRUE.

!------------------------------------------------------------------------------
  END SUBROUTINE UpdateSolverMesh
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!> Split a mesh equally to smaller pieces by performing a uniform split.
!> Also known as mesh multiplication. A 2D element splits into 4 elements of
!> same form, and 3D element into 8 elements. 
!> Currently works only for linear elements.
!------------------------------------------------------------------------------
  FUNCTION SplitMeshEqual(Mesh,h) RESULT( NewMesh )
!------------------------------------------------------------------------------
    REAL(KIND=dp), OPTIONAL :: h(:)
    TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: u(:),v(:),w(:),x(:),y(:),z(:),xh(:)
    INTEGER :: i, j, k, n, NewElCnt, NodeCnt, EdgeCnt, FaceCnt, Node, ParentId, Diag, NodeIt
    LOGICAL :: Found, EdgesPresent
    TYPE(Element_t), POINTER :: Enew,Eold,Edge,Eptr,Eparent,Face,Faces(:)
    INTEGER, POINTER :: Child(:,:)
    INTEGER :: n1,n2,n3,EoldNodes(4),FaceNodes(4),EdgeNodes(2) ! Only linears so far
    INTEGER :: FaceNumber,Edge1,Edge2,Edge3,Edge4,Node12,Node23,Node34,Node41,Node31
    REAL(KIND=dp) :: dxyz(3,3),Dist(3),r,s,t,h1,h2
    TYPE(PElementDefs_t), POINTER :: PDefs
    INTEGER :: ierr, ParTmp(6), ParSizes(6)
    INTEGER, ALLOCATABLE :: FacePerm(:), BulkPerm(:)
!------------------------------------------------------------------------------
    IF ( .NOT. ASSOCIATED( Mesh ) ) RETURN

    CALL Info( 'SplitMeshEqual', 'Mesh splitting works for first order elements 303, 404, 504, (706) and 808.', Level = 6 )

    DO i=1,Mesh % NumberOfBulkElements
      SELECT CASE(Mesh % Elements(i) % TYPE % ElementCode/100)
      CASE(6)
        CALL Fatal('SplitMeshEqual','Pyramids not supported, sorry.')
      END SELECT
    END DO

    NewMesh => AllocateMesh()

    EdgesPresent = ASSOCIATED(Mesh % Edges)
    IF(.NOT.EdgesPresent) CALL FindMeshEdges( Mesh )

    CALL ResetTimer('SplitMeshEqual')

    CALL Info( 'SplitMeshEqual', '******** Old mesh ********', Level = 6 )
    WRITE( Message, * ) 'Nodes             : ',Mesh % NumberOfNodes
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',Mesh % NumberOfBulkElements
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',Mesh % NumberOfBoundaryElements
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Edges             : ',Mesh % NumberOfEdges
    CALL info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Faces             : ',Mesh % NumberOfFaces
    CALL info( 'SplitMeshEqual', Message, Level=6 )
!
!   Update nodal coordinates:
!   -------------------------
    NodeCnt = Mesh % NumberOfNodes + Mesh % NumberOfEdges
!
!   For quad faces add one node in the center:
!   ------------------------
    ALLOCATE(FacePerm(Mesh % NumberOfFaces)); FacePerm = 0
    FaceCnt = 0
    DO i = 1, Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       IF( Face % TYPE % NumberOfNodes == 4 ) THEN
         NodeCnt = NodeCnt+1
         FaceCnt = FaceCnt+1
         FacePerm(i) = NodeCnt
       END IF
    END DO
    
    WRITE( Message, * ) 'Added nodes in the center of faces : ', FaceCnt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )
!
!   For quads and bricks, count centerpoints:
!   -----------------------------------------
    NodeIt = 0
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(4,8)
          NodeCnt = NodeCnt + 1
          NodeIt = NodeIt + 1
       END SELECT
    END DO
    
    WRITE( Message, * ) 'Added nodes in the center of bulks : ', NodeIt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )
!
!   new mesh nodecoordinate arrays:
!   -------------------------------
    CALL AllocateVector( NewMesh % Nodes % x, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % y, NodeCnt )
    CALL AllocateVector( NewMesh % Nodes % z, NodeCnt )

!   shortcuts (u,v,w) old mesh  nodes,
!   (x,y,z) new mesh nodes:
!   ----------------------------------
    u => Mesh % Nodes % x
    v => Mesh % Nodes % y
    w => Mesh % Nodes % z

    x => NewMesh % Nodes % x
    y => NewMesh % Nodes % y
    z => NewMesh % Nodes % z
!
!   new mesh includes old mesh nodes:
!   ----------------------------------
    x(1:Mesh % NumberOfNodes) = u
    y(1:Mesh % NumberOfNodes) = v
    z(1:Mesh % NumberOfNodes) = w

! what is h? - pointer to nodal element size
    IF (PRESENT(h)) THEN
      ALLOCATE(xh(SIZE(x)))
      xh(1:SIZE(h)) = h
    END IF
!
!   add edge centers:
!   -----------------
    j =  Mesh % NumberOfNodes
    DO i=1,Mesh % NumberOfEdges
       j = j + 1
       Edge => Mesh % Edges(i)
       k = Edge % TYPE % NumberOfNodes
       IF (PRESENT(h)) THEN
         h1=h(Edge % NodeIndexes(1))
         h2=h(Edge % NodeIndexes(2))
         r=1._dp/(1+h1/h2)
         x(j) = r*u(Edge%NodeIndexes(1))+(1-r)*u(Edge%NodeIndexes(2))
         y(j) = r*v(Edge%NodeIndexes(1))+(1-r)*v(Edge%NodeIndexes(2))
         z(j) = r*w(Edge%NodeIndexes(1))+(1-r)*w(Edge%NodeIndexes(2))
         xh(j)=r*h1+(1-r)*h2
       ELSE
         x(j) = SUM(u(Edge % NodeIndexes))/k
         y(j) = SUM(v(Edge % NodeIndexes))/k
         z(j) = SUM(w(Edge % NodeIndexes))/k
       END IF
    END DO
    
    CALL Info('SplitMeshEqual','Added edge centers to the nodes list.', Level=10 )  
!
!   add quad face centers for bricks and prisms(wedges):
!   ----------------------------
    j = Mesh % NumberOfNodes + Mesh % NumberOfEdges
    DO i=1,Mesh % NumberOfFaces
       Face => Mesh % Faces(i)
       k = Face % TYPE % NumberOfNodes
       IF( k == 4 ) THEN
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Face % EdgeIndexes(2))
            h2=xh(n+Face % EdgeIndexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Face % EdgeIndexes(3))
            h2=xh(n+Face % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Face,u(Face % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Face,v(Face % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Face,w(Face % NodeIndexes),r,s)
            xh(j) = InterpolateInElement2D(Face,h(Face % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Face % NodeIndexes))/k
            y(j) = SUM(v(Face % NodeIndexes))/k
            z(j) = SUM(w(Face % NodeIndexes))/k
          END IF
       END IF
    END DO
    
    CALL Info('SplitMeshEqual','Added face centers to the nodes list.', Level=10 )
!
!   add centerpoint for quads & bricks:
!   -----------------------------------
    DO i=1,Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(i)
       k = Eold % TYPE % NumberOfNodes
       SELECT CASE( Eold % TYPE % ElementCode / 100 )

       CASE(4)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes
            h1=xh(n+Eold % Edgeindexes(2))
            h2=xh(n+Eold % Edgeindexes(4))
            r=2._dp/(1+h1/h2)-1
            h1=xh(n+Eold % EdgeIndexes(3))
            h2=xh(n+Eold % EdgeIndexes(1))
            s=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement2D(Eold,u(Eold % NodeIndexes),r,s)
            y(j) = InterpolateInElement2D(Eold,v(Eold % NodeIndexes),r,s)
            z(j) = InterpolateInElement2D(Eold,w(Eold % NodeIndexes),r,s)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       CASE(8)
          j = j + 1
          IF (PRESENT(h)) THEN
            n=Mesh % NumberOfNodes+Mesh % NumberOfEdges
            h1=xh(n+Eold % FaceIndexes(4))
            h2=xh(n+Eold % FaceIndexes(6))
            r=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(5))
            h2=xh(n+Eold % FaceIndexes(3))
            s=2._dp/(1+h1/h2)-1

            h1=xh(n+Eold % FaceIndexes(2))
            h2=xh(n+Eold % FaceIndexes(1))
            t=2._dp/(1+h1/h2)-1
            x(j) = InterpolateInElement3D(Eold,u(Eold % NodeIndexes),r,s,t)
            y(j) = InterpolateInElement3D(Eold,v(Eold % NodeIndexes),r,s,t)
            z(j) = InterpolateInElement3D(Eold,w(Eold % NodeIndexes),r,s,t)
          ELSE
            x(j) = SUM(u(Eold % NodeIndexes))/k
            y(j) = SUM(v(Eold % NodeIndexes))/k
            z(j) = SUM(w(Eold % NodeIndexes))/k
          END IF
       END SELECT
    END DO
!
!   Update new mesh node count:
!   ---------------------------
    NewMesh % NumberOfEdges = 0
    NewMesh % NumberOfFaces = 0
    NewMesh % MaxBDOFs = Mesh % MaxBDOFs
    NewMesh % MinEdgeDOFs = Mesh % MinEdgeDOFs
    NewMesh % MinFaceDOFs = Mesh % MinFaceDOFs
    NewMesh % MaxEdgeDOFs = Mesh % MaxEdgeDOFs
    NewMesh % MaxFaceDOFs = Mesh % MaxFaceDOFs
    NewMesh % MaxElementDOFs = Mesh % MaxElementDOFs
    NewMesh % MeshDim = Mesh % MeshDim

    NewMesh % NumberOfNodes = NodeCnt
    NewMesh % Nodes % NumberOfNodes = NodeCnt
!
!   Update bulk elements:
!   =====================
!
!   First count new elements:
!   -------------------------
    NewElCnt = 0
    DO i=1, Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
       Eold => Mesh % Elements(i)
       SELECT CASE( Eold % TYPE % ElementCode/100 )

!      Each element will be divided into 2**Dim new elements:
!      ------------------------------------------------------
       CASE(2)
          NewElCnt = NewElCnt + 2 ! lines
       CASE(3)
          NewElCnt = NewElCnt + 4 ! trias
       CASE(4)
          NewElCnt = NewElCnt + 4 ! quads
       CASE(5)
          NewElCnt = NewElCnt + 8 ! tetras
       CASE(7)
          NewElCnt = NewElCnt + 8 ! prisms (wedges)
       CASE(8)
          NewElCnt = NewElCnt + 8 ! hexas
       END SELECT
    END DO

    WRITE( Message, * ) 'Count of new elements : ', NewElCnt
    CALL Info( 'SplitMeshEqual', Message, Level=10 )

    CALL AllocateVector( NewMesh % Elements, NewElCnt )
    CALL Info('SplitMeshEqual','New mesh allocated.', Level=10 )

    CALL AllocateArray( Child, Mesh % NumberOfBulkElements, 8 )
    CALL Info('SplitMeshEqual','Array for bulk elements allocated.', Level=10 )
    
    NewElCnt = 0
    NodeCnt = Mesh % NumberOfNodes
    EdgeCnt = Mesh % NumberOfEdges

!
!   Index to old quad/hexa centerpoint node in the new mesh nodal arrays:
!   ---------------------------------------------------------------------
    Node = NodeCnt + EdgeCnt + FaceCnt
!
!   Now update all new mesh elements:
!   ---------------------------------
    DO i=1,Mesh % NumberOfBulkElements

       Eold => Mesh % Elements(i)

       SELECT CASE( Eold % TYPE % ElementCode )
       CASE(303)
!
!         Split triangle to four triangles from
!         edge centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,1) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,2) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,3) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Child(i,4) = NewElCnt
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 3)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt

       CASE(404)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split quad to four new quads from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)


       CASE(504)
!
!         Split tetra to 8 new elements from
!         corners and edge centerpoints:
!         ----------------------------------
!
!         1st new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
!
!         2nd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         3rd new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         4th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt

!         Then the annoying part; we still have to split the
!         remaining octahedron into four elements. This can
!         be done in three ways of which only one preserves
!         the minimum angle condition (Delaunay splitting):
!         --------------------------------------------------
          dxyz(1,1) = x(Eold % EdgeIndexes(4) + NodeCnt) &
                    - x(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(2,1) = y(Eold % EdgeIndexes(4) + NodeCnt) &
                    - y(Eold % EdgeIndexes(2) + NodeCnt)
          dxyz(3,1) = z(Eold % EdgeIndexes(4) + NodeCnt) &
                    - z(Eold % EdgeIndexes(2) + NodeCnt)

          dxyz(1,2) = x(Eold % EdgeIndexes(5) + NodeCnt) &
                    - x(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(2,2) = y(Eold % EdgeIndexes(5) + NodeCnt) &
                    - y(Eold % EdgeIndexes(3) + NodeCnt)
          dxyz(3,2) = z(Eold % EdgeIndexes(5) + NodeCnt) &
                    - z(Eold % EdgeIndexes(3) + NodeCnt)

          dxyz(1,3) = x(Eold % EdgeIndexes(6) + NodeCnt) &
                    - x(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(2,3) = y(Eold % EdgeIndexes(6) + NodeCnt) &
                    - y(Eold % EdgeIndexes(1) + NodeCnt)
          dxyz(3,3) = z(Eold % EdgeIndexes(6) + NodeCnt) &
                    - z(Eold % EdgeIndexes(1) + NodeCnt)

          Dist(1) = SQRT( dxyz(1,1)**2 + dxyz(2,1)**2 + dxyz(3,1)**2 )
          Dist(2) = SQRT( dxyz(1,2)**2 + dxyz(2,2)**2 + dxyz(3,2)**2 )
          Dist(3) = SQRT( dxyz(1,3)**2 + dxyz(2,3)**2 + dxyz(3,3)**2 )

          Diag = 1  ! The default diagonal for splitting is between edges 2-4
          IF (Dist(2) < Dist(1) .AND. Dist(2) < Dist(3)) Diag = 2 ! Edges 3-5
          IF (Dist(3) < Dist(1) .AND. Dist(3) < Dist(2)) Diag = 3 ! Edges 1-6

          SELECT CASE( Diag )
          CASE(1)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(2) + NodeCnt
!
          CASE(2)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(5) + NodeCnt
!
          CASE(3)
!
!         5th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         6th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(1) + NodeCnt
!
!         7th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
!
!         8th new element:
!         ----------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 4)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt

          END SELECT


       CASE(706)
!
!         Split prism to 8 new prism from edge
!         centerpoints:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt 
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt 
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt 
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))

!
!         3rd new element (near node 3)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(9) + NodeCnt

!
!         4th new element (bottom center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(5))

!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt

!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt

!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
!
!         8th new element (top half, center)
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 6)
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt



       CASE(808)
!
!         Index to old quad centerpoint node in the
!         new mesh nodal arrays:
!         ------------------------------------------
          Node = Node + 1
!
!         Split brick to 8 new bricks from edge
!         centerpoints and centerpoint of the
!         element:
!         --------------------------------------
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,1) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL  AllocateVector( ENew % NodeIndexes, 8)
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(7) = Node
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(6))
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,2) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(1) + NodeCnt
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(8) = Node
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,3) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(4) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(6) = Node
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(12)+ NodeCnt
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,4) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(1))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(2) + NodeCnt
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(3) + NodeCnt
          Enew % NodeIndexes(5) = Node
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(5))
!
!         5th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,5) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Eold % EdgeIndexes(9) + NodeCnt
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(5) = Eold % NodeIndexes(5)
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(7) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(8) + NodeCnt
!
!         6th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,6) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(3))
          Enew % NodeIndexes(2) = Eold % EdgeIndexes(10)+ NodeCnt
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(4) = Node
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(5) + NodeCnt
          Enew % NodeIndexes(6) = Eold % NodeIndexes(6)
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(8) = FacePerm(Eold % FaceIndexes(2))
!
!         7th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,7) = NewElCnt 
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = FacePerm(Eold % FaceIndexes(6))
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(4) = Eold % EdgeIndexes(12)+ NodeCnt
          Enew % NodeIndexes(5) = Eold % EdgeIndexes(8) + NodeCnt
          Enew % NodeIndexes(6) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(7) = Eold % EdgeIndexes(7) + NodeCnt
          Enew % NodeIndexes(8) = Eold % NodeIndexes(8)
!
!         8th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Child(i,8) = NewElCnt
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( ENew % NodeIndexes, 8 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = FacePerm(Eold % FaceIndexes(4))
          Enew % NodeIndexes(3) = Eold % EdgeIndexes(11)+ NodeCnt
          Enew % NodeIndexes(4) = FacePerm(Eold % FaceIndexes(5))
          Enew % NodeIndexes(5) = FacePerm(Eold % FaceIndexes(2))
          Enew % NodeIndexes(6) = Eold % EdgeIndexes(6) + NodeCnt
          Enew % NodeIndexes(7) = Eold % NodeIndexes(7)
          Enew % NodeIndexes(8) = Eold % EdgeIndexes(7) + NodeCnt

       CASE DEFAULT
          WRITE( Message,* ) 'Element type ', Eold % TYPE % ElementCode, &
              ' not supprted by the multigrid solver.'
          CALL Fatal( 'SplitMeshEqual', Message )
       END SELECT
    END DO

!
!   Update new mesh element counts:
!   -------------------------------
    NewMesh % NumberOfBulkElements = NewElCnt

!
!   Update boundary elements:
!   NOTE: Internal boundaries not taken care of...:!!!!
!   ---------------------------------------------------
    DO i=1,Mesh % NumberOfBoundaryElements

       j = i + Mesh % NumberOfBulkElements
       Eold => Mesh % Elements(j)
!
!      get parent of the boundary element:
!      -----------------------------------
       Eparent => Eold % BoundaryInfo % Left
       IF ( .NOT.ASSOCIATED(Eparent) ) &
          eParent => Eold % BoundaryInfo % Right
       IF ( .NOT. ASSOCIATED( Eparent ) ) CYCLE

       ParentId = Eparent % ElementIndex

       SELECT CASE( Eold % TYPE % ElementCode / 100 )
       CASE(2)
!
!         Line segments:
!         ==============
!
!         which edge of the parent element are we ?
!         -----------------------------------------
          DO Edge1=1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             IF ( Eold % NodeIndexes(1) == Edge % NodeIndexes(1) .AND. &
                  Eold % NodeIndexes(2) == Edge % NodeIndexes(2) .OR.  &
                  Eold % NodeIndexes(2) == Edge % NodeIndexes(1) .AND. &
                  Eold % NodeIndexes(1) == Edge % NodeIndexes(2) ) EXIT
          END DO
!
!         index of the old edge centerpoint in the
!         new mesh nodal arrays:
!         ----------------------------------------
          Node = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,4
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found = .FALSE.
             DO k=1,n-1
                IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k+1) .OR.  &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(1) == Eptr % NodeIndexes(k+1) ) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END DO
             IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(1) .OR.  &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(1) == Eptr % NodeIndexes(1) ) THEN
                Found = .TRUE.
             END IF
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 2 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,4
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             Found = .FALSE.
             DO k=1,n-1
                IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k+1) .OR.  &
                     Enew % NodeIndexes(2) == Eptr % NodeIndexes(k)   .AND. &
                     Enew % NodeIndexes(1) == Eptr % NodeIndexes(k+1) ) THEN
                   Found = .TRUE.
                   EXIT
                END IF
             END DO
             IF ( Enew % NodeIndexes(1) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(1) .OR.  &
                  Enew % NodeIndexes(2) == Eptr % NodeIndexes(n) .AND. &
                  Enew % NodeIndexes(1) == Eptr % NodeIndexes(1) ) THEN
                Found = .TRUE.
             END IF
             IF ( Found ) EXIT
          END DO
          Enew % BoundaryInfo % Left => Eptr

       CASE(3)
!
!         Trias:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:3) = Eold % NodeIndexes(1:3)
          CALL sort( 3, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:3) = Face % NodeIndexes(1:3)
             CALL sort( 3, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) ) EXIT

          END DO
!
!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(1) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node31 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Node31
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 3 )
          Enew % NodeIndexes(1) = Node31
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,3
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr

       CASE(4)
!
!         Quads:
!         ======
!
!         On which face of the parent element are we ?
!         --------------------------------------------
          EoldNodes(1:4) = Eold % NodeIndexes(1:4)
          CALL sort( 4, EoldNodes )

          DO FaceNumber = 1, SIZE( Eparent % FaceIndexes )
             Face => Mesh % Faces( Eparent % FaceIndexes(FaceNumber) )
             FaceNodes(1:4) = Face % NodeIndexes(1:4)
             CALL sort( 4, FaceNodes )

             IF ( EoldNodes(1) == FaceNodes(1) .AND. &
                  EoldNodes(2) == FaceNodes(2) .AND. &
                  EoldNodes(3) == FaceNodes(3) .AND. &
                  EoldNodes(4) == FaceNodes(4) ) EXIT

          END DO

!         Then, what are the edges on this face?
!         --------------------------------------
!
!         First edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(1), Eold % NodeIndexes(2) )
          DO Edge1 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge1) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Second edge:
!         ------------
          EoldNodes(1) = MIN( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(2), Eold % NodeIndexes(3) )
          DO Edge2 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge2) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Third edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(3), Eold % NodeIndexes(4) )
          DO Edge3 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge3) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO

!         Fourth edge:
!         -----------
          EoldNodes(1) = MIN( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          EoldNodes(2) = MAX( Eold % NodeIndexes(4), Eold % NodeIndexes(1) )
          DO Edge4 = 1,SIZE(Eparent % EdgeIndexes)
             Edge => Mesh % Edges( Eparent % EdgeIndexes(Edge4) )
             EdgeNodes(1) = MIN( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             EdgeNodes(2) = MAX( Edge % NodeIndexes(1), Edge % NodeIndexes(2) )
             IF ( EoldNodes(1) == EdgeNodes(1) .AND. &
                  EoldNodes(2) == EdgeNodes(2) ) EXIT
          END DO
!
!         index of the old face and edge centerpoints
!         in the new mesh nodal arrays:
!         ----------------------------------------
          Node = FacePerm(Eparent % FaceIndexes(FaceNumber)) ! faces mid-point
          Node12 = Eparent % EdgeIndexes(Edge1) + Mesh % NumberOfNodes
          Node23 = Eparent % EdgeIndexes(Edge2) + Mesh % NumberOfNodes
          Node34 = Eparent % EdgeIndexes(Edge3) + Mesh % NumberOfNodes
          Node41 = Eparent % EdgeIndexes(Edge4) + Mesh % NumberOfNodes
!
!         1st new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Eold % NodeIndexes(1)
          Enew % NodeIndexes(2) = Node12
          Enew % NodeIndexes(3) = Node
          Enew % NodeIndexes(4) = Node41
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 )  CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         2nd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node12
          Enew % NodeIndexes(2) = Eold % NodeIndexes(2)
          Enew % NodeIndexes(3) = Node23
          Enew % NodeIndexes(4) = Node
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         3rd new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node41
          Enew % NodeIndexes(2) = Node
          Enew % NodeIndexes(3) = Node34
          Enew % NodeIndexes(4) = Eold % NodeIndexes(4)
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
!
!         4th new element
!         ---------------
          NewElCnt = NewElCnt + 1
          Enew => NewMesh % Elements(NewElCnt)
          Enew = Eold
          Enew % ElementIndex = NewElCnt
          CALL AllocateVector( Enew % NodeIndexes, 4 )
          Enew % NodeIndexes(1) = Node
          Enew % NodeIndexes(2) = Node23
          Enew % NodeIndexes(3) = Eold % NodeIndexes(3)
          Enew % NodeIndexes(4) = Node34
          ALLOCATE( Enew % BoundaryInfo )
          Enew % BoundaryInfo = Eold % BoundaryInfo
          NULLIFY( Enew % BoundaryInfo % Left )
          NULLIFY( Enew % BoundaryInfo % Right )
!
!         Search the new mesh parent element among the
!         children of the old mesh parent element:
!         --------------------------------------------
          DO j=1,8
             Eptr => NewMesh % Elements( Child(ParentId,j) )
             n = Eptr % TYPE % NumberOfNodes
             n3 = 0 ! Count matches (metodo stupido)
             DO n1 = 1,4
                DO n2 = 1,SIZE(Eptr % NodeIndexes)
                   IF( Enew % NodeIndexes(n1) == Eptr % NodeIndexes(n2) ) n3 = n3+1
                END DO
             END DO
             IF ( n3 > 2 ) EXIT
          END DO
          IF( n3 < 3 ) CALL Error( 'SplitMeshEqual', 'Parent element not found' )
          Enew % BoundaryInfo % Left => Eptr
       END SELECT
    END DO

!
!   Update new mesh boundary element counts:
!   ----------------------------------------
    NewMesh % NumberOfBoundaryElements = NewElCnt - &
            NewMesh % NumberOfBulkElements
    NewMesh % MaxElementDOFs  = Mesh % MaxElementDOFs
    NewMesh % MaxElementNodes = Mesh % MaxElementNodes

    j = 0
    DO i=1,NewMesh % NumberOfBulkElements+NewMesh % NumberOfBoundaryElements
      Enew => NewMesh % Elements(i)

      IF ( Enew % DGDOFs>0 ) THEN
        ALLOCATE(Enew % DGIndexes(Enew % DGDOFs))
        DO k=1,Enew % DGDOFs
          j = j + 1
          Enew % DGIndexes(k)=j
        END DO
      ELSE
        Enew % DGIndexes=>NULL()
      END IF

      IF (i<=NewMesh % NumberOfBulkElements) THEN
         PDefs => Enew % PDefs

         IF(ASSOCIATED(PDefs)) THEN
           CALL AllocatePDefinitions(Enew)
           Enew % PDefs = PDefs

           ! All elements in actual mesh are not edges
           Enew % PDefs % pyramidQuadEdge = .FALSE.
           Enew % PDefs % isEdge = .FALSE.

           ! If element is of type tetrahedron and is a p element,
           ! do the Ainsworth & Coyle trick
           IF (Enew % TYPE % ElementCode == 504) CALL ConvertToACTetra(Enew)
            CALL GetRefPElementNodes( Enew % Type,  Enew % Type % NodeU, &
                 Enew % Type % NodeV, Enew % Type % NodeW )
         END IF
      ELSE
        Enew % PDefs=>NULL()
      END IF
      Enew % EdgeIndexes => NULL()
      Enew % FaceIndexes => NULL()
      Enew % BubbleIndexes => NULL()
    END DO

    CALL Info( 'SplitMeshEqual', '******** New mesh ********', Level=6 )
    WRITE( Message, * ) 'Nodes             : ',NewMesh % NumberOfNodes
    CALL Info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Bulk elements     : ',NewMesh % NumberOfBulkElements
    CALL Info( 'SplitMeshEqual', Message, Level=6 )
    WRITE( Message, * ) 'Boundary elements : ',NewMesh % NumberOfBoundaryElements
    CALL Info( 'SplitMeshEqual', Message, Level=6 )


    ! Information of the new system size, also in parallel
    !----------------------------------------------------------------------
    ParTmp(1) = Mesh % NumberOfNodes
    ParTmp(2) = Mesh % NumberOfBulkElements
    ParTmp(3) = Mesh % NumberOfBoundaryElements
    ParTmp(4) = NewMesh % NumberOfNodes
    ParTmp(5) = NewMesh % NumberOfBulkElements
    ParTmp(6) = NewMesh % NumberOfBoundaryElements

    IF( .FALSE. .AND. ParEnv % PEs > 1 ) THEN
      CALL MPI_ALLREDUCE(ParTmp,ParSizes,6,MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)

      CALL Info('SplitMeshEqual','Information on parallel mesh sizes')
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(1),' nodes'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(2),' bulk elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'Initial mesh has ',ParSizes(3),' boundary elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(4),' nodes'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(5),' bulk elements'
      CALL Info('SplitMeshEqual',Message)
      WRITE ( Message,'(A,I0,A)') 'New mesh has ',ParSizes(6),' boundary elements'
      CALL Info('SplitMeshEqual',Message)
    END IF


    CALL CheckTimer('SplitMeshEqual',Delete=.TRUE.)

!
!   Update structures needed for parallel execution:
!   ------------------------------------------------
    CALL UpdateParallelMesh( Mesh, NewMesh )
!
!
!   Finalize:
!   ---------
    DEALLOCATE( Child )
    IF(.NOT.EdgesPresent) THEN
      CALL ReleaseMeshEdgeTables( Mesh )
      CALL ReleaseMeshFaceTables( Mesh )
    ELSE
      CALL FindMeshEdges( NewMesh )
    END IF

!call writemeshtodisk( NewMesh, "." )
!stop
CONTAINS

!------------------------------------------------------------------------------
    SUBROUTINE UpdateParallelMesh( Mesh, NewMesh )
!------------------------------------------------------------------------------
       TYPE(Mesh_t), POINTER :: Mesh, NewMesh
!------------------------------------------------------------------------------
       TYPE(Element_t), POINTER :: Edge, Face, Element, BoundaryElement
       INTEGER :: i,j,k,l,m,n,p,q, istat
       INTEGER, POINTER :: IntCnts(:),IntArray(:),Reorder(:)
       INTEGER, ALLOCATABLE :: list1(:), list2(:)
       LOGICAL, ALLOCATABLE :: InterfaceTag(:)

       INTEGER :: jedges
       LOGICAL :: Found
!------------------------------------------------------------------------------

       IF ( ParEnv % PEs <= 1 ) RETURN
!
!      Update mesh interfaces for parallel execution.
!      ==============================================
!
!      Try to get an agreement about the  global numbering
!      of new mesh nodes among set of processes solving
!      this specific eq. Also allocate and generate
!      all other control information needed in parallel
!      execution:
!      ----------------------------------------------------
       n = NewMesh % NumberOfNodes
       ALLOCATE( NewMesh % ParallelInfo % NeighbourList(n), stat=istat )
       IF ( istat /= 0 ) &
         CALL Fatal( 'UpdateParallelMesh', 'Allocate error.' )
       CALL AllocateVector( NewMesh % ParallelInfo % INTERFACE,n  )
       CALL AllocateVector( NewMesh % ParallelInfo % GlobalDOFs,n )

       DO i=1,n
          NULLIFY( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours )
       END DO

       n = Mesh % NumberOfNodes
       NewMesh % ParallelInfo % INTERFACE = .FALSE.
       NewMesh % ParallelInfo % INTERFACE(1:n) = Mesh % ParallelInfo % INTERFACE

       NewMesh % ParallelInfo % GlobalDOFs = 0
       NewMesh % ParallelInfo % GlobalDOFs(1:n) = &
          Mesh % ParallelInfo % GlobalDOFs
!
!      My theory is, that a new node will be an
!      interface node only if all the edge or face
!      nodes which contribute to its existence are
!      interface nodes (the code immediately below
!      will only count sizes):
!      -------------------------------------------
!

       ! New version based on edges and faces (2. March 2007):
       !=====================================================
       SELECT CASE( CoordinateSystemDimension() )
          
       CASE(2)
          !
          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % INTERFACE(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'
          !
          ! Determine possible interface edges:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfEdges ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfEdges
             Edge => Mesh % Edges(i)
             IF( ASSOCIATED(Edge % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Edge % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % INTERFACE( Edge % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          !
          ! Eliminate false positives based on BoundaryElement -data:
          !----------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements( Mesh % NumberOfBulkElements + i )
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED( Element ) ) &
                  Element => BoundaryElement % BoundaryInfo % Right
             IF( .NOT.ASSOCIATED( Element ) ) CYCLE
             IF( .NOT.ASSOCIATED( Element % EdgeIndexes ) ) CYCLE
             
             ALLOCATE( list1( SIZE( BoundaryElement % NodeIndexes )))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort( SIZE(list1), list1 )
             
             DO j = 1,Element % TYPE % NumberOfEdges
                k = Element % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                IF( SIZE( Edge % NodeIndexes ) /= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Edge % NodeIndexes )))
                list2 = Edge % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO

                DEALLOCATE(list2)
                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfEdges
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Edge => Mesh % Edges(i)
             
             ! This is just for the edge count:
             !---------------------------------
             IF( NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + i) ) CYCLE
             
             ! Mark interface nodes and count edges:
             !--------------------------------------
             NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + i) = .TRUE.
             p = p+1

          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 2*p ! check
          
       CASE(3)

          ! Count interface nodes:
          !-----------------------
          p = 0 
          DO i = 1, Mesh % NumberOfNodes
             IF( Mesh % ParallelInfo % INTERFACE(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface nodes'

          ! Determine possible interface faces:
          !------------------------------------
          ALLOCATE( InterfaceTag( Mesh % NumberOfFaces ) )
          InterfaceTag = .FALSE.
          DO i = 1,Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( ASSOCIATED(Face % BoundaryInfo % Left) .AND. &
                  ASSOCIATED(Face % BoundaryInfo % Right) ) CYCLE
             IF( .NOT.ALL( Mesh % ParallelInfo % INTERFACE( Face % NodeIndexes ) )) CYCLE
             InterfaceTag(i) = .TRUE.
          END DO
          
          ! Eliminate false interface faces based on BoundaryElement -data:
          !----------------------------------------------------------------
          DO i = 1,Mesh % NumberOfBoundaryElements
             BoundaryElement => Mesh % Elements(Mesh % NumberOfBulkElements+i)
             Element => BoundaryElement % BoundaryInfo % Left
             IF( .NOT.ASSOCIATED(Element) ) &
                Element => BoundaryElement % BoundaryInfo % Right
              IF( .NOT.ASSOCIATED(Element) ) CYCLE
              IF( .NOT.ASSOCIATED(Element % FaceIndexes) ) CYCLE
             
             ALLOCATE(list1(SIZE(BoundaryElement % NodeIndexes)))
             list1 = BoundaryElement % NodeIndexes
             CALL Sort(SIZE(list1),list1)
             
             DO j = 1,Element % TYPE % NumberOfFaces
                k = Element % FaceIndexes(j)
                Face => Mesh % Faces(k)
                IF(SIZE(Face % NodeIndexes)/= SIZE(list1) ) CYCLE
                
                ALLOCATE( list2( SIZE( Face % NodeIndexes )))
                list2 = Face % NodeIndexes
                CALL Sort( SIZE(list2), list2 )

                Found = .TRUE.
                DO l = 1,SIZE(list2)
                   Found = Found .AND. ( list1(l)==list2(l) )
                END DO
                
                DEALLOCATE(list2)

                IF( Found ) InterfaceTag(k) = .FALSE.
             END DO

             DEALLOCATE(list1)
          END DO
          
          ! Count interface faces:
          !-----------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             Face => Mesh % Faces(i)
             IF( InterfaceTag(i) ) p = p+1
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface faces'
          
          ! Mark all new interface nodes and count interface edges:
          !--------------------------------------------------------
          p = 0
          DO i = 1, Mesh % NumberOfFaces
             IF( .NOT. InterfaceTag(i) ) CYCLE
             Face => Mesh % Faces(i)
             
             DO j = 1,SIZE( Face % EdgeIndexes )
                k = Face % EdgeIndexes(j)
                Edge => Mesh % Edges(k)
                
                ! This is just for the edge count:
                !---------------------------------
                IF( NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + k) ) CYCLE
                
                ! Mark interface nodes and count edges:
                !--------------------------------------
                NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes + k) = .TRUE.
                p = p+1
             END DO
          END DO
!         WRITE(*,'(A,I4,A,I6,A)')'SplitMeshEqual: PE:', &
!              Parenv % MyPE+1, ' Found',p,' interface edges'
          
          DEALLOCATE( InterfaceTag )

          j = p
          k = 3*p ! check
          
       END SELECT

!======================================================================================================
       j = p
       jedges = p

!      For bricks, check also the faces:
!      ---------------------------------
       DO i = 1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i) 
          IF( Face % TYPE % NumberOfNodes == 4 ) THEN
             IF ( ALL( Mesh % ParallelInfo % INTERFACE( Face % NodeIndexes ) ) ) THEN
                NewMesh % ParallelInfo % INTERFACE( Mesh % NumberOfNodes &
                     + Mesh % NumberOfEdges + i ) = .TRUE.
                j = j + 1
                k = k + Face % TYPE % NumberOfNodes
             END IF
          END IF
       END DO

!      CALL AllocateVector( IntCnts,  j )
!      CALL AllocateVector( IntArray, k )
!
!      Old mesh nodes were copied as is...
!      -----------------------------------
       DO i=1,Mesh % NumberOfNodes
          CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours, &
                SIZE( Mesh % ParallelInfo % Neighbourlist(i) % Neighbours) )

          NewMesh % ParallelInfo % NeighbourList(i) % Neighbours = &
             Mesh % ParallelInfo % NeighbourList(i) % Neighbours
       END DO
!
!      Take care of the new mesh internal nodes.
!      Parallel global numbering will take care
!      of the interface nodes:
!      ----------------------------------------
       DO i=Mesh % NumberOfNodes+1, NewMesh % NumberOfNodes
          IF ( .NOT. NewMesh % ParallelInfo % INTERFACE(i) ) THEN
            CALL AllocateVector( NewMesh % ParallelInfo % NeighbourList(i) % Neighbours,1 )
            NewMesh % ParallelInfo % NeighbourList(i) %  Neighbours(1) = ParEnv % MyPE
          END IF
       END DO
!
!      Copy global indices of edge and/or face nodes
!      to temporary work arrays:
!      ---------------------------------------------
!
! check also this:
!      j = 0
!      k = 0
!      DO i = 1,Mesh % NumberOfEdges
!         Edge => Mesh % Edges(i)
!         
!         ! Added check for parent elements 25.2.2007:
!         Found = .NOT.( ASSOCIATED(edge % boundaryinfo % left) &
!              .AND.  ASSOCIATED(edge % boundaryinfo % right) )
!         
!         IF ( ALL(Mesh % ParallelInfo % INTERFACE(Edge % NodeIndexes)) .AND. Found ) THEN
!            j = j + 1
!            IntCnts(j) = Edge % TYPE % NumberOfNodes
!            IntArray( k+1:k+IntCnts(j) ) = &
!                 Mesh % Parallelinfo % GlobalDOFs(Edge % NodeIndexes)
!            CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!            k = k + IntCnts(j)
!         END IF
!      END DO
!      !
!      ! For bricks, check also the faces:
!      ! ---------------------------------
!      DO i = 1,Mesh % NumberOfFaces
!         Face => Mesh % Faces(i)
!         IF( Face % TYPE % NumberOfNodes == 4 ) THEN
!            IF ( ALL( Mesh % ParallelInfo % INTERFACE(Face % NodeIndexes) ) ) THEN
!               j = j + 1
!               IntCnts(j) = Face % TYPE % NumberOfNodes
!               IntArray(k+1:k+IntCnts(j)) = &
!                    Mesh % ParallelInfo % GlobalDOFs(Face % NodeIndexes)
!               CALL Sort( IntCnts(j), IntArray(k+1:k+IntCnts(j)) )
!               k = k + IntCnts(j)
!            END IF
!         END IF
!      END DO
!
!      Finally the beef, do the exchange of new
!      interfaces. The parallel global numbering
!      subroutine will also do reordering of the
!      nodes, hence the reorder array:
!      -------------------------------------------
       CALL AllocateVector( Reorder, NewMesh % NumberOfNodes )
       Reorder = [ (i, i=1,NewMesh % NumberOfNodes) ]

       k = NewMesh % Nodes % NumberOfNodes - Mesh % Nodes % NumberOfNodes

       CALL ParallelGlobalNumbering( NewMesh, Mesh, k, IntCnts, IntArray, Reorder )

!      Account for the reordering of the nodes:
!      ----------------------------------------
       DO i=1,NewMesh % NumberOfBulkElements + &
            NewMesh % NumberOfBoundaryElements
          NewMesh % Elements(i) % NodeIndexes = &
              Reorder( NewMesh % Elements(i) % NodeIndexes )
       END DO

!      DEALLOCATE( IntCnts, IntArray, Reorder )
!      DEALLOCATE( Reorder )
!------------------------------------------------------------------------------
    END SUBROUTINE UpdateParallelMesh
  END FUNCTION SplitMeshEqual
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMesh( Mesh )
!------------------------------------------------------------------------------
     TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
     TYPE(Projector_t), POINTER :: Projector
     TYPE(Projector_t), POINTER :: Projector1
     TYPE(Variable_t), POINTER  :: Var, Var1
     INTEGER :: i,j,k
     LOGICAL :: GotIt
     REAL(KIND=dp), POINTER :: ptr(:)
!------------------------------------------------------------------------------
 
!    Deallocate mesh variables:
!    --------------------------


     CALL Info('ReleaseMesh','Releasing mesh variables',Level=15)
     CALL ReleaseVariableList( Mesh % Variables )
     Mesh % Variables => NULL()

!    Deallocate mesh geometry (nodes,elements and edges):
!    ----------------------------------------------------
     IF ( ASSOCIATED( Mesh % Nodes ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh nodes',Level=15)
       IF ( ASSOCIATED( Mesh % Nodes % x ) ) DEALLOCATE( Mesh % Nodes % x )
       IF ( ASSOCIATED( Mesh % Nodes % y ) ) DEALLOCATE( Mesh % Nodes % y )
       IF ( ASSOCIATED( Mesh % Nodes % z ) ) DEALLOCATE( Mesh % Nodes % z )
       DEALLOCATE( Mesh % Nodes )

       IF ( ASSOCIATED( Mesh % ParallelInfo % GlobalDOFs ) ) &
           DEALLOCATE( Mesh % ParallelInfo % GlobalDOFs )

       IF ( ASSOCIATED( Mesh % ParallelInfo % NeighbourList ) ) THEN 
         DO i=1,Mesh % NumberOfNodes
           IF(ASSOCIATED( Mesh % ParallelInfo % NeighbourList(i) % Neighbours ) ) &
               DEALLOCATE( Mesh % ParallelInfo % NeighbourList(i) % Neighbours )
         END DO
         DEALLOCATE( Mesh % ParallelInfo % NeighbourList )
       END IF

       IF ( ASSOCIATED( Mesh % ParallelInfo % INTERFACE ) ) &
           DEALLOCATE( Mesh % ParallelInfo % INTERFACE )
     END IF

     Mesh % Nodes => NULL()

     IF ( ASSOCIATED( Mesh % Edges ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh edges',Level=15)
       CALL ReleaseMeshEdgeTables( Mesh )
       Mesh % Edges => NULL()
     END IF

     IF ( ASSOCIATED( Mesh % Faces ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh faces',Level=15)
       CALL ReleaseMeshFaceTables( Mesh )
       Mesh % Faces => NULL()
     END IF

     IF (ASSOCIATED(Mesh % ViewFactors) ) THEN
     CALL Info('ReleaseMesh','Releasing mesh view factors',Level=15)
       CALL ReleaseMeshFactorTables( Mesh % ViewFactors )
       Mesh % ViewFactors => NULL()
     END IF


!    Deallocate mesh to mesh projector structures:
!    ---------------------------------------------
     Projector => Mesh % Projector
     DO WHILE( ASSOCIATED( Projector ) )
       CALL Info('ReleaseMesh','Releasing mesh projector',Level=15)
       CALL FreeMatrix( Projector % Matrix )
       CALL FreeMatrix( Projector % TMatrix )
       Projector1 => Projector
       Projector => Projector % Next
       DEALLOCATE( Projector1 )
     END DO
     Mesh % Projector => NULL()


!    Deallocate quadrant tree (used in mesh to mesh interpolation):
!    --------------------------------------------------------------
     IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh quadrant tree',Level=15)
       CALL FreeQuadrantTree( Mesh % RootQuadrant )
       Mesh % RootQuadrant => NULL()
     END IF


     IF ( ASSOCIATED( Mesh % Elements ) ) THEN
       CALL Info('ReleaseMesh','Releasing mesh elements',Level=15)

        DO i=1,Mesh % NumberOfBulkElements+Mesh % NumberOfBoundaryElements

!          Boundaryinfo structure for boundary elements
!          ---------------------------------------------
           IF ( Mesh % Elements(i) % Copy ) CYCLE

           IF ( i > Mesh % NumberOfBulkElements ) THEN
             IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo ) ) THEN
               IF (ASSOCIATED(Mesh % Elements(i) % BoundaryInfo % GebhardtFactors)) THEN
                 IF ( ASSOCIATED( Mesh % Elements(i) % BoundaryInfo % &
                     GebhardtFactors % Elements ) ) THEN
                   DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                       GebhardtFactors % Elements )
                   DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % &
                       GebhardtFactors % Factors )
                 END IF
                 DEALLOCATE( Mesh % Elements(i) % BoundaryInfo % GebhardtFactors )
               END IF
               DEALLOCATE( Mesh % Elements(i) % BoundaryInfo )
             END IF
           END IF

           IF ( ASSOCIATED( Mesh % Elements(i) % NodeIndexes ) ) &
               DEALLOCATE( Mesh % Elements(i) % NodeIndexes )
           Mesh % Elements(i) % NodeIndexes => NULL()
           
           IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) &
              DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
           Mesh % Elements(i) % EdgeIndexes => NULL()

           IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) &
              DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
           Mesh % Elements(i) % FaceIndexes => NULL()

           IF ( ASSOCIATED( Mesh % Elements(i) % DGIndexes ) ) &
              DEALLOCATE( Mesh % Elements(i) % DGIndexes )
           Mesh % Elements(i) % DGIndexes => NULL()

           IF ( ASSOCIATED( Mesh % Elements(i) % BubbleIndexes ) ) &
             DEALLOCATE( Mesh % Elements(i) % BubbleIndexes )
           Mesh % Elements(i) % BubbleIndexes => NULL()

           ! This creates problems later on!!!
           !IF ( ASSOCIATED( Mesh % Elements(i) % PDefs ) ) &
           !   DEALLOCATE( Mesh % Elements(i) % PDefs )

           Mesh % Elements(i) % PDefs => NULL()
 
        END DO
        DEALLOCATE( Mesh % Elements )
        Mesh % Elements => NULL()
      END IF

      CALL Info('ReleaseMesh','Releasing mesh finished',Level=15)
     
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshEdgeTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Edge
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Edges ) ) THEN
       DO i=1,Mesh % NumberOfEdges
          Edge => Mesh % Edges(i)
          IF ( ASSOCIATED( Edge % NodeIndexes ) ) THEN
             DEALLOCATE( Edge % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Edge % BoundaryInfo ) ) THEN
             DEALLOCATE( Edge % BoundaryInfo )
          END IF
       END DO

       DEALLOCATE( Mesh % Edges )
    END IF
    NULLIFY( Mesh % Edges )
    Mesh % NumberOfEdges = 0

    DO i=1,Mesh % NumberOfBulkElements
       IF ( ASSOCIATED( Mesh % Elements(i) % EdgeIndexes ) ) THEN
          DEALLOCATE( Mesh % Elements(i) % EdgeIndexes )
          NULLIFY( Mesh % Elements(i) % EdgeIndexes )
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshEdgeTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFaceTables( Mesh )
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    INTEGER :: i
    TYPE(Element_t), POINTER :: Face
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Mesh % Faces ) ) THEN
       DO i=1,Mesh % NumberOfFaces
          Face => Mesh % Faces(i)
          IF ( ASSOCIATED( Face % NodeIndexes ) ) THEN
             DEALLOCATE( Face % NodeIndexes )
          END IF
          IF ( ASSOCIATED( Face % BoundaryInfo ) ) THEN
             DEALLOCATE( Face % BoundaryInfo )
          END IF
       END DO

       DEALLOCATE( Mesh % Faces )
    END IF
    NULLIFY( Mesh % Faces )
    Mesh % NumberOfFaces = 0

    DO i=1,Mesh % NumberOfBulkElements
       IF ( ASSOCIATED( Mesh % Elements(i) % FaceIndexes ) ) THEN
          DEALLOCATE( Mesh % Elements(i) % FaceIndexes )
          NULLIFY( Mesh % Elements(i) % FaceIndexes )
       END IF
    END DO
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFaceTables
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
  SUBROUTINE ReleaseMeshFactorTables( Factors )
!------------------------------------------------------------------------------
    TYPE(Factors_t), POINTER :: Factors(:)
!------------------------------------------------------------------------------
    INTEGER :: i
!------------------------------------------------------------------------------
    IF ( ASSOCIATED( Factors ) ) THEN
       DO i=1,SIZE( Factors)
          IF (ASSOCIATED(Factors(i) % Factors))  DEALLOCATE(Factors(i) % Factors)
          IF (ASSOCIATED(Factors(i) % Elements)) DEALLOCATE(Factors(i) % Elements)
       END DO
       DEALLOCATE(  Factors )
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE ReleaseMeshFactorTables
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
  SUBROUTINE SetCurrentMesh( Model, Mesh )
!------------------------------------------------------------------------------
    TYPE(Model_t) :: Model
    TYPE(Mesh_t),  POINTER :: Mesh
!------------------------------------------------------------------------------
    Model % Variables => Mesh % Variables

    Model % Mesh  => Mesh
    Model % Nodes => Mesh % Nodes
    Model % NumberOfNodes = Mesh % NumberOfNodes
    Model % Nodes % NumberOfNodes = Mesh % NumberOfNodes

    Model % Elements => Mesh % Elements
    Model % MaxElementNodes = Mesh % MaxElementNodes
    Model % NumberOfBulkElements = Mesh % NumberOfBulkElements
    Model % NumberOfBoundaryElements = Mesh % NumberOfBoundaryElements
!------------------------------------------------------------------------------
  END SUBROUTINE SetCurrentMesh
!------------------------------------------------------------------------------


!----------------------------------------------------------------------------------
  SUBROUTINE DisplaceMesh( Mesh, Update, SIGN, Perm, DOFs, StabRecomp, UpdateDirs )
!----------------------------------------------------------------------------------
    TYPE(Mesh_t) , POINTER :: Mesh 
    REAL(KIND=dp) :: Update(:)
    INTEGER :: DOFs,SIGN,Perm(:)
    LOGICAL, OPTIONAL :: StabRecomp
    INTEGER, OPTIONAL :: UpdateDirs

    INTEGER :: i,k,dim
    LOGICAL :: StabFlag

    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element

    IF ( PRESENT( UpdateDirs ) ) THEN
      dim = UpdateDirs
    ELSE
      dim = DOFs
    END IF

    DO i=1,MIN( SIZE(Perm), SIZE(Mesh % Nodes % x) )
       k = Perm(i)
       IF ( k > 0 ) THEN
         k = DOFs * (k-1)
         Mesh % Nodes % x(i)   = Mesh % Nodes % x(i) + SIGN * Update(k+1)
         IF ( dim > 1 ) &
           Mesh % Nodes % y(i) = Mesh % Nodes % y(i) + SIGN * Update(k+2)
         IF ( dim > 2 ) &
           Mesh % Nodes % z(i) = Mesh % Nodes % z(i) + SIGN * Update(k+3)
        END IF
    END DO

    StabFlag = .TRUE.
    IF ( PRESENT( StabRecomp ) ) StabFlag = StabRecomp

    IF ( SIGN == 1 .AND. StabFlag ) THEN
       k = Mesh % MaxElementDOFs
       CALL AllocateVector( ElementNodes % x,k )
       CALL AllocateVector( ElementNodes % y,k )
       CALL AllocateVector( ElementNodes % z,k )

       DO i=1,Mesh % NumberOfBulkElements
          Element => Mesh % Elements(i)
          IF ( ANY( Perm( Element % NodeIndexes ) == 0 ) ) CYCLE

          k = Element % TYPE % NumberOfNodes
          ElementNodes % x(1:k) = Mesh % Nodes % x(Element % NodeIndexes)
          ElementNodes % y(1:k) = Mesh % Nodes % y(Element % NodeIndexes)
          ElementNodes % z(1:k) = Mesh % Nodes % z(Element % NodeIndexes)
          IF ( Mesh % Stabilize ) THEN
             CALL StabParam( Element,ElementNodes,k, &
                          Element % StabilizationMk, Element % Hk )
          ELSE
             Element % hK = ElementDiameter( Element, ElementNodes )
          END IF
       END DO

       DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z)
    END IF
!------------------------------------------------------------------------------
  END SUBROUTINE DisplaceMesh
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>  Convert tetrahedral element to Ainsworth & Coyle type tetrahedron.
!------------------------------------------------------------------------------
  SUBROUTINE ConvertToACTetra( Tetra )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getTetraEdgeMap, getTetraFaceMap
    IMPLICIT NONE
    
    TYPE(Element_t), POINTER :: Tetra  !< Tetrahedral element to convert
!------------------------------------------------------------------------------
    INTEGER :: i, globalMin, globalMax, globalMinI
    INTEGER, DIMENSION(3) :: face, globalFace
    INTRINSIC MIN, MAX, CSHIFT

    ! Sanity check
    IF (Tetra % TYPE % ElementCode /= 504 .OR. &
         .NOT. ASSOCIATED(Tetra % PDefs)) THEN
       CALL Warn('MeshUtils::ConvertToACTetra','Element to convert not p tetrahedron!')
       RETURN
    END IF    
   
    ! Find global min and max vertices
    globalMin = Tetra % NodeIndexes(1)
    globalMinI = 1
    globalMax = Tetra % NodeIndexes(1)
    DO i=2,4
       ! Find min
       IF (globalMin > Tetra % NodeIndexes(i)) THEN
          globalMin = Tetra % NodeIndexes(i)
          globalMinI = i
       ELSE IF (globalMax < Tetra % NodeIndexes(i)) THEN
          globalMax = Tetra % NodeIndexes(i)
       END IF
    END DO
    
    ! Get face containing global min (either face 1 or 2)
    IF (globalMinI == 4) THEN
       face = getTetraFaceMap(2)
    ELSE
       face = getTetraFaceMap(1)
    END IF
    globalFace(1:3) = Tetra % NodeIndexes(face)

    ! Rotate face until first local index is min global
    DO 
       ! Check if first node matches global min node
       IF (globalMin == globalFace(1)) EXIT
       
       globalFace(1:3) = CSHIFT(globalFace,1)
    END DO
    ! Assign new local numbering
    Tetra % NodeIndexes(face) = globalFace(1:3)

    ! Face 3 now contains global max
    face = getTetraFaceMap(3)
    globalFace(1:3) = Tetra % NodeIndexes(face)
    ! Rotate face until last local index is max global
    DO
       ! Check if last node matches global max node
       IF (globalMax == globalFace(3)) EXIT

       globalFace(1:3) = CSHIFT(globalFace,1)
    END DO
    ! Assign new local numbering
    Tetra % NodeIndexes(face) = globalFace(1:3)

    ! Set AC tetra type
    IF (Tetra % NodeIndexes(2) < Tetra % NodeIndexes(3)) THEN
       Tetra % PDefs % TetraType = 1
    ELSE IF (Tetra % NodeIndexes(3) < Tetra % NodeIndexes(2)) THEN
       Tetra % PDefs % TetraType = 2
    ELSE 
       CALL Fatal('MeshUtils::ConvertToACTetra','Corrupt element type')
    END IF
   
  END SUBROUTINE ConvertToACTetra


!------------------------------------------------------------------------------
!>     Assign local number of edge to given boundary element. Also copies all 
!>     p element attributes from element edge to boundary edge.
!------------------------------------------------------------------------------
  SUBROUTINE AssignLocalNumber( EdgeElement, Element, Mesh,NoPE )
!------------------------------------------------------------------------------
    USE PElementMaps, ONLY : getFaceEdgeMap 
    IMPLICIT NONE

    ! Parameters
    TYPE(Mesh_t) :: Mesh            !< Finite element mesh containing faces and edges.
    TYPE(Element_t), POINTER :: EdgeElement  !< Edge element to which assign local number
    TYPE(Element_t), POINTER :: Element      !< Bulk element with some global numbering to use to assign local number
    LOGICAL, OPTIONAL :: NoPE
!------------------------------------------------------------------------------
    ! Local variables

    INTEGER i,j,n,edgeNumber, numEdges, bMap(4)
    TYPE(Element_t), POINTER :: Edge
    LOGICAL :: EvalPE

    EvalPE = .TRUE.
    IF(PRESENT(NoPE)) EvalPE = .NOT.NoPE
    
    ! Get number of points, edges or faces
    numEdges = 0
    SELECT CASE (Element % TYPE % DIMENSION)
    CASE (1)
      RETURN
    CASE (2)
       numEdges = Element % TYPE % NumberOfEdges
    CASE (3)   
       numEdges = Element % TYPE % NumberOfFaces
    CASE DEFAULT
       WRITE (*,*) 'MeshUtils::AssignLocalNumber, Unsupported dimension:', Element % TYPE % DIMENSION
       RETURN
    END SELECT

    ! For each edge or face in element try to find local number
    DO edgeNumber=1, numEdges
       ! If edges have not been created, stop search. This should not happen, actually.
       IF (.NOT. ASSOCIATED(Element % EdgeIndexes)) THEN
          ! EdgeElement % localNumber = 0
          RETURN
       END IF

       Edge => GetElementEntity(Element,edgeNumber,Mesh)

       ! Edge element not found. This should not be possible, unless there
       ! is an error in the mesh read in process..
       IF (.NOT. ASSOCIATED(Edge)) THEN
          CALL Warn('MeshUtils::AssignLocalNumber','Edge element not found')
          ! EdgeElement % localNumber = 0
          RETURN
       END IF

       n = 0
       ! For each element node
       DO i=1, Edge % TYPE % NumberOfNodes
          ! For each node in edge element
          DO j=1, EdgeElement % TYPE % NumberOfNodes
             ! If edge and edgeelement node match increment counter
             IF (Edge % NodeIndexes(i) == EdgeElement % NodeIndexes(j)) n = n + 1
          END DO
       END DO

       ! If all nodes are on boundary, edge was found
       IF (n == EdgeElement % TYPE % NumberOfNodes) THEN
          IF(EvalPE) &
              EdgeElement % PDefs % localNumber = edgeNumber

          ! Change ordering of global nodes to match that of element
          bMap = getElementBoundaryMap( Element, edgeNumber )
          DO j=1,n
          	EdgeElement % NodeIndexes(j) = Element % NodeIndexes(bMap(j))
	  END DO

          ! Copy attributes of edge element to boundary element
          ! Misc attributes
          IF(EvalPE) THEN
            EdgeElement % PDefs % isEdge = Edge % PDefs % isEdge
          
          ! Gauss points
            EdgeElement % PDefs % GaussPoints = Edge % PDefs % GaussPoints

          ! Element p
            EdgeElement % PDefs % P = Edge % PDefs % P
          END IF
          
          !(and boundary bubble dofs)
          EdgeElement % BDOFs = Edge % BDOFs


          ! If this boundary has edges copy edge indexes
          IF (ASSOCIATED(Edge % EdgeIndexes)) THEN
             ! Allocate element edges to element
             n = Edge % TYPE % NumberOfEdges
             bmap(1:4) = getFaceEdgeMap( Element, edgeNumber )
             
             IF ( ASSOCIATED( EdgeElement % EdgeIndexes) ) THEN
                DEALLOCATE( EdgeElement % EdgeIndexes )
             END IF
             
             CALL AllocateVector( EdgeElement % EdgeIndexes, n )
             ! Copy edges from edge to boundary edge
             DO i=1,n
                EdgeElement % EdgeIndexes(i) = Element % EdgeIndexes(bmap(i))
             !    EdgeElement % EdgeIndexes(i) = Element % EdgeIndexes(i)
             END DO
          END IF
          
          ! Edge fields copied and local edge found so return
          RETURN
       END IF
    END DO

    ! If we are here local number not found
    CALL Warn('MeshUtils::AssignLocalNumber','Unable to find local edge')
    ! EdgeElement % localNumber = 1
  CONTAINS

    FUNCTION GetElementEntity(Element, which, Mesh) RESULT(Entity)
      IMPLICIT NONE

      TYPE(Element_t), POINTER :: Element, Entity 
      INTEGER :: which
      TYPE(Mesh_t) :: Mesh

      NULLIFY(Entity)
      ! Switch by element dimension
      SELECT CASE (Element % TYPE % DIMENSION)
         CASE (2)
            Entity => Mesh % Edges( Element % EdgeIndexes(which))
         CASE (3)
            Entity => Mesh % Faces( Element % FaceIndexes(which))
         CASE DEFAULT
            WRITE (*,*) 'AssignLocalNumber::GetElementEntity: Unsupported dimension'
            RETURN
      END SELECT
    END FUNCTION GetElementEntity
  END SUBROUTINE AssignLocalNumber
    

!------------------------------------------------------------------------------
!>     Based on element degrees of freedom, return the sum of element
!>     degrees of freedom.
!------------------------------------------------------------------------------
  FUNCTION getElementMaxDOFs( Mesh, Element ) RESULT(dofs)
!------------------------------------------------------------------------------
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh        !< Finite element mesh
    TYPE(Element_t), POINTER :: Element  !< Element to get maximum dofs for
    INTEGER :: dofs                      !< maximum number of dofs for Element
!------------------------------------------------------------------------------

    TYPE(ELement_t), POINTER :: Edge, Face
    INTEGER :: i, edgeDofs, faceDofs
    
    ! Get sum of edge dofs if any
    edgeDofs = 0
    IF (ASSOCIATED(Element % EdgeIndexes)) THEN
       DO i=1, Element % TYPE % NumberOfEdges
          Edge => Mesh % Edges(Element % EdgeIndexes(i))
          edgeDofs = edgeDofs + Edge % BDOFs
       END DO
    END IF

    ! Get sum of face dofs if any
    faceDofs = 0
    IF (ASSOCIATED(Element % FaceIndexes)) THEN
       DO i=1, Element % TYPE % NumberOfFaces
          Face => Mesh % Faces(Element % FaceIndexes(i))
          faceDofs = faceDofs + Face % BDOFs
       END DO
    END IF

    ! Get sum of all dofs in element
    dofs = Element % TYPE % NumberOfNodes + &
         edgeDofs + faceDofs + Element % BDOFs
  END FUNCTION getElementMaxDOFs




!------------------------------------------------------------------------------
!> Creates a permutation table for bodies or boundaries using a free chosen string
!> as mask. The resulting permutation is optimized in order, if requested. The
!> subroutine is intended to help in saving boundary data in an ordered manner,
!> but it can find other uses as well. Currently the implementation is limited
!> to normal Lagrangian elements.
!------------------------------------------------------------------------------
  SUBROUTINE MakePermUsingMask( Model,Solver,Mesh,MaskName, &
       OptimizeBW, Perm, LocalNodes, MaskOnBulk, RequireLogical, ParallelComm )
!------------------------------------------------------------------------------
    TYPE(Model_t)  :: Model
    TYPE(Mesh_t)   :: Mesh
    TYPE(SOlver_t) :: Solver
    INTEGER :: LocalNodes
    LOGICAL :: OptimizeBW
    INTEGER, POINTER :: Perm(:)
    CHARACTER(LEN=*) :: MaskName
    LOGICAL, OPTIONAL :: MaskOnBulk
    LOGICAL, OPTIONAL :: RequireLogical
    LOGICAL, OPTIONAL :: ParallelComm
!------------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm(:), Neighbours(:)
    INTEGER, ALLOCATABLE :: s_e(:,:), r_e(:), fneigh(:), ineigh(:)
    TYPE(ListMatrix_t), POINTER :: ListMatrix(:)
    INTEGER :: t,i,j,k,l,m,k1,k2,n,p,q,e1,e2,f1,f2,This,bf_id,nn,ii(ParEnv % PEs)
    INTEGER :: ierr, status(MPI_STATUS_SIZE), NewDofs
    LOGICAL :: Flag, Found, FirstRound, MaskIsLogical, Hit, Parallel
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)
    INTEGER :: Indexes(30), ElemStart, ElemFin, Width
    TYPE(ListMatrixEntry_t), POINTER :: CList, Lptr
    TYPE(Element_t), POINTER :: CurrentElement,Elm
    REAL(KIND=dp) :: MinDist, Dist
!------------------------------------------------------------------------------

    IF(PRESENT(ParallelComm)) THEN
      Parallel = ParallelComm
    ELSE
      Parallel = ParEnv % PEs > 1
    END IF

    ! First check if there are active elements for this mask
    IF( PRESENT( MaskOnBulk ) ) MaskOnBulk = .FALSE.
    IF( PRESENT( RequireLogical ) ) THEN
      MaskIsLogical = RequireLogical
    ELSE
      MaskIsLogical = .FALSE.
    END IF

    IF(.NOT. ASSOCIATED( Perm ) ) THEN
      ALLOCATE( Perm( Mesh % NumberOfNodes ) )
      Perm = 0
    END IF

    ElemStart = HUGE(ElemStart) 
    ElemFin = 0     
    DO l = 1, Model % NumberOfBodyForces
       IF( MaskIsLogical ) THEN
         Hit = ListGetLogical( Model % BodyForces(l) % Values,MaskName,Found) 
       ELSE
         Hit = ListCheckPresent( Model % BodyForces(l) % Values,MaskName)
       END IF 
       IF( Hit ) THEN
          ElemStart = 1
          ElemFin = Mesh % NumberOfBulkElements
          IF( PRESENT( MaskOnBulk ) ) MaskOnBulk = .TRUE.
          EXIT
       END IF
    END DO
    DO l = 1, Model % NumberOfBCs
       IF( MaskIsLogical ) THEN
         Hit = ListGetLogical(Model % BCs(l) % Values,MaskName,Found )
       ELSE
         Hit = ListCheckPresent(Model % BCs(l) % Values,MaskName )
       END IF
       IF( Hit ) THEN
          ElemStart = MIN( ElemStart, Mesh % NumberOfBulkElements + 1)
          ElemFin = Mesh % NumberOfBulkElements + Mesh % NumberOFBoundaryElements
          EXIT
       END IF
    END DO

    IF( ElemFin - ElemStart <= 0) THEN
       LocalNodes = 0
       RETURN
    END IF

    k = 0
    Perm = 0
    FirstRound = .TRUE.

    ! Loop over the active elements
    ! 1st round initial numbering is given
    ! 2nd round a list matrix giving all the connections is created

100 DO t=ElemStart, ElemFin
       
       CurrentElement => Mesh % Elements(t)
       
       Hit = .FALSE.
       IF(t <= Mesh % NumberOfBulkElements) THEN
          l = CurrentElement % BodyId
	  bf_id = ListGetInteger( Model % Bodies(l) % Values, 'Body Force',Found)
	  IF( bf_id>0 ) THEN
            IF( MaskIsLogical ) THEN
              Hit = ListGetLogical( Model % BodyForces(bf_id) % Values, MaskName, Found )
            ELSE
              Hit = ListCheckPresent( Model % BodyForces(bf_id) % Values, MaskName )
            END IF
	  END IF 
       ELSE
          DO l=1, Model % NumberOfBCs
            IF ( Model % BCs(l) % Tag /= CurrentElement % BoundaryInfo % Constraint ) CYCLE
            IF( MaskIsLogical ) THEN
              Hit = ListGetLogical(Model % BCs(l) % Values,MaskName, Found ) 
            ELSE
              Hit = ListCheckPresent(Model % BCs(l) % Values,MaskName ) 
            END IF
            EXIT
          END DO
       END IF       
       IF( .NOT. Hit ) CYCLE       
       
       n = CurrentElement % TYPE % NumberOfNodes
       Indexes(1:n) = CurrentElement % NodeIndexes(1:n)
       
       IF( FirstRound ) THEN
          DO i=1,n
             j = Indexes(i)
             IF ( Perm(j) == 0 ) THEN
                k = k + 1
                Perm(j) = k
             END IF
          END DO
       ELSE
          DO i=1,n
             k1 = Perm(Indexes(i))
             IF ( k1 <= 0 ) CYCLE
             DO j=1,n
                k2 = Perm(Indexes(j))
                IF ( k2 <= 0 ) CYCLE
                Lptr => List_GetMatrixIndex( ListMatrix,k1,k2 )
             END DO
          END DO
       END IF
    END DO
    LocalNodes = k

    !In parallel case, detect nodes which are shared with another partition
    !which may not have an element on this boundary
    !Code borrowed from CommunicateLinearSystemTag
    IF( Parallel ) THEN

      ALLOCATE( IsNeighbour(ParEnv % PEs), fneigh(ParEnv % PEs), ineigh(ParEnv % PEs) )

      nn = MeshNeighbours(Mesh, IsNeighbour)
      nn = 0
      ineigh = 0
      DO i=0, ParEnv % PEs-1
        k = i+1
        IF(i==ParEnv % myPE) CYCLE
        IF(.NOT. IsNeighbour(k) ) CYCLE
        nn = nn + 1
        fneigh(nn) = k
        ineigh(k) = nn
      END DO

      n = COUNT(Perm > 0 .AND. Mesh % ParallelInfo % Interface)
      ALLOCATE( s_e(n, nn ), r_e(n) )

      CALL CheckBuffer( nn*3*n )

      ii = 0
      DO i=1, Mesh % NumberOfNodes
        IF(Perm(i) > 0 .AND. Mesh % ParallelInfo % Interface(i) ) THEN
          DO j=1,SIZE(Mesh % ParallelInfo % Neighbourlist(i) % Neighbours)
            k = Mesh % ParallelInfo % Neighbourlist(i) % Neighbours(j)
            IF ( k == ParEnv % MyPE ) CYCLE
            k = k + 1
            k = ineigh(k)
            IF ( k> 0) THEN
              ii(k) = ii(k) + 1
              s_e(ii(k),k) = Mesh % ParallelInfo % GlobalDOFs(i)
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

      NewDofs = 0

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
            k = SearchNode( Mesh % ParallelInfo, r_e(j), Order=Mesh % ParallelInfo % Gorder )
            IF ( k>0 ) THEN
              IF(.NOT. Perm(k) > 0) THEN
                NewDofs = NewDofs + 1
                Perm(k) = LocalNodes + NewDofs
              END IF
            END IF
          END DO
        END IF
      END DO
      DEALLOCATE(s_e, r_e )

      LocalNodes = LocalNodes + NewDofs
    END IF

    ! Don't optimize bandwidth for parallel cases
    IF( Parallel .OR. .NOT. OptimizeBW ) RETURN

    IF(FirstRound) THEN
       ! Allocate space 
       NULLIFY( ListMatrix )
       ListMatrix => List_AllocateMatrix(LocalNodes)
       FirstRound = .FALSE.

       ! Find the node in the lower left corner at give it the 1st index
       ! since it will probably determine the 1st index
       MinDist = HUGE(MinDist)
       DO i=1,SIZE(Perm)
          IF( Perm(i) <= 0) CYCLE
          Dist = Mesh % Nodes % x(i) + Mesh % Nodes % y(i) + Mesh % Nodes % z(i)
          IF(Dist < MinDist) THEN
             MinDist = Dist
             j = i
          END IF
       END DO

       ! Find the 1st node and swap it with the lower corner
       DO i=1,SIZE(Perm)
          IF( Perm(i) == 1) EXIT
       END DO       
       Perm(i) = Perm(j)
       Perm(j) = 1

       GOTO 100
    END IF

!------------------------------------------------------------------------------

    ALLOCATE( InvPerm(LocalNodes) )
    InvPerm = 0
    DO i=1,SIZE(Perm)
       IF (Perm(i)>0) InvPerm(Perm(i)) = i
    END DO

    ! The bandwidth optimization for lines results to perfectly ordered 
    ! permutations. If there is only one line the 1st node should be the 
    ! lower left corner.

    Flag = .TRUE.
    Width = OptimizeBandwidth( ListMatrix, Perm, InvPerm, &

         LocalNodes, Flag, Flag, MaskName )

    ! We really only need the permutation, as there will be no matrix equation
    ! associated with it.
    DEALLOCATE( InvPerm )
    CALL List_FreeMatrix( LocalNodes, ListMatrix )

!------------------------------------------------------------------------------
  END SUBROUTINE MakePermUsingMask
!------------------------------------------------------------------------------




!------------------------------------------------------------------------
!> Find a point in the mesh structure
!> There are two strategies:
!> 1) Recursive where the same routine is repeated with sloppier criteria
!> 2) One-sweep strategy where the best hit is registered and used if of 
!>    acceptable accuracy. 
!> There are two different epsilons that control the search. One for the 
!> rough test in absolute coordinates and another one for the more accurate
!> test in local coordinates.   
!-------------------------------------------------------------------------
  FUNCTION PointInMesh(Solver, GlobalCoords, LocalCoords, HitElement, &
      CandElement, ExtInitialize ) RESULT ( Hit )
        
    TYPE(Solver_t) :: Solver
    REAL(KIND=dp) :: GlobalCoords(3), LocalCoords(3)
    TYPE(Element_t), POINTER :: HitElement 
    TYPE(Element_t), POINTER, OPTIONAL :: CandElement
    LOGICAL, OPTIONAL :: ExtInitialize
    LOGICAL :: Hit
!-------------------------------------------------------------------------
    LOGICAL :: Initialize, Allocated = .FALSE., Stat, DummySearch, &
        MaskExists, Found, IsRecursive
    INTEGER :: i,j,k,n,bf_id,dim,mini
    REAL(KIND=dp) :: u,v,w,dist,mindist,MinLocalCoords(3)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: CurrentElement
    TYPE(Quadrant_t), POINTER, SAVE :: RootQuadrant =>NULL(), LeafQuadrant
    REAL(kind=dp) :: BoundingBox(6), eps2, eps1 = 1d-3, GlobalEps, LocalEps
    CHARACTER(LEN=MAX_NAME_LEN) :: MaskName


    SAVE :: Allocated, ElementNodes, DummySearch, Mesh, MaskName, MaskExists, &
        GlobalEps, LocalEps, IsRecursive


    IF( PRESENT( ExtInitialize ) ) THEN
      Initialize = ExtInitialize
    ELSE
      Initialize = .NOT. Allocated 
    END IF

    IF( Initialize ) THEN
      Mesh => Solver % Mesh
      n = Mesh % MaxElementNodes
      IF( Allocated ) THEN
        DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
      END IF
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n))
      Allocated = .TRUE.

      IsRecursive = ListGetLogical( CurrentModel % Simulation,&
          'Interpolation Search Recursive',Stat )
!      IF(.NOT. Stat ) IsRecursive = .TRUE.

      LocalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Local Epsilon', Stat )
      IF(.NOT. stat) LocalEps = 1.0d-10

      GlobalEps = ListGetConstReal( CurrentModel % Simulation,  &
          'Interpolation Global Epsilon', Stat ) 
      IF(.NOT. stat) THEN
        IF( IsRecursive ) THEN
          GlobalEps = 2.0d-10
        ELSE
          GlobalEps = 1.0d-4
        END IF
      END IF

      DummySearch = ListGetLogical( CurrentModel % Simulation,&
          'Interpolation Search Dummy',Stat )

      MaskName = ListGetString( CurrentModel % Simulation,&
          'Interpolation Search Mask',MaskExists )

      IF( ASSOCIATED( Mesh % RootQuadrant ) ) THEN
        CALL FreeQuadrantTree( Mesh % RootQuadrant )
        Mesh % RootQuadrant => NULL()
      END IF
    END IF
      

    !-----------------------------------------------
    ! Create the octree search structure, if needed 
    !-----------------------------------------------
    IF ( .NOT. ( DummySearch .OR.  ASSOCIATED( Mesh % RootQuadrant ) ) ) THEN
      BoundingBox(1) = MINVAL( Mesh % Nodes % x )
      BoundingBox(2) = MINVAL( Mesh % Nodes % y )
      BoundingBox(3) = MINVAL( Mesh % Nodes % z )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % x )
      BoundingBox(5) = MAXVAL( Mesh % Nodes % y )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
      
      eps2 = eps1 * MAXVAL( BoundingBox(4:6) - BoundingBox(1:3) )
      BoundingBox(1:3) = BoundingBox(1:3) - eps2
      BoundingBox(4:6) = BoundingBox(4:6) + eps2
      
      CALL BuildQuadrantTree( Mesh,BoundingBox,Mesh % RootQuadrant)
      RootQuadrant => Mesh % RootQuadrant
      IF (.NOT. ASSOCIATED(RootQuadrant) ) THEN
        Hit = .FALSE.
        CALL Warn('PointInMesh','No RootQuadrant associated')
        RETURN
      END IF
    END IF


    Hit = .FALSE.

    ! Check that the previous hit is not hit even now
    !-------------------------------------------------
    IF( PRESENT( CandElement ) ) THEN

      IF( ASSOCIATED(CandElement)) THEN

        CurrentElement => CandElement
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        IF ( PointInElement( CurrentElement, ElementNodes, &
            GlobalCoords, LocalCoords ) ) THEN
          Hit = .TRUE.
          HitElement => CurrentElement
          RETURN
        END IF
      END IF
    END IF


    Eps1 = GlobalEps
    Eps2 = LocalEps


100 IF( DummySearch ) THEN

      mindist = HUGE( mindist ) 
      
      !----------------------------------------------------------
      ! Go through all bulk elements in a dummy search.
      ! This algorithm is mainly here for debugging purposes, or
      ! if just a few nodes need to be searched.
      !----------------------------------------------------------
      DO k=1,Mesh % NumberOfBulkElements
        CurrentElement => Mesh % Elements(k)
        n = CurrentElement % TYPE % NumberOfNodes
        NodeIndexes => CurrentElement % NodeIndexes
        
        IF( MaskExists ) THEN
          bf_id = ListGetInteger( CurrentModel % Bodies(CurrentElement % BodyId) % Values, &
              'Body Force', Found )
          IF( .NOT. Found ) CYCLE
          IF(.NOT. ListCheckPresent( CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
        END IF

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
        
        Hit = PointInElement( CurrentElement, ElementNodes, &
            GlobalCoords, LocalCoords, Eps1, Eps2, LocalDistance = dist )
        IF( dist < mindist ) THEN
          mini = k
          mindist = dist
        END IF
        IF( Hit ) EXIT
      END DO      
    ELSE
      !-----------------------------------------------
      ! Find the right element using an octree search
      ! This is the preferred algorithms of the two.
      !-----------------------------------------------
      NULLIFY(CurrentElement)
      CALL FindLeafElements(GlobalCoords, Mesh % MeshDim, RootQuadrant, LeafQuadrant)
      IF ( ASSOCIATED(LeafQuadrant) ) THEN
        DO j=1, LeafQuadrant % NElemsInQuadrant
          k = LeafQuadrant % Elements(j)
          CurrentElement => Mesh % Elements(k)
          
          IF( MaskExists ) THEN
            bf_id = ListGetInteger( CurrentModel % Bodies(CurrentElement % BodyId) % Values, &
                'Body Force', Found )
            IF( .NOT. Found ) CYCLE
            IF(.NOT. ListCheckPresent( CurrentModel % BodyForces(bf_id) % Values,MaskName) ) CYCLE
          END IF
          
          n = CurrentElement % TYPE % NumberOfNodes
          NodeIndexes => CurrentElement % NodeIndexes
                    
          ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
          ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
          ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
          
          Hit = PointInElement( CurrentElement, ElementNodes, &
              GlobalCoords, LocalCoords, Eps1, Eps2, LocalDistance = dist ) 
          IF( dist < mindist ) THEN
            mini = k
            mindist = dist
            MinLocalCoords = LocalCoords
          END IF
          IF( Hit ) EXIT
        END DO
      END IF      
    END IF

    IF( .NOT. Hit ) THEN
      IF( IsRecursive ) THEN
        Eps1 = 10.0 * Eps1
        Eps2 = 10.0 * Eps2
        IF( Eps1 <= 1.0_dp ) GOTO 100
      ELSE
        IF( mindist < Eps1 ) THEN
          CurrentElement => Mesh % Elements(k)
          LocalCoords = MinLocalCoords
          Hit = .TRUE.
        END IF
      END IF
    END IF

    IF( Hit ) HitElement => CurrentElement
    
  END FUNCTION PointInMesh

  
  !----------------------------------------------------------------
  !> Maps coordinates from the original nodes into a new coordinate
  !> system while optionally maintaining the original coordinates. 
  !> Note that this may be called 
  !---------------------------------------------------------------
  SUBROUTINE CoordinateTransformation( Mesh, CoordTransform, Params, &
      IrreversibleTransformation )
    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform
    TYPE(ValueList_t), POINTER :: Params
    LOGICAL, OPTIONAL :: IrreversibleTransformation
    !---------------------------------------------------------------   
    REAL(KIND=dp) :: R0(3),R1(3),Coeff,Rad0
    LOGICAL :: Irreversible,FirstTime,Reuse,UpdateNodes,Found
    REAL(KIND=dp), POINTER :: x0(:),y0(:),z0(:),x1(:),y1(:),z1(:)
    REAL(KIND=dp), POINTER CONTIG :: NewCoords(:)
    INTEGER :: i,j,k,n,Mode
    TYPE(Variable_t), POINTER :: Var

    ! The coordinate transformation may either be global for all the solvers
    ! and this overrides the original nodes permanently. 
    ! Or it can be a solver specific transformation which saves the initial 
    ! coordinates. 
    CALL Info('CoordinateTransformation','Starting')

    IF(.NOT. ASSOCIATED(Mesh) ) THEN
      CALL Fatal('CoordinateTransformation','Mesh not associated!')
    END IF

    IF( PRESENT( IrreversibleTransformation ) ) THEN
      Irreversible = IrreversibleTransformation
    ELSE
      Irreversible = .FALSE.
    END IF

    n = Mesh % NumberOfNodes 

    x0 => Mesh % Nodes % x
    y0 => Mesh % Nodes % y
    z0 => Mesh % Nodes % z
    
    IF( Irreversible ) THEN
      UpdateNodes = .TRUE.
      ! Map to the same nodes
      x1 => Mesh % Nodes % x
      y1 => Mesh % Nodes % y
      z1 => Mesh % Nodes % z
    ELSE
      ReUse = ListGetLogical(Params,'Coordinate Transformation Reuse',Found ) 
      FirstTime = .NOT. ASSOCIATED( Mesh % NodesMapped )
      IF( FirstTime ) THEN
        ALLOCATE( Mesh % NodesMapped )
        NULLIFY( NewCoords )
        ALLOCATE( NewCoords(3*n) )
        NewCoords = 0.0_dp
        Mesh % NodesMapped % x => NewCoords(1:n)
        Mesh % NodesMapped % y => NewCoords(n+1:2*n)
        Mesh % NodesMapped % z => NewCoords(2*n+1:3*n)
        ! Mesh % NodesMapped % x => NewCoords(1::3)
        ! Mesh % NodesMapped % y => NewCoords(2::3)
        ! Mesh % NodesMapped % z => NewCoords(3::3)
      ELSE
        IF( n /= SIZE(Mesh % NodesMapped % x) ) THEN
          CALL Fatal('CoordinateTransformation','Sizes of original and mapped mesh differ!')
        END IF
      END IF

      IF( CoordTransform == 'previous' ) THEN
        IF( FirstTime ) THEN
          CALL Fatal('CoordinateTransformation','One cannot reuse unexisting transformation!')
        END IF
        ReUse = .TRUE.
      END IF

      ! Note that if many solvers reutilize the same coordinates then they must 
      ! also have the same coordinate mapping. 
      !------------------------------------------------------------------------
      UpdateNodes = FirstTime .OR. .NOT. ReUse 
      ! Map different nodes if the original ones are kept
      x1 => Mesh % NodesMapped % x
      y1 => Mesh % NodesMapped % y
      z1 => Mesh % NodesMapped % z      

      IF( FirstTime ) THEN
        IF( ListGetLogical(Params,'Coordinate Transformation Save',Found ) ) THEN
          CALL Info('CoordinateTranformation',&
              'Creating variables for > Transformed Coordinate < ')
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 1',1,x1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 2',1,y1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate 3',1,z1) 
          CALL VariableAdd( Mesh % Variables,Mesh,CurrentModel % Solver,&
              'Transformed Coordinate',3,NewCoords)
        END IF
      END IF
    END IF
      
    IF( UpdateNodes ) THEN
      IF( ListGetLogical( Params,'Coordinate Transformation Use Degrees',Found) ) THEN
        Coeff = 180.0_dp / PI
        CALL Info('CoordinateTranformation','Using degrees for angles')
      ELSE
        Coeff = 1.0_dp
      END IF

      Rad0 = ListGetConstReal( Params,'Coordinate Transformation Radius',Found )
  
      SELECT CASE ( CoordTransform ) 
        
      CASE('cartesian to polar')
        Mode = 1
      CASE('cartesian to cylindrical')
        Mode = 1
      CASE('polar to cartesian')
        Mode = -1
      CASE('cylindrical to cartesian')
        Mode = -1
        
      CASE DEFAULT
        CALL Fatal('CoordinateTransformation','Unknown transformation: '//TRIM(CoordTransform) )
        
      END SELECT

      DO i=1,n    
        R0(1) = x0(i)
        R0(2) = y0(i)
        R0(3) = z0(i)
        
        IF( Mode == 1 ) THEN
          R1(1) = Rad0 + SQRT( R0(1)**2 + R0(2)**2)
          R1(2) = Coeff * ATAN2( R0(2), R0(1)  ) 
          R1(3) = R0(3)    
       
        ELSE IF( Mode == -1 ) THEN
          R1(1) = COS( R0(2) / Coeff ) * ( R0(1) + Rad0 )
          R1(2) = SIN( R0(2) / Coeff ) * ( R0(1) + Rad0 )
          R1(3) = R0(3)          
        END IF

        x1(i) = R1(1)
        y1(i) = R1(2)
        z1(i) = R1(3)

      END DO
    END IF

    IF( .NOT. Irreversible ) THEN
      Mesh % NodesOrig => Mesh % Nodes
      Mesh % Nodes => Mesh % NodesMapped

      Var => VariableGet( CurrentModel % Variables,'Coordinate 1')
      Var % Values => Mesh % Nodes % x

      Var => VariableGet( CurrentModel % Variables,'Coordinate 2')
      Var % Values => Mesh % Nodes % y

      Var => VariableGet( CurrentModel % Variables,'Coordinate 3')
      Var % Values => Mesh % Nodes % z
    END IF

    CALL Info('CoordinateTransformation','All done',Level=8)

  END SUBROUTINE CoordinateTransformation
!---------------------------------------------------------------


!---------------------------------------------------------------
!> Return back to the original coordinate system. 
!---------------------------------------------------------------
  SUBROUTINE BackCoordinateTransformation( Mesh, DeleteTemporalMesh )
    TYPE(Mesh_t) :: Mesh
    LOGICAL, OPTIONAL :: DeleteTemporalMesh
!---------------------------------------------------------------
    TYPE(Variable_t), POINTER :: Var

    IF( PRESENT( DeleteTemporalMesh ) ) THEN
      IF( DeleteTemporalMesh ) THEN
        DEALLOCATE( Mesh % NodesMapped % x, &
            Mesh % NodesMapped % y, &
            Mesh % NodesMapped % z ) 
        DEALLOCATE( Mesh % NodesMapped )
      END IF
    END IF

    IF( .NOT. ASSOCIATED( Mesh % NodesOrig ) ) THEN
      CALL Fatal('BackCoordinateTransformation','NodesOrig not associated')
    END IF

    Mesh % Nodes => Mesh % NodesOrig

    Var => VariableGet( CurrentModel % Variables,'Coordinate 1')
    Var % Values => Mesh % Nodes % x
    
    Var => VariableGet( CurrentModel % Variables,'Coordinate 2')
    Var % Values => Mesh % Nodes % y

    Var => VariableGet( CurrentModel % Variables,'Coordinate 3')
    Var % Values => Mesh % Nodes % z

  END SUBROUTINE BackCoordinateTransformation
!---------------------------------------------------------------


!---------------------------------------------------------------
!> This partitions the mesh into a given number of partitions in each 
!> direction. It may be used in clustering multigrid or similar, 
!> and also to internal partitioning within ElmerSolver. 
!---------------------------------------------------------------
  SUBROUTINE ClusterNodesByDirection(Params,Mesh,Clustering,MaskActive)
 
    USE GeneralUtils

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, POINTER :: Clustering(:)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,GotIt,Hit
    REAL(KIND=dp), ALLOCATABLE :: Measure(:)
    INTEGER :: i,j,k,k0,l,ind,n,dim,dir,divs,nsize,elemsinpart,clusters
    INTEGER, POINTER :: Iarray(:),Order(:),NodePart(:),NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,clustersize
    REAL(KIND=dp), POINTER :: PArray(:,:), Arrange(:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), Weights(3), &
        avepart,devpart
!---------------------------------------------------------------

    ! CALL Info('ClusterNodesByDirection','')

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = COUNT( MaskActive )
    ELSE
      nsize = Mesh % NumberOfNodes
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal('ClusterNodesByDirection','No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    Parray => ListGetConstRealArray( Params,'Clustering Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal(1) = 1.0
      Normal(2) = 1.0d-2
      IF( dim == 3) Normal(3) = 1.0d-4
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )

    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    

    IF( .FALSE. ) THEN
      PRINT *,'Normal:',Normal
      PRINT *,'Tangent1:',Tangent1
      PRINT *,'Tangent2:',Tangent2
    END IF


    Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',GotIt )
    IF(.NOT. GotIt) Iarray => ListGetIntegerArray( Params,'MG Cluster Divisions',GotIt )
    Divisions = 1
    IF( GotIt ) THEN
      n = MIN( SIZE(Iarray), dim ) 
      Divisions(1:n) = Iarray(1:n)
    ELSE
      clustersize = ListGetInteger( Params,'Partitioning Size',GotIt)
      IF(.NOT. GotIt) clustersize = ListGetInteger( Params,'MG Cluster Size',GotIt)
      IF( GotIt .AND. ClusterSize > 0) THEN
        IF( dim == 2 ) THEN
          Divisions(1) = ( nsize / clustersize ) ** 0.5_dp
          Divisions(2) = ( nsize / ( clustersize * Divisions(1) ) )
        ELSE
          Divisions(1:2) = ( nsize / clustersize ) ** (1.0_dp / 3 )
          Divisions(3) = ( nsize / ( clustersize * Divisions(1) * Divisions(2) ) )
        END IF
      ELSE
        CALL Fatal('ClusterNodesByDirection','Clustering Divisions not given!')
      END IF
    END IF

    Clusters = Divisions(1) * Divisions(2) * Divisions(3)

    IF( .FALSE. ) THEN
      PRINT *,'dim:',dim
      PRINT *,'divisions:',divisions
      PRINT *,'clusters:',clusters
      PRINT *,'nsize:',nsize
    END IF

    ALLOCATE(Order(nsize),Arrange(nsize),NodePart(nsize),NoPart(Clusters))
    

    ! These are needed as an initial value for the loop over dimension
    elemsinpart = nsize
    nodepart = 1
    

    ! Go through each direction and cumulatively add to the clusters
    !-----------------------------------------------------------

    DO dir = 1,dim      
      divs = Divisions(dir)
      IF( divs <= 1 ) CYCLE
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF

      ! Initialize ordering for the current direction
      !----------------------------------------------
      DO i=1,nsize
        Order(i) = i
      END DO
      

      ! Now compute the weights for each node
      !----------------------------------------
      DO i=1,Mesh % NumberOfNodes
        j = i
        IF( MaskExists ) THEN
          IF( .NOT. MaskActive(j) ) CYCLE
        END IF
        
        Coord(1) = Mesh % Nodes % x(i)
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)

        Arrange(j) = SUM( Weights * Coord )
      END DO

      ! Order the nodes for given direction
      !----------------------------------------------
      CALL SortR(nsize,Order,Arrange)

      ! For each direction the number of elements in cluster becomes smaller
      elemsinpart = elemsinpart / divs

      ! initialize the counter partition
      nopart = 0


      ! Go through each node and locate it to a cluster taking into consideration
      ! the previous clustering (for 1st direction all one)
      !------------------------------------------------------------------------
      j = 1
      DO i = 1,nsize
        ind = Order(i)
        
        ! the initial partition offset depends on previous partitioning
        k0 = (nodepart(ind)-1) * divs

        ! Find the correct new partitioning, this loop is just long enough
        DO l=1,divs
          Hit = .FALSE.
          
          ! test for increase of local partition
          IF( j < divs ) THEN
            IF( nopart(k0+j) >= elemsinpart ) THEN
              j = j + 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! test for decrease of local partition
          IF( j > 1 )  THEN            
            IF( nopart(k0+j-1) < elemsinpart ) THEN
              j = j - 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! If either increase or decrease is needed, this must be ok 
          IF(.NOT. Hit) EXIT
        END DO
          
        k = k0 + j
        nopart(k) = nopart(k) + 1
        nodepart(ind) = k
      END DO

    END DO


    minpart = HUGE(minpart)
    maxpart = 0
    avepart = 1.0_dp * nsize / clusters
    devpart = 0.0_dp
    DO i=1,clusters
      minpart = MIN( minpart, nopart(i))
      maxpart = MAX( maxpart, nopart(i))
      devpart = devpart + ABS ( nopart(i) - avepart )
    END DO
    devpart = devpart / clusters

    WRITE(Message,'(A,T25,I10)') 'Min nodes in cluster:',minpart
    CALL Info('ClusterNodesByDirection',Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info('ClusterNodesByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info('ClusterNodesByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info('ClusterNodesByDirection',Message)
    

    IF( ASSOCIATED(Clustering)) THEN
      Clustering = Nodepart 
      DEALLOCATE(Nodepart)
    ELSE
      Clustering => Nodepart
      NULLIFY( Nodepart ) 
    END IF
    
    DEALLOCATE(Order,Arrange,NoPart)


  END SUBROUTINE ClusterNodesByDirection



  SUBROUTINE ClusterElementsByDirection(Params,Mesh,Clustering,MaskActive)
 
    USE GeneralUtils

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, POINTER :: Clustering(:)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,GotIt,Hit
    REAL(KIND=dp), ALLOCATABLE :: Measure(:)
    INTEGER :: i,j,k,k0,l,ind,n,dim,dir,divs,nsize,elemsinpart,clusters
    INTEGER, POINTER :: Iarray(:),Order(:),NodePart(:),NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,clustersize
    REAL(KIND=dp), POINTER :: PArray(:,:), Arrange(:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), Weights(3), &
        avepart,devpart, dist
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
!---------------------------------------------------------------

    ! CALL Info('ClusterElementsByDirection','')

    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = COUNT( MaskActive ) 
    ELSE
      nsize = Mesh % NumberOfBulkElements
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal('ClusterElementsByDirection','No parameter list associated')
    END IF

    dim = Mesh % MeshDim
    Parray => ListGetConstRealArray( Params,'Clustering Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal(1) = 1.0
      Normal(2) = 1.0d-2
      IF( dim == 3) THEN
        Normal(3) = 1.0d-4
      ELSE
        Normal(3) = 0.0_dp
      END IF
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )

    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    
    IF( .FALSE. ) THEN
      PRINT *,'Normal:',Normal
      PRINT *,'Tangent1:',Tangent1
      PRINT *,'Tangent2:',Tangent2
    END IF

    Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',GotIt )
    IF(.NOT. GotIt ) THEN
      Iarray => ListGetIntegerArray( Params,'MG Cluster Divisions',GotIt )
    END IF

    Divisions = 1
    IF( GotIt ) THEN
      n = MIN( SIZE(Iarray), dim ) 
      Divisions(1:n) = Iarray(1:n)
    ELSE
      clustersize = ListGetInteger( Params,'Partitioning Size',GotIt)
      IF(.NOT. GotIt) clustersize = ListGetInteger( Params,'MG Cluster Size',GotIt)
      IF( GotIt .AND. ClusterSize > 0) THEN
        IF( dim == 2 ) THEN
          Divisions(1) = ( nsize / clustersize ) ** 0.5_dp
          Divisions(2) = ( nsize / ( clustersize * Divisions(1) ) )
        ELSE
          Divisions(1:2) = ( nsize / clustersize ) ** (1.0_dp / 3 )
          Divisions(3) = ( nsize / ( clustersize * Divisions(1) * Divisions(2) ) )
        END IF
      ELSE
        CALL Fatal('ClusterNodesByDirection','Clustering Divisions not given!')
      END IF
    END IF

    Clusters = Divisions(1) * Divisions(2) * Divisions(3)

    IF( .FALSE. ) THEN
      PRINT *,'dim:',dim
      PRINT *,'divisions:',divisions
      PRINT *,'clusters:',clusters
      PRINT *,'nsize:',nsize
    END IF

    ALLOCATE(Order(nsize),Arrange(nsize),NodePart(nsize),NoPart(Clusters))
    

    ! These are needed as an initial value for the loop over dimension
    elemsinpart = nsize
    nodepart = 1
    

    ! Go through each direction and cumulatively add to the clusters
    !-----------------------------------------------------------

    DO dir = 1,dim      
      divs = Divisions(dir)
      IF( divs <= 1 ) CYCLE
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF

      ! Now compute the weights for each node
      !----------------------------------------
      j = 0
      DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        IF( MaskExists ) THEN
          IF( .NOT. MaskActive( i ) ) CYCLE
        ELSE 
          IF( i > Mesh % NumberOfBulkElements ) EXIT
        END IF
        
        Element => Mesh % Elements(i)
        NodeIndexes => Element % NodeIndexes 
        n = Element % TYPE % NumberOfNodes

        Coord(1) = SUM( Mesh % Nodes % x( NodeIndexes ) ) / n
        Coord(2) = SUM( Mesh % Nodes % y( NodeIndexes ) ) / n
        Coord(3) = SUM( Mesh % Nodes % z( NodeIndexes ) ) / n

        j = j + 1
        Arrange(j) = SUM( Weights * Coord )

        ! Initialize ordering for the current direction
        Order(j) = j
      END DO

      ! Order the distances for given direction, only the active ones
      !--------------------------------------------------------------
      CALL SortR(nsize,Order,Arrange)

      ! For each direction the number of elements in cluster becomes smaller
      elemsinpart = elemsinpart / divs

      ! initialize the counter partition
      nopart = 0

      ! Go through each node and locate it to a cluster taking into consideration
      ! the previous clustering (for 1st direction all one)
      !------------------------------------------------------------------------
      j = 1
      DO i = 1,nsize
        ind = Order(i)
        
        ! the initial partition offset depends on previous partitioning
        k0 = (nodepart(ind)-1) * divs

        ! Find the correct new partitioning, this loop is just long enough
        DO l=1,divs
          Hit = .FALSE.
          
          ! test for increase of local partition
          IF( j < divs ) THEN
            IF( nopart(k0+j) >= elemsinpart ) THEN
              j = j + 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! test for decrease of local partition
          IF( j > 1 )  THEN            
            IF( nopart(k0+j-1) < elemsinpart ) THEN
              j = j - 1
              Hit = .TRUE.
            END IF
          END IF
          
          ! If either increase or decrease is needed, this must be ok 
          IF(.NOT. Hit) EXIT
        END DO
          
        k = k0 + j
        nopart(k) = nopart(k) + 1

        ! Now set the partition 
        nodepart(ind) = k
      END DO

    END DO


    minpart = HUGE(minpart)
    maxpart = 0
    avepart = 1.0_dp * nsize / clusters
    devpart = 0.0_dp
    DO i=1,clusters
      minpart = MIN( minpart, nopart(i))
      maxpart = MAX( maxpart, nopart(i))
      devpart = devpart + ABS ( nopart(i) - avepart )
    END DO
    devpart = devpart / clusters

    WRITE(Message,'(A,T25,I10)') 'Min nodes in cluster:',minpart
    CALL Info('ClusterElementsByDirection',Message)
    WRITE(Message,'(A,T25,I10)') 'Max nodes in cluster:',maxpart
    CALL Info('ClusterElementsByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Average nodes in cluster:',avepart
    CALL Info('ClusterElementsByDirection',Message)
    WRITE(Message,'(A,T28,F10.2)') 'Deviation of nodes:',devpart
    CALL Info('ClusterElementsByDirection',Message)
    
    
    IF( ASSOCIATED(Clustering)) THEN
      IF( PRESENT( MaskActive ) ) THEN
        j = 0
        DO i=1, SIZE(MaskActive)
          IF( MaskActive(i) ) THEN
            j = j + 1
            Clustering(i) = Nodepart(j)
          END IF
        END DO
      ELSE
        Clustering = Nodepart 
      END IF
      DEALLOCATE(Nodepart)
    ELSE
      Clustering => Nodepart
      NULLIFY( Nodepart ) 
    END IF
    
    DEALLOCATE(Order,Arrange,NoPart)

  END SUBROUTINE ClusterElementsByDirection



  SUBROUTINE ClusterElementsUniform(Params,Mesh,Clustering,MaskActive,PartitionDivisions)
 
    USE GeneralUtils

    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: Clustering(:)
    LOGICAL, OPTIONAL :: MaskActive(:)
    INTEGER, OPTIONAL :: PartitionDivisions(3)
!---------------------------------------------------------------
    LOGICAL :: MaskExists,UseMaskedBoundingBox,Found
    INTEGER :: i,j,k,ind,n,dim,nsize,nmask,clusters
    INTEGER, POINTER :: Iarray(:),ElemPart(:)
    INTEGER, ALLOCATABLE :: NoPart(:)
    INTEGER :: Divisions(3),minpart,maxpart,Inds(3)
    REAL(KIND=dp) :: Coord(3), Weights(3), avepart,devpart
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp) :: BoundingBox(6)
    INTEGER, ALLOCATABLE :: CellCount(:,:,:)
    LOGICAL, ALLOCATABLE :: NodeMask(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: Caller="ClusterElementsUniform"

    CALL Info(Caller,'Clustering elements uniformly in bounding box',Level=6)

    IF( Mesh % NumberOfBulkElements == 0 ) RETURN
    
    MaskExists = PRESENT(MaskActive)
    IF( MaskExists ) THEN
      nsize = SIZE( MaskActive ) 
      nmask = COUNT( MaskActive ) 
      CALL Info(Caller,'Applying division to masked element: '//TRIM(I2S(nmask)),Level=8)
    ELSE
      nsize = Mesh % NumberOfBulkElements 
      nmask = nsize
      CALL Info(Caller,'Applying division to all bulk elements: '//TRIM(I2S(nsize)),Level=8)
    END IF
     
    IF( .NOT. ASSOCIATED( Params ) ) THEN
      CALL Fatal(Caller,'No parameter list associated')
    END IF

    dim = Mesh % MeshDim

    ! We can use the masked bounding box
    UseMaskedBoundingBox = .FALSE.
    IF( MaskExists ) UseMaskedBoundingBox = ListGetLogical( Params,&
        'Partition Masked Bounding Box',Found ) 

    IF( UseMaskedBoundingBox ) THEN
      ALLOCATE( NodeMask( Mesh % NumberOfNodes ) )
      NodeMask = .FALSE.

      ! Add all active nodes to the mask
      DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
        IF( .NOT. MaskActive( i ) ) CYCLE        
        Element => Mesh % Elements(i)
        NodeIndexes => Element % NodeIndexes 
        NodeMask( NodeIndexes ) = .TRUE.
      END DO

      i = COUNT( NodeMask ) 
      CALL Info(Caller,'Masked elements include nodes: '//TRIM(I2S(i)),Level=8)
      
      ! Define the masked bounding box
      BoundingBox(1) = MINVAL( Mesh % Nodes % x, NodeMask )
      BoundingBox(2) = MAXVAL( Mesh % Nodes % x, NodeMask )
      BoundingBox(3) = MINVAL( Mesh % Nodes % y, NodeMask )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % y, NodeMask )
      BoundingBox(5) = MINVAL( Mesh % Nodes % z, NodeMask )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z, NodeMask )

      DEALLOCATE( NodeMask ) 
    ELSE      
      BoundingBox(1) = MINVAL( Mesh % Nodes % x )
      BoundingBox(2) = MAXVAL( Mesh % Nodes % x )
      BoundingBox(3) = MINVAL( Mesh % Nodes % y )
      BoundingBox(4) = MAXVAL( Mesh % Nodes % y )
      BoundingBox(5) = MINVAL( Mesh % Nodes % z )
      BoundingBox(6) = MAXVAL( Mesh % Nodes % z )
    END IF
      
    
    IF( PRESENT( PartitionDivisions ) ) THEN
      Divisions = PartitionDivisions
    ELSE      
      Iarray => ListGetIntegerArray( Params,'Partitioning Divisions',Found)
      IF(.NOT. Found ) THEN
        CALL Fatal(Caller,'> Partitioning Divisions < not given!')
      END IF      
      Divisions = 1
      IF( Found ) THEN
        n = MIN( SIZE(Iarray), dim ) 
        Divisions(1:n) = Iarray(1:n)
      END IF
    END IF
      
    ALLOCATE( CellCount(Divisions(1), Divisions(2), Divisions(3) ) )
    CellCount = 0
    Clusters = 1
    DO i=1,dim
      Clusters = Clusters * Divisions(i)
    END DO

    IF( .FALSE. ) THEN
      PRINT *,'dim:',dim
      PRINT *,'divisions:',divisions
      PRINT *,'clusters:',clusters
      PRINT *,'nsize:',nsize
    END IF

    ALLOCATE(ElemPart(nsize),NoPart(Clusters))
    NoPart = 0
    ElemPart = 0

    !----------------------------------------
    Inds = 1
    Coord = 0.0_dp

    DO i=1,Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
      IF( MaskExists ) THEN
        IF( .NOT. MaskActive( i ) ) CYCLE
      ELSE 
        IF( i > Mesh % NumberOfBulkElements ) EXIT
      END IF
      
      Element => Mesh % Elements(i)
      NodeIndexes => Element % NodeIndexes 
      n = Element % TYPE % NumberOfNodes
      
      ! Find the center of the element
      Coord(1) = SUM( Mesh % Nodes % x( NodeIndexes ) ) / n
      Coord(2) = SUM( Mesh % Nodes % y( NodeIndexes ) ) / n
      IF( dim == 3 ) THEN
        Coord(3) = SUM( Mesh % Nodes % z( NodeIndexes ) ) / n
      END IF

      Inds = 1
      DO j=1,dim
        Inds(j) = CEILING( Divisions(j) * &
            ( Coord(j) - BoundingBox(2*j-1) ) / &
            ( BoundingBox(2*j) - BoundingBox(2*j-1) ) )
      END DO
      Inds = MAX( Inds, 1 ) 

      CellCount(Inds(1),Inds(2),Inds(3)) = &
          CellCount(Inds(1),Inds(2),Inds(3)) + 1

      ind = (Inds(1)-1)*Divisions(2)*Divisions(3) + &
          (Inds(2)-1)*Divisions(3) +  &
          Inds(3)
      ElemPart(i) = ind
      NoPart(ind) = NoPart(ind) + 1
    END DO

    ! Compute statistical information of the partitioning
    n = COUNT( NoPart > 0 )    
    minpart = HUGE(minpart)
    maxpart = 0
    avepart = 1.0_dp * nmask / n
    devpart = 0.0_dp
    DO i=1,clusters
      IF( nopart(i) > 0 ) THEN
        minpart = MIN( minpart, nopart(i))
        maxpart = MAX( maxpart, nopart(i))
        devpart = devpart + ABS ( nopart(i) - avepart )
      END IF
    END DO
    devpart = devpart / n

    CALL Info(Caller,'Number of partitions: '//TRIM(I2S(n)),Level=8)
    CALL Info(Caller,'Min elements in cluster: '//TRIM(I2S(minpart)),Level=8)
    CALL Info(Caller,'Max elements in cluster: '//TRIM(I2S(maxpart)),Level=8)

    WRITE(Message,'(A,F10.2)') 'Average elements in cluster:',avepart
    CALL Info(Caller,Message,Level=8)    
    WRITE(Message,'(A,F10.2)') 'Average deviation in size:',devpart
    CALL Info(Caller,Message,Level=8)

    ! Renumber the partitions using only the active ones
    n = 0
    DO i=1,clusters
      IF( NoPart(i) > 0 ) THEN
        n = n + 1
        NoPart(i) = n
      END IF
    END DO
    
    ! Renumbering only needed if there are empty cells
    IF( n < clusters ) THEN
      DO i=1,nsize
        j = ElemPart(i)
        IF( j > 0 ) ElemPart(i) = NoPart(j) 
      END DO
    END IF

    !DO i=1,clusters
    !  PRINT *,'count in part:',i,COUNT( ElemPart(1:nsize) == i ) 
    !END DO
    
    IF( ASSOCIATED( Clustering ) ) THEN
      WHERE( ElemPart > 0 ) Clustering = ElemPart
      DEALLOCATE( ElemPart ) 
    ELSE
      Clustering => ElemPart
      NULLIFY( ElemPart ) 
    END IF
    
    DEALLOCATE(NoPart,CellCount)

    CALL Info(Caller,'Clustering of elements finished',Level=10)

  END SUBROUTINE ClusterElementsUniform

 
  !> Find the node closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-----------------------------------------------------------------
  FUNCTION ClosestNodeInMesh(Mesh,Coord,MinDist) RESULT ( NodeIndx )
    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coord(3)
    REAL(KIND=dp), OPTIONAL :: MinDist
    INTEGER :: NodeIndx

    REAL(KIND=dp) :: Dist2,MinDist2,NodeCoord(3)
    INTEGER :: i

    MinDist2 = HUGE( MinDist2 ) 

    DO i=1,Mesh % NumberOfNodes
      
      NodeCoord(1) = Mesh % Nodes % x(i)
      NodeCoord(2) = Mesh % Nodes % y(i)
      NodeCoord(3) = Mesh % Nodes % z(i)
    
      Dist2 = SUM( ( Coord - NodeCoord )**2 )
      IF( Dist2 < MinDist2 ) THEN
        MinDist2 = Dist2
        NodeIndx = i  
      END IF
    END DO
    
    IF( PRESENT( MinDist ) ) MinDist = SQRT( MinDist2 ) 

  END FUNCTION ClosestNodeInMesh


  !> Find the element that owns or is closest to the given coordinate. 
  !> The linear search only makes sense for a small number of points. 
  !> Users include saving routines of pointwise information. 
  !-------------------------------------------------------------------
  FUNCTION ClosestElementInMesh(Mesh, Coords) RESULT ( ElemIndx )

    TYPE(Mesh_t) :: Mesh
    REAL(KIND=dp) :: Coords(3)
    INTEGER :: ElemIndx

    REAL(KIND=dp) :: Dist,MinDist,LocalCoords(3)
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Nodes_t) :: ElementNodes
    INTEGER :: k,l,n,istat
    REAL(KIND=dp) :: ParallelHits,ParallelCands
    LOGICAL :: Hit

    n = Mesh % MaxElementNodes
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), STAT=istat)
    IF( istat /= 0 ) CALL Fatal('ClosestElementInMesh','Memory allocation error') 	
    ElemIndx = 0
    MinDist = HUGE( MinDist ) 
    Hit = .FALSE.
    l = 0
    
    ! Go through all bulk elements and look for hit in each element.
    ! Linear search makes only sense for a small number of nodes
    DO k=1,Mesh % NumberOfBulkElements

      Element => Mesh % Elements(k)
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
      
      Hit = PointInElement( Element, ElementNodes, &
          Coords, LocalCoords, LocalDistance = Dist )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        l = k
      END IF
      IF( Hit ) EXIT
    END DO
    
    ! Count the number of parallel hits
    !-----------------------------------------------------------------------
    IF( Hit ) THEN
      ParallelHits = 1.0_dp
    ELSE
      ParallelHits = 0.0_dp
    END IF
    ParallelHits = ParallelReduction( ParallelHits )
    
    ! If there was no proper hit go through the best candidates so far and 
    ! see if they would give a acceptable hit
    !----------------------------------------------------------------------
    IF( ParallelHits < 0.5_dp ) THEN	  

      ! Compute the number of parallel candidates
      !------------------------------------------
      IF( l > 0 ) THEN
        ParallelCands = 1.0_dp
      ELSE
        ParallelCands = 0.0_dp
      END IF
      ParallelCands = ParallelReduction( ParallelCands ) 

      IF( l > 0 ) THEN
        Element => Mesh % Elements(l)
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

        ! If there are more than two competing parallel hits then use more stringent conditions
        ! since afterwords there is no way of deciding which one was closer.
        !--------------------------------------------------------------------------------------
        IF( ParallelCands > 1.5_dp ) THEN
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0d-3, LocalEps=1.0d-4 )	
        ELSE
          Hit = PointInElement( Element, ElementNodes, &
              Coords, LocalCoords, GlobalEps = 1.0_dp, LocalEps=0.1_dp )	
        END IF
      END IF
    END IF

    IF( Hit ) ElemIndx = l

    IF( ParallelHits < 0.5_dp ) THEN
      IF( Hit ) THEN
        ParallelHits = 1.0_dp
      ELSE
        ParallelHits = 0.0_dp
      END IF
      ParallelHits = ParallelReduction( ParallelHits )
      IF( ParallelHits < 0.5_dp ) THEN
        WRITE( Message, * ) 'Coordinate not found in any of the elements!',Coords
        CALL Warn( 'ClosestElementInMesh', Message )
      END IF
    END IF

    DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
 
  END FUNCTION ClosestElementInMesh



!---------------------------------------------------------------
!> This find two fixing nodes for each coordinate direction
!> The indexes are returned in order: x1 x2 y1 y2 z1 z2.
!---------------------------------------------------------------
  SUBROUTINE FindRigidBodyFixingNodes(Solver,FixingDofs,MaskPerm)
!------------------------------------------------------------------------------
    USE GeneralUtils

    TYPE(Solver_t) :: Solver
    INTEGER, OPTIONAL :: FixingDofs(0:)
    INTEGER, OPTIONAL :: MaskPerm(:)

!---------------------------------------------------------------

    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: MaskExists,FixBestDirection,FoundBetter, GotIt
    INTEGER :: i,j,k,l,ind,n,dim,dir,nsize,Sweep,MaxSweep,DirBest
    INTEGER :: PosMeasureIndex, NegMeasureIndex, FixingNodes(0:6)
    LOGICAL, ALLOCATABLE :: ForbiddenNodes(:)
    REAL(KIND=dp), POINTER :: Parray(:,:)
    REAL(KIND=dp) :: Normal(3), Tangent1(3), Tangent2(3), Coord(3), &
        SumCoord(3), AveCoord(3), Weights(3), RefScore, Score, &
        PosMeasure, NegMeasure, OffLineCoeff, DirDistance, &
        InLine, OffLine, Dist, MinDist, InLineMeasure, ScoreLimit
    CHARACTER(LEN=MAX_NAME_LEN) :: Method
!---------------------------------------------------------------

    CALL Info('FindRigidBodyFixingNodes','Starting',Level=6)

    Mesh => Solver % Mesh
    dim = Mesh % MeshDim 
    
    ALLOCATE( ForbiddenNodes(Mesh % NumberOfNodes) )
    CALL DetermineForbiddenNodes( )
    nsize = COUNT(.NOT. ForbiddenNodes) 

!   PRINT *,'Number of allowed Nodes:',nsize

    ! Find the center from the average of node positions
    !-----------------------------------------------------------
    SumCoord = 0.0_dp
    DO i=1,Mesh % NumberOfNodes
      IF( ForbiddenNodes( i ) ) CYCLE
      
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
    
      SumCoord = SumCoord + Coord
    END DO
    AveCoord = SumCoord / nsize


    ! Find the node closest to center and make that the new center
    !--------------------------------------------------------------
    MinDist = HUGE( MinDist ) 

    DO i=1,Mesh % NumberOfNodes
      IF( ForbiddenNodes( i ) ) CYCLE
      
      Coord(1) = Mesh % Nodes % x(i)
      Coord(2) = Mesh % Nodes % y(i)
      Coord(3) = Mesh % Nodes % z(i)
    
      Dist = SUM( ( Coord - AveCoord )**2 )
      IF( Dist < MinDist ) THEN
        MinDist = Dist
        k = i  
      END IF
    END DO

    AveCoord(1) = Mesh % Nodes % x(k)
    AveCoord(2) = Mesh % Nodes % y(k)
    AveCoord(3) = Mesh % Nodes % z(k)
    IF(PRESENT(FixingDOFs)) FixingDOFs(0)=k
    

!   PRINT *,'AveCoord:',AveCoord

    ! Parameters of the search
    !-----------------------------------------------------------

    OffLineCoeff = ListGetConstReal( Solver % Values,'Fixing Nodes Off Line Coefficient',GotIt)
    IF(.NOT. GotIt) OffLineCoeff = 1.0_dp

    ScoreLimit = ListGetConstReal( Solver % Values,'Fixing Nodes Limit Score',GotIt)
    IF(.NOT. GotIt) ScoreLimit = 0.99_dp

    FixBestDirection = ListGetLogical( Solver % Values,'Fixing Nodes Axis Freeze',GotIt)

    Parray => ListGetConstRealArray( Solver % Values,'Fixing Nodes Normal Vector',GotIt )
    IF( GotIt ) THEN
      Normal = Parray(1:3,1)
    ELSE
      Normal = 0.0_dp
      Normal(1) = 1.0
    END IF
    Normal = Normal / SQRT( SUM( Normal ** 2) )      
    CALL TangentDirections( Normal,Tangent1,Tangent2 )
    
    ! Find the fixing nodes by looping over all nodes
    !-----------------------------------------------------------
    DirDistance = 0.0_dp
    DirBest = 0
    DO dir = 1, dim
      
      ! Use the three principal directions as the weight
      !-------------------------------------------------
      IF( dir == 1 ) THEN
        Weights = Normal
      ELSE IF( dir == 2 ) THEN
        Weights = Tangent1
      ELSE 
        Weights = Tangent2
      END IF
      
      PosMeasure = 0.0_dp
      PosMeasureIndex = 0
      NegMeasure = 0.0_dp
      NegMeasureIndex = 0


      ! Choose the nodes within the cones in the given three directions
      !---------------------------------------------------------------
      DO i=1,Mesh % NumberOfNodes
        IF( ForbiddenNodes( i ) ) CYCLE
        
        Coord(1) = Mesh % Nodes % x(i) 
        Coord(2) = Mesh % Nodes % y(i)
        Coord(3) = Mesh % Nodes % z(i)
        
        Coord = Coord - AveCoord
        Dist = SQRT( SUM( Coord ** 2 ) )
 
        ! Signed distance in in-line direction
        InLine = SUM( Coord * Weights )
        
        ! Distance in off-line direction 
        OffLine = SQRT( Dist**2 - InLine**2 )
        
        ! This defines a cone within which nodes are accepted
        InLineMeasure = ABS( InLine ) - OffLineCoeff * OffLine 
        IF( InLineMeasure < 0.0_dp ) CYCLE
        
        IF( InLine < 0.0_dp ) THEN
          IF( InLineMeasure > NegMeasure ) THEN
            NegMeasure = InLineMeasure
            NegMeasureIndex = i
          END IF
        ELSE           
          IF( InLineMeasure > PosMeasure ) THEN
            PosMeasure = InLineMeasure 
            PosMeasureIndex = i
          END IF
        END IF      
      END DO
      
      FixingNodes(2*dir-1) = NegMeasureIndex
      FixingNodes(2*dir) = PosMeasureIndex      

      IF( NegMeasureIndex > 0 .AND. PosMeasureIndex > 0 ) THEN
        IF( PosMeasure + NegMeasure > DirDistance ) THEN
          DirDistance = PosMeasure + NegMeasure
          DirBest = dir
        END IF
      END IF

    END DO


 
    ! To be on the safe side check that no node is used twice
    ! However, do not break the best direction
    !-----------------------------------------------------------------------------------
    DO i=1,2*dim
      DO j=1,2*dim
        IF( FixBestDirection ) THEN
          IF( j == 2*DirBest-1 .OR. j == 2*DirBest ) CYCLE
        END IF        
        IF( FixingNodes(j) == FixingNodes(i) ) FixingNodes(j) = 0
      END DO
    END DO


    ! Go through the fixing nodes one-by-one and set the node so that the harmonic sum
    ! is minimized. This means that small distances are hopefully eliminated. 
    !-----------------------------------------------------------------------------------
    MaxSweep = ListGetInteger( Solver % Values,'Fixing Nodes Search Loops',GotIt)
    DO Sweep = 0,MaxSweep
      FoundBetter = .FALSE.
      DO j=1,2*dim 
        RefScore = FixingNodesScore(j,FixingNodes(j)) 

        ! The first round set the unfixed nodes
        IF( Sweep == 0 ) THEN
!         PRINT *,'Initial Score:',j,RefScore
          IF( FixingNodes(j) /= 0 ) CYCLE
        END IF

        ! Fir the best direction because otherwise there are too 
        ! many moving parts.
        IF( FixBestDirection ) THEN
          IF( j == 2*DirBest-1 .OR. j == 2*DirBest ) CYCLE
        END IF

        RefScore = FixingNodesScore(j,FixingNodes(j)) 

        DO i=1,Mesh % NumberOfNodes
          IF( ForbiddenNodes(i) ) CYCLE
          Score = FixingNodesScore(j,i)
          IF( Score < ScoreLimit * RefScore ) THEN
            RefScore = Score 
            FixingNodes(j) = i            
            FoundBetter = .TRUE.
          END IF
        END DO
      END DO
      IF(.NOT. FoundBetter ) EXIT
    END DO

    DO j=1,2*dim
      RefScore = FixingNodesScore(j,FixingNodes(j)) 
!     PRINT *,'Final Score:',j,RefScore
    END DO

    ! Output the selected nodes
    !-----------------------------------------------------------------------------------
    DO i=1,2*dim
      j = FixingNodes(i)
      WRITE(Message,'(A,I0,3ES10.2)') 'Fixing Node: ',j,&
          Mesh % Nodes % x( j ), &
          Mesh % Nodes % y( j ), &
          Mesh % Nodes % z( j ) 
      CALL Info('FindRigidBodyFixingNodes',Message,Level=6)
      IF( PRESENT( FixingDofs ) ) FixingDofs(i) = j     
    END DO

    DEALLOCATE( ForbiddenNodes )


  CONTAINS

    !> Find the nodes that are either on interface, boundary or do not belong to the field.
    !-----------------------------------------------------------------------------------
    SUBROUTINE DetermineForbiddenNodes()

      TYPE(Element_t), POINTER :: Element
      LOGICAL, POINTER :: ig(:)
      INTEGER :: t
      
      ! Mark all interface nodes as forbidden nodes
      !-----------------------------------------------
      IF( ParEnv % PEs > 1 ) THEN
        ig => Mesh % ParallelInfo % INTERFACE
        ForbiddenNodes = ig(1:Mesh % NumberOfNodes)
      END IF

      ! Mark all nodes on boundary elements as forbidden nodes
      !--------------------------------------------------------
      DO t=Mesh % NumberOfBulkElements + 1, &
          Mesh % NumberOfBoundaryElements + Mesh % NumberOfBulkElements

        Element => Mesh % Elements( t )
        ForbiddenNodes( Element % NodeIndexes ) = .TRUE.
      END DO

      ! If mask exists then add all nodes not in mask to forbidden nodes
      !-----------------------------------------------------------------
      IF( PRESENT( MaskPerm) ) THEN
        DO i=1,Mesh % NumberOfNodes
          IF( MaskPerm(i) == 0 ) ForbiddenNodes(i) = .TRUE.
        END DO
      END IF
      
    END SUBROUTINE DetermineForbiddenNodes


    !> Give a value of goodness to the chosen fixing node.
    !-----------------------------------------------------------------------------------
    FUNCTION FixingNodesScore(direction,cand) RESULT ( Score )

      INTEGER :: direction, cand
      INTEGER :: i,j
      REAL(KIND=dp) :: Score

      REAL(KIND=dp) :: x0(3), x1(3), Dist

      IF( cand == 0 ) THEN
        Score = HUGE( Score ) 
        RETURN
      END IF

      Score = 0.0_dp
      x0(1) = Mesh % Nodes % x( cand )
      x0(2) = Mesh % Nodes % y( cand )
      x0(3) = Mesh % Nodes % z( cand )

      DO i=1,2*dim
        IF( i == direction ) CYCLE
        j = FixingNodes( i )

        ! Do not measure distance to unset nodes!
        IF( j == 0 ) CYCLE

        ! This would lead to division by zero later on
        IF( cand == j ) THEN
          Score = HUGE( Score ) 
          RETURN
        END IF

        x1(1) = Mesh % Nodes % x( j )
        x1(2) = Mesh % Nodes % y( j )
        x1(3) = Mesh % Nodes % z( j )

        Dist = SQRT( SUM( (x0 - x1 ) ** 2 ) )
        Score = Score + 1 / Dist
      END DO

    END FUNCTION FixingNodesScore


!------------------------------------------------------------------------------
  END SUBROUTINE FindRigidBodyFixingNodes
!------------------------------------------------------------------------------



  SUBROUTINE ElmerMeshToDualGraph(Mesh, DualGraph, UseBoundaryMesh)
    IMPLICIT NONE

    TYPE(Mesh_t) :: Mesh
    TYPE(Graph_t) :: DualGraph
    LOGICAL, OPTIONAL :: UseBoundaryMesh

    TYPE(Element_t), POINTER :: Element, Elements(:)

    ! MESH DATA
    ! Mesh (CRS format)
    INTEGER, ALLOCATABLE :: eptr(:), eind(:)
    INTEGER :: nelem
    ! Vertex to element map (CRS format)
    INTEGER, ALLOCATABLE :: vptr(:), vind(:)
    INTEGER :: nvertex

    ! WORK ARRAYS
    ! Pointers to vertex-element maps of the current element
    INTEGER, ALLOCATABLE :: ptrli(:), ptrti(:)
    ! Neighbour indices
    INTEGER, ALLOCATABLE :: neighind(:)
    ! ARRAY MERGE: map for merge
    INTEGER, ALLOCATABLE :: wrkmap(:)

    TYPE :: IntTuple_t
      INTEGER :: i1, i2
    END type IntTuple_t

    TYPE(IntTuple_t), ALLOCATABLE :: wrkheap(:)

    ! OpenMP thread block leads for work division
    INTEGER, ALLOCATABLE :: thrblk(:)
    ! Work indices
    INTEGER, ALLOCATABLE :: wrkind(:), wrkindresize(:)
    INTEGER :: nwrkind

    ! Variables
    INTEGER :: i, dnnz, eid, nl, nli, nti, nn, nv, nthr, &
            te, thrli, thrti, vli, vti, TID, allocstat
    INTEGER :: mapSizePad, maxNodesPad, neighSizePad
    LOGICAL :: Boundary

    INTEGER, PARAMETER :: HEAPALG_THRESHOLD = 24

    CALL Info('ElmerMeshToDualGraph','Creating a dual graph for the mesh',Level=8)

    Boundary = .FALSE.
    IF (Present(UseBoundaryMesh)) Boundary = UseBoundaryMesh

    ! Pointers to mesh data
    IF (.NOT. Boundary) THEN
       nelem = Mesh % NumberOfBulkElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements
    ELSE
       nelem = Mesh % NumberOfBoundaryElements
       nvertex = Mesh % NumberOfNodes
       Elements => Mesh % Elements(&
            Mesh % NumberOfBulkElements+1:Mesh % NumberOfBulkElements+nelem)
    END IF

    ! Initialize dual mesh size and number of nonzeroes
    DualGraph % n = nelem
    dnnz = 0

    ! Copy mesh to CRS structure
    ALLOCATE(eptr(nelem+1), eind(nelem*Mesh % MaxElementNodes), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate mesh structure!')

    eptr(1)=1 ! Fortran numbering
    DO i=1, nelem
      Element => Elements(i)
      nl = Element % TYPE % NumberOfNodes
      nli = eptr(i) ! Fortran numbering
      nti = nli+nl-1
      eind(nli:nti) = Element % NodeIndexes(1:nl) ! Fortran numbering
      eptr(i+1) = nli+nl
    END DO

    ! Construct vertex to element list (in serial!)
    CALL VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)

    ! Allocate pointers to dual mesh
    ALLOCATE(DualGraph % ptr(nelem+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')

    ! Divide work by number of rows in the vertex graph
    nthr = 1 
    !$ nthr = omp_get_max_threads()

    ! Load balance the actual work done by threads (slow)
    ! CALL ThreadLoadBalanceElementNeighbour(nthr, nelem, eptr, eind, vptr, thrblk)
    CALL ThreadStaticWorkShare(nthr, nelem, thrblk)

    !$OMP PARALLEL SHARED(nelem, nvertex, eptr, eind, &
    !$OMP                 vptr, vind, Mesh, DualGraph, &
    !$OMP                 nthr, thrblk, dnnz) &
    !$OMP PRIVATE(i, eid, nli, nti, nn, nv, vli, vti, te, &
    !$OMP         maxNodesPad, neighSizePad, ptrli, ptrti, &
    !$OMP         wrkheap, wrkmap, neighind, &
    !$OMP         wrkind, nwrkind, wrkindresize, allocstat, &
    !$OMP         mapSizePad, thrli, thrti, TID) NUM_THREADS(nthr) &
    !$OMP DEFAULT(NONE)

    TID = 1
    !$ TID = OMP_GET_THREAD_NUM()+1

    ! Ensure that the vertex to element lists are sorted
    !$OMP DO 
    DO i=1,nvertex
      vli = vptr(i)
      vti = vptr(i+1)-1

      CALL Sort(vti-vli+1, vind(vli:vti))
    END DO
    !$OMP END DO NOWAIT

    ! Allocate work array (local to each thread)
    maxNodesPad = IntegerNBytePad(Mesh % MaxElementNodes, 8)
    neighSizePad = IntegerNBytePad(Mesh % MaxElementNodes*20, 8)

    ! Pointers to vertex maps
    ALLOCATE(neighind(neighSizePad), &
            ptrli(maxNodesPad), ptrti(maxNodesPad), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')
    ! Initialize neighbour indices
    neighind = 0

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      ! With multiple threads, use heap based merge
      ALLOCATE(wrkheap(maxNodesPad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
    ELSE
      ! With a small number of threads, use map -based merge
      mapSizePad = IntegerNBytePad(nelem, 8)
      ALLOCATE(wrkmap(mapSizePad), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
              'Unable to allocate local workspace!')
      ! Initialize local map
      wrkmap=0
    END IF

    ! Allocate local list for results
    nwrkind = 0
    ALLOCATE(wrkind(nelem/nthr*20), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate local workspace!')

    ! Ensure that all the threads have finished sorting the vertex indices
    !$OMP BARRIER

    ! Get thread indices
    thrli = thrblk(TID)
    thrti = thrblk(TID+1)

    ! For each element
    DO eid=thrli,thrti-1
      nli = eptr(eid)
      nti = eptr(eid+1)-1
      nv = nti-nli+1

      ! Get pointers to vertices related to the nodes of the element
      te = 0
      DO i=nli,nti
        ptrli(i-nli+1)=vptr(eind(i))
        ptrti(i-nli+1)=vptr(eind(i)+1) ! NOTE: This is to make comparison cheaper
        te = te + ptrti(i-nli+1)-ptrli(i-nli+1)
      END DO

      ! Allocate neighind large enough
      IF (SIZE(neighind)<te) THEN
        DEALLOCATE(neighind)
        neighSizePad = IntegerNBytePad(te,8)
        ALLOCATE(neighind(neighSizePad), STAT=allocstat)
        neighind = 0
      END IF

      ! Merge vertex lists (multi-way merge of ordered lists)
      IF (nthr >= HEAPALG_THRESHOLD) THEN
        CALL kWayMergeHeap(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkheap)
      ELSE
        CALL kWayMergeArray(eid, nv, ptrli, ptrti, &
                te, vind, nn, neighind, wrkmap)
      END IF

      ! Add merged list to final list of vertices
      IF (nn+nwrkind>SIZE(wrkind)) THEN
        ALLOCATE(wrkindresize(MAX(nn+nwrkind,2*SIZE(wrkind))), STAT=allocstat)
        IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
                'Unable to allocate local workspace!')
        wrkindresize(1:nwrkind)=wrkind(1:nwrkind)
        DEALLOCATE(wrkind)
        CALL MOVE_ALLOC(wrkindresize, wrkind)
      END IF
      wrkind(nwrkind+1:nwrkind+nn) = neighind(1:nn)
      nwrkind = nwrkind + nn

      ! Store number of row nonzeroes
      DualGraph % ptr(eid)=nn
    END DO

    ! Get the global size of the dual mesh
    !$OMP DO REDUCTION(+:dnnz)
    DO i=1,nthr
      dnnz = nwrkind
    END DO
    !$OMP END DO

    ! Allocate memory for dual mesh indices
    !$OMP SINGLE
    ALLOCATE(DualGraph % ind(dnnz), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to allocate dual mesh!')
    ! ptr stores row counts, build crs pointers from them
    CALL ComputeCRSIndexes(nelem, DualGraph % ptr)
    !$OMP END SINGLE

    DualGraph % ind(&
            DualGraph % ptr(thrli):DualGraph % ptr(thrti)-1)=wrkind(1:nwrkind)

    IF (nthr >= HEAPALG_THRESHOLD) THEN
      DEALLOCATE(wrkheap, STAT=allocstat)
    ELSE
      DEALLOCATE(wrkmap, STAT=allocstat)
    END IF
    IF (allocstat /= 0) CALL Fatal('ElmerMeshToDualGraph', &
            'Unable to deallocate local workspace!')
    DEALLOCATE(neighind, ptrli, ptrti, wrkind)

    !$OMP END PARALLEL

    ! Deallocate the rest of memory
    DEALLOCATE(eind, eptr, vptr, vind, thrblk)

    CALL Info('ElmerMeshToDualGraph','Dual graph created with size '//TRIM(I2S(dnnz)),Level=8)


  CONTAINS

    SUBROUTINE VertexToElementList(nelem, nvertex, eptr, eind, vptr, vind)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nelem, nvertex
      INTEGER :: eptr(:), eind(:)
      INTEGER, ALLOCATABLE :: vptr(:), vind(:)

      INTEGER :: i, j, v, eli, eti, ind, tmpi, tmpip, allocstat

      ! Initialize vertex structure (enough storage for nvertex vertices
      ! having eptr(nelem+1) elements)
      ALLOCATE(vptr(nvertex+1), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
              'Vertex allocation failed!')
      vptr = 0

      ! For each element

      ! Compute number of elements attached to each vertex (size of lists)
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        DO j=eli, eti
          vptr(eind(j))=vptr(eind(j))+1
        END DO
      END DO

      ! Compute in-place cumulative sum (row pointers!)
      CALL ComputeCRSIndexes(nvertex, vptr)

      ! Allocate vertex to element lists
      ALLOCATE(vind(vptr(nvertex+1)), STAT=allocstat)
      IF (allocstat /= 0) CALL Fatal('VertexToElementList', &
              'Vertex allocation failed!')

      ! Construct element lists for each vertex
      DO i=1,nelem
        eli = eptr(i)
        eti = eptr(i+1)-1

        ! For each vertex in element
        DO j=eli, eti
          ! Add connection to vertex eind(j)
          ind = eind(j)
          vind(vptr(ind))=i
          vptr(ind)=vptr(ind)+1
        END DO
      END DO

      ! Correct row pointers
      DO i=nvertex,2,-1
        vptr(i)=vptr(i-1)
      END DO
      vptr(1)=1
    END SUBROUTINE VertexToElementList

    ! k-way merge with an array
    SUBROUTINE kWayMergeArray(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, map)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      INTEGER :: map(:)

      INTEGER :: i, j, k, vindi

      ! Merge nv lists using a map (i.e. an array)
      nn = 1
      DO i=1,nv
        DO j=ptrli(i), ptrti(i)-1
          vindi = vind(j)
          ! Put element to map if it is not already there
          IF (map(vindi)==0 .AND. vindi /= node) THEN
            neighind(nn)=vindi
            ! Increase counter
            map(vindi)=1
            nn=nn+1
          END IF
        END DO
      END DO
      nn=nn-1

      ! Clear map
      DO i=1,nn
        map(neighind(i)) = 0
      END DO
    END SUBROUTINE kWayMergeArray

    ! k-way merge with an actual heap
    SUBROUTINE kWayMergeHeap(node, nv, ptrli, ptrti, te, vind, &
            nn, neighind, heap)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: node, nv
      INTEGER :: ptrli(:)
      INTEGER, INTENT(IN) ::ptrti(:), te
      INTEGER, INTENT(IN) :: vind(:)
      INTEGER, INTENT(OUT) :: nn
      INTEGER :: neighind(:)
      TYPE(IntTuple_t) :: heap(:)

      TYPE(IntTuple_t) :: tmp
      INTEGER :: ii, l, r, mind, ll, tmpval, tmpind

      ! Local variables
      INTEGER :: i, e, nzheap, vindi, lindi, pind

      ! Put elements to heap
      nzheap = 0
      DO i=1,nv
        IF (ptrli(i)<ptrti(i)) THEN
          heap(i) % i1 = vind(ptrli(i))
          heap(i) % i2= i
          ptrli(i) = ptrli(i)+1
          nzheap = nzheap+1
        END IF
      END DO

      ! Build heap
      DO ii=(nzheap/2), 1, -1
        i = ii
        ! CALL BinaryHeapHeapify(heap, nzheap, i)
        DO 
          ! Find index of the minimum element
          IF (2*i<=nzheap) THEN
            IF (heap(2*i) % i1 < heap(i) % i1) THEN
              mind = 2*i
            ELSE
              mind = i
            END IF
            IF (2*i+1<=nzheap) THEN
              IF (heap(2*i+1) % i1 < heap(mind) % i1) mind = 2*i+1
            END IF
          ELSE
            mind = i
          END IF

          IF (mind == i) EXIT

          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO
      END DO

      pind = -1
      nn = 1
      DO e=1,te
        ! Pick the first element from heap
        vindi = heap(1) % i1
        lindi = heap(1) % i2

        ! Remove duplicates
        IF (vindi /= pind .AND. vindi /= node) THEN
          neighind(nn) = vindi
          pind = vindi
          nn = nn+1
        END IF

        ! Add new element from list (if any)
        IF (ptrli(lindi) < ptrti(lindi)) THEN
          heap(1) % i1 = vind(ptrli(lindi))
          heap(1) % i2 = lindi
          ptrli(lindi) = ptrli(lindi)+1
        ELSE
          heap(1) % i1 = heap(nzheap) % i1
          heap(1) % i2 = heap(nzheap) % i2
          nzheap=nzheap-1
        END IF
        ! CALL BinaryHeapHeapify(heap, nzheap, 1)
        i = 1

        DO 
          ! Find the index of the minimum element
          ii = 2*i
          mind = i
          IF (ii+1<=nzheap) THEN
            ! Elements 2*i and 2*i+1 can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
            IF (heap(ii+1) % i1 < heap(mind) % i1) mind = ii+1
          ELSE IF (ii<=nzheap) THEN
            ! Element ii can be tested
            IF (heap(ii) % i1 < heap(i) % i1) mind = ii
          END IF

          IF (mind == i) EXIT

          ! Bubble down the element
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        END DO

      END DO
      nn=nn-1
    END SUBROUTINE kWayMergeHeap

    SUBROUTINE BinaryHeapHeapify(heap, nelem, sind)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      INTEGER, INTENT(IN) :: sind

      INTEGER :: i, l, r, mind
      TYPE(IntTuple_t) :: tmp

      i = sind
      DO
        l = 2*i
        r = 2*i+1
        ! Find index of the minimum element
        mind = i
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) mind = l
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(mind) % i1) mind = r
        END IF

        IF (mind /= i) THEN
          tmp = heap(i)
          heap(i) = heap(mind)
          heap(mind) = tmp
          i = mind
        ELSE
          EXIT
        END IF
      END DO
    END SUBROUTINE BinaryHeapHeapify

    FUNCTION BinaryHeapIsHeap(heap, nelem) RESULT(heaporder)
      IMPLICIT NONE
      TYPE(IntTuple_t) :: heap(:)
      INTEGER, INTENT(IN) :: nelem
      LOGICAL :: heaporder

      INTEGER :: i, l, r

      heaporder = .TRUE.

      DO i=(nelem/2), 1, -1
        l = 2*i
        r = 2*i+1
        IF (l <= nelem) THEN
          IF (heap(l) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'left: ', l, i
            EXIT
          END IF
        END IF
        IF (r <= nelem) THEN
          IF (heap(r) % i1 < heap(i) % i1) THEN
            heaporder = .FALSE.
            write (*,*) 'right: ', r, i
            EXIT
          END IF
        END IF
      END DO
    END FUNCTION BinaryHeapIsHeap

  END SUBROUTINE ElmerMeshToDualGraph

  SUBROUTINE Graph_Deallocate(Graph)
    IMPLICIT NONE
    TYPE(Graph_t) :: Graph

    DEALLOCATE(Graph % ptr)
    DEALLOCATE(Graph % ind)
    Graph % n = 0
  END SUBROUTINE Graph_Deallocate

  SUBROUTINE ElmerGraphColour(Graph, Colouring, ConsistentColours)
    IMPLICIT NONE

    TYPE(Graph_t), INTENT(IN) :: Graph
    TYPE(Graphcolour_t) :: Colouring
    LOGICAL, OPTIONAL :: ConsistentColours

    INTEGER, ALLOCATABLE :: uncolored(:)
    INTEGER, ALLOCATABLE :: fc(:), ucptr(:), rc(:), rcnew(:)

    INTEGER :: nc, dualmaxdeg, i, v, w, uci, wci, vli, vti, vcol, wcol, &
            nrc, nunc, nthr, TID, allocstat, gn
    INTEGER, ALLOCATABLE :: colours(:)
    INTEGER, PARAMETER :: VERTEX_PER_THREAD = 100
    LOGICAL :: consistent

    ! Iterative parallel greedy algorithm (Alg 2.) from 
    ! U. V. Catalyurek, J. Feo, A.H. Gebremedhin, M. Halappanavar, A. Pothen. 
    ! "Graph coloring algorithms for multi-core and massively multithreaded systems".
    ! Parallel computing, 38, 2012, pp. 576--594. 

    ! Initialize number of colours, maximum degree of graph and number of 
    ! uncolored vertices
    nc = 0
    dualmaxdeg = 0
    gn = Graph % n
    nunc = gn

    ! Check if a reproducible colouring is being requested
    consistent = .FALSE.
    IF (PRESENT(ConsistentColours)) consistent = ConsistentColours

    ! Get maximum vertex degree of the given graph
    !$OMP PARALLEL DO SHARED(Graph) &
    !$OMP PRIVATE(v) REDUCTION(max:dualmaxdeg) DEFAULT(NONE)
    DO v=1,Graph % n
      dualmaxdeg = MAX(dualmaxdeg, Graph % ptr(v+1)- Graph % ptr(v))
    END DO
    !$OMP END PARALLEL DO

    nthr = 1
    ! Ensure that each vertex has at most one thread attached to it
    !$ IF (.NOT. consistent) nthr = MIN(omp_get_max_threads(), gn)

    ! Allocate memory for colours of vertices and thread colour pointers
    ALLOCATE(colours(gn), uncolored(gn), ucptr(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate colour maps!')

    !$OMP PARALLEL SHARED(gn, dualmaxdeg, Graph, colours, nunc, &
    !$OMP                 uncolored, ucptr, nthr) &
    !$OMP PRIVATE(uci, vli, vti, v, w, wci, vcol, wcol, fc, nrc, rc, rcnew, &
    !$OMP         allocstat, TID) &
    !$OMP REDUCTION(max:nc) DEFAULT(NONE) NUM_THREADS(nthr)

    TID=1
    !$ TID=OMP_GET_THREAD_NUM()+1

    ! Greedy algorithm colours a given graph with at 
    ! most max_{v\in V} deg(v)+1 colours
    ALLOCATE(fc(dualmaxdeg+1), rc((gn/nthr)+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
            'Unable to allocate local workspace!')
    ! Initialize forbidden colour array (local to thread)
    fc = 0

    ! Initialize colours and uncolored entries
    !$OMP DO 
    DO v=1,gn
      colours(v)=0
      ! U <- V
      uncolored(v)=v
    END DO
    !$OMP END DO

    DO
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1

        ! For each w\in adj(v) do
        DO w=vli, vti
          ! fc[colour[w]]<-v
          !$OMP ATOMIC READ
          wcol = colours(Graph % ind(w))
          IF (wcol /= 0) fc(wcol) = v
        END DO

        ! Find smallest permissible colour for vertex
        ! c <- min\{i>0: fc[i]/=v \}
        DO i=1,dualmaxdeg+1
          IF (fc(i) /= v) THEN
            !$OMP ATOMIC WRITE 
            colours(v) = i
            ! Maintain maximum colour
            nc = MAX(nc, i)
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO

      nrc = 0
      ! For each v\in U in parallel do
      !$OMP DO
      DO uci=1,nunc
        v = uncolored(uci)
        vli = Graph % ptr(v)
        vti = Graph % ptr(v+1)-1
        vcol = colours(v)

        ! Make sure that recolour array has enough storage for 
        ! the worst case (all elements need to be added)
        IF (SIZE(rc)<nrc+(vti-vli)+1) THEN
          ALLOCATE(rcnew(MAX(SIZE(rc)*2, nrc+(vti-vli)+1)), STAT=allocstat)
          IF (allocstat /= 0) CALL Fatal('ElmerDualGraphColour', &
                  'Unable to allocate local workspace!')
          rcnew(1:nrc)=rc(1:nrc)
          DEALLOCATE(rc)
          CALL MOVE_ALLOC(rcnew, rc)
        END IF

        ! For each w\in adj(v) do
        DO wci=vli,vti
          w = Graph % ind(wci)
          IF (colours(w)==vcol .AND. v>w) THEN
            ! R <- R\bigcup {v} (thread local)
            nrc = nrc + 1
            rc(nrc)=v
            EXIT
          END IF
        END DO
      END DO
      !$OMP END DO NOWAIT

      ucptr(TID)=nrc
      !$OMP BARRIER

      !$OMP SINGLE
      CALL ComputeCRSIndexes(nthr, ucptr)
      nunc = ucptr(nthr+1)-1
      !$OMP END SINGLE

      ! U <- R
      uncolored(ucptr(TID):ucptr(TID+1)-1)=rc(1:nrc)
      !$OMP BARRIER

      ! Colour the remaining vertices sequentially if the 
      ! size of the set of uncoloured vertices is small enough
      IF (nunc < nthr*VERTEX_PER_THREAD) THEN
        !$OMP SINGLE
        DO uci=1,nunc
          v = uncolored(uci)
          vli = Graph % ptr(v)
          vti = Graph % ptr(v+1)-1

          ! For each w\in adj(v) do
          DO w=vli, vti
            ! fc[colour[w]]<-v
            wcol = colours(Graph % ind(w))
            IF (wcol /= 0) fc(wcol) = v
          END DO

          ! Find smallest permissible colour for vertex
          ! c <- min\{i>0: fc[i]/=v \}
          DO i=1,dualmaxdeg+1
            IF (fc(i) /= v) THEN
              ! Single thread, no collisions possible 
              colours(v) = i
              ! Maintain maximum colour
              nc = MAX(nc, i)
              EXIT
            END IF
          END DO
        END DO
        !$OMP END SINGLE NOWAIT

        EXIT
      END IF

    END DO

    ! Deallocate thread local storage
    DEALLOCATE(fc, rc)
    !$OMP END PARALLEL

    DEALLOCATE(uncolored, ucptr)

    ! Set up colouring data structure
    Colouring % nc = nc
    CALL MOVE_ALLOC(colours, Colouring % colours)
  END SUBROUTINE ElmerGraphColour

  SUBROUTINE Colouring_Deallocate(Colours)
    IMPLICIT NONE
    TYPE(GraphColour_t) :: Colours

    DEALLOCATE(Colours % colours)
    Colours % nc = 0
  END SUBROUTINE Colouring_Deallocate

  SUBROUTINE ElmerColouringToGraph(Colours, PackedList)
    IMPLICIT NONE

    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(Graph_t) :: PackedList

    INTEGER, ALLOCATABLE :: cptr(:), cind(:)

    INTEGER :: nc, c, i, n, allocstat

    nc = Colours % nc
    n = size(Colours % colours)
    ALLOCATE(cptr(nc+1), cind(n), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ElmerGatherColourLists','Memory allocation failed.')
    cptr = 0
    ! Count number of elements in each colour
    DO i=1,n
      cptr(Colours % colours(i))=cptr(Colours % colours(i))+1
    END DO

    CALL ComputeCRSIndexes(nc, cptr)

    DO i=1,n
      c=Colours % colours(i)
      cind(cptr(c))=i
      cptr(c)=cptr(c)+1
    END DO

    DO i=nc,2,-1
      cptr(i)=cptr(i-1)
    END DO
    cptr(1)=1

    ! Set up graph data structure
    PackedList % n = nc
    CALL MOVE_ALLOC(cptr, PackedList % ptr)
    CALL MOVE_ALLOC(cind, PackedList % ind)
  END SUBROUTINE ElmerColouringToGraph

  ! Routine constructs colouring for boundary mesh based on colours of main mesh
  SUBROUTINE ElmerBoundaryGraphColour(Mesh, Colours, BoundaryColours)
    IMPLICIT NONE

    TYPE(Mesh_t), INTENT(IN) :: Mesh
    TYPE(GraphColour_t), INTENT(IN) :: Colours
    TYPE(GraphColour_t) :: BoundaryColours

    TYPE(Element_t), POINTER :: Element
    INTEGER :: elem, nelem, nbelem, astat, lcolour, rcolour, nbc
    INTEGER, ALLOCATABLE :: bcolours(:)

    nelem = Mesh % NumberOfBulkElements
    nbelem = Mesh % NumberOfBoundaryElements

    ! Allocate boundary colouring
    ALLOCATE(bcolours(nbelem), STAT=astat)
    IF (astat /= 0) THEN
       CALL Fatal('ElmerBoundaryGraphColour','Unable to allocate boundary colouring')
    END IF
    
    nbc = 0
    ! Loop over boundary mesh
    !$OMP PARALLEL DO &
    !$OMP SHARED(Mesh, nelem, nbelem, Colours, bcolours) &
    !$OMP PRIVATE(Element, lcolour, rcolour) &
    !$OMP REDUCTION(max:nbc) &
    !$OMP DEFAULT(NONE)
    DO elem=1,nbelem       
       Element => Mesh % Elements(nelem+elem)

       ! Try to find colour for boundary element based on left / right parent
       lcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Left)) THEN
          lcolour = Colours % colours(Element % BoundaryInfo % Left % ElementIndex)
       END IF
       rcolour = 0
       IF (ASSOCIATED(Element % BoundaryInfo % Right)) THEN
          rcolour = Colours % colours(Element % BoundaryInfo % Right % ElementIndex)
       END IF

       ! Sanity check for debug
       IF (ASSOCIATED(Element % BoundaryInfo % Left) .AND. & 
          ASSOCIATED(Element % BoundaryInfo % Right) .AND. &
            lcolour /= rcolour) THEN
         CALL Warn('ElmerBoundaryGraphColour','Inconsistent colours for boundary element: ' &
               // TRIM(i2s(elem)) // "=>" &
               // TRIM(i2s(lcolour))// " | "//TRIM(i2s(rcolour)))
         WRITE (*,*) Element % BoundaryInfo % Left % ElementIndex, Element % BoundaryInfo % Right % ElementIndex
       END IF

       bcolours(elem)=MAX(lcolour,rcolour)
       nbc=MAX(nbc,bcolours(elem))
    END DO
    !$OMP END PARALLEL DO

    ! Set up colouring data structure
    BoundaryColours % nc = nbc
    CALL MOVE_ALLOC(bcolours, BoundaryColours % colours)
  END SUBROUTINE ElmerBoundaryGraphColour
  
  ! Given CRS indices, referenced indirectly from graph, 
  ! evenly load balance the work among the nthr threads
  SUBROUTINE ThreadLoadBalanceElementNeighbour(nthr, gn, gptr, gind, &
          rptr, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER :: gptr(:), gind(:), rptr(:)
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, j, k, wrk, gwrk, thrwrk, allocstat

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadLoadBalanceElementNeighbour', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Compute total global work
    gwrk = 0
    DO i=1,gn
      DO j=gptr(i),gptr(i+1)-1
        gwrk = gwrk + (rptr(gind(j)+1)-rptr(gind(j)))
      END DO
    END DO

    ! Amount of work per thread
    thrwrk = CEILING(REAL(gwrk,dp) / nthr)

    ! Find rows for each thread to compute
    blkleads(1)=1
    DO i=1,nthr
      wrk = 0
      ! Acquire enough work for thread i
      DO j=blkleads(i),gn
        DO k=gptr(j),gptr(j+1)-1
          wrk = wrk + (rptr(gind(j)+1)-rptr(gind(j)))
        END DO
        IF (wrk >= thrwrk) EXIT
      END DO

      blkleads(i+1)=j+1
      ! Check if we have run out of rows
      IF (j+1>gn) EXIT
    END DO
    ! Reset number of rows (may be less than or equal to original number)
    nthr = i
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadLoadBalanceElementNeighbour

  SUBROUTINE ThreadStaticWorkShare(nthr, gn, blkleads)
    IMPLICIT NONE

    INTEGER :: nthr
    INTEGER, INTENT(IN) :: gn
    INTEGER, ALLOCATABLE :: blkleads(:)

    INTEGER :: i, rem, thrwrk, allocstat
    INTEGER :: totelem

    ! Compute number of nonzeroes / thread
    !$ nthr = MIN(nthr,gn)

    ALLOCATE(blkleads(nthr+1), STAT=allocstat)
    IF (allocstat /= 0) CALL Fatal('ThreadStaticWorkShare', &
            'Unable to allocate blkleads!')

    ! Special case of just one thread
    IF (nthr == 1) THEN
      blkleads(1)=1
      blkleads(2)=gn+1
      RETURN
    END IF

    ! Assuming even distribution of nodes / element, 
    ! distribute rows for each thread to compute 
    blkleads(1)=1
    thrwrk = gn / nthr
    rem = gn-nthr*thrwrk
    ! totelem = 0
    DO i=1,nthr-1
      IF (i<rem) THEN
        blkleads(i+1)=blkleads(i)+thrwrk+1
      ELSE
        blkleads(i+1)=blkleads(i)+thrwrk
      END IF
    END DO
    ! Assign what is left of the matrix to the final thread
    blkleads(nthr+1)=gn+1
  END SUBROUTINE ThreadStaticWorkShare

  ! Given row counts, in-place compute CRS indices to data
  SUBROUTINE ComputeCRSIndexes(n, arr)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: n
    INTEGER :: arr(:)

    INTEGER :: i, indi, indip

    indi = arr(1)
    arr(1)=1
    DO i=1,n-1
      indip=arr(i+1)
      arr(i+1)=arr(i)+indi
      indi=indip
    END DO
    arr(n+1)=arr(n)+indi
  END SUBROUTINE ComputeCRSIndexes

  !> Calcalate body average for a discontinuous galerkin field.
  !> The intended use is in conjunction of saving the results. 
  !> This tampers the field and therefore may have unwanted side effects
  !> if the solution is to be used for something else too.
  !-------------------------------------------------------------------
  SUBROUTINE CalculateBodyAverage( Mesh, Var, BodySum )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    LOGICAL :: BodySum

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: BodyAverage(:)
    INTEGER, ALLOCATABLE :: BodyCount(:)
    INTEGER :: n,i,j,k,l,nodeind,dgind, Nneighbours
    REAL(KIND=dp) :: AveHits
    LOGICAL, ALLOCATABLE :: IsNeighbour(:)

    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    IF( BodySum ) THEN
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal sum for: '&
          //TRIM(Var % Name), Level=8)
    ELSE
      CALL Info('CalculateBodyAverage','Calculating bodywise nodal average for: '&
          //TRIM(Var % Name), Level=8)
    END IF

    n = Mesh % NumberOfNodes
    ALLOCATE( BodyCount(n), BodyAverage(n), IsNeighbour(Parenv % PEs) )


    DO i=1,CurrentModel % NumberOfBodies

      DO k=1,Var % Dofs
        BodyCount = 0
        BodyAverage = 0.0_dp

        DO j=1,Mesh % NumberOfBulkElements 
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE
          DO l = 1, Element % TYPE % NumberOfNodes
            nodeind = Element % NodeIndexes(l)
            dgind = Var % Perm(Element % DGIndexes(l) )
            IF( dgind > 0 ) THEN
              BodyAverage( nodeind ) = BodyAverage( nodeind ) + &
                  Var % Values( Var % DOFs*( dgind-1)+k )
              BodyCount( nodeind ) = BodyCount( nodeind ) + 1 
            END IF
          END DO
        END DO

        IF( k == 1 ) THEN
          AveHits = 1.0_dp * SUM( BodyCount ) / COUNT( BodyCount > 0 )
          !PRINT *,'AveHits:',i,AveHits
        END IF

        IF(ParEnv % Pes>1) THEN
          Nneighbours = MeshNeighbours(Mesh, IsNeighbour)
          CALL SendInterface(); CALL RecvInterface()
        END IF

        ! Do not average weighted quantities. They should only be summed, I guess... 
        
        IF( .NOT. BodySum ) THEN
          DO j=1,n
            IF( BodyCount(j) > 0 ) BodyAverage(j) = BodyAverage(j) / BodyCount(j)
          END DO
        END IF

        DO j=1,Mesh % NumberOfBulkElements 
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE
          DO l = 1, Element % TYPE % NumberOfNodes
            nodeind = Element % NodeIndexes(l)
            dgind = Var % Perm(Element % DGIndexes(l) )
            IF( dgind > 0 ) THEN
              Var % Values( Var % DOFs*( dgind-1)+k ) = BodyAverage( nodeind ) 
            END IF
          END DO
        END DO
      END DO
    END DO

CONTAINS

     SUBROUTINE SendInterface()
       TYPE buf_t
         REAL(KIND=dp), ALLOCATABLE :: dval(:)
         INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       END TYPE buf_t

       INTEGER, ALLOCATABLE :: cnt(:)
       TYPE(buf_t), ALLOCATABLE :: buf(:)

       INTEGER :: i,j,k,ierr

       ALLOCATE(cnt(ParEnv % PEs), buf(ParEnv % PEs))

       cnt = 0
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT.Mesh % ParallelInfo % Interface(i)) CYCLE
         IF(BodyCount(i) <= 0 ) CYCLE

         DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
           k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)+1
           cnt(k) = cnt(k) + 1
         END DO
       END DO

       DO i=1,ParEnv % PEs
         ALLOCATE(buf(i) % gdof(cnt(i)), buf(i) % ival(cnt(i)), buf(i) % dval(cnt(i)))
       END DO

       cnt = 0
       DO i=1,Mesh % NumberOfNodes
         IF(.NOT.Mesh % ParallelInfo % Interface(i)) CYCLE
         IF(BodyCount(i) <= 0 ) CYCLE

         DO j=1,SIZE(Mesh % ParallelInfo % NeighbourList(i) % Neighbours)
           k = Mesh % ParallelInfo % NeighbourList(i) % Neighbours(j)+1
           cnt(k) = cnt(k) + 1
           buf(k) % gdof(cnt(k)) = Mesh % ParallelInfo % GlobalDOFs(i)
           buf(k) % ival(cnt(k)) = BodyCount(i)
           buf(k) % dval(cnt(k)) = BodyAverage(i)
         END DO
       END DO

       DO i=1,ParEnv % PEs
         IF(.NOT. isNeighbour(i)) CYCLE

         CALL MPI_BSEND( cnt(i),1,MPI_INTEGER,i-1,1310,ELMER_COMM_WORLD,ierr )
         IF(cnt(i)>0) THEN
           CALL MPI_BSEND( buf(i) % gdof,cnt(i),MPI_INTEGER,i-1,1311,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % ival,cnt(i),MPI_INTEGER,i-1,1312,ELMER_COMM_WORLD,ierr )
           CALL MPI_BSEND( buf(i) % dval,cnt(i),MPI_DOUBLE_PRECISION,i-1,1313,ELMER_COMM_WORLD,ierr )
         END IF
       END DO
     END SUBROUTINE SendInterface


     SUBROUTINE RecvInterface()
       INTEGER, ALLOCATABLE :: gdof(:), ival(:)
       REAL(KIND=dp), ALLOCATABLE :: dval(:)
       INTEGER :: i,j,k,ierr, cnt, status(MPI_STATUS_SIZE)

       DO i=1,ParEnv % PEs

         IF(.NOT.isNeighbour(i)) CYCLE

         CALL MPI_RECV( cnt,1,MPI_INTEGER,i-1,1310,ELMER_COMM_WORLD,status,ierr )
         IF(cnt>0) THEN
           ALLOCATE( gdof(cnt), ival(cnt), dval(cnt) )
           CALL MPI_RECV( gdof,cnt,MPI_INTEGER,i-1,1311,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( ival,cnt,MPI_INTEGER,i-1,1312,ELMER_COMM_WORLD,status,ierr )
           CALL MPI_RECV( dval,cnt,MPI_DOUBLE_PRECISION,i-1,1313,ELMER_COMM_WORLD,status,ierr )

           DO j=1,cnt
             k = SearchNode(Mesh % ParallelInfo, gdof(j))
             IF (k>0) THEN
               BodyCount(k) = BodyCount(k) + ival(j)
               BodyAverage(k) = BodyAverage(k)  + dval(j)
             END IF
           END DO 
           DEALLOCATE( gdof, ival, dval )
         END IF
       END DO
       CALL MPI_BARRIER(ELMER_COMM_WORLD,ierr)
     END SUBROUTINE RecvInterface

  END SUBROUTINE CalculateBodyAverage



  !> Given an elemental DG field create a minimal reduced set of it that maintains
  !> the necessary continuities. The continuities may be requested between bodies
  !> or materials. Optionally the user may give a boundary mask which defines the 
  !> potential discontinuous nodes that may be greedy or not. 
  !-------------------------------------------------------------------------------
  FUNCTION MinimalElementalSet( Mesh, JumpMode, VarPerm, BcFlag, &
      NonGreedy ) RESULT ( SetPerm )

    TYPE(Mesh_t), POINTER :: Mesh
    CHARACTER(LEN=*) :: JumpMode
    INTEGER, POINTER, OPTIONAL :: VarPerm(:)
    CHARACTER(LEN=*), OPTIONAL :: BcFlag
    LOGICAL, OPTIONAL :: NonGreedy
    INTEGER, POINTER :: SetPerm(:)

    TYPE(Element_t), POINTER :: Element, Left, Right
    INTEGER :: n,i,j,k,l,bc_id,mat_id,body_id,NoElimNodes,nodeind,JumpModeIndx,&
        LeftI,RightI,NumberOfBlocks
    LOGICAL, ALLOCATABLE :: JumpNodes(:)
    INTEGER, ALLOCATABLE :: NodeVisited(:)
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Found
    

    CALL Info('MinimalDiscontSet','Creating discontinuous subset from DG field',Level=5)

    ! Calculate size of permutation vector
    ALLOCATE( NodeVisited( Mesh % NumberOfNodes ) )
    NodeVisited = 0

    NULLIFY( SetPerm ) 
    k = 0
    DO i=1,Mesh % NumberOfBulkElements         
      Element => Mesh % Elements(i)
      k = k + Element % TYPE % NumberOfNodes
    END DO
    CALL Info('MinimalElementalSet','Maximum number of dofs in DG: '//TRIM(I2S(k)),Level=12)
    ALLOCATE( SetPerm(k) )
    SetPerm = 0
    l = 0
    NoElimNodes = 0

    CALL Info('MinimalElementalSet','Reducing elemental discontinuity with mode: '//TRIM(JumpMode),Level=7)

    SELECT CASE ( JumpMode )

    CASE('db') ! discontinuous bodies
      NumberOfBlocks = CurrentModel % NumberOfBodies
      JumpModeIndx = 1

    CASE('dm') ! discontinuous materials
      NumberOfBlocks = CurrentModel % NumberOfMaterials
      JumpModeIndx = 2

    CASE DEFAULT
      CALL Fatal('MinimalElementalSet','Unknown JumpMode: '//TRIM(JumpMode))

    END SELECT
  

    IF( PRESENT( BcFlag ) ) THEN
      ALLOCATE( JumpNodes( Mesh % NumberOfNodes ) )
    END IF

    
    DO i=1,NumberOfBlocks
      
      ! Before the 1st block no numbers have been given.
      ! Also if we want discontinuous blocks on all sides initialize the whole list to zero. 
      IF( i == 1 .OR. .NOT. PRESENT( BcFlag ) ) THEN
        NodeVisited = 0

      ELSE
        ! Vector indicating the disontinuous nodes
        ! If this is not given all interface nodes are potentially discontinuous
        JumpNodes = .FALSE.
        
        DO j=Mesh % NumberOfBulkElements + 1, &
            Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
          Element => Mesh % Elements(j)

          DO bc_id=1,CurrentModel % NumberOfBCs
            IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
          END DO
          IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE
          IF( .NOT. ListCheckPresent( CurrentModel % BCs(bc_id) % Values, BcFlag ) ) CYCLE

          Left => Element % BoundaryInfo % Left
          Right => Element % BoundaryInfo % Right
          IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) CYCLE

          IF( JumpModeIndx == 1 ) THEN
            LeftI = Left % BodyId
            RightI = Right % BodyId
          ELSE
            LeftI = ListGetInteger( CurrentModel % Bodies(Left % BodyId) % Values,'Material',Found)
            RightI = ListGetInteger( CurrentModel % Bodies(Right % BodyId) % Values,'Material',Found)
          END IF

          IF( LeftI /= i .AND. RightI /= i ) CYCLE
          JumpNodes( Element % NodeIndexes ) = .TRUE.
        END DO

        IF( PRESENT( NonGreedy ) ) THEN
          IF( NonGreedy ) THEN        
            DO j=Mesh % NumberOfBulkElements + 1, &
                Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements
              Element => Mesh % Elements(j)

              DO bc_id=1,CurrentModel % NumberOfBCs
                IF ( Element % BoundaryInfo % Constraint == CurrentModel % BCs(bc_id) % Tag ) EXIT
              END DO
              IF ( bc_id > CurrentModel % NumberOfBCs ) CYCLE

              IF( ListCheckPresent( CurrentModel % BCs(bc_id) % Values, BcFlag ) ) CYCLE

              Left => Element % BoundaryInfo % Left
              Right => Element % BoundaryInfo % Right

              ! External BCs don't have a concept of jump, so no need to treat them
              IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) CYCLE

              JumpNodes( Element % NodeIndexes ) = .FALSE.
            END DO
          END IF
        END IF

        ! Initialize new potential nodes for the block where we found discontinuity
        WHERE( JumpNodes ) NodeVisited = 0
      END IF


      ! Now do the real thing. 
      ! Add new dofs such that minimal discontinuity is maintained 
      DO j=1,Mesh % NumberOfBulkElements         
        Element => Mesh % Elements(j)

        Body_Id = Element % BodyId 
        IF( JumpModeIndx == 1 ) THEN
          IF( Body_id /= i ) CYCLE
        ELSE
          Mat_Id = ListGetInteger( CurrentModel % Bodies(Body_Id) % Values,'Material',Found)
          IF( Mat_Id /= i ) CYCLE
        END IF

        NodeIndexes => Element % NodeIndexes
        
        DO k=1,Element % TYPE % NumberOfNodes         
          nodeind = NodeIndexes(k)
          IF( PRESENT( VarPerm ) ) THEN
            IF( VarPerm( nodeind ) == 0 ) CYCLE
          END IF
          IF( NodeVisited( nodeind ) > 0 ) THEN
            SetPerm( Element % DGIndexes(k) ) = NodeVisited( nodeind )
            NoElimNodes = NoElimNodes + 1
          ELSE
            l = l + 1
            NodeVisited(nodeind) = l
            SetPerm( Element % DGIndexes(k) ) = l
          END IF
        END DO
      END DO
    END DO

    CALL Info('MinimalElementalSet','Independent dofs in elemental field: '//TRIM(I2S(l)),Level=7)
    CALL Info('MinimalElementalSet','Redundant dofs in elemental field: '//TRIM(I2S(NoElimNodes)),Level=7)     

  END FUNCTION MinimalElementalSet


  !> Calculate the reduced DG field given the reduction permutation.
  !> The permutation must be predefined. This may be called repeatedly
  !> for different variables. Optionally one may take average, or 
  !> a plain sum over the shared nodes. 
  !-------------------------------------------------------------------
  SUBROUTINE ReduceElementalVar( Mesh, Var, SetPerm, TakeAverage )

    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: SetPerm(:)
    LOGICAL :: TakeAverage

    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp), ALLOCATABLE :: SetSum(:)
    INTEGER, ALLOCATABLE :: SetCount(:)
    INTEGER :: dof,n,m,i,j,k,l,nodeind,dgind
    REAL(KIND=dp) :: AveHits

    IF(.NOT. ASSOCIATED(var)) THEN
      CALL Warn('ReduceElementalVar','Variable not associated!')
      RETURN
    END IF

    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) THEN
      CALL Warn('ReduceElementalVar','Var % Perm too small!')
      RETURN
    END IF

    IF( TakeAverage ) THEN
      CALL Info('CalculateSetAverage','Calculating reduced set average for: '&
          //TRIM(Var % Name), Level=7)
    ELSE
      CALL Info('CalculateSetAverage','Calculating reduced set sum for: '&
          //TRIM(Var % Name), Level=7)
    END IF

    n = Mesh % NumberOfNodes

    m = MAXVAL( SetPerm )
    ALLOCATE( SetCount(m), SetSum(m) )
    SetCount = 0
    SetSum = 0.0_dp

    ! Take the sum to nodes, and calculate average if requested
    DO dof=1,Var % Dofs
      SetCount = 0
      SetSum = 0.0_dp

      DO i=1,SIZE(SetPerm)
        j = SetPerm(i)
        l = Var % Perm(i)
        SetSum(j) = SetSum(j) + Var % Values( Var % DOFs * (l-1) + dof )
        SetCount(j) = SetCount(j) + 1
      END DO
        
      m = SUM( SetCount ) 
      IF( m == 0 ) RETURN

      IF( TakeAverage ) THEN
        WHERE( SetCount > 0 ) SetSum = SetSum / SetCount
      END IF

      IF( dof == 1 ) THEN
        AveHits = 1.0_dp * SUM( SetCount ) / COUNT( SetCount > 0 )
        WRITE(Message,'(A,ES15.4)') 'Average number of hits: ',AveHits
        CALL Info('ReduceElementalVar',Message,Level=10)
      END IF

      ! Copy the reduced set back to the original elemental field
      DO i=1,SIZE(SetPerm)
        j = SetPerm(i)
        l = Var % Perm(i)
        Var % Values( Var % DOFs * (l-1) + dof ) = SetSum(j)
      END DO
    END DO

  END SUBROUTINE ReduceElementalVar


  !> Given a elemental DG field and a reduction permutation compute the 
  !> body specific lumped sum. The DG field may be either original one
  !> or already summed up. In the latter case only one incident of the 
  !> redundant nodes is set.
  !---------------------------------------------------------------------
  SUBROUTINE LumpedElementalVar( Mesh, Var, SetPerm, AlreadySummed )
    TYPE(Variable_t), POINTER :: Var
    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER, POINTER :: SetPerm(:)
    LOGICAL :: AlreadySummed

    TYPE(Element_t), POINTER :: Element
    LOGICAL, ALLOCATABLE :: NodeVisited(:)
    INTEGER :: dof,n,m,i,j,k,l,nodeind,dgind
    REAL(KIND=dp), ALLOCATABLE :: BodySum(:)

    IF(.NOT. ASSOCIATED(var)) RETURN
    IF( SIZE(Var % Perm) <= Mesh % NumberOfNodes ) RETURN

    CALL Info('LumpedElementalVar','Calculating lumped sum for: '&
        //TRIM(Var % Name), Level=8)

    n = Mesh % NumberOfNodes

    m = MAXVAL( SetPerm )
    IF( AlreadySummed ) THEN
      ALLOCATE( NodeVisited(m) )
    END IF
    ALLOCATE( BodySum( CurrentModel % NumberOfBodies ) )

    ! Take the sum to nodes, and calculate average if requested
    DO dof=1,Var % Dofs

      BodySum = 0.0_dp

      DO i=1,CurrentModel % NumberOfBodies

        IF( AlreadySummed ) THEN
          NodeVisited = .FALSE.
        END IF

        DO j=1,Mesh % NumberOfBulkElements         
          Element => Mesh % Elements(j)
          IF( Element % BodyId /= i ) CYCLE

          DO k=1,Element % TYPE % NumberOfNodes         
            dgind = Element % DGIndexes(k)
            l = SetPerm(dgind)
            IF( l == 0 ) CYCLE

            IF( AlreadySummed ) THEN
              IF( NodeVisited(l) ) CYCLE           
              NodeVisited(l) = .TRUE.
            END IF

            BodySum(i) = BodySum(i) + &
                Var % Values( Var % Dofs * ( Var % Perm( dgind )-1) + dof )
          END DO
        END DO
      END DO

      IF( Var % Dofs > 1 ) THEN
        CALL Info('LumpedElementalVar','Lumped sum for component: '//TRIM(I2S(dof)),Level=6)
      END IF
      DO i=1,CurrentModel % NumberOfBodies
        WRITE(Message,'(A,ES15.4)') 'Body '//TRIM(I2S(i))//' sum:',BodySum(i)
        CALL Info('LumpedElementalVar',Message,Level=10)
      END DO

    END DO

    DEALLOCATE( NodeVisited, BodySum )

  END SUBROUTINE LumpedElementalVar



!------------------------------------------------------------------------------
  SUBROUTINE SaveParallelInfo( Solver )
!------------------------------------------------------------------------------
   TYPE( Solver_t ), POINTER  :: Solver
!------------------------------------------------------------------------------    
   TYPE(ParallelInfo_t), POINTER :: ParInfo=>NULL()
   TYPE(ValueList_t), POINTER :: Params
   CHARACTER(LEN=MAX_NAME_LEN) :: dumpfile
   INTEGER :: i,j,k,n,maxnei
   LOGICAL :: Found, MeshMode, MatrixMode
   CHARACTER(*), PARAMETER :: Caller = "SaveParallelInfo"
   TYPE(Nodes_t), POINTER :: Nodes
   
   Params => Solver % Values 

   MeshMode = ListGetLogical( Params,'Save Parallel Matrix Info',Found ) 
   MatrixMode = ListGetLogical( Params,'Save Parallel Mesh Info',Found ) 

   IF( .NOT. ( MeshMode .OR. MatrixMode ) ) RETURN

10 IF( MeshMode ) THEN
     CALL Info(Caller,'Saving parallel mesh info',Level=8 ) 
   ELSE
     CALL Info(Caller,'Saving parallel matrix info',Level=8 ) 
   END IF

   IF( MeshMode ) THEN
     ParInfo => Solver % Mesh % ParallelInfo
     Nodes => Solver % Mesh % Nodes
     dumpfile = 'parinfo_mesh.dat'
   ELSE
     ParInfo => Solver % Matrix % ParallelInfo
     dumpfile = 'parinfo_mat.dat'      
   END IF

   IF( .NOT. ASSOCIATED( ParInfo ) ) THEN
     CALL Warn(Caller,'Parallel info not associated!')
     RETURN
   END IF

   n = SIZE( ParInfo % GlobalDOFs )
   IF( n <= 0 ) THEN
     CALL Warn(Caller,'Parallel info size is invalid!')
     RETURN
   END IF

   ! memorize the maximum number of parallel neighbours
   maxnei = 0
   IF( ASSOCIATED( ParInfo % NeighbourList ) ) THEN
     DO i=1,n
       IF( ASSOCIATED( ParInfo % NeighbourList(i) % Neighbours ) ) THEN
         j = SIZE( ParInfo % NeighbourList(i) % Neighbours )
         maxnei = MAX( j, maxnei ) 
       END IF
     END DO
   END IF
   CALL Info(Caller,'Maximum number of parallel neighbours:'//TRIM(I2S(maxnei)))

   IF(ParEnv % PEs > 1) dumpfile = TRIM(dumpfile)//'.'//TRIM(I2S(ParEnv % myPE))      
   CALL Info(Caller,'Saving parallel info to: '//TRIM(dumpfile),Level=8)

   OPEN(1,FILE=dumpfile, STATUS='Unknown')  
   DO i=1,n
     j = ParInfo % GlobalDOFs(i)
     IF( ParInfo % INTERFACE(i) ) THEN
       k = 1
     ELSE
       k = 0
     END IF
     WRITE(1,'(3I6)',ADVANCE='NO') i,j,k
     IF( ASSOCIATED( ParInfo % NeighbourList(i) % Neighbours ) ) THEN
       k = SIZE( ParInfo % NeighbourList(i) % Neighbours )
     ELSE
       k = 0
     END IF
     DO j=1,k
       WRITE(1,'(I6)',ADVANCE='NO')  ParInfo % NeighbourList(i) % Neighbours(j)
     END DO
     DO j=k+1,maxnei
       WRITE(1,'(I6)',ADVANCE='NO')  -1 
     END DO
     IF( MeshMode ) THEN
       WRITE(1,'(3ES12.3)',ADVANCE='NO') &
           Nodes % x(i), Nodes % y(i), Nodes % z(i)
     END IF
     WRITE(1,'(A)') ' ' ! finish the line
   END DO
   CLOSE(1)

   ! Redo with matrix if both modes are requested
   IF( MeshMode .AND. MatrixMode ) THEN
     MeshMode = .FALSE.
     GOTO 10
   END IF
   
   CALL Info(Caller,'Finished saving parallel info',Level=10)

!------------------------------------------------------------------------------
 END SUBROUTINE SaveParallelInfo
!------------------------------------------------------------------------------

  
!------------------------------------------------------------------------------
END MODULE MeshUtils
!------------------------------------------------------------------------------

!> \}

