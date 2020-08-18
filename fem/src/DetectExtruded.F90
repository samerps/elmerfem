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
! *  Email:   elmeradm@csc.fi
! *  Web:     http://www.csc.fi/elmer
! *  Address: CSC - IT Center for Science Ltd.
! *           Keilaranta 14
! *           02101 Espoo, Finland 
! *
! *****************************************************************************/
  
!> \ingroup ElmerLib
!> \{
!------------------------------------------------------------------------------
!>  Utilities for mortar and other routines involving boundary meshes.
!------------------------------------------------------------------------------

MODULE DetectExtruded

  USE ElementUtils
  USE ElementDescription
  USE MeshUtils
  USE Types
  IMPLICIT NONE

CONTAINS

!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh even though it is 
!> given in an unstructured format. The routine may be used by some special
!> solvers that employ the special character of the mesh.
!> The extrusion is found for a given direction and for each node the corresponding 
!> up and down, and thereafter top and bottom node is computed.
!-----------------------------------------------------------------------------
  SUBROUTINE DetectExtrudedStructure( Mesh, Solver, ExtVar, &
      TopNodePointer, BotNodePointer, UpNodePointer, DownNodePointer, &
      MidNodePointer, MidLayerExists, NumberOfLayers, NodeLayer )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopNodePointer(:), BotNodePointer(:), &
        UpNodePointer(:), DownNodePointer(:), MidNodePointer(:)
    INTEGER, POINTER, OPTIONAL :: NodeLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
    LOGICAL, OPTIONAL :: MidLayerExists
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, nnodes, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
	UpHit, DownHit, bc_ind, jmin, jmax
    INTEGER, POINTER :: NodeIndexes(:), MaskPerm(:)
    LOGICAL :: MaskExists, UpActive, DownActive, GotIt, Found, DoCoordTransform
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
    REAL(KIND=dp) :: at0, at1, Length, UnitVector(3), Vector(3), Vector2(3), &
        ElemVector(3), DotPro, MaxDotPro, MinDotPro, Eps, MinTop, &
        MaxTop, MinBot, MaxBot
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName, CoordTransform
    CHARACTER(*), PARAMETER :: Caller='DetectExtrudedStructure'

    
    CALL Info(Caller,'Determining extruded structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim
    Params => Solver % Values

    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal('StructuredMeshMapper','Invalid value for Active Coordinate')
    END IF  
    UnitVector = 0.0_dp
    UnitVector(ActiveDirection) = 1.0_dp


    IF( ListGetLogical(Params,'Project To Bottom',GotIt) ) &
        UnitVector = -1.0_dp * UnitVector

    WRITE(Message,'(A,3F8.3)') 'Unit vector of direction:',UnitVector
    CALL Info(Caller,Message,Level=8)

    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0d-4

    VarName = ListGetString(Params,'Mapping Mask Variable',GotIt )
    MaskExists = .FALSE.
    IF(GotIt) THEN
      Var => VariableGet( Mesh % Variables,  VarName )
      IF(ASSOCIATED(Var)) THEN
        MaskExists = ASSOCIATED(Var % Perm)
        IF( MaskExists ) THEN
          ALLOCATE( MaskPerm( SIZE( Var % Perm ) ) )
          MaskPerm = Var % Perm 
          CALL Info(Caller,&
              'Using variable as mask: '//TRIM(VarName),Level=8)
        END IF
      END IF      
    END IF

    nnodes = Mesh % NumberOfNodes
    IF( MaskExists ) THEN
      nsize = MAXVAL( MaskPerm ) 
      CALL Info(Caller,'Applying mask of size: '//TRIM(I2S(nsize)),Level=10)
    ELSE
      nsize = nnodes
      CALL Info(Caller,'Applying mask to the whole mesh',Level=10)
    END IF 

    CoordTransform = ListGetString(Params,'Mapping Coordinate Transformation',DoCoordTransform )
    IF( DoCoordTransform .OR. MaskExists) THEN
      Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
      IF( ASSOCIATED( Var ) ) THEN
        CALL Info(Caller,'Reusing > Extruded Coordinate < variable',Level=12 )
        Values => Var % Values        
      ELSE
        NULLIFY( Values )
        ALLOCATE( Values( nsize ) )
        Values = 0.0_dp
        IF( MaskExists ) THEN
          CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values, MaskPerm)
        ELSE
          CALL VariableAdd( Mesh % Variables, Mesh, Solver,'Extruded Coordinate',1,Values)
        END IF
        Var => VariableGet( Mesh % Variables,'Extruded Coordinate')
      END IF
    ELSE IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF	      

    IF( MaskExists .OR. DoCoordTransform) THEN
      DO i=1,Mesh % NumberOfNodes
        j = i
	IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
        END IF
        Vector(1) = Mesh % Nodes % x(i)
	Vector(2) = Mesh % Nodes % y(i)
	Vector(3) = Mesh % Nodes % z(i)
	IF( DoCoordTransform ) THEN
          CALL CoordinateTransformationNodal( CoordTransform, Vector )
        END IF
        Values(j) = Vector( ActiveDirection )
      END DO
    END IF
    IF( PRESENT( ExtVar ) ) ExtVar => Var
    
    ! Check which direction is active
    !---------------------------------------------------------------------
    UpActive = PRESENT( UpNodePointer) .OR. PRESENT ( TopNodePointer ) 
    DownActive = PRESENT( DownNodePointer) .OR. PRESENT ( BotNodePointer ) 
    
    IF( PRESENT( NumberOfLayers) .OR. PRESENT( NodeLayer ) ) THEN
      UpActive = .TRUE.
      DownActive = .TRUE.
    END IF

    IF(.NOT. (UpActive .OR. DownActive ) ) THEN
      CALL Warn(Caller,'Either up or down direction should be active')
      RETURN
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    IF( UpActive ) THEN
      ALLOCATE(TopPointer(nsize),UpPointer(nsize))
      DO i=1,nnodes
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE 
        END IF
        TopPointer(j) = i
        UpPointer(j) = i
      END DO
    END IF
    IF( DownActive ) THEN
      ALLOCATE(BotPointer(nsize),DownPointer(nsize))
      DO i=1,nnodes        
        j = i
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE 
        END IF
        BotPointer(j) = i
        DownPointer(j) = i
      END DO
    END IF
    
    CALL Info(Caller,'Determine up and down pointers',Level=15)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    DO elem = 1,Mesh % NumberOfBulkElements      
      
      Element => Mesh % Elements(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element
      
      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)
      
      ! This is probably a copy-paste error, I comment it away for time being.   
      ! IF (.NOT. (Element % PartIndex == Parenv % Mype) ) CYCLE

      IF( MaskExists ) THEN
        IF( ANY(MaskPerm(NodeIndexes) == 0) ) CYCLE
      END IF
      
      DO i=1,n
        ii = NodeIndexes(i)
        
        Vector(1) = Nodes % x(i)
	Vector(2) = Nodes % y(i) 
        Vector(3) = Nodes % z(i)
        
 	IF( DoCoordTransform ) THEN
          CALL CoordinateTransformationNodal( CoordTransform, Vector )
        END IF

        MaxDotPro = -1.0_dp
        MinDotPro = 1.0_dp
        
        DO j=i+1,n
          jj = NodeIndexes(j)
          
	  Vector2(1) = Nodes % x(j)
          Vector2(2) = Nodes % y(j)
          Vector2(3) = Nodes % z(j)

	  IF( DoCoordTransform ) THEN
            CALL CoordinateTransformationNodal( CoordTransform, Vector2 )
          END IF
          
          ElemVector = Vector2 - Vector

          Length = SQRT(SUM(ElemVector*ElemVector))
          DotPro = SUM(ElemVector * UnitVector) / Length

          IF( DotPro > MaxDotPro ) THEN
            MaxDotPro = DotPro
            jmax = jj
          END IF
          IF( DotPro < MinDotPro ) THEN
            MinDotPro = DotPro
            jmin = jj
          END IF          
        END DO
          
        IF(MaxDotPro > 1.0_dp - Eps) THEN 
          IF( MaskExists ) THEN
            IF( UpActive ) UpPointer(MaskPerm(ii)) = jmax
            IF( DownActive ) DownPointer(MaskPerm(jmax)) = ii              
          ELSE
            IF( UpActive ) UpPointer(ii) = jmax
            IF( DownActive ) DownPointer(jmax) = ii
          END IF
        END IF
            
        IF(MinDotPro < Eps - 1.0_dp) THEN
          IF( MaskExists ) THEN
            IF( DownActive ) DownPointer(MaskPerm(ii)) = jmin
            IF( UpActive ) UpPointer(MaskPerm(jmin)) = ii
          ELSE
            IF( DownActive ) DownPointer(ii) = jmin
            IF( UpActive ) UpPointer(jmin) = ii              
          END IF
        END IF

      END DO
    END DO
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info(Caller,'determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      
      DO i=1,nnodes
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0) CYCLE
          IF( UpActive ) THEN
            j = UpPointer(MaskPerm(i))
            IF( TopPointer(MaskPerm(i)) /= TopPointer(MaskPerm(j)) ) THEN
              UpHit = UpHit + 1
              TopPointer(MaskPerm(i)) = TopPointer(MaskPerm(j))
            END IF
          END IF
          IF( DownActive ) THEN
            j = DownPointer(MaskPerm(i))
            IF( BotPointer(MaskPerm(i)) /= BotPointer(MaskPerm(j)) ) THEN
              DownHit = DownHit + 1
              BotPointer(MaskPerm(i)) = BotPointer(MaskPerm(j))
            END IF
          END IF
        ELSE
          IF( UpActive ) THEN
            j = UpPointer(i)
            IF( TopPointer(i) /= TopPointer(j) ) THEN
              UpHit = UpHit + 1
              TopPointer(i) = TopPointer( j )
            END IF
          END IF
          IF( DownActive ) THEN
            j = DownPointer(i)
            IF( BotPointer(i) /= BotPointer( j ) ) THEN
              DownHit = DownHit + 1
              BotPointer(i) = BotPointer( j )
            END IF
          END IF
        END IF
      END DO
      
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO

    ! The last round is always a check
    Rounds = Rounds - 1
    
    CALL Info(Caller,'Layered structure detected in '//TRIM(I2S(Rounds))//' cycles',Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info(Caller,'Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal(Caller,'Zero rounds implies unsuccessful operation')
    END IF

    ! Compute the number of layers. The Rounds above may in some cases
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info(Caller,'Compute number of layers',Level=15)    
      DO i=1,nsize
        IF( MaskExists ) THEN
          IF( MaskPerm(i) == 0 ) CYCLE
        END IF
        EXIT
      END DO

      j = BotPointer(1)      
      CALL Info(Caller,'Starting from node: '//TRIM(I2S(j)),Level=15)

      NumberOfLayers = 0
      DO WHILE(.TRUE.)
        jj = j 
        IF( MaskExists ) THEN
          jj = MaskPerm(j)
        END IF
        k = UpPointer(jj)
        IF( k == j ) THEN
          EXIT
        ELSE
          NumberOfLayers = NumberOfLayers + 1
          j = k
        END IF
      END DO

      IF( NumberOfLayers < Rounds ) THEN
        WRITE( Message,'(A,I0,A,I0)') 'There seems to be varying number of layers: ',&
            NumberOfLayers,' vs. ',Rounds
        CALL Warn(Caller, Message )
        NumberOfLayers = Rounds
      END IF
      CALL Info(Caller,&
          'Extruded structure layers: '//TRIM(I2S(NumberOfLayers)),Level=6)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( NodeLayer ) ) THEN
      CALL Info(Caller,'creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      IF( MaskExists ) THEN
        WHERE( MaskPerm == 0 ) Layer = 0
        
        DO i=1,nnodes
          IF( MaskPerm(i) == 0 ) CYCLE
          Rounds = 1
          j = BotPointer(MaskPerm(i))
          Layer(MaskPerm(j)) = Rounds
          DO WHILE(.TRUE.)
            k = UpPointer(MaskPerm(j))
            IF( k == j ) EXIT          
            Rounds = Rounds + 1
            j = k
            Layer(MaskPerm(j)) = Rounds
          END DO
        END DO
      ELSE        
        DO i=1,nsize
          Rounds = 1
          j = BotPointer(i)
          Layer(j) = Rounds
          DO WHILE(.TRUE.)
            k = UpPointer(j)
            IF( k == j ) EXIT          
            Rounds = Rounds + 1
            j = k
            Layer(j) = Rounds
          END DO
        END DO
      END IF
        
      NodeLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info(Caller,Message,Level=6)
      NULLIFY(Layer)
    END IF

    
    IF( PRESENT( MidNodePointer ) ) THEN
      ALLOCATE( MidPointer( nsize ) )
      MidPointer = 0 
      MidLayerExists = .FALSE.

      DO elem = Mesh % NumberOfBulkElements + 1, &       
          Mesh % NumberOfBulkElements + Mesh % NumberOfBoundaryElements  
        
        Element => Mesh % Elements(elem)
        NodeIndexes => Element % NodeIndexes
        
        DO bc_ind = 1, CurrentModel % NumberOfBCs 
          IF( Element % BoundaryInfo % Constraint == &
              CurrentModel % BCs(bc_ind) % Tag ) THEN
            IF( ListCheckPresent( CurrentModel % BCs(bc_ind) % Values,'Mid Surface') ) THEN
              MidPointer( NodeIndexes ) = NodeIndexes
              MidLayerExists = .TRUE.
            END IF
            EXIT
          END IF
        END DO
      END DO

      IF( MidLayerExists ) THEN
        CALL Info(Caller,'determine mid pointers',Level=15)       
                
        DO Rounds = 1, nsize
          DownHit = 0
          UpHit = 0
          DO i=1,nsize
            IF( MaskExists ) THEN
              IF( MaskPerm(i) == 0) CYCLE
            END IF

            ! We can only start from existing mid pointer
            IF( MidPointer(i) == 0 ) CYCLE
            IF( UpActive ) THEN
              j = UpPointer(i)
              IF( MaskExists ) THEN
                IF( MidPointer(MaskPerm(j)) == 0 ) THEN
                  UpHit = UpHit + 1
                  MidPointer(MaskPerm(j)) = MidPointer(MaskPerm(i))
                END IF
              ELSE
                IF( MidPointer(j) == 0 ) THEN
                  UpHit = UpHit + 1
                  MidPointer(j) = MidPointer(i)
                END IF
              END IF
            END IF
            IF( DownActive ) THEN
              j = DownPointer(i)
              IF( MaskExists ) THEN
                IF( MidPointer(MaskPerm(j)) == 0 ) THEN
                  DownHit = DownHit + 1
                  MidPointer(MaskPerm(j)) = MidPointer(MaskPerm(i))
                END IF           
              ELSE
                IF( MidPointer(j) == 0 ) THEN
                  DownHit = DownHit + 1
                  MidPointer(j) = MidPointer(i)
                END IF
              END IF
            END IF
          END DO
          IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
        END DO

        CALL Info(Caller,&
            'Mid layer structure detected in '//TRIM(I2S(Rounds-1))//' cycles',Level=9)
        MidNodePointer => MidPointer
      ELSE
        DEALLOCATE( MidPointer ) 
        MidNodePointer => NULL()
      END IF
    END IF

  
    ! Count the number of top and bottom nodes, for information only
    !---------------------------------------------------------------
    CALL Info(Caller,'Counting top and bottom nodes',Level=15)        
    IF( UpActive ) THEN
      TopNodes = 0
      MinTop = HUGE( MinTop ) 
      MaxTop = -HUGE( MaxTop )
      DO i=1,nnodes
        IF( MaskExists ) THEN
          j = MaskPerm(i) 
          IF( j == 0 ) CYCLE
          IF(TopPointer(j) == i) THEN
            MinTop = MIN( MinTop, Var % Values(j) )
            MaxTop = MAX( MaxTop, Var % Values(j) )
            TopNodes = TopNodes + 1
          END IF
        ELSE
          IF(TopPointer(i) == i) THEN
            MinTop = MIN( MinTop, Var % Values(i) )
            MaxTop = MAX( MaxTop, Var % Values(i) )
            TopNodes = TopNodes + 1
          END IF
        END IF
      END DO
    END IF

    IF( DownActive ) THEN
      BotNodes = 0
      MinBot = HUGE( MinBot ) 
      MaxBot = -HUGE( MaxBot )
      DO i=1,nnodes
        IF( MaskExists ) THEN
          j = MaskPerm(i)
          IF( j == 0 ) CYCLE
          IF( BotPointer(j) == i) THEN
            MinBot = MIN( MinBot, Var % Values(j))
            MaxBot = MAX( MaxBot, Var % Values(j))
            BotNodes = BotNodes + 1
          END IF
        ELSE          
          IF(BotPointer(i) == i) THEN
            MinBot = MIN( MinBot, Var % Values(i))
            MaxBot = MAX( MaxBot, Var % Values(i))
            BotNodes = BotNodes + 1
          END IF
        END IF
      END DO
    END IF


    ! Return the requested pointer structures, otherwise deallocate
    !---------------------------------------------------------------
    CALL Info(Caller,'Setting pointer structures',Level=15)        
    IF( UpActive ) THEN
      IF( PRESENT( TopNodePointer ) ) THEN
        TopNodePointer => TopPointer 
        NULLIFY( TopPointer )
      ELSE
        DEALLOCATE( TopPointer )
      END IF
      IF( PRESENT( UpNodePointer ) ) THEN
        UpNodePointer => UpPointer 
        NULLIFY( UpPointer )
      ELSE
        DEALLOCATE( UpPointer )
      END IF
    END IF
    IF( DownActive ) THEN
      IF( PRESENT( BotNodePointer ) ) THEN
        BotNodePointer => BotPointer 
        NULLIFY( BotPointer ) 
      ELSE
        DEALLOCATE( BotPointer )
      END IF
      IF( PRESENT( DownNodePointer ) ) THEN
        DownNodePointer => DownPointer 
        NULLIFY( DownPointer ) 
      ELSE
        DEALLOCATE( DownPointer )
      END IF
    END IF

    !---------------------------------------------------------------
    at1 = CPUTime()  
    WRITE(Message,* ) 'Top and bottom pointer init time: ',at1-at0
    CALL Info(Caller,Message,Level=6)
    CALL Info(Caller,&
        'Top and bottom pointer init rounds: '//TRIM(I2S(Rounds)),Level=5)
    IF( UpActive ) THEN
      CALL Info(Caller,'Number of nodes at the top: '//TRIM(I2S(TopNodes)),Level=6)
    END IF
    IF( DownActive ) THEN
      CALL Info(Caller,'Number of nodes at the bottom: '//TRIM(I2S(BotNodes)),Level=6)
    END IF
    

  CONTAINS
    
    
    !---------------------------------------------------------------
    SUBROUTINE CoordinateTransformationNodal( CoordTransform, R )
      CHARACTER(LEN=MAX_NAME_LEN) :: CoordTransform
      REAL(KIND=dp) :: R(3)
      !---------------------------------------------------------------
      REAL(KIND=dp) :: Rtmp(3)
      REAL(KIND=dp), SAVE :: Coeff 
      LOGICAL, SAVE :: Visited = .FALSE.
      

      IF( .NOT. Visited ) THEN
        IF( ListGetLogical( Params,'Angles in Degrees') ) THEN
          Coeff = 180.0_dp / PI
        ELSE
          Coeff = 1.0_dp
        END IF
        Visited = .TRUE.
      END IF
      
      SELECT CASE ( CoordTransform )
        
      CASE('cartesian to cylindrical')
        Rtmp(1) = SQRT( R(1)**2 + R(2)**2)
        Rtmp(2) = Coeff * ATAN2( R(2), R(1)  ) 
        Rtmp(3) = R(3) 
        
      CASE('cylindrical to cartesian')
        Rtmp(1) = COS( R(2) / Coeff ) * R(1)
        Rtmp(2) = SIN( R(2) / Coeff ) * R(1)
        Rtmp(3) = R(3)
        
      CASE DEFAULT
        CALL Fatal('CoordinateTransformationNodal','Unknown transformation: '//TRIM(CoordTransform) )
        
      END SELECT
      
      R = Rtmp

    END SUBROUTINE CoordinateTransformationNodal
   

  END SUBROUTINE DetectExtrudedStructure
 !---------------------------------------------------------------



!--------------------------------------------------------------------------
!> This subroutine finds the structure of an extruded mesh for elements.
!> Otherwise very similar as the DetectExtrudedStructure for nodes.
!> Mesh faces may need to be created in order to determine the up and down
!> pointers.
!-----------------------------------------------------------------------------
  SUBROUTINE DetectExtrudedElements( Mesh, Solver, ExtVar, &
      TopElemPointer, BotElemPointer, UpElemPointer, DownElemPointer, &
      NumberOfLayers, ElemLayer )
    
    USE CoordinateSystems
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Solver_t), POINTER :: Solver
    TYPE(Variable_t), POINTER, OPTIONAL :: ExtVar
    INTEGER, POINTER, OPTIONAL :: TopElemPointer(:), BotElemPointer(:), &
        UpElemPointer(:), DownElemPointer(:)
    INTEGER, POINTER, OPTIONAL :: ElemLayer(:)
    INTEGER, OPTIONAL :: NumberOfLayers
!-----------------------------------------------------------------------------
    REAL(KIND=dp) :: Direction(3)
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Variable_t), POINTER :: Var
    REAL(KIND=dp) :: Tolerance
    TYPE(Element_t), POINTER :: Element, Parent
    TYPE(Nodes_t) :: Nodes
    INTEGER :: i,j,k,n,ii,jj,dim, nsize, elem, TopNodes, BotNodes, Rounds, ActiveDirection, &
	UpHit, DownHit, bc_ind
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: UpActive, DownActive, GotIt, Found
    LOGICAL, POINTER :: TopFlag(:), BotFlag(:)
    REAL(KIND=dp) :: at0, at1
    REAL(KIND=dp) :: FaceCenter(3),FaceDx(3),Height(2),Eps, MinTop, MaxTop, MinBot, MaxBot, Diam
    REAL(KIND=dp), POINTER :: Values(:)
    INTEGER, POINTER :: TopPointer(:), BotPointer(:), UpPointer(:), DownPointer(:),Layer(:),MidPointer(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: VarName
    INTEGER :: TestCounter(3),ElementIndex(2)
    CHARACTER(*), PARAMETER :: Caller='DetectExtrudedElements '
        
    CALL Info(Caller,'Determining extruded element structure',Level=6)
    at0 = CPUTime()

    DIM = Mesh % MeshDim

    IF( DIM /= 3 ) THEN
      CALL Fatal(Caller,'Only implemented for 3D cases: '//TRIM(I2S(dim)))
    END IF

    IF( .NOT. ASSOCIATED( Mesh % Faces ) ) THEN
      CALL FindMeshFaces3D( Mesh )
    END IF

    
    Params => Solver % Values
    TestCounter = 0
    
    ActiveDirection = ListGetInteger(Params,'Active Coordinate')
    IF( ActiveDirection < 1 .OR. ActiveDirection > 3 ) THEN
      CALL Fatal(Caller,'Invalid value for Active Coordinate')
    END IF  

    ! Set the dot product tolerance
    !-----------------------------------------------------------------
    Eps = ListGetConstReal( Params,'Dot Product Tolerance',GotIt)
    IF(.NOT. GotIt) Eps = 1.0d-1

    nsize = Mesh % NumberOfBulkElements
    CALL Info(Caller,'Detecting extrusion in the mesh using coordinate: '&
        //TRIM(I2S(ActiveDirection)),Level=8)

    IF( ActiveDirection == 1 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 1')
    ELSE IF( ActiveDirection == 2 ) THEN
      Var => VariableGet( Mesh % Variables,'Coordinate 2')
    ELSE 
      Var => VariableGet( Mesh % Variables,'Coordinate 3')
    END IF	      

    IF( PRESENT( ExtVar ) ) ExtVar => Var

    ! Check which direction is active
    !---------------------------------------------------------------------
    UpActive = PRESENT( UpElemPointer) .OR. PRESENT ( TopElemPointer ) 
    DownActive = PRESENT( DownElemPointer) .OR. PRESENT ( BotElemPointer ) 

    IF( PRESENT( NumberOfLayers) .OR. PRESENT( ElemLayer ) ) THEN
      UpActive = .TRUE.
      DownActive = .TRUE.
    END IF

    IF(.NOT. (UpActive .OR. DownActive ) ) THEN
      CALL Warn(Caller,'Either up or down direction should be active')
      RETURN
    END IF

    ! Allocate pointers to top and bottom, and temporary pointers up and down
    !------------------------------------------------------------------------
    IF( UpActive ) THEN
      ALLOCATE(TopPointer(nsize),UpPointer(nsize))
      DO i=1,nsize
        TopPointer(i) = i
        UpPointer(i) = i
      END DO
    END IF
    IF( DownActive ) THEN
      ALLOCATE(BotPointer(nsize),DownPointer(nsize))
      DO i=1,nsize
        BotPointer(i) = i
        DownPointer(i) = i
      END DO
    END IF

    CALL Info(Caller,'determine up and down pointers',Level=15)

    ! Determine the up and down pointers using dot product as criterion
    !-----------------------------------------------------------------
    n = Mesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n),Nodes % z(n) )
    
    DO elem = 1,Mesh % NumberOfFaces 

      Element => Mesh % Faces(elem)
      NodeIndexes => Element % NodeIndexes
      CurrentModel % CurrentElement => Element

      n = Element % TYPE % NumberOfNodes
      Nodes % x(1:n) = Mesh % Nodes % x(NodeIndexes)
      Nodes % y(1:n) = Mesh % Nodes % y(NodeIndexes)
      Nodes % z(1:n) = Mesh % Nodes % z(NodeIndexes)

      IF( .NOT. ASSOCIATED( Element % BoundaryInfo ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Left ) ) CYCLE
      IF( .NOT. ASSOCIATED( Element % BoundaryInfo % Right ) ) CYCLE
      
      FaceCenter(1) = SUM( Nodes % x(1:n) ) / n
      FaceCenter(2) = SUM( Nodes % y(1:n) ) / n
      FaceCenter(3) = SUM( Nodes % z(1:n) ) / n

      FaceDx(1) = SUM( ABS( Nodes % x(1:n) - FaceCenter(1) ) ) 
      FaceDx(2) = SUM( ABS( Nodes % y(1:n) - FaceCenter(2) ) ) 
      FaceDx(3) = SUM( ABS( Nodes % z(1:n) - FaceCenter(3) ) ) 
      
      Diam = SQRT( SUM( FaceDx**2 ) )

      ! This is not a face that separates extruded elements
      IF( FaceDx(ActiveDirection) > Eps * Diam ) CYCLE      

      TestCounter(1) = TestCounter(1) + 1      
      
      DO k = 1, 2
        IF( k == 1 ) THEN
          Parent => Element % BoundaryInfo % Left
        ELSE
          Parent => Element % BoundaryInfo % Right
        END IF
        IF( .NOT. ASSOCIATED( Parent ) ) CYCLE
               
        n = Parent % TYPE % NumberOfNodes
        NodeIndexes => Parent % NodeIndexes        

        ElementIndex(k) = Parent % ElementIndex
        Height(k) = SUM( Var % Values(NodeIndexes) ) / n
      END DO      

      IF( Height(1) > Height(2) ) THEN
        IF( UpActive ) UpPointer(ElementIndex(2)) = ElementIndex(1)
        IF( DownActive ) DownPointer(ElementIndex(1)) = ElementIndex(2)
      ELSE
        IF( UpActive ) UpPointer(ElementIndex(1)) = ElementIndex(2)
        IF( DownActive ) DownPointer(ElementIndex(2)) = ElementIndex(1)
      END IF
    END DO  
        
    DEALLOCATE( Nodes % x, Nodes % y,Nodes % z )

    
    ! Pointer to top and bottom are found recursively using up and down
    !------------------------------------------------------------------
    CALL Info(Caller,'determine top and bottom pointers',Level=9)

    DO Rounds = 1, nsize
      DownHit = 0
      UpHit = 0
      DO i=1,nsize
        IF( UpActive ) THEN
          j = UpPointer(i)
          IF( TopPointer(i) /= TopPointer( j ) ) THEN
            UpHit = UpHit + 1
            TopPointer(i) = TopPointer( j )
          END IF
        END IF
        IF( DownActive ) THEN
          j = DownPointer(i)
          IF( BotPointer(i) /= BotPointer( j ) ) THEN
	    DownHit = DownHit + 1
            BotPointer(i) = BotPointer( j )
          END IF
        END IF
      END DO
      CALL Info(Caller,'Hits in determining structure: '//TRIM(I2S(UpHit+DownHit)),Level=10)
      IF( UpHit == 0 .AND. DownHit == 0 ) EXIT
    END DO
    ! The last round is always a check
    Rounds = Rounds - 1


    WRITE( Message,'(A,I0,A)') 'Layered elements detected in ',Rounds,' cycles'
    CALL Info(Caller,Message,Level=9)
    IF( Rounds == 0 ) THEN
      CALL Info(Caller,'Try to increase value for > Dot Product Tolerance < ')
      CALL Fatal(Caller,'Zero rounds implies unsuccessful operation')
    END IF


    ! Compute the number of layers. The Rounds above may in some cases 
    ! be too small. Here just one layer is used to determine the number
    ! of layers to save some time.
    !------------------------------------------------------------------
    IF( PRESENT( NumberOfLayers ) ) THEN
      CALL Info(Caller,'Compute number of layers',Level=15)    

      ! We start from any bottom row entry
      j = BotPointer(1)
      
      NumberOfLayers = 0
      DO WHILE(.TRUE.)
        k = UpPointer(j)

        IF( k == j ) THEN
          EXIT
        ELSE
          NumberOfLayers = NumberOfLayers + 1
          j = k
        END IF
      END DO      

      IF( NumberOfLayers < Rounds ) THEN
        WRITE( Message,'(A,I0,A,I0)') 'There seems to be varying number of layers: ',&
            NumberOfLayers,' vs. ',Rounds
        CALL Warn(Caller, Message )
        NumberOfLayers = Rounds
      END IF
      CALL Info(Caller,'Extruded structure layers: '//TRIM(I2S(NumberOfLayers)),Level=6)
    END IF

    
    ! Create layer index if requested
    !------------------------------------------------------------------
    IF( PRESENT( ElemLayer ) ) THEN
      CALL Info(Caller,'creating layer index',Level=9)        

      NULLIFY(Layer)
      ALLOCATE( Layer(nsize) )
      Layer = 1
      
      DO i=1,nsize
        Rounds = 1
        j = BotPointer(i)
        Layer(j) = Rounds
        DO WHILE(.TRUE.)
          k = UpPointer(j)
          IF( k == j ) EXIT          
          Rounds = Rounds + 1
          j = k
          Layer(j) = Rounds
        END DO
      END DO
      
      ElemLayer => Layer
      WRITE(Message,'(A,I0,A,I0,A)') 'Layer range: [',MINVAL(Layer),',',MAXVAL(Layer),']'
      CALL Info(Caller,Message,Level=6)
      NULLIFY(Layer)
    END IF

  
    ! Count the number of top and bottom elements, for information only
    !---------------------------------------------------------------
    CALL Info(Caller,'Counting top and bottom elements',Level=15)        
    IF( UpActive ) THEN
      TopNodes = 0
      MinTop = HUGE( MinTop ) 
      MaxTop = -HUGE( MaxTop )
      DO i=1,nsize
        IF(TopPointer(i) == i) THEN
          MinTop = MIN( MinTop, Var % Values(i) )
          MaxTop = MAX( MaxTop, Var % Values(i) )
          TopNodes = TopNodes + 1
        END IF
      END DO
      CALL Info(Caller,'Number of top elements: '//TRIM(I2S(TopNodes)),Level=9)
    END IF

    IF( DownActive ) THEN
      BotNodes = 0
      MinBot = HUGE( MinBot ) 
      MaxBot = -HUGE( MaxBot )
      DO i=1,nsize
        IF(BotPointer(i) == i) THEN
          MinBot = MIN( MinBot, Var % Values(i))
          MaxBot = MAX( MaxBot, Var % Values(i))
          BotNodes = BotNodes + 1
        END IF
      END DO
    END IF


    ! Return the requested pointer structures, otherwise deallocate
    !---------------------------------------------------------------
    CALL Info(Caller,'Setting pointer structures',Level=15)        
    IF( UpActive ) THEN
      IF( PRESENT( TopElemPointer ) ) THEN
        TopElemPointer => TopPointer 
        NULLIFY( TopPointer )
      ELSE
        DEALLOCATE( TopPointer )
      END IF
      IF( PRESENT( UpElemPointer ) ) THEN
        UpElemPointer => UpPointer 
        NULLIFY( UpPointer )
      ELSE
        DEALLOCATE( UpPointer )
      END IF
    END IF
    IF( DownActive ) THEN
      IF( PRESENT( BotElemPointer ) ) THEN
        BotElemPointer => BotPointer 
        NULLIFY( BotPointer ) 
      ELSE
        DEALLOCATE( BotPointer )
      END IF
      IF( PRESENT( DownElemPointer ) ) THEN
        DownElemPointer => DownPointer 
        NULLIFY( DownPointer ) 
      ELSE
        DEALLOCATE( DownPointer )
      END IF
    END IF

    !---------------------------------------------------------------
    at1 = CPUTime()  
    WRITE(Message,'(A,ES12.3)') 'Top and bottom pointer init time: ',at1-at0
    CALL Info(Caller,Message,Level=6)

    CALL Info(Caller,'Top and bottom pointer init rounds: '//TRIM(I2S(Rounds)),Level=8)

    IF( UpActive ) THEN
      CALL Info(Caller,'Number of elements at the top: '//TRIM(I2S(TopNodes)),Level=8)
    END IF
    IF( DownActive ) THEN
      CALL Info(Caller,'Number of elements at the bottom: '//TRIM(I2S(BotNodes)),Level=8)
    END IF
   

  END SUBROUTINE DetectExtrudedElements
 !---------------------------------------------------------------

END MODULE DetectExtruded
