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
!>  Utilities for mortar and other routines involving boundary meshes.
!------------------------------------------------------------------------------

MODULE MortarProjector

    USE ElementUtils
    USE ElementDescription
    USE MeshUtils
    USE Interpolation
    USE ParallelUtils
    USE Types
    IMPLICIT NONE

CONTAINS


!------------------------------------------------------------------------------
!> Given two interface meshes check the angle between them using the normal
!> vectors of the first element. Also check that all other elements are
!> aligned with the first one. Only then is it possible to determine the angle.
!------------------------------------------------------------------------------
  SUBROUTINE CheckInterfaceMeshAngle(BMesh1, BMesh2, Angles, GotAngles) 
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    REAL(KIND=dp) :: Angles(3)
    LOGICAL :: GotAngles
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n
    REAL(KIND=dp) :: Normal(3), Normal1(3), Normal2(3), Dot1Min, Dot2Min, Alpha
    LOGICAL :: ConstantNormals
    CHARACTER(*), PARAMETER :: Caller='CheckInterfaceMeshAngle'

    ! Currently check of the normal direction is not enforced since at this stage 
    ! CurrentModel % Nodes may not exist!
    ! This means that there may be a 180 error in the directions. 
    ! Therefore an angle smaller than 180 is always chosen.
    !-----------------------------------------------------------------------------
    N = MAX( BMesh1 % MaxElementNodes, BMesh2 % MaxElementNodes )
    ALLOCATE(ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
 
    DO k=1,2
      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      ! we use the Dot2Min and Normal2 temporarily also for first mesh, with k=1
      !-------------------------------------------------------------------------
      DO i=1, PMesh % NumberOfBoundaryElements
        Element => PMesh % Elements(i)
        
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = PMesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = PMesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = PMesh % Nodes % z(NodeIndexes(1:n))           
        
        Normal = NormalVector( Element, ElementNodes, Check = .FALSE. ) 

        ! we use the Dot2Min and Normal2 temporarily also for first mesh, with k=1
        !-------------------------------------------------------------------------       
        IF( i == 1 ) THEN
          Normal2 = Normal
          Dot2Min = 1.0_dp
        ELSE
          Dot2min = MIN( Dot2Min, SUM( Normal * Normal2 ) )
        END IF
      END DO

      IF( k == 1 ) THEN
        Normal1 = Normal2 
        Dot1Min = Dot2Min
      END IF
    END DO

    ConstantNormals = ( 1 - Dot1Min < 1.0d-6 ) .AND. ( 1 - Dot2Min < 1.0d-6 )     
    IF( ConstantNormals ) THEN
      WRITE(Message,'(A,3ES12.3)') 'Master normal: ',Normal1
      CALL Info(Caller,Message,Level=8)    
      
      WRITE(Message,'(A,3ES12.3)') 'Initial Target normal: ',Normal2
      CALL Info(Caller,Message,Level=8)    
            
      ! The full angle between the two normals
      Alpha = ACOS( SUM( Normal1 * Normal2 ) ) * 180.0_dp / PI
      WRITE(Message,'(A,ES12.3)') &
          'Suggested angle between two normals in degs (+/- 180): ',Alpha 
      CALL Info(Caller,Message,Level=8)
    ELSE
      CALL Warn(Caller,'Could not suggest rotation angle')
    END IF


    GotAngles = .FALSE.
    Angles = 0.0_dp
    IF( .NOT. ConstantNormals ) THEN
      CALL Warn(Caller,'Normals are not constant, cannot test for rotation!')
    ELSE IF( Alpha > EPSILON( Alpha ) ) THEN
      ! Rotation should be performed 
      DO i=1,3
        IF( ABS ( Normal1(i) - Normal2(i) ) < EPSILON( Alpha ) ) THEN
          GotAngles = .TRUE.            
          WRITE(Message,'(A,I0,A,ES12.3)') &
              'Rotation around axis ',i,' in degs ',Alpha 
          CALL Info(Caller,Message,Level=8)
          Angles(i) = Alpha
          EXIT
        END IF
      END DO
      IF(.NOT. GotAngles ) THEN
        CALL Warn(Caller,'could not define rotation axis, improve algorithm!')
      END IF
    END IF

    DEALLOCATE(ElementNodes % x, ElementNodes % y, ElementNodes % z )
    
  END SUBROUTINE CheckInterfaceMeshAngle
!------------------------------------------------------------------------------

  
  !------------------------------------------------------------------------------
  !> Find axial, radial or rotational mortar boundary pairs.
  !------------------------------------------------------------------------------
  SUBROUTINE DetectMortarPairs( Model, Mesh, Tol, BCMode, SameCoordinate )
    !------------------------------------------------------------------------------    
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
    REAL(KIND=dp) :: Tol
    INTEGER :: BcMode
    LOGICAL :: SameCoordinate
    !------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n,MinBC,MaxBC,BC,ElemCode
    TYPE(Element_t), POINTER :: Element, Parent, Left, Right, Elements(:)
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Found 
    LOGICAL, ALLOCATABLE :: BCSet(:), BCPos(:), BCNeg(:), BCNot(:)
    INTEGER, ALLOCATABLE :: BCCount(:)
    REAL(KIND=dp) :: x,y,z,f
    REAL(KIND=dp), ALLOCATABLE :: BCVal(:)
    CHARACTER(LEN=MAX_NAME_LEN) :: str
    LOGICAL :: Debug = .FALSE., Hit
    CHARACTER(*), PARAMETER :: Caller='DetectMortarPairs'
   
    ! The code can detect pairs to be glued in different coordinate systems
    SELECT CASE( BCMode )
    CASE( 1 )
      str = 'x-coordinate'
    CASE( 2 )
      str = 'y-coordinate'
    CASE( 3 )
      str = 'z-coordinate'
    CASE( 4 ) 
      str = 'radius'
    CASE( 5 )
      str = 'angle'
    CASE DEFAULT
      CALL Fatal(Caller,'Invalid BCMode: '//TRIM(I2S(BCMode)))
    END SELECT

    CALL Info(Caller,'Trying to find pairs in: '//TRIM(str),Level=6)

    IF(.NOT. ASSOCIATED( Mesh ) ) THEN
      CALL Fatal(Caller,'Mesh not associated!')
    END IF

    IF( ParEnv % PEs > 1 ) THEN
      CALL Warn(Caller,'Not implemented in parallel yet, be careful!')
    END IF
      
    
    ! Interface meshes consist of boundary elements only    
    Elements => Mesh % Elements( Mesh % NumberOfBulkElements+1: )

    ! Find out the min and max constraint 
    MinBC = HUGE( MinBC )
    MaxBC = 0
    DO i=1, Mesh % NumberOfBoundaryElements
      Element => Elements(i)
      ElemCode = Element % Type % ElementCode 
      IF (ElemCode<=200) CYCLE

      BC = Element % BoundaryInfo % Constraint
      MinBC = MIN( MinBC, BC )
      MaxBC = MAX( MaxBC, BC )
    END DO

    CALL Info(Caller,'Minimum Constraint index: '//TRIM(I2S(MinBC)),Level=8)
    CALL Info(Caller,'Maximum Constraint index: '//TRIM(I2S(MaxBC)),Level=8)    
    IF( MaxBC - MinBC < 1 ) THEN
      CALL Warn(Caller,'Needs at least two different BC indexes to create mortar pair!')
      RETURN
    END IF

    ALLOCATE( BCVal( MinBC:MaxBC ) )
    ALLOCATE( BCSet( MinBC:MaxBC ) )
    ALLOCATE( BCNot( MinBC:MaxBC ) )
    ALLOCATE( BCPos( MinBC:MaxBC ) )
    ALLOCATE( BCNeg( MinBC:MaxBC ) )
    ALLOCATE( BCCount( MinBC:MaxBC ) )

    BCVal = 0.0_dp
    BCSet = .FALSE.
    BCNot = .FALSE.
    BCPos = .FALSE.
    BCNeg = .FALSE.
    BCCount = 0
    

    DO i=1, Mesh % NumberOfBoundaryElements
      Element => Elements(i)
      ElemCode = Element % Type % ElementCode 
      IF (ElemCode<=200) CYCLE

      BC = Element % BoundaryInfo % Constraint

      ! This boundary is already deemed not to be a good candidate
      IF( BCNot( BC ) ) CYCLE

      n = Element % Type % NumberOfNodes

      DO j=1,n
        k = Element % NodeIndexes(j)
        x = Mesh % Nodes % x(k)
        y = Mesh % Nodes % y(k)
        z = Mesh % Nodes % z(k)

        ! Here f is a measure: x, y, z, radius, or angle 
        SELECT CASE( BCMode )
        CASE( 1 )
          f = x
        CASE( 2 )
          f = y
        CASE( 3 )
          f = z
        CASE( 4 ) 
          f = SQRT( x**2 + y**2 )
        CASE( 5 )
          f = ATAN2( y, x )
        END SELECT

        ! If the BC is not set then let the first be the one to compare against
        IF( .NOT. BCSet( BC ) ) THEN
          BCVal( BC ) = f
          BCSet( BC ) = .TRUE.
          IF( Debug ) PRINT *,'Compareing BC '//TRIM(I2S(BC))//' against:',f
        ELSE
          ! In consecutive rounds check that the level is consistent
          IF( ABS( f - BCVal(BC) ) > Tol ) THEN
            IF( Debug ) PRINT *,'Failing BC '//TRIM(I2S(BC))//' with:',f-BCVal(BC)
            BCNot( BC ) = .TRUE.
            EXIT
          END IF
        END IF
      END DO

      IF( BCNot( BC ) ) CYCLE

      Parent => Element % BoundaryInfo % Left
      IF( .NOT. ASSOCIATED( Parent ) ) THEN
        Parent => Element % BoundaryInfo % Right
      ELSE
        ! If there are two parents this is an internal BC
        IF( ASSOCIATED( Element % BoundaryInfo % Right ) ) THEN
          IF( Debug ) PRINT *,'Failing internal BC:',BC
          BCNot( BC ) = .TRUE.
          CYCLE
        END IF
      END IF

      ! To define whether the boundar is on positive or negative side of the master element
      ! study the center point of the master element
      n = Parent % TYPE % NumberOfNodes
      x = SUM( Mesh % Nodes % x( Parent % NodeIndexes) ) / n
      y = SUM( Mesh % Nodes % y( Parent % NodeIndexes) ) / n
      z = SUM( Mesh % Nodes % z( Parent % NodeIndexes) ) / n


      SELECT CASE( BCMode )
      CASE( 1 )
        f = x
      CASE( 2 )
        f = y
      CASE( 3 )
        f = z
      CASE( 4 ) 
        f = SQRT( x**2 + y**2 )
      CASE( 5 )
        f = ATAN2( y, x )
      END SELECT

      ! If the parent element is on alternating sides then this cannot be a proper boundary
      IF( f > BCVal( BC ) ) THEN
        IF( BCNeg( BC ) ) THEN
          IF( Debug ) PRINT *,'Failing inconsistent negative BC:',BC
          BCNot( BC ) = .TRUE.
          BCNeg( BC ) = .FALSE.
          CYCLE
        END IF
        BCPos( BC ) = .TRUE.
      ELSE
        IF( BCPos( BC ) ) THEN
          IF( Debug ) PRINT *,'Failing inconsistent positive BC:',BC
          BCNot( BC ) = .TRUE.
          BCPos( BC ) = .FALSE.
          CYCLE
        END IF
        BCNeg( BC ) = .TRUE.
      END IF
    END DO ! Number of boundary elements

    IF( BCMode == 5 ) THEN
      BCVal = 180.0_dp * BCVal / PI
    END IF
    
    j = COUNT( BCPos )
    IF( Debug ) THEN
      IF( j > 0 ) THEN
        IF( Debug ) PRINT *,'Positive constant levels: ',j
        DO i=MinBC,MaxBC
          IF( BCPos(i) ) PRINT *,'BC:',i,BCVal(i)
        END DO
      END IF
    END IF
      
    k = COUNT( BCNeg )
    IF( Debug ) THEN
      IF( k > 0 ) THEN
        PRINT *,'Negative constant levels: ',k
        DO i=MinBC,MaxBC
          IF( BCNeg(i) ) PRINT *,'BC:',i,BCVal(i)
        END DO
      END IF
    END IF
      
    IF( j * k == 0 ) THEN
      PRINT *,'Not enough candidate sides found'
      RETURN
    END IF

    IF( SameCoordinate ) THEN
      DO i=MinBC,MaxBC
        Hit = .FALSE.
        IF( BCPos(i) ) THEN
          DO j=MinBC,MaxBC
            IF ( BCNeg(j) ) THEN
              IF( ABS( BCVal(i) - BCVal(j)) < Tol ) THEN
                Hit = .TRUE.
                EXIT
              END IF
            END IF
          END DO
          IF( .NOT. Hit ) THEN
            BCPos(i) = .FALSE.
            IF( Debug ) PRINT *,'Removing potential positive hit:',i
          END IF
        END IF
        IF( BCNeg(i) ) THEN
          Hit = .FALSE.
          DO j=MinBC,MaxBC
            IF ( BCPos(j) ) THEN
              IF( ABS( BCVal(i) - BCVal(j)) < Tol ) THEN
                Hit = .TRUE.
                EXIT
              END IF
            END IF
          END DO
          IF( .NOT. Hit ) THEN
            BCNeg(i) = .FALSE.
            IF( Debug ) PRINT *,'Removing potential negative hit:',i
          END IF
        END IF
      END DO

      IF( .NOT. ANY( BCPos ) ) THEN 
        PRINT *,'No possible pairs found at same location'
        RETURN
      END IF
    END IF


    k = 0
    DO i=MinBC,MaxBC
      IF( BCPos(i) ) THEN
        Hit = .FALSE.
        DO j=MinBC,i-1
          IF( BCPos(j) ) THEN
            IF( ABS( BCVal(i) - BCVal(j) ) < Tol ) THEN
              Hit = .TRUE.
              EXIT
            END IF
          END IF
        END DO
        IF(Hit ) THEN
          BCCount(i) = BCCount(j)
        ELSE
          k = k + 1
          BCCount(i) = k
        END IF
      END IF
    END DO
    PRINT *,'Found number of positive levels:',k


    k = 0
    DO i=MinBC,MaxBC
      IF( BCNeg(i) ) THEN
        Hit = .FALSE.
        DO j=MinBC,i-1
          IF( BCNeg(j) ) THEN
            IF( ABS( BCVal(i) - BCVal(j) ) < Tol ) THEN
              Hit = .TRUE.
              EXIT
            END IF
          END IF
        END DO
        IF(Hit ) THEN
          BCCount(i) = BCCount(j)
        ELSE
          k = k + 1
          BCCount(i) = -k
        END IF
      END IF
    END DO
    PRINT *,'Found number of negative levels:',k

    PRINT *,'Slave BCs: '
    DO i=MinBC,MaxBC
      IF( BCPos(i) ) PRINT *,'BC:',i,BCVal(i)
    END DO
    PRINT *,'Master BCs: '
    DO i=MinBC,MaxBC
      IF( BCNeg(i) ) PRINT *,'BC:',i,BCVal(i)
    END DO
    
  END SUBROUTINE DetectMortarPairs

  

!------------------------------------------------------------------------------
!> Create master and slave mesh for the interface in order to at a later 
!> stage create projector matrix to implement periodicity or mortar elements.
!> The idea is to use a reduced set of elements and thereby speed up the 
!> mapping process. Also this gives more flexibility in transformation
!> operations since the nodes may be ereased after use. 
!------------------------------------------------------------------------------
  SUBROUTINE CreateInterfaceMeshes( Model, Mesh, This, Trgt, BMesh1, BMesh2, &
      Success ) 
!------------------------------------------------------------------------------    
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    TYPE(Mesh_t), TARGET :: Mesh
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL :: Success
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,m,n,n1,n2,k1,k2,ind,Constraint,DIM,ii,jj,kk
    TYPE(Element_t), POINTER :: Element, Left, Right, Elements(:)
    LOGICAL :: ThisActive, TargetActive
    INTEGER, POINTER :: NodeIndexes(:), Perm1(:), Perm2(:), PPerm(:)
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    LOGICAL :: OnTheFlyBC, CheckForHalo, NarrowHalo, NoHalo, SplitQuadratic, Found

    TYPE(Element_t), POINTER :: Parent,q
    INTEGER :: en, in, HaloCount, ActiveCount, ElemCode, nSplit
    INTEGER :: SplitMap(4), SplitSizes(5)
    LOGICAL, ALLOCATABLE :: ActiveNode(:)

    LOGICAL :: TagNormalFlip, Turn
    TYPE(Nodes_t) :: ElementNodes
    REAL(KIND=dp) :: Normal(3)
    CHARACTER(*), PARAMETER :: Caller='CreateInterfaceMeshes'
   
    CALL Info(Caller,'Making a list of elements at interface',Level=9)

   
    IF ( This <= 0 .OR. Trgt <= 0 ) THEN
      CALL Fatal(Caller,'Invalid target boundaries')
    END IF

    ! Interface meshes consist of boundary elements only    
    Elements => Mesh % Elements( Mesh % NumberOfBulkElements+1: )

    ! We need direction of initial normal if we have a "normal projector"
    TagNormalFlip = ListGetLogical( Model % BCs(This) % Values,'Normal Projector',Found )
    IF( TagNormalFlip ) THEN
      CALL Info(Caller,'Storing initial information on normal directions',Level=12)
      n = Mesh % MaxElementNodes
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    END IF
    
    
    SplitQuadratic = ListGetLogical( Model % Simulation,'Mortar BCs Split Quadratic',Found ) 
    IF( Mesh % NumberOfFaces > 0 .OR. Mesh % NumberOfEdges > 0 ) THEN
      SplitQuadratic = .FALSE.
    END IF
    IF( SplitQuadratic ) CALL Info(Caller,&
        'Quadratic elements will be split',Level=7)



    ! If the target is larger than number of BCs given then
    ! it has probably been created on-the-fly from a discontinuous boundary.
    OnTheFlyBC = ( Trgt > Model % NumberOfBCs )

    ! In parallel we may have some excess halo elements. 
    ! To eliminate them mark the nodes that are associated to elements truly owned. 
    NarrowHalo = .FALSE.
    NoHalo = .FALSE.

    IF( ParEnv % PEs > 1 ) THEN
      ! Account for halo elements that share some nodes for the master boundary
      NarrowHalo = ListGetLogical(Model % Solver % Values,'Projector Narrow Halo',Found)

      ! Do not allow for any halo elements for the master boundary
      IF( .NOT. Found ) THEN
        NoHalo = ListGetLogical(Model % Solver % Values,'Projector No Halo',Found)
      END IF
      
      IF(.NOT. Found ) THEN
        IF( ListGetLogical(Model % Solver % Values, 'Partition Local Constraints',Found) ) THEN
          NarrowHalo = .TRUE.
        ELSE
          NoHalo = .TRUE.
        END IF
      END IF
    END IF

    ! This is just temporarily set to false always until the logic has been tested. 
    CheckForHalo = NarrowHalo .OR. NoHalo

    IF( CheckForHalo ) THEN
      CALL Info(Caller,'Checking for halo elements',Level=15)
      ALLOCATE( ActiveNode( Mesh % NumberOfNodes ) )
      HaloCount = 0
      ActiveNode = .FALSE.
      DO i=1, Mesh % NumberOfBoundaryElements
        Element => Elements(i)
        IF (Element % TYPE % ElementCode<=200) CYCLE

        Left => Element % BoundaryInfo % Left 
        IF( ASSOCIATED( Left ) ) THEN
          IF( Left % PartIndex == ParEnv % MyPe ) THEN
            ActiveNode( Left % NodeIndexes ) = .TRUE.
          ELSE
            HaloCount = HaloCount + 1
          END IF
        END IF

        Right => Element % BoundaryInfo % Right
        IF( ASSOCIATED( Right ) ) THEN
          IF( Right % PartIndex == ParEnv % MyPe ) THEN
            ActiveNode( Right % NodeIndexes ) = .TRUE.
          ELSE
            HaloCount = HaloCount + 1 
          END IF
        END IF
      END DO

      ! No halo element found on the boundary so no need to check them later
      IF( HaloCount == 0 ) THEN
        CALL Info(Caller,'Found no halo elements to eliminate',Level=15)
        DEALLOCATE( ActiveNode ) 
        CheckForHalo = .FALSE.
      ELSE
        CALL Info(Caller,'Number of halo elements to eliminate: '&
            //TRIM(I2S(HaloCount)),Level=12)
      END IF
    END IF


!   Search elements in this boundary and its periodic
!   counterpart:
!   --------------------------------------------------
    n1 = 0
    n2 = 0
    HaloCount = 0
    DO i=1, Mesh % NumberOfBoundaryElements
      Element => Elements(i)
      ElemCode = Element % Type % ElementCode 
      IF (ElemCode<=200) CYCLE

      nSplit = 1
      IF( SplitQuadratic ) THEN
        IF( ElemCode == 306 .OR. ElemCode == 409 ) THEN
          nSplit = 4
        ELSE IF( ElemCode == 408 ) THEN
          nSplit = 5
        END IF
      END IF

      Constraint = Element % BoundaryInfo % Constraint
      IF( Model % BCs(This) % Tag == Constraint ) THEN
        IF( CheckForHalo ) THEN
          IF( NarrowHalo ) THEN
            IF( ANY(ActiveNode(Element % NodeIndexes) ) ) THEN
              n1 = n1 + nSplit
            ELSE
              HaloCount = HaloCount + 1
            END IF
          ELSE IF( NoHalo ) THEN
            ThisActive = .FALSE.
            Left => Element % BoundaryInfo % Left 
            IF( ASSOCIATED( Left ) ) THEN
              ThisActive = ( Left % PartIndex == ParEnv % MyPe )
            END IF
            Right => Element % BoundaryInfo % Right
            IF( ASSOCIATED( Right ) ) THEN
              ThisActive = ThisActive .OR. &
                  ( Right % PartIndex == ParEnv % MyPe ) 
            END IF
            IF( ThisActive ) THEN
              n1 = n1 + nSplit
            ELSE
              HaloCount = HaloCount + 1
            END IF
          END IF
        ELSE
           n1 = n1 + nSplit
        END IF
      END IF

      IF( OnTheFlyBC ) THEN
        IF( Trgt == Constraint ) n2 = n2 + nSplit
      ELSE
        IF ( Model % BCs(Trgt) % Tag == Constraint ) n2 = n2 + nSplit
      END IF
    END DO

    IF( CheckForHalo ) THEN
      CALL Info(Caller,'Number of halo elements eliminated: '&
          //TRIM(I2S(HaloCount)),Level=12)
    END IF

    IF ( n1 <= 0 .OR. n2 <= 0 ) THEN
      ! This is too conservative in parallel
      ! CALL Warn(Caller,'There are no active boundaries!')
      Success = .FALSE.
      RETURN
    END IF


!   Initialize mesh structures for boundaries, this
!   is for getting the mesh projector:
!   ------------------------------------------------
    BMesh1 % Parent => Mesh
    BMesh2 % Parent => Mesh

    WRITE(Message,'(A,I0,A,I0)') 'Number of interface elements: ',n1,', ',n2
    CALL Info(Caller,Message,Level=9)    
    
    CALL AllocateVector( BMesh1 % Elements,n1 )
    CALL AllocateVector( BMesh2 % Elements,n2 )
    CALL AllocateVector( Perm1, Mesh % NumberOfNodes )
    CALL AllocateVector( Perm2, Mesh % NumberOfNodes )

    IF( TagNormalFlip ) THEN
      ALLOCATE( BMesh1 % PeriodicFlip(n1) )
      ALLOCATE( BMesh2 % PeriodicFlip(n2) )
      BMesh1 % PeriodicFlip = .FALSE.
      BMesh2 % PeriodicFlip = .FALSE.      
    END IF
    
 
!   Fill in the mesh element structures with the
!   boundary elements:
!   ---------------------------------------------
    n1 = 0
    n2 = 0
    Perm1 = 0
    Perm2 = 0
    BMesh1 % MaxElementNodes = 0
    BMesh2 % MaxElementNodes = 0


    DO i=1, Mesh % NumberOfBoundaryElements
      Element => Elements(i)
      
      ElemCode = Element % Type % ElementCode 
      IF (ElemCode <= 200) CYCLE

      IF( TagNormalFlip ) THEN            
        n = Element % TYPE % NumberOfNodes
        NodeIndexes => Element % NodeIndexes

        ElementNodes % x(1:n) = Mesh % Nodes % x(NodeIndexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(NodeIndexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(NodeIndexes(1:n))           
        
        Normal = NormalVector( Element,ElementNodes,Check=.TRUE.,&
            Parent = Element % BoundaryInfo % Left, Turn = Turn )        
      END IF
      
      nSplit = 1
      IF( SplitQuadratic ) THEN
        IF( ElemCode == 306 .OR. ElemCode == 409 ) THEN
          nSplit = 4
        ELSE IF( ElemCode == 408 ) THEN
          nSplit = 5
        END IF
      END IF
       
      Constraint = Element % BoundaryInfo % Constraint
      
      ThisActive = ( Model % BCs(This) % Tag == Constraint ) 
      IF( ThisActive .AND. CheckForHalo ) THEN
        IF( NarrowHalo ) THEN
          IF( .NOT. ANY(ActiveNode(Element % NodeIndexes) ) ) THEN
            ThisActive = .FALSE.
          END IF
        ELSE IF( NoHalo ) THEN
          ThisActive = .FALSE.
          Left => Element % BoundaryInfo % Left 
          IF( ASSOCIATED( Left ) ) THEN
            ThisActive = ( Left % PartIndex == ParEnv % MyPe )
          END IF
          Right => Element % BoundaryInfo % Right
          IF( ASSOCIATED( Right ) ) THEN
            ThisActive = ThisActive .OR. &
                ( Right % PartIndex == ParEnv % MyPe ) 
          END IF
        END IF
      END IF

      IF( OnTheFlyBC ) THEN
        TargetActive = ( Trgt == Constraint )
      ELSE
        TargetActive = ( Model % BCs(Trgt) % Tag == Constraint ) 
      END IF

      IF(.NOT. (ThisActive .OR. TargetActive ) ) CYCLE
      
      ! Set the pointers accordingly so we need to code the complex stuff
      ! only once.
      IF ( ThisActive ) THEN
        n1 = n1 + nSplit
        ind = n1
        PMesh => BMesh1
        PPerm => Perm1
      ELSE
        n2 = n2 + nSplit
        ind = n2
        PMesh => BMesh2
        PPerm => Perm2
      END IF

      
      IF( nSplit > 1 ) THEN
        IF( ElemCode == 408 ) THEN
          SplitSizes(1:nSplit) = [ 4,3,3,3,3 ]
          DO ii=1,nSplit
            jj = ind-nSplit+ii
            m = SplitSizes(ii)
            
            SELECT CASE (ii)
            CASE( 1 )
              SplitMap(1:m) = [ 5,6,7,8 ]
            CASE( 2 )
              SplitMap(1:m) = [ 1, 5, 8 ]
            CASE( 3 ) 
              SplitMap(1:m) = [ 2, 6, 5 ]
            CASE( 4 )
              SplitMap(1:m) = [ 3, 7, 6 ]
            CASE( 5 ) 
              SplitMap(1:m) = [ 4, 8, 7 ]
            END SELECT

            CALL AllocateVector(PMesh % Elements(jj) % NodeIndexes, m )
            PMesh % Elements(jj) % NodeIndexes(1:m) = &
                Element % NodeIndexes(SplitMap(1:m))
            PMesh % Elements(jj) % TYPE => GetElementType(101*m)
            IF( ThisActive ) THEN
              PMesh % Elements(jj) % BoundaryInfo => Element % BoundaryInfo 
            END IF
          END DO          
          PMesh % MaxElementNodes = MAX( PMesh % MaxElementNodes, 4 )

        ELSE IF( ElemCode == 409 ) THEN
          SplitSizes(1:n) = [ 4,4,4,4 ]
          DO ii=1,nSplit
            jj = ind-nSplit+ii
            m = SplitSizes(ii)
            
            SELECT CASE (ii)
            CASE( 1 )
              SplitMap(1:m) = [ 1, 5, 9, 8 ]
            CASE( 2 )
              SplitMap(1:m) = [ 2, 6, 9, 5 ]
            CASE( 3 ) 
              SplitMap(1:m) = [ 3, 7, 9, 6 ]
            CASE( 4 ) 
              SplitMap(1:m) = [ 4, 8, 9, 7 ]
            END SELECT

            CALL AllocateVector(PMesh % Elements(jj) % NodeIndexes, m )
            PMesh % Elements(jj) % NodeIndexes(1:m) = &
                Element % NodeIndexes(SplitMap(1:m))
            PMesh % Elements(jj) % TYPE => GetElementType(101*m)
            IF( ThisActive ) THEN
              PMesh % Elements(jj) % BoundaryInfo => Element % BoundaryInfo 
            END IF
          END DO
          PMesh % MaxElementNodes = MAX( PMesh % MaxElementNodes, 4 )
          
        ELSE IF( ElemCode == 306 ) THEN
          SplitSizes(1:n) = [ 3,3,3,3 ]
          DO ii=1,nSplit
            jj = ind-nSplit+ii
            m = SplitSizes(ii)
            
            SELECT CASE (ii)
            CASE( 1 )
              SplitMap(1:m) = [ 1, 4, 6 ]
            CASE( 2 )
              SplitMap(1:m) = [ 2, 5, 4 ]
            CASE( 3 ) 
              SplitMap(1:m) = [ 3, 6, 5 ]
            CASE( 4 ) 
              SplitMap(1:m) = [ 4, 5, 6 ]
            END SELECT

            CALL AllocateVector(PMesh % Elements(j) % NodeIndexes, m )
            PMesh % Elements(jj) % NodeIndexes(1:m) = &
                Element % NodeIndexes(SplitMap(1:m))
            PMesh % Elements(jj) % TYPE => GetElementType(101*m)
            IF( ThisActive ) THEN
              PMesh % Elements(jj) % BoundaryInfo => Element % BoundaryInfo 
            END IF
          END DO
          PMesh % MaxElementNodes = MAX( PMesh % MaxElementNodes, 3 )
        END IF
        n = Element % TYPE % NumberOfNodes             
        PPerm( Element % NodeIndexes(1:n) ) = 1

      ELSE
        n = Element % TYPE % NumberOfNodes             
        PMesh % MaxElementNodes = MAX( PMesh % MaxElementNodes, n )
        PMesh % Elements(ind) = Element

        IF( TagNormalFlip ) THEN
          PMesh % PeriodicFlip(ind) = Turn
        END IF
                  
        CALL AllocateVector(PMesh % Elements(ind) % NodeIndexes,n )
      
        IF( Mesh % NumberOfFaces == 0 .OR. Mesh % NumberOfEdges == 0 ) THEN
          PMesh % Elements(ind) % NodeIndexes(1:n) = Element % NodeIndexes(1:n)
          PPerm( Element % NodeIndexes(1:n) ) = 1
        ELSE
          ! If we have edge dofs we want the face element be associated with the 
          ! face list since that only has properly defined edge indexes.
          Parent => Element % BoundaryInfo % Left
          IF(.NOT. ASSOCIATED( Parent ) ) THEN
            Parent => Element % BoundaryInfo % Right
          END IF

          q => Find_Face(Mesh,Parent,Element)

          PMesh % Elements(ind) % NodeIndexes(1:n) = q % NodeIndexes(1:n)

          ! set the elementindex to be faceindex as it may be needed
          ! for the edge elements.
          PMesh % Elements(ind) % ElementIndex = q % ElementIndex

          IF(ASSOCIATED(q % Pdefs)) THEN
            ALLOCATE(Pmesh % Elements(ind) % Pdefs)
            PMesh % Elements(ind) % PDefs = q % Pdefs
          END IF

          ! Set also the owner partition
          !       PMesh % Elements(ind) % PartIndex = q % PartIndex

          en = q % TYPE % NumberOfEdges
          ALLOCATE(PMesh % Elements(ind) % EdgeIndexes(en))
          Pmesh % Elements(ind) % EdgeIndexes(1:en) = q % EdgeIndexes(1:en)

          PPerm( q % NodeIndexes(1:n) ) = 1
        END IF
      END IF
        

    END DO
  
!   Fill in the mesh node structures with the
!   boundary nodes:
!   -----------------------------------------
    BMesh1 % NumberOfBulkElements = n1
    BMesh2 % NumberOfBulkElements = n2

    BMesh2 % NumberOfNodes = COUNT(Perm2 > 0)
    BMesh1 % NumberOfNodes = COUNT(Perm1 > 0)

    ! As there were some active boundary elements this condition should 
    ! really never be possible   
    IF (BMesh1 % NumberOfNodes==0 .OR. BMesh2 % NumberOfNOdes==0) THEN
      CALL Fatal(Caller,'No active nodes on periodic boundary!')
    END IF

    WRITE(Message,'(A,I0,A,I0)') 'Number of interface nodes: ',&
        BMesh1 % NumberOfNodes, ', ',BMesh2 % NumberOfNOdes
    CALL Info(Caller,Message,Level=9)    
    
    ALLOCATE( BMesh1 % Nodes )
    CALL AllocateVector( BMesh1 % Nodes % x, BMesh1 % NumberOfNodes ) 
    CALL AllocateVector( BMesh1 % Nodes % y, BMesh1 % NumberOfNodes ) 
    CALL AllocateVector( BMesh1 % Nodes % z, BMesh1 % NumberOfNodes )
    
    ALLOCATE( BMesh2 % Nodes )
    CALL AllocateVector( BMesh2 % Nodes % x, BMesh2 % NumberOfNodes ) 
    CALL AllocateVector( BMesh2 % Nodes % y, BMesh2 % NumberOfNodes ) 
    CALL AllocateVector( BMesh2 % Nodes % z, BMesh2 % NumberOfNodes )
    
    CALL AllocateVector( Bmesh1 % InvPerm, BMesh1 % NumberOfNodes )
    CALL AllocateVector( Bmesh2 % InvPerm, BMesh2 % NumberOfNodes )

    ! Now, create the master and target meshes that only include the active elements
    !---------------------------------------------------------------------------
    k1 = 0; k2 = 0
    DO i=1,Mesh % NumberOfNodes

      IF ( Perm1(i) > 0 ) THEN
        k1 = k1 + 1
        Perm1(i) = k1
        BMesh1 % InvPerm(k1) = i

        BMesh1 % Nodes % x(k1) = Mesh % Nodes % x(i)
        BMesh1 % Nodes % y(k1) = Mesh % Nodes % y(i)
        BMesh1 % Nodes % z(k1) = Mesh % Nodes % z(i)
      END IF
      
      IF ( Perm2(i) > 0 ) THEN
        k2 = k2 + 1
        Perm2(i) = k2
        BMesh2 % InvPerm(k2) = i
        
        BMesh2 % Nodes % x(k2) = Mesh % Nodes % x(i)
        BMesh2 % Nodes % y(k2) = Mesh % Nodes % y(i)
        BMesh2 % Nodes % z(k2) = Mesh % Nodes % z(i)
      END IF
    END DO

!   Finally, Renumber the element node pointers to use
!   only boundary nodes:
!   ---------------------------------------------------

    DO i=1,n1
      BMesh1 % Elements(i) % NodeIndexes = Perm1(BMesh1 % Elements(i) % NodeIndexes)
    END DO

    DO i=1,n2
      BMesh2 % Elements(i) % NodeIndexes = Perm2(BMesh2 % Elements(i) % NodeIndexes)
    END DO
    DEALLOCATE( Perm1, Perm2 )

    IF( CheckForHalo ) DEALLOCATE( ActiveNode ) 

    Success = .TRUE.

  END SUBROUTINE CreateInterfaceMeshes
  !---------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !> Given two meshes that should occupy the same domain in space
  !> use rotation, scaling and translation to achieve this goal.
  !---------------------------------------------------------------------------
  SUBROUTINE OverlayIntefaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    LOGICAL :: GotIt, GotRotate
    REAL(KIND=dp) :: x1_min(3),x1_max(3),x2_min(3),x2_max(3),x2r_min(3),x2r_max(3)
    REAL(KIND=dp) :: x(4), RotMatrix(4,4),TrsMatrix(4,4),SclMatrix(4,4), &
           TrfMatrix(4,4),Identity(4,4),Angles(3),Alpha,scl(3),s1,s2
    REAL(KIND=dp), POINTER :: PArray(:,:)
    INTEGER :: i,j,k
    CHARACTER(*), PARAMETER :: Caller='OverlayInterfaceMeshes'
 
    ! First, check the bounding boxes
    !---------------------------------------------------------------------------
    x1_min(1) = MINVAL( BMesh1 % Nodes % x )
    x1_min(2) = MINVAL( BMesh1 % Nodes % y )
    x1_min(3) = MINVAL( BMesh1 % Nodes % z )
    
    x1_max(1) = MAXVAL( BMesh1 % Nodes % x )
    x1_max(2) = MAXVAL( BMesh1 % Nodes % y )
    x1_max(3) = MAXVAL( BMesh1 % Nodes % z )

    WRITE(Message,'(A,3ES15.6)') 'Minimum values for this periodic BC:  ',x1_min
    CALL Info(Caller,Message,Level=8)    
    WRITE(Message,'(A,3ES15.6)') 'Maximum values for this periodic BC:  ',x1_max
    CALL Info(Caller,Message,Level=8)    

    x2_min(1) = MINVAL( BMesh2 % Nodes % x )
    x2_min(2) = MINVAL( BMesh2 % Nodes % y )
    x2_min(3) = MINVAL( BMesh2 % Nodes % z )
    
    x2_max(1) = MAXVAL( BMesh2 % Nodes % x )
    x2_max(2) = MAXVAL( BMesh2 % Nodes % y )
    x2_max(3) = MAXVAL( BMesh2 % Nodes % z )
    
    WRITE(Message,'(A,3ES15.6)') 'Minimum values for target periodic BC:',x2_min
    CALL Info(Caller,Message,Level=8)    
    WRITE(Message,'(A,3ES15.6)') 'Maximum values for target periodic BC:',x2_max
    CALL Info(Caller,Message,Level=8)    

!    If whole transformation matrix given, it will be used directly
!    --------------------------------------------------------------
    Parray => ListGetConstRealArray( BParams,'Periodic BC Matrix', Gotit )
    IF ( GotIt ) THEN
      DO i=1,SIZE(Parray,1)
        DO j=1,SIZE(Parray,2)
          TrfMatrix(i,j) = Parray(j,i)
        END DO
      END DO
    ELSE    
      ! Otherwise check for rotation, scaling and translation
      !------------------------------------------------------

      ! Initialize the mapping matrices
      Identity = 0.0d0
      DO i=1,4
        Identity(i,i) = 1.0d0
      END DO      
      TrsMatrix = Identity
      RotMatrix = Identity
      SclMatrix = Identity
      
      !   Rotations:
      !   These are called first since they are not accounted for in the 
      !   automatic scaling and translation.
      !   ---------------------------------------------------------------      
      Angles = 0.0_dp
      Parray => ListGetConstRealArray( BParams,'Periodic BC Rotate', GotRotate )
      IF( GotRotate ) THEN
        Angles(1:3) = Parray(1:3,1)   
      ELSE
        IF( ListGetLogical( BParams,'Periodic BC Rotate Automatic', GotIt) ) THEN
          CALL CheckInterfaceMeshAngle( BMesh1, BMesh2, Angles, GotRotate ) 
        END IF
      END IF

      IF ( GotRotate ) THEN
        WRITE(Message,'(A,3ES15.6)') 'Rotating target with: ',Angles
        CALL Info(Caller,Message,Level=8)    
        
        DO i=1,3
          Alpha = Angles(i) * PI / 180.0_dp
          IF( ABS(Alpha) < TINY(Alpha) ) CYCLE 
          TrfMatrix = Identity
          
          SELECT CASE(i)
          CASE(1)
            TrfMatrix(2,2) =  COS(Alpha)
            TrfMatrix(2,3) = -SIN(Alpha) 
            TrfMatrix(3,2) =  SIN(Alpha)
            TrfMatrix(3,3) =  COS(Alpha)
          CASE(2)
            TrfMatrix(1,1) =  COS(Alpha)
            TrfMatrix(1,3) = -SIN(Alpha)
            TrfMatrix(3,1) =  SIN(Alpha)
            TrfMatrix(3,3) =  COS(Alpha)
          CASE(3)
            TrfMatrix(1,1) =  COS(Alpha)
            TrfMatrix(1,2) = -SIN(Alpha)
            TrfMatrix(2,1) =  SIN(Alpha)
            TrfMatrix(2,2) =  COS(Alpha)
          END SELECT
          
          RotMatrix = MATMUL( RotMatrix, TrfMatrix )
        END DO
        
        DO i = 1, BMesh2 % NumberOfNodes          
          x(1) = BMesh2 % Nodes % x(i)
          x(2) = BMesh2 % Nodes % y(i)
          x(3) = BMesh2 % Nodes % z(i)
          
          x(4) = 1.0_dp
          x = MATMUL( RotMatrix, x )
          
          BMesh2 % Nodes % x(i) = x(1)
          BMesh2 % Nodes % y(i) = x(2)
          BMesh2 % Nodes % z(i) = x(3)
        END DO
        
        x2r_min(1) = MINVAL( BMesh2 % Nodes % x )
        x2r_min(2) = MINVAL( BMesh2 % Nodes % y )
        x2r_min(3) = MINVAL( BMesh2 % Nodes % z )
        
        x2r_max(1) = MAXVAL( BMesh2 % Nodes % x )
        x2r_max(2) = MAXVAL( BMesh2 % Nodes % y )
        x2r_max(3) = MAXVAL( BMesh2 % Nodes % z )
        
        WRITE(Message,'(A,3ES15.6)') 'Minimum values for rotated target:',x2r_min
        CALL Info(Caller,Message,Level=8)    
        
        WRITE(Message,'(A,3ES15.6)') 'Maximum values for rotated target:',x2r_max
        CALL Info(Caller,Message,Level=8)    
      ELSE
        x2r_min = x2_min
        x2r_max = x2_max
      END IF
   
!   Scaling:
!   This is either given or enforced by requiring bounding boxes to be of the same size 
!   -----------------------------------------------------------------------------------
      Parray => ListGetConstRealArray( BParams,'Periodic BC Scale', Gotit )      
      IF ( GotIt ) THEN
        DO i=1,SIZE(Parray,1)
          SclMatrix(i,i) = Parray(i,1)
        END DO
      ELSE
        ! Define scaling from the bounding boxes
        ! This assumes isotropic scaling since component-wise scaling 
        ! was prone to errors.
        !------------------------------------------------------
        s1 = SUM( ( x1_max(1:3) - x1_min(1:3) ) ** 2 )
        s2 = SUM( ( x2r_max(1:3) - x2r_min(1:3) ) ** 2 )
        IF( s2 > EPSILON( s2 ) ) THEN
          scl(1:3)  = SQRT( s1 / s2 )
        ELSE
          scl(1:3) = 1.0_dp
        END IF
        
        WRITE(Message,'(A,3ES15.6)') 'Scaling with: ',scl(1:3)
        CALL Info(Caller,Message)
        DO i=1,3 
          SclMatrix(i,i) = scl(i)        
        END DO
      END IF
      
!   Translations:
!   And finally define translations
!   -------------
      Parray => ListGetConstRealArray( BParams,'Periodic BC Translate', Gotit )
      IF ( gotit ) THEN
        DO i=1,SIZE(Parray,1)
          TrsMatrix(4,i) = Parray(i,1)
        END DO
      ELSE
        ! Define translations so that the lower left corner is the same
        !-------------------------------------------------------------
        DO i=1,3
          TrsMatrix(4,i) = x1_min(i) - SclMatrix(i,i) * x2r_min(i)
        END DO
      END IF
      WRITE(Message,'(A,3ES15.6)') 'Translation: ',TrsMatrix(4,1:3)
      CALL Info(Caller,Message)
      TrfMatrix = MATMUL( SclMatrix, TrsMatrix )
    END IF

!    Now transform the coordinates:
!    ------------------------------
    DO i=1,BMesh2 % NumberOfNodes
      x(1) = BMesh2 % Nodes % x(i)
      x(2) = BMesh2 % Nodes % y(i)
      x(3) = BMesh2 % Nodes % z(i)
      x(4) = 1.0d0
      x = MATMUL( x, TrfMatrix ) 
      BMesh2 % Nodes % x(i) = x(1) / x(4)
      BMesh2 % Nodes % y(i) = x(2) / x(4) 
      BMesh2 % Nodes % z(i) = x(3) / x(4)
    END DO

    IF(.FALSE.) THEN
      x2r_min(1) = MINVAL( BMesh2 % Nodes % x )
      x2r_min(2) = MINVAL( BMesh2 % Nodes % y )
      x2r_min(3) = MINVAL( BMesh2 % Nodes % z )
      
      x2r_max(1) = MAXVAL( BMesh2 % Nodes % x )
      x2r_max(2) = MAXVAL( BMesh2 % Nodes % y )
      x2r_max(3) = MAXVAL( BMesh2 % Nodes % z )
      
      WRITE(Message,'(A,3ES15.6)') 'Minimum values for transformed target:',x2r_min
      CALL Info(Caller,Message,Level=8)    
      
      WRITE(Message,'(A,3ES15.6)') 'Maximum values for transformed target:',x2r_max
      CALL Info(Caller,Message,Level=8)    
    END IF

  END SUBROUTINE OverlayIntefaceMeshes
  !---------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Given two interface meshes for nonconforming rotating boundaries make 
  !> a coordinate transformation to each node of the slave boundary (BMesh1) so that 
  !> they hit the master boundary (BMesh2). In case of anti-periodic projector 
  !> mark the nodes that need an odd number of periods.
  !---------------------------------------------------------------------------
  SUBROUTINE PreRotationalProjector(BMesh1, BMesh2, MirrorNode )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    !--------------------------------------------------------------------------
    LOGICAL :: AntiPeriodic
    REAL(KIND=dp) :: F2min,F2max,dFii2,Fii
    INTEGER :: i, Nfii, SectorMax
    INTEGER, ALLOCATABLE :: SectorCount(:)
    CHARACTER(*), PARAMETER :: Caller='PreRotationalProjector'
 
    AntiPeriodic = ALLOCATED( MirrorNode )
    IF( AntiPeriodic ) MirrorNode = .FALSE.

    F2Min =  MINVAL( BMesh2 % Nodes % x )
    F2Max =  MAXVAL( BMesh2 % Nodes % x )
    dFii2 = F2Max - F2Min
    SectorMax = CEILING( 360.0_dp / dFii2 ) 

    WRITE( Message,'(A,I0)') 'Maximum number of sectors: ',SectorMax
    CALL Info(Caller,Message,Level=8)

    ALLOCATE( SectorCount(-SectorMax:SectorMax))
    SectorCount = 0

    DO i = 1, BMesh1 % NumberOfNodes
      Fii = BMesh1 % Nodes % x(i)      
      Nfii = FLOOR( (Fii-F2min) / dFii2 )
      BMesh1 % Nodes % x(i) = BMesh1 % Nodes % x(i) - Nfii * dFii2
      SectorCount(Nfii) = SectorCount(Nfii) + 1     
      IF( AntiPeriodic ) THEN
        IF( MODULO(Nfii,2) /= 0 ) THEN
          MirrorNode(i) = .TRUE.
        END IF
      END IF
    END DO

    IF( SectorCount(0) < BMesh1 % NumberOfNodes ) THEN
      CALL Info(Caller,'Number of nodes by sectors',Level=8)
      DO i=-SectorMax,SectorMax
        IF( SectorCount(i) > 0 ) THEN
          WRITE( Message,'(A,I0,A,I0)') 'Sector:',i,'   Nodes:',SectorCount(i)
          CALL Info('MatchInterfaceNodes',Message,Level=8)
        END IF
      END DO
      IF( AntiPeriodic ) THEN
        WRITE( Message,'(A,I0)') 'Number of mirror nodes:',COUNT(MirrorNode)
        CALL Info(Caller,Message,Level=8)
      END IF
    ELSE
      CALL Info(Caller,'No nodes needed mapping')
    END IF

  END SUBROUTINE PreRotationalProjector
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!> Postprocess projector so that it changes the sign of the anti-periodic
!> entries as assigns by the MirrorNode flag.
!------------------------------------------------------------------------------
  SUBROUTINE PostRotationalProjector( Proj, MirrorNode )
!------------------------------------------------------------------------------
    TYPE(Matrix_t) :: Proj                 !< Projection matrix
    LOGICAL, ALLOCATABLE :: MirrorNode(:)  !< Is the node a mirror node or not
!--------------------------------------------------------------------------
    INTEGER, POINTER :: Cols(:),Rows(:)            
    REAL(KIND=dp), POINTER :: Values(:)    
    INTEGER :: i,j,n
!------------------------------------------------------------------------------

    IF( .NOT. ALLOCATED( MirrorNode ) ) RETURN
    IF( COUNT( MirrorNode ) == 0 ) RETURN

    n = Proj % NumberOfRows
    Rows => Proj % Rows
    Cols => Proj % Cols
    Values => Proj % Values

    DO i=1,n
      IF( MirrorNode(i) ) THEN
        DO j = Rows(i),Rows(i+1)-1
          Values(j) = -Values(j)
        END DO
      END IF
    END DO

!------------------------------------------------------------------------------
  END SUBROUTINE PostRotationalProjector
!------------------------------------------------------------------------------

   !----------------------------------------------------------------------------------------
  !> Given a temporal triangle "ElementT", calculate mass matrix contributions for projection
  !> for the slave element "Element" and master element "ElementM".
  !> The nubmering associated to these surface meshes is InvPerm and InvPermM, respectively. 
  !> This is lifted at an outer level in the hope that it would be called by number of
  !> routines in the future.
  !----------------------------------------------------------------------------------------
  SUBROUTINE TemporalTriangleMortarAssembly(ElementT, NodesT, Element, Nodes, ElementM, NodesM, &
      Biorthogonal, DualMaster, DualLCoeff, NoGaussPoints, Projector, NodeScale, &
      NodePerm, InvPerm, InvPermM, SumArea ) 
    !----------------------------------------------------------------------------------------
    TYPE(Element_t) :: ElementT
    TYPE(Element_t), POINTER :: Element, ElementM
    TYPE(Nodes_t) :: NodesT, Nodes, NodesM
    LOGICAL :: Biorthogonal, DualMaster, DualLCoeff
    INTEGER :: NoGaussPoints
    TYPE(Matrix_t) :: Projector
    REAL(KIND=dp) :: NodeScale, SumArea
    INTEGER, POINTER :: NodePerm(:), InvPerm(:), InvPermM(:)
    !----------------------------------------------------------------------------------------

    TYPE(Element_t), POINTER :: ElementP, ElementLin
    TYPE(GaussIntegrationPoints_t) :: IPT
    REAL(KIND=dp) :: area, xt, yt, zt = 0.0_dp, u, v, w, um, vm, wm, &
        detJ, val, val_dual, weight
    REAL(KIND=dp), ALLOCATABLE :: BasisT(:),Basis(:), BasisM(:), MASS(:,:), CoeffBasis(:)
    INTEGER :: i,j,jj,n,ne,nM,neM,ElemCode,LinCode,ElemCodeM,LinCodeM,nip,nrow,AllocStat
    INTEGER, POINTER :: Indexes(:),IndexesM(:)
    LOGICAL :: Stat, AllocationsDone = .FALSE.

    SAVE :: BasisT, Basis, BasisM, CoeffBasis, MASS

    IF(.NOT. AllocationsDone ) THEN
      n = CurrentModel % Mesh % MaxElementNodes
      ALLOCATE( BasisT(3),Basis(n), BasisM(n), CoeffBasis(n), MASS(n,n), STAT = AllocStat )
      IF( AllocStat /= 0 ) CALL Fatal('LocalMortarAssembly','Allocation error!')
      AllocationsDone = .TRUE.
    END IF


    n = Element % TYPE % NumberOfNodes
    ne = Element % TYPE % ElementCode / 100      
    ElemCode = Element % TYPE % ElementCode 
    LinCode = 101 * ne
    Indexes => Element % NodeIndexes

    nM = ElementM % TYPE % NumberOfNodes
    neM = ElementM % TYPE % ElementCode / 100      
    ElemCodeM = Element % TYPE % ElementCode 
    LinCodeM = 101 * neM
    IndexesM => ElementM % NodeIndexes

    IF( NoGaussPoints > 0 ) THEN
      IPT = GaussPoints( ElementT, NoGaussPoints, PreferenceElement = .FALSE. )
    ELSE
      IPT = GaussPoints( ElementT, PreferenceElement = .FALSE. )
    END IF

    IF(BiOrthogonal) THEN
      MASS  = 0
      CoeffBasis = 0
      area = 0._dp
      DO nip=1, IPT % n 
        stat = ElementInfo( ElementT,NodesT,IPT % u(nip),&
            IPT % v(nip),IPT % w(nip),detJ,BasisT)
        IF(.NOT. Stat ) EXIT

        ! We will actually only use the global coordinates and the integration weight 
        ! from the temporal mesh. 

        ! Global coordinates of the integration point
        xt = SUM( BasisT(1:3) * NodesT % x(1:3) )
        yt = SUM( BasisT(1:3) * NodesT % y(1:3) )

        ! Integration weight for current integration point
        Weight = DetJ * IPT % s(nip) 
        area = area + weight

        ! Integration point at the slave element
        IF( ElemCode /= LinCode ) THEN
          ElementLin % TYPE => GetElementType( LinCode, .FALSE. )
          ElementLin % NodeIndexes => Element % NodeIndexes
          ElementP => ElementLin
          CALL GlobalToLocal( u, v, w, xt, yt, zt, ElementP, Nodes )
        ELSE
          CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
        END IF

        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        IF(.NOT. Stat) CYCLE

        DO i=1,n
          DO j=1,n
            MASS(i,j) = MASS(i,j) + weight * Basis(i) * Basis(j)
          END DO
          CoeffBasis(i) = CoeffBasis(i) + Weight * Basis(i)
        END DO
      END DO

      ! Even if there would be multiple ip points, area is still the same...
      IF(Area<1.d-12) RETURN

      CALL InvertMatrix( MASS, n )

      DO i=1,n
        DO j=1,n
          MASS(i,j) = MASS(i,j) * CoeffBasis(i)
        END DO
      END DO
    END IF

    ! Integration over the temporal element using integration points of that element
    DO nip=1, IPT % n 
      stat = ElementInfo( ElementT,NodesT,IPT % u(nip),&
          IPT % v(nip),IPT % w(nip),detJ,BasisT)
      IF(.NOT. Stat) EXIT

      ! We will actually only use the global coordinates and the integration weight 
      ! from the temporal mesh. 

      ! Global coordinates of the integration point
      xt = SUM( BasisT(1:3) * NodesT % x(1:3) )
      yt = SUM( BasisT(1:3) * NodesT % y(1:3) )

      ! Integration weight for current integration point
      weight = DetJ * IPT % s(nip)
      sumarea = sumarea + weight

      ! Integration point at the slave element
      IF( ElemCode /= LinCode ) THEN
        ElementLin % TYPE => GetElementType( LinCode, .FALSE. )
        ElementLin % NodeIndexes => Element % NodeIndexes
        ElementP => ElementLin
        CALL GlobalToLocal( u, v, w, xt, yt, zt, ElementP, Nodes )
      ELSE
        CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
      END IF

      stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

      ! Integration point at the master element
      IF( ElemCodeM /= LinCodeM ) THEN
        ElementLin % TYPE => GetElementType( LinCodeM, .FALSE. )
        ElementLin % NodeIndexes => ElementM % NodeIndexes
        ElementP => ElementLin
        CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementP, NodesM )
      ELSE
        CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
      END IF

      stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
      IF(.NOT. Stat) CYCLE

      ! Add the nodal dofs
      IF(BiOrthogonal) THEN
        CoeffBasis = 0._dp
        DO i=1,n
          DO j=1,n
            CoeffBasis(i) = CoeffBasis(i) + MASS(i,j) * Basis(j)
          END DO
        END DO
      END IF

      DO j=1,n 
        jj = Indexes(j)                                    

        nrow = NodePerm(InvPerm(jj))
        IF( nrow == 0 ) CYCLE

        Projector % InvPerm(nrow) = InvPerm(jj)
        val = Basis(j) * weight
        IF(Biorthogonal) val_dual = CoeffBasis(j) * weight

        DO i=1,n
          IF( ABS( val * Basis(i) ) < 1.0d-10 ) CYCLE

          !Nslave = Nslave + 1
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              InvPerm(Indexes(i)), Basis(i) * val ) 

          IF(BiOrthogonal) THEN
            CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                InvPerm(Indexes(i)), Basis(i) * val_dual ) 
          END IF
        END DO

        DO i=1,nM
          IF( ABS( val * BasisM(i) ) < 1.0d-12 ) CYCLE

          !Nmaster = Nmaster + 1
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              InvPermM(IndexesM(i)), -NodeScale * BasisM(i) * val )                   

          IF(BiOrthogonal) THEN
            IF(DualMaster .OR. DualLCoeff) THEN
              CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                  InvPermM(IndexesM(i)), -NodeScale * BasisM(i) * val_dual ) 
            ELSE
              CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                  InvPermM(IndexesM(i)), -NodeScale * BasisM(i) * val ) 
            END IF
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE TemporalTriangleMortarAssembly

  
  !---------------------------------------------------------------------------
  !> Create a projector for mapping between interfaces using the Galerkin method
  !> A temporal mesh structure with a node for each Gaussian integration point is 
  !> created. Then this projector matrix is transferred to a projector on the nodal
  !> coordinates.   
  !---------------------------------------------------------------------------
   FUNCTION NormalProjector(BMesh2, BMesh1, BC) RESULT ( Projector )
  !---------------------------------------------------------------------------
    USE Lists

    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
    INTEGER, POINTER :: Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:)
    TYPE(Mesh_t), POINTER :: Mesh
    TYPE(Matrix_t), POINTER :: DualProjector
    LOGICAL :: Found, Parallel, BiOrthogonalBasis, &
        CreateDual, DualSlave, DualMaster, DualLCoeff 
    REAL(KIND=dp) :: NodeScale 
    INTEGER, POINTER :: NodePerm(:)
    TYPE(Element_t), POINTER :: Element
    INTEGER :: i,n,m
    CHARACTER(*), PARAMETER :: Caller='NormalProjector'
    
    CALL Info(Caller,'Creating projector between 3D surfaces',Level=7)
    
    Parallel = ( ParEnv % PEs > 1 )
    Mesh => CurrentModel % Mesh
    BMesh1 % Parent => NULL()
    BMesh2 % Parent => NULL()
    
    InvPerm1 => BMesh1 % InvPerm
    InvPerm2 => BMesh2 % InvPerm

    ! Create a list matrix that allows for unspecified entries in the matrix 
    ! structure to be introduced.
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_GALERKIN
    
    CreateDual = ListGetLogical( BC,'Create Dual Projector',Found ) 
    IF( CreateDual ) THEN
      DualProjector => AllocateMatrix()
      DualProjector % FORMAT = MATRIX_LIST
      DualProjector % ProjectorType = PROJECTOR_TYPE_GALERKIN
      Projector % EMatrix => DualProjector
    END IF
    
    ! Check whether biorthogonal basis for projectors requested:
    ! ----------------------------------------------------------
    BiOrthogonalBasis = ListGetLogical( BC, 'Use Biorthogonal Basis', Found)
    ! If we want to eliminate the constraints we have to have a biortgonal basis
    IF(.NOT. Found ) THEN
      BiOrthogonalBasis = ListGetLogical( CurrentModel % Solver % Values, &
          'Eliminate Linear Constraints',Found )
      IF( BiOrthogonalBasis ) THEN
        CALL Info(Caller,&
            'Enforcing > Use Biorthogonal Basis < to True to enable elimination',Level=8)
        CALL ListAddLogical( BC, 'Use Biorthogonal Basis',.TRUE. )
      END IF
    END IF
    
    IF (BiOrthogonalBasis) THEN
      DualSlave  = ListGetLogical(BC, 'Biorthogonal Dual Slave', Found)
      IF(.NOT.Found) DualSlave  = .TRUE.

      DualMaster = ListGetLogical(BC, 'Biorthogonal Dual Master', Found)
      IF(.NOT.Found) DualMaster = .TRUE.

      DualLCoeff = ListGetLogical(BC, 'Biorthogonal Dual Lagrange Coefficients', Found)
      IF(.NOT.Found) DualLCoeff = .FALSE.

      IF(DualLCoeff) THEN
        DualSlave  = .FALSE.
        DualMaster = .FALSE.
        CALL ListAddLogical( CurrentModel % Solver % Values, 'Use Transpose Values',.FALSE.)
      ELSE
        CALL ListAddLogical( CurrentModel % Solver % Values, 'Use Transpose Values',.TRUE.)
      END IF

      Projector % Child => AllocateMatrix()
      Projector % Child % Format = MATRIX_LIST
      CALL Info(Caller,'Using biorthogonal basis, as requested',Level=8)      
    END IF
    

    ALLOCATE( NodePerm( Mesh % NumberOfNodes ) )
    NodePerm = 0
    
    ! in parallel only consider nodes that truly are part of this partition
    DO i=1,BMesh1 % NumberOfBulkElements
      Element => BMesh1 % Elements(i)        
      IF( Parallel ) THEN
        IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE          
      END IF
      NodePerm( InvPerm1( Element % NodeIndexes ) ) = 1
    END DO
    n = 0
    DO i = 1, Mesh % NumberOfNodes
      IF( NodePerm(i) > 0 ) THEN
        n = n + 1
        NodePerm(i) = n
      END IF
    END DO
    CALL Info(Caller,'Initial number of slave nodes '//TRIM(I2S(n))//&
        ' out of '//TRIM(I2S(BMesh1 % NumberOfNodes ) ), Level = 10 )
    
    ALLOCATE( Projector % InvPerm(n) )
    Projector % InvPerm = 0

    DualMaster = ListGetLogical(BC, 'Biorthogonal Dual Master', Found)
    IF(.NOT.Found) DualMaster = .TRUE.

    NodeScale = ListGetConstReal( BC, 'Mortar BC Scaling', Found)
    IF(.NOT. Found ) NodeScale = 1.0_dp
    
    
    ! Here we create the projector
    !--------------------------------------------------------------
    CALL NormalProjectorWeak3D()
    !--------------------------------------------------------------

    
    ! Now change the matrix format to CRS from list matrix
    !--------------------------------------------------------------
    CALL List_toCRSMatrix(Projector)
    CALL CRS_SortMatrix(Projector,.TRUE.)
    CALL Info(Caller,'Number of rows in projector: '&
        //TRIM(I2S(Projector % NumberOfRows)),Level=12)
    CALL Info(Caller,'Number of entries in projector: '&
        //TRIM(I2S(SIZE(Projector % Values))),Level=12)
  
    IF(ASSOCIATED(Projector % Child)) THEN
      CALL List_toCRSMatrix(Projector % Child)
      CALL CRS_SortMatrix(Projector % Child,.TRUE.)
    END IF

    IF( CreateDual ) THEN
      CALL List_toCRSMatrix(DualProjector)
      CALL CRS_SortMatrix(DualProjector,.TRUE.)
    END IF
    
    m = COUNT( Projector % InvPerm  > 0 ) 
    IF( m > 0 ) THEN
      CALL Info(Caller,'Projector % InvPerm set for dofs: '//TRIM(I2S(m)),Level=7)
    END IF
    m = COUNT( Projector % InvPerm  == 0 ) 
    IF( m > 0 ) THEN
      CALL Warn(Caller,'Projector % InvPerm not set in for dofs: '//TRIM(I2S(m)))
    END IF

    CALL Info(Caller,'Projector created',Level=10)


    
  CONTAINS

    
    !----------------------------------------------------------------------
    ! Create weak projector in a generic 3D case using local coordinates.
    ! For each slave element we move into local normal-tangential coordinates
    ! and use the same coordinate system for the candidate master elements
    ! as well. Only the rought 1st selection is made in the original coordinate
    ! system. Using the n-t coordinate system we can again operate in a local
    ! x-y coordinate system.
    !----------------------------------------------------------------------
    SUBROUTINE NormalProjectorWeak3D()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: i,j,n,jj,ii,sgn0,k,kmax,ind,indM,nip,nn,ne,inds(10),nM,neM,iM,i2,i2M
      INTEGER :: ElemCands, TotCands, ElemHits, TotHits, EdgeHits, CornerHits, &
          MaxErrInd, MinErrInd, InitialHits, ActiveHits, TimeStep, Nrange1, NoGaussPoints, &
          AllocStat, NrangeAve, nrow, SubTri
      TYPE(Element_t), POINTER :: Element, ElementM, ElementP
      TYPE(Element_t) :: ElementT
      TYPE(Element_t), TARGET :: ElementLin
      TYPE(GaussIntegrationPoints_t) :: IP, IPT
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: x(10),y(10),xt,yt,zt,xmax,ymax,xmin,ymin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,u,v,w,RefArea,dArea,&
          SumArea,MaxErr,MinErr,Err,Depth,MinDepth,MaxDepth,phi(10),Point(3),uvw(3), &
          val_dual, zmin, zmax, zave, zminm, zmaxm, uq, vq, TolS, &
          MaxNormalDot, ElemdCoord(3), ElemH, MaxElemH(2), MinElemH(2)
      REAL(KIND=dp) :: A(2,2), B(2), C(2), absA, detA, rlen, &
          x1, x2, y1, y2, x1M, x2M, y1M, y2M, x0, y0, dist
      REAL(KIND=dp) :: TotRefArea, TotSumArea
      REAL(KIND=dp), ALLOCATABLE :: Basis(:) 
      LOGICAL :: Stat, CornerFound(4), CornerFoundM(4)
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: TimestepVar
      TYPE(Mesh_t), POINTER :: pMesh      
      TYPE(Nodes_t) :: Center2
      REAL(KIND=dp) :: Center(3), MaxDistance, Normal(3), Tangent(3), Tangent2(3), &
          NormalM(3), r(3)
      
      ! These are used temporarily for debugging purposes
      INTEGER :: SaveInd, MaxSubElem, MaxSubTriangles, DebugInd, iMesh
      LOGICAL :: SaveElem, DebugElem, SaveErr
      CHARACTER(LEN=20) :: FileName

      CHARACTER(*), PARAMETER :: Caller='NormalProjectorWeak3D'

      CALL Info(Caller,'Creating weak constraints using a generic integrator',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      MaxDistance = ListGetCReal( BC,'Projector Max Distance',Found )
            
      SaveInd = ListGetInteger( BC,'Projector Save Element Index',Found )
      DebugInd = ListGetInteger( BC,'Projector Debug Element Index',Found )
      SaveErr = ListGetLogical( BC,'Projector Save Fraction',Found)
      MaxNormalDot = ListGetCReal( BC,'Max Search Normal',Found)
      IF(.NOT. Found ) MaxNormalDot = -0.1
      
      TimestepVar => VariableGet( Mesh % Variables,'Timestep',ThisOnly=.TRUE. )
      Timestep = NINT( TimestepVar % Values(1) )

      IF( SaveErr ) THEN
        FileName = 'frac_'//TRIM(I2S(TimeStep))//'.dat'
        OPEN( 11,FILE=Filename)
      END IF
     
      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), &
          NodesM % x(n), NodesM % y(n), NodesM % z(n), &
          NodesT % x(3), NodesT % y(3), NodesT % z(3), Basis(n), &
          STAT = AllocStat )
      IF( AllocStat /= 0 ) CALL Fatal(Caller,'Allocation error 1')
                      
      MaxErr = 0.0_dp
      MinErr = HUGE( MinErr )
      MinDepth = HUGE( MinDepth )
      MaxDepth = -HUGE( MaxDepth ) 
      MaxErrInd = 0
      MinErrInd = 0
      zt = 0.0_dp
      NodesT % z = 0.0_dp
      
      ! The temporal triangle used in the numerical integration
      ElementT % TYPE => GetElementType( 303, .FALSE. )
      ElementT % NodeIndexes => IndexesT

      ! Use optionally user defined integration rules           
      NoGaussPoints = ListGetInteger( BC,'Mortar BC Gauss Points',Found ) 
      IF( NoGaussPoints > 0 ) THEN
        IPT = GaussPoints( ElementT, NoGaussPoints, PreferenceElement = .FALSE. )
      ELSE
        IPT = GaussPoints( ElementT, PreferenceElement = .FALSE. )
      END IF
      CALL Info(Caller,'Number of integration points for temporal triangle: '&
          //TRIM(I2S(IPT % n)),Level=7)
      
      TotCands = 0
      TotHits = 0
      EdgeHits = 0
      CornerHits = 0
      InitialHits = 0
      ActiveHits = 0
      TotRefArea = 0.0_dp
      TotSumArea = 0.0_dp
      Point = 0.0_dp
      MaxSubTriangles = 0
      MaxSubElem = 0

      ! Save center of elements for master mesh for fast rough test
      n = BMesh2 % NumberOfBulkElements
      ALLOCATE( Center2 % X(n), Center2 % y(n), Center2 % z(n) )

      MaxElemH = 0.0_dp
      MinElemH = HUGE( ElemH ) 

      ! Calculate maximum and minimum elementsize for slave and master mesh
      DO iMesh=1,2
        IF( iMesh == 1 ) THEN
          pMesh => BMesh1
        ELSE
          pMesh => BMesh2
        END IF

        DO ind=1,pMesh % NumberOfBulkElements
          Element => pMesh % Elements(ind)        
          Indexes => Element % NodeIndexes
          n = Element % TYPE % NumberOfNodes
          ne = Element % TYPE % ElementCode / 100
          
          ! Calculate maximum size of element
          ElemdCoord(1) = MAXVAL( pMesh % Nodes % x(Indexes(1:ne)) ) - &
              MINVAL( pMesh % Nodes % x(Indexes(1:ne)) )
          ElemdCoord(2) = MAXVAL( pMesh % Nodes % y(Indexes(1:ne)) ) - &
              MINVAL( pMesh % Nodes % y(Indexes(1:ne)) )
          ElemdCoord(3) = MAXVAL( pMesh % Nodes % z(Indexes(1:ne)) ) - &
              MINVAL( pMesh % Nodes % z(Indexes(1:ne)) )      

          ElemH = SQRT( SUM( ElemdCoord**2 ) )

          MaxElemH(iMesh) = MAX( MaxElemH(iMesh), ElemH ) 
          MinElemH(iMesh) = MIN( MinElemH(iMesh), ElemH ) 
        
          IF( iMesh == 2 ) THEN
            Center2 % x(ind) = SUM( pMesh % Nodes % x(Indexes(1:ne)) ) / ne
            Center2 % y(ind) = SUM( pMesh % Nodes % y(Indexes(1:ne)) ) / ne
            Center2 % z(ind) = SUM( pMesh % Nodes % z(Indexes(1:ne)) ) / ne
          END IF
          
        END DO

        !PRINT *,'Element size range:',MinElemH(iMesh),MaxElemH(iMesh)
      END DO
      
      ! Use tolerances related to minimum elementsize
      TolS = 1.0d-8 * MINVAL( MinElemH ) 

      ! Maximum theoretical distance of centerpoints  
      ElemH = 0.5 * SUM( MaxElemH )
      
      IF( MaxDistance < ElemH ) THEN
        CALL Info(Caller,'Increasing search distance radius')
        !PRINT *,'MaxDistance:',MaxDistance,ElemH
        MaxDistance = 1.2 * ElemH ! some tolerance!
      END IF
            
      DO ind=1,BMesh1 % NumberOfBulkElements        
        
        ! Optionally save the submesh for specified element, for vizualization and debugging
        SaveElem = ( SaveInd == ind )
        DebugElem = ( DebugInd == ind )

        IF( DebugElem ) THEN
          PRINT *,'Debug element turned on: '//TRIM(I2S(ind))
          PRINT *,'Element is p-element:',isActivePElement(element) 
        END IF

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        ne = Element % TYPE % NumberOfEdges 

        ! The coordinates of the boundary element
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))
        Nodes % z(1:n) = BMesh1 % Nodes % z(Indexes(1:n))

        ! Center in the original coordinates
        Center(1) = SUM( Nodes % x(1:ne) ) / ne
        Center(2) = SUM( Nodes % y(1:ne) ) / ne
        Center(3) = SUM( Nodes % z(1:ne) ) / ne
        
        ! Find the new normal-tangential coordinate system for this particular element
        Normal = NormalVector( Element, Nodes, Check = .FALSE. ) 
        IF( BMesh1 % PeriodicFlip(ind) ) Normal = -Normal
        CALL TangentDirections( Normal,Tangent,Tangent2 )
        
        IF( DebugElem ) THEN
          PRINT *,'Center of element:',Center
          PRINT *,'Normal:',Normal,BMesh1 % PeriodicFlip(ind)
          PRINT *,'Tangent:',Tangent
          PRINT *,'Tangent2:',Tangent2
        END IF

        ! Move to local normal-tangential coordinate system for the slave element
        DO i=1,n        
          r(1) = Nodes % x(i)
          r(2) = Nodes % y(i)
          r(3) = Nodes % z(i)
      
          ! Coordinate projected to nt-coordinates
          Nodes % x(i) = SUM( Tangent * r ) 
          Nodes % y(i) = SUM( Tangent2 * r ) 
          Nodes % z(i) = SUM( Normal * r ) 
        END DO
        
        ! Even for quadratic elements only work with corner nodes (n >= ne)        
        xmin = MINVAL(Nodes % x(1:ne))
        xmax = MAXVAL(Nodes % x(1:ne))

        ymin = MINVAL(Nodes % y(1:ne))
        ymax = MAXVAL(Nodes % y(1:ne))

        zmin = MINVAL( Nodes % z(1:ne))
        zmax = MAXVAL( Nodes % z(1:ne))
        zave = SUM( Nodes % z(1:ne) ) / ne 
        
        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;

        IF( DebugElem ) THEN
          PRINT *,'Element n-t range:'
          PRINT *,'xrange:',xmin,xmax
          PRINT *,'yrange:',ymin,ymax
          PRINT *,'zrange:',zmin,zmax
        END IF

        ! Nullify z since we don't need it anymore after registering (zmin,zmax)
        Nodes % z = 0.0_dp
        
        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        
        IP = GaussPoints( Element, PreferenceElement = .FALSE. )
        RefArea = detJ * SUM( IP % s(1:IP % n) ) 
        SumArea = 0.0_dp

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_a.dat'
          OPEN( 10,FILE=Filename)
          DO i=1,ne
            WRITE( 10, * ) Nodes % x(i), Nodes % y(i), Nodes % z(i)
          END DO
          CLOSE( 10 )
        END IF
        
        DO i=1,n
          j = InvPerm1(Indexes(i))
          nrow = NodePerm(j)
          IF( nrow == 0 ) CYCLE
          CALL List_AddMatrixIndex(Projector % ListMatrix, nrow, j ) 
          IF(ASSOCIATED(Projector % Child)) &
              CALL List_AddMatrixIndex(Projector % Child % ListMatrix, nrow, j ) 
        END DO

        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        ElemCands = 0
        ElemHits = 0
        SubTri = 0
        
        DO indM=1,BMesh2 % NumberOfBulkElements

         ! Rough search, note that this cannot be too tight since then
          ! we loose also the contacts.
          IF( ABS( Center(1) - Center2 % x(indM) ) > MaxDistance ) CYCLE
          IF( ABS( Center(2) - Center2 % y(indM) ) > MaxDistance ) CYCLE
          IF( ABS( Center(3) - Center2 % z(indM) ) > MaxDistance ) CYCLE

          IF( DebugElem ) THEN
            PRINT *,'Candidate Elem Center:',indM,Center2 % x(indM),&
                Center2 % y(indM),Center2 % z(indM)           
          END IF
          
          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes

          nM = ElementM % TYPE % NumberOfNodes
          neM = ElementM % TYPE % ElementCode / 100
            
          DO i=1,nM
            j = IndexesM(i)
            r(1) = BMesh2 % Nodes % x(j)
            r(2) = BMesh2 % Nodes % y(j)
            r(3) = BMesh2 % Nodes % z(j)
            
            ! Coordinate projected to nt-coordinates
            NodesM % x(i) = SUM( Tangent * r ) 
            NodesM % y(i) = SUM( Tangent2 * r ) 
            NodesM % z(i) = SUM( Normal * r ) 
          END DO
                    
          ! Now we can make the 2nd quick search in the nt-system.
          ! Now the tangential coordinates can be treated exactly.
          xminm = MINVAL( NodesM % x(1:neM) )
          IF( xminm > xmax ) CYCLE

          xmaxm = MAXVAL( NodesM % x(1:neM) )
          IF( xmaxm < xmin ) CYCLE

          yminm = MINVAL( NodesM % y(1:neM))
          IF( yminm > ymax ) CYCLE
          
          ymaxm = MAXVAL( NodesM % y(1:neM))
          IF( ymaxm < ymin ) CYCLE

          zminm = MINVAL( NodesM % z(1:neM) )
          IF( zminm > zmax + MaxDistance ) CYCLE

          zmaxm = MAXVAL( NodesM % z(1:neM) )
          IF( zmaxm < zmin - MaxDistance ) CYCLE

          NormalM = NormalVector( ElementM, NodesM, Check = .FALSE. ) 
          IF( BMesh2 % PeriodicFlip(indM) ) NormalM = -NormalM

          IF( DebugElem ) THEN
            PRINT *,'ElementM n-t range:'
            PRINT *,'xrange:',xminm,xmaxm
            PRINT *,'yrange:',yminm,ymaxm
            PRINT *,'zrange:',zminm,zmaxm
            PRINT *,'Candidate elem normal:',NormalM, BMesh2 % PeriodicFlip(indM)
          END IF
                    
          ! We must compare this normal to the nt-system where the slave normal is (0,0,1)
          ! Positive normal means that this element is pointing to the same direction!
          IF( NormalM(3) >= MaxNormalDot ) THEN
            IF( DebugElem ) PRINT *,'Normals are not facing!' 
            CYCLE
          END IF
            
          ! Nullify z since we don't need it anymore 
          NodesM % z = 0.0_dp
                  
          k = 0
          ElemCands = ElemCands + 1
          CornerFound = .FALSE.
          CornerFoundM = .FALSE.

          ! Check through the nodes that are created in the intersections of any two edge
          DO i=1,ne
            x1 = Nodes % x(i)
            y1 = Nodes % y(i)
            i2 = i + 1 
            IF( i2 > ne ) i2 = 1  ! check the (ne,1) edge also
            x2 = Nodes % x(i2)
            y2 = Nodes % y(i2)

            DO iM=1,neM
              x1M = NodesM % x(iM)
              y1M = NodesM % y(iM)
              i2M = iM + 1
              IF( i2M > neM ) i2M = 1
              x2M = NodesM % x(i2M)
              y2M = NodesM % y(i2M)

              ! Upon solution this is tampered so it must be initialized 
              ! before each solution. 
              A(1,1) = x2 - x1
              A(2,1) = y2 - y1           
              A(1,2) = x1M - x2M
              A(2,2) = y1M - y2M

              detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
              absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))

              ! Lines are almost parallel => no intersection possible
              ! Check the dist at the end of the line segments.
              IF(ABS(detA) < 1.0d-8 * absA + 1.0d-20 ) CYCLE

              B(1) = x1M - x1
              B(2) = y1M - y1

              CALL InvertMatrix( A,2 )
              C(1:2) = MATMUL(A(1:2,1:2),B(1:2))

              ! Check that the hit is within the line segment
              IF(ANY(C(1:2) < 0.0) .OR. ANY(C(1:2) > 1.0d0)) CYCLE

              ! We have a hit, two line segments can have only one hit
              k = k + 1

              x(k) = x1 + C(1) * (x2-x1)
              y(k) = y1 + C(1) * (y2-y1)

              ! If the point of intersection is at the end of a line-segment it
              ! is also a corner node.
              IF(ABS(C(1)) < 1.0d-6 ) THEN
                CornerFound(i) = .TRUE.
              ELSE IF( ABS(C(1)-1.0_dp ) < 1.0d-6 ) THEN
                CornerFound(i2) = .TRUE.
              END IF

              IF(ABS(C(2)) < 1.0d-6 ) THEN
                CornerFoundM(iM) = .TRUE.
              ELSE IF( ABS(C(2)-1.0_dp ) < 1.0d-6 ) THEN
                CornerFoundM(i2M) = .TRUE.
              END IF

              EdgeHits = EdgeHits + 1
            END DO
          END DO

          IF( DebugElem ) THEN
            PRINT *,'EdgeHits:',k,COUNT(CornerFound),COUNT(CornerFoundM)
          END IF

          ! Check the nodes that are one of the existing nodes i.e. corner nodes
          ! that are located inside in either element. We have to check both combinations. 
          DO i=1,ne
            ! This corner was already determined active as the end of edge 
            IF( CornerFound(i) ) CYCLE

            Point(1) = Nodes % x(i)
            IF( Point(1) < xminm - tolS ) CYCLE
            IF( Point(1) > xmaxm + tolS ) CYCLE

            Point(2) = Nodes % y(i)
            IF( Point(2) < yminm - TolS ) CYCLE
            IF( Point(2) > ymaxm + TolS ) CYCLE

            ! The edge intersections should catch the sharp hits so here we can use hard criteria
            Found = PointInElement( ElementM, NodesM, Point, uvw, LocalEps = 1.0d-8 )
            IF( Found ) THEN
              k = k + 1
              x(k) = Point(1)
              y(k) = Point(2)
              CornerHits = CornerHits + 1
            END IF
          END DO

                    
          ! Possible corner hits for the master element
          DO i=1,neM
            IF( CornerFoundM(i) ) CYCLE

            Point(1) = NodesM % x(i)
            IF( Point(1) < xmin - tols ) CYCLE
            IF( Point(1) > xmax + tols ) CYCLE

            Point(2) = NodesM % y(i)
            IF( Point(2) < ymin - Tols ) CYCLE
            IF( Point(2) > ymax + Tols ) CYCLE

            Found = PointInElement( Element, Nodes, Point, uvw, LocalEps = 1.0d-8 )
            IF( Found ) THEN
              k = k + 1
              x(k) = Point(1)
              y(k) = Point(2)
              CornerHits = CornerHits + 1
            END IF
          END DO

          IF( DebugElem ) THEN
            PRINT *,'Total and corner hits:',k,CornerHits
          END IF
          
          kmax = k          
          IF( kmax < 3 ) CYCLE

          sgn0 = 1

          InitialHits = InitialHits + kmax

          ! The polygon is convex and hence its center lies inside the polygon
          xt = SUM(x(1:kmax)) / kmax
          yt = SUM(y(1:kmax)) / kmax
            
          ! Set the angle from the center and order the nodes so that they 
          ! can be easily triangulated.
          DO k=1,kmax
            phi(k) = ATAN2( y(k)-yt, x(k)-xt )
            inds(k) = k
          END DO

          IF( DebugElem ) THEN            
            PRINT *,'Polygon Coords:',k
            PRINT *,'x:',x(1:k)
            PRINT *,'y:',y(1:k)
            PRINT *,'PolygonArea:',(MAXVAL(x(1:k))-MINVAL(x(1:k)))*(MAXVAL(y(1:k))-MINVAL(y(1:k)))
            PRINT *,'Center:',xt,yt
            PRINT *,'Phi:',phi(1:kmax)
          END IF

          CALL SortR(kmax,inds,phi)
          
          x(1:kmax) = x(inds(1:kmax))
          y(1:kmax) = y(inds(1:kmax))
          
          IF( DebugElem ) THEN
            PRINT *,'Sorted Inds:',inds(1:kmax)
            PRINT *,'Sorted Phi:',phi(1:kmax)
          END IF
   
          ! Eliminate redundant corners from the polygon
          j = 1
          DO k=2,kmax
            dist = (x(j)-x(k))**2 + (y(j)-y(k))**2 
            IF( dist > Tols ) THEN
              j = j + 1
              IF( j /= k ) THEN
                x(j) = x(k)
                y(j) = y(k)
              END IF
            END IF
          END DO
          
          IF( DebugElem ) THEN
            IF( kmax > j ) PRINT *,'Corners reduced to:',j
          END IF
           
          kmax = j
          IF( kmax < 3 ) CYCLE

          ElemHits = ElemHits + 1
          ActiveHits = ActiveHits + kmax

          IF( kmax > MaxSubTriangles ) THEN
            MaxSubTriangles = kmax
            MaxSubElem = ind
          END IF

          IF( SaveElem ) THEN
            FileName = 't'//TRIM(I2S(TimeStep))//'_b'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) NodesM % x(i), NodesM % y(i)
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_c'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            WRITE( 10, * ) xt, yt
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_e'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,kmax
              WRITE( 10, * ) x(i), y(i)
            END DO
            CLOSE( 10 )           
          END IF

          Depth = zave - SUM( NodesM % z(1:neM) )/neM 
          MaxDepth = MAX( Depth, MaxDepth )
          MinDepth = MIN( Depth, MinDepth ) 
          
          ! Deal the case with multiple corners by making 
          ! triangulariation using one corner point.
          ! This should be ok as the polygon is always convex.
          NodesT % x(1) = x(1)
          NodesT % y(1) = y(1)

          DO k=1,kmax-2                         

            ! This check over area also automatically elimiates redundant nodes
            ! that were detected twice.
            dArea = 0.5_dp*ABS( (x(k+1)-x(1))*(y(k+2)-y(1)) -(x(k+2)-x(1))*(y(k+1)-y(1)))
            
            IF( dArea < TolS**2 * RefArea ) CYCLE

            ! Triangle is created by keeping one corner node fixed and rotating through
            ! the other nodes. 
            NodesT % x(2) = x(k+1)
            NodesT % y(2) = y(k+1)
            NodesT % x(3) = x(k+2)
            NodesT % y(3) = y(k+2)

            IF( DebugElem ) THEN
              PRINT *,'Temporal element n-t coordinates',k
              PRINT *,'x:',NodesT % x
              PRINT *,'y:',NodesT % y
            END IF
            
            IF( SaveElem ) THEN
              SubTri = SubTri + 1
              FileName = 't'//TRIM(I2S(TimeStep))//'_s'//TRIM(I2S(SubTri))//'.dat'
              OPEN( 10,FILE=FileName)
              DO i=1,3
                WRITE( 10, * ) NodesT % x(i), NodesT % y(i)
              END DO
              CLOSE( 10 )
            END IF
            
            CALL TemporalTriangleMortarAssembly(ElementT, NodesT, Element, Nodes, ElementM, NodesM, &
                BiorthogonalBasis, DualMaster, DualLCoeff, NoGaussPoints, Projector, NodeScale, &
                NodePerm, InvPerm1, InvPerm2, SumArea ) 
          END DO
                             
          IF( DebugElem ) PRINT *,'Element integrated:',indM,SumArea,RefArea,SumArea / RefArea

          ! If we have integrated enough area we are done!
          IF( SumArea > RefArea*(1.0_dp - 1.0e-6) ) EXIT

        END DO ! indM

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_n.dat'
          OPEN( 10,FILE=Filename)
          OPEN( 10,FILE=FileName)
          WRITE( 10, * ) ElemHits 
          CLOSE( 10 )
        END IF

        TotCands = TotCands + ElemCands
        TotHits = TotHits + ElemHits
        TotSumArea = TotSumArea + SumArea
        TotRefArea = TotRefArea + RefArea

        Err = SumArea / RefArea
        IF( Err > MaxErr ) THEN
          MaxErr = Err
          MaxErrInd = Err
        END IF
        IF( Err < MinErr ) THEN
          MinErr = Err
          MinErrInd = ind
        END IF

        IF( SaveErr ) THEN
          WRITE( 11, * ) ind,SUM( Nodes % x(1:ne))/ne, SUM( Nodes % y(1:ne))/ne, Err
        END IF

      END DO

      IF( SaveErr ) CLOSE(11)
      
        
      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, &
          NodesM % x, NodesM % y, NodesM % z, &
          NodesT % x, NodesT % y, NodesT % z, &
          Center2 % x, Center2 % y, Center2 % z, Basis )
       
      CALL Info(Caller,'Number of integration pair candidates: '&
          //TRIM(I2S(TotCands)),Level=10)
      CALL Info(Caller,'Number of integration pairs: '&
          //TRIM(I2S(TotHits)),Level=10)

      CALL Info(Caller,'Number of edge intersections: '&
          //TRIM(I2S(EdgeHits)),Level=10)
      CALL Info(Caller,'Number of corners inside element: '&
          //TRIM(I2S(EdgeHits)),Level=10)

      CALL Info(Caller,'Number of initial corners: '&
          //TRIM(I2S(InitialHits)),Level=10)
      CALL Info(Caller,'Number of active corners: '&
          //TRIM(I2S(ActiveHits)),Level=10)

      CALL Info(Caller,'Number of most subelement corners: '&
          //TRIM(I2S(MaxSubTriangles)),Level=10)
      CALL Info(Caller,'Element of most subelement corners: '&
          //TRIM(I2S(MaxSubElem)),Level=10)

      WRITE( Message,'(A,ES12.5)') 'Total reference area:',TotRefArea
      CALL Info(Caller,Message,Level=8)
      WRITE( Message,'(A,ES12.5)') 'Total integrated area:',TotSumArea
      CALL Info(Caller,Message,Level=8)

      Err = TotSumArea / TotRefArea
      WRITE( Message,'(A,ES15.6)') 'Average ratio in area integration:',Err 
      CALL Info(Caller,Message,Level=5)

      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Maximum relative discrepancy in areas (element: ',MaxErrInd,'):',MaxErr-1.0_dp 
      CALL Info(Caller,Message,Level=6)
      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Minimum relative discrepancy in areas (element: ',MinErrInd,'):',MinErr-1.0_dp 
      CALL Info(Caller,Message,Level=6)

      WRITE( Message,'(A,ES12.4)') &
          'Minimum depth in normal direction:',MinDepth
      CALL Info(Caller,Message,Level=8)
      WRITE( Message,'(A,ES12.4)') &
          'Maximum depth in normal direction:',MaxDepth
      CALL Info(Caller,Message,Level=8)
      
    END SUBROUTINE NormalProjectorWeak3D

  END FUNCTION NormalProjector

  

  !---------------------------------------------------------------------------
  !> Create a projector for mapping between interfaces using the Galerkin method
  !> A temporal mesh structure with a node for each Gaussian integration point is 
  !> created. Then this projector matrix is transferred to a projector on the nodal
  !> coordinates.   
  !---------------------------------------------------------------------------
   FUNCTION NodalProjector(BMesh2, BMesh1, &
       UseQuadrantTree, Repeating, AntiRepeating ) &
      RESULT ( Projector )
  !---------------------------------------------------------------------------
    USE Lists

    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    LOGICAL :: UseQuadrantTree, Repeating, AntiRepeating
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    INTEGER :: i,j,k,n
    INTEGER, POINTER :: Rows(:),Cols(:)
    REAL(KIND=dp), POINTER :: Values(:)

    BMesh1 % Parent => NULL()
    BMesh2 % Parent => NULL()

    InvPerm1 => BMesh1 % InvPerm
    InvPerm2 => BMesh2 % InvPerm

    ! Set the nodes of Mesh1 to be in the interval defined by Mesh2
    !-----------------------------------------------------------------
    IF( Repeating ) THEN
      IF( AntiRepeating ) THEN
        ALLOCATE( MirrorNode( BMesh1 % NumberOfNodes ) )
        MirrorNode = .FALSE.
      END IF
      CALL PreRotationalProjector(BMesh1, BMesh2, MirrorNode )
    END IF

    ! Create the projector using nodal points 
    ! This corresponds to numerical integration of the collocation method.
    !-----------------------------------------------------------------
    Projector => MeshProjector( BMesh2, BMesh1, UseQuadrantTree )    
    Projector % ProjectorType = PROJECTOR_TYPE_NODAL

    Values => Projector % Values
    Cols => Projector % Cols
    Rows => Projector % Rows

    ! One needs to change the sign of the projector for the mirror nodes
    !-----------------------------------------------------------------------------
    IF( Repeating .AND. AntiRepeating ) THEN
      CALL PostRotationalProjector( Projector, MirrorNode )
      DEALLOCATE( MirrorNode ) 
    END IF

    ! Now return from the indexes of the interface mesh system to the 
    ! original mesh system.
    !-----------------------------------------------------------------
    n = SIZE( InvPerm1 ) 
    ALLOCATE( Projector % InvPerm(n) )
    Projector % InvPerm = InvPerm1

    DO i=1,Projector % NumberOfRows
       DO j = Rows(i), Rows(i+1)-1
         k = Cols(j)    
         IF ( k > 0 ) Cols(j) = InvPerm2(k)
       END DO
    END DO

  END FUNCTION NodalProjector
!------------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !> Create a nodal projector related to discontinuous interface.
  !---------------------------------------------------------------------------
   FUNCTION NodalProjectorDiscont( Mesh, bc ) RESULT ( Projector )
  !---------------------------------------------------------------------------
    USE Lists

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: bc
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    TYPE(Model_t), POINTER :: Model
    INTEGER, POINTER :: NodePerm(:)
    INTEGER :: i,j,n,m
    INTEGER, POINTER :: Rows(:),Cols(:), InvPerm(:)
    REAL(KIND=dp), POINTER :: Values(:)
    LOGICAL :: Found
    CHARACTER(*), PARAMETER :: Caller='NodalProjectorDiscont'

    
    CALL Info(Caller,'Creating nodal projector for discontinuous boundary',Level=7)

    Projector => Null()
    IF( .NOT. Mesh % DisContMesh ) THEN
      CALL Warn(Caller,'Discontinuous mesh not created?')
      RETURN
    END IF

    Model => CurrentModel
    j = 0
    DO i=1,Model % NumberOfBCs
      IF( ListGetLogical(Model % BCs(i) % Values,'Discontinuous Boundary',Found) ) THEN
        j = j + 1
      END IF
    END DO
    ! This is a temporal limitations
    IF( j > 1 ) THEN
      CALL Warn(Caller,'One BC (not '&
          //TRIM(I2S(j))//') only for discontinuous boundary!')
    END IF


    NodePerm => Mesh % DisContPerm
    n = SIZE( NodePerm ) 
    m = COUNT( NodePerm > 0 ) 

    Projector => AllocateMatrix()
    Projector % ProjectorType = PROJECTOR_TYPE_NODAL
    Projector % ProjectorBC = bc

    ALLOCATE( Projector % Cols(m) )
    ALLOCATE( Projector % Values(m) )
    ALLOCATE( Projector % Rows(m+1) )
    ALLOCATE( Projector % InvPerm(m) )

    Cols => Projector % Cols
    Values => Projector % Values
    Rows => Projector % Rows
    InvPerm => Projector % InvPerm
    Projector % NumberOfRows = m

    Values = 1.0_dp
    DO i=1,m+1
      Rows(i) = i
    END DO

    DO i=1,n
      j = NodePerm(i)
      IF( j == 0 ) CYCLE
      Cols(j) = n + j
      InvPerm(j) = i
    END DO

  END FUNCTION NodalProjectorDiscont
!------------------------------------------------------------------------------

  
  !---------------------------------------------------------------------------------
  ! Create a permutation to eliminate edges in a conforming case.
  !---------------------------------------------------------------------------------
  SUBROUTINE ConformingEdgePerm( Mesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
    TYPE(Mesh_t), POINTER :: Mesh, BMesh1, BMesh2
    INTEGER, POINTER :: PerPerm(:)
    LOGICAL, POINTER :: PerFlip(:)
    LOGICAL, OPTIONAL :: AntiPeriodic 
    !---------------------------------------------------------------------------------      
    INTEGER :: n, ind, indm, e, em, eind, eindm, k1, k2, km1, km2, sgn0, sgn, i1, i2, &
        noedges, noedgesm, Nundefined, n0
    TYPE(Element_t), POINTER :: Edge, EdgeM
    INTEGER, POINTER :: Indexes(:), IndexesM(:)
    REAL(KIND=dp) :: xm1, xm2, ym1, ym2, x1, y1, x2, y2, y2m, nrow
    INTEGER, ALLOCATABLE :: PeriodicEdge(:), EdgeInds(:), EdgeIndsM(:)
    REAL(KIND=dp), ALLOCATABLE :: EdgeX(:,:), EdgeY(:,:), EdgeMX(:,:), EdgeMY(:,:)
    REAL(KIND=dp) :: coordprod, indexprod, ss, minss, maxminss
    INTEGER :: minuscount, samecount, mini, doubleusecount
    LOGICAL :: Parallel, AntiPer
    LOGICAL, ALLOCATABLE :: EdgeUsed(:)
    CHARACTER(*), PARAMETER :: Caller='ConformingEdgePerm'
    
    
    CALL Info(Caller,'Creating permutation for elimination of conforming edges',Level=8)

    n = Mesh % NumberOfEdges
    IF( n == 0 ) RETURN

    AntiPer = .FALSE.
    IF( PRESENT( AntiPeriodic ) ) AntiPer = AntiPeriodic

    CALL CreateEdgeCenters( Mesh, BMesh1, noedges, EdgeInds, EdgeX, EdgeY ) 
    CALL Info(Caller,'Number of edges in slave mesh: '//TRIM(I2S(noedges)),Level=10)

    CALL CreateEdgeCenters( Mesh, BMesh2, noedgesm, EdgeIndsM, EdgeMX, EdgeMY )
    CALL Info(Caller,'Number of edges in master mesh: '//TRIM(I2S(noedgesm)),Level=10)

    IF( noedges == 0 ) RETURN
    IF( noedgesm == 0 ) RETURN
    
    ALLOCATE( PeriodicEdge(noedges),EdgeUsed(noedgesm))
    PeriodicEdge = 0
    EdgeUsed = .FALSE.
    maxminss = 0.0_dp
    n0 = Mesh % NumberOfNodes
    Parallel = ( ParEnv % PEs > 1 )
    samecount = 0
    doubleusecount = 0
    
    DO i1=1,noedges
      x1 = EdgeX(3,i1)
      y1 = EdgeY(3,i1)

      IF( PerPerm( EdgeInds(i1) + n0 ) > 0 ) CYCLE

      minss = HUGE(minss)
      mini = 0

      DO i2=1,noedgesm
        x2 = EdgeMX(3,i2)
        y2 = EdgeMY(3,i2)

        ss = (x1-x2)**2 + (y1-y2)**2
        IF( ss < minss ) THEN
          minss = ss
          mini = i2
        END IF
      END DO

      IF( EdgeInds(i1) == EdgeIndsM(mini) ) THEN        
        samecount = samecount + 1        
        CYCLE
      END IF

      IF( EdgeUsed(mini ) ) THEN
        doubleusecount = doubleusecount + 1
      ELSE
        EdgeUsed(mini) = .TRUE.
      END IF
              
      ! we have a hit
      PeriodicEdge(i1) = mini
      maxminss = MAX( maxminss, minss )
    END DO

    WRITE(Message,'(A,ES12.4)') 'Maximum minimum deviation in edge centers:',SQRT(maxminss)
    CALL Info(Caller,Message,Level=8)

    minuscount = 0

    DO e=1,noedges        
      eind = EdgeInds(e)

      ! This has already been set
      IF( PerPerm(eind+n0) > 0 ) CYCLE

      ! Get the conforming counterpart
      em = PeriodicEdge(e)
      IF( em == 0 ) CYCLE
      eindm = EdgeIndsM(em)        

      ! Get the coordinates and indexes of the 1st edge
      Edge => Mesh % Edges(eind)
      k1 = Edge % NodeIndexes( 1 )
      k2 = Edge % NodeIndexes( 2 )
      IF(Parallel) THEN
        k1 = Mesh % ParallelInfo % GlobalDOFs(k1) !BMesh1 % InvPerm(k1))
        k2 = Mesh % ParallelInfo % GlobalDOFs(k2) !BMesh1 % InvPerm(k2))
      END IF

      ! We cannot use the (x,y) coordinates of the full "Mesh" as the boundary meshes
      ! have been mapped such that interpolation is possible. 
      x1 = EdgeX(1,e)
      x2 = EdgeX(2,e)
      y1 = EdgeY(1,e)
      y2 = EdgeY(2,e)

      ! Get the coordinates and indexes of the 2nd edge
      EdgeM => Mesh % Edges(eindm)
      km1 = EdgeM % NodeIndexes( 1 )
      km2 = EdgeM % NodeIndexes( 2 )
      IF(Parallel) THEN
        km1 = Mesh % ParallelInfo % GlobalDOFs(km1) !BMesh2 % InvPerm(km1))
        km2 = Mesh % ParallelInfo % GlobalDOFs(km2) !BMesh2 % InvPerm(km2))
      END IF
      
      xm1 = EdgeMX(1,em)
      xm2 = EdgeMX(2,em)
      ym1 = EdgeMY(1,em)
      ym2 = EdgeMY(2,em)

      coordprod = (x1-x2)*(xm1-xm2) + (y1-y2)*(ym1-ym2) 
      indexprod = (k1-k2)*(km1-km2)

      IF( coordprod * indexprod < 0 ) THEN
        minuscount = minuscount + 1
        PerFlip(eind+n0) = .NOT. AntiPer
        !PRINT *,'prod:',coordprod,indexprod
        !PRINT *,'x:',x1,x2,xm1,xm2
        !PRINT *,'y:',y1,y2,ym1,ym2
        !PRINT *,'k:',k1,k2,km1,km2
      ELSE
        PerFlip(eind+n0) = AntiPer
      END IF

      ! Mark that this is set so it don't need to be set again
      PerPerm(eind+n0) = eindm + n0
    END DO

    DEALLOCATE( EdgeInds, EdgeX, EdgeY ) 
    DEALLOCATE( EdgeIndsM, EdgeMX, EdgeMY )
    DEALLOCATE( PeriodicEdge )

    IF( samecount > 0 ) THEN
      CALL Info(Caller,'Number of edges are the same: '//TRIM(I2S(samecount)),Level=8)
    END IF
        
    IF( minuscount == 0 ) THEN
      CALL Info(Caller,'All edges in conforming projector have consistent sign!',Level=8)
    ELSE
      CALL Info(Caller,'Flipped sign of '//TRIM(I2S(minuscount))//&
          ' (out of '//TRIM(I2S(noedges))//') edge projectors',Level=6)
    END IF

    IF( doubleusecount > 0 ) THEN
      CALL Fatal(Caller,'This is not conforming! Number of edges used twice: '//TRIM(I2S(doubleusecount)))
    END IF

    
  CONTAINS 
    
    ! Create edge centers for the mapping routines.
    !------------------------------------------------------------------------------
    SUBROUTINE CreateEdgeCenters( Mesh, EdgeMesh, noedges, EdgeInds, EdgeX, EdgeY ) 

      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Mesh_t), POINTER :: EdgeMesh
      INTEGER :: noedges
      INTEGER, ALLOCATABLE :: EdgeInds(:)
      REAL(KIND=dp), ALLOCATABLE :: EdgeX(:,:), EdgeY(:,:)

      LOGICAL, ALLOCATABLE :: EdgeDone(:)
      INTEGER :: ind, eind, i, i1, i2, k1, k2, ktmp
      TYPE(Element_t), POINTER :: Element
      INTEGER, POINTER :: EdgeMap(:,:), Indexes(:)
      LOGICAL :: AllocationsDone 


      ALLOCATE( EdgeDone( Mesh % NumberOfEdges ) )
      AllocationsDone = .FALSE.


100   noedges = 0
      EdgeDone = .FALSE.

      DO ind=1,EdgeMesh % NumberOfBulkElements

        Element => EdgeMesh % Elements(ind)        
        EdgeMap => GetEdgeMap( Element % TYPE % ElementCode / 100)

        Indexes => Element % NodeIndexes

        DO i = 1,Element % TYPE % NumberOfEdges          

          eind = Element % EdgeIndexes(i)

          IF( EdgeDone(eind) ) CYCLE

          noedges = noedges + 1
          EdgeDone(eind) = .TRUE.

          IF( ALLOCATED( EdgeInds ) ) THEN            
            ! Get the nodes of the edge
            i1 = EdgeMap(i,1) 
            i2 = EdgeMap(i,2)

            ! These point to the local boundary mesh
            k1 = Indexes( i1 )
            k2 = Indexes( i2 )

            ! Ensure that the order of node is consistent with the global mesh
            ! because this is later used to check the sign of the edge. 
            IF( EdgeMesh % InvPerm(k1) /= Mesh % Edges(eind) % NodeIndexes(1) ) THEN
              IF( EdgeMesh % InvPerm(k1) /= Mesh % Edges(eind) % NodeIndexes(2) ) THEN
                PRINT *,'We have a problem with the edges:',k1,k2
              END IF
              ktmp = k1
              k1 = k2
              k2 = ktmp
            END IF

            EdgeX(1,noedges) = EdgeMesh % Nodes % x(k1)
            EdgeX(2,noedges) = EdgeMesh % Nodes % x(k2)

            EdgeY(1,noedges) = EdgeMesh % Nodes % y(k1)
            EdgeY(2,noedges) = EdgeMesh % Nodes % y(k2)

            ! The center of the edge (note we skip multiplication by 0.5 is it is redundant)
            EdgeX(3,noedges) = EdgeX(1,noedges) + EdgeX(2,noedges)
            EdgeY(3,noedges) = EdgeY(1,noedges) + EdgeY(2,noedges)

            EdgeInds(noedges) = eind
          END IF
        END DO
      END DO

      IF(noedges > 0 .AND. .NOT. AllocationsDone ) THEN
        CALL Info(Caller,'Allocating stuff for edges',Level=20)
        ALLOCATE( EdgeInds(noedges), EdgeX(3,noedges), EdgeY(3,noedges) )
        AllocationsDone = .TRUE.
        GOTO 100
      END IF

      DEALLOCATE( EdgeDone ) 

    END SUBROUTINE CreateEdgeCenters

    
  END SUBROUTINE ConformingEdgePerm



  ! Create a permutation to eliminate nodes in a conforming case.
  !----------------------------------------------------------------------
  SUBROUTINE ConformingNodePerm( Mesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
    TYPE(Mesh_t), POINTER :: Mesh, BMesh1, BMesh2
    INTEGER, POINTER :: PerPerm(:)
    LOGICAL, POINTER, OPTIONAL :: PerFlip(:)
    LOGICAL, OPTIONAL :: AntiPeriodic 
    !----------------------------------------------------------------------
    INTEGER :: n, i1, i2, j1, j2, k1, k2, mini, samecount, doubleusecount
    REAL(KIND=dp) :: x1, y1, z1, x2, y2, z2
    REAL(KIND=dp) :: ss, minss, maxminss
    LOGICAL, ALLOCATABLE :: NodeUsed(:)
    CHARACTER(*), PARAMETER :: Caller='ConformingNodePerm'

    
    CALL Info(Caller,'Creating permutations for conforming nodes',Level=8)

    n = 0
    IF( PRESENT( PerFlip ) ) n = n + 1
    IF( PRESENT( AntiPeriodic ) ) n = n + 1
    IF( n == 1 ) THEN
      CALL Fatal(Caller,'Either have zero or two optional parameters!')
    END IF
      
    n = Mesh % NumberOfNodes
    IF( n == 0 ) RETURN      

    IF( Bmesh1 % NumberOfNodes == 0 ) RETURN
    IF( Bmesh2 % NumberOfNodes == 0 ) RETURN

    maxminss = 0.0_dp
    samecount = 0
    doubleusecount = 0

    ALLOCATE( NodeUsed(BMesh2 % NumberOfNodes) )
    NodeUsed = .FALSE.
    
    DO i1=1,Bmesh1 % NumberOfNodes

      j1 = BMesh1 % InvPerm(i1)
      IF( PerPerm(j1) > 0 ) CYCLE

      x1 = BMesh1 % Nodes % x(i1)
      y1 = BMesh1 % Nodes % y(i1)      
      z1 = BMesh1 % Nodes % z(i1)      

      minss = HUGE(minss)
      mini = 0

      DO i2=1,Bmesh2 % NumberOfNodes
        x2 = BMesh2 % Nodes % x(i2)
        y2 = BMesh2 % Nodes % y(i2)
        z2 = BMesh2 % Nodes % z(i2)

        ss = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
        IF( ss < minss ) THEN
          minss = ss
          mini = i2
        END IF

        ! This should be a hit even in conservative terms.
        IF( minss < EPSILON( minss ) ) EXIT
      END DO

      ! Assume that the closest node is a hit
      IF( j1 == BMesh2 % InvPerm(mini) ) THEN
        samecount = samecount + 1
        CYCLE
      END IF

      IF( NodeUsed(mini ) ) THEN
        doubleusecount = doubleusecount + 1
      ELSE
        NodeUsed(mini) = .TRUE.
      END IF
        
      PerPerm(j1) = BMesh2 % InvPerm(mini)

      maxminss = MAX( maxminss, minss )

      IF( PRESENT( PerFlip ) ) THEN
        IF( AntiPeriodic ) PerFlip(j1) = .TRUE.
      END IF
    END DO

    IF( samecount > 0 ) THEN
      CALL Info(Caller,'Number of nodes are the same: '//TRIM(I2S(samecount)),Level=8)
    END IF

    WRITE(Message,'(A,ES12.4)') 'Maximum minimum deviation in node coords:',SQRT(maxminss)
    CALL Info(Caller,Message,Level=8)

    IF( doubleusecount > 0 ) THEN
      CALL Fatal(Caller,'This is not conforming! Number of nodes used twice: '//TRIM(I2S(doubleusecount)))
    END IF

  END SUBROUTINE ConformingNodePerm
  !----------------------------------------------------------------------


   
  !---------------------------------------------------------------------------
  !> Create a projector for mixed nodal / edge problems assuming constant level
  !> in the 2nd direction. This kind of projector is suitable for 2D meshes where
  !> the mortar line is effectively 1D, or to 3D cases that have been created by
  !> extrusion. 
  !---------------------------------------------------------------------------
  FUNCTION LevelProjector( BMesh1, BMesh2, Repeating, AntiRepeating, &
      FullCircle, Radius, DoNodes, DoEdges, NodeScale, EdgeScale, BC ) &
      RESULT ( Projector )
    !---------------------------------------------------------------------------
    USE Lists
    USE Messages
    USE Types
    USE GeneralUtils
    IMPLICIT NONE

    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2, Mesh
    LOGICAL :: DoNodes, DoEdges
    LOGICAL :: Repeating, AntiRepeating, FullCircle, NotAllQuads, NotAllQuads2
    REAL(KIND=dp) :: Radius, NodeScale, EdgeScale
    TYPE(ValueList_t), POINTER :: BC
    TYPE(Matrix_t), POINTER :: Projector    
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
    LOGICAL ::  StrongNodes, StrongEdges, StrongLevelEdges, StrongExtrudedEdges, &
        StrongSkewEdges, StrongConformingEdges, StrongConformingNodes
    LOGICAL :: Found, Parallel, SelfProject, EliminateUnneeded, SomethingUndone, &
        EdgeBasis, PiolaVersion, GenericIntegrator, Rotational, Cylindrical, WeakProjector, &
        StrongProjector, CreateDual, HaveMaxDistance
    REAL(KIND=dp) :: XmaxAll, XminAll, YminAll, YmaxAll, Xrange, Yrange, &
        RelTolX, RelTolY, XTol, YTol, RadTol, MaxSkew1, MaxSkew2, SkewTol, &
        ArcCoeff, EdgeCoeff, NodeCoeff, MaxDistance
    INTEGER :: NoNodes1, NoNodes2, MeshDim
    INTEGER :: i,j,k,n,m,Nrange,Nrange2, nrow, Naxial
    INTEGER, ALLOCATABLE :: EdgePerm(:),NodePerm(:),DualNodePerm(:)
    INTEGER :: EdgeRow0, FaceRow0, EdgeCol0, FaceCol0, ProjectorRows
    TYPE(Element_t), POINTER :: Element
    INTEGER, POINTER :: NodeIndexes(:)
    REAL(KIND=dp), ALLOCATABLE :: Cond(:)
    TYPE(Matrix_t), POINTER :: DualProjector    
    LOGICAL :: DualMaster, DualSlave, DualLCoeff, BiorthogonalBasis
    LOGICAL :: SecondOrder
    CHARACTER(*), PARAMETER :: Caller='LevelProjector'

    CALL Info(Caller,'Creating projector for a levelized mesh',Level=7)

    IF(.NOT. (DoEdges .OR. DoNodes ) ) THEN
      CALL Warn(Caller,'Nothing to do, no nonodes, no edges!')
      RETURN
    END IF

    EdgeCoeff = ListGetConstReal( BC,'Projector Edge Multiplier',Found )
    IF( .NOT. Found ) EdgeCoeff = ListGetConstReal( CurrentModel % Simulation,&
        'Projector Edge Multiplier',Found )
    IF( .NOT. Found ) EdgeCoeff = 1.0_dp

    NodeCoeff = ListGetConstReal( BC,'Projector Node Multiplier',Found )
    IF( .NOT. Found ) NodeCoeff = ListGetConstReal( CurrentModel % Simulation,&
        'Projector Node Multiplier',Found )
    IF( .NOT. Found ) NodeCoeff = 1.0_dp

    Rotational = ListGetLogical( BC,'Rotational Projector',Found ) .OR. &
        ListGetLogical( BC,'Anti Rotational Projector',Found )
    Cylindrical = ListGetLogical( BC,'Cylindrical Projector',Found ) 
    
    MaxDistance = ListGetCReal( BC,'Projector Max Distance', HaveMaxDistance) 
    IF(.NOT. HaveMaxDistance ) THEN
      MaxDistance = ListGetCReal( CurrentModel % Solver % Values,&
          'Projector Max Distance', HaveMaxDistance)       
    END IF

    Naxial = ListGetInteger( BC,'Axial Projector Periods',Found ) 

    Parallel = ( ParEnv % PEs > 1 )
    Mesh => CurrentModel % Mesh
    BMesh1 % Parent => NULL()
    BMesh2 % Parent => NULL()

    ! Create a projector in style P=I-Q, or rather just P=Q. 
    SelfProject = .TRUE.
    
    ! Range is needed to define tolerances, and to map the angle in case 
    ! the master mesh is treated as a repeating structure. 
    XMaxAll = MAXVAL(BMesh2 % Nodes % x)
    XMinAll = MINVAL(BMesh2 % Nodes % x)
    XRange = XMaxAll - XMinAll

    YMaxAll = MAXVAL(BMesh2 % Nodes % y)
    YMinAll = MINVAL(BMesh2 % Nodes % y)
    YRange = YMaxAll - YMinAll

    ! Fix here the relative tolerance used to define the search tolerance
    RelTolY = 1.0d-4
    ! In the case of infinite target we can have tighter criteria
    IF( FullCircle .OR. Repeating ) THEN
      RelTolX = 1.0d-6
    ELSE
      RelTolX = RelTolY
    END IF
    YTol = RelTolY * YRange
    XTol = RelTolX * XRange

    ! Determine the coefficient that turns possible angles into units of
    ! ach-lenth. If this is not rotational then there are no angles. 
    IF( Rotational .OR. Cylindrical ) THEN
      ArcCoeff = (2*PI*Radius)/360.0_dp
    ELSE
      ArcCoeff = 1.0_dp
    END IF

    ! We have a weak projector if it is requested 
    WeakProjector = ListGetLogical( BC, 'Galerkin Projector', Found )    

    StrongProjector = ListGetLogical( BC,'Level Projector Strong',Found )
    IF( StrongProjector .AND. WeakProjector ) THEN
      CALL Fatal(Caller,'Projector cannot be weak (Galerkin) and strong at the same time!')
    END IF
    
    MeshDim = Mesh % MeshDim
    IF( MeshDim == 3 ) THEN
      Element => BMesh1 % Elements(1)
      IF( Element % TYPE % DIMENSION == 1 ) THEN
        CALL Warn(Caller,'Enforcing 1D integration for 1D boundary elements in 3D mesh!')
        MeshDim = 2
      END IF
    END IF
    
    ! Generic integrator does not make any assumptions on the way the mesh 
    ! is constructured. Otherwise constant strides in y-direction is assumed. 
    ! For weak strategy always use the generic integrator. 
    GenericIntegrator = ListGetLogical( BC,'Level Projector Generic',Found ) 
    IF(.NOT. Found ) GenericIntegrator = WeakProjector

    ! Maximum skew in degrees before treating edges as skewed
    SkewTol = 0.1_dp

    ! Check whether generic integrator should be enforced
    IF( DoEdges .AND. .NOT. GenericIntegrator ) THEN
      IF( Naxial > 0 ) THEN
        GenericIntegrator = .TRUE.
        CALL Info(Caller,'Generic integrator enforced for axial projector',Level=6)
      END IF
      
      ! It is assumed that that the target mesh is always un-skewed 
      ! Make a test here to be able to skip it later. No test is needed
      ! if the generic integrator is enforced. 
      IF(.NOT. GenericIntegrator ) THEN
        MaxSkew1 = CheckMeshSkew( BMesh1, NotAllQuads )
        IF( NotAllQuads ) THEN
          CALL Info(Caller,'This mesh has also triangles',Level=8)
        END IF
        WRITE( Message,'(A,ES12.3)') 'Maximum skew in this mesh: ',MaxSkew1
        CALL Info(Caller,Message,Level=8)
        
        MaxSkew2 = CheckMeshSkew( BMesh2, NotAllQuads2 )
        IF( NotAllQuads2 ) THEN
          CALL Info(Caller,'Target mesh has also triangles',Level=8)
        END IF
        WRITE( Message,'(A,ES12.3)') 'Maximum skew in target mesh: ',MaxSkew2
        CALL Info(Caller,Message,Level=8)
        
        IF( NotAllQuads .OR. NotAllQuads2 .OR. MaxSkew2 > SkewTol ) THEN
          IF( MaxSkew2 > MaxSkew1 .AND. MaxSkew1 < SkewTol ) THEN
            CALL Warn(Caller,'You could try switching the master and target BC!')
          END IF
          CALL Warn(Caller,'Target mesh has too much skew, using generic integrator when needed!')
          GenericIntegrator = .TRUE. 
        END IF
      END IF
      
      IF( GenericIntegrator ) THEN
        CALL Info(Caller,'Edge projection for the BC requires weak projector!',Level=7)
        CALL Fatal(Caller,'We cannot use fully strong projector as wished in this geometry!')
      END IF
    END IF
      
    ! The projectors for nodes and edges can be created either in a strong way 
    ! or weak way in the special case that the nodes are located in extruded layers. 
    ! The strong way results to a sparse projector. For constant 
    ! levels it can be quite optimal, except for the edges with a skew. 
    ! If strong projector is used for all edges then "StrideProjector" should 
    ! be recovered.
               
    IF( DoNodes ) THEN
      StrongNodes = ListGetLogical( BC,'Level Projector Nodes Strong',Found ) 

      StrongConformingNodes = ListGetLogical( BC,'Level Projector Conforming Nodes Strong', Found ) 

      IF(.NOT. Found) StrongNodes = ListGetLogical( BC,'Level Projector Strong',Found ) 
      IF(.NOT. Found) StrongNodes = .NOT. GenericIntegrator
    END IF

    IF( DoEdges ) THEN
      StrongEdges = ListGetLogical( BC,'Level Projector Strong',Found )
      IF(.NOT. Found ) StrongEdges = ListGetLogical( BC,'Level Projector Plane Edges Strong', Found ) 
      IF(.NOT. Found ) StrongEdges = .NOT. GenericIntegrator
      
      StrongLevelEdges = ListGetLogical( BC,'Level Projector Plane Edges Strong', Found ) 
      IF( .NOT. Found ) StrongLevelEdges = StrongEdges
      IF( StrongLevelEdges .AND. GenericIntegrator ) THEN
        CALL Info(Caller,'Using strong level edges with partially weak projector',Level=7)
      END IF

      StrongConformingEdges = ListGetLogical( BC,'Level Projector Conforming Edges Strong', Found ) 
      
      StrongExtrudedEdges = ListGetLogical( BC,'Level Projector Extruded Edges Strong', Found ) 
      IF( .NOT. Found ) StrongExtrudedEdges = StrongEdges
      IF( StrongExtrudedEdges .AND. GenericIntegrator ) THEN
        CALL Info(Caller,'Using strong extruded edges with partially weak projector',Level=7)
      END IF
      
      ! There is no strong strategy for skewed edges currently
      StrongSkewEdges = .FALSE.
    END IF


    ! If the number of periods is enforced use that instead since
    ! the Xrange periodicity might not be correct if the mesh has skew.
    IF( Rotational ) THEN
      IF( FullCircle ) THEN
        Xrange = 360.0_dp
      ELSE 
        i = ListGetInteger( BC,'Rotational Projector Periods',Found,minv=1 ) 
        IF( GenericIntegrator .AND. .NOT. Found ) THEN
          CALL Fatal(Caller,&
              'Generic integrator requires > Rotational Projector Periods <')
        END IF
        Xrange = 360.0_dp / i
      END IF
    END IF

    ! This is the tolerance used to define constant direction in radians
    ! For consistency it should not be sloppier than the SkewTol
    ! but it could be equally sloppy as below.
    RadTol = PI * SkewTol / 180.0_dp

    ! Given the inverse permutation compute the initial number of
    ! nodes in both cases. 
    NoNodes1 = BMesh1 % NumberOfNodes
    NoNodes2 = BMesh2 % NumberOfNodes

    InvPerm1 => BMesh1 % InvPerm
    InvPerm2 => BMesh2 % InvPerm

    ! Create a list matrix that allows for unspecified entries in the matrix 
    ! structure to be introduced.
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_GALERKIN

    CreateDual = ListGetLogical( BC,'Create Dual Projector',Found ) 
    IF( CreateDual ) THEN
      DualProjector => AllocateMatrix()
      DualProjector % FORMAT = MATRIX_LIST
      DualProjector % ProjectorType = PROJECTOR_TYPE_GALERKIN
      Projector % EMatrix => DualProjector
    END IF

    ! Check whether biorthogonal basis for projectors requested:
    ! ----------------------------------------------------------
    BiOrthogonalBasis = ListGetLogical( BC, 'Use Biorthogonal Basis', Found)

    ! If we want to eliminate the constraints we have to have a biortgonal basis
    IF(.NOT. Found ) THEN
      BiOrthogonalBasis = ListGetLogical( CurrentModel % Solver % Values, &
          'Eliminate Linear Constraints',Found )
      IF( BiOrthogonalBasis ) THEN
        CALL Info(Caller,&
            'Enforcing > Use Biorthogonal Basis < to True to enable elimination',Level=8)
        CALL ListAddLogical( BC, 'Use Biorthogonal Basis',.TRUE. )
      END IF
    END IF

    IF (BiOrthogonalBasis) THEN
      IF( DoEdges ) THEN
        CALL Warn(Caller,'Biorthogonal basis cannot be combined with edge elements!')
      END IF

      DualSlave  = ListGetLogical(BC, 'Biorthogonal Dual Slave', Found)
      IF(.NOT.Found) DualSlave  = .TRUE.

      DualMaster = ListGetLogical(BC, 'Biorthogonal Dual Master', Found)
      IF(.NOT.Found) DualMaster = .TRUE.

      DualLCoeff = ListGetLogical(BC, 'Biorthogonal Dual Lagrange Coefficients', Found)
      IF(.NOT.Found) DualLCoeff = .FALSE.

      IF(DualLCoeff) THEN
        DualSlave  = .FALSE.
        DualMaster = .FALSE.
        CALL ListAddLogical( CurrentModel % Solver % Values, 'Use Transpose Values',.FALSE.)
      ELSE
        CALL ListAddLogical( CurrentModel % Solver % Values, 'Use Transpose Values',.TRUE.)
      END IF

      Projector % Child => AllocateMatrix()
      Projector % Child % Format = MATRIX_LIST
      CALL Info(Caller,'Using biorthogonal basis, as requested',Level=8)      
    END IF


    PiolaVersion = ListGetLogical( CurrentModel % Solver % Values, &
        'Use Piola Transform', Found)
    SecondOrder = ListGetLogical( CurrentModel % Solver % Values, &
        'Quadratic Approximation', Found)

    ! At the 1st stage determine the maximum size of the projector
    ! If the strong projector is used then the numbering is done as we go
    ! this way we can eliminate unneeded rows. 
    ! For the weak projector there is no need to eliminate rows. 
    IF( DoNodes ) THEN      
      ALLOCATE( NodePerm( Mesh % NumberOfNodes ) )
      NodePerm = 0

      ! in parallel only consider nodes that truly are part of this partition
      DO i=1,BMesh1 % NumberOfBulkElements
        Element => BMesh1 % Elements(i)        
        IF( Parallel ) THEN
          IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE          
        END IF        
        NodePerm( InvPerm1( Element % NodeIndexes ) ) = 1
      END DO

      n = SUM( NodePerm )
      CALL Info(Caller,'Initial number of slave nodes '//TRIM(I2S(n))//&
          ' out of '//TRIM(I2S(BMesh1 % NumberOfNodes ) ), Level = 10 )

      ! Eliminate the redundant nodes by default. 
      ! These are noded that depend on themselves.
      EliminateUnneeded = ListGetLogical( BC,&
          'Level Projector Eliminate Redundant Nodes',Found ) 
      IF(.NOT. Found ) EliminateUnneeded = .TRUE.

      IF( EliminateUnneeded ) THEN
        m = 0
        n = SUM( NodePerm )
        CALL Info(Caller,&
            'Number of potential nodes in projector: '//TRIM(I2S(n)),Level=10)        
        ! Now eliminate the nodes which also occur in the other mesh
        ! These must be redundant edges
        DO i=1, SIZE(InvPerm2)
          j = InvPerm2(i) 
          IF( NodePerm(j) /= 0 ) THEN
            NodePerm(j) = 0
            !PRINT *,'Removing node:',j,Mesh % Nodes % x(j), Mesh % Nodes % y(j)
            m = m + 1
          END IF
        END DO
        IF( m > 0 ) THEN
          CALL Info(Caller,&
              'Eliminating redundant nodes from projector: '//TRIM(I2S(m)),Level=10)
        END IF
      END IF
      
      IF( CreateDual ) THEN
        ALLOCATE( DualNodePerm( Mesh % NumberOfNodes ) )
        DualNodePerm = 0

        DO i=1,BMesh2 % NumberOfBulkElements
          Element => BMesh2 % Elements(i)        
          IF( Parallel ) THEN
            IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE          
          END IF
          DualNodePerm( InvPerm2( Element % NodeIndexes ) ) = 1
        END DO
                
        IF( EliminateUnneeded ) THEN
          m = 0
          n = SUM( DualNodePerm )
          CALL Info(Caller,&
              'Number of potential nodes in dual projector: '//TRIM(I2S(n)),Level=10)        
          ! Now eliminate the nodes which also occur in the other mesh
          ! These must be redundant edges
          DO i=1, SIZE(InvPerm1)
            j = InvPerm1(i) 
            IF( DualNodePerm(j) /= 0 ) THEN
              DualNodePerm(j) = 0
              PRINT *,'Removing dual node:',j,Mesh % Nodes % x(j), Mesh % Nodes % y(j)
              m = m + 1
            END IF
          END DO
          IF( m > 0 ) THEN
            CALL Info(Caller,&
                'Eliminating redundant dual nodes from projector: '//TRIM(I2S(m)),Level=10)
          END IF
        END IF
      END IF
      
      IF( ListCheckPresent( BC,'Level Projector Condition') ) THEN
        ALLOCATE( Cond( Mesh % MaxElementNodes ) )
        Cond = 1.0_dp
        m = 0
        DO i=1, BMesh1 % NumberOfBulkElements          
          Element => Mesh % Elements( BMesh1 % Elements(i) % ElementIndex )
          CurrentModel % CurrentElement => Element
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          Cond(1:n) = ListGetReal( BC,'Level Projector Condition', n, NodeIndexes )
          DO j=1,n
            k = NodeIndexes(j)
            IF( NodePerm(k) /= 0 ) THEN
              IF( Cond(j) < 0.0 ) THEN
                m = m + 1
                NodePerm(k) = 0 
              END IF
            END IF
          END DO
        END DO
        CALL Info(Caller,'Eliminated nodes with negative condition: '//&
            TRIM(I2S(m)),Level=10)        
        DEALLOCATE( Cond ) 
      END IF
      
      m = 0
      DO i=1,Mesh % NumberOfNodes
        IF( NodePerm(i) > 0 ) THEN
          m = m + 1
          NodePerm(i) = m
        END IF
      END DO
      
      CALL Info(Caller,&
          'Number of active nodes in projector: '//TRIM(I2S(m)),Level=8)
      EdgeRow0 = m
      
      IF( CreateDual ) THEN
        m = 0
        DO i=1,Mesh % NumberOfNodes
          IF( DualNodePerm(i) > 0 ) THEN
            m = m + 1
            DualNodePerm(i) = m
          END IF
        END DO
        ALLOCATE( DualProjector % InvPerm(m) )
        DualProjector % InvPerm = 0

        IF( DoEdges ) THEN
          CALL Fatal(Caller,'Dual projector cannot handle edges!')
        END IF
      END IF
    ELSE
      EdgeRow0 = 0
    END IF
    ProjectorRows = EdgeRow0

    IF( DoEdges ) THEN
      ALLOCATE( EdgePerm( Mesh % NumberOfEdges ) )
      EdgePerm = 0

      ! Mark the edges for which the projector must be created for
      DO i=1, BMesh1 % NumberOfBulkElements

        ! in parallel only consider face elements that truly are part of this partition
        IF( Parallel ) THEN
          IF( BMesh1 % Elements(i) % PartIndex /= ParEnv % MyPe ) CYCLE          
        END IF

        DO j=1, BMesh1 % Elements(i) % TYPE % NumberOfEdges
          EdgePerm( BMesh1 % Elements(i) % EdgeIndexes(j) ) = 1
        END DO
      END DO

      EliminateUnneeded = ListGetLogical( BC,&
          'Level Projector Eliminate Redundant Edges',Found )
      IF(.NOT. Found ) EliminateUnneeded = .TRUE.

      IF( EliminateUnneeded ) THEN
        n = SUM( EdgePerm )
        CALL Info(Caller,&
            'Number of potential edges in projector: '//TRIM(I2S(n)),Level=10)        
        ! Now eliminate the edges which also occur in the other mesh
        ! These must be redundant edges
        DO i=1, BMesh2 % NumberOfBulkElements
          DO j=1, BMesh2 % Elements(i) % TYPE % NumberOfEdges
            EdgePerm( BMesh2 % Elements(i) % EdgeIndexes(j) ) = 0
          END DO
        END DO

        IF( DoNodes ) THEN
          IF( ListGetLogical( BC,'Level Projector Eliminate Edges Greedy',Found ) ) THEN
            DO i=1, BMesh1 % NumberOfBulkElements
              DO j=1, BMesh1 % Elements(i) % TYPE % NumberOfEdges
                k = BMesh1 % Elements(i) % EdgeIndexes(j) 
                IF( ANY( NodePerm( Mesh % Edges(k) %  NodeIndexes ) == 0 ) ) THEN
                  EdgePerm( k ) = 0
                END IF
              END DO
            END DO
          END IF
        END IF
      END IF

      m = 0
      DO i=1,Mesh % NumberOfEdges
        IF( EdgePerm(i) > 0 ) THEN
          m = m + 1
          EdgePerm(i) = m
        END IF
      END DO

      IF( EliminateUnneeded ) THEN
        CALL Info(Caller,&
            'Eliminating redundant edges from projector: '//TRIM(I2S(n-m)),Level=10)
      END IF
      CALL Info(Caller,&
          'Number of active edges in projector: '//TRIM(I2S(m)),Level=8)
      IF (SecondOrder) THEN
        FaceRow0 = EdgeRow0 + 2*m
      ELSE
        FaceRow0 = EdgeRow0 + m
      END IF
      ProjectorRows = FaceRow0
      
      IF( PiolaVersion ) THEN
        ! Note: this might not work in parallel with halo since some of the face elements
        ! do not then belong to the slave boundary. 
        m = 0
        DO i=1,BMesh1 % NumberOfBulkElements
          m = m + BMesh1 % Elements(i) % BDOFs
        END DO
        CALL Info(Caller,&
            'Number of active faces in projector: '//TRIM(I2S(BMesh1 % NumberOfBulkElements)),Level=8)
        CALL Info(Caller,&
            'Number of active face DOFs in projector: '//TRIM(I2S(m)),Level=8)
        ProjectorRows = FaceRow0 + m
      END IF
    END IF

    CALL Info(Caller,&
        'Max number of rows in projector: '//TRIM(I2S(ProjectorRows)),Level=10)
    ALLOCATE( Projector % InvPerm(ProjectorRows) )
    Projector % InvPerm = 0

    ! If after strong projectors there are still something undone they must 
    ! be dealt with the weak projectors. 
    SomethingUndone = .FALSE.

    ! If requested, create strong mapping for node dofs
    !------------------------------------------------------------------   
    IF( DoNodes ) THEN
      IF( StrongConformingNodes ) THEN
        CALL AddNodeProjectorStrongConforming()
      ELSE IF( StrongNodes ) THEN
        IF( GenericIntegrator ) THEN 
          CALL AddNodalProjectorStrongGeneric()
        ELSE
          CALL AddNodalProjectorStrongStrides()
        END IF
      ELSE
        ! If strong projector is applied they can deal with all nodal dofs
        SomethingUndone = .TRUE.
      END IF
    END IF

    ! If requested, create strong mapping for edge dofs
    !-------------------------------------------------------------
    EdgeBasis = .FALSE.
    IF( DoEdges ) THEN
      EdgeCol0 = Mesh % NumberOfNodes
      IF (SecondOrder) THEN
        FaceCol0 = Mesh % NumberOfNodes + 2 * Mesh % NumberOfEdges
      ELSE
        FaceCol0 = Mesh % NumberOfNodes + Mesh % NumberOfEdges
      END IF

      IF( StrongLevelEdges .OR. StrongExtrudedEdges .OR. StrongConformingEdges ) THEN
        IF( StrongConformingEdges ) THEN
          CALL AddEdgeProjectorStrongConforming()
        ELSE
          CALL AddEdgeProjectorStrongStrides()
        END IF
        ! Compute the unset edge dofs. 
        ! Some of the dofs may have been set by the strong projector. 
        m = COUNT( EdgePerm > 0 )
        IF( m > 0 ) THEN
          CALL Info(Caller,&
              'Number of weak edges in projector: '//TRIM(I2S(m)),Level=10)      
        END IF
        IF( m > 0 .OR. PiolaVersion) THEN
          SomethingUndone = .TRUE.
          EdgeBasis = .TRUE.
        END IF
      ELSE
        SomethingUndone = .TRUE.
        EdgeBasis = .TRUE.
      END IF      
    END IF

    ! And the the rest
    !-------------------------------------------------------------
    IF( SomethingUndone ) THEN      
      IF( MeshDim == 2 ) THEN
        CALL Info(Caller,'Initial mesh is 2D, using 1D projectors!',Level=10) 
        CALL AddProjectorWeak1D()
      ELSE IF( GenericIntegrator ) THEN
        CALL AddProjectorWeakGeneric()
      ELSE
        CALL AddProjectorWeakStrides()
      END IF
    END IF

    ! Now change the matrix format to CRS from list matrix
    !--------------------------------------------------------------
    CALL List_toCRSMatrix(Projector)
    CALL CRS_SortMatrix(Projector,.TRUE.)
    CALL Info(Caller,'Number of rows in projector: '&
        //TRIM(I2S(Projector % NumberOfRows)),Level=12)
    CALL Info(Caller,'Number of entries in projector: '&
        //TRIM(I2S(SIZE(Projector % Values))),Level=12)
  

    IF(ASSOCIATED(Projector % Child)) THEN
      CALL List_toCRSMatrix(Projector % Child)
      CALL CRS_SortMatrix(Projector % Child,.TRUE.)
    END IF

    IF( CreateDual ) THEN
      CALL List_toCRSMatrix(DualProjector)
      CALL CRS_SortMatrix(DualProjector,.TRUE.)
    END IF
    
    IF( DoNodes ) DEALLOCATE( NodePerm )
    IF( CreateDual .AND. DoNodes ) DEALLOCATE( DualNodePerm )
    IF( DoEdges ) DEALLOCATE( EdgePerm )

    m = COUNT( Projector % InvPerm  == 0 ) 
    IF( m > 0 ) THEN
      CALL Warn(Caller,'Projector % InvPerm not set in for dofs: '//TRIM(I2S(m)))
    END IF

    CALL Info(Caller,'Projector created',Level=10)
    
  CONTAINS

    ! Currently the target mesh is assumed to be include only cartesian elements
    ! Check the angle in the elements. When we know the target mesh is cartesian
    ! we can reduce the error control in the other parts of the code. 
    !----------------------------------------------------------------------------
    FUNCTION CheckMeshSkew(BMesh, NotAllQuads) RESULT( MaxSkew )

      TYPE(Mesh_t),POINTER :: BMesh
      REAL(KIND=dp) :: MaxSkew
      LOGICAL :: NotAllQuads

      INTEGER :: i,j,n,indM,k,knext,kprev
      TYPE(Element_t), POINTER :: ElementM
      TYPE(Nodes_t) :: NodesM
      REAL(KIND=dp) :: e1(2),e2(2),DotProdM, PhiM
      INTEGER, POINTER :: IndexesM(:)

      CALL Info(Caller,'Checking mesh skew')

      n = 4
      ALLOCATE( NodesM % x(n), NodesM % y(n) )
      MaxSkew = 0.0_dp
      NotAllQuads = .FALSE.
      
      j = 0
      DO indM=1,BMesh % NumberOfBulkElements
        
        ElementM => BMesh % Elements(indM)        
        n = ElementM % TYPE % ElementCode / 100
        IF( n /= 4 ) THEN
          NotAllQuads = .TRUE.
        END IF
        IndexesM => ElementM % NodeIndexes
        NodesM % y(1:n) = BMesh % Nodes % y(IndexesM(1:n))
        NodesM % x(1:n) = BMesh % Nodes % x(IndexesM(1:n))
        
        ! Transfer into real length units instead of angles
        ! This gives right balance between x and y -directions. 
        NodesM % x(1:n) = ArcCoeff * NodesM % x(1:n)
        
        ! Make unit vectors of the edge
        DO k = 1, n
          knext = MODULO(k,n)+1
          kprev = MODULO(n+k-2,n)+1
          
          e1(1) = NodesM % x(knext) - NodesM % x(k) 
          e1(2) = NodesM % y(knext) - NodesM % y(k) 
          
          e2(1) = NodesM % x(kprev) - NodesM % x(k) 
          e2(2) = NodesM % y(kprev) - NodesM % y(k) 
          
          e1 = e1 / SQRT( SUM( e1**2) )
          e2 = e2 / SQRT( SUM( e2**2) )
          
          ! dot product of the unit vectors
          DotProdM = SUM( e1 * e2 )
          
          ! Cosine angle in degrees        
          PhiM = ACOS( DotProdM ) 
          MaxSkew = MAX( MaxSkew, ABS ( ABS( PhiM ) - PI/2 ) )
        END DO
      END DO

      ! Move to degrees and give the tolerance in them
      MaxSkew = MaxSkew * 180.0_dp / PI
        
100   DEALLOCATE( NodesM % x, NodesM % y )

    END FUNCTION CheckMeshSkew
      

    !-------------------------------------------------------------------------------------
    ! Create projector for nodes on the strides directly from a linear 
    ! combination of two nodes. This approach minimizes the size of the projector
    ! and also minimizes the need for parallel communication.
    !-------------------------------------------------------------------------------------
    SUBROUTINE AddNodalProjectorStrongStrides()

      TYPE(Element_t), POINTER :: ElementM
      INTEGER, POINTER :: IndexesM(:)
      INTEGER :: ncoeff, coeffi(2),sgn0, ind, indm, j1, j2, j3, Nundefined
      REAL(KIND=dp) :: x1, y1, x2, y2, xmin, xmax, xminm, xmaxm, Dist, MinDist
      REAL(KIND=dp) :: coeff(2), val, xm1, xm2, xm3
      INTEGER, POINTER :: EdgeMap(:,:)
      TYPE(Nodes_t) :: NodesM
      LOGICAL :: LeftCircle

      CALL Info(Caller,'Creating strong stride projector for nodal dofs',Level=10)

      n = Mesh % MaxElementNodes
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      NodesM % z = 0.0_dp

      ! By construction there is always two components in the projector for the nodes. 
      ncoeff = 2
      coeffi = 0
      sgn0 = 1
      Nundefined = 0

      ! This flag tells if we're working with a full circle and the problematic part of 
      ! the circle with the discontinuity in the angle. 
      LeftCircle = .FALSE.

      DO ind=1,BMesh1 % NumberOfNodes

        nrow = NodePerm( InvPerm1( ind ) )
        IF( nrow == 0 ) CYCLE
        NodePerm( InvPerm1( ind ) ) = 0
        Projector % InvPerm(nrow) = InvPerm1(ind)

        Found = .FALSE.
        x1 = BMesh1 % Nodes % x(ind)
        y1 = BMesh1 % Nodes % y(ind)
        sgn0 = 1
        coeff = 0.0_dp
        MinDist = HUGE( MinDist )

        IF( Repeating ) THEN
          Nrange = FLOOR( (x1-XMinAll) / XRange )
          x1 = x1 - Nrange * XRange
          
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
          END IF
        ELSE IF( FullCircle ) THEN
          LeftCircle = ABS( x1 ) > 90.0_dp
          IF( LeftCircle ) THEN
            IF( x1 < 0.0 ) x1 = x1 + 360.0_dp
          END IF
        END IF

        ! If the projector is of style Px+Qx=0 then
        ! and the negative sign, otherwise let the initial sign be.
        IF( SelfProject ) sgn0 = -sgn0
        
        ! Currently a cheap n^2 loop but it could be improved
        ! Looping over master elements. Look for constant-y strides only. 
        !--------------------------------------------------------------------
        DO indM = 1, BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)
          n = ElementM % TYPE % NumberOfNodes        
          IndexesM => ElementM % NodeIndexes
          
          ! Quick tests to save time
          ! Element must have nodes at the right level
          NodesM % y(1:n) = BMesh2 % Nodes % y(IndexesM(1:n))           
          IF( ALL( ABS( NodesM % y(1:n) - y1 ) > YTol ) ) CYCLE

          ! The x nodes should be in the interval
          NodesM % x(1:n) = BMesh2 % Nodes % x(IndexesM(1:n))

          ! Transform the master element on-the-fly around the problematic angle
          IF( LeftCircle ) THEN
            ! The master nodes are all on right
            IF( ALL( ABS( NodesM % x(1:n) ) - 90.0_dp < Xtol ) ) CYCLE
            DO j=1,n
              IF( NodesM % x(j) < 0.0 ) NodesM % x(j) = NodesM % x(j) + 360.0_dp
            END DO
          END IF
          
          xmaxm = MAXVAL( NodesM % x(1:n) )
          xminm = MINVAL( NodesM % x(1:n) )

          ! Eliminate this special case since it could otherwise give a faulty hit
          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > 180.0_Dp ) CYCLE
          END IF

          Dist = MAX( x1-xmaxm, xminm-x1 ) 

          ! Mark the minimum distance if this would happen to be a problematic node
          MinDist = MIN( Dist, MinDist )

          IF( Dist > Xtol ) CYCLE

          ! Ok, this may be a proper element, now just find the two nodes
          ! needed for the mapping on the same stride. Basically this means 
          ! finding the correct edge but we don't need to use the data structure for that. 
          ! For 1D edge element this is trivial, note however that only 1st degree projection is used!
          j1 = 0; j2 = 0; j3 = 0
          IF( n <= 3 ) THEN
            j1 = 1 
            j2 = 2
            IF( n == 3 ) j3 = 3
          ELSE
            DO j=1,n
              IF( ABS( NodesM % y(j) - y1 ) > YTol ) CYCLE
              IF( j1 == 0 ) THEN
                j1 = j
              ELSE IF( j2 == 0 ) THEN
                j2 = j
              ELSE
                j3 = j
                ! This means that for higher order edges only three nodes are used
                EXIT
              END IF
            END DO
            IF( j2 == 0 ) CALL Warn(Caller,'Could not locate an edge consistently!')
          END IF

          ! The node to map must be in interval, x1 \in [xm1,xm2]
          IF( NodesM % x(j1) > NodesM % x(j2) ) THEN
             j = j2; j2 = j1; j1 = j
          END IF
          xm1 = NodesM % x(j1)
          xm2 = NodesM % x(j2)          

          ! We are at interval [xm1,xm2] now choose either [xm1,xm3] or [xm3,xm2]
          IF( j3 > 0 ) THEN
             xm3 = NodesM % x(j3)          
             IF( x1 > xm3 ) THEN
                j1 = j3; xm1 = xm3
             ELSE 
                j2 = j3; xm2 = xm3
             END IF
          END IF
          
          ! Ok, the last check, this might fail if the element had skew even though the 
          ! quick test is successful! Then the left and right edge may have different range.
          Dist = MAX( x1-xm2, xm1-x1 )
          IF( Dist > Xtol ) CYCLE

          ! When we have the correct edge, the mapping is trivial.
          ! The sum of weights of the projectors is set to one. 
          IF( ABS(xm1-xm2) < TINY(xm1) ) THEN
            CALL Warn(Caller,'Degenerated edge?')
            PRINT *,'ind',ind,x1,y1,xm1,xm2,j1,j2,j3
            PRINT *,'x:',NodesM % x(1:n)
            PRINT *,'y:',NodesM % y(1:n)
            coeff(1) = 0.5_dp
          ELSE
            coeff(1) = (xm2-x1)/(xm2-xm1) 
          END IF
          coeff(2) = 1.0_dp - coeff(1)

          coeffi(1) = IndexesM(j1)
          coeffi(2) = IndexesM(j2)

          Found = .TRUE.
          
          ! If we really exactly between [xm1,xm2] then we may finish the search for good
          IF( Dist < EPSILON( Dist ) ) EXIT
        END DO

        IF(.NOT. Found ) THEN
          Nundefined = Nundefined + 1
          WRITE( Message,'(A,2I8,3ES12.3)') 'Problematic node: ',&
              ind,ParEnv % MyPe,x1,y1,MinDist
          CALL Warn(Caller,Message)
          CYCLE
        END IF

        IF( SelfProject ) THEN
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              InvPerm1(ind), NodeCoeff ) 
        END IF

        ! The scaling of the projector entries is used, for example, 
        ! to allow antiperiodic projectors. 
        Coeff(1:ncoeff) = sgn0 * Coeff(1:ncoeff)

        ! The projection weights
        DO j=1,ncoeff 

          val = Coeff(j)
          ! Skip too small projector entries
          IF( ABS( val ) < 1.0d-12 ) CYCLE

          ! Use the permutation to revert to original dofs
          CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
              InvPerm2(coeffi(j)), NodeScale * NodeCoeff * val ) 
        END DO

      END DO

      IF( Nundefined > 0 ) THEN
        CALL Warn(Caller,&
            'Nodes could not be determined by any edge: '//TRIM(I2S(Nundefined)))          
      END IF

      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )


    END SUBROUTINE AddNodalProjectorStrongStrides
    !---------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------
    ! Adds a nodal projector assuming generic 2D mesh. 
    ! Otherwise should give same results as the one before. 
    !---------------------------------------------------------------------------------
    SUBROUTINE AddNodalProjectorStrongGeneric()

      TYPE(Element_t), POINTER :: ElementM
      INTEGER, POINTER :: IndexesM(:), coeffi(:)
      REAL(KIND=dp), POINTER :: Basis(:),coeff(:)
      INTEGER :: n, nM, ncoeff, sgn0, ind, indm, j1, j2, j3, Nundefined
      REAL(KIND=dp) :: x1, y1, z1, xmin, xmax, xminm, xmaxm, ymaxm, yminm, &
          Dist, MaxMinBasis, detJ, ArcTol, ArcRange
      REAL(KIND=dp) :: val, u, v, w
      TYPE(Nodes_t) :: NodesM
      LOGICAL :: LeftCircle, Found, Stat

      CALL Info(Caller,'Creating strong generic projector for nodal dofs',Level=10)

      n = Mesh % MaxElementNodes
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n), Basis(n), coeff(n), coeffi(n) )
      NodesM % z = 0.0_dp

      ncoeff = 0
      coeffi = 0
      sgn0 = 1
      Nundefined = 0
      z1 = 0.0_dp

      ArcTol = ArcCoeff * Xtol
      ArcRange = ArcCoeff * Xrange 

      ! This flag tells if we're working with a full circle and the problematic part of 
      ! the circle with the discontinuity in the angle. 
      LeftCircle = .FALSE.

      DO ind=1,BMesh1 % NumberOfNodes

        nrow = NodePerm( InvPerm1( ind ) )
        IF( nrow == 0 ) CYCLE
        NodePerm( InvPerm1( ind ) ) = 0
        Projector % InvPerm(nrow) = InvPerm1(ind)

        Found = .FALSE.
        x1 = ArcCoeff * BMesh1 % Nodes % x(ind)
        y1 = BMesh1 % Nodes % y(ind)
        IF( HaveMaxDistance ) THEN
          z1 = BMesh1 % Nodes % z(ind)
        END IF

        sgn0 = 1
        coeff = 0.0_dp
        MaxMinBasis = -HUGE(MaxMinBasis)

        IF( FullCircle ) THEN
          LeftCircle = ABS( x1 ) > ArcCoeff * 90.0_dp
          IF( LeftCircle ) THEN
            IF( x1 < 0.0 ) x1 = x1 + ArcCoeff * 360.0_dp
          END IF
        END IF

        ! If the projector is of style Px+Qx=0 then
        ! and the negative sign, otherwise let the initial sign be.
        IF( SelfProject ) sgn0 = -sgn0
        
        ! Currently a cheap n^2 loop but it could be improved
        ! Looping over master elements. Look for constant-y strides only. 
        !--------------------------------------------------------------------
        DO indM = 1, BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)
          nM = ElementM % TYPE % NumberOfNodes        
          IndexesM => ElementM % NodeIndexes

          IF( HaveMaxDistance ) THEN
            IF( MINVAL( ABS( BMesh2 % Nodes % z(IndexesM(1:nM)) - z1 ) ) > MaxDistance ) CYCLE          
          END IF
          
          ! Quick tests to save time
          NodesM % y(1:nM) = BMesh2 % Nodes % y(IndexesM(1:nM))           
          ymaxm = MAXVAL( NodesM % y(1:nM) )
          yminm = MINVAL( NodesM % y(1:nM) )

          Dist = MAX( y1-ymaxm, yminm-y1 ) 
          IF( Dist > Ytol ) CYCLE

          ! The x nodes should be in the interval
          NodesM % x(1:nM) = BMesh2 % Nodes % x(IndexesM(1:nM))

          ! Transform the master element on-the-fly around the problematic angle
          ! Full 2D circle is never repeating
          IF( LeftCircle ) THEN
            ! The master nodes are all on right
            IF( ALL( ABS( NodesM % x(1:nM) ) - ArcCoeff * 90.0_dp < ArcTol ) ) CYCLE
            DO j=1,nM
              IF( NodesM % x(j) < 0.0 ) NodesM % x(j) = NodesM % x(j) + ArcCoeff * 360.0_dp
            END DO
          END IF
          
          xmaxm = MAXVAL( NodesM % x(1:nM) )
          xminm = MINVAL( NodesM % x(1:nM) )

          ! Eliminate this special case since it could otherwise give a faulty hit
          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > ArcCoeff * 180.0_dp ) CYCLE
          END IF

          IF( Repeating ) THEN
            Nrange = FLOOR( (xmaxm-x1) / XRange )
            IF( Nrange /= 0 ) THEN
              xminm = xminm - Nrange * ArcRange
              xmaxm = xmaxm - Nrange * ArcRange
              NodesM % x(1:nM) = NodesM % x(1:nM) - NRange * ArcRange 
            END IF

            ! Check whether there could be a intersection in an other interval as well
            IF( xminm + ArcRange < x1 + ArcTol ) THEN
              Nrange2 = 1
            ELSE
              Nrange2 = 0
            END IF
          END IF

100       Dist = MAX( x1-xmaxm, xminm-x1 ) 

          IF( Dist < Xtol ) THEN
            ! Integration point at the slave element
            CALL GlobalToLocal( u, v, w, x1, y1, z1, ElementM, NodesM )              
            stat = ElementInfo( ElementM, NodesM, u, v, w, detJ, Basis )
            
            IF( MINVAL( Basis(1:nM) ) > MaxMinBasis ) THEN
              MaxMinBasis = MINVAL( Basis(1:nM) )
              ncoeff = nM
              coeff(1:nM) = Basis(1:nM)
              coeffi(1:nM) = IndexesM(1:nM)
              Found = ( MaxMinBasis >= -1.0d-12 )
            END IF
         
            IF( Found ) EXIT
          END IF
          
          IF( Repeating ) THEN
            IF( NRange2 /= 0 ) THEN
              xminm = xminm + ArcCoeff * Nrange2 * ArcRange
              xmaxm = xmaxm + ArcCoeff * Nrange2 * ArcRange
              NodesM % x(1:n) = NodesM % x(1:n) + NRange2 * ArcRange 
              NRange = NRange + NRange2
              NRange2 = 0
              GOTO 100
            END IF
          END IF         

        END DO

        IF(.NOT. Found ) THEN
          IF( MaxMinBasis > -1.0d-6 ) THEN
            CALL Info(Caller,Message,Level=8)
            Found = .TRUE.
          ELSE
            Nundefined = Nundefined + 1
            IF( .NOT. HaveMaxDistance ) THEN
              WRITE( Message,'(A,2I8,3ES12.3)') 'Problematic node: ',&
                  ind,ParEnv % MyPe,x1,y1,MaxMinBasis
              CALL Warn(Caller,Message )
            END IF
          END IF
        END IF

        IF( Found ) THEN
          IF( SelfProject ) THEN
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                InvPerm1(ind), NodeCoeff ) 
          END IF
          
          ! The scaling of the projector entries is used, for example, 
          ! to allow antiperiodic projectors. 
          Coeff(1:ncoeff) = sgn0 * Coeff(1:ncoeff)
          
          ! Add the projection weights to the matrix
          DO j=1,ncoeff 
            
            val = Coeff(j)
            ! Skip too small projector entries
            ! These really should sum to one we now the limit quite well
            IF( ABS( val ) < 1.0d-8 ) CYCLE
            
            ! Use the permutation to revert to original dofs
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                InvPerm2(coeffi(j)), NodeScale * NodeCoeff * val ) 
          END DO
        END IF
        
      END DO

      IF( Nundefined > 0 ) THEN
        IF( HaveMaxDistance ) THEN
          CALL Info(Caller,&
              'Nodes could not be found in any element: '//TRIM(I2S(Nundefined)))          
        ELSE
          CALL Warn(Caller,&
              'Nodes could not be found in any element: '//TRIM(I2S(Nundefined)))          
        END IF
      END IF

      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z, Basis, coeffi, coeff )


    END SUBROUTINE AddNodalProjectorStrongGeneric
    !---------------------------------------------------------------------------------


    !---------------------------------------------------------------------------------
    ! Create a projector for edges directly. This minmizes the size of the projector 
    ! but may result to numerically inferior projector compared to the weak projector.
    ! It seems to be ok for unskewed geometries where the simplest edge elements work 
    ! well. For skewed geometries the solution does not easily seem to be compatible
    ! with the strong projector. 
    !---------------------------------------------------------------------------------
    SUBROUTINE AddEdgeProjectorStrongStrides()

      INTEGER :: ind, indm, eind, eindm, k1, k2, km1, km2, sgn0, coeffi(100), &
          ncoeff, dncoeff, ncoeff0, i1, i2, j1, j2, Nundefined, NoSkewed, SkewPart
      TYPE(Element_t), POINTER :: Element, ElementM
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      TYPE(Nodes_t) :: NodesM, Nodes 
      INTEGER, POINTER :: EdgeMap(:,:),EdgeMapM(:,:)
      REAL(KIND=dp) :: xm1, xm2, ym1, ym2, coeff(100), signs(100), wsum, minwsum, maxwsum, val, &
          x1o, y1o, x2o, y2o, cskew, sedge
      REAL(KIND=dp) :: x1, y1, x2, y2, xmin, xmax, xminm, xmaxm, ymin, ymax, yminm, ymaxm, xmean, &
          dx,dy,Xeps
      LOGICAL :: YConst, YConstM, XConst, XConstM, EdgeReady, Repeated, LeftCircle, &
          SkewEdge, AtRangeLimit


      CALL Info(Caller,'Creating strong stride projector for edges assuming strides',Level=10)

      n = Mesh % NumberOfEdges
      IF( n == 0 ) RETURN      

      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      Nodes % z = 0.0_dp
      NodesM % z = 0.0_dp

      minwsum = HUGE( minwsum ) 
      maxwsum = 0.0_dp
      NoSkewed = 0
      Nundefined = 0
      LeftCircle = .FALSE.
      Xeps = EPSILON( Xeps )
      AtRangeLimit = .FALSE.

      DO ind=1,BMesh1 % NumberOfBulkElements
        
        Element => BMesh1 % Elements(ind)        
        EdgeMap => GetEdgeMap( Element % TYPE % ElementCode / 100)

        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))

        dx = MAXVAL( Nodes % x(1:n)) - MINVAL(Nodes % x(1:n))
        dy = MAXVAL( Nodes % y(1:n)) - MINVAL(Nodes % y(1:n))

        ! Go through combinations of edges and find the edges for which the 
        ! indexes are the same. 
        DO i = 1,Element % TYPE % NumberOfEdges
          
          eind = Element % EdgeIndexes(i)
          IF( EdgePerm(eind) == 0 ) CYCLE

          nrow = EdgeRow0 + EdgePerm(eind) 
          
          ! Get the nodes of the edge
          i1 = EdgeMap(i,1) 
          i2 = EdgeMap(i,2)

          k1 = Indexes( i1 )
          k2 = Indexes( i2 )

          ! The coordinates of the edge
          x1 = Nodes % x(i1)
          y1 = Nodes % y(i1)

          x2 = Nodes % x(i2)
          y2 = Nodes % y(i2)

          YConst = ( ABS(y2-y1) < RadTol * dy )
          XConst = ( ABS(x2-x1) < RadTol * dx )

          SkewEdge = .FALSE.
          cskew = 1.0_dp
          
          IF( YConst ) THEN
            IF( .NOT. StrongLevelEdges ) CYCLE         
          ELSE IF( XConst ) THEN
            IF( .NOT. StrongExtrudedEdges ) CYCLE
          ELSE
            !print *,'skewed edge: ',ParEnv % MyPe,x1,x2,y1,y2,dx,dy
            !print *,'tol:',ABS(y2-y1)/dy,ABS(x2-x1)/dx,RadTol

            NoSkewed = NoSkewed + 1
            SkewEdge = .TRUE.
            IF(.NOT. StrongSkewEdges) CYCLE
          END IF
          

          ! Numbering of global indexes is needed to ensure correct direction 
          ! of the edge dofs. Basically the InvPerm could be used also in serial
          ! but the order of numbering is maintained when the reduced mesh is created. 
          IF(Parallel) THEN
            k1 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm1(k1))
            k2 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm1(k2))
          END IF
          ncoeff = 0 

          IF( SkewEdge ) THEN
            SkewPart = 0
            sedge = SQRT(ArcCoeff**2*(x1-x2)**2 + (y1-y2)**2)
            x1o = x1
            y1o = y1
            x2o = x2
            y2o = y2
          END IF

          ! This is mainly a test branch for skewed quadrilaters.
          ! It is based on the composition of a skewed edge into 
          ! four cartesian vectors oriented along x or y -axis. 
          ! Unfortunately the resulting projector does not seem to be 
          ! numerically favourable.           
50        IF( SkewEdge ) THEN
            IF( SkewPart < 2 ) THEN
              XConst = .TRUE.
              YConst = .FALSE.
              IF( SkewPart == 1 ) THEN
                x1 = (3.0_dp*x1o + x2o) / 4.0_dp
              ELSE
                x1 = (x1o + 3.0_dp*x2o) / 4.0_dp
              END IF
              x2 = x1
              y1 = y1o
              y2 = y2o
              cskew = 0.5_dp * ABS(y1-y2) / sedge
            ELSE 
              XConst = .FALSE.
              YConst = .TRUE.
              IF( SkewPart == 2 ) THEN
                x1 = x1o
                x2 = (x1o + x2o) / 2.0_dp
                y1 = y1o
                y2 = y1o
              ELSE
                x1 = (x1o + x2o) / 2.0_dp
                x2 = x2o
                y1 = y2o
                y2 = y2o
              END IF
              cskew = ArcCoeff * ABS(x1-x2) / sedge
            END IF
          END IF 

          ncoeff0 = ncoeff
          dncoeff = 0
          Repeated = .FALSE.

          ! If the edge might be treated in two periodic parts 
          ! then here study whether this is the case (Nrange2 /= 0). 
          IF( Repeating ) THEN
            Nrange = FLOOR( (x1-XMinAll) / XRange )
            x1 = x1 - Nrange * XRange
            x2 = x2 - Nrange * XRange
            
            IF( x2 > XMaxAll ) THEN
              Nrange2 = 1
            ELSE IF( x2 < XMinAll ) THEN
              Nrange2 = -1
            ELSE
              Nrange2 = 0
            END IF
          ELSE IF( FullCircle ) THEN
            ! If we have a full circle then treat the left-hand-side
            ! differently in order to circumvent the discontinuity of the
            ! angle at 180 degrees. 
            LeftCircle = ( ABS(x1) > 90.0_dp .AND. ABS(x2) > 90.0_dp )
            IF( LeftCircle ) THEN
              IF( x1 < 0.0_dp ) x1 = x1 + 360.0_dp
              IF( x2 < 0.0_dp ) x2 = x2 + 360.0_dp
            END IF
          END IF

          EdgeReady = .FALSE.
100       sgn0 = 1
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
          END IF
          
          IF( SelfProject ) sgn0 = -sgn0
          
          xmin = MIN(x1,x2)
          xmax = MAX(x1,x2)
          ymin = MIN(y1,y2)
          ymax = MAX(y1,y2)
          xmean = (x1+x2) / 2.0_dp


          ! If the mesh is not repeating there is a risk that we don't exactly hit the start 
          ! or end of the range. Therefore grow the tolerance close to the ends. 
          IF(.NOT. ( Repeating .OR. FullCircle ) ) THEN
            IF ( xmax < XminAll + Xtol .OR. xmin > XmaxAll - Xtol ) THEN
              Xeps = Xtol 
            ELSE
              Xeps = EPSILON( Xeps ) 
            END IF
          END IF

          
          ! Currently a n^2 loop but it could be improved
          !--------------------------------------------------------------------
          DO indm=1,BMesh2 % NumberOfBulkElements
            
            ElementM => BMesh2 % Elements(indm)        
            n = ElementM % TYPE % NumberOfNodes        
            IndexesM => ElementM % NodeIndexes(1:n)
            
            ! Make first some coarse tests to eliminate most of the candidate elements
            ! The y nodes should always have an exact fit
            NodesM % y(1:n) = BMesh2 % Nodes % y(IndexesM(1:n))           
            IF( MINVAL( ABS( ymin - NodesM % y(1:n) ) ) > YTol ) CYCLE
            IF(.NOT. YConst ) THEN
              IF( MINVAL( ABS( ymax - NodesM % y(1:n) ) ) > YTol ) CYCLE
            END IF
            
            NodesM % x(1:n) = BMesh2 % Nodes % x(IndexesM(1:n))
            
            ! If we have a full circle then treat the left part differently
            IF( LeftCircle ) THEN
              IF( ALL( ABS( NodesM % x(1:n) ) - 90.0_dp < Xtol ) ) CYCLE
              DO j=1,n
                IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = NodesM % x(j) + 360.0_dp
              END DO
            END IF
            
            ! The x nodes should be in the interval
            xminm = MINVAL( NodesM % x(1:n) ) 
            xmaxm = MAXVAL( NodesM % x(1:n) ) 
            
            IF( xminm > xmax + Xeps ) CYCLE
            IF( xmaxm < xmin - Xeps ) CYCLE 
            
            ! Eliminate this special case since it could otherwise give a faulty hit
            IF( FullCircle .AND. .NOT. LeftCircle ) THEN
              IF( xmaxm - xminm > 180.0_dp ) CYCLE
            END IF

            yminm = MINVAL( NodesM % y(1:n) ) 
            ymaxm = MAXVAL( NodesM % y(1:n) ) 
            
            ! Ok, we have found a candicate face that will probably have some hits       
            EdgeMapM => GetEdgeMap( ElementM % TYPE % ElementCode / 100)        
            
            ! Go through combinations of edges and find the edges for which the 
            ! indexes are the same. 
            DO j = 1,ElementM % TYPE % NumberOfEdges

              eindm = ElementM % EdgeIndexes(j)
              
              ! Eliminate the possibilitity that the same edge is accounted for twice
              ! in two different boundary elements. 
              IF( ANY( coeffi(ncoeff0+1:ncoeff) == eindm ) ) CYCLE
              
              j1 = EdgeMap(j,1)
              j2 = EdgeMap(j,2)
              
              km1 = IndexesM( j1 )
              km2 = IndexesM( j2 )
              
              ym1 = NodesM % y(j1)
              ym2 = NodesM % y(j2)
              
              xm1 = NodesM % x(j1)
              xm2 = NodesM % x(j2)
              
              ! The target mesh has already been checked that the elements are rectangular so 
              ! the edges must be have either constant y or x.
              YConstM = ( ABS(ym2-ym1) / (ymaxm-yminm) < ABS(xm2-xm1) / (xmaxm-xminm) )
              XConstM = .NOT. YConstM
              
              ! Either both are lateral edges, or both are vertical
              IF( .NOT. ( ( YConst .AND. YConstM ) .OR. ( XConst .AND. XConstM ) ) ) THEN
                CYCLE
              END IF
              
              ! sign depends on the direction and order of global numbering
              IF(Parallel) THEN
                km1 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm2(km1))
                km2 = CurrentModel % Mesh % ParallelInfo % GlobalDOFs(InvPerm2(km2))
              END IF
              
              IF( YConst ) THEN
                IF( ABS( y1 - ym1 ) > YTol ) CYCLE
                
                ! Check whether the range of master x has a union with the slave x
                xmaxm = MAX( xm1, xm2 ) 
                IF( xmaxm < xmin ) CYCLE
                
                xminm = MIN( xm1, xm2 ) 
                IF( xminm > xmax ) CYCLE

                ! Ok, we have a hit register it 
                ncoeff = ncoeff + 1
                coeffi(ncoeff) = eindm

                ! weight depends on the relative fraction of overlapping
                IF( ABS( xmax-xmin) < TINY( xmax ) ) THEN
                  CALL Warn(Caller,'Degenerated edge 2?')
                  coeff(ncoeff) = cskew * 1.0_dp
                ELSE
                  coeff(ncoeff) = cskew * (MIN(xmaxm,xmax)-MAX(xminm,xmin))/(xmax-xmin)
                END IF

                ! this sets the sign which should be consistent 
                IF( (x1-x2)*(xm1-xm2)*(k1-k2)*(km1-km2) > 0.0_dp ) THEN
                  signs(ncoeff) = sgn0
                ELSE
                  signs(ncoeff) = -sgn0
                END IF

                ! There can be only one lateral edge hit for each element
                EXIT 
              ELSE
                dncoeff = dncoeff + 1
                ncoeff = ncoeff + 1

                IF( (y1-y2)*(ym1-ym2)*(k1-k2)*(km1-km2) > 0.0_dp ) THEN
                  signs(ncoeff) = sgn0 
                ELSE
                  signs(ncoeff) = -sgn0
                END IF

                coeffi(ncoeff) = eindm
                ! note: temporarily save the coordinate to the coefficient!
                coeff(ncoeff) = ( xm1 + xm2 ) / 2.0_dp
              END IF
            END DO

            IF( .NOT. SkewEdge ) THEN
              IF( YConst ) THEN
                ! Test whether the sum of coefficients has already reached unity
                wsum = SUM( coeff(1:ncoeff) )
                EdgeReady = ( 1.0_dp - wsum < 1.0d-12 ) 
              ELSE IF( XConst ) THEN                       
                ! If edge was found both on left and right there is no need to continue search
                EdgeReady = ( dncoeff == 2 ) 
              END IF
              IF( EdgeReady ) EXIT
            END IF
          END DO

          IF( YConst ) THEN
            ! For constant y check the 2nd part 
            ! and redo the search if it is active. 
            IF( Repeating ) THEN
              IF( NRange2 /= 0 ) THEN
                x1 = x1 - NRange2 * XRange
                x2 = x2 - NRange2 * XRange
                NRange = NRange + NRange2
                NRange2 = 0
                Repeated = .TRUE.
                GOTO 100
              END IF
            END IF
          ELSE
            ! Here there can be a second part if a proper hit was not found 
            ! due to some epsilon rules.
            IF( SkewEdge ) THEN
              IF( dncoeff == 1 ) THEN
                coeff(ncoeff) = cskew * 1.0_dp
              ELSE IF( dncoeff == 2 ) THEN
                xm1 = coeff(ncoeff-1)
                xm2 = coeff(ncoeff)
                
                IF( ABS( xm2-xm1) < TINY( xm2 ) ) THEN
                  CALL Warn(Caller,'Degenerated edge 3?')
                  coeff(ncoeff-1) = cskew * 0.5_dp
                ELSE
                  coeff(ncoeff-1) = cskew * ABS((xm2-xmean)/(xm2-xm1))
                END IF
                coeff(ncoeff) = cskew * 1.0_dp - coeff(1)
              END IF
            ELSE
              IF( ncoeff == 1 ) THEN
                coeff(1) = 1.0_dp
              ELSE IF( ncoeff >= 2 ) THEN
                IF( ncoeff > 2 ) THEN
                  CALL Warn(Caller,&
                       'There should not be more than two target edges: '//TRIM(I2S(ncoeff))) 
                END IF
                xm1 = coeff(1)
                xm2 = coeff(2)
                IF( ABS( xm2-xm1) < TINY( xm2 ) ) THEN
                  CALL Warn(Caller,'Degenerated edge 3?')
                  coeff(1) = 0.5_dp
                ELSE
                  coeff(1) = ABS((xm2-xmean)/(xm2-xm1))
                END IF
                coeff(2) = 1.0_dp - coeff(1)
              END IF
            END IF

            wsum = SUM( coeff(1:ncoeff) )
          END IF

          ! Skewed edge is treated in four different parts (0,1,2,3)
          ! Go for the next part, if not finished. 
          IF( SkewEdge ) THEN
            IF( SkewPart < 3 ) THEN
              SkewPart = SkewPart + 1
              GOTO 50
            END IF
          END IF
              
          IF( ncoeff == 0 ) THEN
            Nundefined = Nundefined + 1
            WRITE( Message,'(A,2I8,4ES12.3)') 'Problematic edge: ',&
                eind,ParEnv % MyPe,x1,x2,y1,y2
            CALL Warn(Caller, Message )
            WRITE( Message,'(A,I8,3L4,4ES12.3)') 'Bounding box: ',&
                eind,XConst,YConst,Repeating,XminAll,XmaxAll,YminAll,YmaxAll
            CALL Warn(Caller, Message )
            CYCLE
          END IF

          wsum = SUM( ABS( coeff(1:ncoeff) ) )
          minwsum = MIN( minwsum, wsum ) 
          maxwsum = MAX( maxwsum, wsum ) 

          ! In skewed edges the sum of weights may be different from 1 but otherwise
          ! it should be very close to one. 
!          IF( ABS(wsum) < 0.999 .OR. ( ABS(wsum) > 1.001 .AND. .NOT. SkewEdge ) ) THEN
          IF(.FALSE.) THEN
            PRINT *,'*********************'
            PRINT *,'wsum',eind,ncoeff,wsum,Repeated
            PRINT *,'x coords:',x1,x2
            PRINT *,'y coords:',y1,y2
            PRINT *,'xm:',xm1,xm2
            PRINT *,'ym:',ym1,ym2
            PRINT *,'xm coords:',NodesM % x(1:4)
            PRINT *,'ym coords:',NodesM % y(1:4)
            PRINT *,'Const:',XConst,YConst,XConstM,YConstM
            PRINT *,'coeff:',ncoeff,coeff(1:ncoeff),coeffi(1:ncoeff)
          END IF

          ! Mark that this is set so it don't need to be set again
          EdgePerm(eind) = 0

          ! Ok, we found a true projector entry
          Projector % InvPerm(nrow) = EdgeCol0 + eind

          ! The reference to the edge to be projected
          IF( SelfProject ) THEN
            val = 1.0_dp
            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                EdgeCol0 + eind, EdgeCoeff * val ) 
          END IF

          ! The scaling can be used to create antiperiodic projectors, for example. 
          Coeff(1:ncoeff) = signs(1:ncoeff) * Coeff(1:ncoeff)

          ! And finally add the projection weights to the projection matrix
          DO j=1,ncoeff 
            val = Coeff(j)

            IF( ABS( val ) < 1.0d-12 ) CYCLE

            CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                EdgeCol0 + coeffi(j), EdgeScale * EdgeCoeff * val )
          END DO
        END DO
      END DO
         
      IF( Nundefined > 0 ) THEN
        CALL Error(Caller,&
            'Number of edges could not be mapped: '//TRIM(I2S(Nundefined)))          
      END IF

      WRITE( Message,'(A,ES12.5)') 'Minimum absolute sum of edge weights: ',minwsum
      CALL Info(Caller,Message,Level=10)
      
      WRITE( Message,'(A,ES12.5)') 'Maximum absolute sum of edge weights: ',maxwsum
      CALL Info(Caller,Message,Level=10)
      
      IF( NoSkewed > 0 ) THEN
        CALL Info(Caller,'Number of skewed edge mappings: '//TRIM(I2S(NoSkewed)),Level=8)
      END IF
      CALL Info(Caller,'Created strong constraints for edge dofs',Level=8)      

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, &
          NodesM % x, NodesM % y, NodesM % z )

    END SUBROUTINE AddEdgeProjectorStrongStrides
    !----------------------------------------------------------------------

        
    !---------------------------------------------------------------------------------
    ! Create a strong projector for edges in a conforming case.
    ! We create a periodic permutation first instead of creating a matrix directly.
    ! This enables that we can recycle some code. 
    !---------------------------------------------------------------------------------
    SUBROUTINE AddEdgeProjectorStrongConforming()

      INTEGER :: ne, nn, i, nrow, eind, eindm, sgn
      INTEGER, POINTER :: PerPerm(:)
      LOGICAL, POINTER :: PerFlip(:)
      
      CALL Info(Caller,'Creating strong projector for conforming edges',Level=8)

      ne = Mesh % NumberOfEdges
      IF( ne == 0 ) RETURN      

      nn = Mesh % NumberOfNodes            

      ALLOCATE( PerPerm(nn+ne), PerFlip(nn+ne) )
      PerPerm = 0; PerFlip = .FALSE.

      ! Permutation that tells which slave edge depends on which master edge (1-to-1 map)
      CALL ConformingEdgePerm(Mesh, BMesh1, BMesh2, PerPerm, PerFlip )
      
      DO i=nn+1,nn+ne
        IF( PerPerm(i) == 0 ) CYCLE
        eind = i - nn
        eindm = PerPerm(i) - nn
        
        sgn = -1
        IF( PerFlip(i) ) sgn = 1
        
        nrow = EdgeRow0 + EdgePerm(eind)         
        Projector % InvPerm(nrow) = EdgeCol0 + eind
        
        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
            EdgeCol0 + eind, EdgeCoeff ) 
        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
            EdgeCol0 + eindm, sgn * EdgeScale * EdgeCoeff )

        ! Mark that this is now set
        EdgePerm(eind) = 0        
      END DO
      
      DEALLOCATE( PerPerm, PerFlip ) 

      CALL Info(Caller,'Created strong constraints for conforming edge dofs',Level=10)            
      
    END SUBROUTINE AddEdgeProjectorStrongConforming

    !---------------------------------------------------------------------------------
    ! Create a strong projector for edges in a conforming case.
    ! We create a periodic permutation first instead of creating a matrix directly.
    ! This enables that we can recycle some code. 
    !---------------------------------------------------------------------------------
    SUBROUTINE AddNodeProjectorStrongConforming()

      INTEGER :: nn, i, nrow, ind, indm, sgn
      INTEGER, POINTER :: PerPerm(:)
      
      CALL Info(Caller,'Creating strong projector for conforming edges',Level=8)


      nn = Mesh % NumberOfNodes            

      ALLOCATE( PerPerm(nn) )
      PerPerm = 0

      ! Permutation that tells which slave edge depends on which master node (1-to-1 map)
      CALL ConformingNodePerm(Mesh, BMesh1, BMesh2, PerPerm )
      
      DO i=1, nn
        IF( PerPerm(i) == 0 ) CYCLE
        ind = i 
        indm = PerPerm(i) 
        
        sgn = -1
        
        nrow = NodePerm(ind)         
        Projector % InvPerm(nrow) = ind
        
        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
            ind, EdgeCoeff ) 
        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
            indm, sgn * EdgeScale * EdgeCoeff )

        ! Mark that this is now set
        NodePerm(ind) = 0        
      END DO
      
      DEALLOCATE( PerPerm )

      CALL Info(Caller,'Created strong constraints for conforming node dofs',Level=10)            
      
    END SUBROUTINE AddNodeProjectorStrongConforming

    
    !----------------------------------------------------------------------
    ! Create weak projector for the remaining nodes and edges.
    ! This uses the generic way to introduce the weights. The resulting 
    ! matrix is more dense but should be numerically favourable. 
    ! The integration is done by making an on-the-fly triangularization 
    ! into several triangles. This is not generic - it assumes constant
    ! y levels, and cartesian mesh where the search is done.  
    !----------------------------------------------------------------------
    SUBROUTINE AddProjectorWeakStrides()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: j1,j2,j3,j4,jj,ii,sgn0,k,kmax,ind,indM,nip,nn,ne,nf,inds(10),Ninteg,NintegGen
      TYPE(Element_t), POINTER :: Element, ElementM
      TYPE(Element_t) :: ElementT
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: RightSplit, LeftSplit, LeftSplit2, RightSplit2, TopEdge, BottomEdge
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: x(10),y(10),xt,yt,zt,xmax,ymax,xmin,ymin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,ArcTol,u,v,w,um,vm,wm,val,Overlap,RefArea,dArea,&
          SumOverlap,SumArea,qleft, qright, qleft2, qright2, MaxErr,Err,phi(10)
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), BasisM(:)
      REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:),WBasisM(:,:),RotWbasis(:,:),dBasisdx(:,:)
      LOGICAL :: LeftCircle, Stat
      TYPE(Mesh_t), POINTER :: Mesh

      CALL Info(Caller,'Creating weak projector for stride mesh',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      ALLOCATE( NodesT % x(n), NodesT % y(n), NodesT % z(n) )
      ALLOCATE( Basis(n), BasisM(n) )
      ALLOCATE( dBasisdx(n,3), WBasis(n,3), WBasisM(n,3), RotWBasis(n,3) )

      Nodes % z  = 0.0_dp
      NodesM % z = 0.0_dp
      NodesT % z = 0.0_dp

      MaxErr = 0.0_dp
      zt = 0.0_dp
      n = 4
      LeftCircle = .FALSE.

      ArcTol = ArcCoeff * Xtol     
      Ninteg = 0
      NintegGen = 0

      ! The temporal triangle used in the numerical integration
      ElementT % TYPE => GetElementType( 303, .FALSE. )
      ElementT % NodeIndexes => IndexesT

      DO ind=1,BMesh1 % NumberOfBulkElements

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        ne = Element % TYPE % NumberOfEdges
        IF( PiolaVersion ) THEN
          nf = 2
        ELSE
          nf = 0
        END IF
        
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))

        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))
        ymin = MINVAL(Nodes % y(1:n))
        ymax = MAXVAL(Nodes % y(1:n))

        IF( Repeating ) THEN
          Nrange = FLOOR( (xmin-XMinAll) / XRange )
          xmin = xmin - Nrange * XRange
          xmax = xmax - Nrange * XRange
          Nodes % x(1:n) = Nodes % x(1:n) - NRange * XRange 
          IF( xmax > XMaxAll ) THEN
            Nrange2 = 1
          ELSE IF( xmax < XMinAll ) THEN
            Nrange2 = -1
          ELSE
            Nrange2 = 0
          END IF
        ELSE IF( FullCircle ) THEN
          LeftCircle = ( ALL( ABS( Nodes % x(1:n) ) > 90.0_dp ) )
          IF( LeftCircle ) THEN
            DO j=1,n
              IF( Nodes % x(j) < 0.0 ) Nodes % x(j) = Nodes % x(j) + 360.0_dp
            END DO
          END IF
        END IF

        ! Transform the angle to archlength in order to have correct mapping 
        ! of skewed edges.
        Nodes % x(1:n) = ArcCoeff * Nodes % x(1:n)
        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))

        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;
        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        IP = GaussPoints( Element ) 
        RefArea = detJ * SUM( IP % s(1:IP % n) )

        SumArea = 0.0_dp
        SumOverlap = 0.0_dp
        
200     sgn0 = 1
        IF( AntiRepeating ) THEN
          IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
        END IF

        ! find an index offset such that [j1,j2,j3,j4] is ordered the as the standard
        ! nodes in bilinear elements. This could be made generic as well, but it was
        ! easier for me to fix these indexes in this way and I was feeling lazy. 
        j1 = 1; j2 = 1; j3 = 1; j4 = 1
        DO j=2,4
          ! Lower left
          IF( Nodes % x(j) + Nodes % y(j) < Nodes % x(j1) + Nodes % y(j1) ) j1 = j
          ! Lower right
          IF( Nodes % x(j) - Nodes % y(j) > Nodes % x(j2) - Nodes % y(j2) ) j2 = j
          ! Upper right
          IF( Nodes % x(j) + Nodes % y(j) > Nodes % x(j3) + Nodes % y(j3) ) j3 = j
          ! Upper left
          IF( Nodes % x(j) - Nodes % y(j) < Nodes % x(j4) - Nodes % y(j4) ) j4 = j
        END DO

        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        DO indM=1,BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes
          
          NodesM % y(1:n) = BMesh2 % Nodes % y(IndexesM(1:n))
          
          ! Make the quick and dirty search first
          yminm = MINVAL( NodesM % y(1:n))
          IF( ABS( ymin - yminm ) > YTol ) CYCLE
          
          ymaxm = MAXVAL( NodesM % y(1:n))
          IF( ABS( ymax - ymaxm ) > YTol ) CYCLE
          
          NodesM % x(1:n) = BMesh2 % Nodes % x(IndexesM(1:n))

          ! Treat the left circle differently. 
          IF( LeftCircle ) THEN
            ! Omit the element if it is definitely on the right circle
            IF( ALL( ABS( NodesM % x(1:n) ) - 90.0_dp < Xtol ) ) CYCLE
            DO j=1,n
              IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = NodesM % x(j) + 360.0_dp
            END DO
          END IF
          
          ! Transfer into real length units instead of angles
          ! This gives right balance between x and y -directions. 
          NodesM % x(1:n) = ArcCoeff * NodesM % x(1:n)

          xminm = MINVAL( NodesM % x(1:n))
          xmaxm = MAXVAL( NodesM % x(1:n))
                    
          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > ArcCoeff * 180.0_dp ) CYCLE
          END IF
          
          Overlap = (MIN(xmax, xmaxm)- MAX(xmin,xminm))/(xmax-xmin)
          IF( Overlap < RelTolX ) CYCLE 
          
          SumOverlap = SumOverlap + Overlap
          Ninteg = Ninteg + 1
          
          ! Then if this is a possible element create a list of the corner nodes
          ! for a temporal mesh. There will be 3 to 6 corner nodes. 
          ! Check the crossings between the edges of the quadrilaters. These will
          ! be used as new points when creating the virtual triangle mesh. 
          LeftSplit = ( ( Nodes % x(j1) - xminm ) * ( xminm - Nodes % x(j4) ) > 0.0_dp )
          IF(LeftSplit) qleft =  ( Nodes % x(j1) - xminm ) / ( Nodes % x(j1) - Nodes % x(j4) )

          RightSplit = ( ( Nodes % x(j2) - xmaxm ) * ( xmaxm - Nodes % x(j3) ) > 0.0_dp )
          IF(RightSplit) qright = ( Nodes % x(j2) - xmaxm ) / ( Nodes % x(j2) - Nodes % x(j3) )

          LeftSplit2 = ( ( Nodes % x(j2) - xminm ) * ( xminm - Nodes % x(j3) ) > 0.0_dp )
          IF(LeftSplit2) qleft2 =  ( Nodes % x(j2) - xminm ) / ( Nodes % x(j2) - Nodes % x(j3) )

          RightSplit2 = ( ( Nodes % x(j1) - xmaxm ) * ( xmaxm - Nodes % x(j4) ) > 0.0_dp )
          IF(RightSplit2) qright2 = ( Nodes % x(j1) - xmaxm ) / ( Nodes % x(j1) - Nodes % x(j4) )

            ! Mark the splits on the vertical edges aligned with the y-axis
            k = 0
            IF( LeftSplit ) THEN
              k = k + 1
              x(k) = xminm
              qleft = MAX( 0.0, MIN( 1.0, qleft ) )
              y(k) = Nodes % y(j1) + qleft * ( Nodes % y(j4) - Nodes % y(j1))
            END IF
            IF( RightSplit2 ) THEN
              k = k + 1
              x(k) = xmaxm
              qright2 = MAX( 0.0, MIN( 1.0, qright2 ) )
              y(k) = Nodes % y(j1) + qright2 * ( Nodes % y(j4) - Nodes % y(j1))
            END IF
            IF( RightSplit ) THEN
              k = k + 1
              x(k) = xmaxm
              qright = MAX( 0.0, MIN( 1.0, qright ) )
              y(k) = Nodes % y(j2) + qright * ( Nodes % y(j3) - Nodes % y(j2))
            END IF
            IF( LeftSplit2 ) THEN
              k = k + 1
              x(k) = xminm
              qleft2 = MAX( 0.0, MIN( 1.0, qleft2 ) )
              y(k) = Nodes % y(j2) + qleft2 * ( Nodes % y(j3) - Nodes % y(j2))
            END IF

            ! Mark the splits on the horizontal axis
            BottomEdge = .NOT. ( ( Nodes % x(j2) < xminm ) .OR. ( Nodes % x(j1) > xmaxm ) )
            TopEdge    = .NOT. ( ( Nodes % x(j3) < xminm ) .OR. ( Nodes % x(j4) > xmaxm ) )

            IF( BottomEdge ) THEN
              k = k + 1
              x(k) = MAX( xminm, Nodes % x(j1) )
              y(k) = yminm
              k = k + 1
              x(k) = MIN( xmaxm, Nodes % x(j2) )
              y(k) = yminm
            END IF
            IF( TopEdge ) THEN
              k = k + 1
              x(k) = MIN( xmaxm, Nodes % x(j3) )
              y(k) = ymaxm
              k = k + 1
              x(k) = MAX( xminm, Nodes % x(j4) )
              y(k) = ymaxm
            END IF
            kmax = k 

            IF( kmax < 3 ) THEN
              CALL Warn(Caller,'Cannot integrate over '//TRIM(I2S(kmax))//' nodes')
              CYCLE
            END IF
            
            ! The polygon is convex and hence its center lies inside the polygon
            xt = SUM(x(1:kmax)) / kmax
            yt = SUM(y(1:kmax)) / kmax

            ! Set the angle from the center and order the nodes so that they 
            ! can be easily triangulated.
            DO k=1,kmax
              phi(k) = ATAN2( y(k)-yt, x(k)-xt )
              inds(k) = k
            END DO
                       
            CALL SortR(kmax,inds,phi)
            x(1:kmax) = x(inds(1:kmax))
            y(1:kmax) = y(inds(1:kmax))
            !PRINT *,'Polygon: ',ind,indm,LeftSplit, RightSplit, LeftSplit2, RightSplit2, TopEdge, BottomEdge, kmax 

          ! Deal the case with multiple corners by making 
          ! triangulariation using one corner point.
          ! This should be ok as the polygon is always convex.
          NodesT % x(1) = x(1)
          NodesT % y(1) = y(1)

          ! Use somewhat higher integration rules than the default
          IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2 ) 

          DO k=1,kmax-2                         

            ! This check over area also automatically elimiates redundant nodes
            ! that were detected twice.
            dArea = 0.5_dp*ABS( (x(k+1)-x(1))*(y(k+2)-y(1)) -(x(k+2)-x(1))*(y(k+1)-y(1)))
            IF( dArea < RelTolY**2 * RefArea ) CYCLE

            NodesT % x(2) = x(k+1)
            NodesT % y(2) = y(k+1)
            NodesT % x(3) = x(k+2)
            NodesT % y(3) = y(k+2)
            
            ! Integration over the temporal element
            DO nip=1, IP % n 
              stat = ElementInfo( ElementT,NodesT,IP % u(nip),IP % v(nip),IP % w(nip),detJ,Basis)
              
              ! We will actually only use the global coordinates and the integration weight 
              ! from the temporal mesh. 

              ! Global coordinates of the integration point
              xt = SUM( Basis(1:3) * NodesT % x(1:3) )
              yt = SUM( Basis(1:3) * NodesT % y(1:3) )
              zt = 0.0_dp

              ! Integration weight for current integration point
              Wtemp = DetJ * IP % s(nip)
              sumarea = sumarea + Wtemp
              
              ! Integration point at the slave element
              CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx,EdgeBasis=WBasis)
                ELSE
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx )
                  CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
              END IF

              ! Integration point at the master element
              CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, Basis, dBasisdx, EdgeBasis=WBasisM)
                ELSE
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, BasisM, dBasisdx )
                  CALL GetEdgeBasis(ElementM,WBasisM,RotWBasis,BasisM,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
              END IF

              ! Add the nodal dofs
              IF( DoNodes .AND. .NOT. StrongNodes ) THEN
                DO j=1,n 
                  jj = Indexes(j)                                    
                  nrow = NodePerm( InvPerm1(jj) )
                  IF( nrow == 0 ) CYCLE

                  Projector % InvPerm(nrow) = InvPerm1(jj)
                  val = Basis(j) * Wtemp
                  DO i=1,n
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                        InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val ) 

                    IF( ABS( val * BasisM(i) ) < 1.0d-10 ) CYCLE
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                        InvPerm2(IndexesM(i)), -NodeScale * NodeCoeff * BasisM(i) * val )   
                  END DO
                END DO
              END IF

              IF( DoEdges ) THEN
                ! Dofs are numbered as follows:
                ! 1....number of nodes
                ! + ( 1 ... number of edges )
                ! + ( 1 ... 2 x number of faces )
                !-------------------------------------------
                DO j=1,ne+nf
                  
                  IF( j <= ne ) THEN
                    jj = Element % EdgeIndexes(j) 
                    IF( EdgePerm(jj) == 0 ) CYCLE
                    nrow = EdgeRow0 + EdgePerm(jj)
                    jj = jj + EdgeCol0
                    Projector % InvPerm( nrow ) = jj
                  ELSE
                    jj = 2 * ( ind - 1 ) + ( j - 4 )
                    nrow = FaceRow0 + jj
                    jj = 2 * ( Element % ElementIndex - 1) + ( j - 4 ) 
                    Projector % InvPerm( nrow ) = FaceCol0 + jj
                  END IF
                                   
                  DO i=1,ne+nf
                    IF( i <= ne ) THEN
                      ii = Element % EdgeIndexes(i) + EdgeCol0
                    ELSE
                      ii = 2 * ( Element % ElementIndex - 1 ) + ( i - 4 ) + FaceCol0
                    END IF
                    val = Wtemp * SUM( WBasis(j,:) * Wbasis(i,:) ) 
                    IF( ABS( val ) > 1.0d-12 ) THEN
                      CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          ii, EdgeCoeff * val ) 
                    END IF

                    IF( i <= ne ) THEN
                      ii = ElementM % EdgeIndexes(i) + EdgeCol0
                    ELSE
                      ii = 2 * ( ElementM % ElementIndex - 1 ) + ( i - 4 ) + FaceCol0
                    END IF                    
                    val = -Wtemp * SUM( WBasis(j,:) * WBasisM(i,:) ) 
                    IF( ABS( val ) > 1.0d-12 ) THEN
                      CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          ii, EdgeScale * EdgeCoeff * val  ) 
                    END IF
                  END DO
                END DO
              END IF
            END DO
          END DO
        END DO
        
        IF( Repeating ) THEN
          IF( NRange2 /= 0 ) THEN
            xmin = xmin - ArcCoeff * Nrange2 * XRange
            xmax = xmax - ArcCoeff * Nrange2 * XRange
            Nodes % x(1:n) = Nodes % x(1:n) - ArcCoeff * NRange2 * XRange 
            NRange = NRange + NRange2
            NRange2 = 0
            GOTO 200
          END IF
        END IF

        Err = SumArea/RefArea-1.0_dp
        MaxErr = MAX( MaxErr,ABS(Err))
      END DO

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )
      DEALLOCATE( NodesT % x, NodesT % y, NodesT % z )
      DEALLOCATE( Basis, BasisM )
      DEALLOCATE( dBasisdx, WBasis, WBasisM, RotWBasis )

      CALL Info(Caller,'Number of integration pairs: '&
          //TRIM(I2S(Ninteg)),Level=10)

      WRITE( Message,'(A,ES12.3)') 'Maximum error in area integration:',MaxErr 
      CALL Info(Caller,Message,Level=8)


    END SUBROUTINE AddProjectorWeakStrides


    SUBROUTINE LocalEdgeSolutionCoeffs( BC, Element, Nodes, ne, nf, PiolaVersion, SecondOrder, &
        dim, cFact )
      TYPE(ValueList_t), POINTER :: BC
      TYPE(Element_t), POINTER :: Element
      TYPE(Nodes_t) :: Nodes
      INTEGER :: ne, nf, dim
      LOGICAL :: PiolaVersion, SecondOrder            
      REAL(KIND=dp) :: cFact(:)

      TYPE(GaussIntegrationPoints_t) :: IP
      INTEGER :: i,j,m,nip,AllocStat
      REAL(KIND=dp) :: u,v,w,uq,vq,CMass(6,6),CForce(6),detJ,wtemp
      REAL(KIND=dp), POINTER, SAVE :: Basis(:),WBasis(:,:),RotWBasis(:,:), &
          dBasisdx(:,:)
      LOGICAL :: stat, Visited = .FALSE.
      REAL(KIND=dp) :: cvec(2)
      REAL(KIND=dp), POINTER :: pCvec(:,:)
       
      SAVE Visited, cVec 
      
      
      IF( .NOT. Visited ) THEN
        m = 12 
        ALLOCATE( Basis(m), WBasis(m,3), RotWBasis(m,3), dBasisdx(m,3), STAT=AllocStat )
        IF( AllocStat /= 0 ) CALL Fatal('AddProjectorWeakGeneric','Allocation error 3')
        
        pCvec => ListGetConstRealArray( BC,'Level Projector Debug Vector',Found)
        IF( Found ) THEN                  
          Cvec(1:2) = pCvec(1:2,1)
        ELSE
          Cvec = 1.0_dp
        END IF
        Visited = .TRUE.
      END IF

          
      IP = GaussPoints( Element ) 
      CMass = 0.0_dp
      cForce = 0.0_dp                   
      m = ne + nf
      
      DO nip=1, IP % n 
        u = IP % u(nip)
        v = IP % v(nip)
        w = 0.0_dp

        IF (PiolaVersion) THEN
          ! Take into account that the reference elements are different:
          IF ( ne == 3) THEN
            uq = u
            vq = v
            u = -1.0d0 + 2.0d0*uq + vq
            v = SQRT(3.0d0)*vq
          END IF
          IF (SecondOrder) THEN
            stat = EdgeElementInfo( Element, Nodes, u, v, w, &
                DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
                BasisDegree = 2, ApplyPiolaTransform = .TRUE.)
          ELSE
            stat = ElementInfo( Element, Nodes, u, v, w, &
                detJ, Basis, dBasisdx, EdgeBasis=WBasis)
          END IF
        ELSE
          stat = ElementInfo( Element, Nodes, u, v, w, &
              detJ, Basis, dBasisdx )
          CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)              
        END IF

        wtemp = detJ * IP % s(nip)
        DO i=1,m
          DO j=1,m
            CMASS(i,j) = CMASS(i,j) + wtemp * SUM( WBasis(i,1:dim) * WBasis(j,1:dim) )
          END DO
          CFORCE(i) = CFORCE(i) + wtemp * SUM( WBasis(i,1:dim) * cVec(1:dim) )
        END DO
      END DO
      CALL LUSolve(m, CMass(1:m,1:m), cForce(1:m) )
      cFact(1:m) = cForce(1:m)                    
      
    END SUBROUTINE LocalEdgeSolutionCoeffs
    

  
    !----------------------------------------------------------------------
    ! Create weak projector for the remaining nodes and edges
    ! using generic algo that can deal with triangles and quadrilaterals.
    !----------------------------------------------------------------------
    SUBROUTINE AddProjectorWeakGeneric()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: jj,ii,sgn0,k,kmax,ind,indM,nip,nn,ne,nf,inds(10),nM,neM,nfM,iM,i2,i2M
      INTEGER :: edge, edof, fdof
      INTEGER :: ElemCands, TotCands, ElemHits, TotHits, EdgeHits, CornerHits, &
          MaxErrInd, MinErrInd, InitialHits, ActiveHits, TimeStep, Nrange1, NoGaussPoints, &
          Centeri, CenteriM, CenterJ, CenterJM, AllocStat, NrangeAve
      TYPE(Element_t), POINTER :: Element, ElementM, ElementP
      INTEGER :: ElemCode, LinCode, ElemCodeM, LinCodeM
      TYPE(Element_t) :: ElementT
      TYPE(Element_t), TARGET :: ElementLin
      TYPE(GaussIntegrationPoints_t) :: IP
      LOGICAL :: RightSplit, LeftSplit, LeftSplit2, RightSplit2, TopEdge, BottomEdge
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: x(10),y(10),xt,yt,zt,xmax,ymax,xmin,ymin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,ArcTol,u,v,w,um,vm,wm,val,RefArea,dArea,&
          SumArea,MaxErr,MinErr,Err,phi(10),Point(3),uvw(3),ArcRange , &
          val_dual, zmin, zmax, zminm, zmaxm, dAlpha, uq, vq
      REAL(KIND=dp) :: A(2,2), B(2), C(2), absA, detA, rlen, &
          x1, x2, y1, y2, x1M, x2M, y1M, y2M, x0, y0, dist, DistTol, &
          amin, amax, aminM, amaxM, rmin2, rmax2, rmin2M, rmax2M
      REAL(KIND=dp) :: TotRefArea, TotSumArea, Area
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), BasisM(:)
      REAL(KIND=dp), POINTER :: Alpha(:), AlphaM(:)
      REAL(KIND=dp), ALLOCATABLE :: WBasis(:,:),WBasisM(:,:),RotWbasis(:,:),dBasisdx(:,:)
      LOGICAL :: LeftCircle, Stat, CornerFound(4), CornerFoundM(4), PosAngle
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: TimestepVar

      ! These are used temporarily for debugging purposes
      INTEGER :: SaveInd, MaxSubElem, MaxSubTriangles, DebugInd, Nslave, Nmaster
      LOGICAL :: SaveElem, DebugElem, SaveErr, DebugEdge
      REAL(KIND=dp) :: sums, summ, summ2, summabs, EdgeProj(2), EdgeProjM(2), ci, &
          EdgeErr, MaxEdgeErr, cFact(6),cFactM(6)
      CHARACTER(LEN=20) :: FileName
      REAL(KIND=dp), ALLOCATABLE :: CoeffBasis(:), MASS(:,:)

      
      CALL Info(Caller,'Creating weak constraints using a generic integrator',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      SaveInd = ListGetInteger( BC,'Level Projector Save Element Index',Found )
      DebugInd = ListGetInteger( BC,'Level Projector Debug Element Index',Found )
      SaveErr = ListGetLogical( BC,'Level Projector Save Fraction',Found)
      DebugEdge = ListGetLogical( BC,'Level Projector Debug Edge',Found )
      
      TimestepVar => VariableGet( Mesh % Variables,'Timestep',ThisOnly=.TRUE. )
      Timestep = NINT( TimestepVar % Values(1) )

      IF( SaveErr ) THEN
        FileName = 'frac_'//TRIM(I2S(TimeStep))//'.dat'
        OPEN( 11,FILE=Filename)
      END IF
     
      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n), &
          NodesM % x(n), NodesM % y(n), NodesM % z(n), &
          NodesT % x(n), NodesT % y(n), NodesT % z(n), & 
          Basis(n), BasisM(n), dBasisdx(n,3), STAT = AllocStat )
      IF( AllocStat /= 0 ) CALL Fatal('AddProjectorWeakGeneric','Allocation error 1')
      
      IF( Naxial > 1 ) THEN
        ALLOCATE( Alpha(n), AlphaM(n) )
      ELSE
        Alpha => Nodes % x
        AlphaM => NodesM % x
      END IF

      IF(BiOrthogonalBasis) THEN
        ALLOCATE(CoeffBasis(n), MASS(n,n), STAT=AllocStat)
        IF( AllocStat /= 0 ) CALL Fatal('AddProjectorWeakGeneric','Allocation error 2')        
      END IF
        
      IF( EdgeBasis ) THEN 
        n = 12 ! Hard-coded size sufficient for second-order edge elements
        ALLOCATE( WBasis(n,3), WBasisM(n,3), RotWBasis(n,3), STAT=AllocStat )
        IF( AllocStat /= 0 ) CALL Fatal('AddProjectorWeakGeneric','Allocation error 3')
      END IF
        
      Nodes % z  = 0.0_dp
      NodesM % z = 0.0_dp
      NodesT % z = 0.0_dp

      MaxErr = 0.0_dp
      MinErr = HUGE( MinErr )
      MaxErrInd = 0
      MinErrInd = 0
      zt = 0.0_dp
      LeftCircle = .FALSE.

      ArcTol = ArcCoeff * Xtol
      ArcRange = ArcCoeff * Xrange 
     
      DistTol = ArcTol**2 + YTol**2

      ! The temporal triangle used in the numerical integration
      ElementT % TYPE => GetElementType( 303, .FALSE. )
      ElementT % NodeIndexes => IndexesT
      TotCands = 0
      TotHits = 0
      EdgeHits = 0
      CornerHits = 0
      InitialHits = 0
      ActiveHits = 0
      TotRefArea = 0.0_dp
      TotSumArea = 0.0_dp
      Point = 0.0_dp
      MaxSubTriangles = 0
      Nslave = 0
      Nmaster = 0

      IF( DebugEdge ) THEN        
        sums = 0.0_dp; summ = 0.0_dp; summ2 = 0.0_dp; summabs = 0.0_dp
        MaxEdgeErr = 0.0_dp
      END IF
      
      ! Identify center nodes for axial projectors since at the origin the angle
      ! is impossible to determine. Instead for the origin the angle is the average
      ! of the other angles in the element.
      CenterI = 0
      CenterIM = 0
      CenterJ = 0
      CenterJM = 0
      IF( Naxial > 1 ) THEN
        DO i=1,BMesh1 % NumberOfNodes
          IF( BMesh1 % Nodes % x(i)**2 + BMesh1 % Nodes % y(i)**2 < 1.0d-20 ) THEN
            CenterI = i
            CALL Info(Caller,'Found center node in slave: '&
                //TRIM(I2S(CenterI)),Level=10)
            EXIT
          END IF
        END DO
        DO i=1,BMesh2 % NumberOfNodes
          IF( BMesh2 % Nodes % x(i)**2 + BMesh2 % Nodes % y(i)**2 < 1.0d-20 ) THEN
            CenterIM = i
            CALL Info(Caller,'Found center node in master: '&
                //TRIM(I2S(CenterI)),Level=10)
            EXIT
          END IF
        END DO
      END IF
        
              
      DO ind=1,BMesh1 % NumberOfBulkElements

        ! Optionally save the submesh for specified element, for vizualization and debugging
        SaveElem = ( SaveInd == ind )
        DebugElem = ( DebugInd == ind )

        IF( DebugElem ) THEN
          PRINT *,'Debug element turned on:',ind
        END IF

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes

        n = Element % TYPE % NumberOfNodes
        ! We use 'ne' also to indicate number of corners since for triangles and quads these are the same
        ne = Element % TYPE % NumberOfEdges  ! #(SLAVE EDGES)
        nf = Element % BDOFs                 ! #(SLAVE FACE DOFS)

        ElemCode = Element % TYPE % ElementCode 
        LinCode = 101 * ne

        ! Transform the angle to archlength in order to have correct balance between x and y
        Nodes % x(1:n) = ArcCoeff * BMesh1 % Nodes % x(Indexes(1:n))
        Nodes % y(1:n) = BMesh1 % Nodes % y(Indexes(1:n))
        
        ! For axial projector the angle is neither of the coordinates
        IF( Naxial > 1 ) THEN
          ! Calculate the [min,max] range of radius squared for slave element.
          ! We are working with squares because squareroot is a relatively expensive operation. 
          rmax2 = 0.0_dp
          DO j=1,ne
            val = Nodes % x(j)**2 + Nodes % y(j)**2 
            rmax2 = MAX( rmax2, val )
          END DO

          ! The minimum distance in (r,phi) system is not simply minimum of r
          ! We have to find minimum between (0,0) and the line passing (x1,y1) and (x2,y2) 
          rmin2 = HUGE( rmin2 )
          DO j=1,ne
            k = j+1
            IF( k > ne ) k = 1
            val = SegmentOriginDistance2( Nodes % x(j), Nodes % y(j), &
                Nodes % x(k), Nodes % y(k) )
            rmin2 = MIN( rmin2, val )
          END DO

          ! Calculate the angle, and its [-180,180] range
          DO j=1,ne
            alpha(j) = ( 180.0_dp / PI ) * ATAN2( Nodes % y(j), Nodes % x(j)  ) 
          END DO

          ! If we have origin replace it with the average           
          IF( CenterI > 0 ) THEN
            CenterJ = 0
            DO j=1,ne
              IF( Indexes(j) == CenterI ) THEN
                alpha(j) = 0.0_dp
                alpha(j) = SUM( Alpha(1:ne) ) / ( ne - 1 ) 
                CenterJ = j
                EXIT
              END IF
            END DO
          END IF
            
          amin = MINVAL( Alpha(1:ne) )
          amax = MAXVAL( Alpha(1:ne) )
          IF( amax - amin < 180.0_dp ) THEN
            PosAngle = .FALSE.
          ELSE
            PosAngle = .TRUE.
            ! Map the angle to [0,360]
            DO j=1,ne
              IF( Alpha(j) < 0.0 ) Alpha(j) = Alpha(j) + 360.0_dp
            END DO
            IF( CenterJ > 0 ) THEN
              alpha(CenterJ) = 0.0_dp
              alpha(CenterJ) = SUM( Alpha(1:ne) ) / ( ne - 1 ) 
            END IF
            amin = MINVAL( Alpha(1:ne) )
            amax = MAXVAL( Alpha(1:ne) )                        
          END IF
        END IF ! Naxial > 1

        ! If we have full angle eliminate the discontinuity of the angle
        ! since we like to do the mapping using continuous coordinates.
        IF( FullCircle ) THEN
          LeftCircle = ( ALL( ABS( Alpha(1:ne) ) > ArcCoeff * 90.0_dp ) )
          IF( LeftCircle ) THEN
            DO j=1,n
              IF( Alpha(j) < 0.0 ) Alpha(j) = Alpha(j) + ArcCoeff * 360.0_dp
            END DO
          END IF
        END IF
        
        ! Even for quadratic elements only work with corner nodes (n >= ne)        
        xmin = MINVAL(Nodes % x(1:ne))
        xmax = MAXVAL(Nodes % x(1:ne))

        ymin = MINVAL(Nodes % y(1:ne))
        ymax = MAXVAL(Nodes % y(1:ne))
                
        IF( HaveMaxDistance ) THEN
          zmin = MINVAL( BMesh1 % Nodes % z(Indexes(1:ne)) )
          zmax = MAXVAL( BMesh1 % Nodes % z(Indexes(1:ne)) )
        END IF
        
        IF( DebugEdge ) THEN
          CALL LocalEdgeSolutionCoeffs( BC, Element, Nodes, ne, nf, &
              PiolaVersion, SecondOrder, 2, cFact )
          EdgeProj = 0.0_dp; EdgeProjM = 0.0_dp
        END IF
        
        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;

        IF( DebugElem ) THEN
          PRINT *,'inds',n,ne,LinCode,ElemCode
          PRINT *,'x:',Nodes % x(1:n)
          PRINT *,'y:',Nodes % y(1:n)
          PRINT *,'z:',Nodes % z(1:n)
          PRINT *,'xrange:',xmin,xmax
          PRINT *,'yrange:',ymin,ymax
          PRINT *,'zrange:',zmin,zmax
          IF( Naxial > 1 ) PRINT *,'Alpha: ',Alpha(1:n)
        END IF


        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

        IP = GaussPoints( Element ) 
        RefArea = detJ * SUM( IP % s(1:IP % n) )
        SumArea = 0.0_dp

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_a.dat'
          OPEN( 10,FILE=Filename)
          DO i=1,ne
            WRITE( 10, * ) Nodes % x(i), Nodes % y(i)
          END DO
          CLOSE( 10 )
        END IF
        
        IF( DebugElem ) THEN
          PRINT *,'RefArea:',RefArea,detJ
          PRINT *,'Basis:',Basis(1:n)
        END IF


        IF( DoNodes .AND. .NOT. StrongNodes ) THEN
          DO i=1,n
            j = InvPerm1(Indexes(i))
            nrow = NodePerm(j)
            IF( nrow == 0 ) CYCLE
            CALL List_AddMatrixIndex(Projector % ListMatrix, nrow, j ) 
             IF(ASSOCIATED(Projector % Child)) &
               CALL List_AddMatrixIndex(Projector % Child % ListMatrix, nrow, j ) 
          END DO
        END IF


        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        ElemCands = 0
        ElemHits = 0

        
        DO indM=1,BMesh2 % NumberOfBulkElements

          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes

          nM = ElementM % TYPE % NumberOfNodes
          neM = ElementM % TYPE % ElementCode / 100

          ElemCodeM = Element % TYPE % ElementCode 
          LinCodeM = 101 * neM
            
          IF( DebugElem ) THEN
            PRINT *,'Candidate Elem:',indM,nM,NeM, ElemCodeM,LinCodeM
          END IF
 
          IF( HaveMaxDistance ) THEN
            zminm = MINVAL( BMesh2 % Nodes % z(IndexesM(1:neM)) )
            zmaxm = MINVAL( BMesh2 % Nodes % z(IndexesM(1:neM)) )
            IF( zmaxm < zmin - MaxDistance ) CYCLE
            IF( zminm > zmax + MaxDistance ) CYCLE
          END IF
          
          NodesM % y(1:nM) = BMesh2 % Nodes % y(IndexesM(1:nM))
        
          ! Make the quick and dirty search first
          ! This requires some minimal width of the cut
          IF(Naxial <= 1 ) THEN
            yminm = MINVAL( NodesM % y(1:neM))
            IF( yminm > ymax ) CYCLE
            
            ymaxm = MAXVAL( NodesM % y(1:neM))
            IF( ymaxm < ymin ) CYCLE

            NodesM % x(1:nM) = ArcCoeff * BMesh2 % Nodes % x(IndexesM(1:nM))
          ELSE
            NodesM % x(1:nM) = ArcCoeff * BMesh2 % Nodes % x(IndexesM(1:nM))

            ! For axial projector first check the radius since it does not have complications with
            ! periodicity and is therefore cheaper. 
            rmax2M = 0.0_dp
            DO j=1,neM
              val = NodesM % x(j)**2 + NodesM % y(j)**2 
              rmax2M = MAX( rmax2M, val )
            END DO
            IF( rmax2m < rmin2 ) CYCLE
              
            ! The minimum distance in (r,phi) system is not simply minimum of r
            ! We have to find minimum between (0,0) and the line passing (x1,y1) and (x2,y2) 
            rmin2M = HUGE( rmin2M )
            DO j=1,neM
              k = j+1
              IF( k > neM ) k = 1
              val = SegmentOriginDistance2( NodesM % x(j), NodesM % y(j), &
                  NodesM % x(k), NodesM % y(k) )
              rmin2M = MIN( rmin2M, val )
            END DO
            IF( rmin2m > rmax2 ) CYCLE
           
            ! Angle in [-180,180] or [0,360] depending where the slave angle is mapped
            DO j=1,neM
              alphaM(j) = ( 180.0_dp / PI ) * ATAN2( NodesM % y(j), NodesM % x(j)  ) 
            END DO
            
            ! If we have origin replace it with the average 
            IF( CenterIM > 0 ) THEN
              CenterJm = 0
              DO j=1,neM
                IF( IndexesM(j) == CenterIM ) THEN
                  CenterJM = j
                  alphaM(j) = 0.0_dp
                  alphaM(j) = SUM( AlphaM(1:neM) ) / ( neM - 1 ) 
                  EXIT
                END IF
              END DO
            END IF
              
            aminm = MINVAL( AlphaM(1:neM) )
            amaxm = MAXVAL( AlphaM(1:neM) )

            IF( amaxm - aminm > 180.0_dp ) THEN
              ! Map the angle to [0,360]
              DO j=1,neM
                IF( AlphaM(j) < 0.0 ) AlphaM(j) = AlphaM(j) + 360.0_dp
              END DO
              IF( CenterJM > 0 ) THEN
                alphaM(CenterJM) = 0.0_dp
                alphaM(CenterJM) = SUM( AlphaM(1:ne) ) / ( ne - 1 ) 
              END IF
              aminm = MINVAL( AlphaM(1:neM) )
              amaxm = MAXVAL( AlphaM(1:neM) )                        
            END IF
          END IF

          ! Treat the left circle differently. 
          IF( LeftCircle ) THEN
            ! Omit the element if it is definitely on the right circle
            IF( ALL( ABS( AlphaM(1:neM) ) - ArcCoeff * 90.0_dp < ArcTol ) ) CYCLE
            DO j=1,neM
              IF( AlphaM(j) < 0.0_dp ) AlphaM(j) = AlphaM(j) + ArcCoeff * 360.0_dp
            END DO
          END IF

          IF( Repeating ) THEN
            ! Enforce xmaxm to be on the same interval than xmin
            IF( Naxial > 1 ) THEN
              Nrange1 = FLOOR( Naxial * (amaxm-amin+RelTolX) / 360.0_dp )
              Nrange2 = FLOOR( Naxial * (amax-aminm+RelTolX) / 360.0_dp )
              
              ! The two ranges could have just offset of 2*PI, eliminate that
              !Nrange2 = Nrange2 + ((Nrange1 - Nrange2)/Naxial) * Naxial
              !  Nrange2 = Nrange1
              !END IF

              IF( MODULO( Nrange1 - Nrange2, Naxial ) == 0 )  THEN
                Nrange2 = Nrange1
              END IF
              
              IF( MODULO( Nrange1, Naxial) /= 0 ) THEN
                dAlpha = Nrange1 * 2.0_dp * PI / Naxial
                DO i=1,nM
                  x0 = NodesM % x(i)
                  y0 = NodesM % y(i)
                  NodesM % x(i) = COS(dAlpha) * x0 - SIN(dAlpha) * y0
                  NodesM % y(i) = SIN(dAlpha) * x0 + COS(dAlpha) * y0
                END DO
              END IF
                
              !IF( Nrange2 > Nrange1 + Naxial / 2 ) THEN
              !  Nrange2 = Nrange2 - Naxial
              !ELSE IF( Nrange2 < Nrange1 - Naxial / 2 ) THEN
              !  Nrange2 = Nrange2 + Naxial
              !END IF

              IF( DebugElem) THEN
                PRINT *,'axial:',ind,indM,amin,aminm,Nrange1,Nrange2
                PRINT *,'coord:',Nodes % x(1), Nodes % y(1), NodesM % x(1), NodesM % y(1)
                PRINT *,'Alphas:',Alpha(1:n),AlphaM(1:nM)
              END IF
              
            ELSE
              xminm = MINVAL( NodesM % x(1:nM) )
              xmaxm = MAXVAL( NodesM % x(1:nM) )

              Nrange1 = FLOOR( (xmaxm-xmin+ArcTol) / ArcRange )
              Nrange2 = FLOOR( (xmax-xminm+ArcTol) / ArcRange )
              IF( Nrange1 /= 0 ) THEN
                NodesM % x(1:nM) = NodesM % x(1:nM) - NRange1 * ArcRange 
              END IF
            END IF

            Nrange = Nrange1
          END IF

          xminm = MINVAL( NodesM % x(1:neM) )
          xmaxm = MAXVAL( NodesM % x(1:neM) )

          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > ArcCoeff * 180.0_dp ) CYCLE
          END IF

200       IF( xminm > xmax ) GOTO 100
          IF( xmaxm < xmin ) GOTO 100


          ! Rotation alters also the y-coordinate for "axial projector"
          ! Therefore this check is postponed until here.
          IF( Naxial > 1 ) THEN
            yminm = MINVAL( NodesM % y(1:nM) )
            IF( yminm > ymax ) GOTO 100
            
            ymaxm = MAXVAL( NodesM % y(1:nM))
            IF( ymaxm < ymin ) GOTO 100
          END IF

          neM = ElementM % TYPE % NumberOfEdges 
          nfM = ElementM % BDOFs

          k = 0
          ElemCands = ElemCands + 1
          CornerFound = .FALSE.
          CornerFoundM = .FALSE.

          ! Check through the nodes that are created in the intersections of any two edge
          DO i=1,ne
            x1 = Nodes % x(i)
            y1 = Nodes % y(i)
            i2 = i + 1 
            IF( i2 > ne ) i2 = 1  ! check the (ne,1) edge also
            x2 = Nodes % x(i2)
            y2 = Nodes % y(i2)

            DO iM=1,neM
              x1M = NodesM % x(iM)
              y1M = NodesM % y(iM)
              i2M = iM + 1
              IF( i2M > neM ) i2M = 1
              x2M = NodesM % x(i2M)
              y2M = NodesM % y(i2M)
              
              ! Upon solution this is tampered so it must be initialized 
              ! before each solution. 
              A(1,1) = x2 - x1
              A(2,1) = y2 - y1           
              A(1,2) = x1M - x2M
              A(2,2) = y1M - y2M

              detA = A(1,1)*A(2,2)-A(1,2)*A(2,1)
              absA = SUM(ABS(A(1,1:2))) * SUM(ABS(A(2,1:2)))
              
              ! Lines are almost parallel => no intersection possible
              ! Check the dist at the end of the line segments.
              IF(ABS(detA) < 1.0d-8 * absA + 1.0d-20 ) CYCLE

              B(1) = x1M - x1
              B(2) = y1M - y1
              
              CALL InvertMatrix( A,2 )
              C(1:2) = MATMUL(A(1:2,1:2),B(1:2))

              ! Check that the hit is within the line segment
              IF(ANY(C(1:2) < 0.0) .OR. ANY(C(1:2) > 1.0d0)) CYCLE
              
              ! We have a hit, two line segments can have only one hit
              k = k + 1
              
              x(k) = x1 + C(1) * (x2-x1)
              y(k) = y1 + C(1) * (y2-y1)

              ! If the point of intersection is at the end of a line-segment it
              ! is also a corner node.
              IF(ABS(C(1)) < 1.0d-6 ) THEN
                CornerFound(i) = .TRUE.
              ELSE IF( ABS(C(1)-1.0_dp ) < 1.0d-6 ) THEN
                CornerFound(i2) = .TRUE.
              END IF              

              IF(ABS(C(2)) < 1.0d-6 ) THEN
                CornerFoundM(iM) = .TRUE.
              ELSE IF( ABS(C(2)-1.0_dp ) < 1.0d-6 ) THEN
                CornerFoundM(i2M) = .TRUE.
              END IF
         
              EdgeHits = EdgeHits + 1
            END DO
          END DO

          IF( DebugElem ) THEN
            PRINT *,'EdgeHits:',k
          END IF

          ! Check the nodes that are one of the existing nodes i.e. corner nodes
          ! that are located inside in either element. We have to check both combinations. 
          DO i=1,ne
            ! This corner was already determined active as the end of edge 
            IF( CornerFound(i) ) CYCLE

            Point(1) = Nodes % x(i)
            IF( Point(1) < xminm - ArcTol ) CYCLE
            IF( Point(1) > xmaxm + ArcTol ) CYCLE

            Point(2) = Nodes % y(i)
            IF( Point(2) < yminm - YTol ) CYCLE
            IF( Point(2) > ymaxm + YTol ) CYCLE

            ! The edge intersections should catch the sharp hits so here we can use hard criteria
            Found = PointInElement( ElementM, NodesM, Point, uvw, LocalEps = 1.0d-8 )
            IF( Found ) THEN
              k = k + 1
              x(k) = Point(1)
              y(k) = Point(2)
              CornerHits = CornerHits + 1
            END IF
          END DO

          IF( DebugElem ) THEN
            PRINT *,'CornerHits:',k
          END IF

          ! Possible corner hits for the master element
          DO i=1,neM
            IF( CornerFoundM(i) ) CYCLE

            Point(1) = NodesM % x(i)
            IF( Point(1) < xmin - ArcTol ) CYCLE
            IF( Point(1) > xmax + ArcTol ) CYCLE

            Point(2) = NodesM % y(i)
            IF( Point(2) < ymin - YTol ) CYCLE
            IF( Point(2) > ymax + YTol ) CYCLE
         
            Found = PointInElement( Element, Nodes, Point, uvw, LocalEps = 1.0d-8 )
            IF( Found ) THEN
              k = k + 1
              x(k) = Point(1)
              y(k) = Point(2)
              CornerHits = CornerHits + 1
            END IF
          END DO

          IF( DebugElem ) THEN
            PRINT *,'CornerHitsM:',k
          END IF

          kmax = k          
          IF( kmax < 3 ) GOTO 100

          IF( DebugEdge ) THEN          
            CALL LocalEdgeSolutionCoeffs( BC, ElementM, NodesM, neM, nfM, &
                PiolaVersion, SecondOrder, 2, cFactM )
          END IF
          
          sgn0 = 1
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) sgn0 = -1
          END IF
          
          InitialHits = InitialHits + kmax

          ! The polygon is convex and hence its center lies inside the polygon
          xt = SUM(x(1:kmax)) / kmax
          yt = SUM(y(1:kmax)) / kmax
          
          ! Set the angle from the center and order the nodes so that they 
          ! can be easily triangulated.
          DO k=1,kmax
            phi(k) = ATAN2( y(k)-yt, x(k)-xt )
            inds(k) = k
          END DO

          IF( DebugElem ) THEN
            PRINT *,'Phis:',phi(1:kmax)
          END IF

          CALL SortR(kmax,inds,phi)
          x(1:kmax) = x(inds(1:kmax))
          y(1:kmax) = y(inds(1:kmax))

          ! Eliminate redundant corners from the polygon
          j = 1
          DO k=2,kmax
            dist = (x(j)-x(k))**2 + (y(j)-y(k))**2 
            IF( dist > DistTol ) THEN
              j = j + 1
              IF( j /= k ) THEN
                x(j) = x(k)
                y(j) = y(k)
              END IF
            END IF
          END DO
          kmax = j

          IF( DebugElem ) THEN
            PRINT *,'Corners:',kmax
            PRINT *,'Center:',xt,yt
          END IF

          IF( kmax < 3 ) GOTO 100

          ElemHits = ElemHits + 1
          ActiveHits = ActiveHits + kmax

          IF( kmax > MaxSubTriangles ) THEN
            MaxSubTriangles = kmax
            MaxSubElem = ind
          END IF

          IF( SaveElem ) THEN
            FileName = 't'//TRIM(I2S(TimeStep))//'_b'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) NodesM % x(i), NodesM % y(i)
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_d'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) xt, yt
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_e'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,kmax
              WRITE( 10, * ) x(i), y(i)
            END DO
            CLOSE( 10 )           
          END IF

          
          ! Deal the case with multiple corners by making 
          ! triangulariation using one corner point.
          ! This should be ok as the polygon is always convex.
          NodesT % x(1) = x(1)
          NodesT % y(1) = y(1)
          
          ! Use somewhat higher integration rules than the default
          
          NoGaussPoints = ListGetInteger( BC,'Mortar BC Gauss Points',Found ) 
          IF(.NOT. Found ) NoGaussPoints = ElementT % Type % GaussPoints2
          IP = GaussPoints( ElementT, NoGaussPoints )
            

          DO k=1,kmax-2                         
            
            ! This check over area also automatically elimiates redundant nodes
            ! that were detected twice.
            dArea = 0.5_dp*ABS( (x(k+1)-x(1))*(y(k+2)-y(1)) -(x(k+2)-x(1))*(y(k+1)-y(1)))

            IF( DebugElem ) THEN
              PRINT *,'dArea:',dArea,dArea / RefArea
            END IF

            IF( dArea < RelTolY**2 * RefArea ) CYCLE

            ! Triangle is created by keeping one corner node fixed and rotating through
            ! the other nodes. 
            NodesT % x(2) = x(k+1)
            NodesT % y(2) = y(k+1)
            NodesT % x(3) = x(k+2)
            NodesT % y(3) = y(k+2)

            IF(BiOrthogonalBasis) THEN
              MASS  = 0
              CoeffBasis = 0
              area = 0._dp
              DO nip=1, IP % n 
                stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                    IP % v(nip),IP % w(nip),detJ,Basis)
                IF(.NOT. Stat ) EXIT

                ! We will actually only use the global coordinates and the integration weight 
                ! from the temporal mesh. 
              
                ! Global coordinates of the integration point
                xt = SUM( Basis(1:3) * NodesT % x(1:3) )
                yt = SUM( Basis(1:3) * NodesT % y(1:3) )
                zt = 0.0_dp
              
                ! Integration weight for current integration point
                Wtemp = DetJ * IP % s(nip)
                area = area + wtemp
              
                ! Integration point at the slave element
                IF( ElemCode /= LinCode ) THEN
                  ElementLin % TYPE => GetElementType( LinCode, .FALSE. )
                  ElementLin % NodeIndexes => Element % NodeIndexes
                  ElementP => ElementLin
                  CALL GlobalToLocal( u, v, w, xt, yt, zt, ElementP, Nodes )
                ELSE
                  CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
                END IF

                stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
                IF(.NOT. Stat) CYCLE

                DO i=1,n
                  DO j=1,n
                    MASS(i,j) = MASS(i,j) + wTemp * Basis(i) * Basis(j)
                  END DO
                  CoeffBasis(i) = CoeffBasis(i) + wTemp * Basis(i)
                END DO
              END DO

              IF(Area<1.d-12) GOTO 300

              CALL InvertMatrix( MASS, n )

              DO i=1,n
                DO j=1,n
                  MASS(i,j) = MASS(i,j) * CoeffBasis(i)
                END DO
              END DO
            END IF
            
            ! Integration over the temporal element
            DO nip=1, IP % n 
              stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                  IP % v(nip),IP % w(nip),detJ,Basis)
              IF(.NOT. Stat) EXIT

              ! We will actually only use the global coordinates and the integration weight 
              ! from the temporal mesh. 
              
              ! Global coordinates of the integration point
              xt = SUM( Basis(1:3) * NodesT % x(1:3) )
              yt = SUM( Basis(1:3) * NodesT % y(1:3) )
              zt = 0.0_dp
              
              ! Integration weight for current integration point
              Wtemp = DetJ * IP % s(nip)
              sumarea = sumarea + Wtemp
              
              ! Integration point at the slave element
              IF( ElemCode /= LinCode ) THEN
                ElementLin % TYPE => GetElementType( LinCode, .FALSE. )
                ElementLin % NodeIndexes => Element % NodeIndexes
                ElementP => ElementLin
                CALL GlobalToLocal( u, v, w, xt, yt, zt, ElementP, Nodes )
              ELSE
                CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
              END IF


              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  ! Take into account that the reference elements are different:
                  IF ( ne == 3) THEN
                    uq = u
                    vq = v
                    u = -1.0d0 + 2.0d0*uq + vq
                    v = SQRT(3.0d0)*vq
                  END IF
                  IF (SecondOrder) THEN
                    stat = EdgeElementInfo( Element, Nodes, u, v, w, &
                        DetF = DetJ, Basis = Basis, EdgeBasis = WBasis, &
                        BasisDegree = 2, ApplyPiolaTransform = .TRUE.)
                  ELSE
                    stat = ElementInfo( Element, Nodes, u, v, w, &
                        detJ, Basis, dBasisdx,EdgeBasis=WBasis)
                  END IF
                ELSE
                  stat = ElementInfo( Element, Nodes, u, v, w, &
                      detJ, Basis, dBasisdx )
                  CALL GetEdgeBasis(Element,WBasis,RotWBasis,Basis,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
              END IF

              ! Integration point at the master element
              IF( ElemCodeM /= LinCodeM ) THEN
                ElementLin % TYPE => GetElementType( LinCodeM, .FALSE. )
                ElementLin % NodeIndexes => ElementM % NodeIndexes
                ElementP => ElementLin
                CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementP, NodesM )
              ELSE
                CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
              END IF


              IF( EdgeBasis ) THEN
                IF (PiolaVersion) THEN
                  ! Take into account that the reference elements are different:
                  IF ( neM == 3) THEN
                    uq = um
                    vq = vm
                    um = -1.0d0 + 2.0d0*uq + vq
                    vm = SQRT(3.0d0)*vq
                  END IF
                  IF (SecondOrder) THEN
                    stat = EdgeElementInfo( ElementM, NodesM, um, vm, wm, &
                        DetF=detJ, Basis=BasisM, EdgeBasis=WBasisM, &
                        BasisDegree = 2, ApplyPiolaTransform = .TRUE.)                   
                  ELSE
                    stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                        detJ, BasisM, dBasisdx, EdgeBasis=WBasisM)
                  END IF
                ELSE
                  stat = ElementInfo( ElementM, NodesM, um, vm, wm, &
                      detJ, BasisM, dBasisdx )
                  CALL GetEdgeBasis(ElementM,WBasisM,RotWBasis,BasisM,dBasisdx)
                END IF
              ELSE
                stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
              END IF
              IF(.NOT. Stat) CYCLE

              ! Add the nodal dofs
              IF( DoNodes .AND. .NOT. StrongNodes ) THEN
                IF(BiOrthogonalBasis) THEN
                  CoeffBasis = 0._dp
                  DO i=1,n
                    DO j=1,n
                      CoeffBasis(i) = CoeffBasis(i) + MASS(i,j) * Basis(j)
                    END DO
                  END DO
                END IF

                DO j=1,n 
                  jj = Indexes(j)                                    

                  nrow = NodePerm(InvPerm1(jj))
                  IF( nrow == 0 ) CYCLE

                  Projector % InvPerm(nrow) = InvPerm1(jj)
                  val = Basis(j) * Wtemp
                  IF(BiorthogonalBasis) val_dual = CoeffBasis(j) * Wtemp

                  !IF( DebugElem ) PRINT *,'Vals:',val

                  DO i=1,n
                    Nslave = Nslave + 1
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                          InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val ) 

                    IF(BiOrthogonalBasis) THEN
                      CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                            InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val_dual ) 
                    END IF
                  END DO

                  DO i=1,nM
                    IF( ABS( val * BasisM(i) ) < 1.0d-10 ) CYCLE

                    Nmaster = Nmaster + 1
                    CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                        InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val )                   

                    IF(BiOrthogonalBasis) THEN
                      IF(DualMaster.OR.DualLCoeff) THEN
                        CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                              InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val_dual ) 
                      ELSE
                        CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                              InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val ) 
                      END IF
                    END IF
                  END DO
                END DO
              END IF

              IF( DoEdges ) THEN
                IF (SecondOrder) THEN

                  DO j=1,2*ne+nf   ! for all slave dofs
                    IF (j<=2*ne) THEN
                      edge = 1+(j-1)/2    ! The edge to which the dof is associated
                      edof = j-2*(edge-1) ! The edge-wise index of the dof
                      jj = Element % EdgeIndexes(edge) 
                      IF( EdgePerm(jj) == 0 ) CYCLE
                      nrow = EdgeRow0 + 2*(EdgePerm(jj)-1) + edof  ! The row to be written
                      jj = EdgeCol0 + 2*(jj-1) + edof              ! The index of the corresponding DOF
                      Projector % InvPerm( nrow ) = jj
                    ELSE
                      IF( Parallel ) THEN
                        IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
                      END IF
                      fdof = j-2*ne ! The face-wise index of the dof
                      nrow = FaceRow0 + nf * ( ind - 1 ) + fdof
                      jj = FaceCol0 + nf * ( Element % ElementIndex - 1) + fdof
                      Projector % InvPerm( nrow ) = jj
                    END IF

                    DO i=1,2*ne+nf ! for all slave dofs
                      IF( i <= 2*ne ) THEN
                        edge = 1+(i-1)/2    ! The edge to which the dof is associated
                        edof = i-2*(edge-1) ! The edge-wise index of the dof
                        ii = EdgeCol0 + 2*(Element % EdgeIndexes(edge) - 1) + edof
                      ELSE
                        fdof = i-2*ne ! The face-wise index of the dof
                        ii = FaceCol0 + nf * ( Element % ElementIndex - 1) + fdof
                      END IF

                      val = Wtemp * SUM( WBasis(j,:) * Wbasis(i,:) ) 
                      IF( ABS( val ) > 1.0d-12 ) THEN
                        Nslave = Nslave + 1
                        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                            ii, EdgeCoeff * val ) 
                      END IF
                    END DO
                    
                    DO i=1,2*neM+nfM ! for all master dofs
                      IF( i <= 2*neM ) THEN
                        edge = 1+(i-1)/2    ! The edge to which the dof is associated
                        edof = i-2*(edge-1) ! The edge-wise index of the dof
                        ii = EdgeCol0 + 2*(ElementM % EdgeIndexes(edge) - 1) + edof
                      ELSE
                        fdof = i-2*neM ! The face-wise index of the dof
                        ii = FaceCol0 + nfM * ( ElementM % ElementIndex - 1) + fdof
                      END IF

                      val = -Wtemp * sgn0 * SUM( WBasis(j,:) * WBasisM(i,:) ) 
                      IF( ABS( val ) > 1.0d-12 ) THEN
                        Nmaster = Nmaster + 1
                        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                            ii, EdgeScale * EdgeCoeff * val  ) 
                      END IF
                    END DO
                  END DO

                ELSE
                  ! Dofs are numbered as follows:
                  ! 1....number of nodes
                  ! + ( 1 ... number of edges )
                  ! + ( 1 ... 2 x number of faces )
                  !-------------------------------------------
                  DO j=1,ne+nf

                    IF( j <= ne ) THEN
                      jj = Element % EdgeIndexes(j) 
                      IF( EdgePerm(jj) == 0 ) CYCLE
                      nrow = EdgeRow0 + EdgePerm(jj)
                      jj = jj + EdgeCol0
                      Projector % InvPerm( nrow ) = jj
                    ELSE
                      IF( Parallel ) THEN
                        IF( Element % PartIndex /= ParEnv % MyPe ) CYCLE
                      END IF

                      jj = 2 * ( ind - 1 ) + ( j - ne )
                      nrow = FaceRow0 + jj
                      jj = 2 * ( Element % ElementIndex - 1) + ( j - ne ) 
                      Projector % InvPerm( nrow ) = FaceCol0 + jj
                    END IF


                    DO i=1,ne+nf
                      IF( i <= ne ) THEN
                        ii = Element % EdgeIndexes(i) + EdgeCol0
                      ELSE
                        ii = 2 * ( Element % ElementIndex - 1 ) + ( i - ne ) + FaceCol0
                      END IF

                      IF( DebugEdge ) THEN
                        ci = cFact(i)
                        sums = sums + ci * EdgeCoeff * val                         
                        EdgeProj(1:2) = EdgeProj(1:2) + ci * Wtemp * Wbasis(i,1:2)
                      END IF
                        
                      val = Wtemp * SUM( WBasis(j,:) * Wbasis(i,:) ) 
                      IF( ABS( val ) > 1.0d-12 ) THEN
                        Nslave = Nslave + 1                          
                        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                            ii, EdgeCoeff * val ) 
                      END IF
                    END DO
                      
                    DO i=1,neM+nfM
                      IF( i <= neM ) THEN
                        ii = ElementM % EdgeIndexes(i) + EdgeCol0
                      ELSE
                        ii = 2 * ( ElementM % ElementIndex - 1 ) + ( i - neM ) + FaceCol0
                      END IF

                      IF( DebugEdge ) THEN
                        ci = cFactM(i)
                        summ = summ + ci * EdgeScale * EdgeCoeff * val
                        summabs = summabs + ABS( ci * EdgeScale * EdgeCoeff * val )                        
                        IF( NRange /= NRange1 ) THEN
                          summ2 = summ2 + ci * EdgeScale * EdgeCoeff * val
                        END IF                        
                        EdgeProjM(1:2) = EdgeProjM(1:2) + ci * Wtemp * sgn0 * WbasisM(i,1:2)
                      END IF
                        
                      val = -Wtemp * sgn0 * SUM( WBasis(j,:) * WBasisM(i,:) ) 
                      IF( ABS( val ) > 1.0d-12 ) THEN
                        Nmaster = Nmaster + 1
                        CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                            ii, EdgeScale * EdgeCoeff * val  ) 
                      END IF
                    END DO
                  END DO
                END IF
              END IF
            END DO

300         CONTINUE
            
          END DO

100       IF( Repeating ) THEN
            IF( NRange /= NRange2 ) THEN
              ! Rotate the sector to a new position for axial case
              ! Or just some up the angle in the radial/2D case
              IF( Naxial > 1 ) THEN

                IF( Nrange /= Nrange2 ) THEN
                  dAlpha = 2.0_dp * PI * (Nrange2 - Nrange ) / Naxial
                  Nrange = Nrange2
                END IF
             
                DO i=1,nM
                  x0 = NodesM % x(i)
                  y0 = NodesM % y(i)
                  NodesM % x(i) = COS(dAlpha) * x0 - SIN(dAlpha) * y0
                  NodesM % y(i) = SIN(dAlpha) * x0 + COS(dAlpha) * y0
                END DO
              ELSE
                Nrange = Nrange2
                NodesM % x(1:n) = NodesM % x(1:n) + ArcRange  * (Nrange2 - Nrange1)
              END IF
              xminm = MINVAL( NodesM % x(1:neM))
              xmaxm = MAXVAL( NodesM % x(1:neM))
              GOTO 200
            END IF
          END IF

        END DO

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_n.dat'
          OPEN( 10,FILE=Filename)
          OPEN( 10,FILE=FileName)
          WRITE( 10, * ) ElemHits 
          CLOSE( 10 )
        END IF
        
        TotCands = TotCands + ElemCands
        TotHits = TotHits + ElemHits
        TotSumArea = TotSumArea + SumArea
        TotRefArea = TotRefArea + RefArea

        Err = SumArea / RefArea
        IF( Err > MaxErr ) THEN
          MaxErr = Err
          MaxErrInd = Err
        END IF
        IF( Err < MinErr ) THEN
          MinErr = Err
          MinErrInd = ind
        END IF

        IF( SaveErr ) THEN
          WRITE( 11, * ) ind,SUM( Nodes % x(1:ne))/ne, SUM( Nodes % y(1:ne))/ne, Err
        END IF

        IF( DebugEdge ) THEN        
          EdgeErr = SUM( ABS( EdgeProj-EdgeProjM) ) / SUM( ABS(EdgeProj)+ABS(EdgeProjM) )          
          IF( EdgeErr > 1.0e-3 ) THEN
            PRINT *,'EdgeProj:',ind,EdgeErr,EdgeProj,EdgeProjM          
          END IF
          MaxEdgeErr = MAX( MaxEdgeErr, EdgeErr ) 
        END IF
        
      END DO

      IF( SaveErr ) CLOSE(11)
      
      
      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z, &
          NodesM % x, NodesM % y, NodesM % z, &
          NodesT % x, NodesT % y, NodesT % z, &
          Basis, BasisM, dBasisdx )
      IF( EdgeBasis ) THEN
        DEALLOCATE( WBasis, WBasisM, RotWBasis )
      END IF
      IF(BiOrthogonalBasis) THEN
        DEALLOCATE(CoeffBasis, MASS )
      END IF
       
      CALL Info(Caller,'Number of integration pair candidates: '&
          //TRIM(I2S(TotCands)),Level=10)
      CALL Info(Caller,'Number of integration pairs: '&
          //TRIM(I2S(TotHits)),Level=10)

      CALL Info(Caller,'Number of edge intersections: '&
          //TRIM(I2S(EdgeHits)),Level=10)
      CALL Info(Caller,'Number of corners inside element: '&
          //TRIM(I2S(EdgeHits)),Level=10)

      CALL Info(Caller,'Number of initial corners: '&
          //TRIM(I2S(InitialHits)),Level=10)
      CALL Info(Caller,'Number of active corners: '&
          //TRIM(I2S(ActiveHits)),Level=10)

      CALL Info(Caller,'Number of most subelement corners: '&
          //TRIM(I2S(MaxSubTriangles)),Level=10)
      CALL Info(Caller,'Element of most subelement corners: '&
          //TRIM(I2S(MaxSubElem)),Level=10)

      WRITE( Message,'(A,ES12.5)') 'Total reference area:',TotRefArea
      CALL Info(Caller,Message,Level=8)
      WRITE( Message,'(A,ES12.5)') 'Total integrated area:',TotSumArea
      CALL Info(Caller,Message,Level=8)

      Err = TotSumArea / TotRefArea
      WRITE( Message,'(A,ES15.6)') 'Average ratio in area integration:',Err 
      CALL Info(Caller,Message,Level=8)

      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Maximum relative discrepancy in areas (element: ',MaxErrInd,'):',MaxErr-1.0_dp 
      CALL Info(Caller,Message,Level=8)
      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Minimum relative discrepancy in areas (element: ',MinErrInd,'):',MinErr-1.0_dp 
      CALL Info(Caller,Message,Level=8)

      CALL Info(Caller,'Number of slave entries: '&
          //TRIM(I2S(Nslave)),Level=10)
      CALL Info(Caller,'Number of master entries: '&
          //TRIM(I2S(Nmaster)),Level=10)

      IF( DebugEdge ) THEN
        CALL ListAddConstReal( CurrentModel % Simulation,'res: err',err) 

        WRITE( Message,'(A,ES15.6)') 'Slave entries total sum:', sums
        CALL Info(Caller,Message,Level=8)
        WRITE( Message,'(A,ES15.6)') 'Master entries total sum:', summ
        CALL Info(Caller,Message,Level=8)
        WRITE( Message,'(A,ES15.6)') 'Master entries total sum2:', summ2
        CALL Info(Caller,Message,Level=8)
        WRITE( Message,'(A,ES15.6)') 'Maximum edge projection error:', MaxEdgeErr
        CALL Info(Caller,Message,Level=6)

        CALL ListAddConstReal( CurrentModel % Simulation,'res: sums',sums) 
        CALL ListAddConstReal( CurrentModel % Simulation,'res: summ',summ) 
        CALL ListAddConstReal( CurrentModel % Simulation,'res: summ2',summ2) 
        CALL ListAddConstReal( CurrentModel % Simulation,'res: summabs',summabs) 
        CALL ListAddConstReal( CurrentModel % Simulation,'res: maxedgerr',MaxEdgeErr)
      END IF

    END SUBROUTINE AddProjectorWeakGeneric


    
    ! Return shortest distance squared of a point to a line segment.
    ! This is limited to the spacial case when the point lies in origin. 
    FUNCTION SegmentOriginDistance2(x1,y1,x2,y2) RESULT ( r2 )
      REAL(KIND=dp) :: x1,y1,x2,y2,r2
      REAL(KIND=dp) :: q,xc,yc

      q = ( x1*(x1-x2) + y1*(y1-y2) ) / &
          SQRT((x1**2+y1**2) * ((x1-x2)**2+(y1-y2)**2))
      IF( q <= 0.0_dp ) THEN
        r2 = x1**2 + y1**2
      ELSE IF( q >= 1.0_dp ) THEN
        r2 = x2**2 + y2**2
      ELSE
        xc = x1 + q * (x2-x1)
        yc = y1 + q * (y2-y1)
        r2 = xc**2 + yc**2
      END IF
             
    END FUNCTION SegmentOriginDistance2

    
    !----------------------------------------------------------------------
    ! Create weak projector for the nodes in 1D mesh.
    !----------------------------------------------------------------------
    SUBROUTINE AddProjectorWeak1D()

      INTEGER, TARGET :: IndexesT(3)
      INTEGER, POINTER :: Indexes(:), IndexesM(:)
      INTEGER :: jj,ii,sgn0,k,kmax,ind,indM,nip,nn,inds(10),nM,iM,i2,i2M
      INTEGER :: ElemHits, TotHits, MaxErrInd, MinErrInd, TimeStep, AntiPeriodicHits
      TYPE(Element_t), POINTER :: Element, ElementM
      TYPE(Element_t) :: ElementT 
      TYPE(GaussIntegrationPoints_t) :: IP
      TYPE(Nodes_t) :: Nodes, NodesM, NodesT
      REAL(KIND=dp) :: xt,yt,zt,xmax,xmin,xmaxm,ymaxm,&
          xminm,yminm,DetJ,Wtemp,q,u,v,w,um,vm,wm,val,RefArea,dArea,&
          SumArea,MaxErr,MinErr,Err,uvw(3),val_dual,dx,dxcut, &
          zmin,zmax, zminm, zmaxm
      REAL(KIND=dp) :: TotRefArea, TotSumArea
      REAL(KIND=dp), ALLOCATABLE :: Basis(:), BasisM(:)
      LOGICAL :: LeftCircle, Stat
      TYPE(Mesh_t), POINTER :: Mesh
      TYPE(Variable_t), POINTER :: TimestepVar

      ! These are used temporarily for debugging purposes
      INTEGER :: SaveInd
      LOGICAL :: SaveElem
      CHARACTER(LEN=20) :: FileName

      REAL(KIND=dp), ALLOCATABLE :: CoeffBasis(:), MASS(:,:)

      CALL Info(Caller,'Creating weak constraints using a 1D integrator',Level=8)      

      Mesh => CurrentModel % Solver % Mesh 

      SaveInd = ListGetInteger( BC,'Level Projector Save Element Index',Found )
      TimestepVar => VariableGet( Mesh % Variables,'Timestep',ThisOnly=.TRUE. )
      Timestep = NINT( TimestepVar % Values(1) )
 
      n = Mesh % MaxElementNodes
      ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )
      ALLOCATE( NodesM % x(n), NodesM % y(n), NodesM % z(n) )
      ALLOCATE( NodesT % x(n), NodesT % y(n), NodesT % z(n) )
      ALLOCATE( Basis(n), BasisM(n) )

      IF (BiOrthogonalBasis) ALLOCATE(CoeffBasis(n), MASS(n,n))

      Nodes % y  = 0.0_dp
      NodesM % y = 0.0_dp
      NodesT % y = 0.0_dp
      Nodes % z  = 0.0_dp
      NodesM % z = 0.0_dp
      NodesT % z = 0.0_dp
      yt = 0.0_dp
      zt = 0.0_dp

      MaxErr = 0.0_dp
      MinErr = HUGE( MinErr )
      MaxErrInd = 0
      MinErrInd = 0
      zt = 0.0_dp
      LeftCircle = .FALSE.
     
      ! The temporal element segment used in the numerical integration
      ElementT % TYPE => GetElementType( 202, .FALSE. )
      ElementT % NodeIndexes => IndexesT
      IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2  ) 

      TotHits = 0
      AntiPeriodicHits = 0
      TotRefArea = 0.0_dp
      TotSumArea = 0.0_dp


      DO ind=1,BMesh1 % NumberOfBulkElements

        ! Optionally save the submesh for specified element, for vizualization and debugging
        SaveElem = ( SaveInd == ind )

        Element => BMesh1 % Elements(ind)        
        Indexes => Element % NodeIndexes
        
        n = Element % TYPE % NumberOfNodes        
        Nodes % x(1:n) = BMesh1 % Nodes % x(Indexes(1:n))

        ! There is a discontinuity of angle at 180 degs
        ! If we are working on left-hand-side then add 360 degs to the negative angles
        ! to remove this discontinuity.
        IF( FullCircle ) THEN
          LeftCircle = ( ALL( ABS( Nodes % x(1:n) ) > 90.0_dp ) )
          IF( LeftCircle ) THEN
            DO j=1,n
              IF( Nodes % x(j) < 0.0 ) Nodes % x(j) = &
                  Nodes % x(j) + 360.0_dp
            END DO
          END IF
        END IF

        xmin = MINVAL(Nodes % x(1:n))
        xmax = MAXVAL(Nodes % x(1:n))
        dx = xmax - xmin 

        ! The flattened dimension is always the z-component
        IF( HaveMaxDistance ) THEN
          zmin = MINVAL( BMesh1 % Nodes % z(Indexes(1:n)) )
          zmax = MAXVAL( BMesh1 % Nodes % z(Indexes(1:n)) )
        END IF
                        
        ! Compute the reference area
        u = 0.0_dp; v = 0.0_dp; w = 0.0_dp;
        stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )
        RefArea = detJ * ArcCoeff * SUM( IP % s(1:IP % n) )
        SumArea = 0.0_dp
        
        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_a.dat'
          OPEN( 10,FILE=Filename)
          DO i=1,n
            WRITE( 10, * ) Nodes % x(i)
          END DO
          CLOSE( 10 )
        END IF

        ! Set the values to maintain the size of the matrix
        ! The size of the matrix is used when allocating for utility vectors of contact algo.
        ! This does not set the Projector % InvPerm to nonzero value that is used to 
        ! determine whether there really is a projector. 
        DO i=1,n
          j = InvPerm1(Indexes(i))
          nrow = NodePerm(j)
          IF( nrow == 0 ) CYCLE
          CALL List_AddMatrixIndex(Projector % ListMatrix, nrow, j ) 
        END DO

        ! Currently a n^2 loop but it could be improved
        !--------------------------------------------------------------------
        ElemHits = 0
        DO indM=1,BMesh2 % NumberOfBulkElements
          
          ElementM => BMesh2 % Elements(indM)        
          IndexesM => ElementM % NodeIndexes

          nM = ElementM % TYPE % NumberOfNodes

 
          NodesM % x(1:nM) = BMesh2 % Nodes % x(IndexesM(1:nM))

          ! Treat the left circle differently. 
          IF( LeftCircle ) THEN
            ! Omit the element if it is definitely on the right circle
            IF( ALL( ABS( NodesM % x(1:nM) ) - 90.0_dp < XTol ) ) CYCLE
            DO j=1,nM
              IF( NodesM % x(j) < 0.0_dp ) NodesM % x(j) = &
                  NodesM % x(j) + 360.0_dp
            END DO
          END IF
          
          xminm = MINVAL( NodesM % x(1:nM))
          xmaxm = MAXVAL( NodesM % x(1:nM))

          IF( Repeating ) THEN
            ! Enforce xmaxm to be on the same interval than xmin
            Nrange = FLOOR( (xmaxm-xmin+XTol) / XRange )
            IF( Nrange /= 0 ) THEN
              xminm = xminm - Nrange * XRange
              xmaxm = xmaxm - Nrange * XRange
              NodesM % x(1:nM) = NodesM % x(1:nM) - NRange * XRange 
            END IF

            ! Check whether there could be a intersection in an other interval as well
            IF( xminm + XRange < xmax + XTol ) THEN
              Nrange2 = 1
            ELSE
              Nrange2 = 0
            END IF
          END IF

          IF( FullCircle .AND. .NOT. LeftCircle ) THEN
            IF( xmaxm - xminm > 180.0_dp ) CYCLE
          END IF          

200       IF( xminm >= xmax ) GOTO 100
          IF( xmaxm <= xmin ) GOTO 100

          
          ! This is a cheap test so perform that first, if requested
          IF( HaveMaxDistance ) THEN
            zminm = MINVAL( BMesh2 % Nodes % z(IndexesM(1:nM)) )
            zmaxm = MAXVAL( BMesh2 % Nodes % z(IndexesM(1:nM)) )
            IF( zmaxm < zmin - MaxDistance ) GOTO 100 
            IF( zminm > zmax + MaxDistance ) GOTO 100
          END IF
          

          NodesT % x(1) = MAX( xmin, xminm ) 
          NodesT % x(2) = MIN( xmax, xmaxm ) 
          dxcut = ABS( NodesT % x(1)-NodesT % x(2) )

          ! Too small absolute values may result to problems when inverting matrix
          IF( dxcut < 1.0d-12 ) GOTO 100

          ! Too small relative value is irrelevant
          IF( dxcut < 1.0d-8 * dx ) GOTO 100

          sgn0 = 1
          IF( AntiRepeating ) THEN
            IF ( MODULO(Nrange,2) /= 0 ) THEN
              sgn0 = -1
              AntiPeriodicHits = AntiPeriodicHits + 1
            END IF
          END IF
          
          ElemHits = ElemHits + 1

          IF( SaveElem ) THEN
            FileName = 't'//TRIM(I2S(TimeStep))//'_b'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,nM
              WRITE( 10, * ) NodesM % x(i)
            END DO
            CLOSE( 10 )

            FileName = 't'//TRIM(I2S(TimeStep))//'_e'//TRIM(I2S(ElemHits))//'.dat'
            OPEN( 10,FILE=FileName)
            DO i=1,2
              WRITE( 10, * ) NodesT % x(i)
            END DO
            CLOSE( 10 )           
          END IF
                   
          ! Use somewhat higher integration rules than the default
          IP = GaussPoints( ElementT, ElementT % TYPE % GaussPoints2 ) 
          
          IF(BiOrthogonalBasis) THEN
            MASS  = 0
            CoeffBasis = 0
            DO nip=1, IP % n 
              stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                  IP % v(nip),IP % w(nip),detJ,Basis)

              ! Global coordinate of the integration point
              xt = SUM( Basis(1:2) * NodesT % x(1:2) )
            
              ! Integration weight for current integration point
              Wtemp = DetJ * ArcCoeff * IP % s(nip)
            
              ! Integration point at the slave element
              CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
              stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

              DO i=1,n
                DO j=1,n
                  MASS(i,j) = MASS(i,j) + wTemp * Basis(i) * Basis(j)
                END DO
                CoeffBasis(i) = CoeffBasis(i) + wTemp * Basis(i)
              END DO
            END DO

            CALL InvertMatrix( MASS, n )

            DO i=1,n
              DO j=1,n
                MASS(i,j) = MASS(i,j) * CoeffBasis(i)
              END DO
            END DO
          END IF


          DO nip=1, IP % n 
            stat = ElementInfo( ElementT,NodesT,IP % u(nip),&
                IP % v(nip),IP % w(nip),detJ,Basis)
            
            ! We will actually only use the global coordinates and the integration weight 
            ! from the temporal mesh. 
            
            ! Global coordinate of the integration point
            xt = SUM( Basis(1:2) * NodesT % x(1:2) )
            
            ! Integration weight for current integration point
            ! Use the real arc length so that this projector weights correctly 
            ! in rotational case when used with other projectors.
            Wtemp = ArcCoeff * DetJ * IP % s(nip)
            sumarea = sumarea + Wtemp

            ! Integration point at the slave element
            CALL GlobalToLocal( u, v, w, xt, yt, zt, Element, Nodes )              
            stat = ElementInfo( Element, Nodes, u, v, w, detJ, Basis )

            ! Integration point at the master element
            CALL GlobalToLocal( um, vm, wm, xt, yt, zt, ElementM, NodesM )
            stat = ElementInfo( ElementM, NodesM, um, vm, wm, detJ, BasisM )
            
            IF(BiOrthogonalBasis) THEN
              CoeffBasis = 0._dp
              DO i=1,n
                DO j=1,n
                  CoeffBasis(i) = CoeffBasis(i) + MASS(i,j) * Basis(j)
                END DO
              END DO
            END IF

            ! Add the entries to the projector
            DO j=1,n 
              jj = Indexes(j)                                    
              nrow = NodePerm(InvPerm1(jj))
              IF( nrow == 0 ) CYCLE
              
              Projector % InvPerm(nrow) = InvPerm1(jj)
              val = Basis(j) * Wtemp
              IF(BiorthogonalBasis) THEN
                val_dual = CoeffBasis(j) * Wtemp
              END IF

              DO i=1,n
                CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                      InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val )

                IF(BiorthogonalBasis ) THEN
                  CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                        InvPerm1(Indexes(i)), NodeCoeff * Basis(i) * val_dual )
                END IF
              END DO
              
              DO i=1,nM
                CALL List_AddToMatrixElement(Projector % ListMatrix, nrow, &
                    InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val )

                IF(BiorthogonalBasis) THEN
                  IF(DualMaster .OR. DualLCoeff) THEN
                    CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                      InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val_dual )
                  ELSE
                    CALL List_AddToMatrixElement(Projector % Child % ListMatrix, nrow, &
                      InvPerm2(IndexesM(i)), -sgn0 * NodeScale * NodeCoeff * BasisM(i) * val )
                  END IF
                END IF
              END DO
            END DO

            ! Add the entries to the dual projector 
            IF( CreateDual ) THEN
              DO j=1,nM 
                jj = IndexesM(j)                                    
                nrow = DualNodePerm(InvPerm2(jj))
                IF( nrow == 0 ) CYCLE
                
                DualProjector % InvPerm(nrow) = InvPerm2(jj)
                val = BasisM(j) * Wtemp

                DO i=1,nM
                  CALL List_AddToMatrixElement(DualProjector % ListMatrix, nrow, &
                      InvPerm2(IndexesM(i)), sgn0 * NodeCoeff * BasisM(i) * val ) 
                END DO

                DO i=1,n
                  !IF( ABS( val * BasisM(i) ) < 1.0d-10 ) CYCLE
                  CALL List_AddToMatrixElement(DualProjector % ListMatrix, nrow, &
                      InvPerm1(Indexes(i)), -NodeScale * NodeCoeff * Basis(i) * val )                   
                END DO
              END DO
            END IF
          END DO

100       IF( Repeating ) THEN
            IF( NRange2 /= 0 ) THEN
              xminm = xminm + Nrange2 * XRange
              xmaxm = xmaxm + Nrange2 * XRange
              NodesM % x(1:n) = NodesM % x(1:n) + NRange2 * XRange 
              NRange = NRange + NRange2
              NRange2 = 0
              GOTO 200
            END IF
          END IF

        END DO

        IF( SaveElem ) THEN
          FileName = 't'//TRIM(I2S(TimeStep))//'_n.dat'
          OPEN( 10,FILE=Filename)
          WRITE( 10, * ) ElemHits 
          CLOSE( 10 )
        END IF
        
        TotHits = TotHits + ElemHits
        TotSumArea = TotSumArea + SumArea
        TotRefArea = TotRefArea + RefArea

        Err = SumArea / RefArea
        IF( Err > MaxErr ) THEN
          MaxErr = Err
          MaxErrInd = Err
        END IF
        IF( Err < MinErr ) THEN
          MinErr = Err
          MinErrInd = ind
        END IF
      END DO

      DEALLOCATE( Nodes % x, Nodes % y, Nodes % z )
      DEALLOCATE( NodesM % x, NodesM % y, NodesM % z )
      DEALLOCATE( NodesT % x, NodesT % y, NodesT % z )
      DEALLOCATE( Basis, BasisM )

      CALL Info(Caller,'Number of integration pairs: '&
          //TRIM(I2S(TotHits)),Level=10)
      IF( AntiPeriodicHits > 0 ) THEN
        CALL Info(Caller,'Number of antiperiodic pairs: '&
          //TRIM(I2S(AntiPeriodicHits)),Level=10)
      END IF

      WRITE( Message,'(A,ES12.5)') 'Total reference length:',TotRefArea / ArcCoeff
      CALL Info(Caller,Message,Level=8) 
      WRITE( Message,'(A,ES12.5)') 'Total integrated length:',TotSumArea / ArcCoeff
      CALL Info(Caller,Message,Level=8)

      Err = TotSumArea / TotRefArea
      WRITE( Message,'(A,ES12.3)') 'Average ratio in length integration:',Err 
      CALL Info(Caller,Message,Level=8)

      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Maximum relative discrepancy in length (element: ',MaxErrInd,'):',MaxErr-1.0_dp 
      CALL Info(Caller,Message,Level=8)
      WRITE( Message,'(A,I0,A,ES12.4)') &
          'Minimum relative discrepancy in length (element: ',MinErrInd,'):',MinErr-1.0_dp 
      CALL Info(Caller,Message,Level=8)


    END SUBROUTINE AddProjectorWeak1D

  END FUNCTION LevelProjector
  !------------------------------------------------------------------------------


!---------------------------------------------------------------------------
!> Create a Galerkin projector related to discontinuous interface.
!> This uses the information stored when the discontinuous interface
!> was first coined. This enables simple one-to-one mapping. Integration
!> weight is used for the nodel projector to allow physical jump conditions.
!> For the edge dofs there is no such jumps and hence the projector uses
!> weights of one. 
!---------------------------------------------------------------------------
  FUNCTION WeightedProjectorDiscont(Mesh, bc ) RESULT ( Projector )
    !---------------------------------------------------------------------------
    USE Lists
    USE ListMatrix

    TYPE(Mesh_t), POINTER :: Mesh
    INTEGER :: bc
    TYPE(Matrix_t), POINTER :: Projector
    !--------------------------------------------------------------------------
    INTEGER, POINTER :: NodePerm(:)
    TYPE(Model_t), POINTER :: Model
    TYPE(GaussIntegrationPoints_t), TARGET :: IntegStuff
    INTEGER :: p,q,i,j,it,nn,n,m,t,NoOrigNodes, NoDiscontNodes, indp, indq, &
        e1, e2, e12, i1, i2, j1, j2, ParentMissing, ParentFound, PosSides, ActSides, &
        InvPermSize, indpoffset
    INTEGER, POINTER :: Rows(:),Cols(:), InvPerm(:)
    REAL(KIND=dp), POINTER :: Values(:), Basis(:), WBasis(:,:), &
                 Wbasis2(:,:),RotWBasis(:,:),dBasisdx(:,:)
    REAL(KIND=dp) :: u,v,w,val,detJ,Scale,x,weight,Coeff
    INTEGER, ALLOCATABLE :: Indexes(:), DiscontIndexes(:)
    TYPE(Nodes_t) :: ElementNodes
    TYPE(Element_t), POINTER :: Element, Left, Right, OldFace, NewFace, Swap
    LOGICAL :: Stat,DisCont,Found,NodalJump,AxisSym, SetDiag, &
        SetDiagEdges, DoNodes, DoEdges, LocalConstraints, NoHalo
    LOGICAL, ALLOCATABLE :: EdgeDone(:)
    REAL(KIND=dp) :: point(3), uvw(3), DiagEps
    INTEGER, ALLOCATABLE :: EQind(:)
    INTEGER, POINTER :: OldMap(:,:), NewMap(:,:)
    TYPE(ValueList_t), POINTER :: BCParams
    LOGICAL :: CheckHaloNodes
    LOGICAL, POINTER :: HaloNode(:)
    CHARACTER(*), PARAMETER :: Caller = 'WeightedProjectorDiscont'
    
    CALL Info(Caller,'Creating projector for discontinuous boundary '&
         //TRIM(I2S(bc)),Level=7)

    Projector => NULL()
    IF( .NOT. Mesh % DisContMesh ) THEN
      CALL Warn(Caller,'Discontinuous mesh not created?')
      RETURN
    END IF

    Model => CurrentModel

    j = 0
    DO i=1,Model % NumberOfBCs
      IF( ListGetLogical(Model % BCs(i) % Values,'Discontinuous Boundary',Found) ) THEN
        j = j + 1
      END IF
    END DO
    IF( j > 1 ) THEN
      CALL Warn(Caller,'One BC (not '&
          //TRIM(I2S(j))//') only for discontinuous boundary!')
    END IF
 
    BCParams => Model % BCs(bc) % Values

    Scale = ListGetCReal( BCParams,'Mortar BC Scaling',Stat )  
    IF(.NOT. Stat) Scale = -1.0_dp

    NodalJump = ListCheckPrefix( BCParams,'Mortar BC Coefficient')
    IF(.NOT. NodalJump ) THEN
      NodalJump = ListCheckPrefix( BCParams,'Mortar BC Resistivity')
    END IF

    ! Take the full weight when creating the constraints since the values will 
    ! not be communicated
    LocalConstraints = ListGetLogical(Model % Solver % Values, &
        'Partition Local Projector',Found)
    IF(.NOT. Found ) LocalConstraints = ListGetLogical(Model % Solver % Values, &
        'Partition Local Constraints',Found)

    ! Don't consider halo when creating discontinuity
    NoHalo = ListGetLogical(Model % Solver % Values, &
        'Projector No Halo',Found)

    ! Don't consider single halo nodes when creating discontinuity
    CheckHaloNodes = ListGetLogical( Model % Solver % Values,&
        'Projector No Halo Nodes',Found ) 
    IF( CheckHaloNodes ) THEN
      CALL MarkHaloNodes( Mesh, HaloNode, CheckHaloNodes )
    END IF


    IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',Found ) ) THEN
      DoEdges = .FALSE. 
    ELSE IF( ListGetLogical( BCParams,'Projector Skip Edges',Found ) ) THEN
      DoEdges = .FALSE.
    ELSE
      DoEdges = ( Mesh % NumberOfEdges > 0 )
    END IF
    IF( DoEdges .AND. Mesh % NumberOfEdges == 0 ) THEN
      CALL Warn(Caller,'Edge basis requested but mesh has no edges!')
      DoEdges = .FALSE.
    END IF

    IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',Found ) ) THEN
      DoNodes = .FALSE. 
    ELSE IF( ListGetLogical( BCParams,'Projector Skip Nodes',Found ) ) THEN
      DoNodes = .FALSE.
    ELSE
      DoNodes = ( Mesh % NumberOfNodes > 0 )
    END IF

    ! Should the projector be diagonal or mass matrix type 
    SetDiag = ListGetLogical( BCParams,'Mortar BC Diag',Found ) 

    IF(.NOT. Found ) SetDiag = ListGetLogical( BCParams, 'Use Biorthogonal Basis', Found)

    ! If we want to eliminate the constraints we have to have a biortgonal basis
    IF(.NOT. Found ) THEN
      SetDiag = ListGetLogical( CurrentModel % Solver % Values, &
          'Eliminate Linear Constraints',Found )
      IF( SetDiag ) THEN
        CALL Info(Caller,&
            'Setting > Use Biorthogonal Basis < to True to enable elimination',Level=8)
      END IF
    END IF


    SetDiagEdges = ListGetLogical( BCParams,'Mortar BC Diag Edges',Found )
    IF(.NOT. Found ) SetDiagEdges = SetDiag
    DiagEps = ListGetConstReal( BCParams,'Mortar BC Diag Eps',Found ) 

    ! Integration weights should follow the metrics if we want physical nodal jumps. 
    AxisSym = .FALSE.
    IF ( CurrentCoordinateSystem() == AxisSymmetric .OR. &
        CurrentCoordinateSystem() == CylindricSymmetric ) THEN
      IF( NodalJump ) THEN
        AxisSym = .TRUE.
      ELSE IF (ASSOCIATED(CurrentModel % Solver)) THEN
        AxisSym = ListGetLogical(CurrentModel % Solver % Values,'Projector Metrics',Found)
      END IF
      IF( AxisSym ) CALL Info('weightedProjectorDiscont','Projector will be weighted for axi symmetry',Level=7)
    END IF


    n = Mesh % MaxElementDOFs
    ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n) )
    ALLOCATE( Indexes(n), DisContIndexes(n), Basis(n), Wbasis(n,3), &
            Wbasis2(n,3), dBasisdx(n,3), RotWBasis(n,3) )
    Indexes = 0
    Basis = 0.0_dp
    DiscontIndexes = 0

    NodePerm => Mesh % DisContPerm
    NoOrigNodes = SIZE( NodePerm ) 
    NoDiscontNodes = COUNT( NodePerm > 0 ) 

    IF( DoNodes ) THEN
      indpoffset = NoDiscontNodes
    ELSE
      indpoffset = 0
    END IF
    InvPerm => NULL()
    InvPermSize = indpoffset
    
    ! Compute the number of potential edges. This mimics the loop that really creates the projector 
    ! below. 
    IF( DoEdges ) THEN
      ALLOCATE( EdgeDone( Mesh % NumberOfEdges ) )
      EdgeDone = .FALSE.
      indp = indpoffset

      DO t = 1, Mesh % NumberOfBoundaryElements
        
        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )        
        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
        
        Left => Element % BoundaryInfo % Left
        Right => Element % BoundaryInfo % Right 
        
        IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
          CYCLE
        END IF

        ActSides = 0
        IF( ASSOCIATED( Left ) ) THEN
          IF( Left % PartIndex == ParEnv % myPE ) ActSides = ActSides + 1
        END IF
        IF( ASSOCIATED( Right ) ) THEN
          IF( Right % PartIndex == ParEnv % myPe ) ActSides = ActSides + 1
        END IF 
        IF( NoHalo .AND. ActSides == 0 ) CYCLE
        
        ! Consistently choose the face with the old edges 
        IF( ALL( Left % NodeIndexes <= NoOrigNodes ) ) THEN
          OldFace => Left
        ELSE IF( ALL( Right % NodeIndexes <= NoOrigNodes ) ) THEN
          OldFace => Right
        ELSE
          CALL Warn(Caller,'Neither face is purely old!')
          CYCLE
        END IF

        OldMap => GetEdgeMap( OldFace % TYPE % ElementCode / 100)

        DO i = 1,OldFace % TYPE % NumberOfEdges          
          e1 = OldFace % EdgeIndexes(i)
          IF( EdgeDone(e1) ) CYCLE

          i1 = OldFace % NodeIndexes( OldMap(i,1) )
          i2 = OldFace % NodeIndexes( OldMap(i,2) )
                    
          ! i1 and i2 were already checked to be "old" nodes
          IF( NodePerm(i1) == 0 ) CYCLE
          IF( NodePerm(i2) == 0 ) CYCLE

          indp = indp + 1
          EdgeDone(e1) = .TRUE.
        END DO
      END DO
      InvPermSize = indp
      CALL Info(Caller,&
          'Size of InvPerm estimated to be: '//TRIM(I2S(InvPermSize)),Level=8)
    END IF

    ! Ok, nothing to do just go end tidy things up
    IF( InvPermSize == 0 ) GOTO 100

    ! Create a list matrix that allows for unspecified entries in the matrix 
    ! structure to be introduced.
    Projector => AllocateMatrix()
    Projector % FORMAT = MATRIX_LIST
    Projector % ProjectorType = PROJECTOR_TYPE_GALERKIN
    Projector % ProjectorBC = bc
    
    ! Create the inverse permutation needed when the projector matrix is added to the global 
    ! matrix. 
    ALLOCATE( Projector % InvPerm( InvPermSize ) )
    InvPerm => Projector % InvPerm
    InvPerm = 0

    
    ! Projector for the nodal dofs. 
    !------------------------------------------------------------------------
    IF( DoNodes ) THEN

      ParentMissing = 0
      ParentFound = 0
      DO t = 1, Mesh % NumberOfBoundaryElements

        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )
        n = Element % TYPE % NumberOfNodes        
        Indexes(1:n) = Element % NodeIndexes(1:n)

        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
        
        Left => Element % BoundaryInfo % Left
        Right => Element % BoundaryInfo % Right 

        ! Here we really need both sides to be able to continue!
        !IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
        !  ParentMissing = ParentMissing + 1
        !  CYCLE
        !END IF

        PosSides = 0
        ActSides = 0
        IF( ASSOCIATED( Left ) ) THEN
          PosSides = PosSides + 1
          IF( Left % PartIndex == ParEnv % myPE ) ActSides = ActSides + 1
        END IF
        IF( ASSOCIATED( Right ) ) THEN
          PosSides = PosSides + 1
          IF( Right % PartIndex == ParEnv % myPe ) ActSides = ActSides + 1
        END IF
        IF( NoHalo .AND. ActSides == 0 ) CYCLE        

        IF( LocalConstraints ) THEN
          Coeff = 1.0_dp
        ELSE
          Coeff = 1.0_dp * ActSides / PosSides 
        END IF
        IF( ABS( Coeff ) < TINY( 1.0_dp ) ) CYCLE

        ParentFound = ParentFound + 1

        ElementNodes % x(1:n) = Mesh % Nodes % x(Indexes(1:n))
        ElementNodes % y(1:n) = Mesh % Nodes % y(Indexes(1:n))
        ElementNodes % z(1:n) = Mesh % Nodes % z(Indexes(1:n))

        IF( ALL( NodePerm(Indexes(1:n)) == 0 ) ) CYCLE
        
        IF( CheckHaloNodes ) THEN
          IF( ALL( HaloNode(Indexes(1:n)) ) ) CYCLE
        END IF

        ! Get the indexes on the other side of the discontinuous boundary
        DO i=1,n
          j = NodePerm( Indexes(i) ) 
          IF( j == 0 ) THEN
            DiscontIndexes(i) = Indexes(i)
          ELSE
            DiscontIndexes(i) = j + NoOrigNodes
          END IF
        END DO

        IntegStuff = GaussPoints( Element )
        DO j=1,IntegStuff % n
          u = IntegStuff % u(j)
          v = IntegStuff % v(j)
          w = IntegStuff % w(j)

          Stat = ElementInfo(Element, ElementNodes, u, v, w, detJ, Basis)

          weight = Coeff * detJ * IntegStuff % s(j)
          IF( AxisSym ) THEN
            x = SUM( Basis(1:n) * ElementNodes % x(1:n) )
            weight = weight * x
          END IF

          DO p=1,n             
            indp = NodePerm( Indexes(p) )
            IF( indp == 0 ) CYCLE
            IF( CheckHaloNodes ) THEN
              IF( HaloNode( Indexes(p) ) ) CYCLE
            END IF

            val = weight * Basis(p)

            ! Only set for the nodes are are really used
            InvPerm(indp) = Indexes(p)

            IF( SetDiag ) THEN
              CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                  Indexes(p), val ) 

              CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                  DiscontIndexes(p), Scale * val )             
            ELSE
              DO q=1,n

                indq = NodePerm(Indexes(q))
                IF( indq == 0 ) CYCLE

                IF( CheckHaloNodes ) THEN
                  IF( HaloNode( Indexes(p) ) ) CYCLE
                END IF
                
                CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                    Indexes(q), Basis(q) * val ) 
                CALL List_AddToMatrixElement(Projector % ListMatrix, indp, &
                    DiscontIndexes(q), Scale * Basis(q) * val ) 
              END DO
            END IF
          END DO
        END DO
      END DO
      IF( ParentMissing > 0 ) THEN
        CALL Warn(Caller,'Number of half-sided discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentMissing)) )
        CALL Warn(Caller,'Number of proper discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentFound)) )
      END IF
      CALL Info(Caller,'Created projector for '&
          //TRIM(I2S(NoDiscontNodes))//' discontinuous nodes',Level=10)
    END IF


    ! Create the projector also for edge dofs if they exist and are
    ! requested. 
    !----------------------------------------------------------------
    IF( DoEdges ) THEN
      ParentMissing = 0
      ParentFound = 0
      n = Mesh % NumberOfNodes

      val = 1.0_dp
      Scale = 1.0_dp

      indp = indpoffset
      ALLOCATE( Eqind(Mesh % NumberOfEdges) ); EQind = 0

      DO t = 1, Mesh % NumberOfBoundaryElements
        
        Element => Mesh % Elements(Mesh % NumberOfBulkElements + t )
        
        IF ( Element % BoundaryInfo % Constraint /= Model % BCs(bc) % Tag ) CYCLE
        
        Left => Element % BoundaryInfo % Left
        Right => Element % BoundaryInfo % Right 
        
        ! Here we really need both sides to be able to continue!
        IF(.NOT. ASSOCIATED( Left ) .OR. .NOT. ASSOCIATED( Right ) ) THEN
          ParentMissing = ParentMissing + 1
          CYCLE
        END IF

        PosSides = 0
        ActSides = 0
        IF( ASSOCIATED( Left ) ) THEN
          PosSides = PosSides + 1
          IF( Left % PartIndex == ParEnv % myPE ) ActSides = ActSides + 1
        END IF
        IF( ASSOCIATED( Right ) ) THEN
          PosSides = PosSides + 1
          IF( Right % PartIndex == ParEnv % myPe ) ActSides = ActSides + 1
        END IF

        IF( NoHalo .AND. ActSides == 0 ) CYCLE

        IF( LocalConstraints ) THEN
          Coeff = 1.0_dp
        ELSE          
          Coeff = (1.0_dp * ActSides) / (1.0_dp * PosSides)
        END IF

        ! Consistently choose the face with the old edges
        IF( ALL( Left % NodeIndexes <= NoOrigNodes ) ) THEN
        ELSE IF( ALL( Right % NodeIndexes <= NoOrigNodes ) ) THEN
          swap  => Left
          Left  => Right
          Right => swap
        ELSE
          ! We already complained once
          CYCLE
        END IF

        OldFace => Find_Face( Mesh, Left, Element )
        nn = SIZE(Element % NodeIndexes)
        Indexes(1:nn) = Element % NodeIndexes
        Element % NodeIndexes = NodePerm(Indexes(1:nn)) + NoOrigNodes
        NewFace => Find_Face( Mesh, Right, Element )
        Element % NodeIndexes = Indexes(1:nn)
 
        ParentFound = ParentFound + 1

        OldMap => GetEdgeMap( OldFace % TYPE % ElementCode / 100 )
        NewMap => GetEdgeMap( NewFace % TYPE % ElementCode / 100 )

        IntegStuff = GaussPoints( oldface )
        DO it = 1,IntegStuff % n
          u = integstuff % u(it)
          v = integstuff % v(it)
          w = integstuff % w(it)

          nn = OldFace % TYPE % NumberOfNodes
          ElementNodes % x(1:nn) = Mesh % Nodes % x(oldface % NodeIndexes(1:nn))
          ElementNodes % y(1:nn) = Mesh % Nodes % y(oldface % NodeIndexes(1:nn))
          ElementNodes % z(1:nn) = Mesh % Nodes % z(oldface % NodeIndexes(1:nn))

          Stat = ElementInfo( OldFace, ElementNodes,u,v,w, DetJ, Basis,dBasisdx )
          CALL GetEdgeBasis( OldFace, Wbasis, RotWbasis, Basis, dBasisdx )

          Point(1) = SUM(Basis(1:nn) * ElementNodes % x(1:nn))
          Point(2) = SUM(Basis(1:nn) * ElementNodes % y(1:nn))
          Point(3) = SUM(Basis(1:nn) * ElementNodes % z(1:nn))

          nn = NewFace % TYPE % NumberOfNodes
          ElementNodes % x(1:nn) = Mesh % Nodes % x(newface % NodeIndexes(1:nn))
          ElementNodes % y(1:nn) = Mesh % Nodes % y(newface % NodeIndexes(1:nn))
          ElementNodes % z(1:nn) = Mesh % Nodes % z(newface % NodeIndexes(1:nn))

          Found = PointInElement( NewFace, ElementNodes, Point, uvw )
          u = uvw(1); v=uvw(2); w=uvw(3)
          Stat = ElementInfo(NewFace, ElementNodes,u,v,w, detj, Basis,dbasisdx )
          CALL GetEdgeBasis( NewFace, Wbasis2, RotwBasis, Basis, dBasisdx )

          Weight = detJ * IntegStuff % s(it) * Coeff
        
          ! Go through combinations of edges and find the edges for which the 
          ! indexes are the same. 
          DO i = 1,OldFace % TYPE % NumberOfEdges
            e1 = OldFace % EdgeIndexes(i)

            IF ( EQind(e1) == 0 ) THEN
              indp = indp + 1
              EQind(e1) = indp
              InvPerm(indp) = n + e1
            END IF

            IF( SetDiagEdges ) THEN
              i1 = OldFace % NodeIndexes( OldMap(i,1) )
              i1 = NoOrigNodes + NodePerm(i1)
              i2 = OldFace % NodeIndexes( OldMap(i,2) )
              i2 = NoOrigNodes + NodePerm(i2)

              DO j = 1,NewFace % TYPE % NumberOfEdges
                j1 = NewFace % NodeIndexes( NewMap(j,1) )
                j2 = NewFace % NodeIndexes( NewMap(j,2) )
                IF (i1==j1 .AND. i2==j2 .OR. i1==j2 .AND. i2==j1 ) EXIT
              END DO
              val = Weight * SUM(WBasis(i,:) * Wbasis(i,:))
              IF ( ABS(Val)>= 10*AEPS ) &
                  CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e1, Val )
              
              e2  = NewFace % EdgeIndexes(j)
              val = Weight * SUM(WBasis(i,:) * Wbasis2(j,:))
              IF ( ABS(val) >= 10*AEPS ) &
                  CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e2, -Val )              
            ELSE
              DO j = 1,NewFace % TYPE % NumberOfEdges
                e2  = NewFace % EdgeIndexes(j)
                e12 = OldFace % EdgeIndexes(j)
                
                val = Weight * SUM(WBasis(i,:) * Wbasis(j,:))
                IF ( ABS(Val)>= 10*AEPS ) &
                    CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e12, Val )
                
                val = Weight * SUM(WBasis(i,:) * Wbasis2(j,:))
                IF ( ABS(val) >= 10*AEPS ) &
                    CALL List_AddToMatrixElement(Projector % ListMatrix, EQind(e1), n + e2, -Val )
              END DO
            END IF

          END DO
        END DO
      END DO

      DEALLOCATE( EdgeDone )
      IF( .NOT. DoNodes .AND. ParentMissing > 0 ) THEN
        CALL Warn(Caller,'Number of half-sided discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentMissing)) )
        CALL Warn(Caller,'Number of proper discontinuous BC elements in partition '&
           //TRIM(I2S(ParEnv % myPE))//': '//TRIM(I2S(ParentFound)) )
      END IF
      CALL Info(Caller,'Created projector for '&
          //TRIM(I2S(indp-NoDiscontNodes))//' discontinuous edges',Level=10)
    END IF

    ! Convert from list matrix to CRS matrix format
    CALL List_ToCRSMatrix(Projector)

    IF( Projector % NumberOfRows > 0) THEN
      CALL CRS_SortMatrix(Projector,.TRUE.)
      CALL Info(Caller,'Number of entries in projector matrix: '//&
          TRIM(I2S(SIZE(Projector % Cols)) ), Level=9)
    ELSE
      CALL FreeMatrix(Projector); Projector=>NULL()
    END IF

100 DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z )
    DEALLOCATE( Indexes, DisContIndexes, Basis, dBasisdx, WBasis, WBasis2, RotWBasis )
    IF( CheckHaloNodes ) DEALLOCATE( HaloNode )

           
  END FUNCTION WeightedProjectorDiscont
  !------------------------------------------------------------------------------

   !---------------------------------------------------------------------------
  ! Simply fitting of cylinder into a point cloud. This is done in two phases.
  ! 1) The axis of the cylinder is found by minimizing the \sum((n_i*t)^2)
  !    for each component of of t where n_i:s are the surface normals. 
  !    This is fully generic and assumes no positions. 
  ! 2) The radius and center point of the cylinder are found by fitting a circle
  !    in the chosen plane to three representative points. Currently the fitting
  !    can only be done in x-y plane. 
  !---------------------------------------------------------------------------
  SUBROUTINE CylinderFit(PMesh, PParams) 
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Valuelist_t), POINTER :: PParams

    INTEGER :: i,j,k,n,t,AxisI,iter
    INTEGER, POINTER :: NodeIndexes(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(Nodes_t) :: Nodes
    REAL(KIND=dp) :: NiNj(3,3),A(3,3),F(3),M11,M12,M13,M14
    REAL(KIND=dp) :: d1,d2,MinDist,MaxDist,Dist,X0,Y0,Rad
    REAL(KIND=dp) :: Normal(3), AxisNormal(3), Tangent1(3), Tangent2(3), Coord(3), &
        CircleCoord(3,3)
    INTEGER :: CircleInd(3) 
    CHARACTER(*), PARAMETER :: Caller = 'CylinderFit'
    
    CALL Info(Caller,'Trying to fit a cylinder to the surface patch',Level=10)

    NiNj = 0.0_dp

    n = PMesh % MaxElementNodes
    ALLOCATE( Nodes % x(n), Nodes % y(n), Nodes % z(n) )

    ! If the initial mesh is in 2D there is really no need to figure out the 
    ! direction of the rotational axis. It can only be aligned with the z-axis. 
    IF( CurrentModel % Mesh % MeshDim == 2 ) THEN
      AxisNormal = 0.0_dp
      AxisNormal(3) = 1.0_dp
      GOTO 100 
    END IF


    ! Compute the inner product of <N*N> for the elements
    DO t=1, PMesh % NumberOfBulkElements
      Element => PMesh % Elements(t)
      
      n = Element % TYPE % NumberOfNodes
      NodeIndexes => Element % NodeIndexes
      
      Nodes % x(1:n) = PMesh % Nodes % x(NodeIndexes(1:n))
      Nodes % y(1:n) = PMesh % Nodes % y(NodeIndexes(1:n))
      Nodes % z(1:n) = PMesh % Nodes % z(NodeIndexes(1:n))           
      
      Normal = NormalVector( Element, Nodes, Check = .FALSE. ) 

      DO i=1,3
        DO j=1,3
          NiNj(i,j) = NiNj(i,j) + Normal(i) * Normal(j)
        END DO
      END DO      
    END DO

    ! Normalize by the number of boundary elements
    NiNj = NiNj / PMesh % NumberOfBulkElements

    ! The potential direction for the cylinder axis is the direction with 
    ! least hits for the normal.
    AxisI = 1 
    DO i=2,3
      IF( NiNj(i,i) < NiNj(AxisI,AxisI) ) AxisI = i 
    END DO

    CALL Info(Caller,'Axis coordinate set to be: '//TRIM(I2S(AxisI)))

    ! Keep the dominating direction fixed and iteratively solve the two other directions
    AxisNormal = 0.0_dp
    AxisNormal(AxisI) = 1.0_dp

    ! Basically we could solve from equation Ax=0 the tangent but only up to a constant.
    ! Thus we enforce the axis direction to one by manipulation the matrix equation 
    ! thereby can get a unique solution. 
    A = NiNj
    A(AxisI,1:3) = 0.0_dp
    A(AxisI,AxisI) = 1.0_dp
    CALL InvertMatrix( A, 3 )
    AxisNormal = A(1:3,AxisI)

    ! Normalize the axis normal length to one    
    AxisNormal = AxisNormal / SQRT( SUM( AxisNormal ** 2 ) )
    IF( 1.0_dp - ABS( AxisNormal(3) ) > 1.0d-5 ) THEN
      CALL Warn(Caller,'The cylinder axis is not aligned with z-axis!')
    END IF

100 CALL TangentDirections( AxisNormal,Tangent1,Tangent2 )

    IF(.FALSE.) THEN
      PRINT *,'Axis Normal:',AxisNormal
      PRINT *,'Axis Tangent 1:',Tangent1
      PRINT *,'Axis Tangent 2:',Tangent2
    END IF

    ! Finding three points with maximum distance in the tangent directions

    ! First, find the single extremum point in the first tangent direction
    ! Save the local coordinates in the N-T system of the cylinder
    MinDist = HUGE(MinDist) 
    DO i=1, PMesh % NumberOfNodes
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)

      d1 = SUM( Tangent1 * Coord )
      IF( d1 < MinDist ) THEN
        MinDist = d1
        CircleInd(1) = i
      END IF
    END DO

    i = CircleInd(1)
    Coord(1) = PMesh % Nodes % x(i)
    Coord(2) = PMesh % Nodes % y(i)
    Coord(3) = PMesh % Nodes % z(i)
      
    CircleCoord(1,1) = SUM( Tangent1 * Coord ) 
    CircleCoord(1,2) = SUM( Tangent2 * Coord ) 
    CircleCoord(1,3) = SUM( AxisNormal * Coord )
   

    !PRINT *,'MinDist1:',MinDist,CircleInd(1),CircleCoord(1,:)

    ! Find two more points such that their minimum distance to the previous point(s)
    ! is maximized. This takes some time but the further the nodes are apart the more 
    ! accurate it will be to fit the circle to the points. Also if there is just 
    ! a symmetric section of the cylinder it is important to find the points rigorously.
    DO j=2,3
      ! The maximum minimum distance of any node from the previously defined nodes
      MaxDist = 0.0_dp
      DO i=1, PMesh % NumberOfNodes
        Coord(1) = PMesh % Nodes % x(i)
        Coord(2) = PMesh % Nodes % y(i)
        Coord(3) = PMesh % Nodes % z(i)
        
        ! Minimum distance from the previously defined nodes
        MinDist = HUGE(MinDist)
        DO k=1,j-1
          d1 = SUM( Tangent1 * Coord )
          d2 = SUM( Tangent2 * Coord )
          Dist = ( d1 - CircleCoord(k,1) )**2 + ( d2 - CircleCoord(k,2) )**2
          MinDist = MIN( Dist, MinDist )
        END DO
        
        ! If the minimum distance is greater than in any other node, choose this
        IF( MaxDist < MinDist ) THEN
          MaxDist = MinDist 
          CircleInd(j) = i
        END IF
      END DO

      ! Ok, we have found the point now set the circle coordinates 
      i = CircleInd(j)
      Coord(1) = PMesh % Nodes % x(i)
      Coord(2) = PMesh % Nodes % y(i)
      Coord(3) = PMesh % Nodes % z(i)
      
      CircleCoord(j,1) = SUM( Tangent1 * Coord ) 
      CircleCoord(j,2) = SUM( Tangent2 * Coord ) 
      CircleCoord(j,3) = SUM( AxisNormal * Coord )
    END DO
      

    !PRINT *,'Circle Indexes:',CircleInd

    ! Given three nodes it is possible to analytically compute the center point and
    ! radius of the cylinder from a 4x4 determinant equation. The matrices values
    ! m1i are the determinants of the comatrices. 

    A(1:3,1) = CircleCoord(1:3,1)  ! x
    A(1:3,2) = CircleCoord(1:3,2)  ! y
    A(1:3,3) = 1.0_dp
    m11 = Det3x3( a )

    A(1:3,1) = CircleCoord(1:3,1)**2 + CircleCoord(1:3,2)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(1:3,2)  ! y
    A(1:3,3) = 1.0_dp
    m12 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1:3,1)**2 + CircleCoord(1:3,2)**2  ! x^2+y^2
    A(1:3,2) = CircleCoord(1:3,1)  ! x
    A(1:3,3) = 1.0_dp
    m13 = Det3x3( a )
 
    A(1:3,1) = CircleCoord(1:3,1)**2 + CircleCoord(1:3,2)**2 ! x^2+y^2
    A(1:3,2) = CircleCoord(1:3,1)  ! x
    A(1:3,3) = CircleCoord(1:3,2)  ! y
    m14 = Det3x3( a )

    !PRINT *,'determinants:',m11,m12,m13,m14

    IF( ABS( m11 ) < EPSILON( m11 ) ) THEN
      CALL Fatal(Caller,'Points cannot be an a circle')
    END IF

    X0 =  0.5_dp * m12 / m11 
    Y0 = -0.5_dp * m13 / m11
    rad = SQRT( x0**2 + y0**2 + m14/m11 )

    Coord = x0 * Tangent1 + y0 * Tangent2

    !PRINT *,'Center point in cartesian coordinates:',Coord
    
    CALL ListAddConstReal( PParams,'Rotational Projector Center X',Coord(1))
    CALL ListAddConstReal( PParams,'Rotational Projector Center Y',Coord(2))
    CALL ListAddConstReal( PParams,'Rotational Projector Center Z',Coord(3))

    CALL ListAddConstReal( PParams,'Rotational Projector Normal X',AxisNormal(1))
    CALL ListAddConstReal( PParams,'Rotational Projector Normal Y',AxisNormal(2))
    CALL ListAddConstReal( PParams,'Rotational Projector Normal Z',AxisNormal(3))

    
  CONTAINS
    
    ! Compute the value of 3x3 determinant
    !-------------------------------------------
    FUNCTION Det3x3( A ) RESULT ( val ) 
      
      REAL(KIND=dp) :: A(:,:)
      REAL(KIND=dp) :: val

      val = A(1,1) * ( A(2,2) * A(3,3) - A(2,3) * A(3,2) ) &
          - A(1,2) * ( A(2,1) * A(3,3) - A(2,3) * A(3,1) ) &
          + A(1,3) * ( A(2,1) * A(3,2) - A(2,2) * A(3,1) ) 

    END FUNCTION Det3x3

  END SUBROUTINE CylinderFit

  !---------------------------------------------------------------------------
  !> Given two interface meshes for nonconforming rotating boundaries make 
  !> a coordinate transformation to (phi,z) level where the interpolation
  !> accuracy is not limited by the curvilinear coordinates. Also ensure
  !> that the master nodes manipulated so they for sure hit the target nodes.
  !---------------------------------------------------------------------------
  SUBROUTINE RotationalInterfaceMeshes(BMesh1, BMesh2, BParams, Cylindrical, &
      Radius, FullCircle )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    REAL(KIND=dp) :: Radius
    LOGICAL :: FullCircle, Cylindrical
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: x1_min(3),x1_max(3),x2_min(3),x2_max(3),&
        x1r_min(3),x1r_max(3),x2r_min(3),x2r_max(3)
    REAL(KIND=dp) :: x(3), xcyl(3),rad2deg,F1min,F1max,F2min,F2max,dFii1,dFii2,eps_rad,&
        err1,err2,dF,Fii,Fii0,Nsymmetry,fmin,fmax,DegOffset,rad,alpha,x0(3),xtmp(3),&
        Normal(3), Tangent1(3), Tangent2(3) 
    REAL(KIND=dp), POINTER :: TmpCoord(:)
    REAL(KIND=dp),ALLOCATABLE :: Angles(:)
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,ind,Nmax,Nmin,Nfii,Nnodes,MaxElemNodes,NElems
    LOGICAL :: Found, Hit0, Hit90, Hit180, Hit270, SetDegOffset
    LOGICAL :: GotNormal, GotCenter, MoveAngle

    ! We choose degrees as they are more intuitive
    rad2deg = 180.0_dp / PI
    MaxElemNodes = BMesh2 % MaxElementNodes 
    ALLOCATE( Angles(MaxElemNodes) )
    
    Nnodes = BMesh2 % NumberOfNodes
    NElems = BMesh2 % NumberOfBulkElements
    FullCircle = .FALSE.

    ! Cylindrical projector is fitted always and rotational only when requested.
    IF( ListGetLogical( BParams,'Rotational Projector Center Fit',Found ) .OR. &
       Cylindrical ) THEN
      IF( .NOT. ListCheckPresent( BParams,'Rotational Projector Center X') ) THEN
        CALL CylinderFit( BMesh1, BParams ) 
      END IF
    END IF
    
    x0(1) = ListGetCReal( BParams,'Rotational Projector Center X',GotCenter ) 
    x0(2) = ListGetCReal( BParams,'Rotational Projector Center Y',Found ) 
    GotCenter = GotCenter .OR. Found
    x0(3) = ListGetCReal( BParams,'Rotational Projector Center Z',Found ) 
    GotCenter = GotCenter .OR. Found

    Normal(1) = ListGetCReal( BParams,'Rotational Projector Normal X',GotNormal ) 
    Normal(2) = ListGetCReal( BParams,'Rotational Projector Normal Y',Found ) 
    GotNormal = GotNormal .OR. Found
    Normal(3) = ListGetCReal( BParams,'Rotational Projector Normal Z',Found ) 
    GotNormal = GotNormal .OR. Found

    IF( GotNormal ) THEN
      CALL TangentDirections( Normal,Tangent1,Tangent2 )
    END IF

    ! Go through master (k=1) and target mesh (k=2)
    !--------------------------------------------
    DO k=1,2

      ! Potentially the projector may be set to rotate by just adding an offset 
      ! to the angle. This may depende on time etc. 
      IF( k == 1 ) THEN
        DegOffset = ListGetCReal(BParams,'Rotational Projector Angle Offset',SetDegOffset ) 
      ELSE
        SetDegOffset = .FALSE.
      END IF

      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      ! Check the initial bounding boxes
      !---------------------------------------------------------------------------
      x2_min(1) = MINVAL( PMesh % Nodes % x )
      x2_min(2) = MINVAL( PMesh % Nodes % y )
      x2_min(3) = MINVAL( PMesh % Nodes % z )
      
      x2_max(1) = MAXVAL( PMesh % Nodes % x )
      x2_max(2) = MAXVAL( PMesh % Nodes % y )
      x2_max(3) = MAXVAL( PMesh % Nodes % z )
      
      IF( k == 1 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Initial extrema for this boundary (x,y,z)',Level=8)
      ELSE IF( k == 2 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Initial extrema for target boundary (x,y,z)',Level=8)
      END IF
      DO i=1,3
        WRITE(Message,'(A,I0,A,2ES12.3)') 'Coordinate ',i,': ',x2_min(i),x2_max(i)
        CALL Info('RotationalInterfaceMeshes',Message,Level=8)    
      END DO

      ! Memorize the bounding box of the master mesh
      !--------------------------------------------------------------------------
      IF( k == 1 ) THEN
        x1_min = x2_min
        x1_max = x2_max
      END IF

      ! Do the actual coordinate transformation
      !---------------------------------------------------------------------------
      n = PMesh % NumberOfNodes
      DO i=1,n
        x(1) = PMesh % Nodes % x(i)
        x(2) = PMesh % Nodes % y(i)
        x(3) = PMesh % Nodes % z(i)

        ! Subtract the center of axis
        IF( GotCenter ) THEN
          x = x - x0
        END IF

        IF( GotNormal ) THEN
          xtmp = x
          x(1) = SUM( Tangent1 * xtmp ) 
          x(2) = SUM( Tangent2 * xtmp ) 
          x(3) = SUM( Normal * xtmp ) 
        END IF


        ! Set the angle to be the first coordinate as it may sometimes be the 
        ! only nonzero coordinate. Z-coordinate is always unchanged. 
        !------------------------------------------------------------------------
        alpha = rad2deg * ATAN2( x(2), x(1)  ) 
        rad = SQRT( x(1)**2 + x(2)**2)

        ! Set the offset and revert then the angle to range [-180,180] 
        IF( SetDegOffset ) THEN
          alpha = MODULO( alpha + DegOffset, 360.0_dp )            
          IF( alpha > 180.0_dp ) alpha = alpha - 360.0
        END IF

        PMesh % Nodes % x(i) = alpha
        PMesh % Nodes % y(i) = x(3)
        PMesh % Nodes % z(i) = rad      
      END DO
      

      ! For cylindrical projector follow exactly the same logic for slave and master
      !------------------------------------------------------------------------------
      IF( Cylindrical .AND. k == 2 ) THEN
        IF( MoveAngle ) THEN
          CALL Info('RotationalInterfaceMeshes','Moving the 2nd mesh discontinuity to same angle',Level=6)
          DO j=1,PMesh % NumberOfNodes
            IF( PMesh % Nodes % x(j) < Fii0 ) PMesh % Nodes % x(j) = &
                PMesh % Nodes % x(j) + 360.0_dp
          END DO
        END IF
      ELSE
        ! Let's see if we have a full angle to operate or not.
        ! If not, then make the interval continuous. 
        ! Here we check only four critical angles: (0,90,180,270) degs.
        Hit0 = .FALSE.; Hit90 = .FALSE.; Hit180 = .FALSE.; Hit270 = .FALSE.
        MoveAngle = .FALSE.; Fii = 0.0_dp; Fii0 = 0.0_dp
        
        DO i=1, PMesh % NumberOfBulkElements
          Element => PMesh % Elements(i)
          n = Element % TYPE % NumberOfNodes        
          NodeIndexes => Element % NodeIndexes
          Angles(1:n) = PMesh % Nodes % x(NodeIndexes)
          
          fmin = MINVAL( Angles(1:n) ) 
          fmax = MAXVAL( Angles(1:n) )
          
          IF( fmax - fmin > 180.0_dp ) THEN
            Hit180 = .TRUE.
          ELSE
            IF( fmax >= 0.0 .AND. fmin <= 0.0 ) Hit0 = .TRUE.
            IF( fmax >= 90.0_dp .AND. fmin <= 90.0_dp ) Hit90 = .TRUE.
            IF( fmax >= -90.0_dp .AND. fmin <= -90.0_dp ) Hit270 = .TRUE.
          END IF
        END DO
        FullCircle = Hit0 .AND. Hit90 .AND. Hit180 .AND. Hit270
        
        ! Eliminate the problematic discontinuity in case we have no full circle
        ! The discontinuity will be moved to some of angles (-90,0,90).
        IF( FullCircle ) THEN
          CALL Info('RotationalInterfaceMeshes','Cylindrical interface seems to be a full circle',&
              Level=6)
        ELSE IF( Hit180 ) THEN
          MoveAngle = .TRUE.
          IF( .NOT. Hit0 ) THEN
            Fii = 0.0_dp
          ELSE IF( .NOT. Hit270 ) THEN
            Fii = -90.0_dp
          ELSE IF( .NOT. Hit90 ) THEN
            Fii = 90.0_dp
          END IF

          DO j=1,PMesh % NumberOfNodes
            IF( PMesh % Nodes % x(j) < Fii ) PMesh % Nodes % x(j) = &
                PMesh % Nodes % x(j) + 360.0_dp
          END DO
          WRITE( Message,'(A,F8.3)') 'Moving discontinuity of angle to: ',Fii
          Fii0 = Fii
          CALL Info('RotationalInterfaceMesh',Message,Level=6)
        END IF
      END IF


      ! Check the transformed bounding boxes
      !---------------------------------------------------------------------------
      x2r_min(1) = MINVAL( PMesh % Nodes % x )
      x2r_min(2) = MINVAL( PMesh % Nodes % y )
      x2r_min(3) = MINVAL( PMesh % Nodes % z )
      
      x2r_max(1) = MAXVAL( PMesh % Nodes % x )
      x2r_max(2) = MAXVAL( PMesh % Nodes % y )
      x2r_max(3) = MAXVAL( PMesh % Nodes % z )
      
      IF( k == 1 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Transformed extrema for this boundary (phi,z,r)',Level=8)
      ELSE IF( k == 2 ) THEN
        CALL Info('RotationalInterfaceMeshes',&
            'Transformed extrema for target boundary (phi,z,r)',Level=8)
      END IF
      DO i=1,3
        WRITE(Message,'(A,I0,A,2ES12.3)') 'Coordinate ',i,': ',x2r_min(i),x2r_max(i)
        CALL Info('RotationalInterfaceMeshes',Message,Level=8)    
      END DO

      IF( x2r_min(3) < EPSILON( Radius ) ) THEN
        CALL Fatal('RotationalInterfaceMeshes','Radius cannot be almost zero!')
      END IF

      ! Memorize the bounding box for the 1st mesh
      IF( k == 1 ) THEN
        x1r_min = x2r_min
        x1r_max = x2r_max
      END IF
    END DO

    eps_rad = 1.0d-3 

    ! Choose radius to be max radius of this boundary
    Radius = x1r_max(3) 
    
    err1 = ( x1r_max(3) - x1r_min(3) ) / Radius
    err2 = ( x2r_max(3) - x2r_min(3) ) / Radius

    WRITE(Message,'(A,ES12.3)') 'Discrepancy from constant radius:',err1
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    WRITE(Message,'(A,ES12.3)') 'Discrepancy from constant radius:',err2
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    IF( err1 > eps_rad .OR. err2 > eps_rad ) THEN
      CALL Warn('RotationalInterfaceMeshes','Discrepancy of radius is rather large!')
    END IF

    ! Ok, so we have concluded that the interface has constant radius
    ! therefore the constant radius may be removed from the mesh description.
    ! Or perhaps we don't remove to allow more intelligent projector building 
    ! for contact mechanics. 
    !---------------------------------------------------------------------------
    !Bmesh1 % Nodes % z = 0.0_dp
    !BMesh2 % Nodes % z = 0.0_dp

    ! Check whether the z-coordinate is constant or not.
    ! Constant z-coordinate implies 1D system, otherwise 2D system.
    !---------------------------------------------------------------------------
    err1 = ( x1r_max(2) - x1r_min(2) ) / Radius
    err2 = ( x2r_max(2) - x2r_min(2) ) / Radius
    
    IF( err1 < eps_rad .AND. err2 < eps_rad ) THEN
      CALL Info('RotationalInterfaceMeshes','The effective interface meshes are 1D')
      Bmesh1 % Nodes % y = 0.0_dp
      Bmesh2 % Nodes % y = 0.0_dp
    ELSE
      CALL Info('RotationalInterfaceMeshes','The effective interface meshes are 2D')
    END IF

    ! Some pieces of the code cannot work with 1D meshes, this choice is ok for all steps
    Bmesh1 % MeshDim = 2
    Bmesh2 % MeshDim = 2      

    ! Cylindrical interface does not have symmetry as does the rotational!
    IF( Cylindrical .OR. FullCircle ) RETURN

    ! If were are studying a symmetric segment then anylyze further the angle 
    !-------------------------------------------------------------------------
    dFii1 = x1r_max(1)-x1r_min(1)
    dFii2 = x2r_max(1)-x2r_min(1)

    WRITE(Message,'(A,ES12.3)') 'This boundary dfii:  ',dFii1
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    WRITE(Message,'(A,ES12.3)') 'Target boundary dfii:  ',dFii2
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)    

    err1 = 2 * ABS( dFii1 - dFii2 ) / ( dFii1 + dFii2 )
    WRITE(Message,'(A,ES12.3)') 'Discrepancy in dfii:',err1
    CALL Info('RotationalInterfaceMeshes',Message,Level=8)        

    i = ListGetInteger(BParams,'Rotational Projector Periods',Found ) 
    IF( .NOT. Found ) THEN
      Nsymmetry = 360.0_dp / dFii2 
      WRITE(Message,'(A,ES12.3)') 'Suggested sections in target:',Nsymmetry
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
      IF( ABS( Nsymmetry - NINT( Nsymmetry ) ) < 0.01 .OR. Nsymmetry < 1.5 ) THEN          
        CALL Info('RotationalINterfaceMeshes','Assuming number of periods: '&
            //TRIM(I2S(NINT(Nsymmetry))),Level=8)
      ELSE
        IF( dFii1 < dFii2 ) THEN
          CALL Info('RotationalInterfaceMeshes','You might try to switch master and target!',Level=3)
        END IF
        CALL Fatal('RotationalInterfaceMeshes','Check your settings, this cannot be periodic!')
      END IF
      CALL ListAddInteger(BParams,'Rotational Projector Periods', NINT( Nsymmetry ) ) 
    ELSE
      WRITE(Message,'(A,I0)') 'Using enforced number of periods: ',i
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
      Nsymmetry = 360.0_dp / dFii2 
      WRITE(Message,'(A,ES12.3)') 'Suggested number of periods:',Nsymmetry
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
    END IF

  END SUBROUTINE RotationalInterfaceMeshes
!------------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Given axial projectors compute the number of cycles.
  !---------------------------------------------------------------------------
  SUBROUTINE AxialInterfaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    TYPE(Element_t), POINTER :: Element
    REAL(KIND=dp) :: minalpha, maxalpha, minalpha2, maxalpha2
    REAL(KIND=dp) :: x(3), xcyl(3),rad2deg,F1min,F1max,F2min,F2max,dFii, dFii1,dFii2,eps_rad,&
        err1,err2,dF,Nsymmetry,rad,alpha,x0(3),xtmp(3), maxrad, &
        Normal(3), Tangent1(3), Tangent2(3) 
    REAL(KIND=dp), POINTER :: TmpCoord(:)
    REAL(KIND=dp),ALLOCATABLE :: Angles(:)
    INTEGER, POINTER :: NodeIndexes(:)
    INTEGER :: i,j,k,n,ind,Nmax,Nmin,Nfii,Nnodes,MaxElemNodes,sweep
    LOGICAL :: Found, Hit0, Hit90, Hit180, Hit270
    LOGICAL :: GotNormal, GotCenter, FullCircle

    ! We choose degrees as they are more intuitive
    rad2deg = 180.0_dp / PI
    MaxElemNodes = BMesh2 % MaxElementNodes 
    
    x0(1) = ListGetCReal( BParams,'Axial Projector Center X',GotCenter ) 
    x0(2) = ListGetCReal( BParams,'Axial Projector Center Y',Found ) 
    GotCenter = GotCenter .OR. Found
    x0(3) = ListGetCReal( BParams,'Axial Projector Center Z',Found ) 
    GotCenter = GotCenter .OR. Found

    Normal(1) = ListGetCReal( BParams,'Axial Projector Normal X',GotNormal ) 
    Normal(2) = ListGetCReal( BParams,'Axial Projector Normal Y',Found ) 
    GotNormal = GotNormal .OR. Found
    Normal(3) = ListGetCReal( BParams,'Axial Projector Normal Z',Found ) 
    GotNormal = GotNormal .OR. Found

    IF( GotNormal ) THEN
      CALL TangentDirections( Normal,Tangent1,Tangent2 )
    ELSE
      CALL Info('AxialInterfaceMeshes',&
          'Assuming axial interface to have z-axis the normal!',Level=8)
    END IF

    ! Go through master (k=1) and target mesh (k=2)
    !--------------------------------------------
    FullCircle = .FALSE.

    DO k=1,2
     
      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      ! Do the actual coordinate transformation
      !---------------------------------------------------------------------------
      n = PMesh % NumberOfNodes

      ! Register the hit in basic quadrants
      Hit0 = .FALSE.; Hit90 = .FALSE.; Hit180 = .FALSE.; Hit270 = .FALSE.
      maxrad = 0.0_dp
      minalpha = HUGE( minalpha ); maxalpha = -HUGE(maxalpha)
      minalpha2 = HUGE( minalpha2 ); maxalpha2 = -HUGE(maxalpha2)
      
      ! 1st sweep only find max radius, 2nd sweep register the angle range
      DO sweep = 1, 2
        DO i=1,n
          x(1) = PMesh % Nodes % x(i)
          x(2) = PMesh % Nodes % y(i)
          x(3) = PMesh % Nodes % z(i)
          
          ! Subtract the center of axis
          IF( GotCenter ) x = x - x0
          
          IF( GotNormal ) THEN
            xtmp = x
            x(1) = SUM( Tangent1 * xtmp ) 
            x(2) = SUM( Tangent2 * xtmp ) 
            x(3) = SUM( Normal * xtmp ) 
          END IF
          
          ! Compute the angle
          !------------------------------------------------------------------------
          rad = SQRT( x(1)**2 + x(2)**2)
          
          IF( sweep == 1 ) THEN
            maxrad = MAX( maxrad, rad ) 
            CYCLE
          END IF

          ! Do the logic for large enough radius
          IF( rad < 0.5_dp * maxrad ) CYCLE

          IF( x(1) > 0.0 .AND. ABS(x(2)) < ABS(x(1)) ) Hit0 = .TRUE.
          IF( x(2) > 0.0 .AND. ABS(x(1)) < ABS(x(2)) ) Hit90 = .TRUE.
          IF( x(1) < 0.0 .AND. ABS(x(2)) < ABS(x(1)) ) Hit180 = .TRUE.
          IF( x(2) < 0.0 .AND. ABS(x(1)) < ABS(x(2)) ) Hit270 = .TRUE.
          
          ! This can compute the range if there is no nodes close to discontinuity at 180 degs
          alpha = rad2deg * ATAN2( x(2), x(1)  ) 
          minalpha = MIN( alpha, minalpha ) 
          maxalpha = MAX( alpha, maxalpha ) 

          ! This eliminates the discontinuity and moves it to 0 degs
          IF( alpha < 0.0_dp ) alpha = alpha + 360.0_dp          
          minalpha2 = MIN( alpha, minalpha2 ) 
          maxalpha2 = MAX( alpha, maxalpha2 ) 
        END DO
      END DO
      
      FullCircle = Hit0 .AND. Hit90 .AND. Hit180 .AND. Hit270
      IF( FullCircle ) THEN
        CALL Info('RotationalInterfaceMeshes','Axial interface seems to be a full circle',&
            Level=6)
        EXIT
      END IF
      
      dFii = MIN( maxalpha2 - minalpha2, maxalpha - minalpha ) 

      ! memorize the max angle for 1st boundary mesh
      IF( k == 1 ) THEN
        WRITE(Message,'(A,ES12.3)') 'This boundary dfii: ',dFii
        dFii1 = dFii
      ELSE
        WRITE(Message,'(A,ES12.3)') 'Target boundary dfii: ',dFii
        dFii2 = dFii
      END IF
      CALL Info('AxialInterfaceMeshes',Message,Level=8)    
    END DO

    IF( FullCircle ) THEN
      Nsymmetry = 1.0_dp
    ELSE
      err1 = 2 * ABS( dFii1 - dFii2 ) / ( dFii1 + dFii2 )
      WRITE(Message,'(A,ES12.3)') 'Discrepancy in dfii:',err1
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
      Nsymmetry = 360.0_dp / ( MIN( dfii1, dfii2 ) ) 
    END IF
    
    WRITE(Message,'(A,ES12.3)') 'Suggested number of periods:',Nsymmetry
    CALL Info('AxialInterfaceMeshes',Message,Level=8)        

    i = ListGetInteger(BParams,'Axial Projector Periods',Found ) 
    IF( .NOT. Found ) THEN
      CALL ListAddInteger(BParams,'Axial Projector Periods', NINT( Nsymmetry ) ) 
    ELSE
      WRITE(Message,'(A,I0)') 'Using enforced number of periods: ',i
      CALL Info('RotationalInterfaceMeshes',Message,Level=8)        
    END IF

  END SUBROUTINE AxialInterfaceMeshes
!------------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !> Given two interface meshes for nonconforming radial boundaries make 
  !> a coordinate transformation to (r,z) level.
  !> This is always a symmetry condition and can not be a contact condition.
  !---------------------------------------------------------------------------
  SUBROUTINE RadialInterfaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: PMesh
    REAL(KIND=dp) :: x1_min(3),x1_max(3),x2_min(3),x2_max(3), x(3), r, phi, z, &
        err1, err2, phierr, eps_rad, rad, rad2deg
    INTEGER :: i,j,k

    ! We choose degrees as they are more intuitive
    rad2deg = 180.0_dp / PI

    ! Go through master (k=1) and target mesh (k=2)
    !--------------------------------------------
    DO k=1,2

      IF( k == 1 ) THEN
        PMesh => BMesh1
      ELSE
        PMesh => BMesh2
      END IF

      x2_min = HUGE( x2_min )
      x2_max = -HUGE( x2_max )
      
      ! Loop over all nodes
      !----------------------------------------------------------------------------
      DO i=1,PMesh % NumberOfNodes
        x(1) = PMesh % Nodes % x(i)
        x(2) = PMesh % Nodes % y(i)
        x(3) = PMesh % Nodes % z(i)
        
        ! Do the actual coordinate transformation
        !---------------------------------------------------------------------------
        r = SQRT( x(1)**2 + x(2)**2 )
        phi = rad2deg * ATAN2( x(2), x(1)  )
        z = x(3)

        !PRINT *,'interface node:',k,i,r,phi,x(1:2)
        
        PMesh % Nodes % x(i) = r
        PMesh % Nodes % y(i) = z
        PMesh % Nodes % z(i) = 0.0_dp

        ! This is just to check a posteriori that the ranges are ok
        x2_min(1) = MIN(r,x2_min(1))
        IF( r > EPSILON( r ) ) THEN
          x2_min(2) = MIN(phi,x2_min(2))
        END IF
        x2_min(3) = MIN(z,x2_min(3))

        x2_max(1) = MAX(r,x2_max(1))
        IF( r > EPSILON(r) ) THEN
          x2_max(2) = MAX(phi,x2_max(2))
        END IF
        x2_max(3) = MAX(z,x2_max(3))
      END DO

      ! Memorize the bounding box of the master mesh
      !--------------------------------------------------------------------------
      IF( k == 1 ) THEN
        x1_min = x2_min
        x1_max = x2_max
      END IF

      IF( k == 1 ) THEN
        CALL Info('RadialInterfaceMeshes',&
            'Transformed extrema for this boundary (r,phi,z)',Level=8)
      ELSE IF( k == 2 ) THEN
        CALL Info('RadialInterfaceMeshes',&
            'Transformed extrema for target boundary (r,phi,z)',Level=8)
      END IF

      DO i=1,3
        WRITE(Message,'(A,I0,A,2ES12.3)') 'Coordinate ',i,': ',x2_min(i),x2_max(i)
        CALL Info('RadialInterfaceMeshes',Message,Level=8)    
      END DO

      phierr = x2_max(2) - x2_min(2)  
      WRITE(Message,'(A,ES12.3)') 'Discrepancy from constant angle (degs):',phierr
      CALL Info('RadialInterfaceMeshes',Message,Level=8)    
    END DO

    ! Error in radius
    ! Choose radius to be max radius of either boundary
    rad = MAX( x1_max(1), x2_max(1) )    
    err1 = ABS( x1_max(1) - x2_max(1) ) / rad
    err2 = ABS( x1_min(1) - x2_min(1) ) / rad

    WRITE(Message,'(A,ES12.3)') 'Discrepancy in maximum radius:',err1
    CALL Info('RadialInterfaceMeshes',Message,Level=8)    

    WRITE(Message,'(A,ES12.3)') 'Discrepancy in minimum radius:',err2
    CALL Info('RadialInterfaceMeshes',Message,Level=8)    

    eps_rad = 1.0d-3
    IF( err1 > eps_rad .OR. err2 > eps_rad ) THEN
      CALL Warn('RadialInterfaceMeshes','Discrepancy of radius may be too large!')
    END IF

    ! Some pieces of the code cannot work with 1D meshes, this choice is ok for all steps
    Bmesh1 % MeshDim = 2
    Bmesh2 % MeshDim = 2      
    
  END SUBROUTINE RadialInterfaceMeshes
!------------------------------------------------------------------------------

  !---------------------------------------------------------------------------
  !> Given two interface meshes flatten them to (x,y) plane.
  !---------------------------------------------------------------------------
  SUBROUTINE FlatInterfaceMeshes(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Bmesh
    INTEGER :: FlatDim, MeshDim, MinDiffI, i, j
    REAL(KIND=dp), POINTER CONTIG :: Coord(:)
    REAL(KIND=dp) :: Diff, MaxDiff, MinDiff, RelDiff, RelDiff1
    LOGICAL :: Found, ReduceDim

    CALL Info('FlatInterfaceMeshes','Flattening interface meshes to 2D',Level=8)    
    
    MeshDim = CurrentModel % Dimension
    FlatDim = ListGetInteger( BParams,'Flat Projector Coordinate',Found,minv=1,maxv=3) 
    ReduceDim = ListGetLogical( BParams,'Flat Projector Reduce Dimension',Found )

    IF(.NOT. Found ) THEN
      DO j=1, 2
        IF( j == 1 ) THEN
          Bmesh => BMesh1
        ELSE
          BMesh => BMesh2
        END IF
        
        MaxDiff = 0.0
        MinDiff = HUGE( MinDiff ) 
        
        DO i = 1, MeshDim
          IF( i == 1 ) THEN
            Coord => BMesh % Nodes % x 
          ELSE IF( i == 2 ) THEN
            Coord => Bmesh % Nodes % y
          ELSE
            Coord => Bmesh % Nodes % z
          END IF

          Diff = MAXVAL( Coord ) - MINVAL( Coord )
          MaxDiff = MAX( Diff, MaxDiff ) 
          IF( Diff < MinDiff ) THEN
            MinDiff = Diff
            MinDiffI = i
          END IF
        END DO

        RelDiff = MinDiff / MaxDiff
        IF( j == 1 ) THEN
          FlatDim = MinDiffI
          RelDiff1 = RelDiff 
        ELSE IF( j == 2 ) THEN
          IF( RelDiff < RelDiff1 ) FlatDim = MinDiffI
        END IF
      END DO

      CALL Info('FlatInterfaceMeshes','> Flat Projector Coordinate < set to: '//TRIM(I2S(FlatDim)))
      CALL ListAddInteger( BParams,'Flat Projector Coordinate',FlatDim )
    END IF


    DO j=1,2
      ! Some pieces of the code cannot work with 1D meshes, this choice is ok for all steps
      IF( j == 1 ) THEN
        Bmesh => BMesh1
      ELSE
        BMesh => BMesh2
      END IF

      ! Set the 3rd component to be the "distance" in the flat interface      
      IF( FlatDim == 3 ) THEN
        CONTINUE
      ELSE IF( FlatDim == 2 ) THEN
        Coord => BMesh % Nodes % y
        BMesh % Nodes % y => BMesh % Nodes % z
        BMesh % Nodes % z => Coord
        IF( MeshDim == 2 ) BMesh % Nodes % y = 0.0_dp
      ELSE IF( FlatDim == 1 ) THEN
        Coord => BMesh % Nodes % x
        BMesh % Nodes % x => BMesh % Nodes % y
        BMesh % Nodes % y => BMesh % Nodes % z
        Bmesh % Nodes % z => Coord
        IF( MeshDim == 2 ) BMesh % Nodes % y = 0.0_dp
      END IF

      IF( ReduceDim ) BMesh % Nodes % z = 0.0_dp

      Bmesh % MeshDim = 2
    END DO

  END SUBROUTINE FlatInterfaceMeshes
!------------------------------------------------------------------------------


  !---------------------------------------------------------------------------
  !> Given two interface meshes flatten them into the plane that 
  !> best fits either of the meshes. 
  !---------------------------------------------------------------------------
  SUBROUTINE PlaneInterfaceMeshes(BMesh1, BMesh2, BParams )
    !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Bmesh
    INTEGER :: i, j, n, nip, MeshDim
    REAL(KIND=dp) :: Normal(3), NormalSum(3), RefSum, Length, Planeness, &
        PlaneNormal(3,1), PlaneNormal1(3,1), Planeness1, Normal1(3), &
        Tangent(3), Tangent2(3), Coord(3), detJ, Normal0(3)
    REAL(KIND=dp), POINTER :: PNormal(:,:), Basis(:)
    TYPE(Element_t), POINTER :: Element
    TYPE(GaussIntegrationPoints_t) :: IP
    TYPE(Nodes_t) :: ElementNodes
    INTEGER, POINTER :: NodeIndexes(:)
    LOGICAL :: Found, Stat, Normal0Set

    CALL Info('PlaneInterfaceMeshes','Flattening interface meshes to a plane',Level=8)    

    MeshDim = CurrentModel % Dimension
    PNormal => ListGetConstRealArray( BParams,'Plane Projector Normal',Found) 

    ! If the projector normal is not given determine it first 
    IF(.NOT. Found ) THEN     
      CALL Info('PlaneInterfaceMeshes','Could not find > Plane Projector Normal < so determining it now',Level=12)    

      n = MAX_ELEMENT_NODES
      ALLOCATE( ElementNodes % x(n), ElementNodes % y(n), ElementNodes % z(n), Basis(n) )
      ElementNodes % x = 0; ElementNodes % y = 0; ElementNodes % z = 0

      ! Fit a plane to both datasets
      DO j=1, 2
        IF( j == 1 ) THEN
          Bmesh => BMesh1
        ELSE      
          BMesh => BMesh2
        END IF

        NormalSum = 0.0_dp
        RefSum = 0.0_dp
        Normal0Set = .FALSE.

        ! we use the Dot2Min and Normal2 temporarily also for first mesh, with k=1
        !-------------------------------------------------------------------------
        DO i=1, BMesh % NumberOfBulkElements
          Element => BMesh % Elements(i)
          n = Element % TYPE % NumberOfNodes
          NodeIndexes => Element % NodeIndexes
          IP = GaussPoints( Element ) 

          ElementNodes % x(1:n) = BMesh % Nodes % x(NodeIndexes(1:n))
          ElementNodes % y(1:n) = BMesh % Nodes % y(NodeIndexes(1:n))
          ElementNodes % z(1:n) = BMesh % Nodes % z(NodeIndexes(1:n))           

          DO nip=1, IP % n 
            stat = ElementInfo( Element,ElementNodes,&
                IP % u(nip),IP % v(nip),IP % w(nip),detJ,Basis)

            Normal = NormalVector( Element, ElementNodes, &
                IP % u(nip), IP % v(nip), .FALSE. ) 
            IF( .NOT. Normal0Set ) THEN
              Normal0 = Normal
              Normal0Set = .TRUE.
            END IF

            IF( SUM( Normal * Normal0 ) < 0.0 ) Normal = -Normal

            NormalSum = NormalSum + IP % S(nip) * DetJ * Normal
            RefSum = RefSum + IP % S(nip) * DetJ
          END DO
        END DO

        ! Normalize the normal to unity length
        Length = SQRT( SUM( NormalSum ** 2 ) )
        PlaneNormal(:,1) = NormalSum / Length

        ! Planeness is one if all the normals have the same direction
        Planeness = Length / RefSum 
        
        ! Save the key parameters of the first mesh
        IF( j == 1 ) THEN          
          PlaneNormal1 = PlaneNormal
          Planeness1 = Planeness
        END IF
      END DO

      ! Choose the mesh for which is close to a plane 
      IF( Planeness1 > Planeness ) THEN
        PRINT *,'PlaneNormal: Selecting slave normal'
        PlaneNormal = PlaneNormal1
      ELSE
        PRINT *,'PlaneNormal: Selecting master normal'        
        PlaneNormal = -PlaneNormal
      END IF

      PRINT *,'PlaneNormal selected:',PlaneNormal(:,1)

      CALL ListAddConstRealArray( BParams,'Plane Projector Normal',&
          3,1,PlaneNormal )
      DEALLOCATE( ElementNodes % x, ElementNodes % y, ElementNodes % z, Basis )

      PNormal => ListGetConstRealArray( BParams,'Plane Projector Normal',Found) 
    END IF

    Normal = Pnormal(1:3,1)
    CALL TangentDirections( Normal, Tangent, Tangent2 )

    IF(.FALSE.) THEN
      PRINT *,'Normal:',Normal
      PRINT *,'Tangent1:',Tangent
      PRINT *,'Tangent2:',Tangent2
    END IF

    DO j=1,2
      IF( j == 1 ) THEN
        Bmesh => BMesh1
      ELSE
        BMesh => BMesh2
      END IF

      DO i=1,BMesh % NumberOfNodes
        Coord(1) = BMesh % Nodes % x(i)
        Coord(2) = BMesh % Nodes % y(i)
        Coord(3) = BMesh % Nodes % z(i)

        BMesh % Nodes % x(i) = SUM( Coord * Tangent )
        IF( MeshDim == 3 ) THEN
          BMesh % Nodes % y(i) = SUM( Coord * Tangent2 )
        ELSE
          BMesh % Nodes % y(i) = 0.0_dp
        END IF
        BMesh % Nodes % z(i) = SUM( Coord * Normal ) 
      END DO

      IF(.FALSE.) THEN
        PRINT *,'Range for mesh:',j
        PRINT *,'X:',MINVAL(BMesh % Nodes % x),MAXVAL(BMesh % Nodes % x)
        PRINT *,'Y:',MINVAL(BMesh % Nodes % y),MAXVAL(BMesh % Nodes % y)
        PRINT *,'Z:',MINVAL(BMesh % Nodes % z),MAXVAL(BMesh % Nodes % z)
      END IF
    END DO

    Bmesh % MeshDim = 2

  END SUBROUTINE PlaneInterfaceMeshes
  !------------------------------------------------------------------------------



  !---------------------------------------------------------------------------
  !> Given a permutation map the (x,y,z) such that the projector can better 
  !> be applied. E.g. if boundary has constant x, take that as the last coordinate.
  !---------------------------------------------------------------------------
  SUBROUTINE MapInterfaceCoordinate(BMesh1, BMesh2, BParams )
  !---------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
    TYPE(Valuelist_t), POINTER :: BParams
    !--------------------------------------------------------------------------
    LOGICAL :: Found
    REAL(KIND=dp), POINTER CONTIG:: NodesX(:), NodesY(:), NodesZ(:), Wrk(:,:)
    INTEGER, POINTER :: CoordMap(:)
    INTEGER :: MeshNo
    TYPE(Mesh_t), POINTER :: BMesh
    
    ! Perform coordinate mapping
    !------------------------------------------------------------
    CoordMap => ListGetIntegerArray( BParams, & 
        'Projector Coordinate Mapping',Found )
    IF( .NOT. Found ) RETURN

    CALL Info('MapInterfaceCoordinates','Performing coordinate mapping',Level=8)
    
    IF ( SIZE( CoordMap ) /= 3 ) THEN
      WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
      CALL Error( 'MapInterfaceCoordinates', Message )
      WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
      CALL Fatal( 'MapInterfaceCoordinates', Message )
    END IF
    
    IF ( ALL( CoordMap(1:3) /= 1 ) .OR. ALL( CoordMap(1:3) /= 2 ) .OR. ALL( CoordMap(1:3) /= 3 ) ) THEN
      WRITE( Message, * ) 'Inconsistent Coordinate Mapping: ', CoordMap
      CALL Error( 'MapInterfaceCoordinates', Message )
      WRITE( Message, * ) 'Coordinate mapping should be a permutation of 1,2 and 3'
      CALL Fatal( 'MapInterfaceCoordinates', Message )
    END IF

    DO MeshNo = 1,2
      IF( MeshNo == 1 ) THEN
        BMesh => BMesh1
      ELSE
        BMesh => BMesh2 
      END IF

      IF( CoordMap(1) == 1 ) THEN
        NodesX => BMesh % Nodes % x
      ELSE IF( CoordMap(1) == 2 ) THEN
        NodesX => BMesh % Nodes % y
      ELSE
        NodesX => BMesh % Nodes % z
      END IF
    
      IF( CoordMap(2) == 1 ) THEN
        NodesY => BMesh % Nodes % x
      ELSE IF( CoordMap(2) == 2 ) THEN
        NodesY => BMesh % Nodes % y
      ELSE
        NodesY => BMesh % Nodes % z
      END IF
      
      IF( CoordMap(3) == 1 ) THEN
        NodesZ => BMesh % Nodes % x
      ELSE IF( CoordMap(3) == 2 ) THEN
        NodesZ => BMesh % Nodes % y
      ELSE
        NodesZ => BMesh % Nodes % z
      END IF

      BMesh % Nodes % x => NodesX
      BMesh % Nodes % y => NodesY
      BMesh % Nodes % z => NodesZ
    END DO

  END SUBROUTINE MapInterfaceCoordinate

  !------------------------------------------------------------------------------
!> Create a projector between Master and Target boundaries.
!> The projector may be a nodal projector x=Px or a weigted 
!> Galerking projector such that Qx=Px. In the first case the projector 
!> will be P and in the second case [Q-P]. 
!------------------------------------------------------------------------------
  FUNCTION PeriodicProjector( Model, Mesh, This, Trgt, cdim, &
      Galerkin ) RESULT(Projector)
!------------------------------------------------------------------------------   
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    INTEGER, OPTIONAL :: cdim
    TYPE(Mesh_t), TARGET :: Mesh
    TYPE(Matrix_t), POINTER :: Projector
    LOGICAL, OPTIONAL :: Galerkin
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,dim
    LOGICAL :: GotIt, UseQuadrantTree, Success, WeakProjector, &
        Rotational, AntiRotational, Sliding, AntiSliding, Repeating, AntiRepeating, &
        Discontinuous, NodalJump, Radial, AntiRadial, DoNodes, DoEdges, Axial, AntiAxial, &
        Flat, Plane, AntiPlane, LevelProj, FullCircle, Cylindrical, &
        ParallelNumbering, TimestepNumbering, EnforceOverlay, NormalProj
    LOGICAL, ALLOCATABLE :: MirrorNode(:)
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    TYPE(Nodes_t), POINTER :: MeshNodes, GaussNodes
    REAL(KIND=dp) :: NodeScale, EdgeScale, Radius, Coeff
    TYPE(ValueList_t), POINTER :: BC
    CHARACTER(LEN=MAX_NAME_LEN) :: FilePrefix
    TYPE(Variable_t), POINTER :: v

    INTERFACE
      FUNCTION WeightedProjector(BMesh2, BMesh1, InvPerm2, InvPerm1, &
          UseQuadrantTree, Repeating, AntiRepeating, PeriodicScale, &
          NodalJump ) &
         RESULT ( Projector )
        USE Types
        TYPE(Mesh_t), POINTER :: BMesh1, BMesh2
        REAL(KIND=dp) :: PeriodicScale
        INTEGER, POINTER :: InvPerm1(:), InvPerm2(:)
        LOGICAL :: UseQuadrantTree, Repeating, AntiRepeating
        TYPE(Matrix_t), POINTER :: Projector
        LOGICAL :: NodalJump
      END FUNCTION WeightedProjector



    END INTERFACE
!------------------------------------------------------------------------------
    Projector => NULL()
    IF ( This <= 0  ) RETURN    
    CALL Info('PeriodicProjector','Starting projector creation',Level=12)

    DIM = CoordinateSystemDimension()

    CALL ResetTimer('PeriodicProjector')
    
    Projector => NULL()
    BC => Model % BCs(This) % Values
    PMesh => Mesh

    
    ! Whether to choose nodal or Galerkin projector is determined by an optional
    ! flag. The default is the nodal projector.
    !--------------------------------------------------------------------------
    IF( PRESENT( Galerkin) ) THEN
      WeakProjector = Galerkin
    ELSE
      WeakProjector = ListGetLogical( BC, 'Galerkin Projector', GotIt )
    END IF

    ! If the boundary is discontinuous then we have the luxury of creating the projector
    ! very cheaply using the permutation vector. This does not need the target as the 
    ! boundary is self-contained.
    !------------------------------------------------------------------------------------
    IF( ListGetLogical( BC, 'Discontinuous Boundary', GotIt ) .AND. Mesh % DisContMesh )THEN
      IF( WeakProjector ) THEN
        Projector => WeightedProjectorDiscont( PMesh, This )
      ELSE
        Projector => NodalProjectorDiscont( PMesh, This )
      END IF
      
      IF ( .NOT. ASSOCIATED( Projector ) ) RETURN
      GOTO 100
    END IF
    
    IF ( Trgt <= 0 ) RETURN    

    ! Create the mesh projector, and if needed, also eliminate the ghost nodes
    ! There are two choices of projector: a nodal projector P in x=Px, and a 
    ! Galerkin projector [Q-P] in Qx=Px. 
    ! The projector is assumed to be either a rotational projector with no translation
    ! and rotation, or then generic one with possible coordinate mapping.
    !---------------------------------------------------------------------------------
    CALL Info('PeriodicProjector','-----------------------------------------------------',Level=8)
    WRITE( Message,'(A,I0,A,I0)') 'Creating projector between BCs ',This,' and ',Trgt
    CALL Info('PeriodicProjector',Message,Level=8)

    ! Create temporal mesh structures that are utilized when making the 
    ! projector between "This" and "Trgt" boundary.
    !--------------------------------------------------------------------------
    BMesh1 => AllocateMesh()
    BMesh2 => AllocateMesh()
    
    CALL CreateInterfaceMeshes( Model, Mesh, This, Trgt, Bmesh1, BMesh2, &
        Success ) 

    IF(.NOT. Success) THEN
      CALL ReleaseMesh(BMesh1)
      CALL ReleaseMesh(BMesh2)
      RETURN
    END IF

    ! If requested map the interface coordinate from (x,y,z) to any permutation of these. 
    CALL MapInterfaceCoordinate( BMesh1, BMesh2, Model % BCs(This) % Values )

    NormalProj = ListGetLogical( BC,'Normal Projector',GotIt )
    
    ! Check whether to use (anti)rotational projector.
    ! We don't really know on which side the projector was called so 
    ! let's check both sides.
    !--------------------------------------------------------------------------
    Rotational = ListGetLogical( BC,'Rotational Projector',GotIt )
    AntiRotational = ListGetLogical( BC,'Anti Rotational Projector',GotIt )
    IF( AntiRotational ) Rotational = .TRUE.

    Cylindrical =  ListGetLogical( BC,'Cylindrical Projector',GotIt )

    Radial = ListGetLogical( BC,'Radial Projector',GotIt )
    AntiRadial = ListGetLogical( BC,'Anti Radial Projector',GotIt )
    IF( AntiRadial ) Radial = .TRUE.

    Axial = ListGetLogical( BC,'Axial Projector',GotIt )
    AntiAxial = ListGetLogical( BC,'Anti Axial Projector',GotIt )
    IF( AntiAxial ) Axial = .TRUE.

    Sliding = ListGetLogical( BC,'Sliding Projector',GotIt )
    AntiSliding = ListGetLogical( BC,'Anti Sliding Projector',GotIt )
    IF( AntiSliding ) Sliding = .TRUE. 

    Flat = ListGetLogical( BC,'Flat Projector',GotIt )
    Plane = ListGetLogical( BC, 'Plane Projector',GotIt )
    AntiPlane = ListGetLogical( BC,'Anti Plane Projector',GotIt )    
    IF( AntiPlane ) Plane = .TRUE.
    
    IF( Radial ) CALL Info('PeriodicProjector','Enforcing > Radial Projector <',Level=12)
    IF( Axial ) CALL Info('PeriodicProjector','Enforcing > Axial Projector <',Level=12)
    IF( Sliding ) CALL Info('PeriodicProjector','Enforcing > Sliding Projector <',Level=12)
    IF( Cylindrical ) CALL Info('PeriodicProjector','Enforcing > Cylindrical Projector <',Level=12)
    IF( Rotational ) CALL Info('PeriodicProjector','Enforcing > Rotational Projector <',Level=12)
    IF( Flat ) CALL Info('PeriodicProjector','Enforcing > Flat Projector <',Level=12)
    IF( Plane ) CALL Info('PeriodicProjector','Enforcing > Plane Projector <',Level=12)

    NodeScale = ListGetConstReal( BC, 'Mortar BC Scaling',GotIt)
    IF(.NOT.Gotit ) THEN
      IF( AntiRadial .OR. AntiPlane ) THEN
        NodeScale = -1._dp
      ELSE
        NodeScale = 1.0_dp
      END IF
    END IF
    EdgeScale = NodeScale

    NodalJump = ListCheckPrefix( BC,'Mortar BC Coefficient')
    IF(.NOT. NodalJump ) THEN
      NodalJump = ListCheckPrefix( BC,'Mortar BC Resistivity')
    END IF

    ! There are tailored projectors for simplified interfaces
    !-------------------------------------------------------------

    ! Stride projector is obsolete and has been eliminated.
    IF( ListGetLogical( BC,'Stride Projector',GotIt) ) THEN
      CALL ListAddLogical( BC,'Level Projector',.TRUE.)
      CALL ListAddLogical( BC,'Level Projector Strong',.TRUE.)
      CALL Warn('PeriodicProjector','Enforcing > Level Projector < instead of old > Stride Projector <')
    END IF

    LevelProj = ListGetLogical( BC,'Level Projector',GotIt) 
    IF( Rotational .OR. Cylindrical .OR. Radial .OR. Flat .OR. Plane .OR. Axial ) THEN
      IF(.NOT. GotIt ) THEN
        CALL Info('PeriodicProjector','Enforcing > Level Projector = True < with dimensional reduction',&
            Level = 7 )
        LevelProj = .TRUE. 
      ELSE IF(.NOT. LevelProj ) THEN
        ! If we have dimensionally reduced projector but don't use LevelProjector 
        ! to integrate over it, then ensure that the 3rd coordinate is set to zero.
        BMesh1 % Nodes % z = 0.0_dp
        BMesh2 % Nodes % z = 0.0_dp
      END IF
    END IF


    IF( LevelProj ) THEN
      IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',GotIt ) ) THEN
        DoNodes = .FALSE.
      ELSE
        IF( ListGetLogical( BC,'Projector Skip Nodes',GotIt) ) THEN
          DoNodes = .FALSE.
        ELSE
          DoNodes = ( Mesh % NumberOfNodes > 0 ) 
        END IF
      END IF

      IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',GotIt ) ) THEN
        DoEdges = .FALSE.
      ELSE
        IF( ListGetLogical( BC,'Projector Skip Edges',GotIt) ) THEN
          DoEdges = .FALSE.
        ELSE
          ! We are conservative here since there may be edges in 2D which 
          ! still cannot be used for creating the projector
          DoEdges = ( Mesh % NumberOfEdges > 0 .AND. &
              Mesh % MeshDim == 3 .AND. Dim == 3 )

          ! Ensure that there is no p-elements that made us think that we have edges
          ! Here we assume that if there is any p-element then also the 1st element is such
          IF( DoEdges ) THEN
            IF(isPelement(Mesh % Elements(1))) THEN
              DoEdges = .FALSE.
              CALL Info('PeriodicProjector','Edge projector will not be created for p-element mesh',Level=10)
            END IF
          END IF
        END IF
      END IF
    END IF


    ! If the interface is rotational move to (phi,z) plane and alter the phi coordinate
    ! so that the meshes coincide.
    ! Otherwise make the two meshes to coincide using rotation, translation &
    ! scaling.
    !---------------------------------------------------------------------------------
    Radius = 1.0_dp
    FullCircle = .FALSE.
    EnforceOverlay = ListGetLogical( BC, 'Mortar BC enforce overlay', GotIt )

    IF( Rotational .OR. Cylindrical ) THEN
      CALL RotationalInterfaceMeshes( BMesh1, BMesh2, BC, Cylindrical, &
          Radius, FullCircle )
    ELSE IF( Radial ) THEN
      CALL RadialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Flat ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Axial ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )      
      CALL AxialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Plane ) THEN
      CALL PlaneInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( .NOT. ( Sliding .OR. NormalProj ) ) THEN
      IF( .NOT. GotIt ) EnforceOverlay = .TRUE.
    END IF

    IF( EnforceOverlay ) THEN
      CALL OverlayIntefaceMeshes( BMesh1, BMesh2, BC )
    END IF

    Repeating = ( Rotational .OR. Sliding .OR. Axial ) .AND. .NOT. FullCircle 
    AntiRepeating = .FALSE.
    IF( Repeating ) THEN
      AntiRepeating = ListGetLogical( BC,'Antisymmetric BC',GotIt ) 
      IF( .NOT. GotIt ) THEN
        AntiRepeating = ( AntiRotational .OR. AntiSliding .OR. AntiAxial ) .AND. .NOT. FullCircle 
      END IF
    END IF
      
    IF( LevelProj ) THEN 
      Projector => LevelProjector( BMesh1, BMesh2, Repeating, AntiRepeating, &
          FullCircle, Radius, DoNodes, DoEdges, &          
          NodeScale, EdgeScale, BC )
    ELSE
      IF( FullCircle ) THEN
        CALL Fatal('PeriodicProjector','A full circle cannot be dealt with the generic projector!')
      END IF

      UseQuadrantTree = ListGetLogical(Model % Simulation,'Use Quadrant Tree',GotIt)
      IF( .NOT. GotIt ) UseQuadrantTree = .TRUE.
      
      IF( NormalProj ) THEN
        Projector => NormalProjector( BMesh2, BMesh1, BC )
      ELSE IF( WeakProjector ) THEN
        Projector => WeightedProjector( BMesh2, BMesh1, BMesh2 % InvPerm, BMesh1 % InvPerm, &
            UseQuadrantTree, Repeating, AntiRepeating, NodeScale, NodalJump )
      ELSE
        Projector => NodalProjector( BMesh2, BMesh1, &
            UseQuadrantTree, Repeating, AntiRepeating )
      END IF
    END IF


    ! Deallocate mesh structures:
    !---------------------------------------------------------------
    BMesh1 % Projector => NULL()
    BMesh1 % Parent => NULL()
    !DEALLOCATE( BMesh1 % InvPerm ) 
    CALL ReleaseMesh(BMesh1)

    BMesh2 % Projector => NULL()
    BMesh2 % Parent => NULL()
    !DEALLOCATE( BMesh2 % InvPerm ) 
    CALL ReleaseMesh(BMesh2)

100 Projector % ProjectorBC = This

    IF( ListGetLogical( BC,'Projector Set Rowsum',GotIt ) ) THEN
      CALL SetProjectorRowsum( Projector )
    END IF

    Coeff = ListGetConstReal( BC,'Projector Multiplier',GotIt) 
    IF(.NOT. GotIt) Coeff = ListGetConstReal( Model % Simulation,&
        'Projector Multiplier',GotIt) 
    IF( GotIt ) Projector % Values = Coeff * Projector % Values

    IF( ListGetLogical( BC,'Save Projector',GotIt ) ) THEN
      ParallelNumbering = ListGetLogical( BC,'Save Projector Global Numbering',GotIt )

      FilePrefix = 'p'//TRIM(I2S(This))
      
      TimestepNumbering = ListGetLogical( BC,'Save Projector Timestep Numbering',GotIt )
      IF( TimestepNumberIng ) THEN
        i = 0
        v => VariableGet( Mesh % Variables, 'timestep' )
        IF( ASSOCIATED( v ) ) i = NINT( v % Values(1) )
        WRITE( FilePrefix,'(A,I4.4)') TRIM(FilePrefix)//'_',i
      END IF
        
      CALL SaveProjector( Projector, .TRUE.,TRIM(FilePrefix), &
          Parallel = ParallelNumbering) 
      
      ! Dual projector if it exists
      IF( ASSOCIATED( Projector % Ematrix ) ) THEN
        CALL SaveProjector( Projector % Ematrix, .TRUE.,'dual_'//TRIM(FilePrefix),&
            Projector % InvPerm, Parallel = ParallelNumbering) 
      END IF

      ! Biorthogonal projector if it exists
      IF( ASSOCIATED( Projector % Child ) ) THEN
        CALL SaveProjector( Projector % Child, .TRUE.,'biortho_'//TRIM(FilePrefix), &
            Projector % InvPerm, Parallel = ParallelNumbering ) 
      END IF

      IF( ListGetLogical( BC,'Save Projector And Stop',GotIt ) ) STOP EXIT_OK
    END IF    

    CALL CheckTimer('PeriodicProjector',Delete=.TRUE.)
    CALL Info('PeriodicProjector','Projector created, now exiting...',Level=8)

!------------------------------------------------------------------------------
  END FUNCTION PeriodicProjector
!------------------------------------------------------------------------------


  

!------------------------------------------------------------------------------
!> Create a permutation between two meshes such that we can solve a smaller system.
!------------------------------------------------------------------------------
  SUBROUTINE PeriodicPermutation( Model, Mesh, This, Trgt, PerPerm, PerFlip ) 
!------------------------------------------------------------------------------   
    TYPE(Model_t) :: Model
    INTEGER :: This, Trgt
    TYPE(Mesh_t), TARGET :: Mesh
    INTEGER, POINTER :: PerPerm(:)
    LOGICAL, POINTER :: PerFlip(:)
!------------------------------------------------------------------------------
    INTEGER :: i,j,k,n,dim
    LOGICAL :: GotIt, Success, Rotational, AntiRotational, Sliding, AntiSliding, Repeating, &
        Radial, AntiRadial, DoNodes, DoEdges, Axial, AntiAxial, &
        Flat, Plane, AntiPlane, Cylindrical, ParallelNumbering, EnforceOverlay, &
        FullCircle, AntiPeriodic
    REAL(KIND=dp) :: Radius
    TYPE(Mesh_t), POINTER ::  BMesh1, BMesh2, PMesh
    TYPE(ValueList_t), POINTER :: BC
    
!------------------------------------------------------------------------------
    IF ( This <= 0  .OR. Trgt <= 0 ) RETURN    
    CALL Info('PeriodicPermutation','Starting periodic permutation creation',Level=12)

    CALL ResetTimer('PeriodicPermutation')
    
    DIM = CoordinateSystemDimension()
    BC => Model % BCs(This) % Values
    PMesh => Mesh
    
    CALL Info('PeriodicPermutation','-----------------------------------------------------',Level=8)
    WRITE( Message,'(A,I0,A,I0)') 'Creating mapping between BCs ',This,' and ',Trgt
    CALL Info('PeriodicPermutation',Message,Level=8)

    BMesh1 => AllocateMesh()
    BMesh2 => AllocateMesh()
    
    CALL CreateInterfaceMeshes( Model, Mesh, This, Trgt, Bmesh1, BMesh2, Success ) 
    
    IF(.NOT. Success) THEN
      CALL ReleaseMesh(BMesh1)
      CALL ReleaseMesh(BMesh2)
      RETURN
    END IF

    ! If requested map the interface coordinate from (x,y,z) to any permutation of these. 
    CALL MapInterfaceCoordinate( BMesh1, BMesh2, Model % BCs(This) % Values )
    
    ! Lets check what kind of symmetry we have.
    Rotational = ListGetLogical( BC,'Rotational Projector',GotIt )
    AntiRotational = ListGetLogical( BC,'Anti Rotational Projector',GotIt )

    Cylindrical =  ListGetLogical( BC,'Cylindrical Projector',GotIt )

    Radial = ListGetLogical( BC,'Radial Projector',GotIt )
    AntiRadial = ListGetLogical( BC,'Anti Radial Projector',GotIt )
    IF( AntiRadial ) Radial = .TRUE.
    
    Axial = ListGetLogical( BC,'Axial Projector',GotIt )
    AntiAxial = ListGetLogical( BC,'Anti Axial Projector',GotIt )
    IF( AntiAxial ) Axial = .TRUE.
    
    Sliding = ListGetLogical( BC, 'Sliding Projector',GotIt )
    AntiSliding = ListGetLogical( BC, 'Anti Sliding Projector',GotIt )
    IF( AntiSliding ) Sliding = .TRUE.
    
    Flat = ListGetLogical( BC, 'Flat Projector',GotIt )
    Plane = ListGetLogical( BC, 'Plane Projector',GotIt )
    AntiPlane = ListGetLogical( BC,'Anti Plane Projector',GotIt )    
    IF( AntiPlane ) Plane = .TRUE.

    AntiPeriodic = ListGetLogical( BC,'Antisymmetric BC',GotIt )
    IF( .NOT. GotIt ) THEN   
      AntiPeriodic = ( AntiRotational .OR. AntiRadial .OR. AntiAxial .OR. AntiPlane ) 
    END IF
      
    IF( AntiPeriodic ) CALL Info('PeriodicPermutation','Assuming antiperiodic conforming projector',Level=8)
    
    IF( Radial ) CALL Info('PeriodicProjector','Enforcing > Radial Projector <',Level=12)
    IF( Axial ) CALL Info('PeriodicProjector','Enforcing > Axial Projector <',Level=12)
    IF( Sliding ) CALL Info('PeriodicProjector','Enforcing > Sliding Projector <',Level=12)
    IF( Cylindrical ) CALL Info('PeriodicProjector','Enforcing > Cylindrical Projector <',Level=12)
    IF( Rotational ) CALL Info('PeriodicProjector','Enforcing > Rotational Projector <',Level=12)
    IF( Flat ) CALL Info('PeriodicProjector','Enforcing > Flat Projector <',Level=12)
    IF( Plane ) CALL Info('PeriodicProjector','Enforcing > Plane Projector <',Level=12)

    DoNodes = .TRUE.
    IF( ListGetLogical( Model % Solver % Values,'Projector Skip Nodes',GotIt ) ) DoNodes = .FALSE.    
    IF( ListGetLogical( BC,'Projector Skip Nodes',GotIt) ) DoNodes = .FALSE.

    ! We are conservative here since there may be edges in 2D which 
    ! still cannot be used for creating the projector
    DoEdges = ( Mesh % NumberOfEdges > 0 .AND. Mesh % MeshDim == 3 .AND. Dim == 3 )
    
    ! Ensure that there is no p-elements that made us think that we have edges
    ! Here we assume that if there is any p-element then also the 1st element is such
    IF( DoEdges ) THEN
      IF(isPelement(Mesh % Elements(1))) THEN
        DoEdges = .FALSE.
        CALL Info('PeriodicProjector','Edge projector will not be created for p-element mesh',Level=10)
      END IF
    END IF
        
    IF( ListGetLogical( Model % Solver % Values,'Projector Skip Edges',GotIt ) ) DoEdges = .FALSE.
    IF( ListGetLogical( BC,'Projector Skip Edges',GotIt) ) DoEdges = .FALSE.
      
    ! Make the two meshes to coincide using rotation, translation scaling.
    !---------------------------------------------------------------------------------
    Radius = 1.0_dp
    EnforceOverlay = ListGetLogical( BC, 'Mortar BC enforce overlay', GotIt )

    IF( Rotational .OR. Cylindrical ) THEN
      CALL RotationalInterfaceMeshes( BMesh1, BMesh2, BC, Cylindrical, &
          Radius, FullCircle )
      IF( FullCircle ) CALL Fatal('PeriodicPermutation','Cannot deal full circle with permutation')
    ELSE IF( Radial ) THEN
      CALL RadialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Flat ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Axial ) THEN
      CALL FlatInterfaceMeshes( BMesh1, BMesh2, BC )      
      CALL AxialInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( Plane ) THEN
      CALL PlaneInterfaceMeshes( BMesh1, BMesh2, BC )
    ELSE IF( .NOT. Sliding ) THEN
      IF( .NOT. GotIt ) EnforceOverlay = .TRUE.
    END IF
    
    IF( EnforceOverlay ) THEN
      CALL OverlayIntefaceMeshes( BMesh1, BMesh2, BC )
    END IF
    
    IF( DoNodes ) CALL ConformingNodePerm(PMesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
    IF( DoEdges ) CALL ConformingEdgePerm(PMesh, BMesh1, BMesh2, PerPerm, PerFlip, AntiPeriodic )
        
    ! Deallocate mesh structures:
    !---------------------------------------------------------------
    BMesh1 % Projector => NULL()
    BMesh1 % Parent => NULL()
    !DEALLOCATE( BMesh1 % InvPerm ) 
    CALL ReleaseMesh(BMesh1)

    BMesh2 % Projector => NULL()
    BMesh2 % Parent => NULL()
    !DEALLOCATE( BMesh2 % InvPerm ) 
    CALL ReleaseMesh(BMesh2)

    CALL CheckTimer('PeriodicPermutation',Delete=.TRUE.)
           
    CALL Info('PeriodicPermutation','Periodic permutation created, now exiting...',Level=8)
   
    
!------------------------------------------------------------------------------
  END SUBROUTINE PeriodicPermutation
!------------------------------------------------------------------------------


  
  !> If periodic BCs given, compute boundary mesh projector.
  !> If conforming BCs given, create permutation for elimination.
  !------------------------------------------------------
  SUBROUTINE GeneratePeriodicProjectors( Model, Mesh ) 
    TYPE(Model_t) :: Model
    TYPE(Mesh_t), POINTER :: Mesh
  !------------------------------------------------------
    INTEGER :: i,j,k,n,nocyclic,noconf,noflip,mini,maxi
    LOGICAL :: Found
    INTEGER, POINTER :: PerPerm(:)
    LOGICAL, POINTER :: PerFlip(:)
    
    DO i = 1,Model % NumberOfBCs
      k = ListGetInteger( Model % BCs(i) % Values, 'Periodic BC', Found )
      IF( Found ) THEN
        Model % BCs(i) % PMatrix => PeriodicProjector( Model, Mesh, i, k )
      END IF
    END DO

    IF( ListCheckPresentAnyBC( Model,'Conforming BC' ) ) THEN
      IF(.NOT. ASSOCIATED( Mesh % PeriodicPerm ) ) THEN
        n = Mesh % NumberOfNodes + Mesh % NumberOfEdges
        ALLOCATE( Mesh % PeriodicPerm(n) )
        ALLOCATE( Mesh % PeriodicFlip(n) )
      END IF
      PerPerm => Mesh % PeriodicPerm      
      PerPerm = 0
      PerFlip => Mesh % PeriodicFlip
      PerFlip = .FALSE.
      DO i = 1,Model % NumberOfBCs
        k = ListGetInteger( Model % BCs(i) % Values, 'Conforming BC', Found )
        IF( Found ) THEN
          CALL PeriodicPermutation( Model, Mesh, i, k, PerPerm, PerFlip )
        END IF
      END DO
      nocyclic = 0
      noconf = 0
      mini = HUGE(mini)
      maxi = 0
      
      DO i = 1,n
        j = PerPerm(i)
        IF( j > 0 ) THEN
          mini = MIN( mini, i )
          maxi = MAX( maxi, i )
          noconf = noconf + 1
          IF( PerPerm(j) > 0 ) THEN
            PerPerm(i) = PerPerm(j)
            IF( PerFlip(i) ) THEN
              PerFlip(i) = .NOT. PerFlip(j)
            ELSE
              PerFlip(i) = PerFlip(j)
            END IF
            nocyclic = nocyclic + 1
          END IF
        END IF
      END DO
      noflip = COUNT( PerFlip )
            
      CALL Info('GeneratePeriodicProjectors','Number of conforming maps: '//TRIM(I2S(noconf)),Level=8)
      IF(nocyclic>0) CALL Info('GeneratePeriodicProjectors','Number of cyclic maps: '//TRIM(I2S(nocyclic)),Level=8)
      IF(noflip>0) CALL Info('GeneratePeriodicProjectors','Number of periodic flips: '//TRIM(I2S(noflip)),Level=8)
    END IF

    
  END SUBROUTINE GeneratePeriodicProjectors

!------------------------------------------------------------------------------
END MODULE MortarProjector
!------------------------------------------------------------------------------
