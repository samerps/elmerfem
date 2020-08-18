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
!>  Internal mesh generation utilities.
!------------------------------------------------------------------------------


MODULE MeshGeneration

  USE ElementUtils
  USE ElementDescription
  USE MeshUtils 
  USE Interpolation
  USE ParallelUtils
  USE Types
  IMPLICIT NONE

CONTAINS

!------------------------------------------------------------------------------
!> Create node distribution for a unit segment x \in [0,1] with n elements 
!> i.e. n+1 nodes. There are different options for the type of distribution.
!> 1) Even distribution 
!> 2) Geometric distribution
!> 3) Arbitrary distribution determined by a functional dependence
!> Note that the 3rd algorithm involves iterative solution of the nodal
!> positions and is therefore not bullet-proof.
!------------------------------------------------------------------------------
  SUBROUTINE UnitSegmentDivision( w, n, ExtList )
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    INTEGER :: n
    TYPE(ValueList_t), POINTER, OPTIONAL :: ExtList
    !---------------------------------------------------------------
    INTEGER :: i,J,iter,maxiter
    REAL(KIND=dp) :: q,r,h1,hn,minhn,err_eps,err,xn
    REAL(KIND=dp), ALLOCATABLE :: wold(:),h(:)
    LOGICAL :: Found, GotRatio, FunExtruded, Fun1D
    TYPE(Nodes_t) :: Nodes
    TYPE(ValueList_t), POINTER :: ParList
    CHARACTER(*), PARAMETER :: Caller='UnitSegmentDivision'
    
    IF( PRESENT( ExtList ) ) THEN
      ParList => ExtList
    ELSE
      ParList => CurrentModel % Simulation
    END IF

    FunExtruded = ListCheckPresent( ParList,'Extruded Mesh Density')
    Fun1D = ListCheckPresent( ParList,'1D Mesh Density')
    
    ! Geometric division
    !---------------------------------------------------------------
    q = ListGetConstReal( ParList,'Extruded Mesh Ratio',GotRatio)
    IF(.NOT. GotRatio) q = ListGetConstReal( ParList,'1D Mesh Ratio',GotRatio)
    IF( GotRatio ) THEN
      IF( ( ABS(ABS(q)-1.0_dp) < 1.0d-6 ) .OR. (q < 0.0_dp .AND. n <= 2) ) THEN
        CALL Info(Caller,'Assuming linear division as mesh ratio is close to one!')
        GotRatio = .FALSE.
      END IF
    END IF
    
    IF( GotRatio ) THEN
      CALL Info(Caller,'Creating geometric division',Level=5)

      IF( q > 0.0_dp ) THEN      
        r = q**(1.0_dp/(n-1))
        h1 = (1-r)/(1-r**n)
        w(0) = 0.0_dp
        DO i=1,n-1
          w(i) = h1 * (1-r**i)/(1-r)
        END DO
        w(n) = 1.0_dp
      ELSE
        q = -q
        IF(MODULO(n,2) == 0) THEN
          r = q**(1.0_dp/(n/2-1))
          h1 = 0.5_dp*(1-r)/(1-r**(n/2))
        ELSE 
          r = q**(1.0_dp/((n-1)/2))
          h1 = 0.5_dp / ( (1-r**((n+1)/2))/(1-r) - 0.5_dp * r**((n-1)/2))
        END IF
        
        w(0) = 0.0_dp
        DO i=1,n
          IF( i <= n/2 ) THEN
            w(i) = h1 * (1-r**i)/(1-r)
          ELSE
            w(i) = 1.0_dp -  h1 * (1-r**(n-i))/(1-r)
          END IF
        END DO
        w(n) = 1.0_dp
      END IF
            
    ! Generic division given by a function
    !-----------------------------------------------------------------------
    ELSE IF( FunExtruded .OR. Fun1D ) THEN

      CALL Info(Caller,'Creating functional division',Level=5)

      ! Initial guess is an even distribution
      DO i=0,n
        w(i) = i/(1._dp * n)
      END DO

      ALLOCATE( wold(0:n),h(1:n))
      wold = w

      ! parameters that determine the accuracy of the iteration
      maxiter = 10000
      err_eps = 1.0d-6

      ! Iterate to have a density distribution
      !---------------------------------------
      DO iter=1,maxiter
        
        minhn = HUGE(minhn)
        wold = w

        ! Compute the point in the local mesh xn \in [0,1]  
        ! and get the mesh parameter for that element from
        ! external function.
        !---------------------------------------------------
        DO i=1,n
          xn = (w(i)+w(i-1))/2.0_dp
          minhn = MIN( minhn, w(i)-w(i-1) )
          IF( FunExtruded ) THEN
            h(i) = ListGetFun( ParList,'Extruded Mesh Density', xn )
          ELSE
            h(i) = ListGetFun( ParList,'1D Mesh Density', xn )
          END IF
          IF( h(i) < EPSILON( h(i) ) ) THEN
            CALL Fatal(Caller,'Given value for h(i) was negative!')
          END IF
        END DO

        ! Utilize symmetric Gauss-Seidel to compute the new positions, w(i).
        ! from a weigted mean of the desired elemental densities, h(i).
        ! Note that something more clever could be applied here. 
        ! This was just a first implementation...
        !-------------------------------------------------------------
        DO i=1,n-1
          w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
        END DO
        DO i=n-1,1,-1
          w(i) = (w(i-1)*h(i+1)+w(i+1)*h(i))/(h(i)+h(i+1))
        END DO
        
        ! If the maximum error is small compared to the minimum elementsize then exit
        !-----------------------------------------------------------------------------
        err = MAXVAL( ABS(w-wold))/minhn

        IF( err < err_eps ) THEN
          WRITE( Message, '(A,I0,A)') 'Convergence obtained in ',iter,' iterations'
          CALL Info(Caller, Message, Level=9 )
          EXIT
        END IF
      END DO

      IF( iter > maxiter ) THEN
        CALL Warn(Caller,'No convergence obtained for the unit mesh division!')
      END IF

    ! Uniform division 
    !--------------------------------------------------------------
    ELSE
      CALL Info(Caller,'Creating linear division',Level=5)
      DO i=0,n     
        w(i) = i/(1._dp * n)
      END DO
    END IF
    
    CALL Info(Caller,'Mesh division ready',Level=9)
    DO i=0,n
      WRITE( Message, '(A,I0,A,ES12.4)') 'w(',i,') : ',w(i)
      CALL Info(Caller, Message, Level=9 )
    END DO

  END SUBROUTINE UnitSegmentDivision
!------------------------------------------------------------------------------


  
!------------------------------------------------------------------------------
!> Given a 2D mesh extrude it to be 3D. The 3rd coordinate will always
!> be at the interval [0,1]. Therefore the adaptation for different shapes
!> must be done with StructuredMeshMapper, or some similar utility. 
!> The top and bottom surface will be assigned Boundary Condition tags
!> with indexes one larger than the maximum used on by the 2D mesh. 
!------------------------------------------------------------------------------
  FUNCTION MeshExtrude(Mesh_in, in_levels, ExtrudedMeshName) RESULT(Mesh_out)
!------------------------------------------------------------------------------
    TYPE(Mesh_t), POINTER :: Mesh_in, Mesh_out
    INTEGER :: in_levels
    CHARACTER(LEN=MAX_NAME_LEN),INTENT(IN),OPTIONAL :: ExtrudedMeshName

!------------------------------------------------------------------------------
    INTEGER :: i,j,k,l,n,cnt,cnt101,ind(8),max_baseline_bid,max_bid,l_n,max_body,bcid,&
        ExtrudedCoord,dg_n,totalnumberofelements
    TYPE(ParallelInfo_t), POINTER :: PI_in, PI_out
    INTEGER :: nnodes,gnodes,gelements,ierr
    LOGICAL :: isParallel, Found, NeedEdges, PreserveBaseline, PreserveEdges, &
        Rotational, Rotate2Pi
    REAL(KIND=dp)::w,MinCoord,MaxCoord,CurrCoord
    REAL(KIND=dp), POINTER :: ActiveCoord(:)
    REAL(KIND=dp), ALLOCATABLE :: Wtable(:)
    CHARACTER(*), PARAMETER :: Caller='MeshExtrude'
!------------------------------------------------------------------------------

    CALL Info(Caller,'Creating '//TRIM(I2S(in_levels+1))//' extruded element layers',Level=10)

    Mesh_out => AllocateMesh()

    isParallel = ParEnv % PEs>1

    ! Generate volume nodal points:
    ! -----------------------------
    n=Mesh_in % NumberOfNodes
    nnodes=(in_levels+2)*n
    gnodes = nnodes

    ALLOCATE( Mesh_out % Nodes % x(nnodes) )
    ALLOCATE( Mesh_out % Nodes % y(nnodes) )
    ALLOCATE( Mesh_out % Nodes % z(nnodes) )

    gelements = Mesh_in % NumberOfBulkElements

    IF (isParallel) THEN
      PI_in  => Mesh_in % ParallelInfo
      PI_out => Mesh_out % ParallelInfo
    
      IF(.NOT. ASSOCIATED( PI_in ) ) CALL Fatal(Caller,'PI_in not associated!')
      IF(.NOT. ASSOCIATED( PI_out ) ) CALL Fatal(Caller,'PI_out not associated!')
            
      ALLOCATE(PI_out % NeighbourList(nnodes))
      ALLOCATE(PI_out % INTERFACE(nnodes))
      ALLOCATE(PI_out % GlobalDOFs(nnodes))

      IF(.NOT. ASSOCIATED( PI_in % NeighbourList ) ) THEN
        CALL Fatal(Caller,'Neighnours not associated!')
      END IF

      ! For unset neighbours just set the this partition to be the only owner
      DO i=1,Mesh_in % NumberOfNodes
        IF (.NOT.ASSOCIATED(PI_in % NeighbourList(i) % Neighbours)) THEN
          CALL AllocateVector(PI_in % NeighbourList(i) % Neighbours,1)
          PI_in % NeighbourList(i) % Neighbours(1) = ParEnv % Mype
        END IF
      END DO
          
      j=0
      DO i=1,Mesh_in % NumberOfNodes
        IF (PI_in % NeighbourList(i) % &
            Neighbours(1) == ParEnv % MyPE ) j=j+1
      END DO

      CALL MPI_ALLREDUCE(j,gnodes,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
      
      j=0
      DO i=1,Mesh_in % NumberOfBulkElements
        IF (Mesh_in % Elements(i) % PartIndex == ParEnv % MyPE) j=j+1
      END DO
      
      CALL MPI_ALLREDUCE(j,gelements,1, &
           MPI_INTEGER,MPI_SUM,ELMER_COMM_WORLD,ierr)
    END IF

    CALL Info(Caller,'Number of extruded nodes: '//TRIM(I2S(nnodes)),Level=12)
    CALL Info(Caller,'Number of extruded elements: '//TRIM(I2S(gelements)),Level=12)


    ! Create the division for the 1D unit mesh
    !--------------------------------------------
    ALLOCATE( Wtable( 0: in_levels + 1 ) )
    CALL UnitSegmentDivision( Wtable, in_levels + 1 ) 

    ExtrudedCoord = ListGetInteger( CurrentModel % Simulation,'Extruded Coordinate Index', &
        Found, minv=1,maxv=3 )
    IF(.NOT. Found) ExtrudedCoord = 3 

    IF( ExtrudedCoord == 1 ) THEN
      ActiveCoord => Mesh_out % Nodes % x
    ELSE IF( ExtrudedCoord == 2 ) THEN
      ActiveCoord => Mesh_out % Nodes % y
    ELSE IF( ExtrudedCoord == 3 ) THEN
      ActiveCoord => Mesh_out % Nodes % z
    END IF


    PreserveBaseline = ListGetLogical( CurrentModel % Simulation,'Preserve Baseline',Found )
    IF(.NOT. Found) PreserveBaseline = .FALSE.

    PreserveEdges = ListGetLogical( CurrentModel % Simulation,'Preserve Edges',Found )
    IF(.NOT. Found) PreserveEdges = .FALSE.

    MinCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Min Coordinate',Found )
    IF(.NOT. Found) MinCoord = 0.0_dp

    MaxCoord = ListGetConstReal( CurrentModel % Simulation,'Extruded Max Coordinate',Found )
    IF(.NOT. Found) MaxCoord = 1.0_dp

    Rotate2Pi = .FALSE.
    Rotational = ListGetLogical( CurrentModel % Simulation,'Extruded Mesh Rotational',Found )    
    IF( Rotational ) THEN
      Rotate2Pi = ( ABS(ABS( MaxCoord-MinCoord ) - 2*PI) < 1.0d-3*PI )
      IF( Rotate2Pi ) CALL Info(Caller,'Perfoming full 2Pi rotation',Level=6)
    END IF

    
    cnt=0
    DO i=0,in_levels+1

      ! If we rotate full 2Pi then we have natural closure!
      IF( Rotate2Pi ) THEN
        IF( i == in_levels+1) EXIT
      END IF
      
      w = Wtable( i ) 
      CurrCoord = w * MaxCoord + (1-w) * MinCoord      
      
      DO j=1,Mesh_in % NumberOfNodes

        cnt = cnt + 1

        Mesh_out % Nodes % x(cnt) = Mesh_in % Nodes % x(j) 
        Mesh_out % Nodes % y(cnt) = Mesh_in % Nodes % y(j) 
        Mesh_out % Nodes % z(cnt) = Mesh_in % Nodes % z(j) 

        ! Override the coordinate in the extruded direction by the value on the layer.
        ActiveCoord(cnt) = CurrCoord

        IF (isParallel) THEN
          PI_out % INTERFACE(cnt) = PI_in % INTERFACE(j)

          ALLOCATE(PI_out % NeighbourList(cnt) % Neighbours(&
               SIZE(PI_in % NeighbourList(j) % Neighbours)))
          PI_out % NeighbourList(cnt) % Neighbours = &
            PI_in % NeighbourList(j) % Neighbours

          PI_out % GlobalDOFs(cnt) = PI_in % GlobalDOFs(j)+i*gnodes
        END IF

      END DO
    END DO
    Mesh_out % NumberOfNodes=cnt

    
    IF( Rotational ) THEN
      BLOCK
        REAL(KIND=DP) :: x,y,z,r        
        DO i=1,cnt          
          x = Mesh_out % Nodes % x(i)
          y = Mesh_out % Nodes % y(i)
          z = Mesh_out % Nodes % z(i)

          Mesh_out % Nodes % x(i) = COS(z) * x
          Mesh_out % Nodes % y(i) = SIN(z) * x
          Mesh_out % Nodes % z(i) = y
        END DO
      END BLOCK
    END IF
    
    
    ! Count 101 elements:
    ! (these require an extra layer)
    ! -------------------

    cnt101 = 0
    DO i=Mesh_in % NumberOfBulkElements+1, &
         Mesh_in % NumberOfBulkElements+Mesh_in % NumberOfBoundaryElements
       IF(Mesh_in % Elements(i) % TYPE % ElementCode == 101) cnt101 = cnt101+1
    END DO

    n=SIZE(Mesh_in % Elements)

    ! inquire total number of needed 
    IF( Rotate2Pi ) THEN
      totalnumberofelements = n*(in_levels+1) + cnt101
    ELSE
      totalnumberofelements = n*(in_levels+3) + cnt101
    END IF

    IF (PreserveBaseline) &
        totalnumberofelements = totalnumberofelements + Mesh_in % NumberOfBoundaryElements
    ALLOCATE(Mesh_out % Elements(totalnumberofelements))
    
    ! Generate volume bulk elements:
    ! ------------------------------

    Mesh_out % MaxElementNodes = 0

    NeedEdges=.FALSE.
    n=Mesh_in % NumberOfNodes
    cnt=0; dg_n  = 0
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBulkElements

        cnt=cnt+1
        Mesh_out % Elements(cnt) = Mesh_in % Elements(j)

        l_n=0
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+i*n
        END DO
        DO k=1,Mesh_in % Elements(j) % TYPE % NumberOfNodes
          l_n=l_n+1
          IF( Rotate2Pi .AND. i==in_levels ) THEN
            ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)
          ELSE
            ind(l_n) = Mesh_in % Elements(j) % NodeIndexes(k)+(i+1)*n
          END IF
        END DO
        Mesh_out % Elements(cnt) % NDOFs = l_n
        Mesh_out % MaxElementNodes=MAX(Mesh_out % MaxElementNodes,l_n)

        SELECT CASE(l_n)
        CASE(6)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(706)
        CASE(8)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(808)
        END SELECT

        Mesh_out % Elements(cnt) % GElementIndex = &
             Mesh_in % Elements(j) % GelementIndex + gelements*i

        Mesh_out % Elements(cnt) % ElementIndex = cnt
        ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n)) 
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % NodeIndexes = ind(1:l_n)
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    END DO
    Mesh_out % NumberOfBulkElements=cnt

    max_bid=0
    max_baseline_bid=0

    ! include edges (see below)
    NeedEdges =  (NeedEdges .OR. PreserveEdges)
    
    ! -------------------------------------------------------
    IF (PreserveBaseline) THEN
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

        ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
        Mesh_out % Elements(cnt) % BoundaryInfo = &
           Mesh_in % Elements(k) % BoundaryInfo

        max_bid = MAX(max_bid, Mesh_in % Elements(k) % &
                BoundaryInfo % Constraint)

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in %  NumberOfBulkElements*(in_levels+1)+ &
	                   (in_levels+2)*Mesh_in % NumberOfBoundaryElements+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Right => &
              Mesh_out % Elements(Mesh_in % NumberOfBulkElements*(in_levels+1)+ &
	      (in_levels+2)*Mesh_in % NumberOfBoundaryElements+l)
        END IF

        IF(Mesh_in % Elements(k) % TYPE % ElementCode>=200) THEN
          Mesh_out % Elements(cnt) % NDOFs = 2
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(2)) 
          ind(1) = Mesh_in % Elements(k) % NodeIndexes(1)
          ind(2) = Mesh_in % Elements(k) % NodeIndexes(2)
          Mesh_out % Elements(cnt) % NodeIndexes = ind(1:2)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(202)
        ELSE
          Mesh_out % Elements(cnt) % NDOFs = 1
          l=SIZE(Mesh_in % Elements(k) % NodeIndexes)
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l))
          Mesh_out % Elements(cnt) % NodeIndexes = &
            Mesh_in % Elements(k) % NodeIndexes
          Mesh_out % Elements(cnt) % TYPE => &
             Mesh_in % Elements(k) % TYPE
        END IF
        Mesh_out % Elements(cnt) % DGDOFs = 0
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % ElementIndex = cnt
        Mesh_out % Elements(cnt) % PDefs => NULL()
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    
      IF(isParallel) THEN
        j=max_bid
        CALL MPI_ALLREDUCE(j,max_bid,1, &
            MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
      END IF

      max_baseline_bid = max_bid

    END IF


    ! Add side boundaries with the bottom mesh boundary id's:
    ! (or shift ids if preserving the baseline boundary)
    ! -------------------------------------------------------
    DO i=0,in_levels
      DO j=1,Mesh_in % NumberOfBoundaryElements
        k = j + Mesh_in % NumberOfBulkElements

        cnt=cnt+1

        Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

        ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
        Mesh_out % Elements(cnt) % BoundaryInfo = &
           Mesh_in % Elements(k) % BoundaryInfo

        Mesh_out % Elements(cnt) % BoundaryInfo % constraint = &
           Mesh_out % Elements(cnt) % BoundaryInfo % constraint + max_baseline_bid

        max_bid = MAX(max_bid, max_baseline_bid + &
           Mesh_in % Elements(k) % BoundaryInfo % Constraint)

        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Left)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Left % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF
        IF(ASSOCIATED(Mesh_in % Elements(k) % BoundaryInfo % Right)) THEN
          l=Mesh_in % Elements(k) % BoundaryInfo % Right % ElementIndex
          Mesh_out % Elements(cnt) % BoundaryInfo % Right => &
             Mesh_out % Elements(Mesh_in % NumberOfBulkElements*i+l)
        END IF

        IF(Mesh_in % Elements(k) % TYPE % ElementCode>=200) THEN
          Mesh_out % Elements(cnt) % NDOFs = 4
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(4)) 

          ind(1) = Mesh_in % Elements(k) % NodeIndexes(1)+i*n
          ind(2) = Mesh_in % Elements(k) % NodeIndexes(2)+i*n

          IF( Rotate2Pi .AND. i==in_levels ) THEN
            ind(3) = Mesh_in % Elements(k) % NodeIndexes(2)
            ind(4) = Mesh_in % Elements(k) % NodeIndexes(1)
          ELSE
            ind(3) = Mesh_in % Elements(k) % NodeIndexes(2)+(i+1)*n
            ind(4) = Mesh_in % Elements(k) % NodeIndexes(1)+(i+1)*n
          END IF
            Mesh_out % Elements(cnt) % NodeIndexes = ind(1:4)
          Mesh_out % Elements(cnt) % TYPE => GetElementType(404)
        ELSE
          Mesh_out % Elements(cnt) % NDOFs = 1
          l=SIZE(Mesh_in % Elements(k) % NodeIndexes)
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l))
          Mesh_out % Elements(cnt) % NodeIndexes = &
            Mesh_in % Elements(k) % NodeIndexes+i*n
          Mesh_out % Elements(cnt) % TYPE => &
             Mesh_in % Elements(k) % TYPE
        END IF 
        Mesh_out % Elements(cnt) % ElementIndex = cnt
        Mesh_out % Elements(cnt) % DGDOFs = 0
        Mesh_out % Elements(cnt) % DGIndexes => NULL()
        Mesh_out % Elements(cnt) % PDefs => NULL()
        Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
        Mesh_out % Elements(cnt) % FaceIndexes => NULL()
        Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
      END DO
    END DO

    !Take care of extra 101 elements
    !-------------------------------

    IF(cnt101 > 0) THEN
       DO j=1,Mesh_in % NumberOfBoundaryElements
          k = j + Mesh_in % NumberOfBulkElements

          IF(Mesh_in % Elements(k) % TYPE % ElementCode /= 101) CYCLE
          cnt=cnt+1

          Mesh_out % Elements(cnt) = Mesh_in % Elements(k)

          ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
          Mesh_out % Elements(cnt) % BoundaryInfo = &
               Mesh_in % Elements(k) % BoundaryInfo

          Mesh_out % Elements(cnt) % BoundaryInfo % constraint = &
               Mesh_out % Elements(cnt) % BoundaryInfo % constraint + max_baseline_bid

          max_bid = MAX(max_bid, max_baseline_bid + &
               Mesh_in % Elements(k) % BoundaryInfo % Constraint)

          Mesh_out % Elements(cnt) % NDOFs = 1
          ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(1))
          Mesh_out % Elements(cnt) % NodeIndexes = &
               Mesh_in % Elements(k) % NodeIndexes+(in_levels+1)*n
          Mesh_out % Elements(cnt) % TYPE => &
               Mesh_in % Elements(k) % TYPE

          Mesh_out % Elements(cnt) % ElementIndex = cnt
          Mesh_out % Elements(cnt) % DGDOFs = 0
          Mesh_out % Elements(cnt) % DGIndexes => NULL()
          Mesh_out % Elements(cnt) % PDefs => NULL()
          Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
          Mesh_out % Elements(cnt) % FaceIndexes => NULL()
          Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
       END DO
    END IF
    
    IF(isParallel) THEN
      j=max_bid
      CALL MPI_ALLREDUCE(j,max_bid,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'First Extruded BC set to: ',max_bid+1
    CALL Info(Caller,Message,Level=8)

    max_body=0
    DO i=1,Mesh_in % NumberOfBulkElements
      max_body = MAX(max_body,Mesh_in % Elements(i) % Bodyid)
    END DO
    IF(isParallel) THEN
      j=max_body
      CALL MPI_ALLREDUCE(j,max_body,1, &
          MPI_INTEGER,MPI_MAX,ELMER_COMM_WORLD,ierr)
    END IF

    WRITE( Message,'(A,I0)') 'Number of new BCs for layers: ',max_body
    CALL Info(Caller,Message,Level=8)


    ! Add start and finish planes except if we have a full rotational symmetry
    IF( .NOT. Rotate2Pi ) THEN

    ! Add bottom boundary:
    ! --------------------
    DO i=1,Mesh_in % NumberOfBulkElements
      cnt=cnt+1

      Mesh_out % Elements(cnt) = Mesh_in % Elements(i)

      l_n=Mesh_in % Elements(i) % TYPE % NumberOfNodes
      Mesh_out % Elements(cnt) % NDOFs = l_n

      ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
      Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
           Mesh_out % Elements(i)
      Mesh_out % Elements(cnt) % BoundaryInfo % Right => NULL()

      bcid = max_bid + Mesh_out % Elements(cnt) % BodyId
      Mesh_out % Elements(cnt) % BoundaryInfo % Constraint = bcid

      Mesh_out % Elements(cnt) % BodyId = 0
      IF( bcid<=CurrentModel % NumberOfBCs) THEN
        j=ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
        IF(Found) Mesh_out % Elements(cnt) % BodyId=j
      END IF

      ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n))
      Mesh_out % Elements(cnt) % NodeIndexes = &
        Mesh_in % Elements(i) % NodeIndexes
      Mesh_out % Elements(cnt) % ElementIndex = cnt
      Mesh_out % Elements(cnt) % TYPE => &
        Mesh_in % Elements(i) % TYPE
      Mesh_out % Elements(cnt) % DGDOFs = 0
      Mesh_out % Elements(cnt) % DGIndexes => NULL()
      Mesh_out % Elements(cnt) % PDefs => NULL()
      Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
      Mesh_out % Elements(cnt) % FaceIndexes => NULL()
      Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
    END DO

    ! Add top boundary:
    ! -----------------
    DO i=1,Mesh_in % NumberOfBulkElements
      cnt=cnt+1

      Mesh_out % Elements(cnt) = Mesh_in % Elements(i)

      l_n=Mesh_in % Elements(i) % TYPE % NumberOfNodes
      Mesh_out % Elements(cnt) % NDOFs = l_n

      ALLOCATE(Mesh_out % Elements(cnt) % BoundaryInfo)
      Mesh_out % Elements(cnt) % BoundaryInfo % Left => &
           Mesh_out % Elements(in_levels*Mesh_in % NumberOfBulkElements+i)
      Mesh_out % Elements(cnt) % BoundaryInfo % Right => NULL()

      bcid = max_bid + Mesh_out % Elements(cnt) % BodyId + max_body
      Mesh_out % Elements(cnt) % BoundaryInfo % Constraint = bcid

      Mesh_out % Elements(cnt) % BodyId = 0
      IF( bcid<=CurrentModel % NumberOfBCs) THEN
        j=ListGetInteger(CurrentModel % BCs(bcid) % Values,'Body Id',Found)
        IF(Found) Mesh_out % Elements(cnt) % BodyId=j
      END IF

      ALLOCATE(Mesh_out % Elements(cnt) % NodeIndexes(l_n))
      Mesh_out % Elements(cnt) % NodeIndexes = &
        Mesh_in % Elements(i) % NodeIndexes+(in_Levels+1)*n
      Mesh_out % Elements(cnt) % ElementIndex = cnt
      Mesh_out % Elements(cnt) % TYPE => &
        Mesh_in % Elements(i) % TYPE
      Mesh_out % Elements(cnt) % DGDOFs = 0
      Mesh_out % Elements(cnt) % DGIndexes => NULL()
      Mesh_out % Elements(cnt) % PDefs => NULL()
      Mesh_out % Elements(cnt) % EdgeIndexes => NULL()
      Mesh_out % Elements(cnt) % FaceIndexes => NULL()
      Mesh_out % Elements(cnt) % BubbleIndexes => NULL()
    END DO

    END IF ! .NOT. Rotate2Pi
    

    Mesh_out % NumberOfBoundaryElements=cnt-Mesh_out % NumberOfBulkElements

    Mesh_out % Name=Mesh_in % Name
    Mesh_out % DiscontMesh = Mesh_in % DiscontMesh
    Mesh_out % MaxElementDOFs  = Mesh_out % MaxElementNodes
    Mesh_out % Stabilize = Mesh_in % Stabilize
    Mesh_out % MeshDim = 3
    CurrentModel % Dimension = 3

    CALL PrepareMesh( CurrentModel, Mesh_out, isParallel )
    
    IF (PRESENT(ExtrudedMeshName)) THEN
       CALL WriteMeshToDisk(Mesh_out, ExtrudedMeshName)
    END IF

    !------------------------------------------------------------------------------
  END FUNCTION MeshExtrude
!------------------------------------------------------------------------------


!------------------------------------------------------------------------------
!>   Create a 1D mesh, may be used in 1D outlet conditions, for example.
!------------------------------------------------------------------------------
  FUNCTION CreateLineMesh( Params ) RESULT( Mesh )
!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params 
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    INTEGER :: i, j, k, n, NoNodes, NoElements, ActiveDirection, Order, BodyId, ne
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshName
    REAL(KIND=dp), ALLOCATABLE :: w(:)
    CHARACTER(*), PARAMETER :: Caller='CreateLineMesh'
    
!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info(Caller,'Creating 1D mesh on-the-fly')

!   Read in the parameters defining a uniform 1D mesh
!--------------------------------------------------------------    
    Order = ListGetInteger( Params,'1D Element Order',Found,minv=1,maxv=2)
    NoElements = ListGetInteger( Params,'1D Number Of Elements',minv=1)
    Length = ListGetConstReal( Params,'1D Mesh Length',Found)
    IF(.NOT. Found) Length = 1.0_dp
    ActiveDirection = ListGetInteger( Params,'1D Active Direction',Found,minv=-3,maxv=3)
    IF(.NOT.Found) ActiveDirection = 1
    BodyId = ListGetInteger( Params,'1D Body Id',Found,minv=1)
    IF(.NOT. Found) BodyId = 1
    MeshName = ListGetString( Params,'1D Mesh Name',Found)
    IF(.NOT. Found) MeshName = '1d_mesh'
    
    Mesh % Name = MeshName
    Mesh % OutputActive = .FALSE.

!   Compute the resulting mesh parameters
!--------------------------------------------------------------
    ne = Order + 1
    NoNodes = NoElements + 1 + NoElements * (Order - 1)    
    MeshVector = 0.0_dp
    MeshVector( ABS( ActiveDirection ) ) = 1.0_dp
    IF( ActiveDirection < 0 ) MeshVector = -MeshVector
    MeshVector = MeshVector * Length
    
!   Define nodal coordinates
!   -------------------------------
    CALL AllocateVector( Mesh % Nodes % x, NoNodes )
    CALL AllocateVector( Mesh % Nodes % y, NoNodes )
    CALL AllocateVector( Mesh % Nodes % z, NoNodes )

    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z

    ALLOCATE( w(0:NoNodes-1) )
    
    CALL UnitSegmentDivision( w, NoNodes-1, Params )
    
    DO i=1, NoNodes
      Coord = MeshVector * w(i-1)

      x(i) = Coord(1)
      y(i) = Coord(2)
      z(i) = Coord(3)
    END DO
    

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    Elmt => GetElementType( 200 + ne )

    DO i=1,NoElements
      Element => Mesh % Elements(i)      
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()     
      Element % ElementIndex = i

      CALL AllocateVector( Element % NodeIndexes, ne )
      Element % Ndofs = ne

      Element % NodeIndexes(1) = (i-1)*Order + 1
      Element % NodeIndexes(2) = i*Order + 1

      DO j=3,ne
        Element % NodeIndexes(j) = (i-1)*Order + j-1
      END DO
      
      Element % BodyId = BodyId
      Element % PartIndex = ParEnv % myPE
    END DO
    
!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = ne
    Mesh % MaxElementDOFs = ne
    Mesh % MeshDim = 1

    WRITE(Message,'(A,I0)') 'Number of elements created: ',NoElements
    CALL Info(Caller,Message)

    WRITE(Message,'(A,I0)') 'Number of nodes created: ',NoNodes
    CALL Info(Caller,Message)
 
    CALL Info(Caller,'All done')

  END FUNCTION CreateLineMesh

  
  !> Creates a regular 2D mesh of 404 elements
  !> The resulting mesh has no boundary elements etc for now
  !> Should only be used for e.g. mesh to mesh interpolation
  !-------------------------------------------------------------
  FUNCTION CreateRectangularMesh(Params) RESULT(Mesh)

!------------------------------------------------------------------------------
    TYPE(ValueList_t), POINTER :: Params
    TYPE(Mesh_t), POINTER :: Mesh
!------------------------------------------------------------------------------
    REAL(KIND=dp), POINTER :: x(:),y(:),z(:)
    REAL(KIND=dp) :: min_x, max_x, min_y, max_y, dx, dy
    INTEGER :: i, j, k, n, counter, nnx, nny, nex, ney, &
         NoNodes, NoElements, col, row
    LOGICAL :: Found
    TYPE(Element_t), POINTER :: Element
    TYPE(ElementType_t),POINTER :: elmt
    REAL(KIND=dp) :: MeshVector(3), Length, Coord(3)
    CHARACTER(LEN=MAX_NAME_LEN) :: MeshName
    CHARACTER(*), PARAMETER :: Caller='CreateReactangularMesh'

!------------------------------------------------------------------------------
    Mesh => NULL()
    IF ( .NOT. ASSOCIATED( Params ) ) RETURN
    Mesh => AllocateMesh()

    CALL Info(Caller,'Creating 2D mesh on-the-fly')

    !Get parameters from valuelist
    min_x = ListGetConstReal(Params, "Grid Mesh Min X",UnfoundFatal=.TRUE.)
    max_x = ListGetConstReal(Params, "Grid Mesh Max X",UnfoundFatal=.TRUE.)
    min_y = ListGetConstReal(Params, "Grid Mesh Min Y",UnfoundFatal=.TRUE.)
    max_y = ListGetConstReal(Params, "Grid Mesh Max Y",UnfoundFatal=.TRUE.)
    dx    = ListGetConstReal(Params, "Grid Mesh dx",UnfoundFatal=.TRUE.)
    dy    = ListGetConstReal(Params, "Grid Mesh dy",Found)
    IF(.NOT. Found) dy = dx

    IF(max_x <= min_x .OR. max_y <= min_y .OR. dx <= 0.0_dp .OR. dy <= 0.0_dp) &
         CALL Fatal(Caller, "Bad Grid Mesh parameters!")

    !number of nodes in x and y direction (and total)
    nnx = FLOOR((max_x - min_x) / dx) + 1
    nny = FLOOR((max_y - min_y) / dy) + 1
    NoNodes = nnx * nny

    !number of elements in x and y direction (and total)
    nex = nnx - 1
    ney = nny - 1
    NoElements = nex * ney


!   Define nodal coordinates
!   -------------------------------
    CALL AllocateVector( Mesh % Nodes % x, NoNodes )
    CALL AllocateVector( Mesh % Nodes % y, NoNodes )
    CALL AllocateVector( Mesh % Nodes % z, NoNodes )
    x => Mesh % Nodes % x
    y => Mesh % Nodes % y
    z => Mesh % Nodes % z

    z = 0.0_dp !2D

    !Define node positions
    counter = 0
    DO i=1,nnx
      DO j=1,nny
        counter = counter + 1
        x(counter) = min_x + (i-1)*dx
        y(counter) = min_y + (j-1)*dy
      END DO
    END DO

!   Define elements
!   -------------------------------
    CALL AllocateVector( Mesh % Elements, NoElements )

    Elmt => GetElementType( 404 )

    DO i=1,NoElements
      Element => Mesh % Elements(i)
      Element % TYPE => Elmt
      Element % EdgeIndexes => NULL()
      Element % FaceIndexes => NULL()
      Element % ElementIndex = i
      CALL AllocateVector( Element % NodeIndexes, 4 )
      Element % Ndofs = 4

      col = MOD(i-1,ney)
      row = (i-1)/ney

      !THIS HERE NEEDS FIXED!!!!!
      Element % NodeIndexes(1) = (row * nny) + col + 1
      Element % NodeIndexes(2) = (row * nny) + col + 2
      Element % NodeIndexes(4) = ((row+1) * nny) + col + 1
      Element % NodeIndexes(3) = ((row+1) * nny) + col + 2

      Element % BodyId = 1
      Element % PartIndex = ParEnv % myPE
    END DO

!   Update new mesh node count:
!   ---------------------------

    Mesh % NumberOfNodes = NoNodes
    Mesh % Nodes % NumberOfNodes = NoNodes
    Mesh % NumberOfBulkElements = NoElements
    Mesh % MaxElementNodes = 4
    Mesh % MaxElementDOFs = 4
    Mesh % MeshDim = 2

  END FUNCTION CreateRectangularMesh

END MODULE MeshGeneration
  
