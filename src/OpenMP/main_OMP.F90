! ----------- ECE 563 Programming Parallel Machines --------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver ----------- !
! ----------- Author: Sidarth Narayanan ----------- !

Program Main

    Use MPI
    Implicit None

    Integer(kind=8) :: npts,ntime,it,i,npts_local
    Integer(kind=8) :: remeinder,start_indx,end_indx
    Real,Allocatable,Dimension(:) :: T_old,T_new,T_exact
    Real,Allocatable,Dimension(:) :: x,t
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r
    Real(kind=8) :: T_left,T_right,T_initial
    Real(kind=8) :: ghost_left,ghost_right
    Integer :: request_left,request_right,ierr,myid,nprocs
    Integer,Dimension(MPI_STATUS_SIZE) :: status_left,status_right

! Initialize MPI and call rank and size functions
    Call MPI_INIT(ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myid  , ierr)

! Set the problem parameters
    length = 1.0
    Total_time = 10.0
    npts = 10
    ntime = 1
    alpha = 1.0
    T_left = 0.0 
    T_right = 100.0
    T_initial = 25.0
    dx = length/real(npts)
    dt = Total_time/real(ntime)
    r = (alpha*dt)/(dx*dx)

! Split the domain among the processors in a load balanced method 
    npts_local = npts/nprocs
    remeinder = mod(npts,nprocs)
    If(myid .lt. remeinder) then
        npts_local = npts_local + 1
        start_indx = (myid*npts_local) + 1
        end_indx = start_indx + npts_local - 1
    Else
        start_indx = (myid * npts_local) + remeinder + 1
        end_indx = start_indx + npts_local - 1
    Endif

! Allocate the variables
    Allocate(T_old(npts_local))
    Allocate(T_new(npts_local))
    Allocate(T_exact(npts_local))
    Allocate(x(npts_local))
    Allocate(t(ntime))

! Set initial condition
    Do i=start_indx,end_indx
        T_new(i) = T_initial
        x(i) = (i-1)*dx
    Enddo

! Set boundary condition
    If(myid .EQ. 0) T_new(start_indx) = T_left
    If(myid .EQ. (nprocs-1)) T_new(end_indx) = T_right

    Do it = 1,ntime
    
        Call MPI_Barrier(MPI_COMM_WORLD,ierr)
        T_old = T_new
        If(myid .EQ. 0) then
            Call MPI_Isend(T_new(end_indx),1,MPI_DOUBLE,myid+1,myid,MPI_COMM_WORLD,request_right,ierr)
            Call MPI_Irecv(ghost_right,1,MPI_DOUBLE,myid+1,myid+1,MPI_COMM_WORLD,status_right,ierr)
            Call MPI_Wait(request_right,status_right,ierr)
        Elseif(myid .EQ. (nprocs-1)) then
            Call MPI_Isend(T_new(end_indx),1,MPI_DOUBLE,myid-1,myid,MPI_COMM_WORLD,request_left,ierr)
            Call MPI_Irecv(ghost_left,1,MPI_DOUBLE,myid-1,myid-1,MPI_COMM_WORLD,status_left,ierr)
            Call MPI_Wait(request_left,status_left,ierr)
        Else
            Call MPI_Isend(T_new(start_indx),1,MPI_DOUBLE,myid-1,myid,MPI_COMM_WORLD,request_right,ierr)
            Call MPI_Isend(T_new(end_indx),1,MPI_DOUBLE,myid+1,myid,MPI_COMM_WORLD,request_left,ierr)
            Call MPI_Irecv(ghost_right,1,MPI_DOUBLE,myid+1,myid+1,MPI_COMM_WORLD,status_right,ierr)
            Call MPI_Irecv(ghost_left,1,MPI_DOUBLE,myid-1,myid-1,MPI_COMM_WORLD,status_left,ierr)
            Call MPI_Wait(request_right,status_right,ierr)
            Call MPI_Wait(request_left,status_left,ierr)
        Endif
    
        If(myid .EQ. 0) then
            Do i=start_indx+1,end_indx
                If(i .EQ. end_indx) then
                    T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*ghost_right
                Else
                    T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
                Endif
            Enddo
        Elseif(myid .EQ. (nprocs-1)) then
            Do i=start_indx,end_indx-1
                If(i .EQ. start_indx) then
                    T_new(i) = r*ghost_left + (1-(2*r))*T_old(i) + r*T_old(i+1)
                Else
                    T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
                Endif
            Enddo
        Else
            Do i=start_indx,end_indx
                T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
            Enddo
        Endif
    Enddo

    Deallocate(T_old)
    Deallocate(T_new)
    Deallocate(T_exact)
    Deallocate(x)
    Deallocate(t)

    Call MPI_FINALIZE(ierr)

End Program Main
