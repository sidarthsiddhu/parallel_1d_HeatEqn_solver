! ----------- ECE 563 Programming Parallel Machines --------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver ----------- !
! ----------- Author: Sidarth Narayanan ----------- !

Program Main

    Use MPI
    Implicit None

    Integer(kind=8) :: npts,ntime,it,i,npts_local
    Integer(kind=8) :: remeinder,start_indx,end_indx
    Real(kind=8),Allocatable,Dimension(:) :: T_old,T_new,T_exact
    Real(kind=8),Allocatable,Dimension(:) :: x,t
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r
    Real(kind=8) :: T_left,T_right,T_initial
    Real(kind=8) :: ghost_left,ghost_right
    Integer :: request_left,request_right,ierr,myid,nprocs,debug_version
    Integer,Dimension(MPI_STATUS_SIZE) :: status

! Initialize MPI and call rank and size functions
    Call MPI_INIT(ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myid  , ierr)

! Set the problem parameters
    debug_version = 0
    length = 1.0
    Total_time = 1.0
    npts = 10
    ntime = 1
    alpha = 0.0001
    T_left = 0.0 
    T_right = 100.0
    T_initial = 25.0
    dx = length/real(npts)
    dt = Total_time/real(ntime)
    r = (alpha*dt)/(dx*dx)
	
	If(myid .EQ. 0) then
		write(6,'(a75)') " --- Parallel 1-D Heat Equation Solver [dT/dt = alpha*(d2T/dx2)] --- "
		write(6,'(a75)') " --------------------- Problem Parameters -------------------------- "
		write(6,'(a34,F12.4,a)') " Length: ",length," meters"
		write(6,'(a34,F12.4,a)') " Total simulation time: ",Total_time," seconds"
		write(6,'(a34,ES12.5,a)') " dt: ",dt," seconds"
		write(6,'(a34,ES12.5,a)') " dx: ",dx," meters"
		write(6,'(a34,I0,a)') " Number of points (spatial): ",npts," Point(s)"
		write(6,'(a34,I0,a)') " Number of time steps : ",ntime," Step(s)"
		write(6,'(a34,ES12.5,a)') " Thermal Diffusivity: ",alpha," m^2/s"
		write(6,'(a34,ES12.5)') " R(alpha*dt/dx^2) : ",r
		write(6,'(a34,I0,a)') " Number of Processors : ",nprocs," Processor(s)"
		write(6,*) ""
		Call flush(6)
	Endif
	Call MPI_Barrier(MPI_COMM_WORLD,ierr)

! Split the domain among the processors in a load balanced method 
    npts_local = npts/nprocs
    remeinder = mod(npts,nprocs)
    If(myid .LT. remeinder) then
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
        T_new(i) = T_initial + i
        x(i) = (i-1)*dx
    Enddo
	t(1) = 0.0 ! Start with time 0.0 seconds
! Set boundary condition
    If(myid .EQ. 0) T_new(start_indx) = T_left
    If(myid .EQ. (nprocs-1)) T_new(end_indx) = T_right
	
    Do it = 2,ntime
    
        Call MPI_Barrier(MPI_COMM_WORLD,ierr)
		t(it) = t(it-1) + dt
		Do i=start_indx,end_indx
			T_old(i) = T_new(i)
		Enddo
		
		If(nprocs .GT. 1) then
			If(myid .EQ. 0) then
				Call MPI_Sendrecv(T_new(end_indx),1,MPI_DOUBLE_PRECISION,myid+1,myid,ghost_right,1,MPI_DOUBLE_PRECISION,myid+1,myid+1,MPI_COMM_WORLD,status,ierr)
			Elseif(myid .EQ. (nprocs-1)) then
				Call MPI_Sendrecv(T_new(start_indx),1,MPI_DOUBLE_PRECISION,myid-1,myid,ghost_left,1,MPI_DOUBLE_PRECISION,myid-1,myid-1,MPI_COMM_WORLD,status,ierr)
			Else
				Call MPI_Sendrecv(T_new(end_indx),1,MPI_DOUBLE_PRECISION,myid+1,myid,ghost_right,1,MPI_DOUBLE_PRECISION,myid+1,myid+1,MPI_COMM_WORLD,status,ierr)
				Call MPI_Sendrecv(T_new(start_indx),1,MPI_DOUBLE_PRECISION,myid-1,myid,ghost_left,1,MPI_DOUBLE_PRECISION,myid-1,myid-1,MPI_COMM_WORLD,status,ierr)
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
		Else
			Do i=start_indx+1,end_indx-1
				T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
			Enddo
		Endif
			
    Enddo

    Call MPI_Barrier(MPI_COMM_WORLD,ierr)

    If(debug_version .EQ. 1) then
        If(myid .EQ. 0) write(6,'(a75)') " --------------------------- Final Answer -------------------------- "
            Do i=start_indx,end_indx
                write(6,'(a34,I0,a,F12.4)') "T(",i,") = ",T_new(i)
            Enddo
            Call flush(6)
    Endif

    Deallocate(T_old)
    Deallocate(T_new)
    Deallocate(T_exact)
    Deallocate(x)
    Deallocate(t)
	
    Call MPI_FINALIZE(ierr)

End Program Main
