! ----------- ECE 563 Programming Parallel Machines --------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver (MPI)----------- !
! ----------- Author: Sidarth Narayanan ----------- !

Program Main

    Use MPI
    Implicit None

    Integer(kind=8) :: npts,ntime,it,i,npts_local
    Real(kind=8),Allocatable,Dimension(:) :: T_old,T_new
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r,k
    Real(kind=8) :: T_left,T_right,T_initial,iter_time_start,iter_time_end,iter_time
    Real(kind=8) :: ghost_left,ghost_right,time_start,time_end,time_elapsed
    Integer :: ierr,myid,nprocs,remeinder,debug_version
    Integer,Dimension(MPI_STATUS_SIZE) :: status_right,status_left

! Initialize MPI and call rank and size functions
    Call MPI_INIT(ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myid  , ierr)

! Set the problem parameters
	debug_version = 0	
	ntime = 1000
    npts = 10**7
    dx = 1e-03
	dt = 1e-03
	r = 0.5
	alpha = r*(dx**2)/dt
	k = 1e-03
    T_left = 0.0 
    T_right = 0.0
    T_initial = 50.0
	
	length = real(npts)*dx
	Total_time = real(ntime)*dt

! Print out problem parameters
	If(myid .EQ. 0) then
		write(6,'(a)') " "
		write(6,'(a)') " "
		write(6,'(a75)') " --- Parallel 1-D Heat Equation Solver [dT/dt = alpha*(d2T/dx2)] [MPI] --- "
		write(6,'(a75)') " --------------------- Problem Parameters -------------------------- "
		write(6,'(a34,ES12.5,a)') " Length: ",length," meters"
		write(6,'(a34,ES12.5,a)') " Total simulation time: ",Total_time," seconds"
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
	
	!Call MPI_Barrier(MPI_COMM_WORLD,ierr)

! Split the domain among the processors in a load balanced method 
    npts_local = npts/nprocs
    remeinder = mod(npts,nprocs)
    If(myid .LT. remeinder) npts_local = npts_local + 1
        
! Print start and end indx
	If(debug_version .EQ. 1) then
		write(6,'(a,I0,a,I0)') "Myid: ",myid," || npts_local: ",npts_local
		Call Flush(6)
	Endif
	
! Allocate the variables
    Allocate(T_old(npts_local))
    Allocate(T_new(npts_local))

! Set initial condition
    Do i=1,npts_local
        T_new(i) = T_initial
    Enddo
	
! Set boundary condition
    If(myid .EQ. 0) T_new(1) = T_left
    If(myid .EQ. (nprocs-1)) T_new(npts_local) = T_right

! Print starting calculation
	If(myid .EQ. 0) then
		write(6,'(a)') "Starting Calculation"
		write(6,'(a)') " "
		write(6,'(a)') " "
	Endif

! Start the time counter
	Call MPI_Barrier(MPI_COMM_WORLD,ierr)
	time_start = MPI_Wtime()

! Start the time iteration
    Do it = 1,ntime
    
	! Print Iteration number
		if((myid .EQ. 0) .AND. (mod(it,(ntime/10)) .EQ. 0)) write(6,'(a,I0,a)') "Runnnig Iteration ",it," ..."
		
	! Time one of the iteration
		If(it .EQ. ntime/2) iter_time_start = MPI_Wtime()
	
	! Set T_old = T_new
		T_old = T_new
		
	! Pass the details of the ghostg points to the neighbours		
		If(nprocs .GT. 1) then
		
		! Send the left most element
			If(myid .GT. 0)	Call MPI_Send(T_new(1),1,MPI_DOUBLE_PRECISION,myid-1,1,MPI_COMM_WORLD,ierr)
		! Recieve the right ghost element
			If(myid .LT. nprocs-1) Call MPI_Recv(ghost_right,1,MPI_DOUBLE_PRECISION,myid+1,1,MPI_COMM_WORLD,status_right,ierr)
			
		! Send the right most element	
			If(myid .LT. nprocs-1)	Call MPI_Send(T_new(npts_local),1,MPI_DOUBLE_PRECISION,myid+1,1,MPI_COMM_WORLD,ierr)
		! Recieve the left ghost element
			If(myid .GT. 0) Call MPI_Recv(ghost_left,1,MPI_DOUBLE_PRECISION,myid-1,1,MPI_COMM_WORLD,status_left,ierr)
		
		! Do the actual calculation point by point
			If(myid .EQ. 0) then
			! Left end "processor"
				Do i=2,npts_local
					If(i .EQ. npts_local) then
					! End point uses right ghost node	
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*ghost_right
						T_new(i) = T_new(i) + (k/t_new(i))
					Else
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
						T_new(i) = T_new(i) + (k/t_new(i))
					Endif
				Enddo
			Elseif(myid .EQ. (nprocs-1)) then
			! Right end "processor"
				Do i=1,npts_local-1
					If(i .EQ. 1) then
					! Start point uses left ghost node
						T_new(i) = r*ghost_left + (1-(2*r))*T_old(i) + r*T_old(i+1)
						T_new(i) = T_new(i) + (k/t_new(i))
					Else
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
						T_new(i) = T_new(i) + (k/t_new(i))
					Endif
				Enddo
			Else
			! Interior "processors"
				Do i=1,npts_local
					If(i .EQ. 1) then
					! Start point uses left ghost node
						T_new(i) = r*ghost_left + (1-(2*r))*T_old(i) + r*T_old(i+1)
						T_new(i) = T_new(i) + (k/t_new(i))
					Elseif(i .EQ. npts_local) then
					! End point uses right ghost node
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*ghost_right
						T_new(i) = T_new(i) + (k/t_new(i))
					Else
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
						T_new(i) = T_new(i) + (k/t_new(i))
					Endif
				Enddo
			Endif
		Else
		! Basic Calculation for just 1 element without message passing
			Do i=2,npts_local-1
				T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
				T_new(i) = T_new(i) + (k/t_new(i))
			Enddo
		Endif
		
		If(it .EQ. ntime/2) iter_time_end = MPI_Wtime()
			
    Enddo
	
! Stop time count and calculate elapsed time
	Call MPI_Barrier(MPI_COMM_WORLD,ierr)
	time_end = MPI_Wtime()
	time_elapsed = time_end - time_start
	iter_time = iter_time_end - iter_time_start

! Print out elapsed time	
	If(myid .EQ. 0) then
		write(6,'(a)') "Calculation Complete"
		write(6,'(a)') " "
		write(6,'(a,ES12.5,a)') " The time per iteration is: ",iter_time," seconds"
		write(6,'(a,ES12.5,a)') "The total elapsed time is: ",time_elapsed," seconds"
	Endif
	
! Print End Result
	If(debug_version .EQ. 1) then
		if(myid .EQ. 0) write(6,'(a75)') " ---------- Final Answer --------------- "
		Do i=1,npts_local
			write(6,'(a,I0,a,I0,a,ES12.5)') "Myid: ",myid," || T(",i,") = ",T_new(i)
		Enddo
	Endif
! Deallocate Variables	
    Deallocate(T_old)
    Deallocate(T_new)
   
! Finalize MPI
	Call MPI_FINALIZE(ierr)

End Program Main
