! ----------- ECE 563 Programming Parallel Machines --------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver ----------- !
! ----------- Author: Sidarth Narayanan ----------- !

Program Main

    Use MPI
    Implicit None

    Integer(kind=8) :: npts,ntime,it,i,npts_local
    Integer(kind=8) :: start_indx,end_indx
    Real(kind=8),Allocatable,Dimension(:) :: T_old,T_new
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r
    Real(kind=8) :: T_left,T_right,T_initial,T_exact,error,error_max
    Real(kind=8) :: ghost_left,ghost_right,time_start,time_end,time_elapsed
    Integer :: ierr,myid,nprocs,debug_version,remeinder,print_freq,print_result
    Integer,Dimension(MPI_STATUS_SIZE) :: status_right,status_left

! Initialize MPI and call rank and size functions
    Call MPI_INIT(ierr)
    Call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    Call MPI_COMM_RANK(MPI_COMM_WORLD, myid  , ierr)

! Set code debug parameters
    debug_version = 1
	print_result = 0
	print_freq = ntime/10
	
! Set the problem parameters	
	ntime = 10
    npts = 10**9
    dx = 1e-03
	dt = 1e-03
	r = 0.5
	alpha = r*(dx**2)/dt
    T_left = 0.0 
    T_right = 0.0
    T_initial = 50.0
	T_exact = 0.0
	error = 0.0
	
	length = real(npts)*dx
	Total_time = real(ntime)*dt

! Print out problem parameters
	If(myid .EQ. 0) then
		write(6,'(a75)') " --- Parallel 1-D Heat Equation Solver [dT/dt = alpha*(d2T/dx2)] --- "
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
	
! Check npts_local and indices
	If(debug_version .EQ. 1) then
		write(6,'(a,I0,a,I0,a,I0,a,I0)') "Id: ",myid," || Start: ",start_indx," || End: ",end_indx," || npts_local: ",npts_local
	Endif

! Allocate the variables
    Allocate(T_old(start_indx:end_indx))
    Allocate(T_new(start_indx:end_indx))

! Set initial condition
    !Do i=start_indx,end_indx
    !    T_new(i) = T_initial
    !Enddo
    T_new(:) = T_initial
	
! Set boundary condition
    If(myid .EQ. 0) T_new(start_indx) = T_left
    If(myid .EQ. (nprocs-1)) T_new(end_indx) = T_right
	
! Calculate initial RMS error from exact solution
	Do i = start_indx,end_indx
		error = error + ABS(T_new(i) - T_exact)
	Enddo
	error = (error/npts_local)
	Call MPI_Reduce(error,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
	If(myid .EQ. 0) write(6,'(a,ES12.5)') "The initial average Error is: ",error_max
	
! Reset Error value for final calculation
	error = 0.0

	If(myid .EQ. 0) write(6,'(a)') "Starting Calculation ... "

! Start the time counter
	Call MPI_Barrier(MPI_COMM_WORLD,ierr)
	time_start = MPI_Wtime()

! Start the time iteration
    Do it = 2,ntime
    
	! Set T_old = T_new
		!Do i=start_indx,end_indx
			T_old(i) = T_new(i)
		!Enddo
		T_old = T_new
		
	! Pass the details of the ghostg points to the neighbours		
		If(nprocs .GT. 1) then
		
		! Send the left most element
			If(myid .GT. 0)	Call MPI_Send(T_new(start_indx),1,MPI_DOUBLE_PRECISION,myid-1,1,MPI_COMM_WORLD,ierr)
		! Recieve the right ghost element
			If(myid .LT. nprocs-1) Call MPI_Recv(ghost_right,1,MPI_DOUBLE_PRECISION,myid+1,1,MPI_COMM_WORLD,status_right,ierr)
			
		! Send the right most element	
			If(myid .LT. nprocs-1)	Call MPI_Send(T_new(end_indx),1,MPI_DOUBLE_PRECISION,myid+1,1,MPI_COMM_WORLD,ierr)
		! Recieve the left ghost element
			If(myid .GT. 0) Call MPI_Recv(ghost_left,1,MPI_DOUBLE_PRECISION,myid-1,1,MPI_COMM_WORLD,status_left,ierr)
		
		! Do the actual calculation point by point
			If(myid .EQ. 0) then
				Do i=start_indx+1,end_indx
					If(i .EQ. end_indx) then
					! End point uses right ghost node	
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*ghost_right
					Else
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
					Endif
				Enddo
			Elseif(myid .EQ. (nprocs-1)) then
				Do i=start_indx,end_indx-1
					If(i .EQ. start_indx) then
					! Start point uses left ghost node
						T_new(i) = r*ghost_left + (1-(2*r))*T_old(i) + r*T_old(i+1)
					Else
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
					Endif
				Enddo
			Else
				Do i=start_indx,end_indx
					If(I .EQ. start_indx) then
					! Start point uses left ghost node
						T_new(i) = r*ghost_left + (1-(2*r))*T_old(i) + r*T_old(i+1)
					Elseif(i .EQ. end_indx) then
					! End point uses right ghost node
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*ghost_right
					Else
						T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
					Endif
				Enddo
			Endif
		Else
			Do i=start_indx+1,end_indx-1
				T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
			Enddo
		Endif
			
    Enddo
	
! Stop time count and calculate elapsed time
	Call MPI_Barrier(MPI_COMM_WORLD,ierr)
	time_end = MPI_Wtime()
	time_elapsed = time_end - time_start
	
	If(myid .EQ. 0) write(6,'(a)') "Calculation Complete"
	
! Calculate error from exact solution to check the results
	Do i = start_indx,end_indx
		error = error + ABS(T_new(i) - T_exact)
	Enddo
	error = (error/npts_local)
	Call MPI_Reduce(error,error_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

! Print out elapsed time
	If(myid .EQ. 0) then
		write(6,'(a,ES12.5)') "The average Error is: ",error_max
		write(6,'(a,ES12.5,a)') "The Elapsed time is: ",time_elapsed," seconds"
	Endif
	
    If((debug_version .EQ. 1) .AND. (print_result .EQ. 1)) then
        If(myid .EQ. 0) write(6,'(a75)') " --------------------------- Final Answer -------------------------- "
            Do i=start_indx,end_indx
                write(6,'(a34,I0,a,ES12.5)') "T(",i,") = ",T_new(i)
            Enddo
    Endif

! Deallocate Variables	
    Deallocate(T_old)
    Deallocate(T_new)
   
! Finalize MPI
	Call MPI_FINALIZE(ierr)

End Program Main
