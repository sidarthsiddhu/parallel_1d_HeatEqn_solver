! ----------- ECE 563 Programming Parallel Machines --------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver (Sequential Version)----------- !
! ----------- Author: Sidarth Narayanan ----------- !

Program Main

    Implicit None

    Integer(kind=8) :: npts,ntime,it,i
    Real,Allocatable,Dimension(:) :: T_old,T_new
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r,k
    Real(kind=8) :: T_left,T_right,T_initial
    Real(kind=8) :: t1,t2,t3,t4,time_total,time_iter
    Integer :: debug_version

! Initialize MPI and call rank and size functions

! Set the problem parameters
    debug_version = 0
    ntime = 1000
    npts = 1000000
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


! Allocate the variables
    Allocate(T_old(npts))
    Allocate(T_new(npts))
    
! Print out problem parameters
	write(6,'(a)') " "
	write(6,'(a)') " "
	write(6,'(a75)') " -- Parallel 1-D Heat Equation Solver [dT/dt = alpha*(d2T/dx2)] [Sequential] -- "
	write(6,'(a75)') " --------------------- Problem Parameters -------------------------- "
	write(6,'(a34,ES12.5,a)') " Length: ",length," meters"
	write(6,'(a34,ES12.5,a)') " Total simulation time: ",Total_time," seconds"
	write(6,'(a34,ES12.5,a)') " dt: ",dt," seconds"
	write(6,'(a34,ES12.5,a)') " dx: ",dx," meters"
	write(6,'(a34,I0,a)') " Number of points (spatial): ",npts," Point(s)"
	write(6,'(a34,I0,a)') " Number of time steps : ",ntime," Step(s)"
	write(6,'(a34,ES12.5,a)') " Thermal Diffusivity: ",alpha," m^2/s"
	write(6,'(a34,ES12.5)') " R(alpha*dt/dx^2) : ",r
	write(6,*) ""

! Set initial condition
	Do i=1,npts
    		T_new(:) = T_initial
	Enddo

! Set boundary condition
    T_new(1) = T_left
    T_new(npts) = T_right

! Start timing 
	Call CPU_TIME(t1)
! Start Main time iteration
    Do it = 1,ntime
    	
    	If(it .EQ. ntime/2) Call CPU_TIME(t2)
    	if(mod(it,(ntime/10)) .EQ. 0) write(6,'(a,I0,a)') "Runnnig Iteration ",it," ..."
        
        T_old = T_new
	
	Do i=2,npts-1
        	T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
        	T_new(i) = T_new(i) + (k/T_new(i))
	Enddo
	
	If(it .EQ. ntime/2) Call CPU_TIME(t3)
	
    Enddo
    
    Call CPU_TIME(t4)
    time_total = t4 - t1
    time_iter = t3 - t2

! Print out time
	Write(6,'(a,ES12.5,a)') "Total Time taken: ",time_total," seconds"
	Write(6,'(a,ES12.5,a)') "Time taken per iteration: ",time_iter," seconds"
	
    If(debug_version .EQ. 1) then
	write(6,*) "-------------- Final Answer ---------------"
        Do i=1,npts
        	write(6,*) "T(",i,") = ",T_new(i)
        Enddo
    Endif

    Deallocate(T_old)
    Deallocate(T_new)

End Program Main
