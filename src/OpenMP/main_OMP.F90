! ----------------- ECE 563 Programming Parallel Machines ------------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver ----------- !
! ----------------------------- Version : OpenMP ---------------------------!
! ------------------------ Author: Sidarth Narayanan ---------------------- !

Program Main
	
	USE OMP_LIB
    Implicit None

    Integer(kind=8) :: npts,ntime,it,i
    Real,Allocatable,Dimension(:) :: T_old,T_new
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r,k
    Real(kind=8) :: T_left,T_right,T_initial
    Integer :: nprocs,debug_version
	Real(kind=8) :: time_taken,time_taken_per_iteration

! Set the problem parameters
    debug_version = 0
    ntime = 1000
    npts = 10**8
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
	
	!$OMP PARALLEL
		nprocs = OMP_GET_NUM_THREADS()
	!$OMP END PARALLEL
	
! Print out problem parameters
	write(6,'(a)') " "
	write(6,'(a)') " "
	write(6,'(a75)') " -- Parallel 1-D Heat Equation Solver [dT/dt = alpha*(d2T/dx2)] [OpenMP] -- "
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
	
! Allocate the variables
    Allocate(T_old(npts))
    Allocate(T_new(npts))

! Set initial condition
    Do i=1,npts
        T_new(i) = T_initial + i
    Enddo

! Set boundary condition
	T_new(1) = T_left
    T_new(npts) = T_right
	
	If(debug_version .EQ. 1) then
        write(6,*) "-------------- Initial Condition ---------------"
        Do i=1,npts
            write(6,*) "T(",i,") = ",T_new(i)
        Enddo
    Endif
    
	time_taken = -1.0*omp_get_wtime()
    Do it = 1,ntime
	
	if(it .EQ. ntime/2) time_taken_per_iteration = -1.0*omp_get_wtime()
	
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(DYNAMIC,10000)
	Do i=2,npts-1
		T_old(i) = T_new(i)
	Enddo
	!$OMP END PARALLEL DO
		
	!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(DYNAMIC,10000)
        Do i=2,npts-1
            T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
            T_new(i) = T_new(i) + (dt*k/T_new(i))
        Enddo
	
	!$OMP END PARALLEL DO
	
	if(it .EQ. ntime/2) time_taken_per_iteration = time_taken_per_iteration + omp_get_wtime()	
    Enddo
	time_taken = time_taken + omp_get_wtime()
	
	write(6,'(a,F12.5,a)') "Time taken to solve 1 iteration is: ",time_taken_per_iteration," seconds"
	write(6,'(a,F12.5,a)') "Total time taken to solve: ",time_taken," seconds"
	
    If(debug_version .EQ. 1) then
        write(6,*) "-------------- Final Answer ---------------"
        Do i=1,npts
            write(6,*) "T(",i,") = ",T_new(i)
        Enddo
    Endif

    Deallocate(T_old)
    Deallocate(T_new)

End Program Main
