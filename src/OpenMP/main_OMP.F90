! ----------------- ECE 563 Programming Parallel Machines ------------------!
! ----------- Final Project : Parallel 1-D Heat Equation Solver ----------- !
! ----------------------------- Version : OpenMP ---------------------------!
! ------------------------ Author: Sidarth Narayanan ---------------------- !

Program Main
	
	USE OMP_LIB
    Implicit None

    Integer(kind=8) :: npts,ntime,it,i
    Real,Allocatable,Dimension(:) :: T_old,T_new,T_exact
    Real,Allocatable,Dimension(:) :: x,t
    Real(kind=8) :: Total_time,length,dx,dt,alpha,r
    Real(kind=8) :: T_left,T_right,T_initial
    Integer :: request_left,request_right,ierr,myid,nprocs,debug_version
	Real(kind=8) :: time_taken

! Set the problem parameters
    debug_version = 0
    length = 1.0
    Total_time = 1.0
    npts = 1000000
    ntime = 100000
    alpha = 0.001
    T_left = 0.0 
    T_right = 100.0
    T_initial = 25.0
    dx = length/real(npts)
    dt = Total_time/real(ntime)
    r = (alpha*dt)/(dx*dx)

! Allocate the variables
    Allocate(T_old(npts))
    Allocate(T_new(npts))
    Allocate(T_exact(npts))
    Allocate(x(npts))
    Allocate(t(ntime))

! Set initial condition
    Do i=1,npts
        T_new(i) = T_initial + i
        x(i) = (i-1)*dx
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
		
		!$OMP PARALLEL DO PRIVATE(i)
		Do i=2,npts-1
			T_old(i) = T_new(i)
		Enddo
		!$OMP END PARALLEL DO
		
		!$OMP PARALLEL DO PRIVATE(i)
        Do i=2,npts-1
            T_new(i) = r*T_old(i-1) + (1-(2*r))*T_old(i) + r*T_old(i+1)
        Enddo
		!$OMP END PARALLEL DO
		
    Enddo
	time_taken = time_taken + omp_get_wtime()
	
	write(6,'(a,F12.5,a)') "Time taken to solve: ",time_taken," seconds"
	
    If(debug_version .EQ. 1) then
        write(6,*) "-------------- Final Answer ---------------"
        Do i=1,npts
            write(6,*) "T(",i,") = ",T_new(i)
        Enddo
    Endif

    Deallocate(T_old)
    Deallocate(T_new)
    Deallocate(T_exact)
    Deallocate(x)
    Deallocate(t)

End Program Main
