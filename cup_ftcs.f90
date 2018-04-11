!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------------!
!-------------Module containing central variables for the program-------------!
!-----------------------------------------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module declarations

	implicit none
	private
	double precision,public, allocatable, dimension(:) :: x,x_u,y,y_v,er
	double precision,public :: xl,yl
	integer,public :: ipref,jpref
	double precision,public,allocatable,dimension(:,:) :: T
	double precision,public, parameter :: pi=3.1415927
	double precision,public :: radius,perim,h,Tpar
	double precision,public :: cp,rho,mu,k
	double precision,public :: dx,dy
	real ,public :: relax(6)
	integer,public :: iter,last,npi,npj,nsteps
	double precision,public :: cup_width,k_cup,T_amb,g_const,alpha,tend,dt
	
end module declarations

module sub
	contains
	
	SUBROUTINE tick(t)
        INTEGER, INTENT(OUT) :: t

        CALL system_clock(t)
    END SUBROUTINE tick
    
    REAL FUNCTION tock(t)
        INTEGER, INTENT(IN) :: t
        INTEGER :: now, clock_rate

        call system_clock(now,clock_rate)

        tock = real(now - t)/real(clock_rate)
    END FUNCTION tock

	subroutine init()
		
		!----------------!
		! Initialization !
		!----------------!
		use declarations
		implicit none
		
		!-----------------------------------------------------------------!
		! Set number of grid points in each direction and allocate arrays !
		!-----------------------------------------------------------------!
		npi=52
	  	npj=102
	  	last=400	
		allocate(T(npi,npj))
		allocate(er(last))
	  
		!--------------------------------------!
		! Set timestep and simulation end time !
		!--------------------------------------!
		tend=3600
		dt=0.15
		nsteps=ceiling(tend/dt)

		!---------------------------------!
		! Set reference for zero pressure !
		!---------------------------------!
		ipref=npi/2
		jpref=2
	
		!----------------------------------!
		! Initialize simulation parameters !
		!----------------------------------!
		T_amb=273.16+20
		T(1:npi,1:npj)=273.16+83.00
	
		!---------------------!
		! Set constant values !
		!---------------------!
	 	g_const=9.81	! Standard gravity
	 	!alpha=0.000069	! Thermal expansion coeff.
	 	rho=1000		! Density
	 	!mu=0.001   		! Dynamic viscosity
	 	k=0.64  		! Thermal conductivety 
	 	cp=4181   		! Specific heat
	  
	  
	  	!-----------------------------------!
		! Set cup/wall heat loss parameters !
		!-----------------------------------!
	    cup_width=0.01
	    k_cup=77.1	
	
	end subroutine init


	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-----------------Defines grid and computational domain-------------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine grid()
		
		!----------------!
		! Initialization !
		!----------------!
		use declarations
	  	implicit none
	 	integer :: i,j
		
		!-------------------------!
		! Allocate length vectors !
		!-------------------------!
	  	allocate(x(npi),x_u(npi))
	  	allocate(y(npj),y_v(npj))

		!----------------------------------!
		! Set size of computational domain !
		!----------------------------------!
	  	xl=0.05
	  	yl=0.1
		
		!---------------------------!
		! Calculate inner grid size !
		!---------------------------!
	  	dx=xl/real(npi-2)
	  	dy=yl/real(npj-2)

		!------------------------------------------------------!
		! Length variable for scalar points in the x-direction !
		!------------------------------------------------------!
	  	x(1)=0.
	  	x(2)=0.5*dx
	  	do i=3,npi-1
	    	x(i)=x(i-1)+dx
	  	end do
	  	x(npi)=x(npi-1)+0.5*dx

		!------------------------------------------------------!
		! Length variable for scalar points in the y-direction !
		!------------------------------------------------------!
	  	y(1)=0.
	  	y(2)=0.5*dy
	  	do j=3,npj-1
	    	y(j)=y(j-1)+dy
	  	end do
	  	y(npj)=y(npj-1)+0.5*dy 
	  	
	  	
	  	!-------------------------------------------------------------------!
		! Length variable for velocity components u(i,j) in the x-direction !
		!-------------------------------------------------------------------!
	  	x_u(1)=0.
	  	x_u(2)=0.
	  	do i=3,npi
	    	x_u(i)=x_u(i-1)+dx
	  	end do
		
		!-------------------------------------------------------------------!
		! Length variable for velocity components v(i,j) in the x-direction !
		!-------------------------------------------------------------------!
	 	y_v(1)=0.
	  	y_v(2)=0.
	  	do j=3,npj
	    	y_v(j)=y_v(j-1)+dy
	  	end do  
	  	
	end subroutine grid
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-------------Specifies boundary values for every iteration---------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine bound()
	 	
	 	!----------------!
		! Initialization !
		!----------------!
	    use declarations
	    implicit none
	    		
	  	!------------------------------------------------------!
		! von Neumann BC for temperature on east and west wall !
		!------------------------------------------------------!
	    T(npi,2:npj-1)=T(npi-1,2:npj-1)
	    T(1,2:npj-1)=T(2,2:npj-1)
	
	  	!--------------------------------------------------------!
		! von Neumann BC for temperature on north and south wall !
		!--------------------------------------------------------!
	  	T(2:npi-1,1)=T(2:npi-1,2)
	  	T(2:npi-1,npj)=T(2:npi-1,npj-1)
	
	end subroutine bound
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-----Prints out results to output directory for the current timestep-----!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine print(time)
	
		!----------------!
		! Initialization !
		!----------------!
		use declarations
		implicit none
		integer :: i,j
		real :: time
		character(len=30) :: rowfmt
		character(len=30) :: filename
		
		!------------------------!
		! Establish write format !
		!------------------------!
		write(rowfmt,'(A,I4,A)') '(',npj,'(2X,E15.6E2))'
		
		!-------------------!
		! Open output files !
		!-------------------!	    
	    write(filename,'("output/temp/temp_",f7.2,".dat")') time
	   	open(12,file=filename,status='unknown',RECL=(17*npj+120))
	   	

	    !----------------------------------------!
		! Write x- and y-vectors to output files !
		!----------------------------------------!
	    open(110,file='output/x.dat',status='unknown')    
	    open(111,file='output/y.dat',status='unknown')
	    write(110,'(1600F14.7)') x(:)
	    write(111,'(1600F14.7)') y(:)

		!----------------------------------!
		! Write temperature to output file !
		!----------------------------------!
	    do i=1,npi
	    	write(12,FMT=rowfmt) T(i,:)
	    end do
		

		!-------------!
		! Close files !
		!-------------!
	    close(110)
	    close(111)
	    close(12)
	    
	end subroutine print
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!--------Solves linear system of eq.s by line Gauss-Seidel method---------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine solve(fi,b,ae,aw,an,as,ap,istart,iend,jstart,jend)
		
		!----------------!
		! Initialization !
		!----------------!
	  	implicit none
	  	double precision, dimension(:,:), intent(in out) :: fi
	  	double precision, dimension(:,:), intent(in) :: b
	  	double precision, dimension(:,:), intent(in) :: ae,aw,an,as,ap
	  	integer, intent(in) :: istart,iend,jstart,jend
	  	integer :: i,j
	  	double precision, allocatable, dimension(:) :: ath, cmth
	  	double precision :: cth
		
		!-----------------------!
		! Allocate line vectors !
		!-----------------------!	
	  	allocate(ath(iend),cmth(iend))
	  	
		!----------------------------------------!
		! Solving the (e-w) lines from the south !
		!----------------------------------------!
	  	do j=jstart+1,jend-1 
			!-----------------------------!
			! Set values at west boundary !
			!-----------------------------!
	    	ath(istart)=0.   
	    	cmth(istart)=fi(istart,j) 
	    	
	    	!----------------------!
			! Forward substitution !
			!----------------------!
	    	do i=istart+1,iend-1
	      		ath(i)=ae(i,j)/(ap(i,j)-aw(i,j)*ath(i-1))
	      		cth=an(i,j)*fi(i,j+1)+as(i,j)*fi(i,j-1)+b(i,j)   
	      		cmth(i)=(aw(i,j)*cmth(i-1)+cth)/(ap(i,j)-aw(i,j)*ath(i-1))
	    	end do   
	    	
	    	!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do i=iend-1,istart+1,-1
	      		fi(i,j)=ath(i)*fi(i+1,j)+cmth(i)
	    	end do
	  	end do

		!----------------------------------------!
		! Solving the (e-w) lines from the north !
		!----------------------------------------!
	  	do j=jend-2,jstart+1,-1 
			
			!-----------------------------!
			! Set values at west boundary !
			!-----------------------------!
	    	ath(istart)=0. 
	    	cmth(istart)=fi(istart,j)
			
			!----------------------!
			! Forward substitution !
			!----------------------!
	    	do i=istart+1,iend-1
	      		ath(i)=ae(i,j)/(ap(i,j)-aw(i,j)*ath(i-1))
	      		cth=an(i,j)*fi(i,j+1)+as(i,j)*fi(i,j-1)+b(i,j)  
	      		cmth(i)=(aw(i,j)*cmth(i-1)+cth)/(ap(i,j)-aw(i,j)*ath(i-1))
	    	end do   

			!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do i=iend-1,istart+1,-1  
	      		fi(i,j)=ath(i)*fi(i+1,j)+cmth(i)
	    	end do
	  	end do      
		
		!-------------------------!
		! Deallocate line vectors !
		!-------------------------!
	  	deallocate(ath,cmth)
		
		!-----------------------!
		! Allocate line vectors !
		!-----------------------!
	  	allocate(ath(jend),cmth(jend))
		
		!---------------------------------------!
		! Solving the (n-s) lines from the west !
		!---------------------------------------!
	  	do i=istart+1,iend-1
			
			!------------------------------!
			! Set values at south boundary !
			!------------------------------!
	    	ath(jstart)=0.   
	    	cmth(jstart)=fi(i,jstart) 
	    	
	    	!----------------------!
			! Forward substitution !
			!----------------------!
	    	do j=jstart+1,jend-1
	      		ath(j)=an(i,j)/(ap(i,j)-as(i,j)*ath(j-1))
	      		cth=ae(i,j)*fi(i+1,j)+aw(i,j)*fi(i-1,j)+b(i,j)  
	      		cmth(j)=(as(i,j)*cmth(j-1)+cth)/(ap(i,j)-as(i,j)*ath(j-1))
	    	end do   
			
			!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do j=jend-1,jstart+1,-1
	      		fi(i,j)=ath(j)*fi(i,j+1)+cmth(j)
	    	end do   
	  	end do

		!---------------------------------------!
		! Solving the (n-s) lines from the east !
		!---------------------------------------!
	  	do i=iend-2,istart+1,-1
			
			!------------------------------!
			! Set values at south boundary !
			!------------------------------!
	    	ath(jstart)=0. 
	    	cmth(jstart)=fi(i,jstart)
   			
   			!----------------------!
			! Forward substitution !
			!----------------------!      
	    	do j=jstart+1,jend-1
	      		ath(j)=an(i,j)/(ap(i,j)-as(i,j)*ath(j-1))
	      		cth=ae(i,j)*fi(i+1,j)+aw(i,j)*fi(i-1,j)+b(i,j)  
	      		cmth(j)=(as(i,j)*cmth(j-1)+cth)/(ap(i,j)-as(i,j)*ath(j-1))
	    	end do   

			!-----------------------!
			! Backward substitution !
			!-----------------------!
	    	do j=jend-1,jstart+1,-1
	      		fi(i,j)=ath(j)*fi(i,j+1)+cmth(j) 
	    	end do   
	  	end do  
	  	
	  	!-------------------------!
		! Deallocate line vectors !
		!-------------------------!
	  	deallocate(ath,cmth)
	  	
	end subroutine solve
	
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!-------------------------------------------------------------------------!
	!-------------Calculates coefficients for the temp. equation--------------!
	!-------------------------------------------------------------------------!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	subroutine tcoeff(ae,aw,an,as,ap,b)
		
		!----------------!
		! Initialization !
		!----------------!
	  	use declarations
	  	implicit none
	  	double precision, dimension(:,:), intent(out) :: ae,aw,an,as,ap,b
	  	integer :: i,j
	  	real :: fw,fe,fs,fn,dw,de,ds,dn,areaw,areae,areas,arean,sp,su
	  	real :: ap_tilde,ap_zero
		
		!-----------------------!
		! Iterate through array !
		!-----------------------!
	  	do i=2,npi-1
	    	do j=2,npj-1
				
				!-------------------------!
				! Find area of cell faces !
				!-------------------------!
	      		areaw=y_v(j+1)-y_v(j)
	      		areae=areaw
	      		areas=x_u(i+1)-x_u(i)
	      		arean=areas
				
				!---------------------------!
				! Find convective mass flux !
				!---------------------------!
	      		fw=0.
	      		fe=0.
	      		fs=0.
	      		fn=0.

				!---------------------------------------------!
				! Find diffusion conductance by harmonic mean !
				!---------------------------------------------!
!				dw=((k(i-1,j)*k(i,j))/(k(i-1,j)*    &
!	            	(x(i)-x_u(i))+k(i,j)*(x_u(i)-x(i-1))))*areaw
!	        	de=((k(i,j)*k(i+1,j))/(k(i,j)*       &
!	            	(x(i+1)-x_u(i+1))+k(i+1,j)*(x_u(i+1)-x(i))))*areae
!	        	ds=((k(i,j-1)*k(i,j))/(k(i,j-1)*      &
!	            	(y(j)-y_v(j))+k(i,j)*(y_v(j)-y(j-1))))*areas
!	        	dn=((k(i,j)*k(i,j+1))/(k(i,j)*         &
!	            	(y(j+1)-y_v(j+1))+k(i,j+1)*(y_v(j+1)-y(j))))*arean
				dw=0.5*( k+k )/(x(i)-x(i-1))*areaw		
				de=0.5*( k+k )/(x(i+1)-x(i))*areae
				ds=0.5*( k+k )/(y(j)-y(j-1))*areas
				dn=0.5*( k+k )/(y(j+1)-y(j))*arean
	
				!-------------------!
				! Find source terms !
				!-------------------!
				sp=0.
	      		su=0.
	      		! Check if next to wall. If true, add heat loss through wall
	        	if (i==2) then
	          		su=su-(k_cup/cup_width)*areaw*(T(i,j)-T_amb)
	        	elseif (i==npi-1) then
	        		su=su-(k_cup/cup_width)*areae*(T(i,j)-T_amb)
	        	end if
	        	if (j==2) then
	          		su=su-(k_cup/cup_width)*areas*(T(i,j)-T_amb)
	        	elseif (j==npj-1) then
	        		su=su-(k_cup/cup_width)*arean*(T(i,j)-T_amb)
	        	end if
	
	            
	     		!----------------------------------!
				! Establish neighbour coefficients !
				!----------------------------------!
	        	aw(i,j)=dw+fw/2
	        	ae(i,j)=de-fe/2
	        	as(i,j)=ds+fs/2
	        	an(i,j)=dn-fn/2            
	            
	            !------------------------------!
				! Establish center coefficient !
				!------------------------------!
	        	ap_zero= rho*cp*areas*areaw/dt
	        	ap_tilde=aw(i,j)+ae(i,j)+as(i,j)+an(i,j)+fe-fw+fn-fs-sp
	        	ap(i,j)=ap_zero+ap_tilde
	      	
				!-----------------------------!
				! Establish total source term !
				!-----------------------------!
	     	 	b(i,j)=ap_zero*T(i,j)+su
	    	end do
	  	end do   
	  	
	end subroutine tcoeff


end module sub


!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
!----------------------------------------------------------------------------!
! SOLVES TRANSIENT CONVECTION-DIFFUSION PROBLEMS USING THE SIMPLE ALGORITHM  !
!----------------------------------------------------------------------------!
! Calculates velocities, pressure and temperature development in time on a   !
! 2D cartesian grid.                                                         !
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%!
 program CUPHEAT
		
	!----------------!
	! Initialization !
	!----------------!
	use declarations
	use sub
	implicit none
	double precision,allocatable,dimension(:,:) :: ae,aw,an,as,ap,b
	integer :: steps,it,jt,n
	integer :: txtclock, binclock
  	real    :: txttime, bintime
  	
  	!-------------------------------------------!
	! Initialize all parameters and create grid !
	!-------------------------------------------!
 	call init()
 	call grid()
 	 
	!---------------------------!
	! Allocate remaining arrays !
	!---------------------------!
	allocate(ae(npi,npj),aw(npi,npj),an(npi,npj),as(npi,npj))
  	allocate(ap(npi,npj),b(npi,npj))

	!----------------------------------!
	! Set node for convergence history !
	!----------------------------------!
  	it=npi/2
  	jt=npj/2
  	write(*,*) 'Node for convergence history:',it,jt
 	
 	call tick(txtclock)	
 	
 	
 	!---------------!
	! March in time !
	!---------------!
	do n = 1,nsteps 
					
		!---------------------!
		! Solve temp-equation !
		!---------------------!
		call tcoeff(ae,aw,an,as,ap,b)
		call solve(T,b,ae,aw,an,as,ap,1,npi,1,npj)	
			
		!-------------------------!
		! Specify boundary values !
		!-------------------------!
		call bound()			
		
		!---------------------------------------!
		! Print results for every 10th timestep !
		!---------------------------------------!
   		if (mod(n,100) == 0) then
      		write (*,'(f7.2,5g15.5)')  n*dt, &
             	T(it,jt)
            
        	call print(real(dt*n)) 
        	
    	end if    		    	      
	end do
	txttime=tock(txtclock)
    print *, 'System simulation time = ', txttime 
     
end program CUPHEAT
