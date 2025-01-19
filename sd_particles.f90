program formigas2d

! Assigning variables

integer i,j,k,l
integer n,m
integer Nf, Nl
real xmin,xmax,ymin,ymax
real dx, dy, dt, mass, friction
real, allocatable :: S(:,:), F(:,:)
integer, allocatable :: SAUX(:)
real, allocatable :: V(:,:), X(:,:)
real, allocatable :: RAND(:,:)
real, allocatable :: RANDFERO(:,:)
real, allocatable :: FORCA(:,:)
real, allocatable :: FICA(:)
real EVAP(1) 
real, allocatable :: KICKS(:,:)
integer, allocatable :: RANDFEROI(:,:) 
integer npast
real p,q
real time

! i,j,k = integer variables used in loops
! n,m =  number of subdivisions in x and y
! Nf = number of self-driven particles
! Nl = number of cells in the search board
! xmin, xmax, etc.. = vertices of the search board
! dx, dy, dt, mass, friction = space step in x and y, time-step, mass e friction coefficient of the particles
! S(:,:) = number of self-driven particles per cell 
! F(:,:) = marker used to account for the trace left by the self-driven particle per cell
! SAUX(:) = auxiliar vector used to identify which particle is located in each cell
! V(:,:) e X(:,:) = velocity and position of each particle in both space directions
! RAND(:,:) = random matrix to simulate a Brownian motion like movement of the particles
! FORCA(:,:) = random force to move the particles
! RANDFERO(:,:) = random matrix to populate the board with initial traces (real numbers)
! RANDFEROI(:,:) = random matrix to populate the board with initial traces (integer numbers)
! FICA(:) = random number between 0 and 1 that evaluate the odds of a particle to stay in a given position
! KICKS(:,:) = random matrix used to kick the particles out of forbiden zones according to the original algorithm
! npast = number of time-steps
! p = odds of the trace to evaporate
! q = odds of a particle to stay in an empty cell
! time = total simulation time

! Defining initial variables to start a simulation

n=50
m=50
Nl=n*m
Nf=50
xmin=0
ymin=0
xmax=10
ymax=10
dt=0.01
mass=0.1
friction=0.1
p=0.6
q=0.8
time=20
npast=time/dt

dx= (xmax-xmin)/(n-1)
dy= (ymax-ymin)/(m-1)

! Allocating variables

allocate(S(n,m))
allocate(SAUX(n*m))
allocate(F(n,m))
allocate(V(Nf,2))
allocate(X(Nf,2))
allocate(RAND(Nf,2))
allocate(RANDFERO(n,2))
allocate(RANDFEROI(n,2))
allocate(FORCA(Nf,2))
allocate(FICA(Nf))
allocate(KICKS(Nf,2))

! Initial distribution of the self-driven particles

call randomica(xmin,xmax,RAND(:,1),Nf,1)
call randomica(ymin,ymax,RAND(:,2),Nf,2)

do i=1,Nf
X(i,1)= RAND(i,1)
X(i,2)= RAND(i,2)
end do

! Assemblying the S field (number of particles per cell) 

do i=1,Nf

do k=1,n
do j=1,m
if(X(i,1).ge.(k-1)*dx.and.X(i,1).le.k*dx) then
if(X(i,2).ge.(j-1)*dy.and.X(i,2).le.j*dy) then
S(k,j)=S(k,j)+1
SAUX(k*j)=i
end if
end if
end do
end do

end do

write(*,*) 'Field S assemblying ok'

! Checking if there is more than a particle per cell

10 do k=1,n
   do j=1,m
   if(S(k,j).gt.1) then
   ! Kicking the last particle that has entered a populated cell
   X(SAUX(k*j),1)=X(SAUX(k*j),1)+0.5*dx
   go to 10   
   end if
   end do 
   end do

! Initial random distribution of traces to start the simulation

call randomica(1.0,20.0,RANDFERO(:,1),n,3)
call randomica(1.0,20.0,RANDFERO(:,2),n,4)

RANDFEROI=RANDFERO

do k=1,n
do j=1,m
do i=1,NF
if(k.eq.RANDFEROI(i,1).and.j.eq.RANDFEROI(i,2))then
F(k,j)= 1.0
end if
end do
end do
end do


! Opening a data file to write the position of the particles in time

open (2,file='trajectories.plt')
write(2,*) 'Variables="x","y","u","v"'

! Opening a data file to write the traces left by the particles in time

open (3,file='traces.plt')
write(3,*) 'Variables="x","y"'

! Making the main loop (time-evolution)

do k=1,npast

write(*,*) 'Current time-step',k

! Generating the random numbers to compute the forces acting in each particle

 call randomica(-dx,dx,FORCA(:,1),NF,k) ! Direção x (note o intervalo -dx,dx)
 call randomica(-dy,dy,FORCA(:,2),NF,k+1) ! Direção y (note o intervalo -dy,dy)

! Writing in the data file the number of the time zone

write(2,*) 'zone t="',k,'"'
write(3,*) 'zone t="',k,'"'

! Solving the velocity and position of each particle in a given time-step

do i=1,Nf 
 
! First, lets calculate the velocity

 call resvel(V(i,1),dt,FORCA(i,1),friction,mass) ! x direction
 call resvel(V(i,2),dt,FORCA(i,2),friction,mass) ! y direction

! Then, we calculate the position of the self-driven particles
 
 call respos(X(i,1),dt,V(i,1))	! x direction
 call respos(X(i,2),dt,V(i,2))	! y direction

! Avoiding particle scape from the board

if(X(i,1).le.xmin) then
X(i,1)=X(i,1)+0.1*dx
end if
if(X(i,1).ge.xmax) then
X(i,1)=X(i,1)-0.1*dx
end if
if(X(i,2).le.ymin) then
X(i,2)=X(i,2)+0.1*dy
end if
if(X(i,2).ge.ymax) then
X(i,2)=X(i,2)-0.1*dy
end if

end do

! Now, we must link the traces of the particles with their next movements

! Updating the S field (number of particles per cell) 

20 S=0.0

  do i=1,NF
    do l=1,n
     do j=1,m
       if(X(i,1).ge.(l-1)*dx.and.X(i,1).le.l*dx) then
        if(X(i,2).ge.(j-1)*dy.and.X(i,2).le.j*dy) then
         S(l,j)=S(l,j)+1
         SAUX(l*j)=i
        end if
       end if

! Checking if there is more than a particle per cell

         if(S(l,j).gt.1) then

   ! If there is more than a particle per cell we must kick the last invader particle out of the cell 
  
         X(SAUX(l*j),1)=X(SAUX(l*j),1)+ FORCA(i,1)*0.1
         X(SAUX(l*j),2)=X(SAUX(l*j),2)+ FORCA(i,2)*0.1
               go to 20   
         end if 
       end do
     end do
    end do

! Lets check now if the particle will stay on the new cell
! That will depend on the existance of traces left by other self-driven particles in the cell
! But before this procedure, lets pick a random number between 0 and 1 for each particle

call randomica(0.0,1.0,FICA,Nf,npast+2)
  
! Now, we will check for each particle if the new occupied cell has or not a trace left by other particles
    do l=1,n ! Doing a sweep in the x direction of all the cells of the board
     do j=1,m ! Doing a sweep in the y direction of all the cells of the board
       if(S(l,j).eq.1) then ! If the current cell is occupied
        if(F(l,j).eq.0) then ! And if there is not a trace in it    
! Then, the particle has a "s" chance to stay, where s < q, here "q" is the odd to stay if there is a trace on the cell
! In order to do so, lets consider that each particle will have an associated random number between 0 and 1. 
! These number are stored on the vector FICA(i), where i = 1,Nf. 
! If this number is smaller than q (chance = q) the particle stays in the cell
! It this number is bigger thant q (chance = 1-q) the particle must left the cell through a random kick      
	          if(FICA(SAUX(l*j)).gt.q) then
            X(SAUX(l*j),1)=X(SAUX(l*j),1)+ FORCA(SAUX(l*j),1)*0.1
	          X(SAUX(l*j),2)=X(SAUX(l*j),2)+ FORCA(SAUX(l*j),2)*0.1
            else    
	          X(SAUX(l*j),1)=X(SAUX(l*j),1)
	          X(SAUX(l*j),2)=X(SAUX(l*j),2)
            end if
          end if
          end if
       end do
     end do

! Finally, we must left a trace in the current cells occupied by the particles

do l=1,n
do j=1,m
   if(S(l,j).eq.1) then ! If there is a particle on the cell
   F(l,j)=1.0 ! Then, we put a trace on the current cell
   else ! If there is not a particle on the cell we must check if the trace will evaporate or not
   if(F(l,j).eq.1.0) then ! Here is the scenario of an empty cell with a current trace
   call randomica(0.0,1.0,EVAP,1,npast+2*k)  ! We ask for a random number between 0 and 1 
   if(EVAP(1).lt.p) then ! If this number is smaller than p (chance "p" to occur)
   F(l,j)=0.0 ! the trace evaporates
   end if
   end if
   end if
end do
end do
  
! Writting the position and velocity of the particles in a data file

do i=1,Nf
write(2,'(F12.4,F12.4,F12.4,F12.4)') X(i,1),X(i,2),V(i,1),V(i,2)
end do

! Writting the position of traces in a data file
do l=1,n ! Doing a sweep in x
do j=1,m ! Doing a sweep in y
if(F(l,j).eq.1.0) then ! If there is a trace on that cell
write(3,'(F12.4,F12.4)') l*dx+0.5,j*dy+0.5  ! We print the coordinates of the trace at the center of the cell
end if
end do
end do

! Closing the main loop (time-evolution)
end do

end

!*****************************************************************************************!
!                        Subroutines used in the code                                     !
!*****************************************************************************************!

! Subroutine used to solve the velocity of the particles through a 4th order 
! Runge-kutta scheme

subroutine resvel(a,b,c,d,e)
real a          ! velocity
real b			    ! time-step
real c          ! random force
real d			    ! friction
real e			    ! mass
real k1,k2,k3,k4 ! internal variables for the 4th order R.K. scheme

k1 = b*(c-d*a)/e
k2 = b*(c-d*(a+(k1*0.5)))/e
k3 = b*(c-d*(a+(k2*0.5)))/e
k4 = b*(c-d*(a+k3))/e

a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
end subroutine resvel

! Subroutine used to solve the position of each particle through a 4th order 
! Runge-kutta scheme

subroutine respos(a,b,c)
real a          ! position
real b			    ! time-step
real c          ! velocity component
real k1,k2,k3,k4 ! internal variables for the 4th order R.K. scheme

k1=b*(c)
k2=b*((0.5*k1)+c)
k3=b*((0.5*k2)+c)
k4=b*((k3)+c)
a=a+(1.0/6.0)*(k1+2.0*k2+2.0*k3+k4)
end subroutine respos

! Subroutine used to generate random numbers

subroutine randomica(a,b,c,n,d)
real a,b                    ! a,b = range of the random number
integer n, m
real c(n)                   ! c = generated random sequence
integer d,i,e
integer f(8)
integer, allocatable :: seed(:)

 call random_seed(size = m)
allocate (seed(m))

 CALL DATE_AND_TIME(values=f)
 CALL SYSTEM_CLOCK(count=e)

do i = 1,m
seed(i) =  47*d + f(8)*i*d*12000 + e*(3*d+i)
end do

 call random_seed(put = seed)

 call random_number(c)

 c = a+(b-a)*c

end subroutine randomica