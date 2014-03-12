! 3D program to solve the complex GPE
! Can restart from initial condition or any point in real time
! Can additionally imprint a vortex or two or a vortex ring
! Written by Joy Allen, revised last March, 2014

  module param
  
  implicit none

    integer, parameter :: dp = selected_real_kind(15,307)
    integer, parameter :: start = 0, restart = 0
    integer, parameter :: t_start = 0
    integer, parameter :: n_imag = 35000, nmax = 100000, Lr=75, Lz=75
    logical :: plotting = .true.  
    integer, parameter :: save_rate= 300  ! how often to save if plotting
    real(dp), parameter :: tau=0.001_dp, tau2 = 0.0001_dp
    real, parameter :: radr = 2.0_dp, radz = 2.0_dp 
    real(dp), parameter :: pi=4.0_dp*atan(1.0_dp)
    complex(dp), parameter :: eye = (0.0_dp,1.0_dp)

    ! if adding a vortex or vortex ring, set below to true
    ! and specify the objects parameters
    logical :: vortex = .true., vortex_ring = .false.
    real(dp), parameter :: x_0 = 0.0_dp, y_0 = 0.0_dp  ! vortex params
    real(dp), parameter :: x_cent = 0.0_dp, y_cent = 0.0_dp, &  ! ring params
                           ring_rad = 2.0_dp, z_pos = 1.0_dp, &
                           x_cent2 = 0.0_dp, y_cent2 = 0.0_dp, &
                           ring_rad2 = 2.0_dp, z_pos2 = -1.0_dp 
   
    ! set condensate parameters 
    ! some experimental quantities to keep notes 
    !real(dp), parameter :: hbar=1.054571628_dp*10**(-34.0_dp), m1=146.00_dp*10**(-27.0_dp)
    !real(dp), parameter :: omegax1=2.0_dp*pi*150.0_dp, omegay1=2.0_dp*pi*150.0_dp, &
    !omegaz1=2.0_dp*pi*150_dp  ! trapping frequencies
    !real(dp), parameter :: l_perp=sqrt(hbar/(m1*omegax1))  ! h.o. length
    real(dp), parameter :: lambda = 1.0_dp  ! trap aspect ratio
    real(dp), parameter :: g11 = 500.0_dp, mu1 = 8.0_dp
  
  end module param


program CGPE

 use omp_lib
 use param
 
 implicit none

 interface
   ! add any further functions in here.
     
   function ran2(idum)
     INTEGER idum!
      Integer  IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL  AM,EPS,RNMX
      PARAMETER  (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER  idum2,j,k,iv(NTAB),iy
      SAVE  iv,iy,idum2
      DATA  idum2/123456789/, iv/NTAB*0/, iy/0/
   end function ran2

   
     
end interface
!! END OF INTERFACES !!


  integer :: id, nthreads, n_start, radr_i, radz_k, n, i, j, a, k, l
  complex(dp) :: t, dt, dt2
  real(dp):: xL, yL, zL
  real(dp), dimension(0:Lr-1) :: x, y
  real(dp), dimension(0:Lz-1) :: z  
  complex(dp),dimension(0:Lz-1) :: vmult_ring, vmult_ring2
  complex(dp),dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: uold1, unew1
  real(dp) :: oldnorm1, E, ang_mom, norm, energy, mu, damp, &
              xd, yd, rd, rd1, rd2, rd3, xd2, yd2, dx, dy, dz, dens_max
  complex(dp) :: vmult, vmult2, vmultxy, vmultxz, vmultyz
  character (len=1) :: CHA
  character (len=2) :: CHA2
  character (len=3) :: CHA3
  character (len=4) :: CHA4
  logical :: real_time = .false.
 
    ! open some files for data storage
    
    open(12,file= "norm_imag.dat") 
    open(13, file= "norm_real.dat")
    open(17, file= "energy.dat")
    open(18, file = "ang_mom.dat")
    open(19,file = "maxdens.dat")
    open(20, file = "ang_mom_compo.dat")

  
  mu = mu1
  
  ! define the size of the box
  xL = -7.5_dp
  yL = -7.5_dp
  zL = -7.5_dp
  
  ! calculate dx,dy,dz
  dx = 2.0_dp*abs(xL)/(Lr)
  dy = 2.0_dp*abs(yL)/(Lr)
  dz = 2.0_dp*abs(zL)/(Lz)
  
  
  t = cmplx(0.0_dp, 0.0_dp, dp)
  
  !setup grid
  call setup_gridx(x, xL, dx)
  call setup_gridy(y, yL, dy)
  call setup_gridz(z, zL, dz)




   if (start.eq.0) then
     
     !Thomas fermi initial condition and perform imag time
     ! else we ignore this time loop and read in the wavefunction
     ! from file 1 - ics.dat
     call ics1(uold1,x,y,z, real_time,t)
     print*, "initial condition has been set up"
    
     open(10,file="initialcond.dat")
       do k=0,Lz-1
         write(10,*) z(k), abs(uold1(Lr/2,Lr/2,k))**2
       end do
     flush(10)
     close(10)
     
     !set the boundary conditions
     call bcs(uold1)
     call get_simpson(uold1,dx,dy,dz,oldnorm1)
   
     print*, 'The norm is initially = ', oldnorm1
     
     
     !pull condensate towards number of particles wanted
     if(.not.real_time) then
        uold1=uold1/sqrt(oldnorm1)
     endif
     
     open(110,file="initialcond2.dat")
       do i=0,Lr-1
         write(110,*) x(i), abs(uold1(i,Lr/2,Lz/2))**2
       end do
     flush(110)
     close(110)
     
     dt=-tau*eye
      
     t=real(n,dp)*dt
  
   !!!! IMAGINARY TIME PROPAGATION !!!!!
   !!!! with or without a vortex in fixed position !!!

  do n=0, n_imag
  
      t=real(n,dp)*dt
      
      call solve(uold1,unew1,dt,dx,dy,dz,t,x,y,z,mu,real_time)
      call get_simpson(unew1,dx,dy,dz,norm)
  
       unew1=unew1/sqrt(norm)
          write(12,*) -aimag(t), norm
       flush(12)
     
       if (mod(n,10*save_rate)==0) then 
         if (plotting) then
          call snapshot(n,unew1,x,y,z,dx,dz)
         end if
        end if
  
       if ((n.ne.0).and.(abs(norm-oldnorm1) < 1e-12_dp)) then
          print*, 'finished in imaginary time', norm-oldnorm1
          exit
       end if
  
       uold1=unew1
       oldnorm1=norm
     
  end do



 ! writing initial conditions to file
 open(1, file='ics.dat', form = "unformatted")
   do k =0, Lz-1
    do j = 0, Lr-1
     do i = 0, Lr-1
     write(1) unew1(i,j,k)
     end do
     end do
  enddo
 close(1)
 


   else ! read initial condition from file previously generated.
    open(1, file='ics.dat', form = "unformatted")     
    do k = 0, Lz-1
     do j = 0, Lr-1
      do i = 0, Lr -1  
      read(1) uold1(i,j,k)
      end do
     end do
    end do
    close(1)
   unew1=uold1
  end if 

  ! add a vortex?  If so - propagate for short time
  ! in imag time to remove density flucts.

  if (vortex) then
    dt2=-tau2*eye 
  
    do n=0, 750
     
     t=real(n,dp)*dt2
  
      call solve(uold1,unew1,dt2,dx,dy,dz,t,x,y,z,mu,real_time)
      call get_simpson(unew1,dx,dy,dz,norm)
  
    ! pull condensate towards number of particles wanted
      if(.not.real_time) then
       unew1=unew1/sqrt(norm)
      endif
  
  
     do j =0,Lr-2
       do i =0,Lr-2
          xd = x(i) - x_0
          yd = y(j) - y_0
          rd = sqrt(xd**2 + yd**2)
          xd2 = x(i) + x_0
          yd2 = y(j) + y_0
          rd2 = sqrt(xd2**2 + yd2**2)
        if (rd ==0.0_dp) then
          vmult = 1.0_dp
        else 
          vmult = cmplx(xd/rd,yd/rd,dp)
        end if
       if (rd2 ==0.0_dp) then
          vmult2 = 1.0_dp
       else
          vmult2 = 1.0_dp!cmplx(xd2/rd2,yd2/rd2,dp)  ! can add second one
     ! pade approx
     !    sqrt(rd**2*(0.3437_dp + 0.0286_dp*rd**2)/ &
     !               (1.0_dp + 0.3333_dp*rd**2 + 0.0286*rd**4))*
     !
       end if
          do k = 0,Lz -2
            unew1(i+1,j+1,k+1) = vmult*vmult2*abs(unew1(i+1,j+1,k+1))
          end do
       end do
     end do
    if (plotting) then
      call snapshot(n,unew1,x,y,z,dx,dz)
    end if 
   uold1=unew1

   end do
    ! write density with vortex to file
    open(41,file="vortex.dat")
      do i=0,Lr-1
        write(41,*) x(i), abs(unew1(i,Lr/2,Lz/2))**2
      end do
     close(41)
 
  end if

  if (vortex_ring) then
    dt=-tau*eye
  
     do n=0, 500
     
     t=real(n,dp)*dt
  
      call solve(uold1,unew1,dt,dx,dy,dz,t,x,y,z,mu,real_time)
      call get_simpson(unew1,dx,dy,dz,norm)
  
    !pull condensate towards number of particles wanted
    if(.not.real_time) then
     unew1=unew1/sqrt(norm)
    endif
  
  
  
      ! add the vortex ring
      if(n.eq.0) then  
      do k = 0, Lz - 2
       do j = 0, Lr - 2
         do i = 0, Lr - 2
            xd = x(i) - x_0
            yd = y(j) - y_0
           rd = sqrt(xd**2 + yd**2) - ring_rad
           rd2= sqrt((x(i)-x_cent)**2 + (y(j)-y_cent)**2)
          if ((rd2 - ring_rad) ==0.0_dp) then
            vmult_ring(k) = 1.0_dp      
          else 
            vmult_ring(k) = ((rd2 - ring_rad) + eye*(z(k)-z_pos))/ &
                             (sqrt((rd2 - ring_rad)**2 + (z(k)-z_pos)**2))
           end if
           rd = sqrt(xd**2 + yd**2) - ring_rad2
           rd3= sqrt((x(i)-x_cent2)**2 + (y(j)-y_cent2)**2)
          if ((rd3 - ring_rad2) ==0.0_dp) then
            vmult_ring2(k) = 1.0_dp      
          else 
            vmult_ring2(k) = ((rd3 - ring_rad2) + eye*(z(k)-z_pos2))/ &
                             (sqrt((rd3 - ring_rad2)**2 + (z(k)-z_pos2)**2))
           end if
    
              unew1(i+1,j+1,k+1) = vmult_ring2(k)*vmult_ring(k)*abs((unew1(i+1,j+1,k+1)))
          end do
         end do
       end do
      end if
      
      if (plotting) then
       call snapshot(n,unew1,x,y,z,dx,dz)
      end if
    
        
       uold1=unew1
      end do
     ! write density containing vortex ring to file 
      open(41,file="vortex_ring.dat")
        do i=0,Lr-1
          write(41,*) x(i), abs(unew1(i,Lr/2,Lz/2))**2
        end do
       close(41)
 
  end if

  print*, 'end of imag imprinting'

   ! SWITCH TO REAL TIME! 
    real_time = .true.
    dt=cmplx(tau,0.0_dp, dp)
    t = cmplx(0.0_dp, 0.0_dp, dp)
  


  !!!! REAL TIME PROPAGATION  !!!!
  t=real(n,dp)*dt 
  print *, "real time"

    if (restart.eq.1) then
      ! if you want to restart the code from a previous point
      ! switch restart to true and specify the file number in t_start

    n_start = t_start*save_rate
    if (t_start.eq.0) then
      open (1,file = "psi-0000000.dat", form = 'unformatted')
       else if (t_start.lt.10) then
        write (CHA,"(i1)") t_start
       open (1,file = "psi-000000"//CHA//".dat", &
                           form = 'unformatted')
      else if (t_start.lt.100) then
        write (CHA2,"(i2)") t_start
        open (1,file = "psi-00000"//CHA2//".dat",ACCESS='STREAM', &
                             form = 'unformatted')
      else if (t_start.lt.1000) then
        write (CHA3,"(i3)") t_start
       open (1,file = "psi-0000"//CHA3//".dat", ACCESS='STREAM', &
                            form = 'unformatted')
       else if (t_start.lt.10000) then
        write (CHA4,"(i4)") t_start
       open (1,file = "psi-000"//CHA4//".dat", ACCESS='STREAM', &
                            form = 'unformatted')
     end if

      do k = 0, Lz-1
        do j = 0, Lr-1 
          do i = 0, Lr-1 
           read (1) unew1(i,j,k)
          end do
        end do
     end do

       uold1 = unew1
       close(1)
     else
       n_start = 0
     end if

  ! REAL TIME TIME LOOP !

  do n=n_start, nmax

    t=real(n,dp)*dt  
    call solve(uold1,unew1,dt,dx,dy,dz,t,x,y,z,mu,real_time)
     
    ! print out wavefunction every save_rate
    ! for analysis
    call snapshot(n,unew1,x,y,z,dx,dz)
    
    dens_max = maxval(abs(unew1)**2)
    ! print out other useful condensate values
    if(mod(n,save_rate/2)==0) then
     
      call get_simpson(unew1,dx,dy,dz,norm)
      call get_energy(unew1,x,y,z,dx,dy,dz,E,real_time,t)
      call get_angmom(unew1,x,y,z,dx,dy,dz,ang_mom,t)
     
      ! write out parameters to check conservation of
      ! mass, energy, and angular momentum
          write(13,*) real(t), norm
          write(17,*) real(t), E
          write(18,*) real(t), ang_mom
          write(19,*) real(t), dens_max
    end if
    flush(13)
    flush(19)
    uold1=unew1

    
  enddo 

! close all files which are still open!

close(13)
close(17)
close(18)
close(19)
close(20)




!!!!!!!! END OF MAIN PROGRAM !!!!!!

!contains
end program  CGPE

!center the grids at 0
subroutine setup_gridx(x, xL, dx)
  
  use omp_lib
  use param
  
  implicit none
  
  real(dp), dimension(0:Lr-1), intent(out) :: x
  integer :: i
  real(dp), intent(in) :: xL, dx
  
    do i=0,Lr-1
      x(i) =xL+real(i,dp)*dx
      print*, x(i)
    enddo
  
  return
end subroutine setup_gridx

subroutine setup_gridy(y, yL, dy)

  use omp_lib
  use param
  
  implicit none
  
  real(dp), dimension(0:Lr-1), intent(out) :: y
  integer :: i
  real(dp), intent(in) :: yL, dy
  
      do i=0,Lr-1
        y(i) =yL+real(i,dp)*dy
      enddo
  
  return
end subroutine setup_gridy

subroutine setup_gridz(z, zL, dz)

  use omp_lib
  use param
  
  implicit none
  
  real(dp), dimension(0:Lz-1), intent(out) :: z
  integer :: i
  real(dp), intent(in) :: zL, dz
  
  
  
    do i=0,Lz-1
       z(i) =zL+real(i,dp)*dz
    enddo
  
  
  return
end subroutine setup_gridz



! Thomas fermi initial profile

subroutine ics1(u,x,y,z,real_time,t)

  use param
  use omp_lib
  
  implicit none
  
  complex(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1), intent (out) :: u
  real(dp), dimension(0:Lr-1), intent(in) :: x,y
  real(dp), dimension(0:Lz-1), intent(in):: z
  real(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1)::trap
  complex(dp), intent (in) :: t
  real(dp) :: r,h
  integer :: i,j,k
  logical :: real_time
  print*, "in ics "
  
  r = sqrt(2.0_dp*mu1)
  
  call trap1(x,y,z,trap,real_time,t)
  
  do k= 0, Lz-1
    do j= 0, Lr-1
      do i= 0, Lr-1
          if (trap(i,j,k).lt.mu1) then
          u(i,j,k) =cmplx(sqrt((mu1 - trap(i,j,k))/(g11)),0.0_dp,dp)
           else 
          u(i,j,k) = 0.0_dp
        end if    
      enddo
    enddo
  enddo
  
  print*,"finished in ics"
  
  return
end subroutine ics1



subroutine trap1(x,y,z,tra,real_time,t)

 use omp_lib
 use param

 implicit none
 
 real(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: tra
 real(dp), dimension(0:Lr-1), intent (in):: x, y
 real(dp), dimension(0:Lz-1), intent (in):: z
 complex(dp), intent (in) :: t
 integer :: i, j, k
 logical :: real_time

  
   do k=0,Lz-1
    do j=0,Lr-1
        do i=0,Lr-1
          tra(i,j,k)=(0.5_dp)*((x(i))**2.0_dp + (y(j))**2.0_dp+ lambda**2.0_dp*z(k)**2.0_dp)
           
       enddo
    enddo
  enddo

end subroutine trap1

!boundary conditions
subroutine bcs(k)

 use omp_lib
 use param 

 implicit none
 
 complex(dp), intent(inout), dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: k
 integer :: j
 
 ! subroutine to define boundary conditions
 do j=0,Lr-1
    k(0,j,0)=0.0_dp
    k(0,j,Lz-1)=0.0_dp
    k(j,0,0)=0.0_dp
    k(j,Lr-1,0)=0.0_dp
    k(j,0,Lz-1)=0.0_dp
    k(j,Lr-1,Lz-1)=0.0_dp
    k(Lr-1,j,0)=0.0_dp
    k(Lr-1,j,Lz-1)=0.0_dp
 enddo
 
 do j=0,Lz-1
    k(0,0,j)=0.0_dp
    k(0,Lr-1,j)=0.0_dp
    k(Lr-1,0,j)=0.0_dp
    k(Lr-1,Lr-1,j)=0.0_dp
 enddo
return
end subroutine bcs



!solve the wave functions using the fourth order Runge-Kutta method
subroutine solve(uold1,unew1,dt,dx,dy,dz,t,x,y,z,mu,real_time)

 use omp_lib
 use param

 implicit none

 complex(dp), intent(in), dimension(0:Lr-1,0:Lr-1,0:Lz-1):: uold1
 complex(dp), intent(out), dimension(0:Lr-1,0:Lr-1,0:Lz-1):: unew1
 complex(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1),save :: k11,k21,k31,k41
 complex(dp) :: dt
 real(dp), intent(in) :: dx, dy, dz, mu
 real(dp), intent(in), dimension(0:Lr-1) :: x, y
 real(dp), intent(in), dimension(0:Lz-1) :: z
 complex(dp), intent(in) :: t
 logical :: real_time

   call f1(uold1,k11,dx,dy,dz,t,x,y,z,mu,real_time)
   call f1(uold1+0.5_dp*dt*k11,k21,dx,dy,dz,t,x,y,z,mu,real_time)
   call f1(uold1+0.5_dp*dt*k21,k31,dx,dy,dz,t,x,y,z,mu,real_time)
   call f1(uold1+dt*k31,k41,dx,dy,dz,t,x,y,z,mu,real_time)

   unew1=uold1+dt*(k11/6.0_dp+k21/3.0_dp+k31/3.0_dp+k41/6.0_dp)

 return
end subroutine solve



subroutine f1(u1,kappa,dx,dy,dz,t,x,y,z,mu,real_time)

  use omp_lib
  use param

  implicit none
  
  complex(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: kappa
  complex(dp), intent(in), dimension(0:Lr-1,0:Lr-1,0:Lz-1):: u1
  real(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: trap
  real(dp), intent(in) :: dx, dy, dz
  real(dp), intent(in), dimension(0:Lr-1) :: x, y
  real(dp), intent(in), dimension(0:Lz-1) :: z
  complex(dp),dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: laplace_x, laplace_y, &
                                           laplace_z, ux, uy, uz
  integer :: i,j,k
  real(dp) :: damp, mu
  complex(dp), intent(in) :: t
  logical :: real_time


  ! call boundary conditions and trap function
  call bcs(kappa)
  call trap1(x,y,z,trap, real_time,t)
  
  ! put in value of damping if required in real time.
  ! decide whether to include mu in real time prop or not.
  if(real_time) then 
    damp = 0.0_dp
    mu = 0.0_dp
  else 
    damp =  0.0_dp
    mu = mu1
  end if
 
  !$omp parallel do
do k=1, Lz-2
  do j=1, Lr-2
    do i=1, Lr-2
      ux(i,j,k)=(u1(i+1,j,k) - u1(i-1,j,k))/(2.0_dp*dx)
      uy(i,j,k)=(u1(i,j+1,k) - u1(i,j-1,k))/(2.0_dp*dy)
      uz(i,j,k)=(u1(i,j,k+1) - u1(i,j,k-1))/(2.0_dp*dz)
      laplace_x(i,j,k)=(u1(i+1,j,k)-2.0_dp*u1(i,j,k)+u1(i-1,j,k))/(dx*dx)
      laplace_y(i,j,k)=(u1(i,j+1,k)-2.0_dp*u1(i,j,k)+u1(i,j-1,k))/(dy*dy)
      laplace_z(i,j,k)=(u1(i,j,k+1)-2.0_dp*u1(i,j,k)+u1(i,j,k-1))/(dz*dz)
      end do
    end do
  end do
  !$omp end parallel do

  !$omp parallel do
do k=1, Lz-2
  do j=1, Lr-2
    do i=1, Lr-2
      kappa(i,j,k)=-eye/(1.0_dp + eye*damp)* &
                ((trap(i,j,k)+g11*abs(u1(i,j,k))**2-mu)*u1(i,j,k) &
                -0.5_dp*(laplace_x(i,j,k)+laplace_y(i,j,k)+laplace_z(i,j,k)))
    enddo
  enddo
enddo
  !$omp end parallel do


return
end subroutine f1
           

!calculate the norm using simpson

subroutine get_simpson(u,dx,dy,dz,integral)

  use omp_lib
  use param


 implicit none

  integer :: i,j,k,ix,iy,iz
  real(dp),intent(out) ::  integral
  real(dp) :: normsum
  real(dp), dimension(0:Lr-1) :: sumy
  real(dp), dimension(0:Lz-1) :: sumz
  real(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: dens
  complex(dp),dimension(0:Lr-1,0:Lr-1,0:Lz-1),intent(in)::u
  real(dp), intent(in) :: dx, dy, dz


  integral = 0.0_dp
  dens = abs(u)**2

  do k = 0, Lz-1
    do j = 0, Lr -1
    sumy(j) = dens(0,j,k) + dens(Lr-1,j,k)
      do ix = 1, Lr-3, 2
      sumy(j) = sumy(j) + 4.0_dp*dens(ix, j, k) + &
                2.0_dp*dens(ix +1, j, k)
      end do
    end do
    
    sumz(k) = sumy(0) + sumy(Lr -1)
    do iy = 1,Lr-3,2
      sumz(k) = sumz(k) + 4.0_dp*sumy(iy) + 2.0_dp*sumy(iy+1)
    end do
  end do

  normsum = sumz(0) + sumz(Lz-1)
  do iz = 1, Lz-3,2
    normsum = normsum + 4.0_dp*sumz(iz) + 2.0_dp*sumz(iz + 1)
  end do
  integral = dx*dy*dz*normsum/27.0_dp

  end subroutine

  subroutine get_energy(u,x,y,z,dx,dy,dz,E,real_time,t)
  !a subroutine to compute the energy
 
  use omp_lib
  use param
 
  implicit none
 
  integer :: i, j, k, ix, iy, iz
  real(dp), intent(in) :: dx, dy, dz
  real(dp), intent(out) :: E
  complex(dp),dimension(0:Lr-1,0:Lr-1,0:Lz-1),intent(in)::u
  real(dp), intent(in), dimension(0:Lr-1) :: x, y
  real(dp), intent(in), dimension(0:Lz-1) :: z
  complex(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1):: ux, uy, uz
  real(dp), dimension(0:Lr-1,0:Lr-1,0:Lz-1) :: trap
  logical :: real_time
  real(dp) :: damp, normsum
  real(dp), dimension(0:Lr - 2):: sumy
  real(dp), dimension(0:Lz - 2) :: sumz
  complex(dp), intent(in) :: t 
 ! need to convert this to simpson!!!

 ! Check if this is accurate!

  E = 0.0_dp

  if(real_time) then
    damp = 0.0_dp
  else
    damp =  0.0_dp
  end if

  call trap1(x,y,z,trap,real_time,t)

  do k = 1, Lz-2
    do j = 1, Lr-2
      do i = 1, Lr-2
        ux(i,j,k)=(u(i+1,j,k) - u(i-1,j,k))/(2.0_dp*dx)
        uy(i,j,k)=(u(i,j+1,k) - u(i,j-1,k))/(2.0_dp*dy)
        uz(i,j,k)=(u(i,j,k+1) - u(i,j,k-1))/(2.0_dp*dz)
      end do
    end do
  end do

  ! now perform the simpsons integral

  do k = 1, Lz-2
    do j = 1, Lr -2
      sumy(j) = 0.0_dp
      do ix = 1, Lr -3, 2
        sumy(j) = sumy(j) + 2.0_dp*(0.5_dp*(abs(ux(ix,j,k)+uy(ix,j,k)+uz(ix,j,k))**2) + &
                   0.5_dp*g11*abs(u(ix,j,k))**4 &
                  + (trap(ix,j,k))*abs(u(ix,j,k))**2) + &
                  4.0_dp*(0.5_dp*(abs(ux(ix+1,j,k)+uy(ix+1,j,k)+uz(ix+1,j,k))**2) + &
                   0.5_dp*g11*abs(u(ix+1,j,k))**4 &
                  + (trap(ix+1,j,k))*abs(u(ix+1,j,k))**2) 
       end do
    end do
    sumz(k) = 0.0_dp
    do iy = 1, Lr -3, 2
      sumz(k) = sumz(k) + 2.0_dp*sumy(iy) + 4.0_dp*sumy(iy + 1)
    end do
  end do


  normsum = 0.0_dp
  do iz = 1, Lz -3 ,2 
    normsum = normsum + 2.0_dp*sumz(iz) + 4.0_dp*sumz(iz)
  end do

  E = normsum*dx*dy*dz/27.0_dp

end subroutine

subroutine get_angmom(u,x,y,z,dx,dy,dz,ang_mom, t)

use omp_lib
use param

implicit none

 integer :: i, j, k
 real(dp), intent(in) :: dx, dy, dz
 real(dp), intent(out) :: ang_mom
 real(dp) :: ang_mom_x, ang_mom_y, ang_mom_z
 complex(dp),dimension(0:Lr-1,0:Lr-1,0:Lz-1),intent(in)::u
 real(dp), intent(in), dimension(0:Lr-1) :: x, y
 real(dp), intent(in), dimension(0:Lz-1) :: z
 complex(dp) :: ux, uy, uz
 complex(dp), intent(in) :: t

 ! check if this is accurate!

  ang_mom = 0.0_dp

do k = 0, Lz-2
  do j = 0, Lr -2
    do i = 0, Lr -2
      ux=(u(i+1,j,k) - u(i-1,j,k))/(2.0_dp*dx)
      uy=(u(i,j+1,k) - u(i,j-1,k))/(2.0_dp*dy)
      uz=(u(i,j,k+1) - u(i,j,k-1))/(2.0_dp*dz)
      
      ang_mom_x = (z(k)*uy - y(j)*uz)
      ang_mom_y = (x(i)*uz - z(k)*ux)
      ang_mom_z = (y(j)*ux - x(i)*uy)
      
      ang_mom = sqrt(ang_mom_x**2 + ang_mom_y**2 + ang_mom_z**2)
    end do
  end do
end do

   write(20,*) real(t), ang_mom_x, ang_mom_y, ang_mom_z
   flush(20)
   
  end subroutine
 
  subroutine snapshot(n,unew1,x,y,z,dx,dz)

   use param

   implicit none


    integer :: i, j, k
    integer, intent(in) :: n
    real(dp), intent(in), dimension(0:Lr-1) :: x, y
    real(dp), intent(in), dimension(0:Lz-1) :: z
    real(dp), intent(in) :: dx, dz
    complex(dp),dimension(0:Lr-1,0:Lr-1,0:Lz-1),intent(in)::unew1
    real(dp), dimension(0:Lr-1,0:Lr-1) :: PHAS_C
    real(dp), dimension(0:Lr-1,0:Lz-1) :: PHAS_Z
    real(dp) ::  XP, YP, ZP, radr_i, radz_k
    character(15) :: solution='soln*******.dat'
    character(15) :: solutioz='solz*******.dat'
    character(15) :: density = 'dens*******.dat'
    character(15) :: phas = 'phas*******.dat'
    character(16) :: isodens = 'isod*******.dat'
    character(24) :: psi = 'psi-*******.dat'
    character(15) :: phaz = 'phaz*******.dat'

    ! print out densities, phase and wavefunction
   
   ! find index for isosurface plots
    radr_i =NINT(radr/dx)
    radz_k =NINT(radz/dz)
  
  
   if(mod(n,save_rate)==0) then
     write(psi(5:11),'(i7.7)')n/(save_rate)
      open(24,file=psi,ACCESS='STREAM', &
                            form="unformatted")
       do k = 0, Lz - 1
         do j = 0, Lr - 1
           do i= 0, Lr - 1
             write(24)unew1(i,j,k)
           end do
         end do
       end do
      close(24)
    endif
  
  
      ! print out 1d densities to check
     if(mod(n,10*save_rate)==0) then
      write(solution(5:11),'(i7.7)')n/(10*save_rate)
      open(11,file=solution)
        do i=0,Lr-1
          write(11,'(3e20.9E3)')x(i), abs(unew1(i,Lr/2,Lz/2))**2
        enddo
      close(11)
     endif
     if(mod(n,10*save_rate)==0) then
      write(solutioz(5:11),'(i7.7)')n/(10*save_rate)
      open(11,file=solutioz)
        do k=0,Lz-1
          write(11,'(3e20.9E3)')z(k), abs(unew1(Lr/2,Lr/2,k))**2
        enddo
      close(11)
     endif
  
  
     if(mod(n,2*save_rate)==0) then
        write(density(5:11),'(i7.7)')n/(2*save_rate)
        open(14, file=density)
        
        do j = 0, Lr-1
            do i=0,Lr-1
              write(14,'(3e20.9E3)')x(i), y(j), abs(unew1(i,j,Lz/2))**2
              end do
         end do
     end if
  
     if(mod(n,10*save_rate)==0) then
       write(isodens(5:11),'(i7.7)')n/(10*save_rate)
        open(16, file=isodens)
  
        do k = NINT(Lz/2 - radz_k), NINT(Lz/2 + radz_k)
         do j = NINT(Lr/2 - radr_i), NINT(Lr/2 + radr_i)
           do i = NINT(Lr/2 - radr_i), NINT(Lr/2 + radr_i)
  
              write(16,'(e20.9E3)') abs(unew1(i,j,k))**2
           end do
         end do
        end do
     end if
  
      if(mod(n,10*save_rate)==0) then
        write(phas(5:11),'(i7.7)')n/(10*save_rate)
        open(15, file=phas)
        do j = 0, Lr -2
            do i=0,Lr -2
            XP = Real(unew1(i+1,j+1,Lz/2+1),DP)
            YP = AImag(unew1(i+1,j+1,Lz/2+1))
            if ((Abs(XP)<1.E-10_DP).and.(Abs(XP)<1.E-10_DP)) then  
             PHAS_C (i,j)=0.0_DP
           else
             PHAS_C (i,j)=ATan(YP/Abs(XP))
           end if
           if (XP<0.0_DP) then
            if (YP>0.0_DP) then
               PHAS_C (i,j)=PI-PHAS_C (i,j)
             else
               PHAS_C (i,j)=-PI-PHAS_C (i,j)
              end if
            end if
  
              write(15,'(3e20.9E3)')x(i), y(j), PHAS_C(i,j)
              end do
         end do
  
       if(vortex_ring) then
           write(phaz(5:11),'(i7.7)')n/(10*save_rate)
        open(15, file=phaz)
        do k = 0, Lz -2
            do j=0,Lr -2
            ZP = Real(unew1(Lr/2+1,j+1,k+1),DP)
            YP = AImag(unew1(Lr/2+1,j+1,k+1))
            if ((Abs(ZP)<1.E-10_DP).and.(Abs(YP)<1.E-10_DP)) then  
             PHAS_Z (j,k)=0.0_DP
           else
             PHAS_Z (j,k)=ATan(YP/Abs(ZP))
           end if
           if (ZP<0.0_DP) then
            if (YP>0.0_DP) then
               PHAS_Z (j,k)=PI-PHAS_Z (j,k)
             else
               PHAS_Z (j,k)=-PI-PHAS_Z (j,k)
              end if
            end if
  
              write(15,'(3e20.9E3)')y(j), z(k), PHAS_Z(j,k)
              end do
         end do
         end if
  
  
     end if
 

    end subroutine




   function ran2(idum)
    ! a function to calculate random numbers - if ever required.
    ! use omp_lib
    ! use param
    ! implicit none

      INTEGER idum!
      Integer  IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      REAL  AM,EPS,RNMX
      PARAMETER  (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
      IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,&
      NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER  idum2,j,k,iv(NTAB),iy
      SAVE  iv,iy,idum2
      DATA  idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
     return
  end function ran2

   subroutine Seedset (ISE)

     integer, intent (out) :: ISE
     integer :: VALUES (8)
     character (len=10) :: DATE, TIME, ZONE

! Initializes seed of random number generator, setting it to a value dependent
! upon the system clock.

     ISE = 0
     do while (ISE==0)
        call Date_And_Time (DATE, TIME, ZONE, VALUES)
        ISE = -(1000*VALUES(7) + VALUES(8))
     end do

   end subroutine Seedset
