      module parameters
      parameter(nxmax=500200,ncmax=99,nenmax=10000,nfmax=19)
      parameter(ntet=1)   ! This will save a lot of space !!!
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,half=0.5d0,ifak=500)
      end module parameters
      
      module realconstants1
      integer*4 iz,nbf,nx1,nc1,nc2,nt,nx,nprint,nc,ngob,ntfin,nbmax,
     1          ngob1,nerg,mfixed,nnauto,llauto
      end module realconstants1
      
      module extra
      double precision ddel
      end module extra
      
      module cmplxconstants
      complex*16 czero,chalf,cone,ctwo,cthree,cfour,ci 
      complex*16 ch,cst,cst1,cst2,cst3,cst4,cst5,h2,h3,h24,h1d5,hd25,c27
      end module cmplxconstants
      
      module realconstants2
      double precision pi,h,dt,tau,zet,gbr,xtime
      end module realconstants2
      
      module complexarray1
      use parameters
      complex*16, allocatable :: a(:),b(:),c(:),r(:),gam(:),u(:)
      complex*16, allocatable :: cmn(:),cpl(:)   
      end module complexarray1
      
      module complexarray2
      use parameters
      complex*16, allocatable :: a1(:),b1(:),c1(:),r1(:),gam1(:),u1(:)
      end module complexarray2
      
      module complexarray3
      use parameters
      complex*16, allocatable :: q(:,:),v(:,:),ff(:,:),bb(:,:)
      end module complexarray3
      
      module realarray1 
      use parameters
      double precision, allocatable :: x(:)
      double precision, allocatable :: g(:)
      double precision, allocatable :: f(:)
      double precision, allocatable :: fd(:,:)
      double precision, allocatable :: s(:)
      double precision, allocatable :: cek3(:,:,:)      
      end module realarray1 
      
      module keys
      integer key1,key2,key3,key4,key5
      end module keys
      
      module in1
      double precision w1,e1,s1,ton1,toff1,del1,
     1                 w2,e2,s2,ton2,toff2,del2,agbr 
      end module in1
      
      module discr
      use parameters
      integer nf,nn(nfmax),ll(nfmax)
      character *2, target
      end module discr
       
      module in2 
      double precision emin,de,enelec
      end module in2 
      
      module alphas
      double precision alpha1,alpha2
      integer ndelay
      end module alphas
      
      module pgcommon 
      use parameters    
      dimension pg(0:ntet,0:2*ncmax),theta(0:ntet)
      end module pgcommon  
      
      module factorials
      use parameters 
      double precision FACT(IFAK)
      end module factorials
      
      module fgercm
      integer IERR,IERCT
      end module fgercm
      
      module adis
      use parameters
      complex*16 cdst(0:ncmax)
      end module adis
      
      module VDD
      use parameters
      double precision VD(nxmax+1000)
      end module VDD
      
      module pecommon
      use parameters
      double precision dw(nxmax+1000)
      end module pecommon

      module rpcommon
      use parameters
      double precision RP(nxmax+1000)
      end module rpcommon
      
      module DAT1
      double precision Z,DRA
      integer NR,NRP
      end module DAT1
      
      module phcommon
      use parameters
      double precision ph(0:ncmax)
      end module phcommon
      
      module pulseparams
      use parameters
      character*1 shape1up,shape1down,shape2up,shape2down
      integer n1up,n1plat,n1down,n2up,n2plat,n2down,nextra
      integer ntfinal,ntfinpulse
      double precision, allocatable  :: timeval(:)
      double precision, allocatable  :: envelope1(:),envelope2(:)
      double precision, allocatable  :: field1(:),field2(:),fieldtot(:)
      double precision ww1,period1,ee1,alph1,cep1,rr1
      double precision ww2,period2,ee2,alph2,cep2,rr2
      double precision tstart1,tend1,tstart2,tend2,tpulses,tfinal
      double precision tup1,tdown1,cep1rad
      double precision tup2,tdown2,cep2rad
      integer  istart1,iup1,idown1,iend1
      integer  istart2,iup2,idown2,iend2
      double precision gvalstart1,gvalend1
      double precision gvalstart2,gvalend2
      double precision correction
      end module pulseparams

      module fieldFT
C     25000 is hardcoded; should be replaced with a ndim or
C     or some other variable for dynamical allocation
      integer n, nomega
      double precision deltat, wmin, wmax, wdel, gsin, gcos
      double precision dimension, w(25000),greal(25000),gimag(25000)
      end module fieldFT
      
      PROGRAM NE1
C**********************************************************************
C  Nonstationary processes with neon atom in pulsed laser fields      *
C  based on version 1.1  Drake,  May-June 2005 and Li code of 2010    *
C  This is an extension for nonzero magnetic quantum number.          *
C**********************************************************************
      use realconstants1
      use extra
      use cmplxconstants
      use realconstants2
      use complexarray1
      use complexarray2
      use complexarray3
      use realarray1
      use keys
      use in1
      use discr
      use in2
      use alphas
      use pulseparams
      
      implicit double precision (a-h,o-z)
      complex*16 bet,xj
      integer*4 itime_start, itime_end, itime
!      character*2 target
      
      include 'omp_lib.h'
 
      if(key5.gt.1) open(30,file='mathem.scr',status='unknown')
      open(40,file='read.me',status='unknown')
      open(50,file='tdse.inp',status='old')
      open(60,file='wfn.out',status='unknown')
      open(70,file='betas.out',status='unknown')
      open(80,file='pulse.inp',status='old')
      open(90,file='pulse.out',status='unknown')
      open(100,file='fieldFT.out',status='unknown')
      open(230,file='intermediate.out',status='unknown')
  
      call readkb      ! new input subroutine
      if (target == ' H') then
        open(10,file='hyd3.out',status='unknown')
        open(20,file='hyd1.out',status='unknown')
      endif
      if (target == 'Ne') then
        open(10,file='ne3.out',status='unknown')
        open(20,file='ne1.out',status='unknown')
      endif
      
      if(target /= ' H' .and. target /= 'Ne') then
        print *,'This code only works for H or Ne.  Program stopped.'
        write(6,98211) target
98211   format(/,' your target = ',a2) 
        stop
      endif

      call cnstnt      ! definition of constants
      call arrays      ! initiate arrays
      call pulsekb     ! set up the electric field for the laser pulse(s)
!
!  interface to the previous part of the code
      f = fieldtot

      if (target == ' H') call radhyd   ! functions, potentials, gobbler for hydrogen runs
      if (target == 'Ne') call radneon  ! functions, potentials, gobbler for neon runs
        
      call inprintkb
      
C**** SET COEFFICIENTS a,b,c AND INPUT FOR SOLVING DE
      a = -cst2
      c = -cst2
      a(1) = czero
      c(nx-1) = czero
      bb = cone + cst3*v + cst1
      ff = cone - cst3*v - cst1

! llauto has the correct l-value for the autocorrelation function
      q(:,llauto) = dcmplx(g,0.0d0)
      
      q(nx,:) = czero
      q(0,:) = czero
      a1(1) = czero
      c1(nc2) = czero
      b1 = cone

      if(nc.eq.0)    ntfin = -1
      if(ntfin.eq.0) ntfin = -1

      kcount = 0
! printout before we even start
      call output(kcount)

      call system_clock(itime_start,itime)

C*****  INITIAL PROPAGATION FROM DIAGONAL TERM BY tau = dt/2
c$omp parallel do private(jj,j,bet,u,gam)
      do jj=0,nc
        bet=bb(1,jj)
        u(1) = (q(1,jj)*ff(1,jj)+cst2*(q(2,jj)+q(0,jj)))/bet
        do j=2,nx-1
          gam(j)=c(j-1)/bet
          bet = bb(j,jj)-a(j)*gam(j)
          u(j) = (q(j,jj)*ff(j,jj)
     >          +cst2*(q(j+1,jj)+q(j-1,jj))-a(j)*u(j-1))/bet
        end do
        do j=nx-2,1,-1
          u(j) = u(j)-gam(j+1)*u(j+1)
        end do
        do j=1,nx-1
          q(j,jj) = u(j)
        end do
      end do 
c$omp end parallel do
C*****end of propagation by dt/2

!       call output(-1)
!       stop
! reset coefficients to allow for steps of dt by the diagonal term
       a = -cst1
       c = -cst1
       a(1) = czero
       c(nx-1) = czero
       bb = cone + cst*v + cst5
       ff = cone - cst*v - cst5

C  MAIN TIME LOOP
*******************************************************
      do 400 k=0,ntfinpulse
        print *,' time step:', k
        call flush(6)
        kcount = kcount+1

C***** PROPAGATION FROM PULSE BY dt
c$omp parallel do private(j,xj,r1,jj,a1,c1,bet,u1,gam1)
      do 200 j=1,nx-1
        xj = dcmplx(x(j)*f(k),0.d0)
        r1(1) = q(j,0) - cpl(0)*xj*q(j,1)
        r1(nc2) = q(j,nc) - cmn(nc)*xj*q(j,nc1)
        do 120 jj=1,nc1
          r1(jj+1) = q(j,jj)
     >             - xj*(cmn(jj)*q(j,jj-1)+cpl(jj)*q(j,jj+1))
120     continue
        do 130 jj=2,nc2
          a1(jj) = cmn(jj-1)*xj
130     continue
        do 140 jj=1,nc
          c1(jj) = cpl(jj-1)*xj
140     continue

        bet=b1(1)
        u1(1)=r1(1)/bet
        do 150 jj=2,nc2
          gam1(jj)=c1(jj-1)/bet
          bet=b1(jj)-a1(jj)*gam1(jj)
          u1(jj)=(r1(jj)-a1(jj)*u1(jj-1))/bet
150     continue
        do 160 jj=nc,1,-1
          u1(jj)=u1(jj)-gam1(jj+1)*u1(jj+1)
160     continue
        do 170 jj=0,nc
          q(j,jj) = u1(jj+1)
170     continue
200   continue
c$omp end parallel do
C*****end propagation by dt

C***** PROPAGATION FROM DIAGONAL TERM BY dt 
C For the last step it is only tau = dt/2. It should not really
C matter in practice (since the field has died off), but we will
C keep the symmetry for now and jump out of the loop.

      if (k == ntfinpulse) goto 500
c$omp parallel do private(jj,j,bet,u,gam)
      do 300 jj=0,nc
        bet=bb(1,jj)
        u(1)=(q(1,jj)*ff(1,jj) + cst1*(q(2,jj)+q(0,jj)))/bet
        do 310 j=2,nx-1
          gam(j)=c(j-1)/bet
          bet=bb(j,jj)-a(j)*gam(j)
          u(j) = (q(j,jj)*ff(j,jj)
     >           + cst1*(q(j+1,jj)+q(j-1,jj))-a(j)*u(j-1))/bet
310     continue
        do 320 j=nx-2,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
320     continue
        do 330 j=1,nx-1
          q(j,jj) = u(j)
330     continue
300   continue
c$omp end parallel do

C***PRINTOUT FOR TIME LOOP WITH NUMBER nprint*N (N=1,2,...)
      if(kcount.eq.nprint) then
        call output(k)
        kcount = 0
      endif
100   format(10e12.5)
400   continue

! This is just for the last time step of the main loop, where the step is dt/2
500   continue
      a = -cst2
      c = -cst2
      a(1) = czero
      c(nx-1) = czero
      bb = cone + cst3*v + cst1
      ff = cone - cst3*v - cst1  
c$omp parallel do private(jj,j,bet,u,gam)
      do 600 jj=0,nc
        bet=bb(1,jj)
        u(1)=(q(1,jj)*ff(1,jj) + cst2*(q(2,jj)+q(0,jj)))/bet
        do 610 j=2,nx-1
          gam(j)=c(j-1)/bet
          bet=bb(j,jj)-a(j)*gam(j)
          u(j) = (q(j,jj)*ff(j,jj)
     >          + cst2*(q(j+1,jj)+q(j-1,jj))-a(j)*u(j-1))/bet
610     continue
        do 620 j=nx-2,1,-1
          u(j)=u(j)-gam(j+1)*u(j+1)
620     continue
        do 630 j=1,nx-1
          q(j,jj) = u(j)
630     continue
600   continue
c$omp end parallel do

!  print output at this point
      call output(k)
      kcount = 0
      
      call system_clock(itime_end)
      write(*,'(1x,a,f11.3,a)') 'Time elapsed for field propagation: ',
     >      real(itime_end-itime_start)/real(itime), ' seconds.'       
      call flush(6)
      
      if(ntfinpulse.eq.ntfinal) goto 900
C*** CONTINUE WITH 'FREE' PROPAGATION WITH STEP dt AFTER ntfinal
      print *,'continue with free propagation'
      call flush(6)
!  need to reset the arrays again to allow for step dt
      a = -cst1
      c = -cst1
      a(1) = czero
      c(nx-1) = czero
      bb = cone + cst*v + cst5
      ff = cone - cst*v - cst5
             
      call system_clock(itime_start,itime)

      do 800 k = ntfinpulse+1,ntfinal
        print *,' time step:', k
        call flush(6)
        kcount = kcount+1
c$omp parallel do private(jj,j,bet,u,gam)
        do 700 jj=0,nc
          bet=bb(1,jj)
          u(1)=(q(1,jj)*ff(1,jj) + cst1*(q(2,jj)+q(0,jj)))/bet
          do 710 j=2,nx-1
            gam(j)=c(j-1)/bet
            bet=bb(j,jj)-a(j)*gam(j)
            u(j) = (q(j,jj)*ff(j,jj)
     >            + cst1*(q(j+1,jj)+q(j-1,jj))-a(j)*u(j-1))/bet
710       continue
          do 720 j=nx-2,1,-1
            u(j)=u(j)-gam(j+1)*u(j+1)
720       continue
          do 730 j=1,nx-1
            q(j,jj) = u(j)
730       continue
700     continue
c$omp end parallel do
        if(kcount.eq.nprint) then
          call output(k)
          kcount = 0
        endif
800   continue

!     print output at this point
      call output(k)
      
C*****end of propagation

      call system_clock(itime_end)
      write(*,'(1x,a,f11.3,a)') 'Time elapsed for free propagation: ',
     >      real(itime_end-itime_start)/real(itime), ' seconds.'
      call flush(6)
900   continue
! save the wavefunction for a restart 
      do jj = 0,nc
        do i = 0,nx
          if(abs(q(i,jj)).lt.1.0d-30) q(i,jj) = czero
        enddo
        write(60,6000) (q(i,jj),i=0,nx)
      end do
 6000 format(1p10e14.6)
 
      call system_clock(itime_start,itime)

!     testing code output
      call intermediate_output(k)

      call dstrm
      call system_clock(itime_end)
      write(*,'(1x,a,f11.3,a)') 'Time elapsed for dstrm is',
     >      real(itime_end-itime_start)/real(itime), ' seconds.'
      call flush(6)
      
      stop 
      end

C********************************************************
      subroutine output(k)
      use parameters
      use cmplxconstants
      use pgcommon
      use realconstants1
      use realconstants2
      use complexarray3
      use realarray1
      use keys
      use discr
      
      implicit double precision (a-h,o-z)
      character*3 zchar
      character*8 probnm
      character*19 genr
      character*27 flpath
      
      complex*16, allocatable       :: psi(:,:)
      double precision, allocatable :: prob1(:,:),anorm(:),
     >                                 vrlp(:),auto(:),conv(:)
      
      allocate (psi(0:nx+1,0:ntet),prob1(0:nx+1,0:ntet),anorm(nx+1),
     >          auto(nx+1),conv(0:nx+1),vrlp(0:nfmax))
     
      print *,'=== subroutine output: k = ',k 
      call flush(6)
      k1 = k+1

      psi = czero
      prob1 = zero

C*** WAVEFUNCTION (WAVEPACKET)  r * psi (psi = prob1 for saving memory)
      do n=0,nc                
        cek = dsqrt((2.d0*dble(n)+1.d0)/(4.d0*pi))
        do j=0,ntet
          sh = cek * pg(j,n)
          do i=0,ngob1
            psi(i,j) = q(i,n)*dcmplx(sh,0.d0) + psi(i,j)
          end do
        end do
      end do

C*** PROBABILITY DENSITY (r^2*|psi|^2) BY DIRECT SQUARE
      do j=0,ntet
        do i=0,ngob1
          prob1(i,j) = cdabs(psi(i,j))**2 
        end do
      end do

C**** AUTOCORRELATION FUNCTION  
      if(key2 /= 0) then
        do i=1,ngob1
! llauto has the correct l-value for the autocorrelation function
          auto(i) = dreal(q(i-1,llauto))*g(i-1)
        end do
        call arsimd(ngob1,h,auto,are)
        do i=1,ngob1
! llauto has the correct l-value for the autocorrelation function
          auto(i) = dimag(q(i-1,llauto))*g(i-1)
        end do
        call arsimd(ngob1,h,auto,aie)
        vrlp(0) = are**2 + aie**2
        write(10,1000) are,aie,vrlp(0)
      endif
1000  format(/,' Autocorr. function:', 
     >         ' Re =',f11.8,' Im =',f11.8,' Square = ',f11.8)

C**** POPULATION OF DISCRETE STATES FOR key1 neq 0
      if(key1 /= 0) then
        do kd=1,nf
          auto = 0.d0
          do i=1,ngob1                    
            auto(i) = dreal(q(i-1,ll(kd)))*fd(i-1,kd)
          end do
          call arsimd(ngob1,h,auto,are)
          do i=1,ngob1                    
            auto(i) = dimag(q(i-1,ll(kd)))*fd(i-1,kd)
          end do
          call arsimd(ngob1,h,auto,aie)
          vrlp(kd) = are**2 + aie**2
          write(10,1001) kd,nn(kd),ll(kd),are,aie,vrlp(kd)
        end do          
      endif
1001  format(/,' Overlap with state',i2,' : n =',i3,'  l =',i2,
     1 ' Re =',d12.4 ,' Im =',d12.4,'  Square = ',d12.4)

C*** MONITOR CONVERGENCY WITH RESPECT TO PARTIAL WAVES AND NORMALIZATION
      if(key3 /= 0) then
        cnorm = 0.d0
        do j=0,nc
          do i=1,ngob1
            conv(i)=abs(q(i-1,j))**2
          end do
          call arsimd(ngob1,h,conv,awh)
          write(10,1002) j,awh
          cnorm = cnorm + awh
        end do
      endif
1002  format(/,' channel = ',i3,'   contribution = ',e13.4)      

C*** CONTROL OF NORMALIZATION
      if(key3 /= 0) then
        do i=1,ngob1
          anorm(i) = 0.d0
          do j=0,nc
            anorm(i) = anorm(i) + abs(q(i-1,j))**2
          end do
        end do
        call arsimd(ngob1,h,anorm,cnorm)
        write(10,1003) cnorm
      endif
1003  format(/,' norm = ',f11.8)

C*** OPTIONAL PRINTOUT FOR MATHEMATICA
C*** The file names may need to be modified, particularly for non-PC use.
      if (key5.gt.1) then
        iz = iz+1
        rewind(30)
        write(30,9000) iz
        rewind(30)
        read(30,9001) zchar
        genr = 'd:/disk.d/sodium.w/'
        probnm = 'prob.'//zchar
        flpath = genr // probnm
        open(7,file=probnm,status='unknown') 
        do i=0,ngob1
          write(7,9002) (prob1(i,j),j=0,ntet)
        end do
        close(7)
      endif
 9000 format(I3)
 9001 format(A3)
 9002 format(181e11.3)

      call flush(10)
      deallocate(psi,prob1,anorm,auto,conv,vrlp)
      return
      end

C****************************************************
      subroutine readkb
      use realconstants1
      use extra
      use cmplxconstants
      use realconstants2
      use complexarray1
      use complexarray2
      use complexarray3
      use realarray1
      use keys
      use in1
      use discr
      use in2
      use alphas
      use pulseparams
      
      implicit double precision (a-h,o-z)

      namelist /element/  target
      namelist /discrete/ nf,nn,ll,mfixed,nnauto,llauto,nc
      namelist /control/  key1,key2,key3,key4,key5,nprint
      namelist /energies/ emin,de,nerg
      namelist /numerics/ dt,h,nx,gbr,agbr
                  
C----------------------------------------------------------------------------------C
C
C      target    ! character variable for the target element
C
C      nf        ! number of discrete states to project on
C      nn,ll     ! (n,l) values of these states
C      mfixed    ! the m-vlues of the initial state
C      llauto    ! the orbital angular momentum of the initial state
C      nc        ! maximum orbital angular momentum to couple: \ell = 0,1, ... nc
C
C      key1=1    ! calculation of overlap with discrete states
C      key2=1    ! calculation of the autocorrelation function
C      key3=1    ! monitor convergency with respect to partial waves
C      key4=1    ! calculation of photoelectron spectrum and
C      key5=1    ! read energy for angular distribution; 
C                ! if >1, open a file for Mathematica output
C      nprint    ! intermediate results are printed every nprint(+1) time  
C
C      emin      ! lowest energy for which ejected electron spectrum is calculated
C      de        ! difference between energies for which the spectrum is calculated
C      nerg      ! number energies for which the spectrum is calculated
C
C
C      dt        ! time step (in a.u.)
C      h         ! mesh step (in a.u.)     
C      nx        ! number of mesh points (0,1,...,nx), nx must be even
C      gbr       ! gobbler starting radius (au)
C      agbr      ! gobbler strength
C
C      x(0:nx) - radial mesh
C      f(0:nt) - time-dependent laser(s) pulse on the shifted time mesh
C      g(0:nx) - initial state on the radial mesh
C      v(0:nx,0:nc) - local channel potential
C      q(0:nx,0:nc) - solution, 1st index - coordinate, 2nd - channel     
C----------------------------------------------------------------------------------C
C
      read(50,element)
      nf = 6
      nn(1)  = 1
      nn(2)  = 2
      nn(3)  = 2 
      nn(4)  = 3 
      nn(5)  = 3
      nn(6)  = 3
      ll(1)  = 0
      ll(2)  = 0
      ll(3)  = 1
      ll(4)  = 0
      ll(5)  = 1
      ll(6)  = 2
      mfixed = 0
      nnauto = 1
      llauto = 0
      nc     = 8
      read(50,discrete)
      key1   = 1
      key2   = 1
      key3   = 1
      key4   = 1
      key5   = 0
      nprint = 500
      read(50,control)
      dt     = 0.02d0
      h      = 0.02d0
      nx     = 10000
      gbr    = 0.9d0*h*dble(nx)
      agbr   = 5.d0
      read(50,numerics) 
      emin   = 1.d-3
      de     = 1.d-3
      nerg   = 1000
      enelec = -1.0d0
      read(50,energies)
      close(50)

      return
      end
      
C****************************************************
C
C New version of printing routine  (Klaus, Feb. 2014)
C
      subroutine inprintkb
      use parameters
      use pulseparams
      use realarray1
      use realconstants2
      use realconstants1
      use extra
      use complexarray3
      use keys
      use in1
      implicit double precision (a-h,o-z)
      dimension anrm(1:nxmax)
      
      write(40,1001)
1001  format(/,' *********  Laser Parameters:  **********')
      if (abs(alph2*ee2).gt.1.d-10) then
       write(40,1002) ww1,alph1*ee1,period1,ww2,alph2*ee2,period2
      else
       write(40,1102) ww1,alph1*ee1,period1
      endif
1002  format(/,
     >       ' omega1   = ',1p,e11.4,0p,' a.u.',
     >       '   amp1   = ',1p,e11.4,0p,' a.u.',
     >       '     T1   = ',1p,e11.4,0p,' a.u.',/,
     >       ' omega2   = ',1p,e11.4,0p,' a.u.',
     >       '   amp2   = ',1p,e11.4,0p,' a.u.',
     >       '     T2   = ',1p,e11.4,0p,' a.u.')
1102  format(/,
     >       ' omega1   = ',1p,e11.4,0p,' a.u.',
     >       '   amp1   = ',1p,e11.4,0p,' a.u.',
     >       '     T1   = ',1p,e11.4,0p,' a.u.')
      if (abs(alph2*ee2).gt.1.d-10) then 
       write(40,1003) n1up,n1plat,n1down,n2up,n2plat,n2down
      else
       write(40,1103) n1up,n1plat,n1down
      endif
1003  format(/,' n1up, n1plat, n1down:',3i5,
     >       /,' n2up, n2plat, n2down:',3i5)
1103  format(/,' n1up, n1plat, n1down:',3i5)
      if (abs(alph2*ee2).gt.1.d-10) then 
       write(40,1100) shape1up,shape1down,shape2up,shape2down
      else
       write(40,1100) shape1up,shape1down
      endif
1100  format(/' shape for on/off ramp (s = sin^2, g = gaussian,',
     >       ' t = linear (trapezoidal):  ',4a3)  
      if (abs(alph2*ee2).gt.1.d-10) then 
       write(40,1004) CEP1,CEP2
      else
       write(40,1104) CEP1
      endif
1004  format(/,
     >       ' CEP1     = ',1p,e11.4,0p,' deg',/,
     >       ' CEP2     = ',1p,e11.4,0p,' deg')
1104  format(/,' CEP1     = ',1p,e11.4,0p,' deg')
      if (abs(alph2*ee2).gt.1.d-10) then 
       write(40,1005) tstart1,tend1,tstart2,tend2,tpulses,tfinal
      else
       write(40,1105) tstart1,tend1,tpulses,tfinal
      endif

1005  format(/,
     >       ' tstart1  = ',1p,e11.4,0p,' a.u.', 
     >       '   tend1  = ',1p,e11.4,0p,' a.u',/,
     >       ' tstart2  = ',1p,e11.4,0p,' a.u.', 
     >       '   tend2  = ',1p,e11.4,0p,' a.u.',/
     >       ' tpulses  = ',1p,e11.4,0p,' a.u.', 
     >       '  tfinal  = ',1p,e11.4,0p,' a.u.') 
1105  format(/,
     >       ' tstart1  = ',1p,e11.4,0p,' a.u.', 
     >       '   tend1  = ',1p,e11.4,0p,' a.u',/,
     >       ' tpulse   = ',1p,e11.4,0p,' a.u.', 
     >       '  tfinal  = ',1p,e11.4,0p,' a.u.')          
      write(40,1006)
1006  format(//,' *********  Target Parameters:  **********')
      write(40,1007) nnauto,llauto,mfixed
1007  format(/,' initial state:  n = ',i2,',  l = ',i2,',  m = ',i2)
      write(40,1008)
1008  format(/,' *********  Numerical Parameters:  *******') 
      write(40,1009) nx, h, nx*h, ntfinal, dt, tfinal
1009  format(/,' nx    =',i8,';    h  =',1p,e11.4,0p,' a.u.;  rmax = ',
     >           1p,e11.4,0p,' a.u.',//,
     >         ' nt    =',i8,';   dt  =',1p,e11.4,0p,' a.u.;  tfin = ',
     >           1p,e11.4,0p,' a.u.')
      write(40,1010) ngob1,gbr,agbr
1010  format(/,' ngob1 =',i8,';   gbr =',1p,e11.4,0p,' a.u.;  agbr = ',
     >           1p,e11.4,0p,' a.u.')      
      write(40,1011) nc
1011  format(/,' highest L-value of coupled channels (ell = 0,nc):',i4)
      write(40,1012) nprint
1012  format(/,' intermediate printout every',i6,' points',i6,/)
      
      do i=1,ngob1
        anrm(i) = g(i-1)**2
      end do
      call arsimd(ngob1,h,anrm,cn)
      write(40,1013) cn
1013  format(/,' initial norm = ',1p,e16.8,/)
      call flush(40)
      
      return
      end


C****************************************************
      subroutine plg
      use parameters
      use pgcommon
      implicit double precision (a-h,o-z)
      dimension p(0:2*ncmax)
C---- pg = Legendre polynomial Pn   
      dtet = 180.d0/dble(ntet) * dacos(-1.d0)/180.d0
      th = 0.d0
      p(0) = 1.d0
      do 1 i = 0,ntet
        y = dcos(th)
        p(1) = y
        pg(i,0) = 1.d0
        pg(i,1) = y     
        do 2 n=2,2*ncmax
          p(n) = (dble(2*n-1)*y*p(n-1)+dble(1-n)*p(n-2))/dble(n)
          pg(i,n) = p(n)
2       continue
        theta(i) = th
        th = th + dtet
1     continue
      return
      end

************************************************************************
      SUBROUTINE FACTOR
      use factorials
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      FACT(1) = 0.d0
      DO 10 I=1,IFAK-1
       FACT(I+1) = FACT(I) + DLOG(DBLE(I))
10    CONTINUE
      RETURN
      END

************************************************************************
      DOUBLE PRECISION FUNCTION CLB(J1,J2,J3)
      use parameters
      use factorials
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      N = ITRI(J1,J2,J3)
      IF(N.EQ.0) GOTO 10
      L = J1+J2+J3
      M = L/2*2
      IF(M.NE.L) GO TO 10
      CLB = ONE
      B = ONE
      K = (J1+J2-J3)/2
      IF(K.NE.K/2*2) B=-ONE
      AA3 = DBLE(J3)
      KG = (J1+J2+J3)/2
      KG1 = KG+1
      KG2 = 2*KG + 1
      CLB = B * DEXP ( HALF*DLOG(TWO*AA3+ONE) + FACT(KG1) -
     *  FACT(KG1-J1) - FACT(KG1-J2) - FACT(KG1-J3) +
     *  HALF*(FACT(KG2-2*J1) + FACT(KG2-2*J2) + FACT(KG2-2*J3) -
     *  FACT(KG2+1)) )
      RETURN
10    CLB=ZERO
      RETURN
      END

*********************************************
      FUNCTION ITRI(J1,J2,J3)
      ITRI=0
      IF(J1+J2.LT.J3) GOTO 10
      IF(J1+J3.LT.J2) GOTO 10
      IF(J2+J3.LT.J1) GOTO 10
      ITRI=1
10    RETURN
      END

*************************************
      SUBROUTINE ARSIMD(N,DEL,A,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N)
      L=N
      SUM=A(1)-A(L)
      DO 1 I=2,L,2
1     SUM=SUM+4.D0*A(I)+2.D0*A(I+1)
      R=(DEL*SUM)/3.D0
      RETURN
      END

**************************************
      subroutine cnstnt
      use realconstants2
      use cmplxconstants
      use realconstants1
      implicit double precision (a-h,o-z)
      
      pi    = dacos(-1.d0)
      iz    = 100
      nbf   = nx-1
      nx1   = nx+1
      nc1   = nc-1
      nc2   = nc+1
      czero = dcmplx(0.d0,0.d0)
      cone  = dcmplx(1.d0,0.d0)
      ctwo  = dcmplx(2.d0,0.d0)
      cthree= dcmplx(3.d0,0.d0)
      cfour = dcmplx(4.d0,0.d0)
      chalf = dcmplx(.5d0,0.d0)
      ci    = dcmplx(0.d0,1.d0)
      c27   = dcmplx(27.d0,0.d0)
      tau   = dt/2.d0
      zet   = 1.d0
      ch    = dcmplx(h,0.d0)
      cst   = ci*dcmplx(tau,0.d0)
      cst1  = cst/(ctwo*ch*ch)
      cst2  = cst1/ctwo
      cst3  = cst/ctwo
      cst4  = ci*dcmplx(dt,0.d0)
      cst5  = cst/(ch*ch)
      h2    = ctwo*ch
      h3    = cthree*ch
      h24   = cone/(dcmplx(24.d0,0.d0)*ch)
      h1d5  = dcmplx(1.5d0,0.d0)*ch
      hd25  = dcmplx(.25d0,0.d0)*ch

      return
      end      

*************************************
      subroutine arrays
      use parameters
      use realconstants1
      use realconstants2
      use cmplxconstants
      use realconstants1
      use complexarray1
      use complexarray2
      use complexarray3
      use realarray1
      implicit double precision (a-h,o-z)
      
      allocate (x(0:nx+1),g(0:nx+1),fd(0:nx+1,nfmax),s(0:2*nc))
      allocate (cek3(0:nc,0:nc,0:2*nc))
      allocate (a(nx+1),b(nx+1),c(nx+1),r(nx+1),gam(nx+1),u(nx+1))
      allocate (cmn(0:nc),cpl(0:nc))
      allocate (a1(nc+1),b1(nc+1),c1(nc+1),r1(nc+1),gam1(nc+1),u1(nc+1))
      allocate (q(0:nx+1,0:nc),v(0:nx+1,0:nc),ff(0:nx+1,0:nc))
      allocate (bb(0:nx+1,0:nc))

      call LFAK
    
C  SET INITIAL ZEROS
      a = czero
      b = czero
      c = czero
      r = czero
      gam = czero
      u = czero

      x = zero
      g = zero
      
      fd = zero
      
      q  = czero
      v  = czero
      ff = czero
      bb = czero
      
      a1 = czero
      b1 = czero
      c1 = czero
      r1 = czero
      gam1 = czero
      u1 = czero

C*** radial mesh
      do 100 i=1,nx           
        x(i)= h * dble(i)
100   continue

C*** find ngob
      ngob = nx
      if(gbr.lt.1.d-10) goto 200
      do 201 i = 0,nx
        if(abs(x(i)-gbr).lt.1.d-10) then
          ngob = i
          goto 200
        endif
201   continue
200   continue
      if(ngob.eq.0) ngob = nx
      ngob1 = ngob
      if(mod(ngob,2).eq.0) ngob1 = ngob-1  ! makes ngob1 odd

C ANGLE MESH, LEGENDRE POLYNOMIALS, LOG OF FACTORIALS
      call plg
      call factor

C  SET MIXING COEFFICIENTS (multiplied by i*tau)
      cmn(0) = czero
      cpl(0) = czero
      am = dble(mfixed)
      if(mfixed.eq.0) cpl(0) = cone/cdsqrt(cthree) * cst
      do 300 i=1,nc
        nn = i
        ckr = dble(nn)
        ckrlow = dble(nn-1)
        ckrhgh = dble(nn+1)
        intckr = 2*nint(ckr)
        intckrhgh = 2*nint(ckrhgh)
        intckrlow = 2*nint(ckrlow)
        call clegor(intckr,intckrhgh,2,-2*mfixed,2*mfixed,0,cghigh,ier)
        call clegor(intckr,intckrlow,2,-2*mfixed,2*mfixed,0,cglow,ier)
        apl = dsqrt((2.d0*ckr+1.d0)*(2.d0*ckr+3.d0))*clb(nn,nn+1,1)
     >        *cghigh/3.0d0
        amn = dsqrt((2.d0*ckr+1.d0)*(2.d0*ckr-1.d0))*clb(nn,nn-1,1)
     >       *cglow/3.d0
        
        cpl(i) = dcmplx(apl,0.d0)*cst 
        cmn(i) = dcmplx(amn,0.d0)*cst 
300   continue
      cpl(nc) = czero

C**** GEOMETRICAL PART OF MULTIPOLE COEFFICIENTS AND nbmax
      do 400 n1 = 0,nc                 
        cek1 = dsqrt(2.d0*dble(n1)+1.d0)
        do 401 n2 = 0,nc
          cek2 = dsqrt(2.d0*dble(n2)+1.d0)*cek1
          do 402 nbig = iabs(n1-n2), n1+n2, 2
            if(nbig.gt.nbmax) nbmax = nbig
            cek3(n1,n2,nbig) = cek2*clb(n1,n2,nbig)**2/(4.d0*pi)
402        continue
401     continue
400   continue

      return
      end

C********************************************************************
C     written by John Emmons and Sean Buczek 
C     date Feb 2014
C
C**** print data at intermediate time steps (i.e. not just the end
C**** the columns of the output are listed beloew with | separating each entry:
C**** OUTPUT FILE: intermediate.out (file number 230)
C
C**** timestep \t overlap in 1s \t overlap in 2s \t ... \t integral n=1 \t integral n=2 \t ... \t integral betas \t sum of integral
C 

C
C     TODO
C     [X] make output acutally write to file 230
C     [ ] make a format statement
C     [ ] compute betas integral
C     [ ] sum the overaps for vairous n
C     [ ] sum the sum of overlap and integral of betas
C     [ ] print timestep (add input parameter)
C     [ ] make code print at nprint intervals

      subroutine intermediate_output (k)
      use parameters
      use in2
      use discr
      use complexarray3
      use realarray1
      use realconstants2
      use realconstants1
      use adis
      use VDD
      use pecommon
      use rpcommon
      use DAT1
      use phcommon
      use keys
! modules from the output subroutine
      use cmplxconstants
      use pgcommon
      use realarray1

! dstrm variables
      implicit double precision (a-h,o-z)
      complex*16 tint,be,ct1

      double precision, allocatable :: dcrall(:),phall(:,:),betall(:,:)

      dimension tint(0:ncmax,0:nenmax),bet(0:2*ncmax)
      dimension dm(0:nenmax),ener(0:nenmax),wr(nxmax),wi(nxmax) 
      dimension rtp(nxmax+1000),rpe(nxmax+1000)
! output variables
      character*3 zchar
      character*8 probnm
      character*19 genr
      character*27 flpath

      complex*16, allocatable       :: psi(:,:)
      double precision, allocatable :: prob1(:,:),anorm(:),
     >                                 vrlp(:),auto(:),conv(:)

! dstrm memory allocation
      allocate(dcrall(nerg),phall(nerg,0:nc),betall(nerg,0:2*nc))

! output memory allocation
      allocate (psi(0:nx+1,0:ntet),prob1(0:nx+1,0:ntet),anorm(nx+1),
     >          auto(nx+1),conv(0:nx+1),vrlp(0:nfmax))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!dstrm code!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ncu = nc
      nrp1 = nrp+1

      dm = 0.d0
      ener(0) = emin-de
      rpe(1) = 0.d0
      rtp(1) = 0.d0
      
      do 990 i = 2,nrp1
        rtp(i) = rp(i-1)
990       continue        

      do nen = 1,nerg
        ener(nen) = ener(nen-1) + de
      end do

c$omp parallel do private(nen,ep,jj,phse,dw,i,rpe,rg,wr,wi,rtr,rti)
      do 390 nen = 1,nerg
        ep = 2.0d0*ener(nen)
        do 290 jj=0,nc
          call dwavekb(jj,ep,phse,dw)         
          phall(nen,jj) = phse
          do 119 i=2,nrp1
            rpe(i)=dw(i-1)              
119               continue
          do 129 i=1,ngob
            call parinv(x(i),rtp,rpe,nrp1,rg)
            wr(i+1) = rg*dreal(q(i,jj))
            wi(i+1) = rg*dimag(q(i,jj))
129               continue
          call arsimd(ngob1,h,wr,rtr)
          call arsimd(ngob1,h,wi,rti)
          tint(jj,nen) = dcmplx(rtr,rti)
c---- this is for angular distribution:
          if( dabs(ener(nen)-enelec).le.1.d-6 ) cdst(jj)=tint(jj,nen)
c----------------------------------------------------------------------
          dm(nen) = dm(nen) + cdabs(tint(jj,nen))**2
290           continue
9901            format(f9.5,3d15.5)
390                continue
c$omp end parallel do
        
c$omp parallel do private(nen,lam,be,rlam,j1,j2,rj1,rj2,
c$omp>                    cgdd,ier,aw,rph,ct1,bet,dcr)
c---- beta parameters
      do 490 nen = 1,nerg
        am = dble(mfixed)
        if (key4.gt.1) then
         do 259 lam=0,min(2*nc,20)
          be = dcmplx(0.d0,0.d0)
          rlam = dble(lam)
          do 249 j1 = 0,nc
            j2min = abs(j1-lam)
            j2max = min(nc,j1+lam)
            do 239 j2 = j2min,j2max,2
              rj1 = dble(j1)
              rj2 = dble(j2)       
              call clegor(2*j1,2*lam,2*j2,2*mfixed,0,2*mfixed,cgdd,ier)
              aw = dsqrt((2.d0*rj1+1.d0)/(2.d0*rj2+1.d0))
     >            *clb(j1,lam,j2)*cgdd
              rph = phall(nen,j1)-phall(nen,j2)
              ct1 = (0.d0,1.d0)**(j2-j1) * cdexp(dcmplx(0.d0,rph)) 
     >            *tint(j1,nen)*dconjg(tint(j2,nen))*dcmplx(aw,0.d0)
              be = be + ct1
239                   continue
249                          continue
          bet(lam) = dreal(be)*(2.d0*rlam+1.d0)/dm(nen)  
          betall(nen,lam) = bet(lam)
259            continue 
        endif
        dcr = dm(nen)*dsqrt(2.d0/ener(nen))
        dcrall(nen) = dcr
490       continue

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!OUTPUT CODE!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      k1 = k+1

      psi = czero
      prob1 = zero

C*** WAVEFUNCTION (WAVEPACKET)  r * psi (psi = prob1 for saving memory)
      do n=0,nc                
        cek = dsqrt((2.d0*dble(n)+1.d0)/(4.d0*pi))
        do j=0,ntet
          sh = cek * pg(j,n)
          do i=0,ngob1
            psi(i,j) = q(i,n)*dcmplx(sh,0.d0) + psi(i,j)
          end do
        end do
      end do

C*** PROBABILITY DENSITY (r^2*|psi|^2) BY DIRECT SQUARE
      do j=0,ntet
        do i=0,ngob1
          prob1(i,j) = cdabs(psi(i,j))**2 
        end do
      end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!BEGIN OUTPUTING!!!!!!!!!!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! orbital probabilities
C**** POPULATION OF DISCRETE STATES FOR key1 neq 0
!      if(key1 /= 0) then
        do kd=1,nf
          auto = 0.d0
          do i=1,ngob1                    
            auto(i) = dreal(q(i-1,ll(kd)))*fd(i-1,kd)
          end do
          call arsimd(ngob1,h,auto,are)
          do i=1,ngob1                    
            auto(i) = dimag(q(i-1,ll(kd)))*fd(i-1,kd)
          end do
          call arsimd(ngob1,h,auto,aie)
          vrlp(kd) = are**2 + aie**2
          write(230,1001) kd,nn(kd),ll(kd),are,aie,vrlp(kd)
        end do          
!      endif
 1001   format(/,' Overlap with state',i2,' : n =',i3,'  l =',i2,
     1 ' Re =',d12.4 ,' Im =',d12.4,'  Square = ',d12.4)   

! sean's beta integral garbage
      do nen=1, nerg, 1
         if (nen < nerg) then
            result = result + (sqrt(2.d0*ener(nen))*dcrall(nen))
         else
            result =result+(.5*(sqrt(2.d0*ener(nen))*dcrall(nen)))
         endif
      end do

C tie up loose ends of integration
      result = result + ((2*(sqrt(2.d0*ener(1))*dcrall(1)))-
     >    (2.d0*ener(2)*dcrall(2)))
      result= result*ener(1)

      write(230,*) 'Integration is equal to', result

c$omp end parallel do
      do nen = 1,nerg

      end do

      deallocate(psi,prob1,anorm,auto,conv,vrlp)

! test output 
     
      flush(230)
      close(230)

      return
      end  

C*********************************************
C***energy and angular distribution in the final state
      subroutine dstrm
      use parameters
      use in2
      use complexarray3
      use realarray1
      use realconstants2
      use realconstants1
      use adis
      use VDD
      use pecommon
      use rpcommon
      use DAT1
      use phcommon
      use keys

      implicit double precision (a-h,o-z)
      complex*16 tint,be,ct1

      double precision, allocatable :: dcrall(:),phall(:,:),betall(:,:)

      dimension tint(0:ncmax,0:nenmax),bet(0:2*ncmax)
      dimension dm(0:nenmax),ener(0:nenmax),wr(nxmax),wi(nxmax) 
      dimension rtp(nxmax+1000),rpe(nxmax+1000)

      allocate(dcrall(nerg),phall(nerg,0:nc),betall(nerg,0:2*nc))

      ncu = nc
      nrp1 = nrp+1

      dm = 0.d0
      ener(0) = emin-de
      rpe(1) = 0.d0
      rtp(1) = 0.d0
      
      do 100 i = 2,nrp1
        rtp(i) = rp(i-1)
100   continue        

      write(10,1000)
1000  format(/,'#  E(au)',6x,'Spectr    bet0    bet1     bet2     bet3',
     >         '     bet4     bet5     bet6     bet7     bet8')
      call flush(10)
      
      do nen = 1,nerg
        ener(nen) = ener(nen-1) + de
      end do

c$omp parallel do private(nen,ep,jj,phse,dw,i,rpe,rg,wr,wi,rtr,rti)
      do 300 nen = 1,nerg
        ep = 2.0d0*ener(nen)
        do 200 jj=0,nc
          call dwavekb(jj,ep,phse,dw)         
          phall(nen,jj) = phse
          do 110 i=2,nrp1
            rpe(i)=dw(i-1)              
110       continue
          do 120 i=1,ngob
            call parinv(x(i),rtp,rpe,nrp1,rg)
            wr(i+1) = rg*dreal(q(i,jj))
            wi(i+1) = rg*dimag(q(i,jj))
120       continue
          call arsimd(ngob1,h,wr,rtr)
          call arsimd(ngob1,h,wi,rti)
          tint(jj,nen) = dcmplx(rtr,rti)
c---- this is for angular distribution:
          if( dabs(ener(nen)-enelec).le.1.d-6 ) cdst(jj)=tint(jj,nen)
c----------------------------------------------------------------------
          dm(nen) = dm(nen) + cdabs(tint(jj,nen))**2
200     continue
1001  format(f9.5,3d15.5)
300   continue
c$omp end parallel do
        
c$omp parallel do private(nen,lam,be,rlam,j1,j2,rj1,rj2,
c$omp>                    cgdd,ier,aw,rph,ct1,bet,dcr)
c---- beta parameters
      do 400 nen = 1,nerg
        am = dble(mfixed)
        if (key4.gt.1) then
         do 250 lam=0,min(2*nc,20)
          be = dcmplx(0.d0,0.d0)
          rlam = dble(lam)
          do 240 j1 = 0,nc
            j2min = abs(j1-lam)
            j2max = min(nc,j1+lam)
            do 230 j2 = j2min,j2max,2
              rj1 = dble(j1)
              rj2 = dble(j2)       
              call clegor(2*j1,2*lam,2*j2,2*mfixed,0,2*mfixed,cgdd,ier)
              aw = dsqrt((2.d0*rj1+1.d0)/(2.d0*rj2+1.d0))
     >            *clb(j1,lam,j2)*cgdd
              rph = phall(nen,j1)-phall(nen,j2)
              ct1 = (0.d0,1.d0)**(j2-j1) * cdexp(dcmplx(0.d0,rph)) 
     >            *tint(j1,nen)*dconjg(tint(j2,nen))*dcmplx(aw,0.d0)
              be = be + ct1
230         continue
240       continue
          bet(lam) = dreal(be)*(2.d0*rlam+1.d0)/dm(nen)  
          betall(nen,lam) = bet(lam)
250      continue 
        endif
        dcr = dm(nen)*dsqrt(2.d0/ener(nen))
        dcrall(nen) = dcr
400   continue
c$omp end parallel do
      do nen = 1,nerg
        do jj = 0,nc
          write(20,1001) ener(nen),phall(nen,jj),tint(jj,nen)
        end do
        call flush(20)
        if (key4.gt.1) then
          write(10,1002) ener(nen),sqrt(2.d0*ener(nen))*dcrall(nen),
     >                   (betall(nen,lam),lam=0,min(8,2*nc))
        else
           write(10,1002) ener(nen),sqrt(2.d0*ener(nen))*dcrall(nen)
        endif
        call flush(10)
        if (key4.gt.1) then
          write(70,1003) ener(nen),sqrt(2.d0*ener(nen))*dcrall(nen),
     >                   (betall(nen,lam),lam=0,min(2*nc,20))
        else
          write(70,1003) ener(nen),sqrt(2.d0*ener(nen))*dcrall(nen)
        endif
        call flush(70)
      end do
1002  format(f9.5,e10.3,15f9.5)
1003  format(1p200e14.6)
      return
      end  

C***********************************************************
c needs clean-up, but works
      DOUBLE PRECISION FUNCTION FACOUZ(E,L,Z)
      implicit double precision (a-h,o-z)
      GAM=-Z/DSQRT(E)
      M=201-L
      DO 19 K=1,M
      AL=DATAN(GAM/(202-K))
      IF(K.EQ.1) GO TO 18
      FACOUZ=FACOUZ-AL
      GO TO 19
18    BE=DSQRT(GAM*GAM+(202-K)**2)
      FACOUZ=AL*200.5d0+GAM*(DLOG(BE)-1.d0)
     *+(-SIN(AL)/12.0d0+SIN(3.d0*AL)/(360.d0*BE*BE))/BE
19    CONTINUE
      RETURN
      END
C***********************************************************
c needs clean-up, but works
      SUBROUTINE PARINV(X,A,F,N,R)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION A(N),F(N)
      IF(X.LT.A(1)) GO TO 11
      IF(X.GT.A(N)) GO TO 4
      K1=1
      K2=N
2     K3=K2-K1
      IF(K3.LE.1) GO TO 6
      K3=K1+K3/2
      IF(A(K3)-X) 7,8,9
7     K1=K3
      GO TO 2
9     K2=K3
      GO TO 2
8     R=F(K3)
      RETURN
3     B1=A(K1)
      B2=A(K1+1)
      B3=A(K1+2)
      B4=F(K1)
      B5=F(K1+1)
      B6=F(K1+2)
      R=B4*((X-B2)*(X-B3))/((B1-B2)*(B1-B3))+B5*((X-B1)*(X-B3))/
     *((B2-B1)*(B2-B3))+B6*((X-B1)*(X-B2))/((B3-B1)*(B3-B2))
      RETURN
6     IF(K2.NE.N) GO TO 3
      K1=N-2
      GO TO 3
4     C=DABS(X-A(N))
      IF(C.LT.0.1D-7) GO TO 5
      K1=N-2
13    WRITE(*,41) X
41    FORMAT(25H X IS OUT OF THE INTERVAL,3H X=,F13.7)
      GO TO 3
5     R=F(N)
      RETURN
11    C=DABS(X-A(1))
      IF(C.LT.0.1D-7) GO TO 12
      K1=1
      GO TO 13
12    R=F(1)
      RETURN
      END

C----------------------------------------------------------------------
      SUBROUTINE COUL90(X, ETA, XLMIN,LRANGE, FC,GC,FCP,GCP, KFN,IFAIL)
C----------------------------------------------------------------------
C
C  COULOMB & BESSEL FUNCTION PROGRAM-- COUL90 -- USING STEED'S METHOD
C
C  COUL90 RETURNS ARRAYS FC = F, GC = G, FCP = (D/DX) F, GCP = (D/DX) G
C   FOR REAL X .GT. 0. ,REAL ETA (INCLUDING 0.), AND REAL XLMIN .GT.-1.
C   FOR (LRANGE+1) INTEGER-SPACED LAMBDA VALUES.
C   IT HENCE GIVES POSITIVE-ENERGY SOLUTIONS TO THE COULOMB SCHRODINGER
C   EQUATION, TO THE KLEIN-GORDON EQUATION AND TO SUITABLE FORMS OF
C   THE DIRAC EQUATION.    BY SETTING ETA = 0.0 AND RENORMALISING
C   SPHERICAL & CYLINDRICAL BESSEL FUNCTIONS ARE COMPUTED INSTEAD.
C----------------------------------------------------------------------
C   CALLING VARIABLES; ALL REALS ARE DOUBLE PRECISION (REAL*8)
C
C   X       - REAL ARGUMENT FOR COULOMB FUNCTIONS > 0.0 
C             [ X > SQRT(ACCUR) : ACCUR IS TARGET ACCURACY 1.0D-14 ]
C   ETA     - REAL SOMMERFELD PARAMETER, UNRESTRICTED > = < 0.0
C   XLMIN   - REAL MINIMUM LAMBDA-VALUE (L-VALUE OR ORDER),
C             GENERALLY IN RANGE 0.0 - 1.0 AND MOST USUALLY 0.0
C   LRANGE  - INTEGER NUMBER OF ADDITIONAL L-VALUES : RESULTS RETURNED
C             FOR L-VALUES XLMIN TO XLMIN + LRANGE INCLUSIVE
C   FC ,GC  - REAL VECTORS F,G OF REGULAR, IRREGULAR COULOMB FUNCTIONS
C   FCP,GCP - REAL VECTORS FOR THE X-DERIVATIVES OF  F,G
C             THESE VECTORS TO BE OF LENGTH AT LEAST MINL + LRANGE
C             STARTING ELEMENT MINL = MAX0( IDINT(XLMIN+ACCUR),0 )
C   KFN     - INTEGER CHOICE OF FUNCTIONS TO BE COMPUTED :
C           = 0         REAL COULOMB FUNCTIONS AND DERIVATIVES F & G
C           = 1    SPHERICAL BESSEL      "      "     "        j & y
C           = 2  CYLINDRICAL BESSEL      "      "     "        J & Y
C
C   PRECISION:  RESULTS TO WITHIN 2-3 DECIMALS OF "MACHINE ACCURACY"
C   IN OSCILLATING REGION X .GE. [ETA + SQRT{ETA**2 + XLM*(XLM+1)}]
C   I.E. THE TARGET ACCURACY ACCUR SHOULD BE 100 * ACC8 WHERE ACC8 IS
C   THE SMALLEST NUMBER WITH 1.+ACC8.NE.1. FOR OUR WORKING PRECISION.
C   THIS RANGES BETWEEN 4E-15 AND 2D-17 ON CRAY, VAX, SUN, PC FORTRANS
C   SO CHOOSE A SENSIBLE  ACCUR = 1.0D-14
C   IF X IS SMALLER THAN [ ] ABOVE THE ACCURACY BECOMES STEADILY WORSE:
C   THE VARIABLE PACCQ IN COMMON /STEED/ HAS AN ESTIMATE OF ACCURACY.
C----------------------------------------------------------------------
C   ERROR RETURNS                THE USER SHOULD TEST IFAIL ON EXIT
C
C   IFAIL ON INPUT IS SET TO 0                        LIMIT = 20000
C   IFAIL IN OUTPUT =  0 : CALCULATIONS SATISFACTORY
C                   =  1 : CF1 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C                   =  2 : CF2 DID NOT CONVERGE AFTER LIMIT ITERATIONS
C                   = -1 : X < 1D-7 = SQRT(ACCUR)
C                   = -2 : INCONSISTENCY IN ORDER VALUES (L-VALUES) 
C----------------------------------------------------------------------
C  MACHINE-DEPENDENT PARAMETERS:    ACCUR - SEE ABOVE
C           SMALL - OFFSET FOR RECURSION = APPROX SQRT(MIN REAL NO.)
C           IE 1D-30 FOR IBM REAL*8,    1D-150 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
C  PROGRAMMING HISTORY AND BACKGROUND: CPC IS COMPUTER PHYSICS COMMUN.
C  ORIGINAL PROGRAM  RCWFN       IN    CPC  8 (1974) 377-395
C                 +  RCWFF       IN    CPC 11 (1976) 141-142
C  FULL DESCRIPTION OF ALGORITHM IN    CPC 21 (1981) 297-314
C  REVISED STANDARD  COULFG      IN    CPC 27 (1982) 147-166
C  BACKGROUND MATERIAL IN J. COMP. PHYSICS 46 (1982) 171-188         
C  CURRENT PROGRAM   COUL90  (FORTRAN77) SUPERCEDES THESE EARLIER ONES
C  (WHICH WERE WRITTEN IN FORTRAN 4) AND ALSO BY INCORPORATING THE NEW
C  LENTZ-THOMPSON ALGORITHM FOR EVALUATING THE FIRST CONTINUED FRACTION
C  ..SEE ALSO APPENDIX TO J. COMP. PHYSICS 64 (1986) 490-509     1.4.94
C----------------------------------------------------------------------
C  AUTHOR: A. R. BARNETT           MANCHESTER  MARCH   1981/95
C                                  AUCKLAND    MARCH   1991
C----------------------------------------------------------------------
      IMPLICIT         NONE
      INTEGER          LRANGE, KFN, IFAIL
      DOUBLE PRECISION X, ETA, XLMIN
      DOUBLE PRECISION FC (0:*), GC (0:*), FCP(0:*), GCP(0:*)
C----- ARRAYS INDEXED FROM 0 INSIDE SUBROUTINE: STORAGE FROM MINL
      DOUBLE PRECISION ACCUR,ACCH,SMALL, ONE,ZERO,HALF,TWO,TEN2, RT2DPI
      DOUBLE PRECISION XINV,PK,CF1,C,D,PK1,ETAK,RK2,TK,DCF1,DEN,XLM,XLL
      DOUBLE PRECISION EL,XL,RL,SL, F,FCMAXL,FCMINL,GCMINL,OMEGA,WRONSK
      DOUBLE PRECISION WI, A,B, AR,AI,BR,BI,DR,DI,DP,DQ, ALPHA,BETA
      DOUBLE PRECISION E2MM1, FJWKB,GJWKB, P,Q,PACCQ, GAMMA,GAMMAI
      INTEGER          IEXP, NFP, NPQ, L, MINL,MAXL, LIMIT
      LOGICAL          ETANE0, XLTURN
      PARAMETER      ( LIMIT = 30000, SMALL = 1.0D-150 )
      COMMON  /STEED/  PACCQ,NFP,NPQ,IEXP,MINL    !not required in code
      COMMON  /DESET/  CF1,P,Q,F,GAMMA,WRONSK     !information only
C----------------------------------------------------------------------
C     COUL90 HAS CALLS TO: DSQRT,DABS,MAX0,IDINT,DSIGN,DFLOAT,DMIN1
C----------------------------------------------------------------------
      DATA ZERO,ONE,TWO,TEN2,HALF /0.0D0, 1.0D0, 2.0D0, 1.0D2, 0.5D0/
      DATA RT2DPI /0.79788 45608 02865  D0/ 
CQ    DATA RT2DPI /0.79788 45608 02865 35587 98921 19868 76373 Q0/
C-----THIS CONSTANT IS  DSQRT(TWO / PI):
C-----USE Q0 FOR IBM REAL*16: D0 FOR REAL*8 AND DOUBLE PRECISION
C----------------CHANGE ACCUR TO SUIT MACHINE AND PRECISION REQUIRED
      ACCUR = 1.0D-14
      IFAIL = 0
      IEXP  = 1
      NPQ   = 0
      GJWKB = ZERO
      PACCQ = ONE
      IF(KFN .NE. 0) ETA = ZERO
                 ETANE0  = ETA .NE. ZERO                                
      ACCH  = DSQRT(ACCUR)
C-----   TEST RANGE OF X, EXIT IF.LE.DSQRT(ACCUR) OR IF NEGATIVE
      IF( X .LE. ACCH )                GO TO 100
      IF( KFN.EQ.2 )   THEN
         XLM = XLMIN - HALF                                  
        ELSE
         XLM = XLMIN                                                     
        ENDIF
      IF( XLM.LE.-ONE .OR. LRANGE.LT.0 )         GO TO 105 
      E2MM1  = XLM * XLM + XLM
      XLTURN = X * (X -  TWO * ETA) .LT. E2MM1
      E2MM1  = E2MM1  +  ETA * ETA
      XLL    = XLM + DFLOAT(LRANGE)
C-----  LRANGE IS NUMBER OF ADDITIONAL LAMBDA VALUES TO BE COMPUTED
C-----  XLL  IS MAX LAMBDA VALUE [ OR 0.5 SMALLER FOR J,Y BESSELS ]
C-----  DETERMINE STARTING ARRAY ELEMENT (MINL) FROM XLMIN
      MINL  = MAX0( IDINT(XLMIN + ACCUR),0 )     ! index from 0
      MAXL  = MINL + LRANGE
C-----   EVALUATE CF1  =  F   =  DF(L,ETA,X)/DX   /   F(L,ETA,X)
      XINV = ONE / X
      DEN  = ONE                       ! unnormalised F(MAXL,ETA,X)
      PK   = XLL + ONE
      CF1  = ETA / PK  +  PK * XINV                                             
           IF( DABS(CF1).LT.SMALL )    CF1 = SMALL
      RK2  = ONE
         D = ZERO
         C = CF1
C----- BEGIN CF1 LOOP ON PK = K STARTING AT LAMBDA + 1: LENTZ-THOMPSON
      DO 10 L =  1 , LIMIT             ! abort if reach LIMIT (20000)    
          PK1 = PK + ONE
          IF( ETANE0 ) THEN
                ETAK = ETA / PK
                RK2  = ONE + ETAK * ETAK
                 TK  = (PK + PK1) * (XINV + ETAK / PK1)
             ELSE
                 TK  = (PK + PK1) * XINV
             ENDIF
          D   =  TK - RK2 * D          ! direct  ratio of B convergents    
          C   =  TK - RK2 / C          ! inverse ratio of A convergents
            IF( DABS(C).LT.SMALL ) C = SMALL
            IF( DABS(D).LT.SMALL ) D = SMALL
          D   = ONE / D
          DCF1=   D * C
          CF1 = CF1 * DCF1
              IF( D.LT.ZERO )    DEN = -DEN
          PK  = PK1
          IF( DABS(DCF1-ONE).LT.ACCUR )     GO TO  20 ! proper exit
   10 CONTINUE
                                            GO TO 110 ! error exit 
   20       NFP = PK - XLL - 1                        ! number of steps
              F = CF1                                 ! need DEN later
C----DOWNWARD RECURRENCE TO LAMBDA = XLM; ARRAYS GC, GCP STORE RL, SL
      IF( LRANGE.GT.0 )       THEN
          FCMAXL    = SMALL  * DEN 
          FCP(MAXL) = FCMAXL * CF1
          FC (MAXL) = FCMAXL
                    XL = XLL                   
                    RL = ONE
          DO 30 L =  MAXL, MINL+1, -1
             IF( ETANE0 )  THEN
                    EL = ETA / XL                
                    RL = DSQRT( ONE + EL * EL )
                    SL = XL * XINV  + EL
                    GC (L) = RL                  ! storage
                    GCP(L) = SL
                ELSE
                    SL = XL * XINV
                ENDIF
             FC (L-1)  = ( FC(L)   * SL  +  FCP(L) ) / RL
             FCP(L-1)  =   FC(L-1) * SL  -  FC (L) * RL
             XL    =  XL - ONE                   ! end value is XLM
   30     CONTINUE
         IF( DABS(FC(MINL)).LT.ACCUR*SMALL )  FC(MINL) = ACCUR * SMALL
          F   = FCP(MINL) / FC(MINL)             ! F'/F at min L-value
          DEN = FC (MINL)                        ! normalisation
      ENDIF
C---------------------------------------------------------------------
C-----   NOW WE HAVE REACHED LAMBDA = XLMIN = XLM
C-----   EVALUATE CF2 = P + I.Q  USING STEED'S ALGORITHM (NO ZEROS)
C---------------------------------------------------------------------
      IF( XLTURN ) CALL JWKB( X,ETA,DMAX1(XLM,ZERO),FJWKB,GJWKB,IEXP )
      IF( IEXP.GT.1 .OR. GJWKB.GT.(ONE / (ACCH*TEN2)) ) THEN
          OMEGA = FJWKB
          GAMMA = GJWKB * OMEGA
          P     = F
          Q     = ONE
        ELSE                                     ! find cf2     
          XLTURN = .FALSE.
          PK =  ZERO
          WI =  ETA + ETA
          P  =  ZERO
          Q  =  ONE - ETA * XINV
          AR = -E2MM1
          AI =  ETA
          BR =  TWO * (X - ETA)
          BI =  TWO
          DR =  BR / (BR * BR + BI * BI)
          DI = -BI / (BR * BR + BI * BI)
          DP = -XINV * (AR * DI + AI * DR)
          DQ =  XINV * (AR * DR - AI * DI)
          DO 40 L = 1, LIMIT
             P  = P  + DP
             Q  = Q  + DQ
             PK = PK + TWO
             AR = AR + PK
             AI = AI + WI                                                   
             BI = BI + TWO                                                  
             D  = AR * DR - AI * DI + BR                                        
             DI = AI * DR + AR * DI + BI                                        
             C  = ONE / (D * D + DI * DI)                  
             DR =  C * D                                                      
             DI = -C * DI                                                     
             A  = BR * DR - BI * DI - ONE                                       
             B  = BI * DR + BR * DI                                             
             C  = DP * A  - DQ * B
             DQ = DP * B  + DQ * A                                              
             DP = C
      IF( DABS(DP)+DABS(DQ).LT.(DABS(P)+DABS(Q)) * ACCUR ) GO TO 50
   40     CONTINUE
                                              GO TO 120 ! error exit
   50     NPQ   = PK / TWO                              ! proper exit
          PACCQ = HALF * ACCUR / DMIN1( DABS(Q),ONE )
          IF( DABS(P).GT.DABS(Q) ) PACCQ = PACCQ * DABS(P)
C---------------------------------------------------------------------
C    SOLVE FOR FCMINL = F AT LAMBDA = XLM AND NORMALISING FACTOR OMEGA
C---------------------------------------------------------------------
          GAMMA   = (F - P) / Q
          GAMMAI  = ONE / GAMMA
          IF( DABS(GAMMA) .LE. ONE )  THEN 
                 OMEGA  = DSQRT( ONE  +  GAMMA * GAMMA )
            ELSE
                 OMEGA  = DSQRT( ONE  +  GAMMAI* GAMMAI) * DABS(GAMMA)
            ENDIF 
          OMEGA  = ONE / ( OMEGA * DSQRT(Q) )
          WRONSK = OMEGA
        ENDIF   
C--------------------------------------------------------------------- 
C    RENORMALISE IF SPHERICAL OR CYLINDRICAL BESSEL FUNCTIONS
C---------------------------------------------------------------------
      IF( KFN.EQ.1 )       THEN         !   spherical Bessel functions
                 ALPHA = XINV
                 BETA  = XINV
        ELSEIF( KFN.EQ.2 ) THEN         ! cylindrical Bessel functions
                 ALPHA = HALF * XINV
                 BETA  = DSQRT( XINV ) * RT2DPI
        ELSE                            ! kfn = 0,   Coulomb functions
                 ALPHA = ZERO     
                 BETA  = ONE
        ENDIF
      FCMINL = DSIGN( OMEGA,DEN ) * BETA
      IF( XLTURN )   THEN
                        GCMINL =   GJWKB * BETA
        ELSE
                        GCMINL =  FCMINL * GAMMA
        ENDIF
      IF( KFN.NE.0 )    GCMINL = -GCMINL         ! Bessel sign differs
      FC (MINL) = FCMINL
      GC (MINL) = GCMINL
      GCP(MINL) = GCMINL * (P - Q * GAMMAI - ALPHA) 
      FCP(MINL) = FCMINL * (F - ALPHA)
      IF( LRANGE.EQ.0 )                          RETURN
C---------------------------------------------------------------------
C    UPWARD RECURRENCE FROM GC(MINL),GCP(MINL) STORED VALUES ARE RL,SL
C    RENORMALISE FC,FCP AT EACH LAMBDA AND CORRECT REGULAR DERIVATIVE
C      XL   = XLM HERE  AND RL = ONE , EL = ZERO FOR BESSELS
C---------------------------------------------------------------------
      OMEGA = BETA * OMEGA / DABS(DEN)
                 XL = XLM
                 RL = ONE 
      DO 60  L = MINL+1 , MAXL                   ! indexed from 0
                 XL = XL + ONE
          IF( ETANE0 ) THEN
                 RL = GC (L)
                 SL = GCP(L)
            ELSE 
                 SL =  XL * XINV
            ENDIF
          GC (L)  = ( (SL - ALPHA) * GC(L-1) - GCP(L-1) ) / RL
          GCP(L)  =    RL *  GC(L-1)  -  (SL + ALPHA) * GC(L)
          FCP(L)  = OMEGA * ( FCP(L)  -  ALPHA * FC(L) )
          FC (L)  = OMEGA *   FC (L)
   60 CONTINUE
      RETURN
C------------------   ERROR MESSAGES
  100 IFAIL = -1
      WRITE(6,1000) X,ACCH
 1000 FORMAT(' FOR X = ',1PD12.3,'     TRY SMALL-X  SOLUTIONS,',
     *' OR X IS NEGATIVE'/ ,' SQUARE ROOT (ACCURACY) =  ',D12.3/)
                     RETURN
  105 IFAIL = -2                                                        
      WRITE (6,1005) LRANGE,XLMIN,XLM                                    
 1005 FORMAT(/' PROBLEM WITH INPUT ORDER VALUES: LRANGE, XLMIN, XLM = ',    
     *I10,1P2D15.6/)                                                        
                     RETURN                                   
  110 IFAIL =  1                                                        
      WRITE (6,1010) LIMIT, CF1,DCF1, PK,ACCUR                              
 1010 FORMAT(' CF1 HAS FAILED TO CONVERGE AFTER ',I10,' ITERATIONS',/ 
     *' CF1,DCF1,PK,ACCUR =  ',1P4D12.3/)                               
                     RETURN                                       
  120 IFAIL =  2                                                        
      WRITE (6,1020) LIMIT,P,Q,DP,DQ,ACCUR
 1020 FORMAT(' CF2 HAS FAILED TO CONVERGE AFTER ',I7,' ITERATIONS',/  
     *' P,Q,DP,DQ,ACCUR =  ',1P4D17.7,D12.3/)
                     RETURN                                              
      END                                                               
C----------------------------------------------------------------------  
      SUBROUTINE  JWKB   (X,ETA,XL, FJWKB,GJWKB, IEXP)            
      DOUBLE PRECISION    X,ETA,XL, FJWKB,GJWKB, DZERO                      
C----------------------------------------------------------------------
C-----COMPUTES JWKB APPROXIMATIONS TO COULOMB FUNCTIONS  FOR XL .GE. 0.
C-----AS MODIFIED BY BIEDENHARN ET AL. PHYS REV 97 (1955) 542-554
C-----CALCULATED IN SINGLE, RETURNED IN DOUBLE PRECISION VARIABLES
C-----CALLS DMAX1, SQRT, ALOG, EXP, ATAN2, FLOAT, INT     
C     AUTHOR:    A.R.BARNETT   FEB 1981    LAST UPDATE MARCH 1991
C----------------------------------------------------------------------
      REAL    ZERO,HALF,ONE,SIX,TEN,RL35,ALOGE
      REAL    GH2,XLL1,HLL,HL,SL,RL2,GH,PHI,PHI10
      INTEGER IEXP, MAXEXP
      PARAMETER  ( MAXEXP = 300 )
      DATA  ZERO,HALF,ONE,SIX,TEN  /0.0E0, 0.5E0, 1.0E0, 6.0E0, 1.0E1/
      DATA DZERO,RL35,ALOGE /0.0D0, 35.0E0, 0.43429 45 E0 /  
C----------------------------------------------------------------------
CHOOSE MAXEXP NEAR MAX EXPONENT RANGE E.G. 1.D300 FOR DOUBLE PRECISION
C----------------------------------------------------------------------
      GH2   =  X * (ETA + ETA - X)                                         
      XLL1  = DMAX1( XL * XL + XL, DZERO )                                   
      IF( GH2 + XLL1 .LE. ZERO )                 RETURN
      HLL  = XLL1 + SIX / RL35                                           
      HL   = SQRT(HLL)                                                 
      SL   = ETA / HL + HL / X                                             
      RL2  = ONE + ETA * ETA / HLL                                         
      GH   = SQRT(GH2 + HLL) / X                                         
      PHI  = X*GH - HALF*( HL*ALOG((GH + SL)**2 / RL2) - ALOG(GH) )      
      IF ( ETA.NE.ZERO ) PHI = PHI - ETA * ATAN2(X*GH,X - ETA)         
      PHI10 = -PHI * ALOGE                                                
      IEXP  =  INT(PHI10)                                               
      IF ( IEXP.GT.MAXEXP ) THEN
           GJWKB = TEN**(PHI10 - FLOAT(IEXP))               
      ELSE
           GJWKB = EXP(-PHI)                                
           IEXP  = 0                                        
      ENDIF
      FJWKB = HALF / (GH * GJWKB)                                           
C---------------------------------------------------------------------
C     END OF CONTINUED-FRACTION COULOMB & BESSEL PROGRAM  COUL90
C---------------------------------------------------------------------
      RETURN                                                            
      END                                                               

*   ------------------------------------------------------------------
*               H N O R M
*   ------------------------------------------------------------------
*
*   Returns the value of the normalization constant for an (nl)
*   hydrogenic function with nuclear charge ZZ.
*
      DOUBLE PRECISION FUNCTION HNORM(N,L,ZZ)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      M = L + L + 1
      A = N + L
      B = M
      T = A
      D = B
      M = M - 1
      IF (M .EQ. 0) GO TO 2
      DO 1 I = 1,M
      A = A - 1.d0
      B = B - 1.d0
      T = T*A
1     D = D*B
2     HNORM = DSQRT(ZZ*T)/( N*D)
      RETURN
      END

*   ------------------------------------------------------------------
*               H W F
*   ------------------------------------------------------------------
*
*   Returns the value of an unnormalized (nl) hydrogenic function
*   with nuclear charge ZZ and radius r.
*
      DOUBLE PRECISION FUNCTION HWF(N,L,ZZ,R)
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)

      K = N-L-1
      P = 1.d0
      A = 1.d0
      B = K
      C = N+ L
      X = -2.d0*ZZ*R/N
*
*  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0
*
      IF ( X .LT. -150.D0 ) GO TO 5
      IF (K) 1,2,3
3     DO 4 I = 1,K
      P = 1.d0 + A/B*P/C*X
      A = A + 1.d0
      B = B - 1.d0
4     C = C - 1.d0
2     HWF = P*DEXP(X/2.d0)*(-X)**(L+1)
      RETURN
1     WRITE(*,7) N,L,ZZ,R
7     FORMAT(' FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM',/,
     1       ' N = ',I4,'   L = ',I4,'   Z = ',F6.1,'   R = ',F8.4)
      STOP
5     HWF = 0.d0
      RETURN
      END

C**************************************************************
      SUBROUTINE MESHDW(R0,RMAXV,HR)
      use parameters
      use DAT1
      use rpcommon
      implicit double precision (a-h,o-z)
      nrpmax = nxmax+1000
      H = HR
      RP(1) = R0
      NRP = 1
      A = 2.0d0
      DO 2 K=1,nrpmax
        DO 1 J=1,NR
          NRP=NRP+1
          RP(NRP)=RP(NRP-1)+H
          IF(RP(NRP).GT.RMAXV) RETURN
          IF(NRP.GT.(nrpmax-1)) GO TO 3
1       CONTINUE
        H=A*H
        IF(H.GE.DRA) A=1.0d0
2     CONTINUE
3     WRITE(*,100)
100   FORMAT(' LAST POINT LESS THEN RMAXV')
      RETURN
      END

C***************************************************************
! new version of subroutine pulse by Klaus; written January 2014
! fourier transform added by John Emmons; written February 2014

      subroutine pulsekb
      use realconstants1
      use extra
      use realconstants2
      use in1
      use realarray1
      use pulseparams
      use fieldFT

      implicit none

      integer i

      namelist /pulse1/ ee1,alph1,ww1,rr1,cep1,n1up,n1plat,n1down,
     >                  shape1up,shape1down,nextra
      namelist /pulse2/ ee2,alph2,ww2,rr2,cep2,n2up,n2plat,n2down,
     >                  shape2up,shape2down
!
! There will be two pulses, and each pulse is characterized by a number
! of parameters.  Here they are for the first pulse: 
!                   
!   independent parameters
!      rr1:        delay in periods from the start of time
!      n1up:       number of periods for ramp-up
!      n1plat:     number of periods on plateau
!      n1down:     number of periods for ramp-down
!      ww1:        frequency in atomic units
!      ee1:        amplitude in atomic units 
!                  e [a.u.]  =  0.5338(-8) * sqrt ( I [W/cm^2] ) 
!      alph1:      factor by which e1 is multiplied
!      shape1up:   character variable for shape of ramp-up
!      shape1down: character variable for shape of ramp-up
!      cep1:       carrier envelope phase relative to sin(ww1*t)
!
!  set defaults for first pulse 
      ee1        = 0.05338d0    ! 10^14 W/cm^2
      alph1      = 1.0d0        ! full strength
      ww1        = 0.5d0        ! 13.6 eV
      rr1        = 0.0d0
      cep1       = 0.0d0
      n1up       = 20           ! 40-cycle sin^2 pulse
      n1plat     = 0
      n1down     = 20
      shape1up   = 's'
      shape1down = 's'
      nextra     = 0            ! no propagation beyond the end of the pulse
      
      read(80,pulse1)
!
!  set defaults for second pulse 
      ee2        = ee1 
      alph2      = 1.0d0        ! full strength
      ww2        = 2.d0*ww1     ! second harmonic
      rr2        = 0.0d0
      cep2       = 0.0d0
      n2up       = 2*n1up       ! consistent with second harmonic
      n2plat     = 2*n1plat
      n2down     = 2*n1down
      shape2up   = 's'
      shape2down = 's'
      read(80,pulse2)
!
!   now get the dependent parameters        
      period1    = 2.d0*pi/ww1
      tstart1    = rr1*period1
      tup1       = tstart1 + n1up*period1
      tdown1     = tup1 + n1plat*period1
      tend1      = tdown1 + n1down*period1
      period2    = 2.d0*pi/ww2
      tstart2    = rr2*period2
      tup2       = tstart2 + n2up*period2
      tdown2     = tup2 + n2plat*period2
      tend2      = tdown2 + n2down*period2
      tpulses    = max(tend1,tend2)
      ntfinpulse = nint(tpulses/dt)
      tfinal     = tpulses + nextra*period1
      ntfinal    = nint(tfinal/dt)
      cep1rad    = cep1*pi/180.0d0
      cep2rad    = cep2*pi/180.0d0
!
! allocate the arrays
      allocate(timeval(0:ntfinal),envelope1(0:ntfinal),
     >         envelope2(0:ntfinal),field1(0:ntfinal),
     >         field2(0:ntfinal),fieldtot(0:ntfinal),
     >         f(0:ntfinal))
!
!   define the time array and the zeros of the envelope function
      envelope1 = 0.d0
      envelope2 = 0.d0
      field1    = 0.d0
      field2    = 0.d0
      fieldtot  = 0.d0
      f         = 0.d0
      
      do i = 0,ntfinal
        timeval(i) = (dble(i)+0.5d0)*dt
      end do
!
!   now get the characteristic pointers
      istart1 = nint(tstart1/dt)
      iup1    = nint(tup1/dt)
      idown1  = nint(tdown1/dt)
      iend1   = nint(tend1/dt)
      istart2 = nint(tstart2/dt)
      iup2    = nint(tup2/dt)
      idown2  = nint(tdown2/dt)
      iend2   = nint(tend2/dt)
!
!   set the envelope function on the plateaus
      do i = iup1,idown1
        envelope1(i) = 1.0d0
      end do
      do i = iup2,idown2
        envelope2(i) = 1.0d0
      end do
!   
!   what happens now depends on the details of the ramp-on
      do i = istart1,iup1
       if(shape1up == 't') 
     >   envelope1(i) = (timeval(i)-tstart1)/(tup1-tstart1)
       if(shape1up == 's') 
     >   envelope1(i) = 
     >     sin(0.5d0*pi*(timeval(i)-tstart1)/(tup1-tstart1))**2
       if(shape1up == 'g')
     >   envelope1(i) = 
     >     exp(-4.0*log(2.d0)*((tup1-timeval(i))/(tup1-tstart1))**2)
      end do
      do i = istart2,iup2
       if(shape2up == 't') 
     >   envelope2(i) = (timeval(i)-tstart2)/(tup2-tstart2)
       if(shape2up == 's') 
     >   envelope2(i) = 
     >     sin(0.5d0*pi*(timeval(i)-tstart2)/(tup2-tstart2))**2
       if(shape2up == 'g')
     >   envelope2(i) = 
     >     exp(-4.0*log(2.d0)*((tup2-timeval(i))/(tup2-tstart2))**2)
      end do
!   
!   similar for the ramp-down, except that we start from 1 and go to 0
      do i = idown1,iend1
       if(shape1down == 't') 
     >   envelope1(i) = (tend1-timeval(i))/(tend1-tdown1)
       if(shape1down == 's') 
     >   envelope1(i) = 
     >     sin(0.5d0*pi*(tend1-timeval(i))/(tend1-tdown1))**2
       if(shape1down == 'g')
     >   envelope1(i) = 
     >   exp(-4.0*log(2.d0)*((tdown1-timeval(i))/(tdown1-tend1))**2)
      end do
      do i = idown2,iend2
       if(shape2down == 't') 
     >   envelope2(i) = (tend2-timeval(i))/(tend2-tdown2)
       if(shape2down == 's') 
     >   envelope2(i) = 
     >     sin(0.5d0*pi*(tend2-timeval(i))/(tend2-tdown2))**2
       if(shape2down == 'g')
     >   envelope2(i) = 
     >   exp(-4.0*log(2.d0)*((tdown2-timeval(i))/(tdown2-tend2))**2)
      end do
!
!   for Gaussians, we need to make a small correction, so that the
!   envelope function actually goes to zero.
      if(shape1up == 'g') then 
        gvalstart1 = envelope1(istart1)
        do i = istart1,iup1
          correction = gvalstart1*(tup1-timeval(i))/(tup1-tstart1)
          envelope1(i) = envelope1(i) - correction
        end do
      endif 
      if(shape1down == 'g') then 
        gvalend1 = envelope1(iend1)
        do i = idown1,iend1
          correction = gvalend1*(timeval(i)-tdown1)/(tend1-tdown1)
          envelope1(i) = envelope1(i) - correction
        end do
      endif 
      if(shape2up == 'g') then 
        gvalstart2 = envelope2(istart2)
        do i = istart2,iup2
          correction = gvalstart2*(tup2-timeval(i))/(tup2-tstart2)
          envelope2(i) = envelope2(i) - correction
        end do
      endif 
      if(shape2down == 'g') then 
        gvalend2 = envelope2(iend2)
        do i = idown2,iend2
          correction = gvalend2*(timeval(i)-tdown2)/(tend2-tdown2)
          envelope2(i) = envelope2(i) - correction
        end do
      endif
!
!     now we finally have the envelope functions and can calculate the fields
      write(90,1000)
1000  format('#   time (a.u.)   time(fs)      ',
     >       'envelope1     envelope2      ',
     >       'field1        field2        field')
      do i = 0,ntfinal
        field1(i) = alph1*ee1*envelope1(i)*sin(ww1*timeval(i)+cep1rad)
        field2(i) = alph2*ee2*envelope2(i)*sin(ww2*timeval(i)+cep2rad)
        fieldtot(i) = field1(i)+field2(i)
        write(90,1001) timeval(i),0.02418884d0*timeval(i),
     >   envelope1(i),envelope2(i),field1(i),field2(i),fieldtot(i)
       end do
1001  format(1p10e14.6)

!     compute the fourier transform of the total field
!     added by John Emmons; written February 2014

       wmin = 0.0001d0
       wmax = 2.000d0
       wdel = 0.0001d0

       nomega = nint((wmax-wmin)/wdel)
       deltat = timeval(2)-timeval(1)

      do 31 n=1,nomega
       w(n) = wmin+dble(n-1)*wdel
       gcos = 0.0d0
       gsin = 0.0d0

       do 41 i=0,ntfinal
        gcos = gcos + cos(w(n)*timeval(i))*fieldtot(i)
        gsin = gsin + sin(w(n)*timeval(i))*fieldtot(i)
!      old debug printing
       !  if (n.eq.380.and.mod(i,100).eq.0) then
       !  print *,'debug for n = ',n,'   i = ',i
       !  print *,i,timeval(i),fieldtot(i) 
       !  print *,n,w(n)
       !  print *,gcos,gsin
       ! endif
41      continue

       greal(n) = gcos*deltat
       gimag(n) = gsin*deltat
31     continue

!     write the fourier transform to a file
      do 51 n=1,nomega
        write(1030,1011) n,w(n),greal(n),gimag(n),
     >                sqrt(greal(n)**2+gimag(n)**2)
51    continue

1011  format(i10,1p10e16.8)

      call flush(90)
      return
      end

C**************************************************************************
C
C This will do the hydrogen case with analytic wavefunctions and potentials
C
      subroutine radhyd
      use parameters
      use realconstants1
      use realarray1
      use realconstants2
      use discr
      use complexarray3
      use keys
      use in1
      use cmplxconstants
      use VDD
      use rpcommon
      use DAT1
      
      implicit double precision (a-h,o-z)

C*** INITIAL RADIAL FUNCTION AND CHANNEL POTENTIALS
      do 100 i=0,nx
!       g(i) = 2.d0*x(i)*exp(-x(i))      ! initial radial function
        xx=x(i)
        if(i.eq.0) xx = 1.d-40            
        do 101 j=0,nc                     ! channel potentials
          cm = dble(j)
          v(i,j) = dcmplx(cm*(cm+1.d0)/(2.d0*xx*xx) - zet/xx,0.d0) 
C--- Gobbler:
          if(x(i).gt.gbr) then
            v(i,j) = v(i,j) - ci*dcmplx(agbr*(x(i)/gbr - 1.d0)**3,0.d0)
          endif            
101     continue
100   continue

C**** radial functions for discrete states if key1 = 1
      if(key1 == 1) then
        do 110 kd = 1, nf
          un = hnorm(nn(kd),ll(kd),1.d0)
          do 111 i = 0,nx1-1
            fd(i,kd)   = un * hwf(nn(kd),ll(kd),1.d0,x(i))
*
            if(nn(kd) == nnauto .and. ll(kd) == llauto) then
              g(i) = fd(i,kd)             ! initial radial function
            endif
*      
111       continue
110     continue
!       STOP 'after 110'
      endif
C***** mesh and screening potential for continuum if key4 = 1
      if(key4.ge.1) then
        z = 1.0d0
        hr = 0.005d0
        nr = 100
        dra = h - 0.001d0
        call meshdw(1.d-4,x(ngob),hr)
        do 130 i=1,nxmax+1000
          vd(i) = 0.d0
130     continue
      endif

c-- initial wavefunction to complex and boundary conditions 
!      do 150 i = 0,nx
!        q(i,0) = dcmplx(g(i),0.d0)
!150   continue
!      do 160 jj=0,nc
!        q(nx,jj) = czero
!        q(0,jj) = czero
!160   continue
!
      return
      end
       
C***************************************************************************
C
C This will do the Neon case with our 2013/2014 wavefunctions and potentials
C
      subroutine radneon
      use parameters
      use realconstants1
      use realarray1
      use realconstants2
      use discr
      use complexarray3
      use keys
      use in1
      use cmplxconstants
      use VDD
      use rpcommon
      use DAT1

      implicit double precision (a-h,o-z)
      
      double precision, allocatable :: v1(:)
      double precision, allocatable :: rr(:),vv(:),
     >    g1s(:),g2s(:),g3s(:),g4s(:),g5s(:),
     >           g2p(:),g3p(:),g4p(:),g5p(:),
     >                  g3d(:),g4d(:),g5d(:)
      
      allocate(v1(0:nxmax))
      np = 464
      allocate(rr(np),vv(np),
     >    g1s(np),g2s(np),g3s(np),g4s(np),g5s(np),
     >            g2p(np),g3p(np),g4p(np),g5p(np),
     >                    g3d(np),g4d(np),g5d(np))

C*** INITIAL RADIAL FUNCTION AND CHANNEL POTENTIALS
      open(65,file='ne.wfn',status='unknown')
      do 100 jn=1,np       
        read(65,1000) rr(jn),vv(jn),
     >                g1s(jn),g2s(jn),g3s(jn),g4s(jn),g5s(jn),
     >                        g2p(jn),g3p(jn),g4p(jn),g5p(jn),
     >                                g3d(jn),g4d(jn),g5d(jn)

100   continue 
1000  format(14d14.6)           
      close(65)
    
      do 110 i=0,nx
        if(x(i).gt.rr(np)) then
          g(i) = 0.d0
          v1(i) = -zet
          goto 120
        endif
        call parinv(x(i),rr,g2p,np,g(i))
        call parinv(x(i),rr,vv,np,v1(i))
120     continue
        xx = x(i)
        if(i.eq.0) xx = 1.d-40            
        do 130 j=0,nc                     ! channel potentials
          cm = dble(j)
          v(i,j) = dcmplx(cm*(cm+1.d0)/(2.d0*xx*xx)+v1(i)/xx,0.d0) 
C--- Gobbler:
          if(x(i).gt.gbr) then
            v(i,j) = v(i,j) - ci*dcmplx(agbr*(x(i)/gbr - 1.d0)**3,0.d0)
          endif            
130     continue
110   continue

C**** RADIAL FUNCTIONS FOR DISCRETE STATES IF KEY1 = 1
      if(key1 == 1) then
        do 140 i = 0,nx1
          if(x(i).gt.rr(np)) goto 140
          call parinv(x(i),rr,g1s,np,fd(i,1))
          call parinv(x(i),rr,g2s,np,fd(i,2))
          call parinv(x(i),rr,g3s,np,fd(i,3))
          call parinv(x(i),rr,g4s,np,fd(i,4))
          call parinv(x(i),rr,g5s,np,fd(i,5))
          call parinv(x(i),rr,g2p,np,fd(i,6))
          call parinv(x(i),rr,g3p,np,fd(i,7))
          call parinv(x(i),rr,g4p,np,fd(i,8))
          call parinv(x(i),rr,g5p,np,fd(i,9))
          call parinv(x(i),rr,g3d,np,fd(i,10))
          call parinv(x(i),rr,g4d,np,fd(i,11))
          call parinv(x(i),rr,g5d,np,fd(i,12))
          g(i) = fd(i,6)             ! initial radial function is 2p
140     continue
      endif
C*****mesh and screening potential for continuum if key4=1
      if(key4.ge.1) then
        z = 10.0d0
        hr = 0.005d0
        nr = 100
        dra = h - 0.001d0
        call meshdw(1.d-4,x(ngob),hr)
        do 160 i = 1,nrp
          vd(i)= z - 1.d0
160     continue   
        do 170 i=1,nrp
          if(rp(i).gt.rr(np)) goto 170
          call parinv(rp(i),rr,vv,np,vk)
          vd(i) = z+vk
170     continue
      endif

      return
      end
************************************************************************
*                                                                      *
*  THE FOLLOWING SUBROUTINE CALCULATES CLEBSCH-GORDAN COEFFICIENTS.    *
*  NOTE: THE INTEGER ARGUMENTS ARE TWICE THE ACTUAL ANGULAR MOMENTA    *
*        AND THEIR Z-COMPONENTS, RESPECTIVELY.                         *
*                                                                      *
************************************************************************
      SUBROUTINE CLEGOR (J1,J2,J3,M1,M2,M3,CG,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (IFAK = 500)

      COMMON / LOF / FAK (IFAK)

      CG = ZERO
      IER = 0

      J  = (J1+J2+J3)/2
      M  = (M1+M2-M3)/2
      JX = (J2-J1+J3)/2
      JY = (J1-J2+J3)/2
      JZ = (J1+J2-J3)/2
*
*  CHECK OF THE COUPLING RULES
*
      IF ( (J1.LT.0).OR.(J2.LT.0).OR.(J3.LT.0).OR.
     1     (JX.LT.0).OR.(JY.LT.0).OR.(JZ.LT.0).OR.
     2     ((J1-IABS(M1)).LT.0).OR.((J2-IABS(M2)).LT.0).OR.
     3     ((J3-IABS(M3)).LT.0).OR.
     4     (MOD(J1,2).NE.MOD(IABS(M1),2)).OR.
     5     (MOD(J2,2).NE.MOD(IABS(M2),2)).OR.
     6     (MOD(J3,2).NE.MOD(IABS(M3),2)).OR.
     7     (MOD(J1+J2,2).NE.MOD(J3,2)).OR.(M.NE.0) )  THEN
       IER = 1
       RETURN
      ENDIF
*
*  CALCULATION OF SPECIAL VALUES
*
      IF ((J1.EQ.0).OR.(J2.EQ.0)) THEN
       CG = ONE
       RETURN
      ENDIF
      IF (J3.EQ.0) THEN
       CG = ONE/SQRT(J1+ONE)
       N = (J1-M1)/2
       IF (MOD(N,2).EQ.1) CG = -CG
       RETURN
      ENDIF
      IF ((M1.EQ.0).AND.(M2.EQ.0)) THEN
       IF (MOD(J,2).EQ.1) RETURN
       N10 = J/2
       L1  = J1/2
       L2  = J2/2
       L3  = J3/2
       N11 = (L1+L2-L3)/2
       CG = (J3+ONE)*EXP(FAK(JX+1)+FAK(JY+1)+FAK(JZ+1)-
     1      FAK(J+2)+TWO*(FAK(N10+1)-FAK(N10-L1+1)-
     2      FAK(N10-L2+1)-FAK(N10-L3+1)))
       CG = SQRT (CG)
       IF (MOD(N11,2).EQ.1) CG = -CG
       RETURN
      ENDIF
*
*  CALCULATION OF THE GENERAL CASE
*
      N1 = (J1+M1)/2
      N2 = (J1-M1)/2
      N3 = (J2+M2)/2
      N4 = (J2-M2)/2
      N5 = (J3+M3)/2
      N6 = (J3-M3)/2
      N7 = (J2-J3-M1)/2
      N8 = (J1-J3+M2)/2
      KMIN = MAX(0,N7,N8)
      KMAX = MIN(JZ,N2,N3)
      L = KMIN
      CONST = FAK(N1+1)+FAK(N2+1)+FAK(N3+1)+FAK(N4+1)
     1       +FAK(N5+1)+FAK(N6+1)+FAK(JX+1)+FAK(JY+1)
     2       +FAK(JZ+1)-FAK(J+2)-TWO*(FAK(JZ-L+1)+FAK(N2-L+1)
     3       +FAK(N3-L+1)+FAK(L-N7+1)+FAK(L-N8+1)+FAK(L+1))
      N = KMAX-KMIN+1
      I = N-1
      SUM = ONE
      IA = JZ-L
      IB = N2-L
      IC = N3-L
      ID = L+1-N7
      IE = L+1-N8
      IG = L+1

   10 AI   = (IA-I)*(IB-I)*(IC-I)
      BIP1 = (ID+I)*(IE+I)*(IG+I)
      SUM  = ONE-AI*SUM/BIP1
      I    = I-1
      IF (I.GE.0) GOTO 10
      CG = SQRT((J3+ONE)*EXP(CONST)) * SUM
      IF (MOD(L,2).EQ.1) CG = -CG

      RETURN
      END
************************************************************************
*                                                                      *
*  THE FOLLOWING SUBROUTINE CALCULATES LOGARITHMS FACTORIALS. THESE    *
*  ARE NEEDED IN THE CLEBSCH-GORDAN PACKAGE.                           *
*                                                                      *
************************************************************************
      SUBROUTINE LFAK
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      PARAMETER (ZERO=0.0D0,HALF=0.5D0,ONE=1.0D0,TWO=2.0D0)
      PARAMETER (IFAK=500)

      COMMON / LOF / FAK(IFAK)

      FAK(1) = ZERO
      DO 1 I=1,IFAK-1
       FAK(I+1) = FAK(I)+LOG(REAL(I))
1     CONTINUE
      RETURN
      END

C**************************************************************
c needs clean-up, but works
      SUBROUTINE DWAVEkb(L,E,PHASE,dw)
      use parameters
      use VDD
      use rpcommon
      use DAT1
      implicit double precision (a-h,o-z)
      dimension dw(nxmax+1000)
      DIMENSION F1(500),F2(500),G1(500),G2(500)
      DIMENSION FC1(500),FC2(500),GC1(500),GC2(500)
      dimension fc(0:500),fcp(0:500),gc(0:500),gcp(0:500)
      pi=dacos(-1.d0)
      L1=L+1
      PHASE=0.0d0
      ANORM=0.0d0
      SE=DSQRT(E)
      AL=L
      PHI=0.0d0
      M=1
      DO 5 I=1,L1
5     PHI=PHI+DLOG(dble(2*I+1))
      R1=DEXP((PHI-35.0d0)/(AL+1.0d0))/SE
      NRP2=NRP-2
      DO 1 J=1,NRP2
      IF(RP(J).GE.R1) GO TO 2
1     DW(J)=0.0d0
      RETURN
2     J1=(J-1-(J-1)/NR*NR)/(NR-1)
      ALL=AL*(AL+1.0d0)
      NR1=NR+1
      NR2=NR+2
      NR3=NR2+1
      N1=J-J1
      N2=N1+1
      N3=N2+1
      X=RP(N1)*SE
      AN=(VD(N1)-Z)/SE
      AL1=AL+1.0d0
6     B1=1.d0
      XX=X*X
      B2=AN/AL1*X
      FL=B1+B2
      S=0.0d0
7     S=S+1
      BS=B2
      B2=(2.d0*AN*B2*X-B1*XX)/((2.d0*AL1+S)*(S+1.d0))
      B1=BS
      FL=FL+B2
      IF(B1*1.d-15.NE.0..OR.B2*1.d-15.NE.0) GO TO 7
      DO 8 I=1,L1
8     FL=FL*X/(dble(2*I-1))
      GO TO (11,12) , M
11    DW(N1)=FL
      AN=(VD(N2)-Z)/SE
      X=RP(N2)*SE
      M=2
      GO TO 6
12    DW(N2)= FL
      H=RP(N2)-RP(N1)
      HH=H*H
      C1=5.d0/6.d0
      C2=1.d0/12.d0
      DO 3 J=N3,NRP
      J1=J-1
      RJ=RP(J)
      RJ1=RP(J1)
      H1=H
      H=RJ-RJ1
      HH=H*H
      JD=(H1+H)/(2.9d0*H1)
      J2=J-2-JD
      RJ2=RP(J2)
      Y=(2.d0+C1*HH*(ALL/RJ1**2-2.d0*(Z-VD(J1))/RJ1-E))*DW(J1)
     *-(1.d0-C2*HH*(ALL/RJ2**2-2.d0*(Z-VD(J2))/RJ2-E))*DW(J2)
      DW(J)=Y/(1.d0-C2*HH*(ALL/RJ**2-2.d0*(Z-VD(J))/RJ-E))
3     CONTINUE
      DW1=DW(NRP-1)
      DW2=DW(NRP)
      ZA=Z-VD(NRP)
      AN=-ZA/SE
      X1=RP(NRP-1)*SE
      X2=RP(NRP)*SE
      CALL COUL90(X1,AN,0.d0,L,FC,GC,FCP,GCP,0,ifail)
      F1(L1)=FC(L)
      FC1(L1)=FCP(L)
      G1(L1)=GC(L)
      GC1(L1)=GCP(L)
      CALL COUL90(X2,AN,0.d0,L,FC,GC,FCP,GCP,0,ifail)
      F2(L1)=FC(L)
      FC2(L1)=FCP(L)
      G2(L1)=GC(L)
      GC2(L1)=GCP(L)
c      CALL COUL(X1,AN,L,L,F1,FC1,G1,GC1,1.0E-6,100.0)
c      CALL COUL(X2,AN,L,L,F2,FC2,G2,GC2,1.0E-6,100.0)
13    DET=F1(L1)*G2(L1)-F2(L1)*G1(L1)
      IF(DET.EQ.0.d0) GO TO 10
      AC=(DW1*G2(L1)-DW2*G1(L1))/DET
      AS=(DW2*F1(L1)-DW1*F2(L1))/DET
      ANORM=AC*AC+AS*AS
      PHASE=pi/2.d0
      IF(AC.EQ.0.d0) GO TO 9
      PHASE=PHASE+DATAN(AS/AC)+FACOUZ(E,L,ZA)
      IF(AC.LT.0.d0) PHASE=PHASE+pi
9     ANORM=1./DSQRT(ANORM*pi*SE)
10    DO 4 J=1,NRP
4     DW(J)=DW(J)*ANORM  
      RETURN
      END
