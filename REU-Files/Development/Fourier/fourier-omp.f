*
      program foutrap
      implicit double precision (a-h,o-z)
      parameter(ndim=25000)
      dimension f(0:ndim),t(0:ndim),w(ndim),greal(ndim),gimag(ndim)
*
      include 'omp_lib.h'
*
      pi = acos(-1.0d0)
      wmin = 0.0001d0
      wmax = 2.000d0
      wdel = 0.0001d0

      nomega = nint((wmax-wmin)/wdel)
      print *, 'nomega = ',nomega
*
* read the function
*
      i = 0
      read(1,*)
10    read(1,*,end=20) t(i),c2,c3,c4,c5,c6,f(i)
      i = i+1
      goto 10 
20    icount = i
      deltat = t(2)-t(1)
      print *,' number of time points:',icount
      print *,' delta t = ',deltat
*
* now calculate the Fourier transform using simple trapezoidal rule;
* no need to worry about the end points, since f is basically zero.
*
c$omp parallel do private(n,gcos,gsin,i)
      do 30 n=1,nomega
       w(n) = wmin+dble(n-1)*wdel
       gcos = 0.0d0
       gsin = 0.0d0
       do 40 i=0,icount
        gcos = gcos + cos(w(n)*t(i))*f(i)
        gsin = gsin + sin(w(n)*t(i))*f(i)
        if (n.eq.380.and.mod(i,100).eq.0) then
         print *,'debug for n = ',n,'   i = ',i
         print *,i,t(i),f(i) 
         print *,n,w(n)
         print *,gcos,gsin
        endif
40     continue
       greal(n) = gcos*deltat
       gimag(n) = gsin*deltat
30    continue
c$omp end parallel do
*
      do 50 n=1,nomega
        write(2,1000) n,w(n),greal(n),gimag(n),
     >                sqrt(greal(n)**2+gimag(n)**2)
50    continue

*
1000  format(i10,1p10e16.8)
      stop
      end
