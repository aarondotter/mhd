      subroutine fscan
      implicit double precision (a-h,o-z), integer (i-n)
c
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c***********************************************************************
c
c     scan hydrogen and helium ionization fractions to find minimum of c
c     free energy                                                      c
c                                                                      c
c -------------------------------------------------------------------- c
c    modifications (May-July 1990): indroduction of molecular scan and c
c    other improvements, such as asking for a helium, or hydrogen scan c
c    only (flags ihesc and ihsc [common/keep/]). furthermore,          c
c    previously the library routine isamax was used to find the        c
c    minimum. In case of equal free energy, however, isamax did not    c
c    return the appropriate point. this minimum search is now 'hard-   c
c    coded' here.
c -------------------------------------------------------------------- c
c                                                                      c
c***********************************************************************
c
      include 'types'
      include 'parms'
      include 'coms'
c
      data niter
     .    /   30/
c
      dimension xx(5), xf(5)
c
      kh = kz(1)
c
      if ( kh .le. 0 .or. ihesc.eq.1 ) go to 100
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     scan hydrogen-proton balance
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     set indices and charges for the "neutral" and "ion"
c
      isn = ish
      isp = ispr
      if(wwt(isn).lt.0.1d0 .and. wwt(isp).lt.0.1d0) goto 90
      cn = 0.0d0
      ccn= 1.0d0
      rcn= 1.0d0
      cp = 1.0d0
      iscan = 1
      go to 200
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     scan hydrogen molecule-atom balance
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     set indices and charges for the "neutral" and "ion"
c
 90   isn = ish2
      isp = ish
      if(wwt(isn).lt.0.1d0 .or. rholog.gt.-1.d0) goto 100
      cn = 0.0d0
      ccn= 2.0d0
      rcn= 0.5d0
      cp = 0.0d0
      iscan = 2
      go to 200
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     scan he-he+ balance
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  100 khe = kz(2)
      if ( ihsc.eq. 1 ) return
      if ( khe .le. 0 ) return
      isn = ichm1(khe)
      isp = isn + 1
      if(wwt(isn).lt.0.1d0 .and. wwt(isp).lt.0.1d0) return
      cn = 0.0d0
      ccn= 1.0d0
      rcn= 1.0d0
      cp = 1.0d0
      iscan = 3
      go to 200
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     scan he+ - he++ balance
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
  110 continue 
      isp = isp + 1
      isn = isn + 1
      if(wwt(isp).lt.0.1d0) return
      cn = 1.0d0
      ccn= 1.0d0
      rcn= 1.0d0
      cp = 2.0d0
      iscan = 4
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     scan procedure
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     save original neutral, ion, and excess electron occupation numbers
c
  200 xne = sn(ise) - cp * sn(isp) - cn * sn(isn)
      xnn = sn(isn)*ccn
      xnp = sn(isp)
      xnt = xnn + xnp
c
c     set trial range of ionization fraction
c
      x1    = -20.d0
      x5    = - 1.d-10
      xx(1) = 10.0d0**x1
      xx(5) = 10.0d0**x5
      xx(3) = dsqrt( xx(1) * xx(5) )
c
c     calculate free energy at endpoints and midpoints
c
      do 210 i = 1, 5, 2
      sn(isn) = ( 1.0d0 - xx(i) ) * xnt * rcn
      sn(isp) = xx(i) * xnt
      sn(ise) = cp * sn(isp) + cn * sn(isn) + xne
      call ftot
      xf(i) = f(1) + f(2) + f(3) + f(4)
  210 continue
      iter = 0
c
c     insert new midpoints on range
c
c---------------------- basic scanning iteration ------------------
  250 iter = iter + 1
      do 270 i = 2, 4, 2
      xx(i) = dsqrt( xx(i-1) * xx(i+1) )
      sn(isn) = ( 1.0d0 - xx(i) ) * xnt * rcn
      sn(isp) = xx(i) * xnt
      sn(ise) = cp * sn(isp) + cn * sn(isn) + xne
      call ftot
      xf(i) = f(1) + f(2) + f(3) + f(4)
 270  continue
c
c     find minimum, reset endpoints and midpoint of range
c
      xfmi = xf(1)
      imin = 1
      do 280 i=2,5
      if(xf(i).le.xfmi) then
      xfmi = xf(i)
      imin = i
      end if
 280  continue
c
      if ( imin .le. 2 ) then
         xx(5) = xx(3)
         xx(3) = xx(2)
         xf(5) = xf(3)
         xf(3) = xf(2)
      else if ( imin .eq. 3 ) then
         xx(1) = xx(2)
         xx(5) = xx(4)
         xf(1) = xf(2)
         xf(5) = xf(4)
      else if ( imin .ge. 4 ) then
         xx(1) = xx(3)
         xx(3) = xx(4)
         xf(1) = xf(3)
         xf(3) = xf(4)
      end if
c
c     check for convergence
c
      xdif = log10( xx(3) / xx(1) )
      if ( xdif .gt. 1.d-3 .and. iter .le. niter ) go to 250
c
c------------------ end basic scanning iteration ------------------
c
      write(6,9001) iscan,imin,iter,xx(imin),xf(imin)
c
      sn(isn) = ( 1.0d0 - xx(imin) ) * xnt * rcn
      sn(isp) = xx(imin) * xnt
      sn(ise) = cp * sn(isp) + cn * sn(isn) + xne
      if(iscan.ne.2) call new eta
c
      go to ( 90, 100, 110, 300 ) iscan
c
300   return
c
 9001 format(' scan: iscan,imin,iter,x,f ',3i3,1pg12.3,e25.15)
      end
