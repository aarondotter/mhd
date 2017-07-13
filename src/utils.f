c
      subroutine smm(l,m,n,a,ia,b,ib,c,ic)
      implicit double precision (a-h,o-z), integer (i-n)
c
      dimension a(ia,1),b(ib,1),c(ic,1)
      do 30 ll=1,l
      do 20 nn=1,n
      c(ll,nn)=0.d0
      do 10 mm=1,m
10    c(ll,nn)=c(ll,nn)+a(ll,mm)*b(mm,nn)
20    continue
30    continue
      return
      end
c
      subroutine smv(m,n,a,ia,x,y)
      implicit double precision (a-h,o-z), integer (i-n)
c
      dimension a(ia,1),x(1),y(1)
      do 20 mm=1,m
      y(mm)=0.d0
      do 10 nn=1,n
10    y(mm)=y(mm)+a(mm,nn)*x(nn)
20    continue
      return
      end
c
      subroutine smmtr(l,m,n,a,ia,b,ib,c,ic)
      implicit double precision (a-h,o-z), integer (i-n)
c
      dimension a(ia,1),b(ib,1),c(ic,1)
      do 30 ll=1,l
      do 20 nn=1,n
      c(ll,nn)=0.d0
      do 10 mm=1,m
10    c(ll,nn)=c(ll,nn)+a(ll,mm)*b(nn,mm)
20    continue
30    continue
      return
      end
c
      double precision function wgt(t,t1,t2,t3,t4,wmin,ex)
      implicit double precision (a-h,o-z), integer (i-n)
c
      if (t.lt.t1) then
         wgt = 1.d0
      else if (t.lt.t2) then
         trel = (t - t1)/(t2 - t1)
         wgt  = 1.d0 - (1.d0 - wmin)*flip(trel,ex)
      else if (t.lt.t3) then
         wgt  = wmin
      else if (t.lt.t4) then
         trel = (t - t3)/(t4 - t3)
         wgt  = wmin + (1.d0 - wmin)*flip(trel,ex)
      else
         wgt  = 1.d0
      end if
c
      return
      end
c
      double precision function flip(x,ex)
      implicit double precision (a-h,o-z), integer (i-n)
c
      c = 0.5d0**(1.d0-ex)
c
      if (x.lt.0.d0) then
          flip = 0.d0
      else if (x.le.0.5d0) then
          flip = c*x**ex
      else if (x.lt.1.d0) then
          flip = 1.d0 - c*(1.d0-x)**ex
      else
          flip = 1.d0
      end if
c
      return
      end
