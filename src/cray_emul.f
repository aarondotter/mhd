      double precision function cvmgp(x1,x2,x3)
      implicit double precision (a-h,o-z), integer (i-n)
      if(x3)1,2,2
    1 cvmgp=x2
      return
    2 cvmgp=x1
      return
      end
c 
      double precision function cvmgz(x1,x2,x3)
      implicit double precision (a-h,o-z), integer (i-n)
      if(x3)1,2,1
    1 cvmgz=x2
      return
    2 cvmgz=x1
      return
      end
c 
      double precision function cvmgn(x1,x2,x3)
      implicit double precision (a-h,o-z), integer (i-n)
      if(x3)2,1,2
    1 cvmgn=x2
      return
    2 cvmgn=x1
      return
      end
c 
      double precision function cvmgt(x1,x2,x3)
      implicit double precision (a-h,o-z), integer (i-n)
      logical x3
      if(x3) then
        cvmgt=x1
      else
        cvmgt=x2
      end if
      return
      end
c 
      double precision function cvmgm(x1,x2,x3)
      implicit double precision (a-h,o-z), integer (i-n)
      if(x3)2,1,1
    1 cvmgm=x2
      return
    2 cvmgm=x1
      return
      end
c 
      double precision function sdot(n,a,i1,b,i2)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1),b(1)
      sdot=0.d0
      j1=1
      j2=1
      do 1 i=1,n
      sdot=sdot+a(j1)*b(j2)
      j1=j1+i1
      j2=j2+i2
    1 continue
      return
      end 
c
      double precision function ssum(n,a,i1)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1)
      ssum=0.d0
      j1=1
      do 1 i=1,n
      ssum=ssum+a(j1)
      j1=j1+i1
    1 continue
      return
      end
c 
      subroutine sscal(n,fak,a,i1)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1)
      j1=1
      do 1 i=1,n
      a(j1)=fak*a(j1)
    1 j1=j1+i1
      return
      end
c 
      integer function ismin(n,a,nstep) 
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1)
      ndim=n*nstep
      ismin=1
      k=1 
      x=a(1)
      do 1 i=1,ndim,nstep
      if(a(i).ge.x)go to 1
      ismin=k
      x=a(i)
    1 k=k+1
      return
      end
c 
      subroutine scopy(n,a,na,b,nb)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=a(ia)
      ia=ia+na
      ib=ib+nb
    1 continue
      return
      end 
c
      integer function ismax(n,a,nstep) 
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1)
c
      ndim=n*nstep
      x=a(1)
      ismax=1
      k=1 
      do 1 i=1,ndim,nstep
      if(a(i).le.x)go to 1
      ismax=k
      x=a(i)
    1 k=k+1
      return
      end
c 
      integer function isamax(n,a,nstep) 
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1)
c
      ndim=n*nstep
      x=dabs(a(1))
      isamax=1
      k=1 
      do 1 i=1,ndim,nstep
      aa=dabs(a(i))
      if(aa.le.x)go to 1
      isamax=k
      x=aa
    1 k=k+1
      return
      end
c 
      subroutine saxpy(n,c,a,na,b,nb)
      implicit double precision (a-h,o-z), integer (i-n)
      dimension a(1),b(1)
      ia=1
      ib=1
      do 1 i=1,n
      b(ib)=c*a(ia)+b(ib)
      ia=ia+na
    1 ib=ib+nb
      return
      end 
c
      subroutine leq(a,b,nn,mm,ia,ib,err)
      implicit double precision (a-h,o-z), integer (i-n)
c......... modified version for s/r solvfi 
c......... (called by finddn, called by fmin):
c......... err does not yield the actual determinant, because that
c......... could in some cases be zero by underflow, when several
c......... a(i,i)=1.d-1000 occur. thus final scaling of err with detsc
c......... is suppressed, and the numerical value of err can not be
c......... used except as a test of succes in leq.
      dimension a(1),b(1)
c
c     this routine will find the inverse of a matrix by the method of
c     partial pivoting and gaussian elimination if b is set equal to
c     the identity matrix
c
c     n - dimension of segment of a to be used
c     m - number of right hand columns of b to be used
c     ia - the total number of rows in large array a
c     ib - the total number of rows in large array b
c     the matrix equation is    ax=b
c     err = det(a)
c     if mm = 0 leq calculates err = det(a)
c
c  note on modification on 23/1 1985:
c
c  previously, the routine contained the statements
c
c      do 14 k=i2,m1,ib
c      b(k)=b(k)+b(i1)*r
c   14 i1=i1+ib
c
c  this caused problems, on some computers, for m = ib = 1.
c  then m1 = 1 and the loop was skipped when i2 .gt. 1.
c  this has been changed.
c
      common/prtmat/iprmat
      if(iprmat.gt.0) then
      mmmm = ia*ia
      write(6,7789) (a(jj),jj=1,mmmm)
 7789 format(' a(ii) in leq',/(1x,1p5g15.6))
      end if
c
      n=nn
      m=mm
      err=0.0d0
      detsc=1.0d0
c
c     treat n .le. 1 separately
      if(n-1) 280,281,285
c     no equation
  280 write(6,295) n
      return
c     n = 1
  281 err=a(1)
      if(m.le.0) return
      ai=1.d0/err
      m1=ib*m
      do 282 j=1,m1,ib
  282 b(j)=b(j)*ai
      return
c
c
c     find maximum element in each row and divide the row by it
  285 n1=ia*n
      ia1=ia+1
      m1=ib*m
      do 1 i=1,n
      r= dabs(a(i))
      do 2 j=  1,n1,ia
      ij=j+i-1
    2 r=max(r, dabs(a(ij)))
      if(r)31,30,31
   30 write(6,298)i
      return
   31 do 3 j=1,n1,ia
      ij=j+i-1
    3 a(ij)=a(ij)/r
      if(m.eq.0) go to 1
      do 4 j=1,m1,ib
      ij=j+i-1
    4 b(ij)=b(ij)/r
    1 detsc=detsc*r
c
c
c     find maximum element in the i'th column
      n2=n-1
      do 5 i=1,n2
      ialow=(i-1)*ia+i
      iaup=(i-1)*ia+n
      r= dabs(a(ialow))
      ialow1=ialow+1
      imax=ialow
      do 6 j=ialow1,iaup
      if(r- dabs(a(j)))7,6,6
    7 imax=j
      r= dabs(a(j))
    6 continue
      if(imax-ialow)8,8,9
c     replace the i'th row with the row that has the maximum element in
c         the respective column and put the i'th row in its place
    9 im=imax
   72 if(im-ia)70,70,71
   71 im=im-ia
      go to 72
   70 do 10 j=1,n1,ia
      jj=i+j-1
      ji=im+j-1
      r=a(jj)
      a(jj)=a(ji)
   10 a(ji)=r
c     change sign of determinant
      detsc=-detsc
c
      if(m.eq.0) go to 8
      do 11 j=1,m1,ib
      jj=i+j-1
      ji=im+j-1
      r=b(jj)
      b(jj)=b(ji)
   11 b(ji)=r
c     multiply the i'th row by (the negative of each i'th column element
c       below the diagonal divided by the diagonal element) and add the
c     resulting row to the respective row of the element used
    8 iaup1=iaup-1
c
c
      do 12 j=ialow,iaup1
      if(a(ialow))32,33,32
   33 joy=i
      if(a(ialow1))81,82,81
   82 write(6,299)joy,joy
      return
   81 write(6,297)joy,joy
      do 34 k=1,n1,ia
      jj=joy+k-1
      ji=joy+k
      if(joy+1-n)35,36,36
   35 write(6,296)
      return
   36 r=a(jj)
      a(jj)=a(ji)
   34 a(ji)=r
c     change sign of determinant
      detsc=-detsc
c
      if(m.eq.0) go to 8
      do 37 k=1,m1,ib
      jj=joy+k-1
      ji=joy+k
      r=b(jj)
      b(jj)=b(ji)
   37 b(ji)=r
      go to 8
   32 j1=j+1
      r=-a(j1)/a(ialow)
      i1=ialow
      do 13 k=j1,n1,ia
      a(k)=a(k)+a(i1)*r
   13 i1=i1+ia
c
c  loop to reset b has been modified, 25/1/1985.
c
      if(m.eq.0) go to 12
      i1=i
      i2=j-ialow+i+1
      do 14 k=1,m1,ib
      b(i2)=b(i2)+b(i1)*r
      i1=i1+ib
   14 i2=i2+ib
   12 continue
c
c
    5 continue
c
c
c     the matrix is now in triangular form
c     first set err=1.0d0
      err=1.0d0
c     if(any diagonal element of a is zero x cannot be solved for
      do 15 i=1,n
      idiag=(i-1)*ia+i
      err=err*a(idiag)
      if(err) 15,16,15
   16 write(6,299)i,i
      return
   15 continue
c     scale determinant
c....                   err=err*detsc ...... suppressed in this version
c
      if(m.eq.0) return
c     find solution to ax=b
      do 18 k=1,m
      ka=(n-1)*ia+n
      kb=(k-1)*ib+n
      b(kb)=b(kb)/a(ka)
      do 19 l=1,n2
      i=n-l
      r=0.0d0
      imax=i+1
      do 20 j=imax,n
      jj=i+n+1-j
      ja=(jj-1)*ia+i
      jb=(k-1)*ib+jj
   20 r=r+a(ja)*b(jb)
      la=(i-1)*ia+i
      lb=(k-1)*ib+i
   19 b(lb)=(b(lb)-r)/a(la)
   18 continue
c
      return
c
  295 format(///20h leq called with n =,i4)
  296 format(///48h the row cannot be changed with the row below it  ,
     .40h because it it the last row in the array   /
     .55h the solution, matrix x, cannot be found by this method  )
  297 format(5h1  a(,i4,1h,,i4,13h) equals zero  /
     .47h   try switching this row with the row below it   ,
     .48h and go back to statement number 8 and try again  )
  298 format(26h1  all the elements in row,i5,20h  are zero therefore,
     .55h the solution, matrix x, cannot be found by this method )
  299 format(50h1  the solution, matrix x, cannot be found by this  ,
     .57h method because there is a zero array element in the main ,
     .9h diagonal  / 30x,2ha(,i4,1h,,i4,8h) = zero  )
      end

      subroutine iscopy(n,ka,na,kb,nb)
      implicit double precision (a-h,o-z), integer (i-n)
c
c----- on the cray, this routine would be unnecessary and its task is
c----- performed by the cray assembly language routine 'scopy'.
c----- however, if one wanted to use this fact, all programs
c----- calling 'iscopy' would have to be maintained in two versions,
c----- for crays, and for all other machines. this would be too clumsy
c----- for a neglibible gain in cpu time.
c-----                                                    [wd, 7/8/00]
c
      dimension ka(1),kb(1)
      ia=1
      ib=1
      do 1 i=1,n
      kb(ib)=ka(ia)
      ia=ia+na
      ib=ib+nb
    1 continue
      return
      end 
      
