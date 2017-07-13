      subroutine matgen
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     assemble matrix of linearized stoichiometric relations,          c
c     number conservation, and charge conservation equations           c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
c
      do 11 kchem = 1, nchem
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
      iw1 = iwin1(kchem)
      iw2 = iwin2(kchem)
c
      if ( nucz(kchem) .eq. 1 ) then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hydrogen                                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     initialize all equations as dummy equations
c
      do 1 is = is1, is2
      a(is, is) = 1.0d0
    1       continue
c
c           hydrogen molecule dissociation
c
      if ( wwt(ish2) .gt. 0.1d0 ) then
               dn(ish2) = 2.0d0 * dfdn(ish) - dfdn(ish2)
               do 2 js = 1, nspes
               a(ish2, js) = d2fdn2(ish2, js) - 2.0d0*d2fdn2(ish, js)
    2          continue
      end if
c
c           h2+ molecule dissociation
c
      if ( wwt(ish2p) .gt. 0.1d0 ) then
               dn(ish2p) = dfdn(ispr) + dfdn(ish) - dfdn(ish2p)
               do 3 js = 1, nspes
               a(ish2p, js) = d2fdn2(ish2p, js) - d2fdn2(ish , js)
     .                                          - d2fdn2(ispr, js)
    3          continue
      end if
c
c           h minus dissociation
c
      if ( wwt(ishm) .gt. 0.1d0 ) then
               dn(ishm) = dfdn(ise) + dfdn(ish) - dfdn(ishm)
               do 4 js = 1, nspes
               a(ishm, js) = d2fdn2(ishm, js) - d2fdn2(ish, js)
     .                                        - d2fdn2(ise, js)
    4          continue
      end if
c
c           h ionization
c
      isnc = ish
      if ( wwt(ispr) .lt. 0.1d0 ) go to 55
      isnc = ispr
      if ( wwt(ish ) .lt. 0.1d0 ) go to 55
            dn(ish) = dfdn(ise) + dfdn(ispr) - dfdn(ish)
            do 5 js = 1, nspes
            a(ish, js) = d2fdn2(ish, js) - d2fdn2(ispr, js)
     .                                   - d2fdn2(ise , js)
    5       continue
c
c           hydrogen number conservation
c
   55       dn(isnc) = abun(kchem) * totn
     .             - 2.0d0*( wwt(ish2)*sn(ish2) + wwt(ish2p)*sn(ish2p) )
     .             - wwt(ishm)*sn(ishm) - wwt(ish)*sn(ish)
     .             - wwt(ispr)*sn(ispr)
      a(isnc, ish2 ) = 2.0d0 * wwt(ish2 )
      a(isnc, ish2p) = 2.0d0 * wwt(ish2p)
      a(isnc, ishm ) = 1.0d0 * wwt(ishm )
      a(isnc, ish  ) = 1.0d0 * wwt(ish  )
      a(isnc, ispr ) = 1.0d0 * wwt(ispr )
      else
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     all elements other than hydrogen                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c           ionization balance equations
c
            do 7 is = iw1, iw2 - 1
            dn(is) = dfdn(ise) + dfdn(is+1) - dfdn(is)
            do 6 js = 1, nspes
            a(is, js) = d2fdn2(is, js) - d2fdn2(is+1, js)
     .                                 - d2fdn2(ise , js)
    6       continue
    7       continue
c
c           number conservation equation
c
            dn(iw2) = abun(kchem) * totn
     .              - sdot( nion(kchem), sn(is1), 1, wwt(is1), 1 )
            do 8 js = iw1, iw2
            a(iw2, js) = 1.0d0
    8       continue
c
c           dummy equations
c
            do 9 is = is1, iw1 - 1
            a(is, is) = 1.0d0
    9       continue
            do 10 is = iw2 + 1, is2
            a(is, is) = 1.0d0
   10       continue
c
       end if
   11  continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     total charge conservation                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dn(ise) = sn(ise) - sdot( nspes - 1, zs, 1, sn, 1 )
      do 12 js = 1, nspes
      a(ise, js) = wwt(js) * zs(js)
   12 continue
c
c      if(rholog.gt.1.0) then
c
c      write(6,7789) (dn(jj),jj=1,nspes)
c 7789 format(' dn:in (1-ns)',/(1x,1p5g15.6))
c      write(6,7749) (dfdn(jj),jj=1,nspes)
c 7749 format(' dfdn  (1-ns)',/(1x,1p5g15.6))
c
c      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     transform matrix to system yielding fractional changes (dn/n)    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 20 js = 1, nspes
      call sscal( nspes, sn(js), a(1, js), 1 )
   20 continue
c
      return
      end
