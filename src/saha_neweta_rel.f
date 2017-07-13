      subroutine saha1
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     saha equilibrium for mixture with hydrogen                       c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      data    nlow,     nc,  niter
     .    /     15,     20,     20/
      data    tolm,   dmax,   conv
     .    /  0.4d0,  0.8d0,  1.d-9/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     evaluate saha phis. for the present ignore electron degeneracy   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call set phi
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     estimate upper and lower bounds on electron number               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     upper bound: complete ionization
      elr = 0.d0
      do 100 kchem = 1, nchem
      elr = elr + abun(kchem) * dfloat( nucz(kchem) )
  100 continue
      elr = totn * elr
c
c     lower bound: ionization of "representative metal" and/or hydrogen
      itlow = 0
c          "representative metal"
      rat = 4.0d0 * abunm * totn / ( phim * vol )
      if ( rat .ge. 0.001d0 ) then
      elm = 0.5d0 * phim * vol * ( dsqrt( 1.0d0 + rat ) - 1.0d0 )
      else
      elm = abunm * totn * ( 1.0d0 - 0.250d0 * rat * ( 1.0d0
     .            - 0.5d0 * rat  * ( 1.0d0 - 0.625d0 * rat )))
      end if
c          hydrogen
      kh = kz(1)
      rat = 4.0d0 * abun(kh) * totn / ( sphi(4, kh) * vol )
      if ( rat .ge. 0.001d0 ) then
      elh = 0.5d0 * sphi(4, kh) * vol * ( dsqrt( 1.0d0 + rat ) - 1.0d0 )
      else
      elh = abun(kh) * totn * ( 1.0d0 - 0.250d0 * rat * ( 1.0d0
     .               - 0.5d0 * rat  * ( 1.0d0 - 0.625d0 * rat )))
      end if
c          choose larger of two contributions
      ell = max( elm, elh, 1.0d-70 )
c
c     given electron number, compute hydrogen molecule/atom/ion fraction
    1 c1  = 1.0d0+ell / ( sphi(3, kh) * vol ) + sphi(4, kh) * vol / ell
      c2  = 2.0d0*abun(kh) * totn * ( 1.0d0/( sphi(1, kh) * vol )
     .                              + sphi(4, kh) / (ell * sphi(2, kh)))
      c3  = 0.5d0 * c1 / c2
      rat = 4.0d0 * c2 / c1**2
      if ( rat .ge. 0.001d0 ) then
      sf(4, kh) = c3 * ( dsqrt( 1.0d0 + rat ) - 1.0d0 )
      else
      sf(4, kh) = ( 1.0d0 - 0.250d0 * rat * ( 1.0d0 - 0.5d0 * rat *
     .                  ( 1.0d0 - 0.625d0 * rat ) ) ) / c1
      end if
      sf(5, kh) = sf(4, kh) * sphi(4, kh) * vol / ell
      sf(3, kh) = sf(4, kh) * ell / ( sphi(3, kh) * vol )
      sf(2, kh) = sf(4, kh) * sf(5, kh)
     .                      * abun(kh) * totn / ( sphi(2, kh) * vol )
      sf(1, kh) = sf(4, kh) * sf(4, kh)
     .                      * abun(kh) * totn / ( sphi(1, kh) * vol )
c     get new estimate of electron number from charge conservation
      call fjk( ell )
      sum = 0.0d0
      do 2 kchem = 1, nchem
      is1 = ichm1(kchem)
      sum = sum
     .    + abun(kchem) *sdot( nion(kchem), sf(1, kchem), 1, zs(is1), 1)
    2 continue
      elec = totn * sum
c
c     check for fortuitous convergence
      if ( elec .eq. ell ) go to 12
c
c     check that "lower bound" really is. if not, try again.
      if ( elec .gt. ell ) go to 4
      elr = ell
      ell = 0.1d0 * ell
      itlow = itlow + 1
      if ( itlow .le. nlow ) go to 1
c
c     failure to get bound
      write  ( 0, 3 ) tlog, rholog, ell, elr, elec
    3 format ('saha1 unable to obtain lower bound',/ 2f10.3, 1p3e10.3 )
      stop 'saha1:1'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     narrow range by geometric averaging                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    4 itc = 0
c
c     take geometric average of upper and lower bounds
    5 elc = dsqrt( ell * elr )
c
c     given electron number, compute hydrogen molecule/atom/ion fraction
      c1  = 1.0d0+elc / ( sphi(3, kh) * vol ) + sphi(4, kh) * vol / elc
      c2  = 2.0d0* abun(kh) * totn * ( 1.0d0/ ( sphi(1, kh) * vol )
     .                              + sphi(4, kh) / (elc * sphi(2, kh)))
      c3  = 0.5d0 * c1 / c2
      rat = 4.0d0 * c2 / c1**2
      if ( rat .ge. 0.001d0 ) then
      sf(4, kh) = c3 * ( dsqrt( 1.0d0 + rat ) - 1.0d0 )
      else
      sf(4, kh) = ( 1.0d0 - 0.250d0 * rat * ( 1.0d0 - 0.5d0 * rat *
     .                  ( 1.0d0 - 0.625d0 * rat ) ) ) / c1
      end if
      sf(5, kh) = sf(4, kh) * sphi(4, kh) * vol / elc
      sf(3, kh) = sf(4, kh) * elc / ( sphi(3, kh) * vol )
      sf(2, kh) = sf(4, kh) * sf(5, kh)
     .                      * abun(kh) * totn / ( sphi(2, kh) * vol )
      sf(1, kh) = sf(4, kh) * sf(4, kh)
     .                      * abun(kh) * totn / ( sphi(1, kh) * vol )
c
c     get new estimate of electron number from charge conservation
      call fjk( elc )
      sum = 0.0d0
      do 6 kchem = 1, nchem
      is1 = ichm1(kchem)
      sum = sum
     .    + abun(kchem) *sdot( nion(kchem), sf(1, kchem), 1, zs(is1), 1)
    6 continue
      elec = totn * sum
c
c     check for fortuitous convergence
      if ( elec .eq. elc ) go to 12
c
c     check if bounds are close enough to switch to newton-raphson
      if ( 2.0d0*dabs( elec - elc )/( elec + elc ) .le. tolm ) go to 8
c
c     replace a bound
      if ( elec .gt. elc ) ell = max( elc, 1.0d-70 )
      if ( elec .lt. elc ) elr = elc
      itc = itc + 1
      if ( itc .le. nc ) go to 5
c
c     failure to narrow range sufficiently
      write  ( 0, 7 ) tlog, rholog, ell, elr, elc, elec
    7 format ('saha1 unable to narrow range by geometric averaging'/
     .        2f10.3, 1p4e10.3 )
      stop 'saha1:2'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration on charge conservation and hydrogen     c
c     number conservation equations.                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    8 iter = 0
c
c     initialize derivatives of hydrogen occupation numbers wrt electron
c     number
      dsfde(1, kh) =    0.d0
      dsfde(2, kh) = -sf(2, kh) / elc
      dsfde(3, kh) = -sf(3, kh) / elc
      dsfde(4, kh) =    0.d0
      dsfde(5, kh) = -sf(5, kh) / elc
c
c     linearized hydrogen number conservation equation
    9 a11 =   4.0d0*(sf(1,kh)+sf(2,kh)) + sf(3,kh) + sf(4,kh) + sf(5,kh)
      a12 =  -2.0d0*          sf(2,kh)  + sf(3,kh)            - sf(5,kh)
      b1  =   1.0d0
     .      - 2.0d0*(sf(1,kh)+sf(2,kh)) - sf(3,kh) - sf(4,kh) - sf(5,kh)
c
c     linearized charge conservation equation
      sum  = 0.0d0
      sum1 = 0.0d0
      do 10 kchem = 1, nchem
      is1  = ichm1(kchem)
      sum  = sum
     .     + abun(kchem)*sdot(nion(kchem), sf   (1,kchem),1, zs(is1),1)
      sum1 = sum1
     .     + abun(kchem)*sdot(nion(kchem), dsfde(1,kchem),1, zs(is1),1)
   10 continue
      sum1 = elc * sum1
      a21 =  totn * abun(kh)*(2.0d0*sf(2, kh) - sf(3, kh) + sf(5, kh) )
      a22 = - elc + totn * sum1
      b2  =   elc - totn * sum
c
c     solve for fractional change in neutral hydrogen occupation fractio
c     and in electron number
      det =  a11 * a22 - a12 * a21
      dfh = ( b1 * a22 - b2  * a12 ) / det
      de  = ( b2 * a11 - b1  * a21 ) / det
c
c     limit changes
      admx = max( dabs( dfh ), dabs( de ) )
      if ( admx .gt. dmax ) then
      dfh = dfh * dmax / admx
      de  = de  * dmax / admx
      end if
c
c     apply them and update charge conservation
      sf(4, kh) = ( 1.0d0 + dfh ) * sf(4, kh)
      elc       = ( 1.0d0 + de  ) * elc
c
c     update occupation fractions and their derivatives wrt electron num
      call fjk( elc )
      sf(5, kh) = sf(4, kh) * sphi(4, kh) * vol / elc
      sf(3, kh) = sf(4, kh) * elc / ( sphi(3, kh) * vol )
      sf(2, kh) = sf(4, kh) * sf(5, kh)
     .                      * abun(kh) * totn / ( sphi(2, kh) * vol )
      sf(1, kh) = sf(4, kh) * sf(4, kh)
     .                      * abun(kh) * totn / ( sphi(1, kh) * vol )
      dsfde(2, kh) = - sf(2, kh) / elc
      dsfde(3, kh) = - sf(3, kh) / elc
      dsfde(5, kh) = - sf(5, kh) / elc
c
c     check for convergence
      elec = elc
      if ( admx .le. conv ) go to 12
      iter = iter + 1
      if ( iter .le. niter ) go to 9
c
c     failure to converge
      write  ( 0, 11 ) tlog, rholog, sf(4, kh), dfh, elc, de
   11 format ('nonconvergence of nondegenerate newton-raphson in saha1'
     .        / 2f10.3, 1p4e10.3 )
      stop 'saha1:3'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     now allow for electron degeneracy                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   12 t32  = t**1.5d0
      conf = vol * t32 / ( ceta * totn )
ccccccccccccccccccccccccccccccccccc temp ccccccccccccccccccccccccccccc
ccc      if(.true.) goto 20
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      iter = 0
c
c     use asymptotic expression for first guess for eta
      eta = log ( elec / ( 2.0d0 * gs(ise) * vol * t32 ) )
c     now find accurate value
      sn(ise) = elec
ccc      write(6,*)'saha1,guessed eta,elec: ',eta,elec
      call new eta
ccc      write(6,*)'saha1,exact eta: ',eta
c
c     convert phis into form to use electron degeneracy parameter as bas
c     variable instead of electron number
      sphi(3, kh) = sphi(3, kh) / ( 2.0d0 * gs(ise) * t32 )
      sphi(4, kh) = sphi(4, kh) / ( 2.0d0 * gs(ise) * t32 )
      do 14 kchem = 1, nchem
      if ( kchem .eq. kh ) go to 14
      do 13 jion  = 1, nion(kchem) - 1
      sphi(jion, kchem) = sphi(jion, kchem) / ( 2.0d0*gs(ise) * t32 )
   13 continue
   14 continue
c
c     initialize occupation fractions and their derivatives wrt eta
      call fjkd
      sf(5, kh) = sf(4, kh) *   sphi(4, kh) * exeta
      sf(3, kh) = sf(4, kh) / ( sphi(3, kh) * exeta )
      sf(2, kh) = sf(4, kh) * sf(5, kh)
     .                      * abun(kh) * totn / ( sphi(2, kh) * vol )
      sf(1, kh) = sf(4, kh) * sf(4, kh)
     .                      * abun(kh) * totn / ( sphi(1, kh) * vol )
      dsfde(2, kh) = - sf(2, kh)
      dsfde(3, kh) =   sf(3, kh)
      dsfde(5, kh) = - sf(5, kh)
c
c     linearized hydrogen number conservation equation
   15 a11 =   4.0d0*(sf(1,kh)+sf(2,kh)) + sf(3,kh) + sf(4,kh) + sf(5,kh)
      a12 = - 2.0d0*          sf(2,kh)  + sf(3,kh)            - sf(5,kh)
      b1  =   1.0d0
     .      - 2.0d0*(sf(1,kh)+sf(2,kh)) - sf(3,kh) - sf(4,kh) - sf(5,kh)
c
c     linearized charge conservation equation
      sum  = 0.0d0
      sum1 = 0.0d0
      do 16 kchem = 1, nchem
      is1 = ichm1(kchem)
      sum = sum
     .     + abun(kchem)*sdot(nion(kchem), sf   (1, kchem),1, zs(is1),1)
      sum1 = sum1
     .     + abun(kchem)*sdot(nion(kchem), dsfde(1, kchem),1, zs(is1),1)
   16 continue
c
c ---- switch nonrel/rel
cccc     beta=0.d0
      beta  = ck*t/cme/cc**2
c
      f12 = r_fd(0,0.5d0,eta,beta)
      f32 = r_fd(0,1.5d0,eta,beta)
      fe12= r_fd(1,0.5d0,eta,beta)
      fe32= r_fd(1,1.5d0,eta,beta)

      a21 = abun(kh) * (2.0d0*sf(2, kh) - sf(3, kh) + sf(5, kh) )
      a22 = - conf * ( fe12 + beta * fe32 ) + sum1
      b2  =   conf * ( f12  + beta * f32  ) - sum
c
c  Non-relativistic
c
c      a22 = - 0.5d0* conf * fd(1) + sum1
c      b2  =          conf * fd(2) - sum
c
ccc      write(6,*)'degen: conf,fd2,sum,b2 ',conf,fd(2),sum,b2
c
c     solve for fractional change in neutral hydrogen occupation fractio
c     and in electron degeneracy parameter
      det   =  a11 * a22 - a12 * a21
      dfh   = ( b1 * a22 - b2  * a12 ) / det
      deta  = ( b2 * a11 - b1  * a21 ) / det
ccc      write(6,*)'det,dfh,deta,a11,a12,a21,a22,b1,b2',
ccc     . det,dfh,deta,a11,a12,a21,a22,b1,b2
c
c     limit changes
      admx  = max( dabs( dfh ), dabs( deta ) )
      if ( admx .gt. dmax ) then
      dfh   = dfh   * dmax / admx
      deta  = deta  * dmax / admx
      end if
c
c     apply them and update electron number
      sf(4, kh) = ( 1.0d0 + dfh ) * sf(4, kh)
      eta = eta + deta
      exeta = exp( - eta )
      f12 = r_fd(0,0.5d0,eta,beta)
      f32 = r_fd(0,1.5d0,eta,beta)
      fe12= r_fd(1,0.5d0,eta,beta)
      fe32= r_fd(1,1.5d0,eta,beta)
      elec = totn * conf * ( f12 + beta * f32 )

c
c  Non-relativistic
c
c      call ferdir( eta, fd )
c      elec = totn * conf * fd(2)
c............... temporary ................
ccc      write(6,*) 'iter,saha1: eta,elec,deta',eta,elec
c..........................................
c
c     update occupation fractions and their derivatives wrt eta
      call fjkd
      sf(5, kh) = sf(4, kh) *   sphi(4, kh) * exeta
      sf(3, kh) = sf(4, kh) / ( sphi(3, kh) * exeta )
      sf(2, kh) = sf(4, kh) * sf(5, kh)
     .                      * abun(kh) * totn / ( sphi(2, kh) * vol )
      sf(1, kh) = sf(4, kh) * sf(4, kh)
     .                      * abun(kh) * totn / ( sphi(1, kh) * vol )
      dsfde(2, kh) = - sf(2, kh)
      dsfde(3, kh) =   sf(3, kh)
      dsfde(5, kh) = - sf(5, kh)
c
c     check for convergence
      if ( admx .le. conv  ) go to 20
      iter = iter + 1
      if ( iter .le. niter ) go to 15
c
c     failure to converge
      write  ( 0, 17 ) tlog, rholog, eta, deta, sf(4, kh), dfh, elec
   17 format ('nonconvergence of degen. charge conservation in saha1'
     .        / 3f10.3, 1p4e10.3 )
      stop 'saha1:4'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convergence. calculate all occupation numbers.                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 do 22 kchem = 1, nchem
      is1 = ichm1(kchem)
      do 21 jion = 1, nion(kchem)
      sn( jion + is1 - 1 ) =max(totn * abun(kchem) * sf(jion, kchem),
     .                            1.0d-70 )
   21 continue
   22 continue
      sn(ise) = max( elec, 1.0d-70 )
c
c     update all variables associated with degeneracy parameter
      call new eta
c
      return
      end
      subroutine saha2
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     saha equilibrium for mixture without hydrogen                    c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      data    nlow,     nc,  niter
     .    /     15,     20,     20/
      data    tolm,   dmax,   conv
     .    /  0.4d0,  0.8d0,  1.d-9/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     evaluate saha phis. for the present ignore electron degeneracy   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call set phi
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     estimate upper and lower bounds on electron number               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     upper bound: complete ionization
      elr = 0.d0
      do 100 kchem = 1, nchem
      elr = elr + abun(kchem) * dfloat( nucz(kchem) )
  100 continue
      elr = totn * elr
c
c     lower bound: ionization of "representative metal" only
      itlow = 0
      rat = 4.0d0 * abunm * totn / ( phim * vol )
      if ( rat .ge. 0.001d0 ) then
      ell = 0.5d0 * phim * vol * ( dsqrt( 1.0d0 + rat ) - 1.0d0 )
      else
      ell = abunm * totn * ( 1.0d0-0.250d0*rat * ( 1.0d0-0.5d0*rat
     .                         * ( 1.0d0-0.625d0*rat )))
      end if
      ell = max( ell, 1.0d-70 )
c
c     get new estimate of electron number from charge conservation
    1 call fjk( ell )
      sum = 0.0d0
      do 2 kchem = 1, nchem
      is1 = ichm1(kchem)
      sum = sum
     .    + abun(kchem)*sdot( nion(kchem), sf(1, kchem), 1, zs(is1) , 1)
    2 continue
      elec = totn * sum
c
c     check for fortuitous convergence
      if ( elec .eq. ell ) go to 12
c
c     check that "lower bound" really is. if not, try again.
      if ( elec .gt. ell ) go to 4
      elr = ell
      ell = 0.1d0 * ell
      itlow = itlow + 1
      if ( itlow .le. nlow ) go to 1
c
c     failure to get bound
      write  ( 0, 3 ) tlog, rholog, ell, elr, elec
    3 format ('saha2 unable to obtain lower bound'/ 2f10.3, 1p3e10.3 )
      stop 'saha2:1'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     narrow range by geometric averaging                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    4 itc = 0
c
c     take geometric average of upper and lower bounds
    5 elc = dsqrt( ell * elr )
c
c     get new estimate of electron number from charge conservation
      call fjk( elc )
      sum = 0.0d0
      do 6 kchem = 1, nchem
      is1 = ichm1(kchem)
      sum = sum
     .    + abun(kchem) * sdot(nion(kchem), sf(1, kchem), 1, zs(is1), 1)
    6 continue
      elec = totn * sum
c
c     check for fortuitous convergence
      if ( elec .eq. elc ) go to 12
c
c     check if bounds are close enough to switch to newton-raphson
      if ( 2.0d0*dabs(elec - elc)/( elec + elc ) .le. tolm ) go to 8
c
c     replace a bound
      if ( elec .gt. elc ) ell = max( elc, 1.0d-70 )
      if ( elec .lt. elc ) elr = elc
      itc = itc + 1
      if ( itc .le. nc ) go to 5
c
c     failure to narrow range sufficiently
      write  ( 0, 7 ) tlog, rholog, ell, elr, elc, elec
    7 format ('saha2 unable to narrow range by geometric averaging'/
     .        2f10.3, 1p4e10.3 )
      stop 'saha2:2'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration on charge conservation equation.        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    8 iter = 0
c     linearized charge conservation equation
    9 sum  = 0.0d0
      sum1 = 0.0d0
      do 10 kchem = 1, nchem
      is1  = ichm1(kchem)
      sum  = sum
     .     + abun(kchem)*sdot(nion(kchem), sf   (1,kchem), 1, zs(is1),1)
      sum1 = sum1
     .     + abun(kchem)*sdot(nion(kchem), dsfde(1,kchem), 1, zs(is1),1)
   10 continue
      sum1 = elc * sum1
      rhs = elc - totn * sum
      de  = rhs / ( totn * sum1 - elc )
c
c     limit change
      ade = dabs( de )
      if ( ade .gt. dmax ) de = de * dmax / ade
c
c     apply it and update charge occupation fractions and their derivati
c     wrt electron number
      elc = ( 1.0d0 + de ) * elc
      elec = elc
      call fjk( elc )
c
c     check for convergence
      if ( ade .le. conv ) go to 12
      iter = iter + 1
      if ( iter .le. niter ) go to 9
c
c     failure to converge
      write  ( 0, 11 ) tlog, rholog, elc, de
   11 format ('nonconvergence of nondegen. charge conservation in saha2'
     .        / 2f10.3, 1p2e10.3 )
      stop 'saha2:3'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     now allow for electron degeneracy                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   12 t32  = t**1.5d0
      conf = vol * t32 / ( ceta * totn )
      iter = 0
c
c     use asymptotic expression for first guess for eta
      eta = log ( elec / ( 2.0d0 * gs(ise) * vol * t32 ) )
c
c     find accurate value for eta
      sn(ise) = elec
      call new eta
c
c     convert phis into form to use electron degeneracy parameter as bas
c     variable instead of electron number
      do 14 kchem = 1, nchem
      do 13 jion  = 1, nion(kchem) - 1
      sphi(jion, kchem) = sphi(jion, kchem) / ( 2.0d0*gs(ise) * t32 )
   13 continue
   14 continue
c
c     initialize occupation fractions and their derivatives wrt electron
      call fjkd
c
c     apply newton-raphson to charge conservation equation, now allowing
c     degeneracy
   15 sum  = 0.0d0
      sum1 = 0.0d0
      do 16 kchem = 1, nchem
      is1  = ichm1(kchem)
      sum  = sum
     .     + abun(kchem)*sdot(nion(kchem), sf   (1,kchem), 1, zs(is1),1)
      sum1 = sum1
     .     + abun(kchem)*sdot(nion(kchem), dsfde(1,kchem), 1, zs(is1),1)
   16 continue
c
c  Non-relativistic
c
c      rhs   = conf * fd(2) - sum
c      deta  = rhs / ( sum1 - 0.5d0 * conf * fd(1) )
c
c ---- switch nonrel/rel
cccc     beta=0.d0
      beta = ck*t/cme/cc**2
c
      f12   = r_fd(0,0.5d0,eta,beta)
      f32   = r_fd(0,1.5d0,eta,beta)
      fe12  = r_fd(1,0.5d0,eta,beta)
      fe32  = r_fd(1,1.5d0,eta,beta)

      rhs   = conf * ( f12 + beta * f32 ) - sum
      deta  = rhs / ( sum1 - conf * ( fe12 + beta * fe32 ) )
c
c     limit change
      adeta = dabs( deta )
      if ( adeta .gt. dmax ) deta = deta * dmax / adeta
c
c     apply it and update electron number
      eta = eta + deta
      exeta = exp( - eta )
      f12   = r_fd(0,0.5d0,eta,beta)
      f32   = r_fd(0,1.5d0,eta,beta)
      fe12  = r_fd(1,0.5d0,eta,beta)
      fe32  = r_fd(1,1.5d0,eta,beta)

c
c  Non-relativistic
c
c      call ferdir( eta, fd )
c      elec = totn * conf * fd(2)

      elec = totn * conf * ( f12 + beta * f32 )
c
c     update occupation fractions and their derivatives wrt electron num
      call fjkd
c
c     check for convergence
      if ( adeta .le. conv ) go to 20
      iter = iter + 1
      if ( iter .le. niter ) go to 15
c
c     failure to converge
      write  ( 0, 17 ) tlog, rholog, eta, deta, elec
   17 format ('nonconvergence of degen. charge conservation in saha2'
     .        / 3f10.3, 1p2e10.3 )
      stop 'saha2:4'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convergence. calculate all occupation numbers.                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 do 22 kchem = 1, nchem
      is1 = ichm1(kchem)
      do 21 jion = 1, nion(kchem)
      sn( jion + is1 - 1 ) =max(totn * abun(kchem) * sf(jion, kchem),
     .                            1.0d-70 )
   21 continue
   22 continue
      sn(ise) = elec
c
c     update all variables associated with degeneracy parameter
      call new eta
c
      return
      end
c
      subroutine new eta
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     from given ( t, ne) calculate degeneracy parameter and its       c
c     derivatives. evaluate theta and its derivatives                  c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      data   niter,   conv
     .    /    20, 1.d-10 /
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      snise = max(sn(ise),1.d-70)
      rhs = ceta * snise / ( t**1.5d0 * vol )
c
c........... initialization to a few per cent accuracy
c........... (see ref. in dappen, astron. astrphys., 1980)
c
      if(rhs.lt.1.d-6) then
            eta = log(rhs)
      else if(rhs.lt.100.d0) then
            g   =rhs
            g2  =rhs*rhs
            et0 =g+(g2/(4.45+2.257*g+g2))*exp((1.32934*g)**0.66666667)
            eta = log(et0)
      else
            eta = (1.32934*rhs)**0.66666667
      end if
      eta = eta + 0.120782
c
ccc      write(6,*)'eta,rhs ini: ',eta,rhs
c
c ---- switch nonrel/rel
cccc     beta=0.d0
      beta  = ck*t/cme/cc**2
c
      f12   = r_fd(0,0.5d0,eta,beta)
      f32   = r_fd(0,1.5d0,eta,beta)
      fe12  = r_fd(1,0.5d0,eta,beta)
      fe32  = r_fd(1,1.5d0,eta,beta)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iter  = 0
c
    1 deta  = ( rhs - f12 - beta * f32 ) / ( fe12 + beta * fe32 )
      eta = eta + deta
      f12   = r_fd(0,0.5d0,eta,beta)
      f32   = r_fd(0,1.5d0,eta,beta)
      fe12  = r_fd(1,0.5d0,eta,beta)
      fe32  = r_fd(1,1.5d0,eta,beta)
      if( dabs(deta) .le. conv ) go to 3
      iter = iter + 1
      if( iter .le. niter ) go to 1
c
c     failure to converge
      write  ( 0, 2 ) tlog, rholog, snise, eta
    2 format ( ' nonconvergence of degeneracy parameter ' /
     .         ' log t =',f5.2,' log rho =',f6.2,' ne =',1pe10.3,
     .         ' eta   =',e10.3 )
      stop 'neweta'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     convergence                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c -------------- note: exeta passed through common to be used elsewhere!
c
    3 exeta = exp( - eta )
c
      f52    = r_fd(0,2.5d0,eta,beta)
      fe52   = r_fd(1,2.5d0,eta,beta)
      f2e12  = r_fd(3,0.5d0,eta,beta)
      f2e32  = r_fd(3,1.5d0,eta,beta)
      f2e52  = r_fd(3,2.5d0,eta,beta)
      f2e72  = r_fd(3,3.5d0,eta,beta)
      f3e12  = r_fd(6,0.5d0,eta,beta)
      f3e32  = r_fd(6,1.5d0,eta,beta)
      f3e52  = r_fd(6,2.5d0,eta,beta)

***************************************************************
*  perfect cancellation appears in dthet etc. when eta is     *
*  small, but not harmful for final results up to 12th digits *
*                                                 z.g. (11/99)*
***************************************************************

      thet   = ( fe12 + beta * fe32 ) / ( f12 + beta * f32 )
      dthet  = thet*((f2e12+beta*f2e32)/(fe12+beta*fe32)-thet)
      d2thet = thet*((f3e12+beta*f3e32)/(fe12+beta*fe32)-3.0d0*dthet
     .         -thet**2)

      game    = ( fe32 + beta * fe52 ) / ( f32 + beta * f52 )
      dgamde  = game*((f2e32+beta*f2e52)/(fe32+beta*fe52)-game)
      d2gamde2= game*((f3e32+beta*f3e52)/(fe32+beta*fe52)-3.d0*dgamde
     .          -game**2)
c
      detdv   = - 1.0d0 / ( thet * vol )
      detdt   = - (fe32+beta*fe52)/(fe12+beta*fe32) / t
      detdn   =   1.0d0 / ( thet * snise )
      d2etdn2 = - detdn**2 * (f2e12+beta*f2e32)/(fe12+beta*fe32)
      d2etdnt = - detdn * detdt *(thet+dthet/thet-game-dgamde/game)
      d2etdnv = - dthet * detdn * detdv / thet
      d2etdtv = - detdt * detdv *(thet+dthet/thet-game-dgamde/game)
      d2etdt2 = - (1.0d0/t+detdt*(thet+dthet/thet-game-dgamde/game))
     .          * detdt+beta*((fe52-f2e72)/t-f2e52*detdt)/(f12+beta
     .          * f32)/2.d0/t/thet
      d2etdv2 =   ( thet - dthet / thet ) * detdv**2

      if(iter.gt.5) write (6,*)' bad convergence in neweta: iter,eta',
     . ',tlog,rholog= ',iter,eta,tlog,rholog
c
ccc      write(6,*) 'eta,fd ',eta,fd
      return
      end
