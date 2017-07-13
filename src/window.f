      subroutine window
      implicit double precision (a-h,o-z), integer (i-n)
c
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c***********************************************************************
c     select windows containing species with occupation fractions      c
c     above threshold.                                                 c
c                                                                      c
c --- modification July 1990:                                          c
c     in very degenerate matter (exeta close to machine zero) the      c
c     sphi(i,.) could become 0 and cause crash of the programme;       c
c     however, under these condition it is safe to use the result of   c
c     the previous call of window, and thus nothing is done here.      c
c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      if ( exeta .lt. 1.d-70 ) return
c
      if ( irho .gt. 1 ) go to 20
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set windows for first point on isotherm                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     loop over all chemical elements
c
      do 10 kchem = 1, nchem
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
      sum = ssum( nion(kchem), sn(is1), 1 )
      kk  = nucz(kchem)
c
      if ( kk .eq. 1 ) then
c
c-----------------------------------------------------------------------
c     hydrogen                                                         c
c-----------------------------------------------------------------------
c
c     account for hydrogen molecules
c
      sum = sum + sn(ish2) + sn(ish2p)
      do 1 is = is1, is2
      xfrac(is) = sn(is) / sum
      wwt(is) = cvmgt( 1.0d0, 0.0d0, xfrac(is) .ge. floor )
    1 continue
c
      if ( wwt(ispr) .lt. 0.1d0 ) wwt(ish2p) = 0.0d0
      iw1 = is1
      iw2 = is2
      sum = 2.0d0*(wwt(ish2)*xfrac(ish2) + wwt(ish2p)*xfrac(ish2p))
     .           + wwt(ishm)*xfrac(ishm) + wwt(ish  )*xfrac(ish  )
     .           + wwt(ispr)*xfrac(ispr)
c
      else
c
c-----------------------------------------------------------------------
c     all elements other than hydrogen                                 c
c-----------------------------------------------------------------------
c
      do 2 is = is1, is2
      xfrac(is) = sn(is) / sum
      wwt(is) = cvmgt( 1.0d0, 0.0d0, xfrac(is) .ge. floor )
    2 continue
c
c     find first species above threshold
      do 3 is = is1, is2
      iw1 = is
      if ( wwt(is) .gt. 0.0d0 ) go to 4
    3 continue
c
c     find last species above threshold
    4 do 5 is = is2, is1, -1
      iw2 = is
      if ( wwt(is) .gt. 0.0d0 ) go to 6
    5 continue
    6 sum = sdot( nion(kchem), wwt(is1), 1, xfrac(is1), 1 )
c
      end if
c
c-----------------------------------------------------------------------
c     set populations, renormalizing to current window                 c
c-----------------------------------------------------------------------
c
      iwin1(kchem) = iw1
      iwin2(kchem) = iw2
      sum = abun(kchem) * totn / sum
      do 7 is = is1, is2
      sn(is) = max( wwt(is)*xfrac(is)*sum, 1.0d-70 )
    7 continue
c
   10 continue
c
c     set unit weight for electrons
      wwt(ise) = 1.0d0
c
      go to 70
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set windows for subsequent points on isotherm                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 ckt = ck * t
      t32 = t**1.5d0
c
c     loop over all chemical elements
c
      do 60 kchem = 1, nchem
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
      kk  = nucz (kchem)
      sum = ssum( nion(kchem), sn(is1), 1)
c
c     initialize all weights for each element to zero
      do 21 is = is1, is2
      wwt(is) = 0.0d0
   21 continue
c
      if ( kk .eq. 1 ) then
c-----------------------------------------------------------------------
c     hydrogen                                                         c
c-----------------------------------------------------------------------
c
      sum = sum + sn(ish2) + sn(ish2p)
      do 30 is = is1, is2
      xfrac(is) = sn(is) / sum
   30 continue
c
c     set saha ratios (including degeneracy effects)
c
c     h2 dissociation
      sphi(1, kchem) =     ( gs(ish)**2 / gs(ish2) ) * t32
     .                   * ( z (ish)**2 / z (ish2) )
     .                   * exp( - ( 2.0d0*e0(ish)  - e0(ish2) ) / ckt )
c
c     h2+ dissociation
      sphi(2, kchem) =
     .            ( gs(ispr) * gs(ish) / gs(ish2p) ) * t32
     .         *  ( z (ispr) * z (ish) / z (ish2p) )
     .         *    exp( - ( e0(ispr) + e0(ish) - e0(ish2p) ) / ckt )
c
c     h- dissociation
      sphi(3, kchem) =  exeta * ( z (ish) / z (ishm) )
     .                * exp( - ( e0(ish) - e0(ishm) ) / ckt )
c
c     h atom ionization
      sphi(4, kchem) =  exeta * ( z (ispr) / z (ish) )
     .                * exp( - ( e0(ispr) - e0(ish) ) / ckt )
c
c     select species currently above threshold
      if ( xfrac(ish2 ) .ge. floor ) wwt(ish2 ) = 1.0d0
      if ( xfrac(ish2p) .ge. floor ) wwt(ish2p) = 1.0d0
      if ( xfrac(ishm ) .ge. floor ) wwt(ishm ) = 1.0d0
      if ( xfrac(ish  ) .ge. floor ) wwt(ish  ) = 1.0d0
      if ( xfrac(ispr ) .ge. floor ) wwt(ispr ) = 1.0d0
      if ( wwt  (ispr ) .lt. 0.1d0 ) wwt(ish2p) = 0.0d0
c
c     check whether h atoms should be added
      if ( wwt(ish) .gt. 0.1d0 ) go to 31
      frs = xfrac(ispr) / sphi(4, kchem)
      if ( frs .lt. floor ) go to 31
      xfrac(ish) = floor
      wwt  (ish) = 1.0d0
c
c     check whether protons should be added
   31 if ( wwt(ispr) .gt. 0.1d0 ) go to 32
      frs = xfrac(ish) * sphi(4, kchem)
      if ( frs .lt. floor .and. rholog .lt. -1.0d0 ) go to 32
c     at high densities add protons if neutral h is significantly
c     populated
      if ( xfrac(ish) .le. 1.0d-4 ) go to 32
      xfrac(ispr) = floor
      wwt  (ispr) = 1.0d0
c
c     check whether h- should be added
   32 if ( wwt(ishm) .gt. 0.1d0 ) go to 33
      frs = xfrac(ish) / sphi(3, kchem)
      if ( frs .lt. floor ) go to 33
      xfrac(ishm) = floor
      wwt  (ishm) = 1.0d0
c
c     check whether h2+ should be added
   33 if ( wwt(ish2p) .gt. 0.1d0 ) go to 34
      frs = sn(ish) * sphi(4, kchem) / ( sphi(2, kchem) * vol )
      if ( frs .lt. floor ) go to 34
      xfrac(ish2p) = floor
      wwt  (ish2p) = 1.0d0
      if ( wwt(ispr)   .lt. 0.1d0 ) wwt(ish2p) = 0.0d0
c
c     check whether h2 should be added
   34 if ( wwt(ish2) .gt. 0.1d0 ) go to 35
      frs = sn(ish) / ( sphi(1, kchem) * vol )
      if ( frs .lt. floor ) go to 35
      xfrac(ish2) = floor
      wwt  (ish2) = 1.0d0
c
   35 iw1 = is1
      iw2 = is2
      sum = 2.0d0*(wwt(ish2)*xfrac(ish2)+wwt(ish2p)*xfrac(ish2p))
     .           + wwt(ishm)*xfrac(ishm) + wwt(ish  )*xfrac(ish  )
     .           + wwt(ispr)*xfrac(ispr)
c
      else
c
c-----------------------------------------------------------------------
c     all elements other than hydrogen                                 c
c-----------------------------------------------------------------------
c
      do 40 is = is1, is2
      xfrac(is) = sn(is) / sum
   40 continue
c
c     set saha ratios (including degeneracy effects)
      do 41 jion = 1, nion(kchem) - 1
      is = jion + is1 - 1
      sr(jion, kchem) = exeta * ( z(is+1) /  z(is) )
     .                     * exp( - ( e0(is+1) - e0(is) ) / ckt )
   41 continue
c
c     find first species currently above threshold
      do 42 is = is1, is2
      iw1  = is
      if ( xfrac(is) .ge. floor ) go to 43
   42 continue
c
c     find last species currently above threshold
   43 do 44 is = is2, is1, -1
      iw2  = is
      if ( xfrac(is) .ge. floor ) go to 45
   44 continue
c
c     check whether a lower species should be added
   45 if ( iw1 .eq. is1 ) go to 46
      jion = iw1 - is1 + 1
      frs = xfrac(iw1) / sr(jion-1, kchem)
      if ( frs.lt.floor .and. xfrac(iw1).lt.1.0d-4 ) go to 46
      iw1 = iw1 - 1
      xfrac(iw1) = floor
c
c     check whether a higher species should be added
   46 if ( iw2 .eq. is2 ) go to 47
      jion = iw2 - is1 + 1
      frs = xfrac(iw2) * sr(jion, kchem)
      if ( frs.lt.floor .and. xfrac(iw2).lt.1.0d-4 ) go to 47
      iw2 = iw2 + 1
      xfrac(iw2) = floor
c
c     set window
   47 do 48 is = iw1, iw2
      wwt(is) = 1.0d0
   48 continue
c
      sum = sdot( nion(kchem), wwt(is1), 1, xfrac(is1), 1 )
c
      end if
c
c-----------------------------------------------------------------------
c     set populations, renormalizing to current window
c-----------------------------------------------------------------------
c
      iwin1(kchem) = iw1
      iwin2(kchem) = iw2
      do 50 is = is1, is2
      xfrac(is) = xfrac(is) / sum
      sn(is) = max( wwt(is)*xfrac(is)*abun(kchem)*totn, 1.0d-70 )
   50 continue
c
   60 continue
c
c     set unit weight for electrons
      wwt(ise) = 1.0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set up transformation matrix between reaction parameters and     c
c     occupation numbers                                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   70 call lambda
c
      return
      end


