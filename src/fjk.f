      subroutine fjk( xelec )
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c    calculate r(j, k), p(j, k), s(k), f(j, k) and their derivatives   c
c    wrt electron number for all elements other than hydrogen, using   c
c    electron number as the basic variable                             c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over all chemical elements                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 20 kchem = 1, nchem
      if ( nucz(kchem) .eq. 1 ) go to 20
      nr = nion(kchem) - 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     evaluate r(jion, kchem) and jmax(kchem)                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 1 jion = 1, nr
      sr(jion, kchem) = sphi(jion, kchem) * vol / xelec
    1 continue
c
c     now figure out ion stage with maximum occupation number
c
      jmax = 0
c
c     hypothesis (1): neutral atom
      if ( sr(   1, kchem) .gt. 1.0 d0 ) go to 3
      do 2 jion = 2, nr
      if ( sr(jion, kchem) .ge. 1.0 d0 ) go to 3
    2 continue
c     accept neutral atom as dominant ion state
      jmax = 1
      go to 7
c
c     hypothesis (2): highest ion
    3 if ( sr(nr  , kchem) .lt. 1.0 d0) go to 5
      do 4 jion = 1, nr - 1
      if ( sr(jion, kchem) .le. 1.0 d0) go to 5
    4 continue
c     accept highest ion as dominant ion state
      jmax = nr + 1
      go to 7
c
c     hypothesis (3): some intermediate ion state
    5 do 6 jion = 1, nr - 1
      if ( sr(jion    , kchem) .ge. 1.0d0 .and.
     .     sr(jion + 1, kchem) .le. 1.0d0 ) jmax = jion + 1
    6 continue
c
c     if the process failed, arbitrarily choose the middle ion
      if ( jmax .eq. 0 ) jmax = ( nr + 1 ) / 2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate p(jion, kchem), s(kchem), f(jion, kchem) and their     c
c     derivatives wrt to electron number                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    7 sp(jmax, kchem) = 1.0d0
      do 8 jion = jmax + 1, nion(kchem)
      sp(jion, kchem) = sp(jion - 1, kchem) * sr(jion - 1, kchem)
    8 continue
      do 9 jion = jmax - 1, 1, -1
      sp(jion, kchem) = sp(jion + 1, kchem) / sr(jion, kchem)
    9 continue
      do 10 jion = 1, nion(kchem)
      dspde(jion, kchem) = dfloat(jmax - jion ) * sp(jion, kchem)/xelec
   10 continue
c
      ss   (kchem) = ssum( nion(kchem), sp   (1, kchem), 1 )
      dssde(kchem) = ssum( nion(kchem), dspde(1, kchem), 1 )
c
      do 11 jion = 1, nion(kchem)
      sf   (jion, kchem) = sp(jion, kchem) / ss(kchem)
      dsfde(jion, kchem) =
     .              ( dspde(jion, kchem) * ss   (kchem)
     .              - sp   (jion, kchem) * dssde(kchem) ) / ss(kchem)**2
   11 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close loop over chemical elements                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 continue
c
      return
      end
c
      subroutine fjkd
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c    calculate r(j, k), p(j, k), s(k), f(j, k) and their derivatives   c
c    wrt electron degeneracy parameter for all elements other than     c
c    hydrogen                                                          c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over all chemical elements                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 20 kchem = 1, nchem
      if ( nucz(kchem) .eq. 1 ) go to 20
      nr = nion(kchem) - 1
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     evaluate r(jion, kchem) and jmax(kchem)                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 1 jion = 1, nr
      sr(jion, kchem) = sphi(jion, kchem) * exeta
    1 continue
c
c     now figure out ion stage with maximum occupation number
c
      jmax = 0
c
c     hypothesis (1): neutral atom
      if ( sr(   1, kchem) .gt. 1.0d0 ) go to 3
      do 2 jion = 2, nr
      if ( sr(jion, kchem) .ge. 1.0d0 ) go to 3
    2 continue
c     accept neutral atom as dominant ion state
      jmax = 1
      go to 7
c
c     hypothesis (2): highest ion
    3 if ( sr(nr  , kchem) .lt. 1.0d0 ) go to 5
      do 4 jion = 1, nr - 1
      if ( sr(jion, kchem) .le. 1.0d0 ) go to 5
    4 continue
c     accept highest ion as dominant ion state
      jmax = nr + 1
      go to 7
c
c     hypothesis (3): some intermediate ion state
    5 do 6 jion = 1, nr - 1
      if ( sr(jion    , kchem) .ge. 1.0d0 .and.
     .     sr(jion + 1, kchem) .le. 1.0d0 ) jmax = jion + 1
    6 continue
c
c     if the process failed, arbitrarily choose the middle ion
      if ( jmax .eq. 0 ) jmax = ( nr + 1 ) / 2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate p(jion, kchem), s(kchem), f(jion, kchem) and their     c
c     derivatives wrt to eta                                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    7 sp(jmax, kchem) = 1.0d0
      do 8 jion = jmax + 1, nion(kchem)
      sp(jion, kchem) = sp(jion - 1, kchem) * sr(jion - 1, kchem)
    8 continue
      do 9 jion = jmax - 1, 1, -1
      sp(jion, kchem) = sp(jion + 1, kchem) / sr(jion, kchem)
    9 continue
      do 10 jion = 1, nion(kchem)
      dspde(jion, kchem) = dfloat(  jmax - jion ) * sp(jion, kchem)
   10 continue
c
      ss   (kchem) = ssum( nion(kchem), sp   (1, kchem), 1 )
      dssde(kchem) = ssum( nion(kchem), dspde(1, kchem), 1 )
c
      do 11 jion = 1, nion(kchem)
      sf   (jion, kchem) = sp(jion, kchem) / ss(kchem)
      dsfde(jion, kchem) =
     .              ( dspde(jion, kchem) * ss   (kchem)
     .              - sp   (jion, kchem) * dssde(kchem) ) / ss(kchem)**2
   11 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close loop over chemical elements                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   20 continue
c
      return
      end
