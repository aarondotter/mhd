      subroutine lambda
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     generate transformation matrix between reaction parameters and   c
c     occupation numbers                                               c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize matrix to zero
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      b = 0.0d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over all chemical elements
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nlam = 0
      do 80 kchem = 1, nchem
      iw1 = iwin1(kchem)
      iw2 = iwin2(kchem)
c
      kk  = nucz(kchem)
      if ( kk .eq. 1 ) then
c
c-----------------------------------------------------------------------
c     hydrogen                                                         c
c-----------------------------------------------------------------------
c
c           no contribution if only one molecule/atom/ion stage
c           in window
      if ( ssum(5, wwt(ish2), 1) .lt. 1.1d0 ) go to 77
c
c           h2 molecules
            if ( wwt(ish2 ) .lt. 0.1d0 ) go to 73
            nlam = nlam + 1
            b(nlam, ish2 ) =  1.0d0
            b(nlam, ish  ) = -2.0d0
c
c           h2+ molecules
   73       if ( wwt(ish2p) .lt. 0.1d0 ) go to 74
            nlam = nlam + 1
            b(nlam, ish2p) =  1.0d0
            b(nlam, ish  ) = -1.0d0
            b(nlam, ispr ) = -1.0d0
c
c           h- ions
   74       if ( wwt(ishm ) .lt. 0.1d0 ) go to 75
            nlam = nlam + 1
            b(nlam, ishm ) =  1.0d0
            b(nlam, ish  ) = -1.0d0
            b(nlam, ise  ) = -1.0d0
c
c           protons
   75       if ( wwt(ispr ) .lt. 0.1d0 ) go to 77
            nlam = nlam + 1
            b(nlam, ispr ) =  1.0d0
            b(nlam, ish  ) = -1.0d0
            b(nlam, ise  ) =  1.0d0
c
      else
c
c-----------------------------------------------------------------------
c     all elements other than hydrogen                                 c
c-----------------------------------------------------------------------
c
c           no contribution if only one ion stage in window
      if ( iw1 .eq. iw2 ) go to 77
      do 76 is = iw1, iw2 - 1
      nlam = nlam + 1
            b(nlam, is   ) = -1.0d0
            b(nlam, is+1 ) =  1.0d0
            b(nlam, ise  ) =  1.0d0
   76       continue
      end if
   77 continue
c
   80 continue
c
      return
      end
