      subroutine setphi
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     evaluate saha phi's                                              c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
c
      ckt = ck * t
      t32 = t ** 1.5d0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     phi for the "representative metal"                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      phim = 2.0d0 * gs(ise) * t32 * exp( - chim / ckt )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over all ions of all elements                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do kchem = 1, nchem
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hydrogen                                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
         if ( nucz(kchem) .eq. 1 ) then
c     
c     h2  molecule dissociation
c     
            sphi(1, kchem) =           ( gs(ish)**2 / gs(ish2) ) * t32
     .           * ( z (ish)**2 / z (ish2) )
     .           * exp( - ( 2.0d0*e0(ish)  - e0(ish2) ) / ckt )
c     
c     h2+ molecule dissociation
c     
            sphi(2, kchem) =  ( gs(ispr) * gs(ish) / gs(ish2p) ) * t32
     .           * ( z (ispr) * z (ish) / z (ish2p) )
     .           * exp( - ( e0(ispr) + e0(ish) - e0(ish2p) ) / ckt )
c     
c     h- negative ion dissociation
c     
            sphi(3, kchem) =    ( gs(ise) * gs(ish) / gs(ishm) ) * t32
     .           * ( 2.0d0  * z (ish) / z (ishm) )
     .           * exp( - ( e0(ish) - e0(ishm) ) / ckt )
c     
c     h atom ionization
c     
            sphi(4, kchem) =    ( gs(ise) * gs(ispr) / gs(ish) ) * t32
     .           * ( 2.0d0  * z (ispr) / z (ish) )
     .           * exp( - ( e0(ispr) - e0(ish) ) / ckt )
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     all other elements                                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
         else
            is1 = ichm1(kchem)
            do jion = 1, nion(kchem) - 1
               is = jion + is1 - 1
               sphi(jion, kchem) = ( gs(ise) * gs(is+1) / gs(is) ) * t32
     .              * ( 2.0d0  * z (is+1) / z (is) )
     .              * exp( - ( e0(is+1) - e0(is) ) / ckt )
            enddo
         end if
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close loop over elements                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      enddo
c     
      return
      end
