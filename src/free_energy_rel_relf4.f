c----------- f1 and f2 same as nonrelativistic. f3,f4,ftot different
      subroutine f1
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     translational free energy and derivatives                        c
c***********************************************************************
c
      include 'types'
      include 'parms'
      include 'coms'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivatives of f1 wrt to occupation numbers                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ckt = ck * t
      ii = - mspes
c
      do is = 1, nspes - 1
         ii = ii + mspes + 1
c     use weight wwt to supress species below threshold
         fscr(is) = - wwt(is) * ckt * ( 1.5d0*log(  t  ) - log( sn(is) )
     .              + log( vol ) + log( gs(is) ) )
         dfdn  (is) = dfdn  (is) + fscr(is)
         d2fdnt(is) = d2fdnt(is) + wwt(is) * ( dfdn(is)/t - 1.5d0 * ck )
         d2fdnv(is) = d2fdnv(is) - wwt(is) * ckt / vol
         d2fdn2(is,is) = d2fdn2(is,is) + wwt(is) * ckt / sn(is)
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy, pressure, internal energy                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sum  =   ckt * sdot( nspes - 1, sn, 1, wwt , 1 )
      f(1) = - sum + sdot( nspes - 1, sn, 1, fscr, 1 )
      p(1) = sum / vol
      e(1) = 1.5 d0 * sum
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     second derivatives of f1 wrt to t and v                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      d2fdt2 = d2fdt2 - e(1) / t**2
      d2fdtv = d2fdtv - p(1) / t
      d2fdv2 = d2fdv2 + p(1) / vol
c
      return
      end
c
      subroutine f2
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     internal free energy                                             c
c***********************************************************************
c
c  ---------------------------------------------------------------------
c  ---- In this version, the treatment of H- and H2+ molecules is 
c  ---- exceptional. They are acting in s/r f2 as if they were neutrals.
c  ---- This is a formal decision (something has to be done) and it 
c  ---- is not justified from physics. This version is equivalent to
c  ---- the s/r f2 before the speed-up (see legacy version of f2 for
c  ---- more info).                [ 7/10/00 - WD ]
c  ---------------------------------------------------------------------
c  includes some streamlining by Johann Reiter (right-hand-side
c  comments: ). incorporated 26/11/95.
c  =====================================================================
c  modified 16/9/91: speed-up by Dimitri and Werner (Boulder June 91)
c  (WD - Aarhus)     is incorporated (exploiting simplicity of derivatives
c                    w.r.t. charged particles)
c  =====================================================================
c  modified 10/3/88: only zero part of dlnwdn actually used.
c                    loop 108 replaced by saxpy
c                    modify expression in loop 22
c
      include 'types'
      include 'parms'
      include 'coms'
c
      dimension dlnwdnz(mlev), z32(mspes)
c
      data  alph, bet
     .    / 10.0d0, 2.0d0 /
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      write(6,*) ' f2 entered'
      ckt    = ck  * t
      ckt2   = ckt * t                                                  
      vcon   = 4.0d0 * cpi / ( 3.0d0 * vol )
c
      do 101 is = 1, nspes
      z32(is)   = max( zs(is), 0.0d0 )**1.5d0
  101 continue
      z32n      = sdot( nspes - 1, sn, 1, z32, 1 )
c
      f(2) = 0.0d0
      e(2) = 0.0d0
      p(2) = 0.0d0
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     loop over all chemical elements                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 40 kchem = 1, nchem
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contributions from bare ions                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      is       = ispes ( nion(kchem), kchem )
      f(2)     = f(2)     + wwt(is) * e0(is) * sn(is)
      e(2)     = e(2)     + wwt(is) * e0(is) * sn(is)
      dfdn(is) = dfdn(is) + wwt(is) * e0(is)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     contributions from  molecule/atom/ion species with bound states  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 30 jion = 1, nion(kchem) - 1
c
c     species number, partition function data table number, and number
c     of bound states
c
      is = ispes(jion, kchem)
      ip = ipf  (jion, kchem)
c
      if ( ip .gt. 0 ) then                                             
         nl = nint(nlev(ip))                                            
      else                                                              
         nl = 1                                                         
      endif                                                             
c
c     zero dz/dt, dz/dv, dz/dn(i), d2z/dtdn(i), d2z/dvdn(i),d2z/dn(i)dn(
c
      dzdt   = 0.0d0
      dzdv   = 0.0d0
      d2zdt2 = 0.0d0
      d2zdtv = 0.0d0
      d2zdv2 = 0.0d0
      do 2 ips    = 1, nspes
      dzdn  (ips) = 0.0d0
      d2zdnt(ips) = 0.0d0
      d2zdnv(ips) = 0.0d0
      dlnpdn(ips) = 0.0d0
      do 1 jps = 1, nspes
      d2zdn2(jps, ips) = 0.0d0                                          
    1 continue
    2 continue
      xi     = 0.0d0
      psi    = 1.0d0
c........... some small psiln is needed, because there
c........... are (not very clever) division of type psiln/psiln
      psiln  = 1.d-20
      dlnpdv = 0.0d0
c
c     zero w, ln w, dlnw/dn(i), scratch
c
      do 102 il   = 1, mlev
      w   (il)    = 0.0d0
      wln (il)    = 0.0d0
      wscr(il)    = 0.0d0
      dlnwdnz(il) = 0.0d0
  102 continue
      do 4 ips = 1, nspes
      do 3 il  = 1, nl
      dlnwdn(il, ips) = 0.0d0
    3 continue
    4 continue
c
      if ( nl .le. 1 ) go to 14
c
c----------------------------------------------------------------------c
c     w, dw/dn(i), dz/dt, dz/dv, dz/dn(i), d2z/dtdn(i), d2z/dvdn(i),   c
c     d2z/dn(i)dn(j), d2z/dt2, d2z/dv2 for species with multiple bound c
c     states                                                           c
c----------------------------------------------------------------------c
c
c     d(ln w)/dn(i), d(ln psi)/dn(i), xi, psi
c
c     charged perturbers
      do 5 il = 1, nl
      dlnwdnz(il) = - vcon * chi(il, ip)
      wln    (il) = z32n * dlnwdnz(il)
    5 continue
c
c     neutral perturbers
      do 7 ips = 1, nspes - 1
      if( zs(ips) .gt. 0.0d0 ) go to 7
c     ignore perturbers below threshold
      if ( wwt(ips) .lt. 0.1d0 ) go to 7
c     ignore perturber if perturbed species is positively
c     charged, unless it happens to be h2+
      if ( zs(is) .le. 0.0d0 .or. is .eq. ish2p ) then
            do 6 il = 1, nl
            dlnwdn(il, ips) = - vcon * ( rad(il, ip) + rad0(ips) )**3
            wln(il) = wln(il) + sn(ips) * dlnwdn(il, ips)
    6       continue
            dlnpdn(ips) = vcon * ( rad0(is) + rad0(ips) )**3
            xi = xi + sn(ips) * dlnpdn(ips)
      end if
    7 continue
      if ( xi .gt. 0.0d0 ) then
            psiln = - alph * xi**bet
            psi   = max( exp( psiln ), 1.d-70 )                      
            call sscal(nspes-1,bet*psiln/xi,dlnpdn,1)                   
            dlnpdv = - bet * psiln / vol
      end if
c
c     w, z
c
      do 9 il = 1, nl
      w(il) = exp( wln(il) ) * ekt(il, ip)
    9 continue
      z(is) = max( ssum( nl, w, 1 ), 1.0d-70 )
c
c---->supress all contributions to f2 and its derivatives from species
c---->below threshold
      if ( wwt(is) .lt. 0.1d0 ) then
      z(is) = max( psi * z(is), 1.0d-70 )
            go to 30
      end if
c
c     dz/dv, dz/dt, d2z/dt2, d2z/dtv, d2z/dv2
c
      dzdv = - sdot( nl, w, 1, wln        , 1) / vol
      dzdt =   sdot( nl, w, 1, elev(1, ip), 1) / ckt2
c
      z1 = 2.0d0 * ckt                                                  
      do 10 il = 1, nl
      wscr(il) = elev(il, ip) * ( elev(il, ip) - z1 )                   
   10 continue
      d2zdt2 = sdot( nl, w, 1, wscr, 1) / ckt2**2
c
      do 110 il = 1, nl
      wscr(il) = elev(il, ip) * wln(il)
  110 continue
      d2zdtv = - sdot( nl, w, 1, wscr, 1) / ( ckt2 * vol )
c
      do 111 il = 1, nl
      wscr(il) = wln(il) * ( 2.0d0 + wln(il) )
  111 continue
      d2zdv2 = sdot( nl, w, 1, wscr, 1) / vol**2
c
c     dz/dn(i), d2z/dtdn(i), d2z/dvdn(i), d2z/dn(i)dn(j)
c
c     derivatives wrt total effective charge density
      do 112 il = 1, nl
      wscr(il) = w(il) * dlnwdnz(il)
  112 continue
      dzdnz   =     ssum(nl, wscr, 1)
      d2zdnzt =     sdot(nl, wscr, 1, elev(1, ip),  1) / ckt2
      d2zdnzv = - ( sdot(nl, wscr, 1, wln, 1) + dzdnz ) / vol
      d2zdnz2 =     sdot(nl, wscr, 1, dlnwdnz, 1)
c
c ---------- the evaluation of charged-charged derivatives
c ---------- happens symmetrically. only half of them computed.
c
      do 120 ips = 1, nspes - 1
      if ( zs (ips) .le. 0.0d0 ) go to 120
c     ignore perturbers below threshold
      if ( wwt(ips) .lt. 0.1d0 ) go to 120
c
c     first derivatives wrt charged perturbers
      dzdn  (ips) = z32(ips) * dzdnz
      d2zdnt(ips) = z32(ips) * d2zdnzt
      d2zdnv(ips) = z32(ips) * d2zdnzv
c     second derivatives wrt charged perturbers
      do 113 jps = ips, nspes - 1
      if(zs(jps) .eq. 0.0d0) goto 113
c     ignore perturbers below threshold
      if(wwt(jps) .lt. 0.1d0) goto 113
      d2zdn2(ips, jps) = z32(ips) * z32(jps) * d2zdnz2
      d2zdn2(jps, ips) = d2zdn2(ips, jps)
  113 continue
c
 120  continue
c
c ---------- the evaluation of neutral-charged derivatives
c ---------- is not done symmetrically. (neutral-neutral is,
c ---------- but it is not worth a third do-loop).
c
c     ignore if perturbed species is positively
c     charged, unless it happens to be h2+
c
      if ( zs(is) .le. 0.0d0 .or. is .eq. ish2p ) then
c
      do 125 ips = 1, nspes - 1
      if ( zs(ips ) .gt. 0.0d0 ) go to 125
c     ignore perturbers below threshold
      if ( wwt(ips) .lt. 0.1d0 ) go to 125
c
            do 11 il  = 1, nl
            wscr(il) = w(il) * dlnwdn(il, ips)
   11       continue
c           first derivatives wrt neutral perturbers
            dzdn  (ips) =     ssum(nl, wscr, 1)
            d2zdnt(ips) =     sdot(nl, wscr, 1 , elev(1, ip), 1) / ckt2
            d2zdnv(ips) = - ( sdot(nl, wscr, 1 , wln, 1)
     .                                              + dzdn(ips) ) / vol
            d2zdnznu    =     sdot(nl, wscr, 1, dlnwdnz, 1)
            do 114 jps = 1, nspes - 1
c           ignore perturbers below threshold
            if ( wwt(jps) .lt. 0.1d0 ) go to 114
c           second derivatives wrt neutral perturbers
            if( zs(jps) .eq. 0.0d0 ) then
                  d2zdn2(ips, jps) =
     .                              sdot(nl, wscr, 1, dlnwdn(1, jps), 1)
c           second mixed derivatives wrt neutral + charged perturbers
            else if( zs(jps) .gt. 0.0d0 ) then
                  d2zdn2(ips, jps) =  z32(jps) * d2zdnznu
            end if
            d2zdn2(jps, ips) = d2zdn2(ips, jps)
  114       continue
125   continue
c
      end if
c
      go to 20
c
c----------------------------------------------------------------------c
c     w, dw/dn(i), dz/dt, dz/dv, dz/dn(i), d2z/dtdn(i), d2z/dvdn(i),   c
c     d2z/dn(i)dn(j), d2z/dt2, d2z/dv2 for species with only one bound c
c     state                                                            c
c----------------------------------------------------------------------c
c
c     d(ln w)/dn(i), d(ln psi)/dn(i), xi, psi
c
   14 do 15 ips = 1, nspes - 1
c     ignore perturbers below threshold
      if ( wwt(ips) .lt. 0.1d0 ) go to 15
c
c     charged perturber
      if ( zs(ips) .gt. 0.d0 ) then
      dlnwdn(1, ips) = - vcon * chi0(is) * z32(ips)
c
c     neutral perturber. ignore if perturbed species is positively
c     charged, unless it happens to be h2+
      else if ( zs(is) .le. 0.0d0 .or. is .eq. ish2p ) then
            dlnwdn(1, ips) = - vcon * ( rad0(is) + rad0(ips) )**3
      dlnpdn(ips) = - dlnwdn(1, ips)
      xi = xi + sn(ips) * dlnpdn(ips)
      end if
   15 continue
      if ( xi .gt. 0.0d0 ) then
        psiln = - alph * xi**bet
        psi   = max( exp( psiln ), 1.d-70 )
        call sscal(nspes-1,bet*psiln/xi,dlnpdn,1)                       
        dlnpdv = - bet * psiln / vol
      end if
c
c     ln w
c
      wln(1) = sdot( nspes - 1, sn, 1, dlnwdn, mlev )
c
c     w, z
c
      w(1)  = g0(is) * exp( wln(1) )
      z(is) = max( w(1), 1.0d-70 )
c
c---->supress all contributions to f2 and its derivatives from species
c---->below threshold
c
      if ( wwt(is) .lt. 0.1d0 ) then
      z(is) = max( psi * z(is), 1.0d-70 )
            go to 30
      end if
c
c     dz/dv, dz/dt, d2z/dt2, d2z/dtv, d2z/dv2
c
      dzdv   =  - w(1) * wln(1) / vol
      d2zdv2 = wln(1) * ( 2.0d0 + wln(1) ) * w(1) / vol**2
c
c     dz/dn(i), d2z/dtdn(i), d2z/dvdn(i), d2z/dn(i)dn(j)
c
      z1 = (1.0d0 + wln(1)) * w(1) / vol                                
      do 17 jps = 1, nspes -1                                           
c     ignore perturbers below threshold
      if ( wwt(jps) .lt. 0.1d0 ) go to 17                               
c
      dzdn  (jps) =  w(1) * dlnwdn(1, jps)                              
      d2zdnv(jps) = - z1 * dlnwdn(1, jps)                               
c
      z2 = wwt(jps) * dzdn(jps)                                         
      do 16 ips = 1, nspes - 1                                          
c     use weight wwt to supress perturbers below threshold
      d2zdn2(ips, jps) = z2 * dlnwdn(1, ips)                            
   16 continue
   17 continue
c
c----------------------------------------------------------------------c
c     f2, p2, e2, df/dn(i), d2f/dtdn(i), d2f/dvdn(i), d2f/dt2,         c
c     d2f/dtv, d2f/dv2, d2f2/dn(i)dn(j)                                c
c----------------------------------------------------------------------c
c
c     contribution to free energy
   20 f(2) = f(2) + sn(is) * ( e0(is) - ckt * log( psi*z(is) ) )
c
c     contribution to pressure
      p(2) = p(2) + ckt * sn(is) * ( dzdv / z(is) + dlnpdv )
c
c     contribution to internal energy
      e(2) = e(2) + sn(is) * ( e0(is) + ckt2 * dzdt/z(is) )
c
c     df/dn(i), d2f/dtdn(i), d2f/dvdn(i)
c
      dfdn  (is) = dfdn  (is) + e0(is) - ckt * ( log(z(is)) + psiln )
      d2fdnt(is) = d2fdnt(is) - ckt * (dzdt/z(is) + log(z(is))/t)
     .                        - ck  * psiln
      d2fdnv(is) = d2fdnv(is) - ckt * ( dzdv/z(is) + dlnpdv )
c
      cktsn = ckt * sn(is)                                              
      betvolr = bet / vol                                               
      zisr = 1.0d0 / z(is)                                              
      trez = 1.0d0 / t                                                  
      z1 = dzdt * zisr                                                  
      z2 = dzdv * zisr                                                  
c
      do 21 ips = 1, nspes - 1
      dfdn  (ips) = dfdn  (ips) - cktsn *                               
     .                            ( dzdn(ips) * zisr + dlnpdn(ips) )    
      d2fdnt(ips) = d2fdnt(ips) - cktsn * (                             
     .             (d2zdnt(ips) + (trez - z1)*dzdn(ips)) * zisr         
     .             + dlnpdn(ips) * trez )                               
      d2fdnv(ips) = d2fdnv(ips) - cktsn * (                             
     .                       (d2zdnv(ips) - z2*dzdn(ips)) * zisr        
     .                       - betvolr * dlnpdn(ips) )                  
   21 continue
c
c     d2f/dn(i)dn(j)
c
c---  cktsn = ckt * sn(is)                                              
c---  bet1  = bet - 1.0d0                                               
c---  betpsi = bet*psiln                                                
      z2 = (bet - 1.0d0) / (bet*psiln)                                  
      do 23 jps = 1, nspes - 1                                          
      dzdnzi = dzdn(jps) * zisr                                         
      z1 = z2 * dlnpdn(jps)                                             
      do 22 ips = 1, nspes - 1                                          
      d2fdn2(ips, jps) = d2fdn2(ips, jps) - cktsn *
     .        ( ( d2zdn2(ips, jps) - dzdnzi*dzdn(ips) ) * zisr          
     .          + z1 * dlnpdn(ips))                                     
   22 continue
   23 continue
c
      do 24 jps = 1, nspes - 1
      d2fdn2(is, jps) = d2fdn2(is, jps)
     .                - ckt * ( dzdn(jps) * zisr + dlnpdn(jps) )        
   24 continue
      do 25 ips = 1, nspes - 1
      d2fdn2(ips, is) = d2fdn2(ips, is)
     .                - ckt * ( dzdn(ips) * zisr + dlnpdn(ips) )        
   25 continue
c
c     d2fdt2, d2fdtv, d2fdv2
c
      z1 = dzdt * zisr                                                  
      d2fdt2 = d2fdt2 - cktsn *                                         
     .                     ( d2zdt2 + (2.0d0/t - z1)*dzdt ) * zisr      
      d2fdtv = d2fdtv - cktsn * (                                       
     .                     ( d2zdtv + (1.0d0/t - z1)*dzdv ) * zisr      
     .                     + dlnpdv / t )
      d2fdv2 = d2fdv2 - cktsn * (
     .                     ( d2zdv2 - dzdv**2*zisr ) *zisr              
     .                     - (bet + 1.0d0) * dlnpdv / vol )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close loop over ion species                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   30 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close loop over elements                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   40 continue
c
      return
      end
c
      subroutine f3
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     free energy of degenerate electron gas and derivatives           c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
c..      write(6,*) ' f3 entered'
c
      ckt   = ck * t
c
c ---- switch nonrel/rel------------------------------------------------
cccc      beta=0.d0
       beta = ckt/cme/cc**2
c
      f12   = r_fd(0,0.5d0,eta,beta)
      f32   = r_fd(0,1.5d0,eta,beta)
      f52   = r_fd(0,2.5d0,eta,beta)
      fe12  = r_fd(1,0.5d0,eta,beta)
      fe32  = r_fd(1,1.5d0,eta,beta)
      fe52  = r_fd(1,2.5d0,eta,beta)
      fe72  = r_fd(1,3.5d0,eta,beta)
      f1232 = f12+beta*f32
      fe1232= fe12+beta*fe32
c
      frat1 = f1232/fe1232
      frat2 = (f32+0.5d0*beta*f52)/f1232
      frat3 = (f32+beta*f52)/f1232
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy, pressure, internal energy                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      f(3) = - ckt * sn(ise) * ( 0.666666666666667d0 * frat2 - eta )
      p(3) =   ckt * sn(ise) *   0.666666666666667d0 * frat2 / vol
      e(3) =   ckt * sn(ise) * frat3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivatives of f3 wrt electron number, t, and v                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dfdn  (ise) = dfdn  (ise) + ckt * eta
      d2fdnv(ise) = d2fdnv(ise) + ckt * detdv
      d2fdnt(ise) = d2fdnt(ise) + ck  * ( eta + t * detdt )
      d2fdv2 = d2fdv2 - ckt * sn(ise) * detdv / vol
      d2fdtv = d2fdtv+ck*sn(ise)*( 3.0d0*frat1-5.0d0*frat3/3.d0+
     .         (beta*f52/3.d0-fe32*frat1)/f1232)/vol
      d2fdt2 = d2fdt2+ck*sn(ise)*(2.25d0*frat1-2.5d0*frat3-0.5d0*beta/
     .         f1232*(fe72-f52-frat1*fe52*(3.d0+0.5d0*beta*fe52/f1232)
     .         ))/t
      d2fdn2(ise, ise) = d2fdn2(ise, ise) + ckt * frat1 / sn(ise)
c
      return
      end
c
      subroutine f4
c---------------- both for lack of a good theory for f4 (relativity
c---------------- is only one issue among exchange/diffraction),
c---------------- and for the fact that one would require substantially
c---------------- more elaborate derivatives with respect to T in f4
c---------------- due to the functional form of theta(T,eta), we adopt a
c---------------- nonrelativistic f4 (and the f4 part of ftot).
c---------------- however, eta must be computed non-relativistically here.
c
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     free energy of coulomb interactions and derivatives              c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      dimension   dxdn(mspes), d2xdnt(mspes), d2xdnv(mspes),
     .          d2xdn2(mspes, mspes)
      DIMENSION AB(20),WT(20)
c
      equivalence (dxdn  , dzdn  ), (d2xdnt, d2zdnt), (d2xdnv, d2zdnv),
     .            (d2xdn2, d2zdn2)
C
C     20 points Gauss Lagendre abscissas and weight
C
      DATA AB(  1) /   7.6526521133 4973337546 4040939883 82110 D -2 /
      DATA AB(  2) /   2.2778585114 1645078080 4961953685 74624 D -1 /
      DATA AB(  3) /   3.7370608871 5419560672 5481770249 27237 D -1 /
      DATA AB(  4) /   5.1086700195 0827098004 3640509552 50998 D -1 /
      DATA AB(  5) /   6.3605368072 6515025452 8366962262 85936 D -1 /
      DATA AB(  6) /   7.4633190646 0150792614 3050703556 41590 D -1 /
      DATA AB(  7) /   8.3911697182 2218823394 5290617015 20685 D -1 /
      DATA AB(  8) /   9.1223442825 1325905867 7524412032 98113 D -1 /
      DATA AB(  9) /   9.6397192727 7913791267 6661311972 77221 D -1 /
      DATA AB( 10) /   9.9312859918 5094924786 1223884713 20278 D -1 /
C
      DATA WT(  1) /   1.5275338713 0725850698 0843319550 97593 D -1 /
      DATA WT(  2) /   1.4917298647 2603746787 8287370019 69436 D -1 /
      DATA WT(  3) /   1.4209610931 8382051329 2983250671 64933 D -1 /
      DATA WT(  4) /   1.3168863844 9176626898 4944997481 63134 D -1 /
      DATA WT(  5) /   1.1819453196 1518417312 3773777113 82287 D -1 /
      DATA WT(  6) /   1.0193011981 7240435036 7501354803 49876 D -1 /
      DATA WT(  7) /   8.3276741576 7047487247 5814322204 62061 D -2 /
      DATA WT(  8) /   6.2672048334 1090635695 0653518704 16063 D -2 /
      DATA WT(  9) /   4.0601429800 3869413310 3995227493 21098 D -2 /
      DATA WT( 10) /   1.7614007139 1521183118 6196235185 28163 D -2 /
C
      data   niter,   conv
     .    /    20, 1.d-10 /
c
ccc      write(6,*) ' f4 entered'
c
      DO K=1,10
        AB(K+10) = -AB(11-K)
        WT(K+10) =  WT(11-K)
      END DO
c
      ckt = ck * t
c
c     ============================================================
c     ============= itau = 1 for usual tau
c     ============= itau = 0 for switched-off tau (tau=1, that is)
c     ============================================================
      itau = 1
c     ============================================================
c
c ---- switch nonrel/rel (leave it beta=0 except when testing) -------
c>>>>>>>>>>>>>>>>... test version for legacy (with wrong T derivatives)
cccc  beta=0.d0
      beta  = ckt/cme/cc**2
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> c
c                      separate new eta part for f4
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> c
c
c --------- cannibalized from s/r new eta to create 'local' eta (=etaf4)
      snise = max(sn(ise),1.d-70)
      if(.true.) goto 5555
      rhs = ceta * snise / ( t**1.5d0 * vol )
c
c........... initialization to a few per cent accuracy
c........... (see ref. in dappen, astron. astrphys., 1980)
c
      if(rhs.lt.1.d-6) then
            etaf4 = log(rhs)
      else if(rhs.lt.100.d0) then
            g   =rhs
            g2  =rhs*rhs
            et0 =g+(g2/(4.45+2.257*g+g2))*exp((1.32934*g)**0.66666667)
            etaf4 = log(et0)
      else
            etaf4 = (1.32934*rhs)**0.66666667
      end if
      etaf4 = etaf4 + 0.120782
c
      f12   = r_fd(0,0.5d0,etaf4,beta)
      f32   = r_fd(0,1.5d0,etaf4,beta)
      fe12  = r_fd(1,0.5d0,etaf4,beta)
      fe32  = r_fd(1,1.5d0,etaf4,beta)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iter  = 0
   31 detaf4  = ( rhs - f12 - beta * f32 ) / ( fe12 + beta * fe32 )
      etaf4 = etaf4 + detaf4
      f12   = r_fd(0,0.5d0,etaf4,beta)
      f32   = r_fd(0,1.5d0,etaf4,beta)
      fe12  = r_fd(1,0.5d0,etaf4,beta)
      fe32  = r_fd(1,1.5d0,etaf4,beta)
      if( dabs(detaf4) .le. conv ) go to 33
      iter = iter + 1
      if( iter .le. niter ) go to 31
c
c     failure to converge
      write  ( iout, 22 ) tlog, rholog, snise, etaf4
   22 format ( ' nonconvergence of degeneracy parameter in f4' /
     .         ' log t =',f5.2,' log rho =',f6.2,' ne =',1pe10.3,
     .         ' eta   =',e10.3 )
      stop 'newef4'
c
   33 continue
ccc      write(6,*) ' etaf4,eta = ',etaf4,eta
c
5555  etaf4 = eta
      f12   = r_fd(0,0.5d0,etaf4,beta)
      f32   = r_fd(0,1.5d0,etaf4,beta)
      fe12  = r_fd(1,0.5d0,etaf4,beta)
      fe32  = r_fd(1,1.5d0,etaf4,beta)
      f52    = r_fd(0,2.5d0,etaf4,beta)
      fe52   = r_fd(1,2.5d0,etaf4,beta)
      f2e12  = r_fd(3,0.5d0,etaf4,beta)
      f2e32  = r_fd(3,1.5d0,etaf4,beta)
      f2e52  = r_fd(3,2.5d0,etaf4,beta)
      f2e72  = r_fd(3,3.5d0,etaf4,beta)
      f3e12  = r_fd(6,0.5d0,etaf4,beta)
      f3e32  = r_fd(6,1.5d0,etaf4,beta)
      f3e52  = r_fd(6,2.5d0,etaf4,beta)
c
      thf4   = ( fe12 + beta * fe32 ) / ( f12 + beta * f32 )
      dthf4  = thf4*((f2e12+beta*f2e32)/(fe12+beta*fe32)-thf4)
      d2thf4 = thf4*((f3e12+beta*f3e32)/(fe12+beta*fe32)-3.0d0*dthf4
     .         -thf4**2)
c ------------------ game,dgamde,d2gamde2 are local
c ------------------ and need not be given special names in s/r f4
      game    = ( fe32 + beta * fe52 ) / ( f32 + beta * f52 )
      dgamde  = game*((f2e32+beta*f2e52)/(fe32+beta*fe52)-game)
      d2gamde2= game*((f3e32+beta*f3e52)/(fe32+beta*fe52)-3.d0*dgamde
     .          -game**2)
c
      de4dv   = - 1.0d0 / ( thf4 * vol )
      de4dt   = - (fe32+beta*fe52)/(fe12+beta*fe32) / t
      de4dn   =   1.0d0 / ( thf4 * snise )
      d2e4dn2 = - de4dn**2 * (f2e12+beta*f2e32)/(fe12+beta*fe32)
      d2e4dnt = - de4dn * de4dt *(thf4+dthf4/thf4-game-dgamde/game)
      d2e4dnv = - dthf4 * de4dn * de4dv / thf4
      d2e4dtv = - de4dt * de4dv *(thf4+dthf4/thf4-game-dgamde/game)
      d2e4dt2 = - (1.0d0/t+de4dt*(thf4+dthf4/thf4-game-dgamde/game))
     .          * de4dt+beta*((fe52-f2e72)/t-f2e52*de4dt)/(f12+beta
     .          * f32)/2.d0/t/thf4
      d2e4dv2 =   ( thf4 - dthf4 / thf4 ) * de4dv**2
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> c
c                        end separate new eta part
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sums over charges                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sumn = max( sdot(nspes - 1, zmask, 1, sn, 1), 1.d-70 )
      zn   = max( sdot(nspes - 1, zs   , 1, sn, 1), 1.d-70 )
      znt  = max( sdot(nspes - 1, zsq  , 1, sn, 1) + sn(ise)*thf4,
     .                                              1.d-70 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     parameter x and its derivatives                                  c
c     use weight wwt to supress species below threshold                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      f1232  = f12+beta*f32
      f3252  = f32+beta*f52
      frf4   = f1232 / f3252
      x      = cx * frf4 * (zn/sumn) * dsqrt( znt ) / dsqrt( vol*t**3 )
c
c-----------------------------------------------------------------------
c     zero everything
c-----------------------------------------------------------------------
c
      do 2 is = 1, nspes
      fscr  (is) = 0.0d0
      dxdn  (is) = 0.0d0
      d2xdnt(is) = 0.0d0
      d2xdnv(is) = 0.0d0
      do 1 js = 1, nspes
      d2xdn2(is, js) = 0.0d0
    1 continue
    2 continue
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      do 3 is = 1, nspes - 1
      dxdn(is) = wwt(is) *
     .              x * ( zs(is)/zn - zmask(is)/sumn+0.5d0*zsq(is)/znt )
    3 continue
c
      dxdn(ise) = x * ( 0.5d0* (thf4 + sn(ise)*dthf4*de4dn)/znt
     .                - 1.5d0* frf4 * de4dn
     .                + 1.0d0 / sn(ise) )
c
      dxdt = x * ( de4dt*(thf4 - 1.5d0*frf4 + 0.5d0*sn(ise)*dthf4/znt)
     .           - 1.5d0/t   )
      dxdv = x * ( de4dv*(thf4 - 1.5d0*frf4 + 0.5d0*sn(ise)*dthf4/znt)
     .           - 0.5d0/vol )
c
c-----------------------------------------------------------------------
c     second derivatives
c-----------------------------------------------------------------------
c
      do 5 js = 1, nspes - 1
      if ( wwt(js) .lt. 0.1d0 ) go to 5
c
      if ( zs (js) .ne. 0.0d0 ) then
            do 4 is = 1, js
            d2xdn2(is, js) = wwt(is) *
     .                   ( dxdn(is) * dxdn(js) / x
     .                   + x * (       zmask(is) * zmask(js) / sumn**2
     .                         -       zs   (is) * zs   (js) / zn**2
     .                       - 0.5d0 * zsq  (is) * zsq  (js) / znt**2) )
    4       continue
c
            if( js .gt. 1 ) call scopy( js - 1, d2xdn2( 1,js), 1,
     .                                          d2xdn2(js, 1), mspes )
      d2xdnt(js) = dxdt * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthf4*de4dt/znt**2
      d2xdnv(js) = dxdv * dxdn(js) / x
     .               - 0.5d0 * x * zsq(js)*sn(ise) * dthf4*de4dv/znt**2
      end if
    5 continue
c
      do 6 is = 1, nspes - 1
      d2xdn2(is, ise) = wwt(is) *
     .       ( dxdn(is) * dxdn(ise) / x
     .   - 0.5d0 * x * zsq(is) * (thf4 + sn(ise)*dthf4*de4dn)/znt**2 )
    6 continue
      call scopy( nspes - 1, d2xdn2(1,ise), 1, d2xdn2(ise,1), mspes )
c
      d2xdn2(ise, ise) =
     .    dxdn(ise)**2/x
     .    - x * ( 1.0d0/sn(ise)**2
     .       + 1.5d0 * frf4 * (d2e4dn2 + (thf4 - 1.5d0*frf4)*de4dn**2)
     .         - 0.5d0 * ( d2thf4 * sn(ise) * de4dn**2
     .                  + dthf4 * (2.d0*de4dn + sn(ise)*d2e4dn2)
     .                - (thf4 + sn(ise)*dthf4*de4dn)**2/znt ) / znt  )
c
      d2xdnt(ise) = dxdt*dxdn(ise)/x + x * (
     .  -1.5d0*frf4* (d2e4dnt + (thf4 - 1.5d0*frf4)*de4dn*de4dt)
     .  +0.5d0*( d2thf4 * sn(ise)*de4dn*de4dt
     .       + dthf4 * (de4dt + sn(ise)*d2e4dnt)
     .     - (thf4+sn(ise)*dthf4*de4dn)*sn(ise)*dthf4*de4dt/znt)/znt )
c
      d2xdnv(ise) = dxdv*dxdn(ise)/x + x * (
     .  -1.5d0*frf4* (d2e4dnv + (thf4 - 1.5d0*frf4)*de4dn*de4dv)
     .  +0.5d0*( d2thf4 * sn(ise)*de4dn*de4dv
     .       + dthf4 * (de4dv + sn(ise)*d2e4dnv)
     .     - (thf4+sn(ise)*dthf4*de4dn)*sn(ise)*dthf4*de4dv/znt)/znt )
c
      d2xdt2 = dxdt**2/x + x * ( 1.5d0/t**2
     .     + d2e4dt2 * ( thf4 - 1.5d0*frf4 + 0.5d0*sn(ise)*dthf4/znt )
     .       + de4dt**2 * ( dthf4 - 1.5d0*frf4*(thf4 - 1.5d0 * frf4)
     .         + 0.5d0*sn(ise)*(d2thf4 - sn(ise)*dthf4**2/znt)/znt ) )
c
      d2xdtv = dxdt*dxdv/x + x * (
     .       d2e4dtv * ( thf4 - 1.5d0*frf4 + 0.5d0*sn(ise)*dthf4/znt )
     .    + de4dt * de4dv * ( dthf4 - 1.5d0*frf4*(thf4 - 1.5d0 * frf4)
     .         + 0.5d0*sn(ise)*(d2thf4 - sn(ise)*dthf4**2/znt)/znt ) )
c
      d2xdv2 = dxdv**2/x + x * ( 0.5d0/vol**2
     .     + d2e4dv2 * ( thf4 - 1.5d0*frf4 + 0.5d0*sn(ise)*dthf4/znt )
     .     + de4dv**2 * ( dthf4 - 1.5d0*frf4*(thf4 - 1.5d0 * frf4)
     .         + 0.5d0*sn(ise)*(d2thf4 - sn(ise)*dthf4**2/znt)/znt ) )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     tau and its derivatives                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ( itau . eq . 1 .or.
     .     (rholog.gt.-1.95 .and. tlog.lt.4.6 ) .or.
     .     (rholog.gt.-1.55 .and. tlog.lt.4.8 ) ) then
        if ( itau .eq. 0 ) then
            write (6,*) ' *** F4: tau enforced. log rho, log t = ',
     .                    rholog,tlog
        end if
c
        IF (X.LE.1.D0) THEN
          F1 = 0.D0
          F2 = 0.D0
          F3 = 0.D0
          DO J=1,20
             F1 = F1+(1.D0+AB(J))**2*WT(J)/(1.D0+X*(1.D0+AB(J))/2.D0)
             F2 = F2+(1.D0+AB(J))**3*WT(J)/(1.D0+X*(1.D0+AB(J))/2.D0)**2
             F3 = F3+(1.D0+AB(J))**4*WT(J)/(1.D0+X*(1.D0+AB(J))/2.D0)**3
          END DO
          TAU  = F1*3.D0/8.D0
          DTAU =-F2*3.D0/16.D0
          D2TAU= F3*3.D0/16.D0
        ELSE
          DP1 =   3.D0 * ( log(1.D0 + X)
     .          - X*(1.D0 - 0.5D0*X) ) / X**3
          DP2 = - 3.D0 * ( DP1 - 1.D0/(1.D0 + X)    ) / X
          DP3 =   - ( 4.D0*DP2 + 3.D0/(1.D0 + X)**2 ) / X
          TAU   = DP1
          DTAU  = DP2
          D2TAU = DP3
        END IF
      else
        tau   = 1.d0
        dtau  = 0.d0
        d2tau = 0.d0
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy, pressure, internal energy                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     free energy
      f(4) = - cf4 * tau * dsqrt( znt**3 / (t * vol) )
c
c     pressure
      p(4) = f(4) * ( 0.5d0/vol-1.5d0* sn(ise) * dthf4 * de4dv / znt
     .                          - dxdv * dtau / tau )
c
c     internal energy
      e(4) = f(4) * ( 1.5d0/t -1.5d0* sn(ise) * dthf4 * de4dt / znt
     .                          - dxdt * dtau / tau ) * t
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     derivatives of f4                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----------------------------------------------------------------------
c     first derivatives                                                c
c-----------------------------------------------------------------------
c
      df4dt =   f(4) * ( -0.5d0/ t + 1.5d0*sn(ise) * dthf4 * de4dt / znt
     .                            + dxdt * dtau / tau )
      df4dv = - p(4)
c
      do 7 is = 1, nspes - 1
      fscr(is) = f(4)*wwt(is)*(1.5d0*zsq(is)/znt + dxdn(is)*dtau/tau)
      dfdn(is) = dfdn(is) + fscr(is)
    7 continue
      fscr(ise) = f(4)*( 1.5d0* ( thf4 + sn(ise)*dthf4*de4dn )/znt
     .                   + dxdn(ise) * dtau/tau )
      dfdn(ise) = dfdn(ise) + fscr(ise)
c
c-----------------------------------------------------------------------
c     second derivatives                                               c
c-----------------------------------------------------------------------
c
      do 9 js = 1, nspes - 1
      if ( wwt(js) .lt. 0.1d0 ) go to 9
c
      if ( zs(js)  .ne. 0.d0 ) then
            do 8 is = 1, js
            d2fdn2(is, js) = d2fdn2(is, js) + wwt(is) *
     .             ( fscr(is) * fscr(js) / f(4)
     .             + ( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(js)
     .               + d2xdn2(is, js) * dtau / tau
     .               - 1.5d0* zsq(is) * zsq(js) / znt**2 ) * f(4) )
    8       continue
            if( js .gt. 1 ) call scopy( js - 1, d2fdn2(1 ,js), 1,
     .                                          d2fdn2(js, 1), mspes )
c
            d2fdnt(js) = d2fdnt(js)
     .                 + fscr(js) * df4dt / f(4)
     .                 + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdt
     .                   + d2xdnt(js) * dtau/tau
     .                   -1.5d0*zsq(js)*sn(ise)*dthf4*de4dt/znt**2)*f(4)
c
            d2fdnv(js) = d2fdnv(js)
     .                 + fscr(js) * df4dv / f(4)
     .                 + ( (d2tau/tau - (dtau/tau)**2) * dxdn(js) * dxdv
     .                   + d2xdnv(js) * dtau/tau
     .                   -1.5d0*zsq(js)*sn(ise)*dthf4*de4dv/znt**2)*f(4)
      end if
    9 continue
c
      do 10 is = 1, nspes - 1
      d2fdn2(is, ise) = d2fdn2(is, ise) + wwt(is) *
     .      ( fscr(is) * fscr(ise) / f(4)
     .      + f(4)*( (d2tau/tau - (dtau/tau)**2) * dxdn(is) * dxdn(ise)
     .             + d2xdn2(is, ise) * dtau / tau
     .             -1.5d0*zsq(is)*(thf4+sn(ise)*dthf4*de4dn)/znt**2))
   10 continue
      call scopy( nspes - 1, d2fdn2(1,ise), 1, d2fdn2(ise,1), mspes )
c
      d2fdn2(ise, ise) = d2fdn2(ise, ise)
     .   + fscr(ise)**2 / f(4)
     .   + f(4) * ( (d2tau/tau - (dtau/tau)**2) * dxdn(ise)**2
     .            + d2xdn2(ise, ise) * dtau / tau
     .            + 1.5d0*( d2thf4 * sn(ise) * de4dn**2
     .                    + dthf4 * ( 2.d0*de4dn + sn(ise) * d2e4dn2 )
     .                    - (thf4 + sn(ise)*dthf4*de4dn)**2/znt) / znt )
c
      d2fdnt(ise) = d2fdnt(ise) + fscr(ise)*df4dt/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdn(ise)
     . + d2xdnt(ise) * dtau/tau
     . + 1.5d0 *( d2thf4*sn(ise)*de4dt*de4dn
     .        +dthf4*(de4dt + sn(ise)*d2e4dnt)
     .        -(thf4+sn(ise)*dthf4*de4dn)*sn(ise)*dthf4*de4dt/znt)/znt )
c
      d2fdnv(ise) = d2fdnv(ise) + fscr(ise)*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdv*dxdn(ise)
     . + d2xdnv(ise) * dtau/tau
     . + 1.5d0 *( d2thf4*sn(ise)*de4dv*de4dn
     .        +dthf4*(de4dv + sn(ise)*d2e4dnv)
     .        -(thf4+sn(ise)*dthf4*de4dn)*sn(ise)*dthf4*de4dv/znt)/znt )
c
      d2fdtv = d2fdtv + df4dt*df4dv/f(4) + f(4) *
     . ( (d2tau/tau - (dtau/tau)**2)*dxdt*dxdv
     . + d2xdtv * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thf4*de4dt*de4dv
     .            + dthf4*(d2e4dtv-sn(ise)*dthf4*de4dt*de4dv/znt))/znt )
c
      d2fdt2 = d2fdt2 + df4dt**2/f(4) + f(4) *
     . ( 0.5d0 / t**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdt**2
     . + d2xdt2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thf4*de4dt**2
     .               + dthf4*(d2e4dt2-sn(ise)*dthf4*de4dt**2/znt))/znt )
c
      d2fdv2 = d2fdv2 + df4dv**2/f(4) + f(4) *
     . ( 0.5d0 / vol**2
     . + (d2tau/tau - (dtau/tau)**2)*dxdv**2
     . + d2xdv2 * dtau/tau
     . + 1.5d0 * sn(ise) * ( d2thf4*de4dv**2
     .               + dthf4*(d2e4dv2-sn(ise)*dthf4*de4dv**2/znt))/znt )
c
      return
      end
c
      subroutine ftot
c---------------- both for lack of a good theory for f4 (relativity
c---------------- is only one issue among exchange/diffraction),
c---------------- and for the fact that one would require substantially
c---------------- more elaborate derivatives with respect to T in f4
c---------------- due to the functional form of theta(T,eta), we adopt a
c---------------- nonrelativistic f4 (and the f4 part of ftot).
c---------------- however, eta must be computed non-relativistically here.
c
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     calculate free energy for scan procedure                         c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      DIMENSION AB(20),WT(20)
c
      data  alph, bet
     .    / 10.0d0, 2.0d0 /
c
      data   niter,   conv
     .    /    20, 1.d-10 /
C
C     20 points Gauss Lagendre abscissas and weight
C
      DATA AB(  1) /   7.6526521133 4973337546 4040939883 82110 D -2 /
      DATA AB(  2) /   2.2778585114 1645078080 4961953685 74624 D -1 /
      DATA AB(  3) /   3.7370608871 5419560672 5481770249 27237 D -1 /
      DATA AB(  4) /   5.1086700195 0827098004 3640509552 50998 D -1 /
      DATA AB(  5) /   6.3605368072 6515025452 8366962262 85936 D -1 /
      DATA AB(  6) /   7.4633190646 0150792614 3050703556 41590 D -1 /
      DATA AB(  7) /   8.3911697182 2218823394 5290617015 20685 D -1 /
      DATA AB(  8) /   9.1223442825 1325905867 7524412032 98113 D -1 /
      DATA AB(  9) /   9.6397192727 7913791267 6661311972 77221 D -1 /
      DATA AB( 10) /   9.9312859918 5094924786 1223884713 20278 D -1 /
C
      DATA WT(  1) /   1.5275338713 0725850698 0843319550 97593 D -1 /
      DATA WT(  2) /   1.4917298647 2603746787 8287370019 69436 D -1 /
      DATA WT(  3) /   1.4209610931 8382051329 2983250671 64933 D -1 /
      DATA WT(  4) /   1.3168863844 9176626898 4944997481 63134 D -1 /
      DATA WT(  5) /   1.1819453196 1518417312 3773777113 82287 D -1 /
      DATA WT(  6) /   1.0193011981 7240435036 7501354803 49876 D -1 /
      DATA WT(  7) /   8.3276741576 7047487247 5814322204 62061 D -2 /
      DATA WT(  8) /   6.2672048334 1090635695 0653518704 16063 D -2 /
      DATA WT(  9) /   4.0601429800 3869413310 3995227493 21098 D -2 /
      DATA WT( 10) /   1.7614007139 1521183118 6196235185 28163 D -2 /
C
      DO K=1,10
        AB(K+10) = -AB(11-K)
        WT(K+10) =  WT(11-K)
      END DO
c
c     ============================================================
c     ============= itau = 1 for usual tau
c     ============= itau = 0 for switched-off tau (tau=1, that is)
c     ============================================================
      itau = 1
c     ============================================================
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     translational free energy                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ckt = ck * t
      ii = - mspes
      do 1 is = 1, nspes - 1
      ii = ii + mspes + 1
      fscr(is) = - wwt(is) * ckt * ( 1.5d0*log(  t  ) - log( sn(is) )
     .                                 + log( vol ) + log( gs(is) ) )
    1 continue
      sum  =   ckt * sdot( nspes - 1, sn, 1, wwt , 1 )
c
      f(1) = - sum + sdot( nspes - 1, sn, 1, fscr, 1 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     internal free energy                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      vcon = 4.0d0 * cpi / ( 3.0d0 * vol )
      f(2) = 0.0d0
c
      do 40 kchem = 1, nchem
c
c     bare ions
c
      is = ispes( nion(kchem), kchem )
      f(2) = f(2) + wwt(is) * e0(is) * sn(is)
c
c     molecule/atom/ion species with bound states
c
      do 30 jion = 1, nion(kchem) - 1
      is = ispes(jion, kchem)
      if ( wwt(is) .lt. 0.1d0 ) go to 30
c
      xi     = 0.0d0
      psi    = 1.0d0
      psiln  = 0.0d0
c
      do 2 il = 1, mlev
      w  (il) = 0.0d0
      wln(il) = 0.0d0
    2 continue
      do 4 ips = 1, nspes
      do 3 il  = 1, mlev
      dlnwdn(il, ips) = 0.0d0
    3 continue
    4 continue
c
      ip = ipf(jion, kchem)
      nl = 1
      if ( ip .le. 0 ) go to 14
      nl = nint(nlev(ip))
c
c     species with multiple bound states
c
      do 7 ips = 1, nspes - 1
      if ( wwt(ips) .lt. 0.1d0 ) go to 7
      if ( zs(ips) .gt. 0.d0 ) then
      zcon = vcon * dsqrt( zs(ips)**3 )
      do 5 il = 1, nl
      dlnwdn(il, ips) = - zcon * chi(il, ip)
    5       continue
      else  if ( zs(is) .le. 0.0d0 .or. is .eq. ish2p ) then
      do 6 il = 1, nl
      dlnwdn(il, ips) = - vcon * ( rad(il, ip) + rad0(ips) )**3
    6       continue
      xi = xi + vcon * sn(ips) * ( rad0(is) + rad0(ips) )**3
      end if
    7 continue
c
      if ( xi .gt. 0.0d0 ) then
      psiln = - alph * xi**bet
      psi   = max( exp( psiln ), 1.d-70 )
      end if
c
      do  9 ips = 1, nspes - 1
      if ( wwt(ips) .lt. 0.1d0 ) go to 9
      do 8 il  = 1, nl
      wln(il) = wln(il) + sn(ips) * dlnwdn(il,ips)
    8 continue
    9 continue
      do 10 il = 1, nl
      w(il) = exp( wln(il) ) * ekt(il, ip)
   10 continue
      z(is) = max( psi * ssum( nl, w, 1 ), 1.0d-70 )
c
      go to 20
c
c     species with only one bound state
c
   14 do 15 ips = 1, nspes - 1
      if ( wwt(ips) .lt. 0.1d0 ) go to 15
      if ( zs(ips) .gt. 0.d0 ) then
      dlnwdn(1, ips) = - vcon * chi0(is) * dsqrt( zs(ips)**3 )
      else if ( zs(is) .le. 0.0d0 .or. is .eq. ish2p ) then
            dlnwdn(1, ips) =  - vcon * ( rad0(is) + rad0(ips) )**3
      xi = xi + vcon * sn(ips) * ( rad0(is) + rad0(ips) )**3
      end if
   15 continue
c
      if ( xi .gt. 0.0d0 ) then
      psiln = - alph * xi**bet
      psi   = max( exp( psiln ), 1.d-70 )
      end if
c
      wln(1) = sdot( nspes - 1, sn, 1, dlnwdn, mlev )
      w  (1) = g0(is) * exp( wln(1) )
      z(is) = max( psi * w(1), 1.0d-70 )
c
   20 f(2) = f(2) + sn(is) * ( e0(is) - ckt * log(z(is)) )
c
   30 continue
   40 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy of degenerate electron gas                           c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call new eta
c
c ---- switch nonrel/rel------------------------------------------------
cccc     beta=0.d0
      beta  = ckt/cme/cc**2
c
      f12   = r_fd(0,0.5d0,eta,beta)
      f32   = r_fd(0,1.5d0,eta,beta)
      f52   = r_fd(0,2.5d0,eta,beta)
      frat  = (f32+0.5d0*beta*f52)/(f12+beta*f32)
c
      f(3) = - ckt * sn(ise) * ( 0.666666666666667 d0 * frat - eta )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     free energy of coulomb interactions                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c ---- switch nonrel/rel (leave it beta=0 except when testing) -------
c>>>>>>>>>>>>>>>>... test version for legacy (with wrong T derivatives)
cccc  beta=0.d0
      beta  = ckt/cme/cc**2
c
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> c
c                      separate new eta part for ftot
c >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> c
c
c --------- cannibalized from s/r new eta to create 'local' eta (=etaf4)
      snise = max(sn(ise),1.d-70)
      if(.true.) goto 5555
      rhs = ceta * snise / ( t**1.5d0 * vol )
c
c........... initialization to a few per cent accuracy
c........... (see ref. in dappen, astron. astrphys., 1980)
c
      if(rhs.lt.1.d-6) then
            etaf4 = log(rhs)
      else if(rhs.lt.100.d0) then
            g   =rhs
            g2  =rhs*rhs
            et0 =g+(g2/(4.45+2.257*g+g2))*exp((1.32934*g)**0.66666667)
            etaf4 = log(et0)
      else
            etaf4 = (1.32934*rhs)**0.66666667
      end if
      etaf4 = etaf4 + 0.120782
c
      f12   = r_fd(0,0.5d0,etaf4,beta)
      f32   = r_fd(0,1.5d0,etaf4,beta)
      fe12  = r_fd(1,0.5d0,etaf4,beta)
      fe32  = r_fd(1,1.5d0,etaf4,beta)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     newton-raphson iteration
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      iter  = 0
   31 detaf4  = ( rhs - f12 - beta * f32 ) / ( fe12 + beta * fe32 )
      etaf4 = etaf4 + detaf4
      f12   = r_fd(0,0.5d0,etaf4,beta)
      f32   = r_fd(0,1.5d0,etaf4,beta)
      fe12  = r_fd(1,0.5d0,etaf4,beta)
      fe32  = r_fd(1,1.5d0,etaf4,beta)
      if( dabs(detaf4) .le. conv ) go to 33
      iter = iter + 1
      if( iter .le. niter ) go to 31
c
c     failure to converge
      write  ( iout, 22 ) tlog, rholog, snise, etaf4
   22 format ( ' nonconvergence of degeneracy parameter in ftot' /
     .         ' log t =',f5.2,' log rho =',f6.2,' ne =',1pe10.3,
     .         ' eta   =',e10.3 )
      stop 'neweft'
c
   33 continue
ccc      write(6,*) ' etaftot,eta = ',etaf4,eta
c
5555  etaf4 = eta
      f12   = r_fd(0,0.5d0,etaf4,beta)
      f32   = r_fd(0,1.5d0,etaf4,beta)
      fe12  = r_fd(1,0.5d0,etaf4,beta)
      fe32  = r_fd(1,1.5d0,etaf4,beta)
      f52    = r_fd(0,2.5d0,etaf4,beta)
c
      thf4   = ( fe12 + beta * fe32 ) / ( f12 + beta * f32 )
      frf4   = (f32+0.5d0*beta*f52)/(f12+beta*f32)
c
      sumn = max( sdot(nspes - 1, zmask, 1, sn, 1), 1.d-70 )
      zn   = max( sdot(nspes - 1, zs   , 1, sn, 1), 1.d-70 )
      znt  = max( sdot(nspes - 1, zsq  , 1, sn, 1) +sn(ise)*thf4,
     .                                              1.d-70 )
c
      x = cx * frf4 * (zn/sumn) * dsqrt( znt ) / dsqrt( vol*t**3 )
c
      if ( itau . eq . 1 .or.
     .     (rholog.gt.-1.95 .and. tlog.lt.4.6 ) .or.
     .     (rholog.gt.-1.55 .and. tlog.lt.4.8 ) ) then
c
        IF (X.LE.1.D0) THEN
          F1 = 0.D0
          DO J=1,20
             F1 = F1+(1.D0+AB(J))**2*WT(J)/(1.D0+X*(1.D0+AB(J))/2.D0)
          END DO
          TAU  = F1*3.D0/8.D0
        ELSE
          DP1 =   3.D0 * ( log(1.D0 + X)
     .          - X*(1.D0 - 0.5D0*X) ) / X**3
          TAU   = DP1
        END IF
c
      else
        tau   = 1.d0
        dtau  = 0.d0
        d2tau = 0.d0
      end if
c
      f(4) = - cf4 * tau * dsqrt( znt**3 / (t * vol) )
c
      return
      end
