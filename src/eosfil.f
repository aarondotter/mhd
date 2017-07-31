c---------------------------------------------------------------------
c---------------------------------------------------------------------
c consolidated non-relativistic version (one-in-all). WD 8/8/2003
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
      program eosfil
      implicit real*8 (a-h,o-z), integer*4 (i-n)
      real*8 nlev,ionh2,ionatm,nlvls8
      character*8 hdr,dsna,dsne

c
c***********************************************************************
c     create binary mass store file containing physical and partition  c
c     function data required by equation of state code                 c
c***********************************************************************
c
c---------------- for coordination and remarks se file 'parms'
      parameter ( mchem =  15, mspes = 205, mz = 26, mion = mz + 1 )
      parameter ( mpf   = 200, mlev  = 500 )
c
      parameter ( mmap  = mchem * ( 2*mion + 2 ) + 2 * ( mspes + mpf )
     .                                + mz + 9 )
      parameter ( mpfdat = mpf * (4*mlev + 1) )
      parameter ( mspdat = 10 * mspes )
      parameter ( mi = 30, mspv = 50 )
c
      character*4 name(mchem)
      integer nucz(mchem), nion(mchem), atwt(mchem)
      
c      common /atoms /
c     .                         name  (mchem ),        nucz  (mchem ),
c     .                         nion  (mchem ),        atwt  (mchem )
      common /const /
     .                                camu   ,               cc     ,
     .                                ce     ,               ch     ,
     .                                ck     ,               cme    ,
     .                                cpi    ,               cevw   ,
     .                                ca0
      common /data  /
     .                         dat   (5*mlev)
      common /head  /
     .                         hdr   (  10  )
      common /io    /
     .                   iin, iout, itty, iatm, ieos, jobdat, jobtim
      common /names /
     .                  dsna, dsne
      common /map   /
     .                                nchem  ,               nspes  ,
     .                                npf    ,               ish2   ,
     .                                ish2p  ,               ishm   ,
     .                                ish    ,               ispr   ,
     .                                ise    ,        kz    (mz    ),
     .                         ichm1 (mchem ),        ichm2 (mchem ),
     .                         jspes (mspes ),        jpf   (mpf   ),
     .                         kspes (mspes ),        kpf   (mpf   ),
     .                  ispes (mion  ,mchem ), ipf   (mion  ,mchem )
      common /pfdat /
     .                         nlev  (mpf   ),
     .                  elev  (mlev  ,mpf   ), stwt  (mlev  ,mpf   ),
     .                  rad   (mlev  ,mpf   ), chi   (mlev  ,mpf   )
      common /spes  /
     .                         smass (mspes ),        gs    (mspes ),
     .                         zs    (mspes ),        zsq   (mspes ),
     .                         zmask (mspes ),        sn    (mspes ),
     .                         e0    (mspes ),        g0    (mspes ),
     .                         rad0  (mspes ),        chi0  (mspes )
c
      dimension imap(mmap), spdata(mspdat), pfdata(mpfdat), qn0(mspes)
      equivalence ( imap, nchem ), ( spdata, smass ), ( pfdata, nlev ),
     .            ( chi0, qn0   )
c
      dimension h2dat (4*mlev), atmdat(5*mlev), ei(mi)
      dimension eh2   (mlev  ), ionh2 (mlev  ), gh2   (mlev  ),
     .          rh2   (mlev  )
      dimension eatm  (mlev  ), ionatm(mlev  ), gatm  (mlev  ),
     .          qnatm (mlev  ), qlatm (mlev  )
      equivalence ( atmdat, dat           ), ( h2dat , dat           )
      equivalence ( eh2   , dat           ), ( ionh2 , dat(  mlev+1) ),
     .            ( gh2   , dat(2*mlev+1) ), ( rh2   , dat(3*mlev+1) )
      equivalence ( eatm  , dat           ), ( ionatm, dat(  mlev+1) ),
     .            ( gatm  , dat(2*mlev+1) ), ( qnatm , dat(3*mlev+1) ),
     .            ( qlatm , dat(4*mlev+1) )
c
      dimension kspv(mspv), ispv(mspv), rspv(mspv)
c
      data kspv
     .         /     1,     2,   6*6,   4*7,   3*8,    10,    12,  2*12,
     .            2*13,  7*14,  4*16,    18,  3*20,  3*20, 11*26/
      data rspv
     .         / 0.529d0,   0.5d0, 6*0.9d0, 4*.75d0, 3*0.6d0,   0.5d0,
     .           1.7d0, 2*2.1d0,  2*1.6d0, 7*1.4d0, 4*1.0d0,
     .           0.8d0, 3*2.1d0, 3*2.6d0, 11*1.3d0/
      data ispv
     .         /   1,   1,   1,   2,   3, 123, 124, 125,   1,   2,
     .             3, 123,   1,   2,   3,   1,   1,  68,  69,   1,
     .            37,   1,   2,   3, 115, 116, 117, 118,   1,   2,
     .             3,  76,   1,   1,  62,  63,  73,  74,  75,   1,
     .            90, 100, 113, 123, 133, 137, 147, 151, 152, 153/
c
      data           camu,            cc,            ce,            ch
     .    / 1.6605655d-24, 2.9979246d+10, 4.8032426d-10, 6.6261764d-27/
      data             ck,           cme,           cpi,          cevw
     .    / 1.3806624d-16, 9.1095345d-28, 3.141592654d0, 8.0654653d+03/
      data            ca0
     .    / 5.2917706d-09/
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     get date and time of job                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccccc      jobdat = date ( )
ccccc      jobtim = clock( )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     activate files.                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      call opnfil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read header                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      read   ( iin, 4 ) hdr
    4 format ( 10a8 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     zero indices in map                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do 5 im = 1, mmap
      imap(im) = 0
    5 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read atomic data and start construction of index map             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read number of chemical elements
      read   ( iin, 6 ) nchem
    6 format ( i5 )
c
c     loop over all chemical elements
      nspes = 1
      do 14 kchem = 1, nchem
c
c     read name, nuclear charge, number of molecule/atom/ion species,
c     atomic weight for element kchem.
      read   ( iin, 7 ) name(kchem), nucz(kchem), nion(kchem),
     .                  atwt(kchem)
    7 format ( 2x, a4, 2i5, f10.2 )
c
c     accumulate total number of particle species (including electrons).
      nspes = nspes + nion(kchem)
c
c     associate chemical element number kchem with charge number nucz
      kz( nucz(kchem) ) = kchem
c
c     define first and last particle species number for chemical element kchem
      if ( kchem .eq. 1 ) then
            ichm1(kchem) = 1
            ichm2(kchem) = nion(kchem)
      else
            ichm1(kchem) = ichm2(kchem - 1) + 1
            ichm2(kchem) = ichm2(kchem - 1) + nion(kchem)
      end if
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
c
c     define ion stage jspes and chemical element number kspes for particle
c     species is
      do 8 is = is1, is2
      jspes(is) = is - is1 + 1
      kspes(is) = kchem
c
c     define particle species number ispes for each molecule/atom/ion stage
c     of chemical element kchem
      ispes( is - is1 + 1, kchem ) = is
    8 continue
c
c     read ground-state statistical weights of all molecule/atom/ion stages
c     of element kchem
      read   ( iin, 9 ) ( g0(is), is = is1, is2 )
    9 format ( 15f5.1 )
c
c     read dissociation/detachment/ionization energies (in wavenumbers) of
c     all molecule/atom/ion species of chemical element kchem
      read   ( iin, 10 ) ( e0(is), is = is1, is2 - 1 )
   10 format ( 8f10.1 )
c     convert to electron volts
      do 11 is = is1, is2 -1
      e0(is) = e0(is) / cevw
   11 continue
c
c     read ground-state effective quantum numbers of all molecule/atom/ion
c     stages of element kchem
      read   ( iin, 12 ) ( qn0(is), is = is1, is2 -1 )
   12 format ( 8f10.1 )
c
c     read ground-state radii (angstroms) of neutral atoms and h2, h2+, and h-
      if ( nucz(kchem) .eq. 1 ) then
c     hydrogen
            read   ( iin, 13 ) ( rad0(is), is = is1, is2 - 1 )
   13       format ( 8f10.1 )
      else
c     all other elements
            read   ( iin, 13 ) rad0(is1)
      end if
   14 continue
c
c     define species numbers of electrons and of hydrogen molecules/atoms/ions
      ise = nspes
      kh  = kz(1)
      if ( kh .gt. 0 ) then
            ish2  = ichm1(kh)
            ish2p = ish2  + 1
            ishm  = ish2p + 1
            ish   = ishm  + 1
            ispr  = ish   + 1
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     print input data                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write  ( iout, 20 )
   20 format ( '1' 15x 'data preparation for equation of state' )
      write  ( iout, 21 )
   21 format ( ' ' )
      write  ( iout, 22 ) hdr, jobdat, jobtim
   22 format ( 10a8/ 2a10 )
      write  ( iout, 23 ) dsna
   23 format ( ' atomdata: ' a8 )
      write  ( iout, 24 ) dsne
   24 format ( '  eosdata: ' a8 )
      write  ( iout, 21 )
      write  ( iout, 25 ) nchem
   25 format ( ' gas consists of' i3 ' chemical elements' )
c
      do 35 kchem = 1, nchem
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
      write  ( iout, 21 )
      write  ( iout, 21 )
      write  ( iout, 26 ) name(kchem), nucz(kchem), atwt(kchem)
   26 format ( ' ' a10 ', z =' i3 ', atomic weight =' f8.4 )
      write  ( iout, 21 )
      write  ( iout, 27 )
   27 format ( 5x ' g0:' )
      write  ( iout, 28 ) ( g0(is), is = is1, is2 )
   28 format ( 5x,10f7.1 )
      write  ( iout, 29 )
   29 format ( 5x ' e0:' )
      write  ( iout, 30 ) ( e0(is), is = is1, is2 - 1 )
   30 format ( 5x, 10f7.1 )
      write  ( iout, 31 )
   31 format ( 5x ' qn0:' )
      write  ( iout, 32 ) ( qn0(is), is = is1, is2 - 1 )
   32 format ( 8x,10f7.4 )
      write  ( iout, 33 )
   33 format ( 5x ' rad0:' )
      if ( nucz(kchem) .eq. 1 ) then
c     hydrogen
            write  ( iout, 34 ) ( rad0(is), is = is1, is2 - 1 )
   34       format ( 5x,10f7.2 )
      else
c     all other elements
            write  ( iout, 34 ) rad0(is1)
      end if
   35 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute all physical data required by eos code for each species  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      gcon = ( 2.0d0 * cpi * ck / ch**2 ) ** 1.5d0
c
c     free electrons
      zs   (ise) = - 1.0d0
      zsq  (ise) =   1.0d0
      zmask(ise) =   1.0d0
      smass(ise) =   cme
      gs   (ise) =   gcon * cme**1.5d0
c
c     all other species
      do 46 kchem = 1, nchem
c
c     mass, charge, and translational statistical weight factor gs
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
      if ( nucz(kchem) .eq. 1 ) then
c           hydrogen
            zs(ish2 ) =   0.0d0
            zs(ish2p) =   1.0d0
            zs(ishm ) = - 1.0d0
            zs(ish  ) =   0.0d0
            zs(ispr ) =   1.0d0
            smass(ish  ) = atwt(kchem) * camu
            smass(ispr ) = smass(ish) - cme
            smass(ishm ) = smass(ish) + cme
            smass(ish2 ) = 2.0d0 * smass(ish)
            smass(ish2p) = smass(ish2) - cme
      else
c           all other elements
            do 40 is = is1, is2
            zs   (is) = dfloat( is - is1 )
            smass(is) = atwt(kchem) * camu - zs(is) * cme
   40       continue
      end if
      do 41 is = is1, is2
      zsq  (is) = zs(is)**2
      zmask(is) = cvmgz( 0.0d0, 1.0d0, zs(is) )
      gs   (is) = gcon * smass(is)**1.5d0
   41 continue
c
c     stark-ionization theory parameter for ground states
      do 42 is = is1, is2 - 1
      fac = cvmgt(
     .     16.0d0 *     qn0(is) **2 * (qn0(is) + 1.16666666666666d0) /
     .  ( 3.0d0 * (1.0d0 + qn0(is))**2*(qn0(is)**2 + qn0(is) + 0.5d0) ),
     .        1.0d0,  qn0(is) .ge. 3.0d0 )
      zz  = dmax1( zs(is) + 1.0d0, 1.0d0 )
      e0(is) = cc * ch * ( cevw * e0(is) )
      chi0(is) = 16.0d0*( ce**2 * dsqrt( zz/fac) / e0(is) )**3
   42 continue
c
c     absolute ground-state energies (ergs) of all molecule/atom/ion
c     stages of element kchem relative to ground state of lowest stage
      if ( nucz(kchem) .eq. 1 ) then
c     hydrogen
            dh2  = e0(is1    )
            dh2p = e0(is1 + 1)
            dhm  = e0(is1 + 2)
            dh   = e0(is1 + 3)
            e0(is1    ) = 0.0d0
            e0(is1 + 1) =         dh2 + dh - dh2p
            e0(is1 + 2) = 0.5d0 * dh2      - dhm
            e0(is1 + 3) = 0.5d0 * dh2
            e0(is1 + 4) = 0.5d0 * dh2 + dh
      else
c     all other elements
            do 43 is = is2, is1 + 1, -1
            e0(is) = e0(is - 1)
   43       continue
            e0(is1) = 0.0d0
            do 44 is = is1 + 1, is2
            e0(is) = e0(is - 1) + e0(is)
   44       continue
      end if
c
c     convert ground-state radii to cm
      if ( nucz(kchem) .eq. 1 ) then
c     hydrogen
            do 45 is = is1, is2 - 1
            rad0(is) = 1.0d-8 * rad0(is)
   45       continue
      else
c     all other elements
            rad0(is1) = 1.0d-8 * rad0(is1)
      end if
   46 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read partition function table information; finish index map      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     number of tables
      read   ( iin, 50 ) npf
   50 format ( i5 )
      write  ( iout, 51 )
   51 format ( ' ' )
      write  ( iout, 52 ) npf
   52 format (' multistate partition functions used for' i4 ' species' )
c
c     loop over all tables
      do 56 ip = 1, npf
c
c     read nuclear charge and ion stage for table number ip
      read   ( iin, 53 ) kp, jp
   53 format ( 2i5 )
      write  ( iout, 54 ) ip, kp, jp
   54 format ( i5, 5x ' z =' i3 ', ion stage =' i3 )
c
c     check that we have the element in mix
      kp = kz(kp)
      if ( kp .le. 0 ) then
      write  ( iout, 55 )
   55 format (' partition function designated for element not in mix' )
      stop 'eosfil3'
      end if
c
c     define ion number jpf and chemical element number kpf for partition
c     function number ip
      jpf(ip) = jp
      kpf(ip) = kp
c
c     define partition function table number ipf for ion species jp of
c     chemical element kp
      ipf(jp, kp) = ip
   56 continue
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write parameters to binary file                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write  ( ieos ) mchem, mz, mion, mspes, mpf, mlev 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write index map to binary file                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write  ( ieos ) imap
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c    now deal with partition function data                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     initialize count
      npff = 0
c
c     read next species
   60 read   ( iatm ) kp, jp, nlvls8
      nlvls = nint(nlvls8)
c
c     check for end of data
      if ( kp .eq. 0 ) go to 70
c
c     check that data for this species is needed. if not, skip it.
      kp = kz(kp)
      is = ispes(jp, kp)
      ip = ipf  (jp, kp)
      if ( ip .le. 0 ) then
            read   ( iatm ) null
            read   ( iatm ) null
            go to 60
      end if
c
c     increment count
      npff = npff + 1
c
c     hydrogen molecule
      if ( is .eq. ish2 ) then
            read   ( iatm ) h2dat
            read   ( iatm ) ei
            nlev(ip) = nlvls8
c           energy levels (ergs), statistical weights, radii
            do 61 il = 1, nlvls
            elev(il, ip) = eh2(il) * cc * ch
c---------------------- output to show place of energy levels
c     write(6,*)'H2: ip,il,elev(il,ip) ',ip,il,elev(il,ip) 
c------------------------------------------------------------
            stwt(il, ip) = gh2(il)
            rad (il, ip) = rh2(il)
   61       continue
            rad0(ish2) = rad(1, ip)
c           stark ionization parameter
            zz  = 1.0d0
            do 62 il = 1, nlvls
            j = nint(ionh2(il))
            qnm = ionh2(il)
            eh2(il) = cc * ch * ei(j)
            fac =
     .      cvmgt(16.0d0*  qnm **2   * (qnm + 1.16666666666666d0) /
     .           ( 3.0d0 * (1.0d0 + qnm)**2 * ( qnm**2 + qnm + 0.5d0) ),
     .                1.0d0,  qnm .ge. 3.0d0 )
            chi(il, ip) = 16.0d0*( ce**2*dsqrt( zz/fac ) / eh2(il) )**3
   62       continue
            chi0(ish2) = chi(1, ip)
c
c     hydrogen molecule ion
      else if ( is .eq. ish2p ) then
            read   ( iatm ) h2dat
            read   ( iatm ) ei
            nlev(ip) = nlvls8
            zz  = 2.0d0
            fac = 1.0d0
            con = 16.0d0*(ce**2*dsqrt(zz/fac) / ( cc*ch*ei(1) ) )**3
c           energy levels (ergs), statistical weights, radii,
c           and stark ionization parameter
            do 63 il = 1, nlvls
            elev(il, ip) = eh2(il) * cc * ch
c---------------------- output to show place of energy levels
c     write(6,*)'H2+: ip,il,elev(il,ip) ',ip,il,elev(il,ip) 
c------------------------------------------------------------
            stwt(il, ip) = gh2(il)
            rad (il, ip) = rh2(il)
            chi (il, ip) = con
   63       continue
            rad0(ish2p) = rad(1, ip)
            chi0(ish2p) = chi(1, ip)
c
c     an atom or ion
      else
            read   ( iatm ) atmdat
            read   ( iatm ) ei
            nlev(ip)= nlvls8
c           energy levels (ergs), statistical weights
            do 64 il = 1, nlvls
            elev(il, ip) = eatm(il) * cc * ch
c---------------------- output to show place of energy levels
c     write(6,*)'ion/atom: ip,il,elev(il,ip) ',ip,il,elev(il,ip) 
c------------------------------------------------------------
            stwt(il, ip) = gatm(il)
   64       continue
c           stark ionization parameter
            zz  = dfloat( jp )
            if ( is .eq. ish ) zz = 1.0d0
            do 65 il = 1, nlvls
            j = nint(ionatm(il))
            eatm(il) = ( ei(j) - eatm(il) ) * cc * ch
            fac = cvmgt(
     .      16.0d0* qnatm(il) **2 * (qnatm(il) + 1.16666666666666d0) /
     .     (3.0d0*(1.d0+qnatm(il))**2*(qnatm(il)**2+qnatm(il)+.5d0)),
     .       1.0d0,  qnatm(il) .ge. 3.0d0 )
            chi(il, ip) = 16.0d0*( ce**2 *dsqrt(zz/fac) / eatm(il) )**3
   65       continue
            chi0(is) = chi(1, ip)
c           radii (neutral species only)
            if ( zs(is) .eq. 0.0d0 ) then
c                 first use hydrogenic formula
                  do 66 il = 1, nlvls
                  qlatm(il) =
     .                    cvmgp( qlatm(il), qnatm(il)-1.0d0, qlatm(il) )
                  rad(il, ip) = 0.5d0 * ca0 * ( 3.0d0 * qnatm(il)**2
     .                        - qlatm(il) * ( qlatm(il) + 1.0d0) ) / zz
   66             continue
c                 then scan for special cases
                  do 67 iv = 1, mspv
                  if ( kspv(iv) .eq. nucz(kp) )
     .                              rad(ispv(iv), ip) = 1.d-8 * rspv(iv)
   67             continue
                  rad0(is) = rad(1, ip)
            end if
      end if
c
c     check for completion of partition function data
      if ( npff .ne. npf ) go to 60
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     check that we got all the data requested                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   70 if ( npff .ne. npf ) then
            write  ( iout, 71 ) npf, npff
   71       format (' number of partition functions requested (' i3
     .              ') disagrees with number found (' i3 ')' )
            stop'eosfil4'
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write physical data for all species                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write  ( ieos ) spdata
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write partition function data                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write  ( ieos ) pfdata
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close and save files                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call clzfil
c
      stop'normal termination'
      end
      subroutine clzfil
      implicit real*8 (a-h,o-z), integer*4 (i-n)
c ----- open and close statements require plain-vanilla integer variables
      integer iera,iere
c***********************************************************************
c     close and save files                                             c
c***********************************************************************
c
      common /io    /
     .                   iin, iout, itty, iatm, ieos, jobdat, jobtim
      common /names /
     .                  dsna, dsne
c
      character*8 dsna, dsne, dsnn
c
      nerror = 0
      close( unit = iatm, iostat = iera )
      nerror = nerror + iera
      close( unit = ieos, iostat = iere )
      nerror = nerror + iere
      if ( nerror .gt. 0 ) then
            write  ( itty, 62 ) iera, iere
            write  ( iout, 62 ) iera, iere
   62       format ( ' dataset closing errors. iera =' i3 ' iere =' i3 )
      end if
      dsnn=dsne
      ier=0
ccccc call save( ier, 'dn' , dsnn , 'pdn' , dsnn )
      if ( ier .gt. 0 ) then
            write  ( itty, 63 ) ier, dsne
            write  ( iout, 63 ) ier, dsne
   63       format ( ' error number' i3 ' in storing dataset ' a8
     .               ' in mass store' )
      end if
c
      return
      end
      subroutine opnfil
      implicit real*8 (a-h,o-z), integer*4 (i-n)
c ----- open and close statements require plain-vanilla integer variables
      integer ier
c***********************************************************************
c     open all files                                                   c
c***********************************************************************
c
      common /io    /
     .                   iin, iout, itty, iatm, ieos, jobdat, jobtim
      common /names /
     .                  dsna, dsne
c
      character*8 dsn, dsna, dsne
c
c     declare unit numbers
c
      itty  = 4
      iin   = 5
      open( unit = iin, file = 'eosfil.dat', status='old')
      iout  = 6
      iatm  = 7
      ieos  = 8
c
c     open terminal
c
      istep = 1
      dsn   = 'tty'
      open( unit = itty, iostat = ier, file = dsn)
      if ( ier .gt. 0 ) go to 5
c
c     read file names
c
      write  ( itty, 1 ) jobdat, jobtim
      read   ( iin , 2 ) dsna
      write  ( itty, 3 ) dsna
      read   ( iin , 2 ) dsne
      write  ( itty, 4 ) dsne
    1 format ( '1' 2a10 )
    2 format ( 10x, a8 )
    3 format ( ' atomdata: ' a8 )
    4 format ( '  eosdata: ' a8 )
c
c     open atomic data file
c
      dsn  = dsna
      istep = 5
ccccc call access(ier,'dn',dsn,'pdn',dsn)
      if(ier .gt. 0) goto 5
c
      istep = 6
      open( unit = iatm , iostat = ier, file = dsn, status = 'old',
     .      form = 'unformatted' )
      if ( ier .gt. 0 ) go to 5
c
c     open eos data file
c
      istep = 7
      dsn = dsne
      open( unit = ieos, iostat = ier, file = dsn, status = 'unknown',
     .      form = 'unformatted' )
      if ( ier .le. 0 ) return
c
c     error returns
c
    5 write  ( itty, 6 ) istep, ier, ier, ier, dsn
      write  ( iout, 6 ) istep, ier, ier, ier, dsn
    6 format (' error return in step' i3 ' of file activation.'/
     .        ' ier =' i20, o22, a8, ' dsn =' a8 )
      stop 'eosfil1'
c
c    7 write  ( itty, 8 ) dsn, istep
c      write  ( iout, 8 ) dsn, istep
c    8 format (' dataset ' a8 ' fails to exist. istep =' i3 )
c      stop 'eosfil2'
c
      end
