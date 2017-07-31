c---------------------------------------------------------------------
c---------------------------------------------------------------------
c consolidated non-relativistic version (one-in-all). WD 7/10/00
c---------------------------------------------------------------------
c---------------------------------------------------------------------
c
      program lvlfil
c
      implicit real*8 (a-h,o-z), integer*4 (i-n)
      real*8 nlev,ionh2,ionatm
c
c***********************************************************************
c     assemble molecule/atom/ion level data from ascii files to a      *
c     binary file. *
c***********************************************************************
c
      parameter ( mlev = 500, mi = 30 )
c
      common /data  /
     .               dat   (5*mlev )
      common /io    /
     .               iin, iout, itty, iold, inew, jobdat, jobtim
      common /names /
     .               dsno, dsnn, dsn
      character*8 dsno, dsnn, dsn
      character*40 name
c
      dimension h2dat(4*mlev), atmdat(5*mlev)
      dimension eh2   (mlev  ), ionh2 (mlev  ), gh2   (mlev  ),
     .          rh2   (mlev  )
      dimension eatm  (mlev  ), ionatm(mlev  ), gatm  (mlev  ),
     .          qnatm (mlev  ), qlatm (mlev  )
c
      equivalence ( atmdat, dat           ), ( h2dat , dat           )
      equivalence ( eh2   , dat           ), ( ionh2 , dat(  mlev+1) ),
     .            ( gh2   , dat(2*mlev+1) ), ( rh2   , dat(3*mlev+1) )
      equivalence ( eatm  , dat           ), ( ionatm, dat(  mlev+1) ),
     .            ( gatm  , dat(2*mlev+1) ), ( qnatm , dat(3*mlev+1) ),
     .            ( qlatm , dat(4*mlev+1) )
c
      dimension ei(mi), n1(mi), n2(mi), gi(mi),ionh24(mlev),ionat4(mlev)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     get date and time of job                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
ccc      jobdat = date ( )
ccc      jobtim = clock( )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     activate files.                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call opnfil
      xindf = -1.d200
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialize count of total number of levels                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ntot = 0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     first copy data from previous binary file (if any)               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     check for old file
      if ( iold .le. 0 ) go to 30
c
c     read specification of species
   11 read   ( iold )     nucz, jion, nlev, nion, ryd, name
c     write it on output file and binary file
      write  ( iout, 12 ) nucz, jion, nlev, nion, ryd, name
      write  ( inew )     nucz, jion, nlev, nion, ryd, name
   12 format ( i3, i2, 2i5, f10.2, a40 'data from previous binary file')
c
c     check for end of data
      if ( nucz .le. 0 ) go to 30
c
      if ( nucz .eq. 1 .and. (jion .eq. 1 .or. jion .eq. 2) ) then
c
c     h2 or h2+ molecule
c
c           read the data in
            read  ( iold ) h2dat
c
c           read energies of parent states
            read  ( iold ) ei
c
c           write all data to the new binary file
            write ( inew ) h2dat
            write ( inew ) ei
c
      else
c
c     atom or ion
c
c           read the data in
            read  ( iold ) atmdat
c
c           read energies of parent states
            read  ( iold ) ei
c
c           write all data to the new binary file
            write ( inew ) atmdat
            write ( inew ) ei
      end if
c
c     tally total number of levels
      ntot = ntot + int(nlev)
c
c     loop back for next set of data
      go to 11
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     now add new data from ascii files                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     read nuclear charge, ion stage, number of levels, number of
c     parent states, name of species
c
   30 read   ( iin , 31 ) nucz, jion, nlev4, nion, ryd, name
      nlev = dfloat(nlev4)
   31 format ( i3, i2, 2i5, f10.1, a40 )
c
c     check for end of data
      if ( nucz .le. 0 ) go to 60
c
c     h2 or h2+ molecule
      if ( nucz .eq. 1 .and. (jion .eq. 1 .or. jion .eq. 2) ) then
c
c           initialize all data to indefinites
            do 32 i = 1, 4 * mlev
            h2dat(i) = xindf
   32       continue
            do 33 i = 1, mi
            ei(i) = xindf
   33       continue
c
c           read the data in
            read   ( iin, 34 ) ( eh2(i), ionh24(i), gh2(i), rh2(i),
     .                            i = 1, nlev4 )
   34       format ( 5x, f10.1, i3, 6x, f7.2, e9.2 )
c
c           convert radii from angstroms to cm
            do 35 i = 1, nlev4
            ionh2(i) = dfloat(ionh24(i))
            rh2(i) = 1.0d-8 * rh2(i)
   35       continue
c
c           read energies of parent states
            read   (iin, 36 ) ( ei(i), i = 1, nion )
   36       format ( 5x, f10.1 )
c
c           write header to the output file
            write  ( iout, 31 ) nucz, jion, nlev4, nion, ryd, name
c
c           write all data to the binary file
            write  ( inew ) nucz, jion, nlev, nion, ryd, name
            write  ( inew ) h2dat
            write  ( inew ) ei
c
      else
c
c     atom or ion
c
c           initialize all data to indefinites
            do 50 i = 1, 5 * mlev
            atmdat(i) = xindf
   50       continue
            do 51 i = 1, mi
            ei(i) = xindf
   51       continue
c
c           read the data in
            read   ( iin, 52 ) ( eatm(i), qnatm(i), qlatm(i), gatm(i),
     .                         ionat4(i), i = 1, nlev4 )
   52       format ( 4x, f11.1, f8.4, f6.0, f7.0, i4 )
            do 5222 i=1,nlev4
 5222       ionatm(i) = ionat4(i)  
c
c           read energies of parent states and number of hydrogenic
c           states converging on them
            read   (iin, 53 ) ( ei(i), n1(i), n2(i), gi(i),
     .                          i = 1, nion )
   53       format ( 4x, f11.1, 2i5, f5.0  )
c
c           fill in upper hydrogenic energy levels converging on
c           parent states
            zsq = dfloat(jion * jion)
c           remember that charge of neutral hydrogen is a special case
            if ( nucz .eq. 1 .and. jion .eq. 4 ) zsq = 1.0d0
            do 55 ii = 1, nion
            if ( n1(ii) .le. 0 ) go to 55
            nmin = nlev4 + 1
            nmax = nmin + n2(ii) - n1(ii)
            do 54 i = nmin, nmax
            j = i - nmin + n1(ii)
            sq = j * j
            eatm  (i) = ei(ii) - zsq * ryd / sq
            qnatm (i) = j
            qlatm (i) = -1.0d0
            gatm  (i) = 2.0d0 * sq * gi(ii)
            ionatm(i) = dfloat(ii)
   54       continue
            nlev4 = nmax
            nlev = dfloat(nlev4)
   55       continue
c
c           write header to the output file
            write  ( iout, 31 ) nucz, jion, nlev4, nion, ryd, name
c
c           write all data to the binary file
            write  ( inew ) nucz, jion, nlev, nion, ryd, name
            write  ( inew ) atmdat
            write  ( inew ) ei
      end if
c
c     tally total number of levels
      ntot = ntot + nlev4
c
c     loop back for next set of data
      go to 30
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end of data                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     write end of data flag and end of file on binary file
   60 write  ( iout, 31 ) nucz, jion, nlev4, nion, ryd, name
      write  ( inew )     nucz, jion, nlev, nion, ryd, name
      end file inew
      write  ( itty, 61 ) ntot
      write  ( iout, 61 ) ntot
   61 format ( ' finished. ntot =' i6 )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     close and save files                                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call clzfil
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     terminate                                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      stop
      end
      subroutine clzfil
      implicit real*8 (a-h,o-z)
c***********************************************************************
c     close and save files                                             c
c***********************************************************************
c
      common /io    /
     .               iin, iout, itty, iold, inew, jobdat, jobtim
      common /names /
     .               dsno, dsnn, dsn
c
      character*8 dsno, dsnn, dsn
c
      nerror = 0
      iero   = 0
      iern   = 0
      if ( iold .gt. 0 ) then
            close( unit = iold, iostat = iero )
            nerror = nerror + iero
      end if
      close( unit = inew, iostat = iern )
      nerror = nerror + iern
      if ( nerror .gt. 0 ) then
            write  ( itty, 62 ) iero, iern
            write  ( iout, 62 ) iero, iern
   62       format ( ' dataset closing errors. iero =' i3 ' iern=' i3 )
      write(6,1111)iern,iern,iern
1111  format(1x,i8,2x,o8,2x,a8)
      end if
cccc  call mass( 'store', 0, 0, '\eos/atomdat/'//dsnn//'\', ier)
      write(6,*) 'here we would save ',dsnn,' to mass storage'
c
      return
      end
      subroutine opnfil
      implicit real*8 (a-h,o-z)
c***********************************************************************
c     open all files                                                   c
c***********************************************************************
c
      common /io    /
     .               iin, iout, itty, iold, inew, jobdat, jobtim
      common /names /
     .               dsno, dsnn, dsn
c
      character*8 dsno, dsnn, dsn
      logical lxst
c
c     declare unit numbers
c
      itty  = 4
      iin   = 5
      open( unit = iin, file = 'lvlfil.dat', status='old' )
      iout  = 6
c
c     open terminal
c
      istep = 1
      dsn   = 'tty'
      open( unit = itty, iostat = ier, file = dsn, status='unknown' )
      if ( ier .gt. 0 ) go to 5
c     read file names
c
ccccc      write  ( iout, 1 ) jobdat, jobtim
      read   ( iin , 2 ) iold, dsno
      write  ( iout, 3 ) iold, dsno
      read   ( iin , 2 ) inew, dsnn
      write  ( iout, 4 ) inew, dsnn
ccccc    1 format ( '1' 2a10 )
    2 format ( 5x, i5, 5x, a8 )
    3 format ( ' old'  i5, 5x, a8 )
    4 format ( ' new'  i5, 5x, a8 )
c
c     open old data file (if any)
c
      if ( iold .gt. 0 ) then
            dsn   = dsno
            inquire( file = dsn, exist = lxst )
            if ( .not. lxst ) then
                 istep = 5
ccccccccccc call mass ......
      write(6,*) 'call mass() attempted. force stop'
      ier=100
                 if ( ier .gt. 0 ) go to 5
            end if
            istep = 6
            open(   unit = iold , iostat = ier, file = dsn,
     .            status = 'old', form = 'unformatted' )
            if ( ier .gt. 0 ) go to 5
      end if
c
c     open new data file
c
      istep = 7
      dsn   = dsnn
      open( unit = inew, iostat = ier, file = dsn, status = 'unknown',
     .      form = 'unformatted' )
      if ( ier .le. 0 ) return
c
c     error returns
c
    5 write  ( itty, 6 ) istep, ier, ier, ier, dsn
      write  ( iout, 6 ) istep, ier, ier, ier, dsn
    6 format (' error return in step' i3 ' of file activation.'/
     .        ' ier =' i20, o22, a8, ' dsn =' a8 )
      stop 'lvlfil1'
c
c    7 write  ( itty, 8 ) dsn, istep
      write  ( iout, 8 ) dsn, istep
    8 format (' dataset ' a8 ' fails to exist. istep =' i3 )
      stop 'lvlfil2'
c
      end
