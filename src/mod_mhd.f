c---------------------------------------------------------------------
c consolidated relativistic version (one-in-all). WD 7/10/00
c---------------------------------------------------------------------
      module MHDeos

      implicit double precision (a-h,o-z), integer (i-n)

      integer, parameter :: nres=23

      private

      public :: mhd_init, eosDT_get, nres, write_result
      
      contains

      subroutine mhd_init(datafile,abunfile)
      character(len=128) :: datafile, abunfile
      call openfile(datafile,abunfile)

      call setup
      
      end subroutine mhd_init

      subroutine eosDT_get(logT,logRho,result)
      double precision, intent(in) :: logT(:), logRho(:)
      double precision, intent(out) :: result(:,:)

      call track(1,logT,logRho)
      call tabdef

      call tdfill

      call outtab(result)
          
      end subroutine eosDT_get
      
      subroutine write_result(io,tl,res)
      integer, intent(in) :: io
      double precision, intent(in) :: tl(:), res(:,:)
      integer :: i,n
      
      write(io,'(99a16)') 'logRho', 'logPgas', 'logE', 'logS',
     >     'dlnPgas_dlnT', 'dlnPgas_dlnRho', 'd2lnPgas_dlnT_dlnRho',
     >     'dlnE_dlnT', 'dlnE_dlnRho', 'd2lnE_dlnT_dlnRho', 'dlnS_dlnT',
     >     'dlnS_dlnRho', 'd2lnS_dlnT_dlnRho', 'mu', 'log_free_e', 
     >     'eta', 'f_H+', 'f_He+', 'f_He++', 'f_H2', 'dse', 'dpe', 'dsp'

      n=size(tl)
      do i=1,n
         write(io,'(1p99e16.8)') tl(i), res(i,:)
      enddo
      end subroutine write_result

      
      subroutine track(nrho0,logT,logRho)
      integer, intent(in) :: nrho0
      double precision, intent(in) :: logT(:), logRho(:)
c  sets up temperature-density part of mhd driver
c  the input-model-specific part is handled in s/r track1
c
c  original version: 3/3/88(aarhus) and 31/3/88 (circe,paris)
c  upgraded for gamma1-dependent flexible meshing: 6/6/88 (circe paris)
c
      common/tr/tlin0(1000),rlin0(1000),gamin0(1000),lines0
      common/gr/tl(1000),rhmin(1000),nrhot,isocom
      dimension rl(1000)

      lines0 = size(logT)
c
      tl(1:lines0) = logT
      rl(1:lines0) = logRho
      
      iso = lines0
      isocom = iso
      
      rhmin = rl
      nrhot = nrho0
      
      return
      end
      

      subroutine tabdef
      include 'tabparms'
      common/gr/tl(1000),rhmin(1000),nrhot,iso
      common/tab/tl2(nt2m),rl2(nt2m,nr2m),nt2,nr2

      nt2 = iso
      nr2 = nrhot

      do it=1,nt2
         tl2(it)=tl(it)
         rl2(it,nr2) = rhmin(it)
      enddo

      return
      end

      
      subroutine openfile(datafile,abunfile)
      character(len=128) :: datafile, abunfile
c***********************************************************************
c     opens i/o files. acquires requested datasets                     c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
      iin    = 5
      ipunch = 1
      idat   = 2

      open(unit=iin,file=trim(abunfile),iostat=ier,status='old')
      if (ier>0) go to 4
      
      open(unit=ipunch,iostat=ier,status='unknown',form='unformatted')
      if (ier>0) go to 4

      open(unit=idat,file=trim(datafile),iostat=ier,status='old',
     >     form='unformatted')
      if(ier>0) go to 4
      
      
      return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     error returns                                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 4    write  ( iout, 5 ) istep, ier, ier, ier, dsn, pdn
 5    format ( ' error return in step',i3,' of file activation'/
     .         ' ier =',i20, o22, a8/ ' dsn = ',a8,' pdn = ',a8 )
      stop 'opnfil1'

      end

      
      subroutine tdfill
      include 'types'
      include 'parms'
      include 'tabparms'
      include 'coms'

      common/tab/tl2(nt2m),rl2(nt2m,nr2m),nt2,nr2
      common/eq/var(ivar)
      common/out/var2(nt2m,nr2m,ivar)

      ddt = 1.0d-4
      ddr = 1.0d-4

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nt     = nt2
      nrho   = nr2
c>>>>>>>>>>>>>>>>>> nt,drho,nrho transfered via common
c
      icase=1
      call tdtab(icase)
c
      return
      end

      

      subroutine outtab(result)
      integer ier
      double precision :: result(:,:)
c
c
c....... create table on (formatted) output file
c
      include 'parms'
      include 'tabparms'

      common/tab/tl2(nt2m),rl2(nt2m,nr2m),nt2,nr2
      common/out/var2(nt2m,nr2m,ivar)
c
c>>>>>>>>>>>>>>>>>>>>>> only part of mhd common needed <<<<<<<<<<<<<<<<
      common /io    /
     .                                iin    ,               iout   ,
     .                                ipunch ,               ifracs ,
     .                                ipops  ,               idat

c>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>><<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
c
      dimension atwt (mchem),abun (mchem)
c
c....................................................................
c
      rewind ipunch
c
      read(ipunch) nchem,(atwt(ic),abun(ic),ic=1,nchem),gasmu
      
      read ( ipunch, iostat = ier ) nt, nr, drho, ddt, ddr
c======================================================================
      call rpun(nt2,nr2,drho,ddt,ddr)
c======================================================================
c     
      do i=1,nt2
         do j=1,nr2
            result(i,:) = var2(i,j,:)
         enddo
      enddo
      
      return
      end

      
      subroutine rpun(nt,nrho,drho,ddt,ddr)
      include 'parms'
      include 'tabparms'

      common/tab/tl2(nt2m),rl2(nt2m,nr2m),nt2,nr2
      common/out/var2(nt2m,nr2m,ivar)
      common /io/iin,iout,ipunch,ifracs,ipops,idat
      common/rp1/ tlg(5),rhomin(5),elg(mrho,5),
     & pglg(mrho,5),pelg(mrho,5),ptlg(mrho,5),etlg(mrho,5),csbv(mrho,5),
     & csbp(mrho,5),chirh(mrho,5),cht(mrho,5),gm1(mrho,5),gm2(mrho,5),
     & gm3(mrho,5),qdb(mrho,5),frp1(mrho,5),frp2(mrho,5),frp3(mrho,5),
     & frp4(mrho,5),elpun(mrho,5),etapun(mrho,5),
     & eltpun(mrho,5),elrpun(mrho,5),tnpun(mrho,5),varspc(mvspc,mrho,5),
     & iloopi(5)
c ------ varspc for various quantities (with the possibility
c ------ to obtain derivatives, therefore here and not in s/r outp0-3).
c
      do n=1,nt
         do j=1,5
            read( ipunch )      nrhtab,iloopi(j),tlg(j),rhomin(j),
     .             (elg  (i,j), i= 1, nrho), (pglg (i,j), i= 1, nrho),
     .             (pelg (i,j), i= 1, nrho), (ptlg (i,j), i= 1, nrho),
     .             (etlg (i,j), i= 1, nrho), (csbv (i,j), i= 1, nrho),
     .             (csbp (i,j), i= 1, nrho), (chirh(i,j), i= 1, nrho),
     .             (cht  (i,j), i= 1, nrho), (gm1  (i,j), i= 1, nrho),
     .             (gm2  (i,j), i= 1, nrho), (gm3  (i,j), i= 1, nrho),
     .             (qdb  (i,j), i= 1, nrho),
     .             (frp1 (i,j), i= 1, nrho), (frp2 (i,j), i= 1, nrho),
     .             (frp3 (i,j), i= 1, nrho), (frp4 (i,j), i= 1, nrho),
     .             (elpun(i,j), i= 1, nrho), (etapun(i,j),i= 1, nrho),
     .             (eltpun(i,j),i= 1, nrho), (elrpun(i,j),i= 1, nrho),
     .             (tnpun (i,j),i= 1, nrho),
     .             ( (varspc(ivspc,i,j), ivspc=1,mvspc),  i= 1, nrho)

         enddo

         if(nrhtab.ne.nrho) stop ' mismatch of density points in rpun '
         
         call rpun1(n,nrho,var2,tl2,nt2m,nr2m,ivar,drho,ddt,ddr)
      enddo

      return
      end

      
      subroutine rpun1(n,nrho,var,tl,ntm,nrm,ivar1,drho,ddt,ddr)
      include 'parms'
      include 'tabparms'

      double precision :: mu, muinv, mu_e, mu_i, ne, nion, deltaT,
     >     deltaR
      
c ----------- ivar from tabparms is dummy variable. effective is ivar1.
c
      dimension var(ntm,nrm,ivar1),tl(ntm)
c
      common/rp1/ tlg(5),rhomin(5),elg(mrho,5),
     & pglg(mrho,5),pelg(mrho,5),ptlg(mrho,5),etlg(mrho,5),csbv(mrho,5),
     & csbp(mrho,5),chirh(mrho,5),cht(mrho,5),gm1(mrho,5),gm2(mrho,5),
     & gm3(mrho,5),qdb(mrho,5),frp1(mrho,5),frp2(mrho,5),frp3(mrho,5),
     & frp4(mrho,5),elpun(mrho,5),etapun(mrho,5),
     & eltpun(mrho,5),elrpun(mrho,5),tnpun(mrho,5),varspc(mvspc,mrho,5),
     & iloopi(5)
c ------ varspc for various quantities (with the possibility
c ------ to obtain derivatives, therefore here and not in s/r outp0-3).
c
      umod   = log(10.d0)
      carad = 7.56567d-15
      camu = 1.6605655d-24
c
      tlog = tlg(2)
      T = exp(umod*tlog)
c
      if(dabs(tlog-tl(n)).gt.1.d-03)then
         write(*,*) 'tlog=', tlog
         write(*,*) 'tlg(2)=', tlg(2)
         write(*,*) 'tl(n)=', tl(n)
         stop 'rpun1: mismatch in temps'
      endif

c
      deltaT = 2.0d0*ddt
      deltaR = 2.0d0*ddr
      
      do m=1,nrho
c     
c............ quantities for table ................
c
         etot        = varspc(1,m,2)
         ftot        = varspc(2,m,2)
         stot        = varspc(5,m,2) !entropy
         sl          = log10(max(stot,1e-20))
         rhol = rhomin(2) + dfloat(m-1)*drho
         rho = exp(umod*rhol)
         
         ne = varspc(3,m,2)
         nion = varspc(4,m,2)
         mu_e = rho/camu/ne
         mu_i = rho/camu/nion
         muinv = (1.0d0/mu_i + 1.0d0/mu_e)
         mu = 1.0d0/muinv
         pgas = exp(umod*pglg(m,2))
         prad = carad*T*T*T*T/3.0d0
         P = Pgas + Prad

         chiT = cht(m,2)
         chiRho = chirh(m,2)

         dlnP_dlnT = (ptlg(m,5) - ptlg(m,4)) / deltaT
         dlnP_dlnRho = (ptlg(m,3) - ptlg(m,1)) / deltaR
         dlnS_dlnT = (log10(varspc(5,m,5))-log10(varspc(5,m,4)))/deltaT
         dlnS_dlnRho=(log10(varspc(5,m,3))-log10(varspc(5,m,1)))/deltaR

         !chiT = dlnP_dlnT
         !chiRho = dlnP_dlnRho

         !print *, cht(m,2), chiT, chirh(m,2), chiRho
         
         dlnE_dlnT = (etlg(m,5) - etlg(m,4))/deltaT
         dlnE_dlnRho = (etlg(m,3) - etlg(m,1))/deltaR

         dS_dT_Rho = (stot/T)*dlnS_dlnT
         dS_dRho_T = (stot/rho)*dlnS_dlnRho
         dE_dT_Rho = (etot/T)*dlnE_dlnT
         dE_dRho_T = (etot/rho) * dlnE_dlnRho
         dP_dT_Rho = (P/T)*dlnP_dlnT

         !Frank's Maxwell relations
         dse = T*dS_dT_Rho/dE_dT_Rho -1.0d0
         
         dpe = (rho*rho*dE_dRho_T + T*dP_dT_Rho)/P - 1.0d0        

         dsp = -(dS_dRho_T*rho*rho)/dP_dT_Rho - 1.0d0
         
         var(n,m, 1) = rhol
         var(n,m, 2) = pglg(m,2)
         var(n,m, 3) = etlg(m,2)
         var(n,m, 4) = sl
         var(n,m, 5) = chiT
         var(n,m, 6) = chiRho
         var(n,m, 7) = 0d0 !d2lnPgas_dlnT_dlnRho
         var(n,m, 8) = dlnE_dlnT
         var(n,m, 9) = dlnE_dlnRho
         var(n,m,10) = 0d0 !d2lnE_dlnT_dlnRho
         var(n,m,11) = dlnS_dlnT
         var(n,m,12) = dlnS_dlnRho
         var(n,m,13) = 0d0 !d2lnS_dlnT_dlnRho
         var(n,m,14) = mu
         var(n,m,15) = 0d0 !log_free_e
         var(n,m,16) = etapun(m,2)
         var(n,m,17) = frp1(m,2)
         var(n,m,18) = frp2(m,2)
         var(n,m,19) = frp3(m,2)
         var(n,m,20) = frp4(m,2)
         var(n,m,21) = dse
         var(n,m,22) = dpe
         var(n,m,23) = dsp
      enddo

      return
      end


c======================================================================
c
c follow those s/r that were significantly modified for the
c high-energy fixup
c
c======================================================================
c
      subroutine setup
c***********************************************************************
c     initialization.  read and validate data.                         c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
      dimension imap(mmap), jmap(mmap), spdata(mspdat), pfdata(mpfdat)

      equivalence ( imap, nchem ), ( spdata, smass ), ( pfdata, nlev )

      data           camu,            cc,            ce,            ch
     .    / 1.6605655d-24, 2.9979246d+10, 4.8032426d-10, 6.6261764d-27/
      data             ck,           cme,           cpi,          cevw
     .    / 1.3806624d-16, 9.1095345d-28, 3.141592654d0, 8.0654653d+03/

      data          carad
     .    / 7.56567  d-15/
ccc   data          carad
ccc  .    / 0.0 d0 /
c
      ifailt = 0
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read header                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      read(iin,*)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     zero all indices in map                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      imap = 0
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read atomic data                                                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     number of chemical elements
      read   ( iin, 2 ) nchem
    2 format ( i5 )
c
c     read name, nuclear charge, number of molecule/atom/ion species,
c     atomic weight, abundance.  compute total number of particle
c     species (including electrons).
      nspes = 1
      do kchem = 1, nchem
         read   ( iin, '(2x,a4,2i5,f10.2,e10.2)' )
     >        name(kchem), nucz(kchem), nion(kchem),
     >        atwt(kchem), abun(kchem)

         write(6,*) 'kchem, atwt, abun =',kchem, atwt(kchem),
     *        abun(kchem)
         nspes = nspes + nion(kchem)
      enddo
      
c     normalize abundances
      abun = abun/sum(abun(1:nchem))

c     mean molecular weight
      gasmu = dot_product(abun(1:nchem),atwt(1:nchem))

c
c     set ionization potential (ergs) and abundance for fictitious
c     representative metal
      chim  = cc * ch * ( 5.5d0 * cevw )
      abunm = 1.0d-10
      do kchem = 1, nchem
         if ( nucz(kchem) .eq. 11 ) abunm = abunm + abun(kchem)
         if ( nucz(kchem) .eq. 13 ) abunm = abunm + abun(kchem)
         if ( nucz(kchem) .eq. 19 ) abunm = abunm + abun(kchem)
         if ( nucz(kchem) .eq. 20 ) abunm = abunm + abun(kchem)
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     construct index map                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     associate chemical element number kchem with charge number nucz
c
c ..... initialize to zero to prevent printouts of nonexisiting
c ..... elements in s/r outp1 .... never assume automatic
c ..... preinitialization .....                             [7/7/00, WD]
      do jz=1,mz
         kz (jz) = 0
      enddo
c
      do kchem = 1, nchem
         kz( nucz(kchem) ) = kchem
      enddo
c
      do kchem = 1, nchem
c     first and last particle species number for chemical element kchem
         if ( kchem .eq. 1 ) then
            ichm1(kchem) = 1
            ichm2(kchem) = nion(kchem)
         else
            ichm1(kchem) = ichm2(kchem-1) + 1
            ichm2(kchem) = ichm2(kchem-1) + nion(kchem)
         end if
c     
         do is = ichm1(kchem), ichm2(kchem)
c     
c     ion number jspes and chemical element number kspes for particle
c     species is
            jspes(is) = is - ichm1(kchem) + 1
            kspes(is) = kchem
c     
c     particle species number ispes for ion species j of chemical elemen
            ispes( is - ichm1(kchem) + 1, kchem ) = is
         enddo
      enddo
c
c     species numbers of electrons and of hydrogen molecules/atoms/ions
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
c     read partition function table information                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     number of tables
      read   ( iin, '(i5)' ) npf
      write  ( iout,  '(a1)' ) '0'
      write  ( iout, 17 ) npf
 17   format (' partition function tables used for',i4,' species:' )
c
c     loop over all tables
      do ip = 1, npf
c
c     nuclear charge and ion stage for table ip
         read   ( iin, '(2i5)') kp, jp

c     convert to chemical element number
         kp = kz(kp)

c     check that we have that element in mix
         if ( kp < 0 ) stop 'setup1'

c     ion number jpf and chemical element number kpf for partition
c     function number ip
         jpf(ip) = jp
         kpf(ip) = kp
         
c     partition function table number ipf for ion species jp of chemical
c     element kp
         ipf(jp, kp) = ip
      enddo      
c

      close(iin)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     validate parameters and index map                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     validate parameters
      read   ( idat ) mmchem, mmz , mmion, mmspes, mmpf, mmlev
      if ( mmchem .eq. mchem ) go to 31
      write  ( iout, 30 ) mmchem, mchem
 30   format(2i5)
      stop 'setup2'
 31   if ( mmz    .eq. mz    ) go to 32
      write  ( iout, 30 ) mmz   , mz
      stop 'setup3'
 32   if ( mmion  .eq. mion  ) go to 33
      write  ( iout, 30 ) mmion , mion
      stop 'setup4'
 33   if ( mmspes .eq. mspes ) go to 34
      write  ( iout, 30 ) mmspes, mspes
      stop 'setup5'
 34   if ( mmpf   .eq. mpf   ) go to 35
      write  ( iout, 30 ) mmpf  , mpf
      stop 'setup6'
 35   if ( mmlev  .eq. mlev  ) go to 40
      write  ( iout, 30 ) mmlev , mlev
      stop 'setup7'
c
c     validate index map
 40   read(idat) jmap
      do im = 1, mmap
         if ( imap(im) /= jmap(im) ) then
            write  ( iout, '(3i5)' ) im, imap(im), jmap(im)
         endif
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read physical data for all species                               c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      read   ( idat ) spdata
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read partition function data                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      read   ( idat ) pfdata
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     specify test volume and some useful constants                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      vol  = 1.0d0
      ceta = dsqrt( cpi ) / ( 4.0d0 * gs(ise) )
      cf4  = 2.0d0 * dsqrt( cpi ) * ce**3 / ( 3.0d0 * dsqrt( ck ) )
      cx   = 3.0d0 * cf4 / ck
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     set threshold occupation fraction                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      floor = 1.0d-14
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     read specification of temperature-density grid                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c ------ add Ladislav Zejda's new indexing of species with 
c ------ partition fcts in preparation for symbolic-algebra-generated
c ------ f2 .....
c ----------------------------------------------------------------------
      call lsetup
c ----------------------------------------------------------------------
c
c     close data file
      close(idat)
      
      return
      end
c
      subroutine thermo(isave)
c---------------------------------------------------------------------
c
c version that relegates the machine-dependent part (that is, the
c linear-equation solving with, or without, the CRAY 
c assembly-language s/r) to a new s/r (solvth).
c                                WD, 7/7/00
c---------------------------------------------------------------------
c
c...... isave=1: save 'frac' into 'fracx'. useful for starting values
c......          of isotherms beginning at densities too high for
c......          usual starting procedure.
c...... isave=0: do not save into 'fracx'.
c
c...... caution: thermo has to be used in a restricted way:
c...... i.e. through stepping on isotherms.
c
c***********************************************************************
c     evaluate all thermodynamic quantities                            c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      dimension d2fdlt(mlam), d2fdlv(mlam), d2fdl2(mlam, mlam)
      dimension dnedni(mspes), dnedl(mlam)
      dimension res2(mlam,2)
      equivalence (d2fdlt, d2zdnt), (d2fdlv, d2zdnv), (d2fdl2, d2zdn2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     gas pressure, internal energy, free energy, electron pressure    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     gas pressure
      pgas = sum(p)
c     ===================================================================
      if ( pgas .le. 0.d0 ) then
        write(6,*) 'pgas not positive: pgas,logt,logrho = ',
     .              pgas,tlog,rholog
        pgas = 1.d-40
      end if
c     ===================================================================
      pglog(irho) = log10( pgas )
c
c     internal energy (per cm**3)
      egas = sum(e) / vol
c     ===================================================================
      if ( egas .le. 0.d0 ) then
        write(6,*) 'egas not positive: egas,logt,logrho = ',
     .              egas,tlog,rholog
        egas = 1.d-40
      end if
c     ===================================================================
      elog(irho) = log10( egas )
c
c     free energy
      fgas = sum(f)
c
c     electron pressure
      pe = p(3)
      pelog(irho) = log10( pe )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     number of electrons, number of massive particles                 c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sne(irho) = sn(ise)
      snm(irho) = totn
c     account for hydrogen molecules
      if ( kz(1) .gt. 0 ) snm(irho) = totn - sn(ish2) - sn(ish2p)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save molecule/atom/ion occupation fractions                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     save old values from previous point on isotherm
      if ( irho .gt. 1 ) call scopy( mion*mz, frac, 1, fraco, 1 )
c
c     generate current values
      do kchem = 1, nchem
         is1 = ichm1(kchem)
         is2 = ichm2(kchem)
         kk  = nucz (kchem)
         sumsn = ssum( nion(kchem), sn(is1), 1 )
c     allow for hydrogen molecules
         if ( kk .eq. 1 ) sumsn = sumsn + sn(ish2) + sn(ish2p)
c
         do is = is1, is2
            jion = is - is1 + 1
            frac(jion, kk) = sn(is) / sumsn
         enddo
c
      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     save electron fraction, total number of nuclei, log temperature, c
c     log density                                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      frac(mion - 3, 1) = sn(ise) / totn
      frac(mion - 2, 1) = totn
      frac(mion - 1, 1) = tlog
      frac(mion    , 1) = rholog
c
c     at first point on isotherm force old values to be current values
      if ( irho .eq. 1 ) call scopy( mion*mz, frac, 1, fraco, 1 )
c
c     save values needed to initialize next isotherm
      if ( isave.eq.1 ) call scopy( mion*mz, frac, 1, fracx, 1 )
c
      if ( isave.eq.4 ) then
           call  scopy( nspes  , sn  , 1, snis(1,irho ), 1 )
           call  scopy( nspes  , wwt , 1, wwis(1,irho ), 1 )
c -------- avoid putting constant into argument if type is nonstandard
           ni1 = 1
           call iscopy( nchem  ,iwin1, ni1, iwis1(1,irho), ni1 )
           call iscopy( nchem  ,iwin2, ni1, iwis2(1,irho), ni1 )
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     thermodynamic quantities                                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c-----------------------------------------------------------------------
c     add radiation terms where needed                                 c
c-----------------------------------------------------------------------
c
      d2fdt2 = d2fdt2 - 4.00000000000000 d0 * carad * t*t * vol
      d2fdtv = d2fdtv - 1.33333333333333 d0 * carad * t*t*t
      ptot   = pgas   + 0.33333333333333 d0 * carad * t*t*t*t
      etot   = egas   + 1.00000000000000 d0 * carad * t*t*t*t * vol
      ftot   = fgas   - 0.33333333333333 d0 * carad * t*t*t*t * vol
      ptlog(irho) = log10( ptot )
      etlog(irho) = log10( etot )
c
c-----------------------------------------------------------------------
c     convert derivatives of f wrt populations to derivatives wrt      c
c     reaction parameters                                              c
c-----------------------------------------------------------------------
c
c     first zero storage
      fscr = 0.0d0
      d2fdlt = 0.0d0
      d2fdlv = 0.0d0
      d2fdl2 = 0.0d0
      dnedni = 0.0d0
      a = 0.0d0
      
c
c-----------------------------------------------------------------------
c     ignore d2f/dlam dt, etc for case of no reaction parameters       c
c-----------------------------------------------------------------------
c
      if ( nlam > 0 ) then
c
c     d2f/dlam dt = b * d2f/dn dt
         call smv( nlam, nspes, b, mlam, d2fdnt, d2fdlt )
c
c     d2f/dlam dv = b * d2f/dn dv
         call smv( nlam, nspes, b, mlam, d2fdnv, d2fdlv )
c
c     dne/dlam = b * dne/dn
         dnedni(ise) = 1.0d0
         call smv( nlam, nspes, b, mlam, dnedni, dnedl )
c
c     form d2f/dn2 * b(transpose); store temporarily in matrix a
         call smmtr(nspes,nspes,nlam,d2fdn2,mspes,b,mlam,a,mspes)
c
c     d2f/dlam2 = b * d2f/dn2 * b(transpose)
         call smm( nlam, nspes, nlam, b, mlam, a, mspes, d2fdl2, mlam)
c
c======================================================================
c     MACHINE-DEPENDENT LINEAR EQUATION SOLVER
c======================================================================
         mdim = mlam
         call solvth(d2fdlt,d2fdlv,d2fdl2,res2,mdim)
c======================================================================
c
c-----------------------------------------------------------------------
c     c sub v                                                          c
c-----------------------------------------------------------------------
c
         do ilam=1,nlam
            fscr(ilam)=res2(ilam,1)
         enddo
c
      endif
      
      csubv(irho) = -t * ( d2fdt2 - sdot(nlam, d2fdlt, 1, fscr, 1) )
     .                 / ( rho * vol )
c
c-----------------------------------------------------------------------
c     chi sub t                                                        c
c-----------------------------------------------------------------------
c
      chit(irho) = -t * ( d2fdtv - sdot(nlam, d2fdlv, 1, fscr, 1) )
     .                / ptot
c
c-----------------------------------------------------------------------
c     dne/dlnt
c-----------------------------------------------------------------------
c
      eltpun(irho) = -t * sdot(nlam, dnedl, 1, fscr, 1)
c
c-----------------------------------------------------------------------
c     dne/dln rho
c-----------------------------------------------------------------------
c
      elrpun(irho) = vol * sdot(nlam, dnedl, 1, fscr, 1)
c
c-----------------------------------------------------------------------
c     chi sub rho                                                      c
c-----------------------------------------------------------------------
c
      if ( nlam > 0 ) fscr(1:nlam) = res2(1:nlam,2)
      chirho(irho) = vol * ( d2fdv2 - sdot(nlam, d2fdlv, 1, fscr, 1) )
     .                   / ptot
c
c-----------------------------------------------------------------------
c     c sub p                                                          c
c-----------------------------------------------------------------------
c
      csubp(irho) = csubv(irho)
     .            + ptot * chit(irho)**2 / ( rho * t * chirho(irho) )
c
c-----------------------------------------------------------------------
c     q adiabatic                                                      c
c-----------------------------------------------------------------------
c
      qadb(irho) = chit(irho) / chirho(irho)
c
c-----------------------------------------------------------------------
c     gamma 3 (minus 1)                                                c
c-----------------------------------------------------------------------
c
      gam3(irho) = ptot * chit(irho) / ( rho * t * csubv(irho) )
c
c-----------------------------------------------------------------------
c     gamma 1
c-----------------------------------------------------------------------
c
      gam1(irho) = chit(irho) * gam3(irho) + chirho(irho)
c
c-----------------------------------------------------------------------
c     gamma 2,3
c-----------------------------------------------------------------------
c
      rat = gam1(irho) / gam3(irho)
      gam2(irho) = rat / ( rat - 1.0d0 )
      gam3(irho) = gam3(irho) + 1.0d0
c
c-----------------------------------------------------------------------
c     h,he ionization
c-----------------------------------------------------------------------
      frapun(1,irho) = frac(5,1)
      frapun(2,irho) = frac(2,2)
      frapun(3,irho) = frac(3,2)
      frapun(4,irho) = frac(1,1)
c
c-----------------------------------------------------------------------
c     eta, total number of nuclei and electrons
c-----------------------------------------------------------------------
c
      etapun(irho)   = eta
      tnpun (irho)   = totn
      elpun (irho)   = frac( mion - 3 , 1 )
c ----- additional variables for entropy (or cv,cp) to be put in varspc
      varspc(1,irho) = etot
      varspc(2,irho) = ftot

      varspc(3,irho) = sne(irho)
      varspc(4,irho) = snm(irho)

      varspc(5,irho) = (etot - ftot)/t ! entropy
    
      
cccc  cvspc = rho * csubv(irho) / ( ck * (snm(irho) + sne(irho)) )
cccc  cpspc = rho * csubp(irho) / ( ck * (snm(irho) + sne(irho)) )
      return
      end
c
      subroutine fmin(isave)
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c    the principal modifications with respect to the usual version are: 
c    1) possibility to transfer the flag 'isave' into s/r thermo,
c       and to announce (via 'isave'=-999) a failure of convergence to
c       the calling programme, 
c    2) in the computations of the cluster of 5 isotherms, the sn(.)'s
c       (and the window switches) are stored during the 1-st pass; these
c       values are recalled during the 4 subsequent passes. the gain in
c       computing time is tremendous, especially in difficult (i.e.
c       pressure ionization regions). the flags to set and read are
c       isave=4 and istart=4, respectively (s/r thermo and start),
c    3) checking for oscilliatory (bad) behaviour
c       of the h2-h-p-he-hep balance, which signals that we are way off
c       the equilibrium: a h2-h-p-he-hep scan is then performed (s/r scan),
c    4) input of a flag 'ifail' that allows, if set unequal to zero, to
c       skip the Newton-raphson iteration and to proceed directly to the
c       computation of the thermodynamic quantities, which will of course
c       be meaningless, but allow to continue long table computation jobs.
c       a clear warning is issued in such cases. typically 'ifail' is set
c       unequal to zero if several attempts (like reductions of steps in
c       density) fail to cause convergence,
c    5) periodically, the quality of convergence is checked. if it remains, 
c       or becomes anew, bad, premature exit (isave = -999) is forced. 
c       this helps to cut down time some of the time-consuming long iterations 
c       (iter --> niter).
c    6) finally, the number of occurences with ifail>0 is counted,
c       and the maximum number can be specified, beyond which the
c       program is aborted. this facility prevents (unplanned, but who
c       knows possible) infinite loops in the recovery attempts.
c
c***********************************************************************
c     free energy minimization                                         c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      dimension sn1st(mspes),wwt1st(mspes),iw11st(mchem),iw21st(mchem),
     .          snvar(5,50),smav(5),s2mav(5),sig(5)
c
c     ============= 1st dim: number of species checked: h2,h,p,he,hep.
c     ============= 2nd dim: max number of iterations.
c     ============= dims at least as much as nchk and niter.
c
      data  niter,  nchk,  nsctst, nsct1, nohop, dnohop, conv, convs
     .    / 50,     5,     12,     5,     7,     0.1d0,  1d-12, 1d-4/
      data  ifailm, nwarn, iqstp,  convq, facq , convw
     .    / 25,     25,    5,      1d-7 , 10.d0, 1d-10  /
c
      iter  = 0
      it1st = 1
      iscan = 0
      iquad = 0
      rth = 0.0d0
      rthe= 0.0d0
c ----------------------------------------------------------------------
c ---- quantities for testing (oscillatory) convergence of h-he species
c ----------------------------------------------------------------------
      snvar = 0.0d0
      kh  = kz(1)
      if ( kh  .gt. 0 ) rth  = 1.d0/(abun( kh)*totn)
c
      khe = kz(2)
      if ( khe .gt. 0 ) rthe = 1.d0/(abun(khe)*totn)
c
c ----------------------------------------------------------------------
c
      call  scopy(nspes,    sn, 1,  sn1st, 1)
      call  scopy(nspes,   wwt, 1, wwt1st, 1)
c -------- avoid putting constant into argument if type is nonstandard
      ni1 = 1
      call iscopy(nchem, iwin1, ni1, iw11st, ni1)
      call iscopy(nchem, iwin2, ni1, iw21st, ni1)
c
c--   write(6,*) 'fmin: t,r ',tlog,rholog
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate changes in the occupation numbers                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      dfmaxm = dnohop*1.0001d0
      dfmold = dfmaxm
c
 100  iter = iter + 1
c
      if(iter.ge.2) it1st = 0
c
      if(ifail.ne.0) goto 600
c
      call finddn
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     update occupation numbers                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call update
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     check for convergence                                            c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c===  write(6,*) iter,dfmax,dnmax,idfmax,idnmax
c
      if ( dfmax .le. conv   ) go to 600
c
      if (iquad.eq.0 .and. dfmax.lt.convq) iquad = iter
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     check for hopeless case                                          c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if (iter.eq.nohop .and. dfmaxm.gt.dnohop .or.
     .   (mod(iter,nohop).eq.0 .and. iter.gt.nsctst .and.
     .    dfmax.gt.dnohop .and. iscan.eq.0) .or.
     .   (iscan.eq.1 .and. iter.gt.2*nsctst .and. dfmax.gt.dnohop)) then
c
         isave = -999
         goto 1500
c
      end if
c
c ----------------------------------------------------------------------
c -------- test for bad (oscillatory) convergence of h-he species
c ----------------------------------------------------------------------
c
      if ( kh  .gt. 0 )  then
         snvar(1,iter) = sn(ish2)
         snvar(2,iter) = sn(ish )
         snvar(3,iter) = sn(ispr)
      end if
c
      if ( khe .gt. 0 )  then
         inhe          = ichm1(khe)
         inhe2p        = inhe + 2
         snvar(4,iter) = sn(inhe  )
         snvar(5,iter) = sn(inhe2p)
      end if
c
      if ( iter .ge. nsctst .and. iscan. eq. 0 ) then
c
         do kk=1,nchk
            smav (kk) = 0.d0
            s2mav(kk) = 0.d0
            do jj=1,nsct1
               ind       = iter   -   jj   +   1
               smav (kk) = smav (kk) + snvar(kk,ind)
               s2mav(kk) = s2mav(kk) + snvar(kk,ind)*snvar(kk,ind)
            enddo
         enddo
c
          facav   = 1.d0/dfloat(nsct1)
c
          do kk=1,nchk
             smav (kk) = facav*smav (kk)
             s2mav(kk) = facav*s2mav(kk)
             sig  (kk) = dsqrt(max(0.d0 , s2mav(kk)-smav(kk)*smav(kk)))
          enddo
c
          sig(1) = sig(1)*rth
          sig(2) = sig(2)*rth
          sig(3) = sig(3)*rthe
          sig(4) = sig(4)*rthe
          sig(5) = sig(5)*rthe
c
          sigth  = sig(1) + sig(2) + sig(3)
          sigthe = sig(4)
          sigthp = sig(5)
          sigtot = sigth  + sigthe + sigthp
          sigt01 = 0.01d0 * sigtot
          sigt99 = 0.99d0 * sigtot
c
          if ( sigtot .gt. convs ) then
c===          write(6,9007) tlog,rholog,sigtot,(sig(ll),ll=1,nchk)
              iscan = 1
              call  scopy(nspes, sn1st,1,   sn,1)
              call  scopy(nspes,wwt1st,1,  wwt,1)
c -------- avoid putting constant into argument if type is nonstandard
              ni1 = 1
              call iscopy(nchem,iw11st,ni1,iwin1,ni1)
              call iscopy(nchem,iw21st,ni1,iwin2,ni1)
              ihsc = 0
              if(sigth.gt.sigt99) ihsc  = 1
              ihesc = 0
              if(sigth.lt.sigt01) ihesc = 1
              if(ihesc.eq.1 .and. sigthe.lt.sigt01) ihesc=2
              dfmaxm = dfmold
              call fscan
          end if
      end if
c
c ----------------------------------------------------------------------
c -------------------- end test ----------------------------------------
c ----------------------------------------------------------------------
c
      iqua1 = iquad + iqstp
      iqua2 = iqua1 + iqstp
      if ( iqua2.ge.niter-2 ) iqua2 = niter - 2
      facqm = facq*dfmaxm
c
      if ( iquad.ne.0 .and. iter .ge. iqua1
     .                .and.(dfmax.lt.facqm .or. iter.ge.iqua2)
     .                .and. dfmax.lt. convq   .and.  iscan.eq.0) then
         if(dfmax.gt.convw)
     .      write(6,9011) iquad,iter,dfmax,dfmaxm,tlog,rholog
         goto 600
      end if
c
      if(iscan.ne.1 .and. dfmax.lt.dfmaxm) dfmaxm = dfmax
c
      if ( iter .lt. niter ) go to 100
c
      isave=-999
      goto 1000
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute thermodynamic quantitites                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
 600  call thermo(isave)
c
c  signal bad (though successful) convergence
c
 1000 if(iter.ge.nwarn .and. isave.ne.-999)
     .    write(6,9009) iter,dfmax,tlog,rholog
c
 1500 if ( ifail .ne. 0 ) then
          ifailt = ifailt + 1
          if(iloop.eq.1) write(6,9013) ifail,ifailt,tlog,rholog
      end if
c
      if (ifailt .ge. ifailm) then
          write(6,*) 'too many ifails. stop'
          stop 'ifail'
      end if
c
      return
c
 9009 format (' slow convergence: iter,dfmax,tlog,rholog: ',
     .         1x,i5,1p3e16.6)
 9011 format (' bad convergence: iquad,iter,dfmax,dfmaxm,tlog,rholog: ',
     .        /1x,2i5,1p4e16.6)
 9013 format (1x,76(1H*),/' all recovery attempts failed,',
     .      ' exit from fmin allowed: if this is a tabulated',
     .     /' point, the results will be completely wrong!',
     .      ' ifail,ifailt,tlog,rholog',
     .       /20x,2i5,1p2g16.6,/1x,76(1H*))
      end
c
      subroutine tdtab(icase)
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c ----- modifications of July 1990: the analagous recovery parts of the
c       pulling and isotherm calculations has been put into s/r addstp.
c -----
c
c...... icase=0: isotherms begin at sufficiently low density to permit
c...... a 'saha start' at irho=1. in this case, istart=0 for irho=1, and
c...... istart=1 for irho.gt.1 to be used in s/r start.
c
c...... icase=1: isotherms do not begin at such low densities. in this case
c...... a preparatory build up from low densities has to precede the first call
c...... of s/r 'fmin' in order to provide the values of 'fracx'. this buildup
c...... is done here. then istart=2 for irho=1, and istart=1 for irho.gt.1, and
c...... istart is to be used in s/r start.
c
c...... itab=1,2,.. choice from the tables in common/tab*/...
c
c
      include 'types'
      include 'parms'
      include 'tabparms'
      include 'coms'
c
      common/tab/tl2(nt2m),rl2(nt2m,nr2m),nt2,nr2
c
c------------------ consts. for pulling and nonconvergence recovery
      rpull0 = -15.
      drpull = 1.0
      drpx   = drho
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     principal loop over isotherms                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do it=1,nt
c     
         tlog0  = tl2(it)
         rhmin0 = rl2(it,1)
c     
c============== pull out from low density if needed (icase=1) ======
c               ****************************************************
c
         if(icase.eq.1) then
            npull = int((rhmin0 - rpull0)/drpull) + 1
            tlog = tlog0
            t    = 10.d0**tlog
c>>>>>>>>>>>>>>>>>>>>>>>>>>>>> pulling do loop <<<<<<<<<<<<<<<<<<<<<<<<
            do irho=1,npull
               rholog = rhmin0 - dfloat(npull-irho)*drpull
               rho    = 10.d0**rholog
               drpx   = drpull
c
               istart = 1
               if(irho.eq.1) istart=0
               istrtx = istart
               call start(istart)
c..................................... to save fractions into fracx
               isave = 0
               if(irho.eq.npull) isave = 1
               isavx = isave
c
c              minimize free energy
               ifail = 0
               call fmin(isave)
c
c------------- save last successful data
c
               if(isave.ne.-999) then
                   call scopy ( mion*mz, frac, 1, fraco,  1 )
                   call scopy ( nspes  , sn  , 1, snkeep, 1 )
                   call scopy ( nspes  , wwt , 1, wwkeep, 1 )
c -------- avoid putting constant into argument if type is nonstandard
                   ni1 = 1
                   call iscopy( nchem  ,iwin1, ni1, iwin1k, ni1 )
                   call iscopy( nchem  ,iwin2, ni1, iwin2k, ni1 )
c
c---------------------  measure taken if nonconvergence ----------------
c
               else
c
                   if(irho.eq.1) then
c         ------------  exit because of hopeless choice of rpull0 ------
                        write(6,9901) tlog,rholog
                        stop 'tdtab1'
                   else
                        call addstp(drpull,isavx,istrtx)
                   end if
c
                end if
c
c---------------------- end measure for nonconvergence -----------------
c
             enddo
c
c>>>>>>>>>>>>>>>>>>>>>>> end pulling do loop <<<<<<<<<<<<<<<<<<<<<<<<<<
c
             write(6,9004) tlog,rholog
c
          end if
c
c============== end pull out from low densities ====================
c               *******************************
c--------------------------------------------------------------------
c     loop over perturbed isotherms (to prepare numerical derivatives)
c---------------------------------------------------------------------
          do iloop = 1, 5
             if ( iloop .eq. 1 ) then
                tlog = tlog0
                rhomin = rhmin0 - ddr
             else if ( iloop .eq. 2 ) then
                tlog = tlog0
                rhomin = rhmin0
             else if ( iloop .eq. 3 ) then
                tlog = tlog0
                rhomin = rhmin0 + ddr
             else if ( iloop .eq. 4 ) then
                tlog = tlog0 - ddt
                rhomin = rhmin0
             else if ( iloop .eq. 5 ) then
                tlog = tlog0 + ddt
                rhomin = rhmin0
             end if
             t    = 10.0d0**tlog
c     
c++++++++++++++++++++++++
c     loop over densities
c++++++++++++++++++++++++
c     
             do irho = 1, nrho
                rholog = rhomin + drho * dfloat( irho - 1 )
                rho    = 10.0d0**rholog
c     
                if(iloop.eq.2) then
c     
c..........test for correct density construction
c     
                   rl10 = rl2(it,irho)
c     
                   deltr  = rholog - rl10
c     
                   if(abs(deltr).gt.1.d-8) then
                      print 9001,it,irho,rl10,rholog,deltr
                      stop'tdtab2'
                   end if
                end if
c     
                istart = 1
                if(icase.eq.0.and.irho.eq.1) istart=0
                if(icase.ne.0.and.irho.eq.1) istart=2
c     
                if( iloop.gt.1 ) istart=4
                istrtx = istart
c     
                call start(istart)
c     
c     minimize free energy

                isave=0
                if( iloop.eq.1 ) isave = 4
                isavx = isave
                ifail = 0
                call fmin(isave)
c     
c-------------save last successful data
c     
                if(isave.ne.-999) then
                   call  scopy( mion*mz, frac, 1, fraco,  1 )
                   call  scopy( nspes  , sn  , 1, snkeep, 1 )
                   call  scopy( nspes  , wwt , 1, wwkeep, 1 )
c     -------- avoid putting constant into argument if type is nonstandard
                   ni1 = 1
                   call iscopy( nchem  ,iwin1, ni1, iwin1k, ni1 )
                   call iscopy( nchem  ,iwin2, ni1, iwin2k, ni1 )
c     
c=====================measure taken if nonconvergence =============
c     
                else
c     
                   if (iloop.eq.1) then
                      if (irho .eq.1) then
c     --------  exit because of hopeless choice of end-pull point
                         write(6,9902) tlog,rholog
                         stop 'tdtab3'
                      else
                         call addstp(drho,isavx,istrtx)
                      end if
                   else
                      ifail = 1
                      call start(istart)
                      call fmin (isave )
                      write(6,9905) tlog,rholog
                   end if
c     
                end if
c     
c=====================end measure for nonconvergence ==============
c     
c     output
c$$$  call outp0
c$$$  call outp1
c     
             enddo              ! rho-loop
c++++++++++++++++++++++++++++
c     end loop over densities
c++++++++++++++++++++++++++++
c     
c     finish output
             call outp3
c     
          enddo
c--------------------------------------
c     end loop over perturbed isotherms
c--------------------------------------
c     
       enddo                    !T-loop
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     end principal loop over isotherms                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
      return
c
 9001 format(' error in s/r tdtab: wrong density on isotherm, ',
     . 'it,irho,r10,rholog,deltr = ',/1x,2i10,1p3g15.6)
 9004 format(' end pulling to intermediate rho: tlog,rholog= ',
     . 2f10.6)
 9901 format(' error in s/r tdtab: starting density to high, ',
     . 'in pulling: tlog,rholog = ',/1x,2f10.6)
 9902 format(' error in s/r tdtab: impossible starting point ',
     . 'from pulling: tlog,rholog = ',/1x,2f10.6)
 9905 format(1x,76(1H*),/' starting point from iloop=1 was unusable. ',
     . /' ifail was put to 1 to allow exit: tlog,rholog = ',
     .    2f10.6,/' the results will be wrong!',/1x,76(1H*))
      end
c
      subroutine addstp(drh,isvx,istx)
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c......  this s/r adds intermediate steps in tabulation
c......  on an isotherm, in case the previous solution does not
c......  lead to convergence of fmin. rholog is the density to reach,
c......  rholog - drh the last successful density.
c
      include 'types'
      include 'parms'
      include 'coms'
c
c......  choose ifsc at least .ge. 2. otherwise, no oscillating
c......  between 'normal' and scanning starts possible
c
      data ifsc,   epstp
     .    /   5,   1.d-4  /
c
c ......... addstp should only be called if iloop=1. therefore
c ......... istx must not be 4.
c
      if(istx.eq.4) then
         write(6,*) 'wrong parameters in addstp. >>>> stop <<<<'
         write(6,*) 'tlog,rholog,isvx,istx = ',tlog,rholog,isvx,istx
         stop 'addstp'
      end if
c
      rlprv  = rholog - drh
      rltrg  = rholog
c
      rhol2  = rlprv
      drfix  = drh
      drl    = drfix
      fac    = 0.5
c
      ifix   = 0
      lstsuc = ifix
      irl    = 0
c
 11   drl    = drl*fac
      drl0   = drl
c
      rhol1  = rhol2 + drl
      rholog = rhol1
      icut   = 0
      if(rhol1.ge.rltrg) then
         rhol1  = rltrg
         rholog = rltrg
         drl    = rltrg - rhol2
         icut   = 1
      end if
      rho    = 10.d0**rholog
      ifix   = ifix + 1
c
c............................... try scan procedure (each ifsc time)
c
      if(mod(ifix-lstsuc,ifsc).eq.0 .and. icut.eq.0) then
         isctr  = 1
         istart = 3
      else
         isctr  = 0
         istart = 5
      end if
c
      isave  = isvx
c
c ---------------- allow an erroneous exit if no recovery possible ----
      if(dabs(drl0).lt.epstp .and. isctr.ne.1 
     .                       .and. ifix.gt.lstsuc+ifsc) then
          ifail  = 1
          rholog = rltrg
          rho    = 10.d0**rholog
      else
          ifail  = 0
      end if
c ---------------------------------------------------------------------
c
      call start(istart)
      call fmin (isave )
c
      if(isave.eq.isvx) then
c
c----------------- save last successful composition (incl. window)
c
         if (ifail.eq.0) then
              call  scopy( mion*mz, frac, 1, fraco,  1 )
              call  scopy( nspes  , sn  , 1, snkeep, 1 )
              call  scopy( nspes  , wwt , 1, wwkeep, 1 )
c -------- avoid putting constant into argument if type is nonstandard
              ni1 = 1
              call iscopy( nchem  ,iwin1, ni1, iwin1k, ni1 )
              call iscopy( nchem  ,iwin2, ni1, iwin2k, ni1 )
         end if
c
         rhol2=rhol1
         fac=1.0
         if(irl.eq.1) fac=1.5
         irl=1
ccc---         write(6,*)' suc rec step:',
ccc---     .             ' ifix,lstsuc,isctr,irl,drl,tlog,rholog '
ccc---         write(6,*)  ifix,lstsuc,isctr,irl,drl,tlog,rholog
         lstsuc = ifix
c
         if(dabs(rholog-rltrg).lt.1.d-8) goto 20
         goto 11
      else
         irl=0
         fac = 0.7
         goto 11
      end if
c
 20   write(6,9019) tlog,rholog
c
      return
 9019 format(' successful recovery (addstp). tlog,rholog= ',2f10.6)
      end
c
      subroutine start(istart)
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c        modification from old s/r start: flag 'istart' to
c        control options. istart=0: compute saha equilibrium;
c        istart=1: use previous value; istart=2: use fracx to begin
c        at a point close to where it was saved in a previous call
c        of thermo; istart=3: call s/r scan independent of tlog and
c        rholog; istart=4: use sn(.)'s that were found when
c        iloop was 1; istart=5: use sn(.)'s of last converging
c        solution of fmin.
c
c    caution: these modifications do not allow 'free' use of fmin;
c             it is still needed to step on isotherms.
c........
      include 'types'
      include 'parms'
      include 'coms'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate total particle number (nuclei)                         c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      totn = rho * vol / ( gasmu * camu )
c
      if ( istart.eq.1 .or. istart.eq.3 ) goto 5
c
c     evaluate boltzmann factors
c
         ckt = ck * t
         do 2 ip = 1, npf
         do 1 il = 1, nint(nlev(ip))
         ekt(il, ip) = stwt(il, ip) * exp( - elev(il, ip) / ckt )
    1    continue
    2    continue
c
c     compute partition functions
c
         do 3 is = 1, nspes - 1
         z(is) = g0(is)
    3    continue
         do 4 ip = 1, npf
         is = ispes( jpf(ip), kpf(ip) )
         z(is) = ssum( nint(nlev(ip)), ekt(1, ip), 1 )
    4    continue
c
         if (istart.eq.0) then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        first point on isotherm: calculate occupation numbers from    c
c        saha equilibrium                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        saha equations: use saha1 if mixture contains hydrogen
c                        use saha2 otherwise
c
            kh = kz(1)
            if ( kh .gt. 0 ) then
               call saha1
            else
               call saha2
            end if
c
         else if( istart.eq.2) then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        use occupation fractions from previous isotherm at            c
c        current starting density                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            do 101 kchem = 1, nchem
            is1 = ichm1(kchem)
            is2 = ichm2(kchem)
            kk  = nucz (kchem)
            do 100 is = is1, is2
            jion = is - is1 + 1
            sn(is) = max( fracx(jion, kk)*abun(kchem)*totn, 1.d-70)
 100        continue
 101        continue
c
         else if( istart.eq.4 ) then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        use starting value from the same point of the first passage   c
c        through the isotherm cluster (iloop=1)                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            call  scopy( nspes, snis(1,irho ), 1  , sn  , 1 )
            call  scopy( nspes, wwis(1,irho ), 1  , wwt , 1 )
c -------- avoid putting constant into argument if type is nonstandard
            ni1 = 1
            call iscopy( nchem, iwis1(1,irho), ni1  ,iwin1, ni1 )
            call iscopy( nchem, iwis2(1,irho), ni1  ,iwin2, ni1 )
c
         else if( istart.eq.5 ) then
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        use starting value from the last coonvergent solution of fmin c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
            call  scopy(nspes,snkeep,1,   sn,1)
            call  scopy(nspes,wwkeep,1,  wwt,1)
c -------- avoid putting constant into argument if type is nonstandard
            ni1 = 1
            call iscopy(nchem,iwin1k,ni1,iwin1,ni1)
            call iscopy(nchem,iwin2k,ni1,iwin2,ni1)
c
      else
c
c........we should never be here, but if....
c........-------------------------------....
         write(6,*) 'error of parm in s/r start; istart,it,irho = ',
     .   istart,it,irho
         stop
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     in case we have used some stored starting values (istart=4,5)    c
c     we must here compute the reaction parameters anew, because       c
c     the windowing could have changed since the last successful       c
c     point (where the composition was stored). if this is negelcted,  c
c     one risks a crash of the next call of thermo. also reset         c
c     degeneracy-related quantities.                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if ( istart.eq.4 .or. istart.eq.5) then
c
           call lambda
           call new eta
           return
c
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c for istart=0,2 which do not fix electron number or specify window    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sn(ise) = max( sdot( nspes - 1, zs, 1, sn, 1 ), 1.0d-70 )
      call new eta
c
c     determine window containing species above threshold
c
      call window
c
      return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     all other cases                                                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
   5  if (istart.eq.3) go to 40
c
c=======================================================================
c     in normal temperature-density regime estimate occupation numbers c
c     from previous values of molecule/atom/ion occupation fractions   c
c=======================================================================
c
      do 30 kchem = 1, nchem
      is1 = ichm1(kchem)
      is2 = ichm2(kchem)
      kk  = nucz (kchem)
c-----------------------------------------------------------------------
c     hydrogen                                                         c
c-----------------------------------------------------------------------
c
      if ( kk .eq. 1 ) then
c     extrapolate if old value is above threshold
      do 10 jion = 1, nion(kchem)
      is = jion + is1 - 1
      xfrac(is) = cvmgt( frac(jion,kk)**2/fraco(jion,kk), 0.0d0,
     .                   fraco(jion,kk) .ge. floor )
   10 continue
c     otherwise rescale current value
      sum = 2.0d0*(wwt(ish2)*xfrac(ish2)+wwt(ish2p)*xfrac(ish2p))
     .           + wwt(ishm)*xfrac(ishm) + wwt(ish  )*xfrac(ish  )
     .           + wwt(ispr)*xfrac(ispr)
      do 11 jion = 1, nion(kchem)
      is = jion + is1 - 1
      xfrac(is) = cvmgt( frac(jion,kk)/sum, xfrac(is),
     .                   fraco(jion,kk) .lt. floor )
   11 continue
c
c
c-----------------------------------------------------------------------
c     all elements other than hydrogen                                 c
c-----------------------------------------------------------------------
c
      else
c     extrapolate if we have two previous values above threshold
      do 20 jion = 1, nion(kchem)
      is = jion + is1 - 1
      xfrac(is) = cvmgt( frac (jion, kk)**2 / fraco(jion, kk),
     .                   0.0d0,
     .                   frac (jion, kk) .ge. floor .and.
     .                   fraco(jion, kk) .ge. floor )
   20 continue
c     upward recursion if current value is above threshold
c     and previous value isn't
      do 21 jion = 2, nion(kchem)
      is = jion + is1 - 1
      xfrac(is) =
     .        cvmgt( xfrac(is-1) * frac(jion, kk) / frac(jion-1, kk),
     .        xfrac(is),
     .        frac (jion  , kk) .ge. floor .and.
     .        fraco(jion  , kk) .lt. floor .and.
     .        frac (jion-1, kk) .ge. floor )
   21 continue
c     downward recursion if current value is above threshold and
c     previous value isn't, and value not already set by upward
c     recursion
      do 22 jion = nion(kchem)-1, 1, -1
      is = jion + is1 - 1
      xfrac(is) =
     .           cvmgt( xfrac(is+1) * frac(jion, kk) / frac(jion+1, kk),
     .                  xfrac(is),
     .                  frac (jion  , kk) .ge. floor .and.
     .                  fraco(jion  , kk) .lt. floor .and.
     .                  frac (jion+1, kk) .ge. floor .and.
     .                  xfrac(is)         .le. 0.0d0 )
   22 continue
c     extrapolate if current value below threshold and previous
c     value isn't
      do 23 jion = 1, nion(kchem)
      is = jion + is1 - 1
      xfrac(is) = cvmgt( frac(jion, kk)**2 / fraco(jion, kk),
     .                   xfrac(is),
     .                   frac (jion, kk) .lt. floor .and.
     .                   fraco(jion, kk) .ge. floor )
   23 continue
c     rescale current value if both are below threshold
      sum = sdot( nion(kchem), wwt(is1), 1, xfrac(is1), 1 )
      do 24 jion = 1, nion(kchem)
      is = jion + is1 - 1
      xfrac(is) = cvmgt( frac(jion, kk) / sum,
     .                   xfrac(is),
     .                   frac (jion, kk) .lt. floor .and.
     .                   fraco(jion, kk) .lt. floor )
   24 continue
      end if
c
c-----------------------------------------------------------------------
c      now compute occupation numbers, renormalizing                   c
c      estimated occupation fractions to unity                         c
c-----------------------------------------------------------------------
c
      sum = sdot( nion(kchem), wwt(is1), 1, xfrac(is1), 1 )
c     account for hydrogen molecules
      if ( kk .eq. 1 ) sum = sum + wwt(ish2 ) * xfrac(ish2 )
     .                           + wwt(ish2p) * xfrac(ish2p)
      do 25 is = is1, is2
      xfrac(is) = xfrac(is) / sum
      sn(is) = max( wwt(is)*xfrac(is)*abun(kchem)*totn, 1.0d-70 )
   25 continue
c
c-----------------------------------------------------------------------
c     close loop over all chemical elements                            c
c-----------------------------------------------------------------------
c
   30 continue
c
c-----------------------------------------------------------------------
c     estimate electron number and degeneracy parameter                c
c-----------------------------------------------------------------------
c
      sn(ise) = max( sdot( nspes - 1, zs, 1, sn, 1 ), 1.0d-70 )
      call new eta
      go to 50
c
c=======================================================================
c     we arrive here if istart=3, which is triggered if the            c
c     stepping on the isotherm is going rough.                         c
c=======================================================================
c
  40  ihsc  = 0
      ihesc = 0
      call fscan
c
c=======================================================================
c     determine window containing species considered to be above       c
c     threshold; recompute electron number and degeneracy parameter    c
c     for consistency                                                  c
c=======================================================================
c
c     set window
c
   50 call window
c
c     update electron number and degeneracy parameter for consistency
c     with ions actually in window
c
      sn(ise) = max( sdot( nspes - 1, zs, 1, sn, 1 ), 1.0d-70 )
      call new eta
c
      return
c
      end
c
      subroutine update
c---------------------------------------------------------------------
c
c version of 30 Aug 1990 (ESTEC/CIRCE/CRAY-2)
c
c---------------------------------------------------------------------
c
c***********************************************************************
c     update occupation numbers                                        c
c
c .............. modification July 1990 ............................
c
c .............. fudge for very small electron fractions (flag ieldir):
c .............. ignore electrons resulting from iteration step,
c .............. do not include relative electron change in dfmax,
c .............. and compute new electron value directly.
c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      data   dmax  ,  ethrsh
     .    /  0.25d0,  1.d-5  /
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     compute and limit fractional change                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c -------------choose typical value around sqrt(floor)
c     
      if(it1st.eq.1) then
          ieldir = 0
          if(sn(ise)/totn .lt. ethrsh) ieldir=1
      end if
c     
      if(ieldir .eq. 1) then
          nscale = nspes - 1
          snelec = max( sdot( nscale, zs, 1, sn, 1 ), 1.0d-70 )
      else
          nscale = nspes
          snelec = -1.0d-70
      end if
c
c ........................................................................
c ........................................................................
c
c     compute change in occupation fractions from fractional change
c     in occupation numbers
      do kchem = 1, nchem
         is1 = ichm1(kchem)
         do jion = 1, nion(kchem)
            df(jion+is1-1) = dn(jion+is1-1) * sn(jion+is1-1)
     .           /(abun(kchem)*totn)
         enddo
      enddo
      df(ise) = dn(ise)
c
c     index of maximum change
      idfmax = isamax( nscale, df, 1 )
      dfmax  = dabs( df(idfmax) )
c
c     find scale factor and scale changes if necessary
      if ( dfmax .gt. dmax ) then
            sclfac = dmax / dfmax
            call sscal( nscale, sclfac, dn, 1 )
      end if
c
c     limit maximum fractional change
      do is = 1, nscale
         dn(is) = cvmgp( dn(is), -0.7d0, 0.7d0 + dn(is) )
         dn(is) = cvmgp( dn(is),  2.0d0, 2.0d0 - dn(is) )
      enddo
      idnmax = isamax( nscale, dn, 1 )
      dnmax  = dabs( dn(idfmax) )
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     apply the changes                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      do is = 1, nscale
         sn(is) = max( sn(is)*(1.0d0 + dn(is)), 1.0d-70 )
      enddo
c
      if (ieldir.eq.1) sn(ise) = snelec
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     update degeneracy parameter                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call new eta
c
      return
      end
c
      SUBROUTINE LSETUP
C-----------------------------------------------------------------------
C Purpose: This s/r provides the commons used in the Mupad-generated
C          s/r F2 of the MHD program.
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C Version: 0.0.0
C-----------------------------------------------------------------------

C-----------------------------------------------------------------------
C Authors:

C  Chia-Hsien Lin
C  chlin@physics.usc.edu
C  University of Southern California
C  1999

C  Ladislav Zejda
C  lada@centauri.usc.edu
C  University of Southern California
C  1999
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C References:
C (1) Hummer, D. G., Mihalas, D.: The Equation of State for Stellar  En-
C    velopes. I. An Occupation Probability Formalism for the  Truncation
C    of Internal Partition Functions; The  Astrophysical  Journal,  331:
C    794-814, 1988 August 15
C (2) Mihalas, D., Dappen, W., Hummer, D. G.: The Equation of State  for
C    Stellar Envelopes. II. Algorithm and selected results;  The  Astro-
C    physical Journal, 331: 815-825, 1988 August 15
C (3) Dappen, W., Mihalas, D., Hummer, D. G., Mihalas, B. W.: The  Equa-
C    tion of State for Stellar Envelopes. III. Thermodynamic Quantities;
C    The Astrophysical Journal, 332: 261-270, 1988 September 1
C-----------------------------------------------------------------------


C-----------------------------------------------------------------------
C Notation:
C NLEVSP(i)   = the number of energy levels for a species  i , for  bare
C              nuclei like proton this is equal to  0   even  they  have
C              ground state energy and CHISP factor for  level = 1 .
C CHISP (l,i) = the factor  16*((ZS(i)+1)^(1/2)*e^2/K(l,i)^(1/2)/chi(l,i
C              ))^3  for the energy level  l  of a species  i . See  the
C              equation  (8)  of the reference  (2)
C ELEVSP(l,i) = the  l-th energy level  of a species  i .
C RADSP (l,i) = the radius of a species  i  in the energy level  l .
C STWTSP(l,i) = the statistical weight of the energy  level   l   for  a
C              species  i .
C-----------------------------------------------------------------------

      IMPLICIT NONE
      
      include 'types'

      integer I,IP,IS,J,L

      include 'parms'
      include 'coms'

      DO I=1,NCHEM
         IS=ISPES(NION(I),I)
         NLEVSP(IS)=0
         CHISP(1,IS)=CHI0(IS)
         ELEVSP(1,IS)=E0(IS)
         RADSP(1,IS)=RAD0(IS)
         STWTSP(1,IS)=G0(IS)
         DO J=1,NION(I)-1
            IS=ISPES(J,I)
            IP=IPF(J,I)
            IF (IP.GT.0) THEN
               NLEVSP(IS)=NINT(NLEV(IP))
               DO L=1,NLEVSP(IS)
                  CHISP(L,IS)=CHI(L,IP)
                  ELEVSP(L,IS)=ELEV(L,IP)+E0(IS)
                  RADSP(L,IS)=RAD(L,IP)
                  STWTSP(L,IS)=STWT(L,IP)
               END DO
            ELSE
               NLEVSP(IS)=1
               CHISP(1,IS)=CHI0(IS)
               ELEVSP(1,IS)=E0(IS)
               RADSP(1,IS)=RAD0(IS)
               STWTSP(1,IS)=G0(IS)
            END IF
         END DO
      END DO
      
      RETURN
      END

      end module MHDeos
