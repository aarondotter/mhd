      subroutine outp0
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     write out basic results                                          c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
      if ( iout .le. 0 ) return
      if ( iloop .ne. 2 ) return
c
c     write header for first point on isotherm
c
      if ( irho .eq. 1 ) then
            write  ( iout, 1 ) tlog, hdr, jobdat, jobtim
    1       format ( '1',25x,'isotherm for log t =',f5.2 / 10a8/ 3a10 )
      end if
c
c     density, total atom number, electron number, degeneracy parameter
c
      write  ( iout, 2 )
    2 format (' ')
      write  ( iout, 3 ) rholog, totn, sn(ise), eta
    3 format ( ' log rho =',f7.2,' natom =',1pe12.5,' ne =',e12.5,
     .         ' eta =',0pf10.5 )
c
c     pressure
c
      write  ( iout, 4 ) ( jfe, p(jfe), jfe = 1, mfe ), pgas
    4 format ( 4(:' p',i1,' =',1pe11.4), ' pg =',e11.4 )
c
c     internal energy
c
      write  ( iout, 5 ) ( jfe, e(jfe), jfe = 1, mfe ), egas
    5 format ( 4(:' e',i1,' =',1pe11.4), ' eg =',e11.4 )
c
c     free energy
c
      write  ( iout, 6 ) ( jfe, f(jfe), jfe = 1, mfe ), fgas
    6 format ( 4(:' f',i1,' =',1pe11.4), ' fg =',e11.4 )
c
c     gas pressure, internal energy, and electron pressure in units of
c     k * t per particle
c
      ckt = ck * t
      pgas = pgas / ( ckt * (snm(irho) + sne(irho)) )
      egas = egas / ( ckt * (snm(irho) + sne(irho)) )
      pe   = pe   / ( ckt * sn(ise) )
      write  ( iout, 7 ) pgas, egas, pe
    7 format ( ' pg/nkt =',1pe11.4,'  e/nkt =',e11.4,'  pe/nekt =',
     .        e11.4 )
c
c     specific heats in units of k per particle;
c     chi sub t, chi sub rho
c
      cv = rho * csubv(irho) / ( ck * (snm(irho) + sne(irho)) )
      cp = rho * csubp(irho) / ( ck * (snm(irho) + sne(irho)) )
      write  ( iout, 8 ) cv, cp, chit(irho), chirho(irho)
    8 format ( '     cv =',1pe11.4,'     cp =',e11.4,'    chi t =',
     .         e11.4, ' chi rho =',e11.4 )
c
c     adiabatic gammas
c
      write  ( iout, 9 ) gam1(irho), gam2(irho), gam3(irho), qadb(irho)
    9 format ( ' gamma1 =',1pe11.4,' gamma2 =',e11.4,'   gamma3 =',
     .         e11.4, '    qadb =',e11.4 )
c
      return
      end
c
      subroutine outp1
      implicit double precision (a-h,o-z), integer (i-n)
c
c***********************************************************************
c     print atom/ion/molecule occupation fractions for h, he, c, n, o, c
c     ne, fe, and electrons                                            c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
      dimension abfrcs(mchem)
c-----------------------------------------------------------------------
c---- iprtz = 0: output fractions only, iprtz = 1: fractions and z -----
c-----------------------------------------------------------------------
      iprtz = 0
c-----------------------------------------------------------------------
      melec = mion - 3
c
      if ( ifracs .le. 0 ) return
      if ( iloop  .ne. 2 ) return
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     log t and log rho and composition                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      sum=0.d0
      do 9010 i=1,nchem
      sum = sum + abun(i)*atwt(i)
 9010 continue
      do 9020 i=1,nchem
      abfrcs(i) = abun(i)*atwt(i)/sum
 9020 continue
c     write(ifracs,9099) (atwt(ic),abun(ic)/abun(1),abfrcs(ic),
c    .                    ic=1,nchem)
c9099 format(/' Composition [atwt,abun,abfrcs]',/(1x,1p3g15.7))
c....................................................................
      write  ( ifracs, 1 ) tlog, rholog, abfrcs(1), abfrcs(2)
    1 format ( ' ' / ' log t = ',f12.8,' log rho = ',f12.8,
     .         '      X,Y = ',2f10.6,
     .         ' ' / 78('-'))
c....................................................................
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     hydrogen atom/ion/molecule occupation fractions                  c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kh = kz(1)
      if ( kh .gt. 0 ) then
         write  ( ifracs, 100 )
  100    format ( ' ' )
         write  ( ifracs, 2 ) ( frac(jion, 1), jion = 1, nion(kh) )
    2    format ( ' fh2    =',1pe11.4,'  fh2+  =',e11.4,
     .            '   fh-  =',  e11.4,'  fhyd  =',e11.4,
     .            ' ' / ' fprot  =',  e11.4 )
         is1 = ichm1(kh)
         is2 = ichm2(kh)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then
            write ( ifracs, 102 ) ( z(is), is = is1, is2 )
  102       format ( ' zh2 =',1pe19.12,' zh2+ =',e19.12,
     .               ' zh- =',  e19.12, 
     .         ' ' / ' zhyd =',e19.12,' zprot =',  e19.12 )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     helium atom/ion occupation fractions                             c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      khe = kz(2)
      if ( khe .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 3 )( jion, frac(jion, 2), jion = 1, nion(khe))
    3    format ( 4(: ' fhe(',i2,')=',1pe11.4) )
         is1 = ichm1(khe)
         is2 = ichm2(khe)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then
            write ( ifracs, 103 ) ( jspes(is), z(is), is = is1, is2 )
  103       format ( 3(: ' zhe(',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     carbon occupation fractions                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kc = kz(6)
      if ( kc .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 4 ) ( jion, frac(jion, 6), jion = 1, nion(kc))
    4    format ( 4(: ' fc (',i2,')=',1pe11.4) )
         is1 = ichm1(kc)
         is2 = ichm2(kc)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then 
            write ( ifracs, 104 ) ( jspes(is), z(is), is = is1, is2 )
  104       format ( 3(: ' zc (',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     nitrogen occupation fractions                                    c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kn = kz(7)
      if ( kn .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 5 ) ( jion, frac(jion, 7), jion = 1, nion(kn))
    5    format ( 4(: ' fn (',i2,')=',1pe11.4) )
         is1 = ichm1(kn)
         is2 = ichm2(kn)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then 
            write ( ifracs, 105 ) ( jspes(is), z(is), is = is1, is2 )
  105       format ( 3(: ' zn (',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     oxygen occupation fractions                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      ko = kz(8)
      if ( ko .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 6 ) ( jion, frac(jion, 8), jion = 1, nion(ko))
    6    format ( 4(: ' fo (',i2,')=',1pe11.4) )
         is1 = ichm1(ko)
         is2 = ichm2(ko)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then 
            write ( ifracs, 106 ) ( jspes(is), z(is), is = is1, is2 )
  106       format ( 3(: ' zo (',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     neon occupation fractions                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kne= kz(10)
      if ( kne .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 7 )( jion, frac(jion,10), jion = 1, nion(kne))
    7    format ( 4(: ' fne(',i2,')=',1pe11.4) )
         is1 = ichm1(kne)
         is2 = ichm2(kne)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then 
            write ( ifracs, 107 ) ( jspes(is), z(is), is = is1, is2 )
  107       format ( 3(: ' zne(',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     sodium occupation fractions                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kna = kz(11)
      if ( kna .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 8 )( jion, frac(jion,11), jion = 1, nion(kna))
    8    format ( 4(: ' fna(',i2,')=',1pe11.4) )
         is1 = ichm1(kna)
         is2 = ichm2(kna)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then 
            write  ( ifracs, 108 ) ( jspes(is), z(is), is = is1, is2 )
  108       format ( 3(: ' zna(',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     iron occupation fractions                                        c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kfe = kz(26)
      if ( kfe .gt. 0 ) then
         write  ( ifracs, 100 )
         write  ( ifracs, 9 )( jion, frac(jion,26), jion = 1, nion(kfe))
    9    format ( 4(: ' ffe(',i2,')=',1pe11.4) )
         is1 = ichm1(kfe)
         is2 = ichm2(kfe)
c---------------------------- optional output of z --------------------
         if(iprtz.eq.1) then 
            write  ( ifracs, 109 ) ( jspes(is), z(is), is = is1, is2 )
  109       format ( 3(: ' zfe(',i2,')=',1pe19.12) )
         end if
c----------------------------------------------------------------------
      end if
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     electron fraction                                                c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      write  ( ifracs, 100 )
      felec = frac(melec,1)
      write  ( ifracs, 10  ) felec
   10 format ( '  felec =',1pe11.4 )
c
      return
      end
c
      subroutine outp3
      implicit double precision (a-h,o-z), integer (i-n)
c ----- open and close statements require plain-vanilla integer variables
      integer ier
c
c***********************************************************************
c     write all thermodynamic quantities for complete isotherm.        c
c     unformatted.                                                     c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
      common/rp1/ tlg_array(5),rhomin_array(5),elg_array(mrho,5),
     & pglg_array(mrho,5),pelg_array(mrho,5),ptlg_array(mrho,5),
     & etlg_array(mrho,5),csbv_array(mrho,5),
     & csbp_array(mrho,5),chirh_array(mrho,5),cht_array(mrho,5),
     & gm1_array(mrho,5),gm2_array(mrho,5),
     & gm3_array(mrho,5),qdb_array(mrho,5),frp1_array(mrho,5),
     & frp2_array(mrho,5),frp3_array(mrho,5),
     & frp4_array(mrho,5),elpun_array(mrho,5),etapun_array(mrho,5),
     & eltpun_array(mrho,5),elrpun_array(mrho,5),
     & tnpun_array(mrho,5),eglg_array(mrho,5),sglg_array(mrho,5),
     & slg_array(mrho,5),iloopi_array(5),nrho_array(5)
c
      if(nchem.gt.mchem .or. nrho.gt.mrho) then
            write(6,1009) nchem,mchem,nrho,mrho
            stop
      end if
c
      ier = 0

c
      
!      if ( ipunch .le. 0 ) return

      
      if ( it .eq. 1 .and. iloop .eq. 1 ) then
c  don't need this when communicating output by common data
c            write ( ipunch, *, iostat = ier ) nchem,(atwt(ic),abun(ic),
c     .      ic=1,nchem),gasmu
c            write ( ipunch, *, iostat = ier ) nt, nrho, drho, ddt, ddr
      end if
      if ( ier .gt. 0 ) go to 2
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     write file                                                       c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      k = iloop
      tlg_array(k) = tlog
      rhomin_array(k) = rhomin
      iloopi_array(k) = iloop
      nrho_array(k) = nrho
      do i=1,nrho
         elg_array(i,k) = elog (i  )
         pglg_array(i,k) = pglog(i  )
         pelg_array(i,k) = pelog(i  )
         ptlg_array(i,k) = ptlog(i  )
         etlg_array(i,k) = etlog(i  )
         csbv_array(i,k) = csubv(i  )
         csbp_array(i,k) = csubp(i  )
         chirh_array(i,k) = chirho(i )
         cht_array(i,k) = chit (i  )
         gm1_array(i,k) = gam1 (i  )
         gm2_array(i,k) = gam2 (i  )
         gm3_array(i,k) = gam3 (i  )
         qdb_array(i,k) = qadb (i  )
         frp1_array(i,k) = frapun(1,i )
         frp2_array(i,k) = frapun(2,i )
         frp3_array(i,k) = frapun(3,i )
         frp4_array(i,k) = frapun(4,i )
         elpun_array(i,k) = elpun(i  )
         etapun_array(i,k) = etapun(i  )
         eltpun_array(i,k) = eltpun(i )
         elrpun_array(i,k) = elrpun(i  )
         tnpun_array(i,k) = tnpun (i )
         eglg_array(i,k) = eglog(i  )
         sglg_array(i,k) = sglog(i  )
         slg_array(i,k) = stlog(i  )
      end do

!      write( ipunch, *, iostat=ier ) nrho,iloop,tlog,rhomin,
!     .             (elog (i  ), i= 1, nrho), (pglog(i  ), i= 1, nrho),
!     .             (pelog(i  ), i= 1, nrho), (ptlog(i  ), i= 1, nrho),
!     .             (etlog(i  ), i= 1, nrho), (csubv(i  ), i= 1, nrho),
!     .             (csubp(i  ), i= 1, nrho), (chirho(i ), i= 1, nrho),
!     .             (chit (i  ), i= 1, nrho), (gam1 (i  ), i= 1, nrho),
!     .             (gam2 (i  ), i= 1, nrho), (gam3 (i  ), i= 1, nrho),
!     .             (qadb (i  ), i= 1, nrho),
!     .             ((frapun(ii,i ), i=1,nrho),ii=1,4),
!     .             (elpun(i  ), i= 1, nrho), (etapun(i  ),i= 1, nrho),
!     .             (eltpun(i ), i= 1, nrho), (elrpun(i  ),i= 1, nrho),
!     .             (tnpun (i ), i= 1, nrho),
!     .             ( (varspc(ivspc,i  ), ivspc=1,mvspc),  i= 1, nrho)
     
     
      if ( ier .le. 0 ) return
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     error return                                                     c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
    2 write  ( iout, 1 ) ier
    1 format ( ' error',i3,' on dataset punch.' )
      stop'outp3'
c
 1009 format(' error in s/r outp3: too small arrays: nchem,mchem,',
     .'nrho,mrho = ',/1x,4i10)
c
      end
