      subroutine solvfi
      implicit double precision (a-h,o-z), integer (i-n)
c======================================================================
c     solve with LEQ (contained in file cray_emul.f)
c======================================================================
c     special linear-equation solver for the s/r finddn.
c     the machine-dependent part is in this s/r, the idea being
c     that finddn becomes machine independent.
c                         WD, 7/7/00
c======================================================================
c     input/output transferred via common
c======================================================================
      include 'types'
      include 'parms'
      include 'coms'
c
      call leq(a,dn,nspes,1,mspes,mspes,err)
      if(err.ne.0.d0) goto 4
      write ( iout, 3 ) err
    3 format(' error return from leq in solvfi, err =',1pg12.5)
      stop'solvfi'
 4    continue
c
      return
      end
c
      subroutine solvth(d2fdlt,d2fdlv,d2fdl2,res2,mdim)
      implicit double precision (a-h,o-z), integer (i-n)
c======================================================================
c     solve with LEQ (contained in file cray_emul.f)
c======================================================================
c     special linear-equation solver for the s/r thermo.
c     the machine-dependent part is in this s/r, the idea being
c     that thermo becomes machine independent.
c                         WD, 7/7/00
c======================================================================
      include 'types'
      include 'parms'
      include 'coms'
c
      dimension d2fdlt(mdim), d2fdlv(mdim), d2fdl2(mdim, mdim)
      dimension res2(mdim,2)
c
      do 10 ilam=1,nlam
      res2(ilam,1)=d2fdlt(ilam)
      res2(ilam,2)=d2fdlv(ilam)
 10   continue
      call leq(d2fdl2,res2,nlam,2,mlam,mlam,err)
      if ( err .eq. 0.d0) then
c           error return
            write ( iout, 77)
   77       format(' error return from leq in solvth')
            stop'solvth'
      end if
c
      return
      end
