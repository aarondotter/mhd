      subroutine finddn
      implicit double precision (a-h,o-z), integer (i-n)
c
c---------------------------------------------------------------------
c
c version that relegates the machine-dependent part (that is, the
c linear-equation solving with or without the CRAY 
c assembly-language s/r) to a new s/r (solvfi).
c                                WD, 7/7/00
c---------------------------------------------------------------------
c
c***********************************************************************
c     solve for changes in occupation numbers                          c
c***********************************************************************
      include 'types'
      include 'parms'
      include 'coms'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     initialization                                                   c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      d2fdt2 = 0.0d0
      d2fdtv = 0.0d0
      d2fdv2 = 0.0d0
c
      do is = 1, nspes
         dfdn  (is) = 0.0d0
         d2fdnt(is) = 0.0d0
         d2fdnv(is) = 0.0d0
         fscr  (is) = 0.0d0
         dn    (is) = 0.0d0

         do js = 1, nspes
            a     (is, js) = 0.0d0
            d2fdn2(is, js) = 0.0d0
         enddo

      enddo
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     calculate f, df/dn, d2f/dn2                                      c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     translational free energy
      call f1
c     internal free energy
      call f2
c     free energy of degenerate electron gas
      call f3
c     free energy of coulomb interactions
      call f4
c..      write(6,*) ' end calling the f1-4'
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     assemble the system                                              c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      call matgen
c
c======================================================================
c     MACHINE-DEPENDENT LINEAR EQUATION SOLVER
c======================================================================
      call solvfi
c======================================================================
c
      return
      end
