C
C  PROGRAM:    GFD_D3
C              Calculate the generalized Fermi-Dirac function and its
C              derivatives up to the 3rd order.
C
C  REFERENFE:
C      Gong Z.G., Zejda L., Dappen W., Aparicio J.M., Comp. Phys. Commun. 
C               in preparation, 2000.
C      Aparicio J.M., ApJS, 117, 627-632, 1998.
C
C  PARAMETERS
C      ON ENTRY
C           IE        - Integer
C                       Index for function to be evaluated
C                       IE = 0 -- F_k
C                       IE = 1 -- d F_k / d (eta)
C                       IE = 2 -- d F_k / d (beta)
C                       IE = 3 -- d^2 F_k / d (eta)^2
C                       IE = 4 -- d^2 F_k / d (eta) d (beta)
C                       IE = 5 -- d^2 F_k / d (beta)^2
C                       IE = 6 -- d^3 F_k / d (eta)^3
C                       IE = 7 -- d^3 F_k / d (eta)^2 d (beta)
C                       IE = 8 -- d^3 F_k / d (eta) d (beta)^2
C                       IE = 9 -- d^3 F_k / d (beta)^3
C
C           DK        - Double precision
C                       Index of the Fermi-Dirac function.
C
C           ETA       - Double precision
C                       Degeneracy parameter.
C                       ETA   = mu / k T
C
C           BETA      - Double precision
C                       Relativity parameter.
C                       BETA = k T / m C**2
C
C  NOTES:
C     20 points Gauss quadrature is used in this example.
C     Choice of other number of points of Gauss quadratures 
C     can be made by calling different function call 
C     DFERMI***.F with its corresponding Gauss quadratures 
C     DQLAG***.F and DQLEG***.F.
C
      FUNCTION R_FD(IE, DK, ETA, BETA)
      DOUBLE PRECISION R_FD,DFERMI200,DK,ETA,BETA
      INTEGER IE
C
      IF (IE.LT.0 .OR. IE.GE.10) THEN
        WRITE(6,*) 'ERROR: INDEX OUT OF RANGE!'
        STOP
      ELSE IF (BETA.LT.0.D0) THEN
        WRITE(6,*) 'ERROR: BETA MUST BE GREATER THAN OR EQUALS TO 0!'
        STOP
      ELSE
        R_FD=DFERMI200(IE,DK,ETA,BETA)
      END IF
C
      RETURN
      END


      FUNCTION DFERMI_1(IB, X, PAR, N)
C
C  To compute the generalized Fermi-Dirac or its derivatives integrand.
C
C  PARAMETERS
C      ON ENTRY
C              IB        - Integer
C                          Index for function to be evaluated
C
C              X         - Double precision
C                          Integration variable for the FD function.
C
C              PAR(1)    - Double precision
C                          Index of the Fermi-Dirac function.
C
C              PAR(2)    - Double precision
C                          Degeneracy parameter.
C
C              PAR(3)    - Double precision
C                          Relativity parameter.
C
C      ON RETURN
C              DFERMI_1  - Double precision
C                          Approximation to the GFD integrand.
C
      DOUBLE PRECISION DFERMI_1,DK,ETA,BETA,X,PAR(N),DXT
C
C***FIRST EXECUTABLE STATEMENT  DFERMI_1
      DK   =PAR(1)
      ETA  =PAR(2)
      BETA =PAR(3)
      DXT=DSQRT(1.D0+0.5D0*X*BETA)
      DFERMI_1=0.0D0
C
      IF (IB.EQ.0) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**DK*DXT/(exp(X-ETA)+1.D0)
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.1) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**DK*DXT/(exp(X-ETA)+2.D0+exp(ETA-X))
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.2) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**(DK+1.D0)/(exp(X-ETA)+1.D0)/4.D0/DXT
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.3) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**DK*DXT/(exp(X-ETA)+2.D0+exp(ETA-X))
     .             *((1.D0-exp(ETA-X))/(1.D0+exp(ETA-X)))
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.4) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**(DK+1.D0)/(exp(X-ETA)+2.D0+exp(ETA-X))
     .             /4.D0/DXT
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.5) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=-X**(DK+2.D0)/(exp(X-ETA)+1.D0)/16.D0/DXT**3
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.6) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**DK*DXT/(exp(X-ETA)+2.D0+exp(ETA-X))
     .             *((1.D0-4.D0*exp(ETA-X)+exp(2.D0*(ETA-X)))/
     .             (1.D0+exp(ETA-X))**2)
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.7) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**(DK+1.D0)/(exp(X-ETA)+2.D0+exp(ETA-X))
     .             *((1.D0-exp(ETA-X))/(1.D0+exp(ETA-X)))
     .             /4.0D0/DXT
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.8) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=-X**(DK+2.D0)/(exp(X-ETA)+2.D0+exp(ETA-X))
     .             /16.D0/DXT**3
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      ELSE IF (IB.EQ.9) THEN
        IF ((X-ETA).LT.7.0D2) THEN
          DFERMI_1=X**(DK+3.D0)/(exp(X-ETA)+1.D0)*3.D0/64.D0
     .             /DXT**5
        ELSE
          DFERMI_1=0.0D0
        ENDIF
      END IF
C       
      RETURN
      END


      FUNCTION DFERMI_2(IB, X, PAR, N)
C
C  To compute generalized Fermi-Dirac or its derivatives integrand.
C
C  PARAMETERS
C      ON ENTRY
C              IB        - Integer
C                          Index for function to be evaluated
C
C              X         - Double precision
C                          Integration variable for the FD function.
C
C              PAR(1)    - Double precision
C                          Index of the Fermi-Dirac function.
C
C              PAR(2)    - Double precision
C                          Degeneracy parameter.
C
C              PAR(3)    - Double precision
C                          Relativity parameter.
C
C       ON RETURN
C              DFERMI_2  - Double precision
C                          Approximation to the GFD integrand when the
C                          Z=DSQRT(X) variable change has been made.
C
      DOUBLE PRECISION DFERMI_2,DK,DXT,ETA,BETA,X,PAR(N)
C
C***FIRST EXECUTABLE STATEMENT  DFERMI_2
      DK   =PAR(1)
      ETA  =PAR(2)
      BETA =PAR(3)
      DXT=DSQRT(1.D0+0.5D0*X*X*BETA)
      DFERMI_2=0.D0
C
      IF (IB.EQ.0) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=2.D0*X**(2.D0*DK+1.0D0)*DXT
        ELSE
          DFERMI_2=2.D0*X**(2.D0*DK+1.0D0)*DXT/(exp(X*X-ETA)+1.D0)
        ENDIF
      ELSE IF (IB.EQ.1) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=0.D0
        ELSE
          DFERMI_2=2.D0*X**(2.D0*DK+1.0D0)*DXT/(exp(X*X-ETA)+2.D0+
     .             exp(ETA-X*X))
        ENDIF
      ELSE IF (IB.EQ.2) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=X**(2.D0*DK+3.D0)/2.D0/DXT
        ELSE
          DFERMI_2=X**(2.D0*DK+3.D0)/(exp(X*X-ETA)+1.D0)/2.D0/DXT
        ENDIF
      ELSE IF (IB.EQ.3) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=0.0D0
        ELSE
          DFERMI_2=2.D0*X**(2.D0*DK+1.0D0)*DXT/(exp(X*X-ETA)+2.D0+
     .             exp(ETA-X*X))*((exp(X*X-ETA)-1.D0)/(exp(X*X-ETA)
     .             +1.D0))
        ENDIF
      ELSE IF (IB.EQ.4) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=0.0D0
        ELSE
          DFERMI_2=X**(2.D0*DK+3.D0)/(exp(X*X-ETA)+2.D0+
     .             exp(ETA-X*X))/2.D0/DXT
        ENDIF
      ELSE IF (IB.EQ.5) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=-X**(2.D0*DK+5.D0)/8.D0/DXT**3
        ELSE
          DFERMI_2=-X**(2.D0*DK+5.D0)/(exp(X*X-ETA)+1.D0)/
     .             8.D0/DXT**3
        ENDIF
      ELSE IF (IB.EQ.6) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=0.0D0
        ELSE
          DFERMI_2=2.D0*X**(2.D0*DK+1)*DXT/(exp(X*X-ETA)+2.D0+
     .             exp(ETA-X*X))*((exp(X*X-ETA)-4.D0+exp
     .             (ETA-X*X))/(exp(X*X-ETA)+2.D0+exp(ETA-X*X)))
        ENDIF
      ELSE IF (IB.EQ.7) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=0.0D0
        ELSE
          DFERMI_2=X**(2.D0*DK+3.D0)/(exp(X*X-ETA)+2.D0+exp(ETA-X*X))
     .           *((exp(X*X-ETA)-1.D0)/(exp(X*X-ETA)+1.D0))/2.0D0/DXT
        ENDIF
      ELSE IF (IB.EQ.8) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=0.0D0
        ELSE
          DFERMI_2=-X**(2.D0*DK+5.D0)/(exp(X*X-ETA)+2.D0+
     .             exp(ETA-X*X))/8.D0/DXT**3
        ENDIF
      ELSE IF (IB.EQ.9) THEN
        IF ((X*X-ETA).LT.-7.0D2) THEN
          DFERMI_2=X**(2.D0*DK+7.D0)*3.D0/32.D0/DXT**5
        ELSE
          DFERMI_2=X**(2.D0*DK+7.D0)/(exp(X*X-ETA)+1.D0)*3.D0/32.D0
     .             /DXT**5
        ENDIF
      END IF
C
      RETURN
      END


