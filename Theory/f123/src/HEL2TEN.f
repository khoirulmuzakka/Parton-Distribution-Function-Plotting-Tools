      SUBROUTINE HEL2TEN(XBJ,Q,HMASS,XHLR,XH123)
C-----------------------------------------------------------------------------
C      Computes F123 HADRON Helicity Amps from CAOT paper
C
C
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)
      data isave /0/
      save isave

cv==========================
cv set hevay mass thresholds
C      IF(isave.eq.0) THEN
C         isave=1
C         write(6,*) " PATCH FOR TESTING: HMASS=0 "
C      ENDIF
C         hmass=0.0d0
cv==========================


      XPLUS=XHLR( 1)
      XZERO=XHLR( 0)
      XMINU=XHLR(-1)

      RHO=1.0
      IF(HMASS.GT.0.0) THEN
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)
      ENDIF

      XH123(1)=(1./2.)*(XPLUS+XMINU)
      XH123(2)=(XBJ/RHO**2)*(XPLUS+XMINU+2.0*XZERO)
      XH123(3)=(1./RHO)*(-XPLUS+XMINU)

      RETURN
      END
