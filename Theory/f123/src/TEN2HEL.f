      SUBROUTINE TEN2HEL(XBJ,Q,HMASS,XHLR,XH123)
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


      F1=XH123(1)
      F2=XH123(2)
      F3=XH123(3)

      RHO=1.0
      IF(HMASS.GT.0.0) THEN
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)
      ENDIF

      XPLUS= F1-(RHO/2.0)*F3
      XMINU= F1+(RHO/2.0)*F3
      XZERO=-F1+(RHO**2/(2.0*XBJ))*F2

      XHLR( 1)=XPLUS
      XHLR( 0)=XZERO
      XHLR(-1)=XMINU

      RETURN
      END
