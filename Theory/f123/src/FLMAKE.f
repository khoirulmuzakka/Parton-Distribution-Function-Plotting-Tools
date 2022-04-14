      FUNCTION FLMAKE(XBJ,Q,HMASS,XH123)
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

      CALL TEN2HEL(XBJ,Q,HMASS,XHLR,XH123)
      FLMAKE= XHLR( 0)

      RETURN
      END
