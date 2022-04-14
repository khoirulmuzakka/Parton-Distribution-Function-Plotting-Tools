      FUNCTION FMAKER(XBJ,Q,HMASS,F123)
C-----------------------------------------------------------------------------
C      Computes R-RATIO given F123 HADRON Helicity Amps
C
C
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123(3)
      PARAMETER (PI=3.14159265359)

      FMAKER=0.0

      RHO=1.0
      IF(HMASS.GT.0.0) THEN
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)
      ENDIF
      RHO2=RHO*RHO


      RATIO=0.0
      IF( F123(1).NE.0.0 )
     >     RATIO=  F123(2) / ( 2.0*XBJ* F123(1) ) * RHO2-1.0

      FMAKER=RATIO

      RETURN
      END
