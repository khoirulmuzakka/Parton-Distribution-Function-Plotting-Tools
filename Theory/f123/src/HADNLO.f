      SUBROUTINE HADNLO(Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XHLR)
C-----------------------------------------------------------------------------
C      Computes F(LR0) HADRON Helicity Amps from CAOT paper
C
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension XH123(3), XHLR(-1:1)
      PARAMETER (PI=3.14159265359)

      XHLR( 1)=FHADHEL( 1,Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)
      XHLR( 0)=FHADHEL( 0,Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)
      XHLR(-1)=FHADHEL(-1,Isch,Ihad,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)

      RETURN
      END
