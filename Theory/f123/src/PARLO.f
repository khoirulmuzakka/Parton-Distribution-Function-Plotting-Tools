      SUBROUTINE PARLO (Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
C-----------------------------------------------------------------------------
C      Computes LO Parton Helicity Amps from AOT paper
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension OmGLO(-1:1), OHelLO(-1:1,-1:1)
      PARAMETER (PI=3.14159265359)
      Data iLeft, iLong, iRight / -1, 0, 1 /
     >     igL2, igRgL, igR2 / -1, 0, 1 /

      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))
C-----------------------------------------------------------------------------
      F1m2    = F1m**2
      F2m2    = F2m**2
      Q2      = Q**2
      DEL   = DELTA(-Q2,F1M2,F2M2)
C-----------------------------------------------------------------------------
      OHelLO(iRight, igR2)  =   (Q2+F1m2+F2m2 + DEL)/DEL
      OHelLO(iRight, igRgL) =   -2*F1m*F2m/DEL
      OHelLO(iRight, igL2)  =   (Q2+F1m2+F2m2 - DEL)/DEL

      OHelLO(iLeft, igR2)  =  OHelLO(iRight, igL2)
      OHelLO(iLeft, igRgL) =  OHelLO(iRight, igRgL)
      OHelLO(iLeft, igL2)  =  OHelLO(iRight, igR2)

      OHelLO(iLong, igR2)  = ((F1m2+F2m2) + ((F1m2-F2m2)**2/Q2))/DEL
      OHelLO(iLong, igRgL) = +2*F1m*F2m/DEL
      OHelLO(iLong, igL2)  = OHelLO(iLong, igR2)
C-----------------------------------------------------------------------------
C      WR=     GRQ**2 *OHelLO(iRight, igR2)  +
C     >    2.0*GRQ*GLQ*OHelLO(iRight, igRgL) +
C     >        GLQ**2 *OHelLO(iRight, igL2)
C
C      WZ=     GRQ**2 *OHelLO(iLong,  igR2)  +
C     >    2.0*GRQ*GLQ*OHelLO(iLong,  igRgL) +
C     >        GLQ**2 *OHelLO(iLong,  igL2)
C
C      WL=     GRQ**2 *OHelLO(iLeft,  igR2)  +
C     >    2.0*GRQ*GLQ*OHelLO(iLeft,  igRgL) +
C     >        GLQ**2 *OHelLO(iLeft,  igL2)
C-----------------------------------------------------------------------------
c    MODIFY TO ALLOW FOR Z-GAMMA TERMS   5/7/07
C-----------------------------------------------------------------------------
      WR=     GRQ1*GRQ2         *OHelLO(iRight, igR2)  +
     >    (GRQ1*GLQ2+GRQ2*GLQ1) *OHelLO(iRight, igRgL) +
     >        GLQ1*GLQ2         *OHelLO(iRight, igL2)

      WZ=    GRQ1*GRQ2          *OHelLO(iLong,  igR2)  +
     >   (GRQ1*GLQ2+GRQ2*GLQ1)  *OHelLO(iLong,  igRgL) +
     >       GLQ1*GLQ2          *OHelLO(iLong,  igL2)

      WL=    GRQ1*GRQ2          *OHelLO(iLeft,  igR2)  +
     >   (GRQ1*GLQ2+GRQ2*GLQ1)  *OHelLO(iLeft,  igRgL) +
     >       GLQ1*GLQ2          *OHelLO(iLeft,  igL2)

      OmGLO( 1)=WR
      OmGLO( 0)=WZ
      OmGLO(-1)=WL

C-----------------------------------------------------------------------------
      Return
      End
