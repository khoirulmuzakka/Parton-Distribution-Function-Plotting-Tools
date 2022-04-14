      SUBROUTINE PARNLO (X,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU,OmGNLO)
C-----------------------------------------------------------------------------
C      Computes Parton Helicity Amps from CAOT paper
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension OmGNLO(-1:1)
      PARAMETER (PI=3.14159265359)
C !*** THIS IS A PATCH TO PASS ALQCD THE PDF SET: FIO 22SEP99
      COMMON /CXINTNLOTMP/ ISET

      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

      OmGNLO( 1)=0.0
      OmGNLO( 0)=0.0
      OmGNLO(-1)=0.0

      ALPHAS= AlQCD(xMu,iset)

      Q2= Q**2
      F1M2 = F1M ** 2
      F2M2 = F2M ** 2

CCCC   THIS CHECK IS THE SAME AS THE Q CHECK BELOW
      XMAX=Q2/((F1M+F2M)**2+Q2)
      IF((X.GE.1.0).OR.(X.GE.XMAX)) RETURN

CCCC     HMASS=0  !*** TEMP PATCH: THIS IS PARTON LEVEL
CCCC     SMIN=((F1M+F2M+HMASS)**2-HMASS**2)
      SMIN=(F1M+F2M)**2
      QMIN=SQRT( SMIN*X/(1.0-X) )
      IF(Q.LE.QMIN) RETURN

C     Cm-Energy for the hard process for this SHAT=S
      S = Q2 * (1./X - 1.)
      IF(S.LE.SMIN) RETURN
      RS= SQRT(S)

      DEL = DELTA(S, F1M2, F2M2)
      TLOG = + Log(4*F1M2*s/(s+F1M2-F2M2+Del)**2)
      ULOG = + Log(4*F2M2*s/(s-F1M2+F2M2+Del)**2)

      XLAM =  ( Q2 + S)          /(2.0*RS)
      BET =  DEL                /(2.0*RS)
      E1  =  ( F1M2 - F2M2 + S) /(2.0*RS)
      E2  =  (-F1M2 + F2M2 + S) /(2.0*RS)
      EQ  =  (-Q2 + S)          /(2.0*RS)

      GSPLUS = -((1./2 + e1*(-1. + e1/XLAM)/XLAM)*TLOG) -
     >   (1./2 + e2*(-1. + e2/XLAM)/XLAM)*ULOG
     >   - 2*bet*eq**2/(XLAM**2*rs)

      GXPLUS = -(F1M*F2M*(4*bet*rs + (TLOG + ULOG)*(-F1M2 - F2M2 + s)))/
     >   (4*XLAM**2*s)

      GXPLUS = (2.0)* GXPLUS        !*** ERROR IN ORIGINAL FIO (4/19/93)

      GAPLUS = bet*(-F1M2 + F2M2)/(XLAM**2*rs) -
     >    TLOG*(1./2 + e1*(-1. + e1/XLAM)/XLAM +
     >    F1M2*(-F1M2 + F2M2)/(2*XLAM**2*s)) +
     >    ULOG*(1./2 + e2*(-1. + e2/XLAM)/XLAM -
     >    F2M2*(-F1M2 + F2M2)/(2*XLAM**2*s))

      GSZERO = eq*(-F1M2 + F2M2)*(TLOG - ULOG)/(XLAM**2*rs) +
     >   bet*((-F1M2 + F2M2)**2 +
     >   Q2*(-F1M2 - F2M2 + 2*Q2))/(XLAM**2*Q2*rs) +
     >   (TLOG + ULOG)*(-(eq*(F1M2 + F2M2 - (-F1M2 + F2M2)**2/Q2))/
     >   (2*XLAM**2*rs) + F1M2*F2M2/(XLAM**2*s) +
     >   (F1M2 + F2M2)*(-(-F1M2 + F2M2)**2 + Q2**2 - 2*XLAM**2*s)/
     >   (4*XLAM**2*Q2*s))

      GXZERO = bet*F1M*F2M/(XLAM**2*rs) + F1M*F2M*(TLOG + ULOG)*
     >      (1/Q2 + (-F1M2 - F2M2 + s)/(2*XLAM**2*s))/2

      GXZERO = (2.0)* GXZERO        !*** ERROR IN ORIGINAL FIO (4/19/93)

C-----------------------------------------------------------------------------
C                                                               The Amplitudes
C                                                               --------------
C                                                    The RR helicity amplitude
       ORRGR2 = GSPLUS + GAPLUS
       ORRGLR = GXPLUS
       ORRGL2 = GSPLUS - GAPLUS
C                                             The LONG-LONG helicity amplitude
       OZZGR2 = GSZERO
       OZZGLR = GXZERO
       OZZGL2 = GSZERO
C                                                    The LL helicity amplitude
       OLLGR2 = GSPLUS - GAPLUS
       OLLGLR = GXPLUS
       OLLGL2 = GSPLUS + GAPLUS
C-----------------------------------------------------------------------------
C                                                       Assemble the integrand
C                                                       ----------------------
c     XFAC=ALPHAS/PI/2.0
c     OmGNLO( 1)=XFAC*(ORRGR2*GRQ**2 +2.0*ORRGLR*GLQ*GRQ +ORRGL2*GLQ**2)
c     OmGNLO( 0)=XFAC*(OZZGR2*GRQ**2 +2.0*OZZGLR*GLQ*GRQ +OZZGL2*GLQ**2)
c     OmGNLO(-1)=XFAC*(OLLGR2*GRQ**2 +2.0*OLLGLR*GLQ*GRQ +OLLGL2*GLQ**2)
C-----------------------------------------------------------------------------
C                                                       Assemble the integrand
C                                                       ----------------------
      XFAC=ALPHAS/PI/2.0
      OmGNLO( 1)=XFAC*(ORRGR2*GRQ1*GRQ2+
     >   ORRGLR*(GLQ1*GRQ2+GLQ2*GRQ1) +ORRGL2*GLQ1*GLQ2)
      OmGNLO( 0)=XFAC*(OZZGR2*GRQ1*GRQ2+
     >   OZZGLR*(GLQ1*GRQ2+GLQ2*GRQ1) +OZZGL2*GLQ1*GLQ2)
      OmGNLO(-1)=XFAC*(OLLGR2*GRQ1*GRQ2+
     >   OLLGLR*(GLQ1*GRQ2+GLQ2*GRQ1) +OLLGL2*GLQ1*GLQ2)
C
      Return
      End
