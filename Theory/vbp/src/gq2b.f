      FUNCTION GQ2B(TTB)
C                                                   -=-=- gq2b

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / F0COM / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / GQ2COM / TA, TB
     1 / E866xf/ ixfx
      TB = TTB
cjfo  E866 xf mods
      xffac=1.
      if(ixfx.eq.1)xffac=1./(xa+xb)
      TAU = XA * XB
      FB = VBP_PDF(LSET, NTG, IB, TB, AMU, IRET)
      GQ2B = 1./ (TB - XB) *
     >         (   GCSTOT(TA, TB, XA, XB, TAU) * FB
     >          -  (XA**2 + (XA-TA)**2) / TA**3 / 2 * FB0(IB)*xffac )
C
     >      + HCSTOT( TA, TB, XA, XB, TAU) * FB

      RETURN
C                                            ****************************
      ENTRY QG2B(TTA)
C
      TA = TTA
cjfo E866 xf mods
      xffac=1.
      if(ixfx.eq.1)xffac=1./(xa+xb)
      TAU = XA * XB
      FA = VBP_PDF( LSET, NBM, IA, TA, AMU, IRET)
      QG2B = 1./ (TA - XA) *
     >         (   GCSTOT(TB, TA, XB, XA, TAU) * FA
     >          -  (XB**2 + (XB-TB)**2) / TB**3 / 2 * FA0(IA)*xffac )
C
     >      + HCSTOT( TB, TA, XB, XA, TAU) * FA

      RETURN
C                        ****************************
      END
