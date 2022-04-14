c ... 03.05.2007 (IJS/JYY)
c ... gq2b_nuc.f

      FUNCTION GQ2B_nuc(TTB)
C                                                   -=-=- gq2b

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / GQ2COM_NUC / TA, TB
     1 / E866xf/ ixfx

      TB = TTB
cjfo  E866 xf mods
      xffac=1.
      if(ixfx.eq.1)xffac=1./(xa+xb)
      TAU = XA * XB
      FB = VBP_PDF(LSET2, NTG, IB, TB, AMU, IRET)
      GQ2B_nuc = 1./ (TB - XB) *
     >         (   GCSTOT(TA, TB, XA, XB, TAU) * FB
     >          -  (XA**2 + (XA-TA)**2) / TA**3 / 2 * FB0(IB)*xffac )
C
     >      + HCSTOT( TA, TB, XA, XB, TAU) * FB

      RETURN
C                                            ****************************
      ENTRY QG2B_nuc(TTA)
C
      TA = TTA
cjfo E866 xf mods
      xffac=1.
      if(ixfx.eq.1)xffac=1./(xa+xb)
      TAU = XA * XB
      FA = VBP_PDF( LSET1, NBM, IA, TA, AMU, IRET)
      QG2B_nuc = 1./ (TA - XA) *
     >         (   GCSTOT(TB, TA, XB, XA, TAU) * FA
     >          -  (XB**2 + (XB-TB)**2) / TB**3 / 2 * FA0(IA)*xffac )
C
     >      + HCSTOT( TB, TA, XB, XA, TAU) * FA

      RETURN
C                        ****************************
      END
