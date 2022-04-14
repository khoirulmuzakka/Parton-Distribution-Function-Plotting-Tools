      FUNCTION QQ4B(TTB)
C                                                   -=-=- qq4b

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / F0COM / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM / AER, RER, ZRO, ONE
     1 / QQ4COM / TA, FA
     1 /E866xf/ ixfx

      TB = TTB

      TAU = XA*XB
      FB = VBP_PDF( LSET, NTG, IB, TB, AMU, IRET)
C
cjfo  E866 xf mods
      xffac=1.
      if(ixfx.eq.1) xffac=1./(xa+xb)
cjfo
C      II = II + 1
C	If (II.eq.1) Print *, 'Tau, Xa, Ta, Xb, Tb, Fa, Fb  in QQ4B:'
C      If (II.le.200) Print '(1pE12.5, 6E12.5)', Tau,Xa,Ta,Xb,Tb,Fa,Fb      
	QQ4B = 1/(TA-XA)/(TB-XB) *
     >         ( GASTOT( TA, TB, XA, XB, TAU) * FA * FB
     >          - ( (XA/TA)**2 + 1) / 2 * FA * FB0(IB)*xffac
     >          - ( (XB/TB)**2 + 1) / 2 * FA0(IA) * FB*xffac
     >          + FA0(IA) * FB0(IB)*xffac)
     >        + HASTOT( TA, TB, XA, XB, TAU) * FA * FB

C      If (II.le.200) Print *, 'QQ4B = ', QQ4B
      RETURN
C                        ****************************
      END
