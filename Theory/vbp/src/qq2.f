      FUNCTION QQ2(TA)
C                                                   -=-=- qq2

C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: QqbarA
C Description: q - qbar annihilation modules
C ListOfFiles: qq2 qq4a qq4b gastot
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/3vbp/RCS/QqbarA.f,v 1.1 98/08/05 23:53:42 wkt Exp $ 
C #Log:	QqbarA.f,v $
c Revision 1.1  98/08/05  23:53:42  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / F0COM / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 /E866xf/ixfx
      QQ2 = 0.
      FA = VBP_PDF(LSET,  NBM, IA, TA, AMU, IRET)

cjfo E866 xf mods
      alfx=2.*xa/(ta+xa)
      xffac=1.
      if(ixfx.eq.1)then
         alfx=(xa+xb)/(ta+xb)
         xffac=1./(xa+xb)
      endif
      F1  = (TA**2 + XA**2) / TA**2
     >    * LOG( alfx *(1.-XB) / XB ) * FA
cjfo     >    * LOG( 2. * XA *(1.-XB) / XB / (TA + XA) ) * FA
      F1S = 2. * LOG( (1.-XB) / XB  )  * FA0(IA)

      QQ2 = QQ2 + ( F1 - F1S) / ( TA - XA )
C
      F2  = (3./2.) * FA
      F2S = (3./2.) * FA0 (IA)

      QQ2 = QQ2 + ( F2 - F2S) / (TA - XA)
C
      F3 = ( - 2./TA - 3 * XA / TA**2 ) * FA
      QQ2 = QQ2 + F3
cjfo
      qq2=qq2*xffac
C                            Add terms due to Rnorm scheme & scale corrections
      IF (LM .OR. LL)
     1   QQ2 = QQ2 - (3./2.) * (1./TA) * FMSQ(XA/TA)
     >       * (FA - FA0(IA) * XA / TA)
cjfo  xffac adjusted as needed in fmsq
C
      RETURN
C                                       ****************************
      ENTRY QQ3(TB)

      QQ3 = 0.
      FB = VBP_PDF(LSET, NTG, IB, TB, AMU,  IRET)
      
cjfo E866 xf mods
      alfx=2.*xb/(tb+xb)
      xffac=1.
      if(ixfx.eq.1)then
         alfx=(xa+xb)/(tb+xa)
         xffac=1./(xa+xb)
      endif
      F1  = (TB**2 + XB**2) / TB**2
     >     * LOG( alfx *(1.-XA) / XA ) * FB
cjfo     >     * LOG( 2. * XB *(1.-XA) / XA / (TB + XB) ) * FB
      F1S = 2.
     >     * LOG( (1.-XA) / XA ) * FB0(IB)
      QQ3 = QQ3 + ( F1 - F1S) / ( TB - XB )
C
      F2  = (3./2.) * FB
      F2S = (3./2.) * FB0(IB)
      QQ3 = QQ3 + ( F2 - F2S) / (TB - XB)
C
      F3 = ( - 2./TB - 3 * XB / TB**2 ) * FB
      QQ3 = QQ3 + F3
cjfo
      qq3=qq3*xffac

C                   Factor 3/2 to compensate relative normalization  with above
      IF (LM .OR. LL) 
     1   QQ3 = QQ3 - (3./2.) * (1./TB) * FMSQ(XB/TB) 
     >       * (FB - FB0(IB) * XB/TB)
cjfo  xffac adusted as needed in fmsq
      RETURN
C                        ****************************
      END
