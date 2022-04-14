c ... 03.05.2007 (IJS/JYY)
c ... qq2_nuc.f

      FUNCTION QQ2_nuc(TA)
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
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 /E866xf/ixfx
      QQ2_nuc = 0.
      FA = VBP_PDF(LSET1,  NBM, IA, TA, AMU, IRET)

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

      QQ2_nuc = QQ2_nuc + ( F1 - F1S) / ( TA - XA )
C
      F2  = (3./2.) * FA
      F2S = (3./2.) * FA0 (IA)

      QQ2_nuc = QQ2_nuc + ( F2 - F2S) / (TA - XA)
C
      F3 = ( - 2./TA - 3 * XA / TA**2 ) * FA
      QQ2_nuc = QQ2_nuc + F3
cjfo
      qq2_nuc=qq2_nuc*xffac
C                            Add terms due to Rnorm scheme & scale corrections
      IF (LM .OR. LL)
     1   QQ2_nuc = QQ2_nuc - (3./2.) * (1./TA) * FMSQ_nuc(XA/TA)
     >       * (FA - FA0(IA) * XA / TA)
cjfo  xffac adjusted as needed in fmsq
C
      RETURN
C                                       ****************************
      ENTRY QQ3_nuc(TB)

      QQ3_nuc = 0.
      FB = VBP_PDF(LSET2, NTG, IB, TB, AMU,  IRET)
      
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
      QQ3_nuc = QQ3_nuc + ( F1 - F1S) / ( TB - XB )
C
      F2  = (3./2.) * FB
      F2S = (3./2.) * FB0(IB)
      QQ3_nuc = QQ3_nuc + ( F2 - F2S) / (TB - XB)
C
      F3 = ( - 2./TB - 3 * XB / TB**2 ) * FB
      QQ3_nuc = QQ3_nuc + F3
cjfo
      qq3_nuc=qq3_nuc*xffac

C                   Factor 3/2 to compensate relative normalization  with above
      IF (LM .OR. LL) 
     1   QQ3_nuc = QQ3_nuc - (3./2.) * (1./TB) * FMSQ_nuc(XB/TB) 
     >       * (FB - FB0(IB) * XB/TB)
cjfo  xffac adusted as needed in fmsq
      RETURN
C                        ****************************
      END
