c ... 03.05.2007 (IJS/JYY)
c ... gq1_nuc.f

      FUNCTION GQ1_nuc(TA)
C                                                   -=-=- gq1_nuc

C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: GluQrk
C Description: gluon - quark (Compton) scattering modules
C ListOfFiles: gq1 gq2a gq2b gcstot
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/3vbp/RCS/GluQrk.f,v 1.1 98/08/05 23:53:33 wkt Exp $ 
C #Log:	GluQrk.f,v $
c Revision 1.1  98/08/05  23:53:33  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 /E866xf/ ixfx
cjfo
C
      IF (IA .NE. 0) THEN
        PRINT *, 'IA .NE. 0 in GQ1_nuc; IA =', IA
        call PAUSE('in gq1_nuc.f') 
      ENDIF
cjfo E866 xf mods
      alfx=2.*xa/(ta+xa)
      xffac=1.
      if(ixfx.eq.1)then
         alfx=(xa+xb)/(ta+xb)
         xffac=1./(xa+xb)
      endif
C                                                   gluon distr.
      FA = VBP_PDF(LSET1, NBM, IA, TA, AMU, IRET)
      FACTOR = ( XA**2 + (TA - XA)**2 ) / 2. / TA**3
     >         * LOG( alfx * (1. - XB ) / XB )
     >      + 1./ 2. / TA
     >      - 3. * XA * ( TA - XA ) / TA**3
      factor=factor*xffac
cjfo     >         * LOG( 2.* XA * (1. - XB ) / XB / ( TA + XA ) )
cjfo     >      + 1./ 2. / TA
cjfo     >      - 3. * XA * ( TA - XA ) / TA**3
      GQ1_nuc = FACTOR * FA

C                     Factor 2 to compensate relative normalization  with above
cjfo xffac adjusted in fmsg as needed
      IF (LM .OR. LL) GQ1_nuc = GQ1_nuc - 2. * (1/TA) * 
     #                          FMSG_nuc(XA/TA) * FA
C
      RETURN
C                                            ****************************
      ENTRY QG1_nuc(TB)

      IF (IB .NE. 0) THEN
        PRINT *, 'IB .NE. 0 in QG1_nuc; IB =', IB
        call PAUSE('in gq1_nuc.f') 
      ENDIF
cjfo E866 xf mods
      alfx=2.*xb/(tb+xb)
      xffac=1.
      if(ixfx.eq.1)then
         alfx=(xa+xb)/(tb+xa)
         xffac=1./(xa+xb)
      endif

      FB = VBP_PDF(LSET2, NTG, IB, TB, AMU, IRET)
      FACTOR = ( XB**2 + (TB - XB)**2 ) / 2. / TB**3
     >         * LOG( alfx * (1. - XA ) / XA )
     >      + 1./ 2. / TB
     >      - 3. * XB * ( TB - XB ) / TB**3
      factor=factor*xffac
cjfo     >         * LOG( 2.* XB * (1. - XA ) / XA / ( TB + XB ) )
cjfo     >      + 1./ 2. / TB
cjfo     >      - 3. * XB * ( TB - XB ) / TB**3
      QG1_nuc = FACTOR * FB

C                     Factor 2 to compensate relative normalization  with above
cjfo  xffac adjusted in fmsg as needed
      IF (LM .OR. LL) QG1_nuc = QG1_nuc - 2. * (1/TB) * 
     #                          FMSG_nuc(XB/TB) * FB
C
      RETURN
C                        ****************************
      END
