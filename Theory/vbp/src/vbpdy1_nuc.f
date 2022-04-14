c ... vbpdy1_nuc.f (IJS/JYY)
c ... 03.05.2007

      FUNCTION VBPDY1_NUC (Y)
C                                                   -=-=- vbpdy1_nuc

C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: VbpXsc
C Description: cross-section calculations and some auxilaries
C ListOfFiles: vbpdy1_nuc wlepasymk gq11 dist
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/3vbp/RCS/VbpXsc.f,v 1.2 98/08/16 21:21:07 wkt Exp $ 
C #Log:	VbpXsc.f,v $
c Revision 1.2  98/08/16  21:21:07  wkt
c Modifications induced by changes from EwkPac5 to EwkPac6.
c   SetEwk changed; Setup and Call of EW coupling constants changed.
c 
c Revision 1.1  98/08/05  23:53:48  wkt
c Initial revision
c 
C			      Input parameters are in the common block /STCOM_NUC/
C
C			      Output crossection is normalized to the simple
C			      parton formula in the tree approximation.
C
C			      1-loop corrections are normalized to Born term;
C			      Formulas are from Kubar et.al. (Nucl.Phys.B175)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
C
      PARAMETER (PI = 3.1415927, EULER = 0.57721566, ALPHEM = 1./137.)
      PARAMETER (NFL = 6, NBN = 4)

      LOGICAL LM, LL
      COMMON
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM_NUC / AER, RER, ZRO, ONE
     1 / E866xf/ ixfx

      integer VBP_NFLTOT

      Dimension za(1000), qq2a(1000)
      EXTERNAL QQ2_nuc,QQ3_nuc,QQ4A_nuc,GQ1_nuc,GQ2A_nuc,QG1_nuc,
     #         QG2A_nuc,FMSQ_nuc
c
      external EwCpl2An
C
      DATA SML / 1.E-6 /
C
C                      Ignore heavy quarks  -- good approx & avoid wasting cpu
c      nflt = 4
c  Not ignoring heavy quarks anymore.  Also, if nflt set to 4, then
c  three flavor fits are confused.
c
      NFLT = VBP_NFLTOT()

      TAU = Q / RTS
      XA  = EXP( Y) * TAU
      XB  = EXP(-Y) * TAU
cjfo  mods for E866 xf
      xffac=1.
      if(ixfx.eq.1) xffac=1./(xa+xb)
cjfo
C                         Leading Order and 1-loop QQ and GQ terms 
      SIGLO = 0.
      SG1QQ = 0.
      SG1GQ = 0.

      DO 10, I = - NFLT, NFLT

          FA0 (I) = VBP_PDF(LSET1, NBM, I, XA, AMU,  IRET)
          FB0 (I) = VBP_PDF(LSET2, NTG, I, XB, AMU,  IRET)

 10   CONTINUE
C                                           ...............Quark-Antiquark Part
      DO 20, IA = - NFLT, NFLT
          IF (IA.EQ.0) GOTO 20
          DO 25 IB = -NFLT, NFLT
             CPL2 = EwCpl2An (IA, IBSN, IB)
             IF (CPL2 .GT. SML) THEN
C                                                            Born Term
              PART1 =   FA0 (IA)
              PART2 =   FB0 (IB)
              PARTON  = PART1*PART2
cjfo
              TBORN  = CPL2 * PARTON*xffac
              SIGLO = SIGLO + TBORN
              IF (Iord .NE. 1) THEN
C                                                                  Q - Q Terms
C                                                   Term 1
                QQ1 = 1. + (5./3.) * PI**2
     >            - (3./2.) * LOG( XA * XB /(1.-XA)/(1.-XB) )
     >            + 2. * LOG( XA/(1.-XA) ) * LOG( XB/(1.-XB) )
                TERM1  = (2./3.) * CPL2 * QQ1 * PARTON
cjfo
                SG1QQ = SG1QQ + TERM1*xffac
         IF (LM .OR. LL) THEN
           TMP1 = ADZINT(FMSQ_nuc,ZRO,XA,AER,RER,ER1,IER1,1,1)
           TMP2 = ADZINT(FMSQ_nuc,ZRO,XB,AER,RER,ER2,IER2,1,1)      
           TERM1M = ( TMP1+TMP2 ) * CPL2 * PARTON 
           SG1QQ = SG1QQ + TERM1M
         ENDIF
C                                                                Term 2
        TMP = ADZINT(QQ2_nuc,XA,ONE,AER,RER,ER,IER,1,1)
C	     NN = Intusz(za, qq2a)
C         TMP = GausInt(QQ2_nuc,XA,ONE,AER,RER,ER,IER)
         TERM2 = (2./3.) * CPL2 * TMP * FB0(IB)
         SG1QQ = SG1QQ + TERM2
C                                                                Term 3
         TMP = ADZINT(QQ3_nuc,XB,ONE,AER,RER,ER,IER,1,1)
         TERM3 = (2./3.) * CPL2 * FA0(IA) * TMP
         SG1QQ = SG1QQ + TERM3
C                                                                Term 4
         TMP = ADZINT(QQ4A_nuc,XA,ONE,AER,RER,ER,IER,1,1)
         TERM4 = (4./3.) * CPL2 * TMP
         SG1QQ = SG1QQ + TERM4
      ENDIF
      ENDIF
 25   CONTINUE
 20   CONTINUE
C                                               ...............Gluon-Quark Part
C
      IF (Iord.NE. 1) THEN
        IA = 0
        DO 30 IB = - NFLT, NFLT
          IF (IB.NE.0) Then
             CPL2 = EwCpl2Cn (IB, IBSN)
             IF (CPL2 .GT. SML) Then
C                                                              Term 1
                TMP = ADZINT(GQ1_nuc,XA,ONE,AER,RER,ER,IER,2,1)
                TERM5 = (1./2.) * CPL2 * TMP * FB0(IB)
cjfo
                SG1GQ = SG1GQ + TERM5
C                                                               Term 2
                TMP = ADZINT(GQ2A_nuc,XA,ONE,AER,RER,ER,IER,2,1)
                TERM6  = (1./2.) * CPL2 * TMP
                SG1GQ = SG1GQ + TERM6
             End If
          End If
C
  30    CONTINUE
C                                            ..................Quark-Gluon Part
        IB = 0
        DO 40, IA = - NFLT, NFLT
          IF (IA.NE.0) Then
             CPL2 = EwCpl2Cn (IA, IBSN)
C                                                              Term 1
             TMP = ADZINT(QG1_nuc,XB,ONE,AER,RER,ER,IER,2,1)
             TERM7  = .5 * CPL2 * TMP * FA0(IA)
cjfo
             SG1GQ = SG1GQ + TERM7
C                                                                Term 2
             TMP = ADZINT(QG2A_nuc,XB,ONE,AER,RER,ER,IER,2,1)
             TERM8  = .5 * CPL2 * TMP
             SG1GQ = SG1GQ + TERM8
          End If
  40    CONTINUE
      END IF
      CONTINUE
C                                            ...........................Done
      
      SIGNLO = SIGLO + VBP_ALPI(AMU) * (SG1QQ + SG1GQ)

      IF (Iord .EQ. 1) THEN
        VbpDy1_nuc = SIGLO
      ELSE
        VbpDy1_nuc = SIGNLO
      ENDIF
      
      RETURN
C                        ****************************
      END
