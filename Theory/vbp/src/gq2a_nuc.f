c ... 03.05.2007 (IJS/JYY)
c ... gq2a_nuc.f

      FUNCTION GQ2A_nuc (TTA)
C                                                   -=-=- gq2a

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM_NUC / AER, RER, ZRO, ONE
     1 / GQ2COM_NUC / TA, TB
C
      EXTERNAL QG2B_NUC, GQ2B_NUC
C
      TA = TTA

      TMP = ADZ2NT( GQ2B_NUC, XB, ONE, AER, RER, ER, IER, 2,1)

      GQ2A_nuc = VBP_PDF( LSET1, NBM, IA, TA, AMU, IRET)* TMP

      RETURN
C                                            ****************************
      ENTRY QG2A_nuc(TTB)

      TB = TTB

      TMP = ADZ2NT( QG2B_NUC, XA, ONE, AER, RER, ER, IER, 2,1)
      
      QG2A_nuc = VBP_PDF( LSET2, NTG, IB, TB, AMU, IRET)*TMP

      RETURN
C                        ****************************
      END
