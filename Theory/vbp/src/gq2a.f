      FUNCTION GQ2A (TTA)
C                                                   -=-=- gq2a

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / F0COM / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM / AER, RER, ZRO, ONE
     1 / GQ2COM / TA, TB
C
      EXTERNAL QG2B, GQ2B
C
      TA = TTA

      TMP = ADZ2NT( GQ2B, XB, ONE, AER, RER, ER, IER, 2,1)

      GQ2A = VBP_PDF( LSET, NBM, IA, TA, AMU, IRET)* TMP

      RETURN
C                                            ****************************
      ENTRY QG2A(TTB)

      TB = TTB

      TMP = ADZ2NT( QG2B, XA, ONE, AER, RER, ER, IER, 2,1)
      
      QG2A = VBP_PDF( LSET, NTG, IB, TB, AMU, IRET)*TMP

      RETURN
C                        ****************************
      END
