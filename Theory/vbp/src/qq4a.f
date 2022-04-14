      FUNCTION QQ4A(TTA)
C                                                   -=-=- qq4a

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / F0COM / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM / AER, RER, ZRO, ONE
     1 / QQ4COM / TA, FA

      EXTERNAL QQ4B
C
      TA = TTA
C
      FA  = VBP_PDF( LSET, NBM, IA, TA, AMU, IRET)
      
      QQ4A = ADZ2NT( QQ4B, XB, ONE, AER, RER, ER, IER, 2,1)

      RETURN
C                        ****************************
      END
