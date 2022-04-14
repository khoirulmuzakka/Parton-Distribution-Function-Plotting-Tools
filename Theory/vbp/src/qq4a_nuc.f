c ... 03.05.2007 (IJS/JYY)
c ... qq4a_nuc.f

      FUNCTION QQ4A_nuc(TTA)
C                                                   -=-=- qq4a

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      LOGICAL LM, LL

      COMMON
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM_NUC / AER, RER, ZRO, ONE
     1 / QQ4COM_NUC / TA, FA

      EXTERNAL QQ4B_nuc
C
      TA = TTA
C
      FA  = VBP_PDF( LSET1, NBM, IA, TA, AMU, IRET)
      
      QQ4A_nuc = ADZ2NT( QQ4B_nuc, XB, ONE, AER, RER, ER, IER, 2,1)

      RETURN
C                        ****************************
      END
