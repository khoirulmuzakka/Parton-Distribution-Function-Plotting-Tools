      Subroutine SetVbpPar (Ihd1, Ihd2, Jbsn, Iset, Jord, Jschm, Jscal)
      
C		Sunday, April 1, 2001 : comment out Setpdf. This call should be made at the front-end
C		program module

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Logical LM, LL

      Common
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, IOrd, IBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / INTCM / AER, RER, ZRO, ONE

      Nbm  = Ihd1
      Ntg  = Ihd2
      IBsn = JBsn
      Lset = Iset
      IOrd = JOrd
C                       iSchm = 0/1  corresponds to msbar/DIS schemes
      schm = JSchm
C                       iScal = 0  : scale Amu = Q
C                             = 1  : scale Amu is treated as indep. variable.
      scal = JScal
C							Use ParVbp to set LL and LM as well.
      Call ParVbp(1, 'Ischm', schm, Ir)
      Call ParVbp(1, 'Iscal', scal, Ir)

C      Call Setpdf(Lset)

      Return
C                        ----------------------------
      Entry GetVbpPar (Ihd1, Ihd2, Jbsn, Iset, Jord, Jschm, Jscal)

      Ihd1 = Nbm
      Ihd2 = Ntg
      Jbsn = Ibsn
      Iset = Lset
      Jord = Iord
      Jschm= Ischm
      Jscal= Iscal

      Return
C                        ----------------------------
      Entry SetVbpVar (Rs, Qq, Scle, Rerr)

      Rts = Rs
      Q = Qq
C                              Set Amu only if Iscal = 0 (.Not. LL)
      If (.Not. LL) Amu = Scle
      Rer = Rerr

      Return
C                        -----------------------------
      Entry GetVbpVar (Rs, Qq, Scle, Rerr)

      Rs = Rts
      Qq = Q
      Scle = Amu 
      Rerr = Rer
      
      Return
C                        ****************************
      END

      Block Data DatVbp

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      Logical LL, LM

      COMMON
     1 / STCOM / RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / LOVBP / LM, LL, Ischm, Iscal
     1 / F0COM / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / INTCM / AER, RER, ZRO, ONE

       Data 
     1    RTS, Q, AMU, NBM, NTG, LSET, iOrd, iBsn
     1 / 1800.0, 10.0, 10.0, -1, 1, 1301, 2, 2 /
     1    AER, RER, ZRO, ONE,  LM, LL, Ischm, Iscal
     1 /  0D0, 2D-2, 0D0, 1D0, .True., .False., 0, 0 /

C                     *******************
      End
