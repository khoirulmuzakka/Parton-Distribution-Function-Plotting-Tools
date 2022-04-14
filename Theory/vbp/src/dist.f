      Subroutine Dist (Ist, x, q, pd)
C                                                   -=-=- dist

      Implicit Double Precision(A-H, O-Z)

      dimension pd(8)

      pd (1) = vbp_pdf (ist, 1, 1, x, q, Ir)
      pd (2) = vbp_pdf (ist, 1, 2, x, q, Ir)
      pd (3) = vbp_pdf (ist, 1, 3, x, q, Ir)
      pd (4) = vbp_pdf (ist, 1, -1, x, q, Ir)
      pd (5) = vbp_pdf (ist, 1, -2, x, q, Ir)
      pd (6) = vbp_pdf (ist, 1, -3, x, q, Ir)
      pd (7) = vbp_pdf (ist, 1, 0, x, q, Ir)
      pd (8) = vbp_pdf (ist, 1, 4, x, q, Ir)

      return
      end
C                                                          =-=-= QqbarA
