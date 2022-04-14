      Function XINTSUB(X)
C-----------------------------------------------------------------------------
C The splitting function for G->Q
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
        Common  / CINTSUB /  XMU, XiTH, iSet, IPARTON , Ihad

      iGluon=0
      GluPdf = Pdf(iSet,iHad,iGluon,XiTH/X, XMU, iRet)
      SpltFn = ( X**2 + (1.0-X)**2 )/2.0
C  ***  moved 2.0 into XINTSUB  FIO 12/7/95
      XINTSUB = SpltFn * GluPdf / X

      Return
      End
