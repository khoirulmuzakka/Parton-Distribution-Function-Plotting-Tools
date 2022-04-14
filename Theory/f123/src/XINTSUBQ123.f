      Function XINTSUBQ123(Xin)
C-----------------------------------------------------------------------------
C The splitting function for Q->Q
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      DATA SMALL /1.E-14/
        Common /CINTSUBQ/
     >   F1Min,XBJin, XMUIN, XiTHIN, iSetIN, IPARTONin,IhadIn

      F1M=F1Min
      XBJ=XBJin
      XMU=XMUin
      XiTH=XiTHin
      iSet=iSetin
      IPARTON=IPARTONin
      Ihad=IhadIN

      X=Xin
      IF(X.GE.1.0-SMALL) X=1.0-SMALL

       Z =X
      Z0=XiTH

C-----------------------------------------------------------------------------
C --- PDF's

      iHadr=Ihad          !*** 2/28/05 FIO Fixed
      QrkPdf = Pdf(iSet,iHadr,IPARTON,XiTH/X, XMU, iRet)
      QrkPdf0= Pdf(iSet,iHadr,IPARTON,XiTH  , XMU, iRet)

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS:


        XLOG1= LOG( (F1M/XMU)**2 )
C ***MODIFIED     \/
         delSUB= (+1.)*(  2. - 3./2. *XLOG1  - 2.*XLOG1*Log(1. - z0)
     >          - 2.*Log(1. - z0)
     >          - 2.*Log(1. - z0)**2
     >          )

         XMEASURE= (1.-XiTH)   !*** INTEGRATION MEASURE

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS:

        regSUB=
     > (
     >  (-XLOG1 -1.)*( (  (1.+ z**2)*QrkPdf/X
     >                   -(1.+1.**2)*QrkPdf0/1.0) * 1./(1.-z) )
     >  +      (-2.)*( (  (1.+ z**2)*QrkPdf/X
     >                  - (1.+1.**2)*QrkPdf0/1.0) *
     >                                       Log(1.-z)/(1.-z) )  !*** modified
     > )


      TOTAL = regSUB + delSUB/XMEASURE* QrkPdf0 /1.0

      TOTAL = (4./3.) * TOTAL   !*** COLOR FACTOR

      XINTSUBQ123 = TOTAL

      Return
      End
