        Function XINTNLO(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADHEL PROGRAM
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       Dimension OmGNLO(-1:1)
       PARAMETER (PI=3.14159265359)
       COMMON /CXINTNLO/
     >   ETA,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,IRET,IHEL,Ihad
       COMMON /CXINTNLOTMP/ JSET

        JSET=ISET  !*** THIS IS A PATCH TO PASS ALQCD THE PDF SET: FIO 22SEP99

        IHADRON=1
        IGLUON=0

CCC ******* NOTE, WE USE ETA HERE INSTEAD OF XIMIN BECAUSE THE SLOW-RESCALING
CCC ******* MASS CORRECTION IS BUILT INTO THE MATRIX ELEMENTS AND NEED NOT BE
CCC ******* PUT IN BY HAND

        XHAT = ETA/XI
        CALL PARNLO (XHAT,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU,OmGNLO)
        FTEMP = OmGNLO(IHEL)
        PdfTMP  = PDF(iSet, Ihad, IGLUON, XI, XMU, iRet)
        TEMP = FTEMP * PDFTMP / XI
        XINTNLO= TEMP

        RETURN
        END
