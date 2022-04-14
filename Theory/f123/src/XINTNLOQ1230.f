        Function XINTNLOQ1230(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADNLOQ1230 PROGRAM
C
C      25 JUNE 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C
C
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DOUBLE COMPLEX  DiLog
       DOUBLE COMPLEX CARG,CSPENCE
       PARAMETER (PI=3.141592653589793)
       DATA SMALL /1.E-14/
       COMMON /CXINTNLOQ1230/
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad

        XINTNLOQ1230=0.0

        IHADRON=Ihad           !*** 2/28/05 FIO Fixed
        IGLUON=0

        XHAT = XIMIN/Xi
        XHAT0= XIMIN

        Z = XI
        z0= XIMIN

        IF(Z.GE.1.0-SMALL) Z=(1.0-SMALL)

        ALPHAS= AlQCD(XMU,iset)

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS:

        CARG= DCMPLX( 1.-z0 )
        CSPENCE= DILOG(CARG)
        SPENCE= DBLE(CSPENCE)

         FDEL2=(-1.0)*(
     >    PI**2/6. + 9./2.*Z0 + 5./4.*Z0**2 + 3.*LOG(1.-Z0)
     >   -(1./2.)* Z0 * (2.+Z0)*LOG(1.-Z0) - (LOG(1.-Z0))**2 - SPENCE
     >    )

         FDEL2=(-1.0)*(
     >   PI*PI/3. + 7.*z0/2. + z0**2 + 3.*Log(1.-z0) - z0*Log(1.-z0) -
     >   z0**2*Log(1. - z0)/2. - ( Log(1. - z0) )**2 + z0*Log(z0) +
     >   z0**2*Log(z0)/2. - 2.*SPENCE
     >   )

         FDEL1= (FDEL2 - 0.0)/2.
         FDEL3= (FDEL2 - 0.0)

         XMEASURE= (1.-z0)   !*** INTEGRATION MEASURE

C-----------------------------------------------------------------------------
C --- REGULAR TERMS:

         FREG2  = 0.0
         FREG21 = (2.0*Z)
         FREG23 = (1.0+Z)

         FREG1= (FREG2 - FREG21)/2.
         FREG3= (FREG2 - FREG23)

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS:

         FPLS2X= ((1.+ Z**2)/(1.-Z) * ( LOG((1.-Z)/Z) -3./4. ) +
     >   (9.+5.*Z)/4. )

         Pdfz  = PDF(iSet, IHADRON, IPART, XHAT , XMU, iRet)
         Pdfz0 = PDF(iSet, IHADRON, IPART, XHAT0, XMU, iRet)

         FPLS2  = FPLS2X*(Pdfz/Z  - Pdfz0/1.)

         FPLS21 = 0.0
         FPLS23 = 0.0

         FPLS1= (FPLS2 - FPLS21)/2.
         FPLS3= (FPLS2 - FPLS23)
C-----------------------------------------------------------------------------
C --- SUM TERMS:   (Standard F2 is defined with an overall XBJ factor)

          FTEMP1= FDEL1/XMEASURE*Pdfz0 + FREG1*Pdfz/Z + FPLS1
          FTEMP2=(FDEL2/XMEASURE*Pdfz0 + FREG2*Pdfz/Z + FPLS2) * XBJ
          FTEMP3= FDEL3/XMEASURE*Pdfz0 + FREG3*Pdfz/Z + FPLS3

        IF    (INDEX.EQ. 1)  THEN
          FTEMP=  FTEMP1 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
        ELSEIF(INDEX.EQ. 2)  THEN
          FTEMP=  FTEMP2 * ( GLQ1*GLQ2 + GRQ1*GRQ2 )
        ELSEIF(INDEX.EQ. 3)  THEN
          FTEMP=  FTEMP3 * ( GLQ1*GLQ2 - GRQ1*GRQ2 )
        ELSE
          WRITE(6,*) ' BAD INDEX ',INDEX
          STOP
        ENDIF

        FTEMP= FTEMP * ALPHAS/(2.0*PI) * (4./3.)
        FTEMP= FTEMP *2.0   !*** TO MATCH OUR NORMALIZATION

        XINTNLOQ1230= FTEMP

        RETURN
        END
