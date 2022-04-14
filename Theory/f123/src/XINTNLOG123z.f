        Function XINTNLOG123z(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADNLOG1230 PROGRAM
C
C      14 DECEMBER 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DOUBLE COMPLEX  DiLog
       DOUBLE COMPLEX CARG,CSPENCE
       PARAMETER (PI=3.141592653589793)
       DATA SMALL /1.E-14/
       COMMON /CXINTNLOG123z/
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,Ihad

        XINTNLOG123z=0.0

        IHADRON=Ihad             !*** 2/28/05 FIO Fixed
        IGLUON=0

        Z = XI
        Z0 = XIMIN
        ZHAT = Z0/Z

        IF(Z.GE.1.0-SMALL) Z=(1.0-SMALL)

        ALPHAS= AlQCD(XMU,iset)

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS:

C   there are none for the gluon case

C-----------------------------------------------------------------------------
C --- REGULAR TERMS:

         FREG2  = 2 * 1/2. *
     >          ((z**2 + (1-z)**2)*LOG((1-z)/z) - 1 + 8*z*(1-z))

         FREG21 = + 2 * 1/2. * (4 * z * (1-z))


         FREG1= (FREG2 - FREG21)/2.0
         FREG3= 0.0

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS:

C   there are none for the gluon case

C-----------------------------------------------------------------------------
C --- SUM TERMS:   (Standard F2 is defined with an overall XBJ factor)

        igluon=0
         Pdfz  = PDF(iSet, IHADRON, igluon, ZHAT , XMU, iRet)

          FTEMP1= FREG1*Pdfz/Z
          FTEMP2= FREG2*Pdfz/Z * XBJ
          FTEMP3= FREG3*Pdfz/Z

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

        FTEMP= FTEMP * ALPHAS/(2.0*PI)
        FTEMP= FTEMP *2.0   !*** TO MATCH OUR NORMALIZATION

        XINTNLOG123z= FTEMP

        RETURN
        END
