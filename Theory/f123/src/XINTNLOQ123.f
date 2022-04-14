        Function XINTNLOQ123(XI)
C-----------------------------------------------------------------------------
C      COMPUTES CONVOLUTION OF PARTON HELICITIES WITH PDF'S
C      USED BY HADHEL PROGRAM
C
C      FIO 1998
C      ALA GOTTSCHALK, VAN DER BIJ ..., KRAMER ...
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
      DOUBLE COMPLEX  DiLog
      DOUBLE COMPLEX CARG,CSPENCE
       Dimension OmGNLO(-1:1)
       PARAMETER (PI=3.141592653589793)
       DATA SMALL /1.E-14/
       COMMON /CXINTNLOQ123/
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad

        XINTNLOQ123=0.0

        IHADRON=Ihad              !*** 2/28/05 FIO Fixed
        IGLUON=0

        XHAT = XIMIN/Xi
        XHAT0= XIMIN

        Z = XI

        IF(Z.GE.1.0-SMALL) Z=(1.0-SMALL)
C       IF (XHAT.GE.1.0-SMALL) XHAT=(1.0-SMALL)
C       IF (XHAT.LE.    SMALL) XHAT=(    SMALL)

        ALPHAS= AlQCD(XMU,iset)

        Q2=Q*Q
        X0=Q2/(Q2+F2M**2)
        IF(X0.GE.1.0-SMALL) X0=(1.0-SMALL)

        XLOG =  Log( Q2/( F1M**2 * X0 * Z * (1-X0*Z)) )
        XLOG1=  Log( Q2/( F1M**2 * X0     * (1-X0  )) )

C-----------------------------------------------------------------------------
C --- DELTA FUNCTION TERMS:

        CARG= DCMPLX( -X0/(1.-X0) )
        CSPENCE= DILOG(CARG)
        SPENCE= DBLE(CSPENCE)

         FDEL2= (3./2. + 2.* LOG(1.-XIMIN)) * (XLOG1 - 2.)
     >  + 1. - PI*PI/3. - 2.* SPENCE  +2.*LOG(1.-X0)

         FDEL21=  (1.-X0)/X0 * LOG(1.-X0)
         FDEL23=  (1.-X0)/X0 * LOG(1.-X0)

         FDEL1= (FDEL2 - FDEL21)/2.
         FDEL3= (FDEL2 - FDEL23)

         XMEASURE= (1.-XIMIN)   !*** INTEGRATION MEASURE

C-----------------------------------------------------------------------------
C --- REGULAR TERMS:

         FREG2= 2. + Z - (1.-X0)*Z/( 2.*(1.-X0*Z)**2 )
     >         -(1.-X0)*(Z+2)/(1.-X0*Z) + 1./(2.*(1.-X0*Z))
         FREG21= (1.+X0)*Z - 2.*(1.-X0)/(1.-X0*Z)
     >         + (1.-X0)**2 * Z**2 /( 1.-X0*Z)
         FREG23= 1.+Z-2.*(1.-X0)/(1.-X0*Z)

         FREG1= (FREG2 - FREG21)/2.
         FREG3= (FREG2 - FREG23)

C-----------------------------------------------------------------------------
C --- PLUS FUNCTION TERMS:

         FPLSZ= (1.+ Z**2)/Z  * (XLOG  -2.)
         FPLS1= (1.+1.**2)/1. * (XLOG1 -2.)

         Pdfz  = PDF(iSet, IHADRON, IPART, XHAT , XMU, iRet)
         Pdfz0 = PDF(iSet, IHADRON, IPART, XHAT0, XMU, iRet)

         FPLS2 = (FPLSZ*Pdfz - FPLS1*Pdfz0)/(1.-Z)

         FPLS1= (FPLS2 - 0.0)/2.
         FPLS3= (FPLS2 - 0.0)
C-----------------------------------------------------------------------------
C --- SUM TERMS:   (F2 is defined with an overall XBJ factor)

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

        XINTNLOQ123= FTEMP

        RETURN
        END
