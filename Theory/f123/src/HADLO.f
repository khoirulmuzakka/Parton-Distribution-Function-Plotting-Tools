      SUBROUTINE HADLO (IMASS,Isch,Ihad,IPARTON,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS,XHADLO)
C-----------------------------------------------------------------------------
C      Computes LO Hadron Helicities
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO(-1:1), XHADLO(-1:1)
        PARAMETER (PI=3.14159265359)
        integer Ftarget
        Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 

        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        Q2=Q**2

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        if(Ftarget.eq.1) then
          ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
        else 
          ETA = XBJ
        endif

        SHATTH= F2M**2
        XITH=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)


         if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
         ZERO=0.0
         IF(IMASS.EQ.0) THEN
            XITH=XBJ   !*** MASSLESS KINEMATICS
            if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
             CALL  PARLO (Q,ZERO,ZERO,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ELSE
             CALL  PARLO (Q,F1M ,F2M ,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ENDIF

        IHADRON= Ihad              !*** 2/28/05 FIO Fixed
        PdfTMP  = PDF(iSet, IHADRON, IPARTON, XITH, XMU, iRet)

        XHADLO( 1) = OmGLO( 1) * PDFTMP
        XHADLO( 0) = OmGLO( 0) * PDFTMP
        XHADLO(-1) = OmGLO(-1) * PDFTMP

        RETURN
        END
