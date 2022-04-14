      SUBROUTINE HADSUB (IMASS,Isch,IhadIn,IPARTON,XBJ,Q,F1Min,F2Min,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMUin,ISETin,HMASS,XHADSUB)
C-----------------------------------------------------------------------------
C      Computes SUB Hadron Helicities:
C           Gluon initiated only; can not do quark with parallel structure: FIO 2/17/99
C
C      25 June 98: Moved XHADSUBQ into separate module
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO(-1:1), XHADSUB(-1:1)
        EXTERNAL XINTSUB
        PARAMETER (PI=3.14159265359)
        Common  / ActInt /  AERR, RERR, iActL, iActU
        Common  / CINTSUB /  XMU, XiTH, iSet, IPARTONout, Ihad
        integer Ftarget
        Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 
 
        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        XHADSUB( 1) = 0.0
        XHADSUB( 0) = 0.0
        XHADSUB(-1) = 0.0

C        write(6,*) 'xmuin, f1min, XMUin.LE.F1Min',
C     >   xmuin, f1min, XMUin.LE.F1Min
        IF(XMUin.LE.F1Min) RETURN

        F1M    =F1Min
        F2M    =F2Min
        IPARTONout=IPARTON
        iSet =ISETin
        XMU  =XMUin
        Ihad =IhadIn

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

       ALPIX = AlQCD(XMU,iset)/pi
       DivFct = ALPIX * Log((XMU/F1m)**2) / 2.0
C  *** changed from 4.0 to 2.0 and move 2.0 into XINTSUB  FIO 12/7/95

         if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)

         ZERO=0.0
         IF(IMASS.EQ.0) THEN
             XITH=XBJ   !*** MASSLESS KINEMATICS
             if(isch.eq.9) xith=Min(eta*(1.0d0+(f1m+f2m)**2/Q2),1.0d0) !*** For S-ACOT(CHI)
             CALL  PARLO (Q,ZERO,ZERO,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ELSE
             CALL  PARLO (Q,F1M ,F2M ,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)
         ENDIF

       XIMAX=1.0
       IF(XITH.LT.XIMAX) THEN
          Convol = AdzInt (XINTSUB, XITH, XIMAX,
     >                 AERR, RERR, ErrEst, iErr, iActL, iActU)
       ELSE
          Convol  = 0.0
       ENDIF

        PDFSUB = DIVFCT * CONVOL


        XHADSUB( 1) = OmGLO( 1) * PDFSUB
        XHADSUB( 0) = OmGLO( 0) * PDFSUB
        XHADSUB(-1) = OmGLO(-1) * PDFSUB

        RETURN
        END
