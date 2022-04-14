      SUBROUTINE HADSUBQ123(Isch,IhadIn,IPARTONin,XBJin,Q,F1Min,F2Min,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMUin,ISETin,HMASS,XHADSUBQ123)
C-----------------------------------------------------------------------------
C      Computes SUB Hadron Helicities
C
C
C
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
        Dimension OmGLO(-1:1),OmGLO123(3),XHADSUBQ123(3)
        EXTERNAL XINTSUBQ123
        PARAMETER (PI=3.14159265359)
        Common  / ActInt /  AERR, RERR, iActL, iActU
        Common /CINTSUBQ/ F1M,XBJ, XMU, XiTH, iSet, IPARTON,Ihad
        integer Ftarget
        Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 
        DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        XHADSUBQ123(1) = 0.0
        XHADSUBQ123(2) = 0.0
        XHADSUBQ123(3) = 0.0

        IF(XMUin.LE.F1Min) RETURN

        F1M    =F1Min
        F2M    =F2Min
        XBJ    =XBJin
        XMU    =XMUin
        iSet   =ISETin
        IPARTON=IPARTONin
        Ihad   =Ihadin

        Q2=Q**2

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        if (Ftarget.eq.1) then
          ETA  = 2.0* XBJ/(1.0 + Sqrt(Discr))
        else
          ETA = XBJ
        endif

        SHATTH= F2M**2
        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        XIMIN=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)
        XIMAX= 1.0

        xith=XIMIN

       ALPHAS = AlQCD(XMU,iset)
       IF(ALPHAS.GT.0.5) ALPIX=0.5
       ComFct = ALPHAS/(2.0*PI)

       XIMAX=1.0
       XIMIN=XITH

C----------------------------------------------------------------------
        Conv = 0.0
        IF(XITH.LT.XIMAX)
     >    Conv = AdzInt (XINTSUBQ123, XIMIN, XIMAX,
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)

        Conv = Conv * ComFct
C----------------------------------------------------------------------

        CALL  PARLO (Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,OmGLO)

        XFAC= 1./(1.+(F2M/Q)**2)

        OmGLO123(1)= ( OmGLO(+1) + OmGLO(-1) )/2.
        OmGLO123(2)= ( OmGLO(+1) + OmGLO(-1) + 2.0*OmGLO(0) )*XFAC
        OmGLO123(3)=  -OmGLO(+1) + OmGLO(-1)
C----------------------------------------------------------------------

      RHO=1.0
      IF(HMASS.GT.0.0) THEN
        XNU=Q**2/(2.0*HMASS*XBJ)
        RHO=SQRT(1.0+(Q/XNU)**2)
      ENDIF


        XHADSUBQ123(1) = OmGLO123(1) * Conv
        XHADSUBQ123(2) = OmGLO123(2) * Conv * XBJ    !*** Def of F2
        XHADSUBQ123(3) = OmGLO123(3) * Conv


        RETURN
        END
