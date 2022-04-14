      Function FHADHEL (IHEL,Isch,IhadIn,XBJ,Q,F1M,F2M,
     >   GLQ1,GRQ1,GLQ2,GRQ2,XMU,ISET,HMASS)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLO to do convolution
C
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
       PARAMETER (PI=3.14159265359)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLO/ ETA,Qout,F1Mout,F2Mout,
     >   GLQout1,GRQout1,GLQout2,GRQout2,
     >   XMUout,ISETout,IRETout,IHELout,Ihad
       EXTERNAL XINTNLO
       integer Ftarget
       Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 

        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        if (Ftarget.eq.0d0) then
          ETA = XBJ
        else 
          ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
        endif

        IHELout=IHEL
        Qout   =Q
        F1Mout =F1M
        F2Mout =F2M
        GLQout1  =GLQ1
        GRQout1  =GRQ1
        GLQout2  =GLQ2
        GRQout2  =GRQ2
        XMUout =XMU
        ISETout=ISET
        Ihad   =Ihadin

        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        XIMAX= 1.0

        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)
     >  Conv = AdzInt (XINTNLO, XIMIN, XIMAX,
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)

        FHADHEL = CONV

        RETURN
        END
