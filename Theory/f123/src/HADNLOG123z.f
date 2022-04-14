      SUBROUTINE HADNLOG123z(Isch,IhadIn,XBJin,Qin,F1Min,F2Min,
     >   GLQin1,GRQin1,GLQin2,GRQin2,XMUin,ISETin,HMASS,FHAD123G2z)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLOG123z to do convolution
C
C      14 DECEMBER 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DIMENSION FHAD123G2z(3)
       PARAMETER (PI=3.141592653589793)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLOG123z/
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,Ihad
       integer Ftarget
       Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 

       EXTERNAL XINTNLOG123z
       DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))

        F1M = F1Min
        F2M = F2Min
        GLQ1 = GLQin1
        GRQ1 = GRQin1
        GLQ2 = GLQin2
        GRQ2 = GRQin2
        XBJ = XBJin
        Q   = Qin
        XMU = XMUin
        ISET= ISETin
        Ihad=IhadIn


        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        if (Ftarget .eq. 1) then
          ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
        else 
          ETA = XBJ
        endif

        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA   !*** Massless Parton Case:
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2  !*** Full Case:
        SHATTH= (F1M+F2M)**2
        Q2=Q*Q
        XIMIN=ETA*(Q2-0.0**2+SHATTH +
     >           DELTA(-Q2,0.0d0**2,SHATTH))/(2.0*Q2)

        XIMIN= xbj   !*** use massless kinematics
        XIMAX= 1.0

        DO INDEX=1,3,1
        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)
     >    Conv = AdzInt (XINTNLOG123z, XIMIN, XIMAX,
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)
          FHAD123G2z(INDEX) = CONV
        ENDDO

        RETURN
        END
