      SUBROUTINE HADNLOQ1230(Isch,IhadIn,XBJin,Qin,F1Min,F2Min,
     >  GLQin1,GRQin1,GLQin2,GRQin2,XMUin,HMASS,IPARTin,ISETin,FHAD123Q)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLOQ1230 to do convolution
C
C      25 JUNE 1998: IMPLEMENT THE MASSLESS CASE ALA FURMANSKI & PETRONZIO
C
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DIMENSION FHAD123Q(3)
       PARAMETER (PI=3.141592653589793)
       Common  / ActInt /  AERR, RERR, iActL, iActU
       COMMON /CXINTNLOQ1230/
     >   XBJ,XIMIN,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   Q,XMU,ISET,IRET,INDEX,IPART,Ihad
       integer Ftarget
       Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 

       EXTERNAL XINTNLOQ1230
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
        IPART= IPARTin
        ISET= ISETin
        Ihad =IhadIn


        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        if(Ftarget .eq. 1) then
          ETA   = 2.0* XBJ/(1.0 + Sqrt(Discr))
        else 
          ETA = XBJ
        endif 

        XIMIN= xbj   !*** debugging patch (Remove later)
        XIMIN= ETA   !*** Massless Parton Case:
        XIMIN= ETA*((F1M+F2M)**2+Q**2)/Q**2
        SHATTH= F2M**2
        Q2=Q*Q
        XIMIN=ETA*(Q2-F1M**2+SHATTH +
     >           DELTA(-Q2,F1M**2,SHATTH))/(2.0*Q2)

        XIMIN= xbj   !*** use massless kinematics
        if(isch.eq.9) XIMIN= Min(ETA*((F1M+F2M)**2+Q**2)/Q**2,1.0d0)
         XIMAX= 1.0

        DO INDEX=1,3,1
        Conv  = 0.0
        IF(XIMIN.LT.XIMAX)
     >    Conv = AdzInt (XINTNLOQ1230, XIMIN, XIMAX,
     >                   AERR, RERR, ErrEst, iErr, iActL, iActU)
          FHAD123Q(INDEX) = CONV
        ENDDO

        RETURN
        END
