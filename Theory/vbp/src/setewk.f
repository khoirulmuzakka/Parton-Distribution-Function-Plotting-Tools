      Subroutine SetEwk
C                                                   -=-=- setewk

C  These comments are enclosed in the lead subprogram to survive forsplit

C ====================================================================
C GroupName: SetEwk
C Description: initial setup of the EwkPac.
C ListOfFiles: setewk filcpl filkmm
C ====================================================================

C #Header: /Net/d2a/wkt/1hep/2ewk/RCS/SetEwk.f,v 6.1 98/08/16 17:24:43 wkt Exp $
C #Log:	SetEwk.f,v $
c Revision 6.1  98/08/16  17:24:43  wkt
c Re-organization; rationalization; initialization for DIS & DY
c 
C   In DisPac: QCPLN ---> EwCplSc
C   In VbpPac: SQCPAN --> EwCpl2An
C              SQCPCM --> EwCpl2Cn
c 
C Revision 6.0  96/12/25  14:38:28  wkt
C Synchronize with version 6 of all the other pac's

c Revision 1.1  94/02/22  11:10:59  lai
c Initial revision
c 
c Revision 2.2  92/03/08  15:21:24  wkt
c Setxxx.f uniformly added to force blockdata linking and 
c to perform other initiation Functions.
c 
C                 Setup the Ewkpac, which contains the following modules

C      Subroutine SetEwk
C      Subroutine FILCPL (CH, GV, GA, GL, GR)
C      Subroutine FILKMM 
C
C      Function ALFEwk (IBOSON)
C      Entry VBNMAS (IBOSON)
C      Entry SWG2F ()
C      Entry GEWLT (IT1, IBS, IBT, IT2)
C      Entry GEWLH (IT1, IBS, IBH, IT2)
C      Entry GEWQT (IQ1, IBS, IBT, IQ2)
C      Entry GEWQH (IQ1, IBS, IBH, IQ2)
C
C      Subroutine SetEwCpl2
C      Function SQCPAN (JP1, JBN, JP2)
C      Entry SQCPCM (JP, JBN)
C      Function EWCPL0 (IBSN, IP1, IP2, IRT)   (Inactive)

C     Set up the basic Electro-weak coupling matrices for the Boson-Fermion
C     Yukawa coupling term in the Effective Lagrangian

C     Boson label: (IBN)   1,   2,   3,   4
C                       gamma   W+   W-   Z    
  
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN)
 
      DIMENSION IP(NSP, NGN), CHL(NSP), CHQ(NSP)
      DIMENSION GGQA(NSP,NGN,NSP,NGN), GGQV(NSP,NGN,NSP,NGN),
     >          GGQR(NSP,NGN,NSP,NGN), GGQL(NSP,NGN,NSP,NGN)
 
      External DatEwk

      COMMON / IOUNIT / NIN, NOUT, NWRT
      COMMON / VBPEwkPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / KMATRX / VKM (NGN, NGN)

C            Left, Right, Vector, and Axial-vector couplings of Leptons & Quarks

      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 
      DATA IP, Lprt   / 1, 2, 4, 3, 6, 5, 0 /
      DATA CHL,CHQ / 1.0, 0.0, 0.666666667, -0.333333333 /
 
C     Adopt the convention:

C                |  IBSN      coupling       current           coupling**2/4pi
C                |
C   ------------ |------------------------------------------------------------
C      Photon    |    1          e             Gmu            ALFE = Alpha(QED) 
C                |
C      W+/W-     |   2/3      g/2Sqrt2      Gmu(1-G5)          ALFE / SWG2 / 8
C                |
C        Z       |    4        g/2CosW  Gmu((T3-2QSWG2)-T3G5)   ALFE / S2WG2
C                |
C   ------------ |-------------------------------------------------------------
C     where Gmu = Gamma-super-mu; G5 = Gamma-5; Cos W = Cos(Theta-Weinberg);
C     SWG2 = (Sin W)**2; & S2WG2 = (Sin 2W)**2.
C   ---------------------------------------------------------------------------
C                                       Fill the overall coupling**2 array ALFEW
      ALFEW (1) = ALFE
      ALFEW (2) = ALFE / SWG2 / 8.
      ALFEW (3) = ALFEW (2)
      ALFEW (4) = ALFE / SWG2 / (1.-SWG2) / 4.
C                                                   ----------------------------
C                                                       For Each Generation:
C                                     Fill the relative coupling in the currents
C                           GXL/R = (GXV -/+ GXA)/2.      (X = Lepton / Quark) ;
C                             hence the currents are:   GXL/R * Gmu * (1-/+G5)

      CALL FILCPL (CHL, GLV, GLA, GLL, GLR)
      CALL FILCPL (CHQ, GQV, GQA, GQL, GQR)
C                                                  ----------------------------
C                          For Quarks, put in Generations and generation-mixing

C                                 fill in the KM matrix array from data in the
C                           MIXing PARameter common block / MIXPAR / consisting
C                                 of Cosines of mixing angles and Phase angles.
      CALL FILKMM

C                     For W+ / W- , Multiply the KM matrix into the Down quarks
      DO 10 IG1 = 1, NGN
      DO 10 IG2 = 1, NGN
      DO 10 IS1 = 1, NSP
      DO 10 IS2 = 1, NSP
        GGQA (IS1, IG1, IS2, IG2) = GQA(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQV (IS1, IG1, IS2, IG2) = GQV(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQR (IS1, IG1, IS2, IG2) = GQR(IS1, 2, IS2) * VKM(IG1, IG2)
        GGQL (IS1, IG1, IS2, IG2) = GQL(IS1, 2, IS2) * VKM(IG1, IG2)
   10 CONTINUE
C                                                  ----------------------------
C                         (Generation, ISPin) labels  ----> parton flavor label
      DO 20 IG1 = 1, NGN
      DO 20 IG2 = 1, NGN
      DO 20 IS1 = 1, NSP
      DO 20 IS2 = 1, NSP
        IP1 = IP(IS1, IG1)
        IP2 = IP(IS2, IG2)

        DO 21 IBN = 1, NBN
        IF (IBN .EQ. 1 .OR. IBN .EQ. 4) THEN
         IF (IG1 .NE. IG2) THEN
            HQA (IP1, IBN, IP2) = 0.0
            HQV (IP1, IBN, IP2) = 0.0
            HQL (IP1, IBN, IP2) = 0.0
            HQR (IP1, IBN, IP2) = 0.0
         ELSE
            HQA (IP1, IBN, IP2) = GQA(IS1, IBN, IS2)
            HQV (IP1, IBN, IP2) = GQV(IS1, IBN, IS2)
            HQL (IP1, IBN, IP2) = GQL(IS1, IBN, IS2)
            HQR (IP1, IBN, IP2) = GQR(IS1, IBN, IS2)
         ENDIF
        ENDIF
   21   CONTINUE

            HQA (IP1, 2, IP2) = GGQA(IS1, IG1, IS2, IG2)
            HQV (IP1, 2, IP2) = GGQV(IS1, IG1, IS2, IG2)
            HQL (IP1, 2, IP2) = GGQL(IS1, IG1, IS2, IG2)
            HQR (IP1, 2, IP2) = GGQR(IS1, IG1, IS2, IG2)

            HQA (IP1, 3, IP2) = GGQA(IS2, IG2, IS1, IG1)
            HQV (IP1, 3, IP2) = GGQV(IS2, IG2, IS1, IG1)
            HQL (IP1, 3, IP2) = GGQL(IS2, IG2, IS1, IG1)
            HQR (IP1, 3, IP2) = GGQR(IS2, IG2, IS1, IG1)

   20 CONTINUE
C                                                  ----------------------------
C                                                                      Finished
C                          Set Lprt>0 to print Quark couplings for confirmation
      If (Lprt .ge. 1) Then
        Print '(A)', 
     >  'Flavor-dependent part of Quark-VectorBoson Vertex:'
        Do 41 Ibn = 1, 4
	Print '(/A, I4)', 'Iboson = ', Ibn
	Print '(/A/(4F12.5))', 'Axial-v: Iflv=1-4:'
     >, ((HqA(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
	Print '(/A/(4F12.5))', 'Vector: Iflv=1-4:'
     >, ((HqV(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
	Print '(/A/(4F12.5))', 'Left-hand: Iflv=1-4:'
     >, ((HqL(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
	Print '(/A/(4F12.5))', 'Right-hand: Iflv=1-4:'
     >, ((HqR(Ip1,Ibn,Ip2), Ip1=1,4), Ip2=1,4)
 41     Continue
      Endif

      Return
C                        ****************************
      End

      BLOCK DATA DATEwk
C                       
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
      PARAMETER (NSP = 2, NGN = 3, NBN = 4, NPOL = 3)
      PARAMETER (NFL = NSP * NGN, Ntt1=NSP*NBN*NSP, Ntt2=NFL*NBN*NFL)
C      PARAMETER (NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)
      PARAMETER (NANG= NGN*(NGN-1)/2,NGNm1=(NGN-1),NPHS=NGNm1*(NGN-2)/2)
 
      COMMON / VBPEwkPAR / WMS, ZMS, SWG2, ALFE, ALFEW(NBN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)
      COMMON / EW1LCP / GLL(NSP,NBN,NSP), GLR(NSP,NBN,NSP), 
     >                  GLV(NSP,NBN,NSP), GLA(NSP,NBN,NSP)
      COMMON / EW1QCP / GQL(NSP,NBN,NSP), GQR(NSP,NBN,NSP), 
     >                  GQV(NSP,NBN,NSP), GQA(NSP,NBN,NSP)
      COMMON / EW2QCP / HQL(NFL,NBN,NFL), HQR(NFL,NBN,NFL), 
     >                  HQV(NFL,NBN,NFL), HQA(NFL,NBN,NFL)
 
c      DATA WMS, ZMS, SWG2, ALFE / 80.41, 91.187, 0.23124, 7.297353E-3 /
      Data (((GLL(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /
     >     (((GLR(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GLV(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GLA(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 / 
C
     >     (((GQL(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /
     >     (((GQR(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GQV(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /    
     >     (((GQA(i,j,k), i=1,nsp), j=1,nbn), k=1,nsp) / Ntt1*0D0 /   
C
     >     (((HQL(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQR(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQV(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
     >     (((HQA(i,j,k), i=1,nfl), j=1,nbn), k=1,nfl) / Ntt2*0D0 /    
 
C                        ****************************
      END

