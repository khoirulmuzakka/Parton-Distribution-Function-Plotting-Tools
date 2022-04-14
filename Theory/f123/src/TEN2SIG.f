      Double Precision Function TEN2SIG(XSIGN,ICHARGED,
     >   XBJ,Q,GLLEP,GRLEP,XLEPOL,
     >   HMASS,SBIG,F123)
C-----------------------------------------------------------------------------
C     Compute dSig/dx/dy
C     Arguments modified: FIO 8/15/02
C      05/07/2007  Include Z-Z, and G-Z terms
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       Dimension F123(3)
       PARAMETER (PI=3.14159265359)

         TEN2SIG= 0.0

C       ICHARGED     =  0 FOR NC;  NON-ZERO FOR CC
C       XSIGN= +1.0  !*** THIS IS FOR      LEPTON SCATTERING.
C       XSIGN= -1.0  !*** THIS IS FOR ANTI-LEPTON SCATTERING.
C  **** ANTI-LEPTON SCATTERING IS OBTAINED BY GLLEP <=> GRLEP _OR_ BY FLIPPING XSIGN
C
C  ***** Z TERMS NOT YET IMPLEMENT HERE !!!!!!!!!!!!!
C

         E1=(SBIG-HMASS**2)/(2.0*HMASS)
         XNU=Q**2/(2.0*HMASS*XBJ)
         Y=XNU/E1
         E2=E1*(1.0-Y)

         IF(Y.GE.1.0) RETURN

         ALPHA=1.0/137.
         EEM2=ALPHA*4.0*PI
         G1PHOTON= EEM2/Q**2

         VEV = 246.0
         WMASS=80.0
         G= (2.0*WMASS)/VEV
         GBW= G/(2.0*SQRT(2.0))
         G1W= GBW**2/(Q**2+WMASS**2)


      IF(ICHARGED.LE.0.0) THEN
C        ******** FOR NEUTRAL CURRENT ONLY *********
C  ***** Z TERMS NOT YET IMPLEMENT HERE !!!!!!!!!!!!!
         XFAC= 2.0 * HMASS * E1/(Pi) * G1PHOTON**2 /XLEPOL
       ELSE
C        ******** FOR CHARGED CURRENT ONLY *********
         XFAC= 2.0 * HMASS * E1/(Pi) * G1W**2 /XLEPOL
C        XFAC= XFAC /(1+Q**2/WMASS**2)**2  !*** remove: FIO 18 mar 2002
       ENDIF

         GPLUS  = GLLEP**2 + GRLEP**2
         GMINUS = GLLEP**2 - GRLEP**2

         TERM1= XBJ * F123(1) * Y*Y
         TERM2= F123(2) *( (1.-Y)-(HMASS*XBJ*Y)/(2.0*E1) )
         TERM3= XBJ * F123(3) * Y * (1.-Y/2.)
         SUM= GPLUS * (TERM1+TERM2) + XSIGN * GMINUS * TERM3

         XSEC= XFAC * SUM
         TEN2SIG= XSEC

        RETURN
        END
