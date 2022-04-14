C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       SUBROUTINE Fcc123(icharge,Mode, XBJ, Q,XMU, F123)
C-----------------------------------------------------------------------------
C      Program to COMPUTE K FACTORS
C      
C      
C      05/07/2007  Include Z-Z, and G-Z terms
C      05/03/2005  FIO update
C      02/03/2005  FIO update
C      02/05/01  FIO. 
C      23 Jan 2008: Include Short-Cut-2: 
C      
C----------------------------------------------------
C----------------------------------------------------
C      Ischeme= 0  !*** Massless MS-Bar
C      Ischeme= 1  !*** Full ACOT Scheme 
C      Ischeme= 2  !*** FFS 
C      Ischeme= 3  !*** Simplified ACOT Scheme 
C      Ischeme= 4  !*** Test Full ACOT Scheme (no NLO Q)
C      Ischeme= 5  !*** LO
C      Ischeme= 6  !*** Massless LO
c      Ischeme= 7  !*** Short-cut2: ACOT w/ Massless NLO-Q
C----------------------------------------------------
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       CHARACTER HEADER*78
       Dimension xout123(3), XXc(3), XXb(3), XMARRAY(6), CHARGE(6)
       Dimension F123(3), Xtot(3)
       PARAMETER (PI=3.14159265359)
       Common /Ischeme/ Isch, Iset, Iflg, Ihad
       Common  / ActInt /  AERR, RERR, iActL, iActU
       common /fred/ xmc,xmb,HMASS,SINW2, XMW, xmz

C-----------------------------------------------------------------------------
         INTEGER UPTYPE(3), DNTYPE(3)
         DIMENSION CKM(3,3) 

       DATA  (UPTYPE(I), I=1,3,1) /1,4,6/
       DATA  (DNTYPE(I), I=1,3,1) /2,3,5/
       DATA ((CKM(I,J), J=1,3,1), I=1,3,1)    !*** ABS OF CKM MATRIX
     1   /0.975268, 0.2209986, 0.00350000,
     2    0.2209541, 0.974422, 0.0409997, 
     3    0.00565041, 0.0407591, 0.999153/

C-----------------------------------------------------------------------------
C--- FOR CHARGED CURRENT
       DATA GLQ,GRQ
     >   /  1.0,0.0/
       DATA IPRINT /0/  !*** NO debugging
C       DATA IPRINT /1/ !*** for debugging
C-----------------------------------------------------------------------------
C                          U      D      S      C      B      T
       DATA XMARRAY  /    0.1,   0.1,   0.2,   1.3,   4.5, 175.0/
       DATA CHARGE   /   +2.0,  -1.0,  -1.0,  +2.0,  -1.0,  +2.0/

C F: nice different values from Fnc
       XMARRAY(4)=xmc
       XMARRAY(5)=xmb

C----------------------------------------------------------------------
C INITIALIZATION
C----------------------------------------------------------------------
          DO I=1,3,1
             F123( I)  = 0.0d0
             xtot( I)  = 0.0d0
             xxc(  I)  = 0.0d0
             xxb(  I)  = 0.0d0
          ENDDO
C----------------------------------------------------------------------
C CHOOSE W+ or W- SCATTERING
C----------------------------------------------------------------------
      IUPDOWN  =-1    !*** W- SCATTERING
      IUPDOWN  =+1    !*** W+ SCATTERING
      IUPDOWN = ICHARGE
      UPDOWN=float(ICHARGE)  !*** CHANGE TO FLOAT
C----------------------------------------------------
C   *************************************************
C----------------------------------------------------
C----------------------------------------------------
C   ***  LOOP OVER PARTON FLAVORS: 
C----------------------------------------------------
      DO IQUARKIN  =1,3,1
      DO IQUARKOUT =1,3,1

C----------------------------------------------------
      IF(UPDOWN.GE.0)  THEN  !*** THIS IS FOR DOWN -> UP
           IPARTIN = +DNTYPE(IQUARKIN )
           IPARTOUT= +UPTYPE(IQUARKOUT)

           IUP= IQUARKOUT   !***  USED TO PICK UP CORRECT CKM TERM
           IDN= IQUARKIN

      ELSE                 !*** THIS IS FOR UP -> DOWN
           IPARTIN = +UPTYPE(IQUARKIN )
           IPARTOUT= +DNTYPE(IQUARKOUT)

           IUP= IQUARKIN    !***  USED TO PICK UP CORRECT CKM TERM
           IDN= IQUARKOUT

      ENDIF
      COUPLING = CKM(IUP,IDN)**2
C----------------------------------------------------------------------
C SETUP PARTON MASSES: 
C----------------------------------------------------------------------
      IF(IPRINT.NE.0) WRITE(6,*) ' PARTON ',IPARTIN,' => ',IPARTOUT

C----------------------------------------------------------------------
C SETUP PARTON MASSES: 
C----------------------------------------------------------------------
       F1M=XMARRAY(ABS(IPARTIN))
       F2M=XMARRAY(ABS(IPARTOUT))

C----------------------------------------------------
C NOTE: ASSUME TOP IS ZERO (NOT IMPLEMENTED IN ALL PDF SETS)
C----------------------------------------------------------------------
      if((abs(IPARTIN).NE.6).and.(abs(IPARTOUT).NE.6)) then
         call TOTF(XBJ,Q,F1M,F2M,GLQ,GRQ,GLQ,GRQ,IPARTIN,IPARTOUT,
     >   XMU,ISET,Ihad,HMASS,Isch,Iflg,Icharge,Xout123) 
      else
           DO I=1,3,1
             Xout123( I)  =0.0
           ENDDO
      endif
C----------------------------------------------------

C----------------------------------------------------
C   ***  ADD PARTON CONTRIBUTION TO TOTAL STRUCTURE FUNCTIONS
C----------------------------------------------------
           DO I=1,3,1
             xtot( I)  = xtot( I)  + COUPLING * Xout123( I) 
           ENDDO

C----------------------------------------------------
C   ***  PICK OUT CHARM AND SAVE THIS
C   ***  don't forget: s->c sums over cbar->sbar
C----------------------------------------------------
       IF((ABS(IPARTin).EQ.4).or.(ABS(IPARTOUT).EQ.4)) THEN    
           DO I=1,3,1
             XXc(I)  = XXc(I) + COUPLING * Xout123(I)    
           ENDDO
        ENDIF
C----------------------------------------------------
C   ***  PICK OUT BOTTOM AND SAVE THIS
C   ***  don't forget: s->c sums over cbar->sbar
C----------------------------------------------------
       IF((ABS(IPARTin).EQ.5).or.(ABS(IPARTOUT).EQ.5)) THEN    
           DO I=1,3,1
             XXb(I)  = XXb(I) + COUPLING * Xout123(I)    
           ENDDO
        ENDIF
C----------------------------------------------------
C----------------------------------------------------
      ENDDO  !***  IQUARKOUT
      ENDDO  !***  IQUARKIN
C----------------------------------------------------
C   ***  END LOOP OVER PARTON FLAVORS: 
C----------------------------------------------------

C----------------------------------------------------
C   ***  pass up proper result
C----------------------------------------------------
           DO I=1,3,1
             IF(MODE.EQ.1) f123( I)  = xtot(I) 
             IF(MODE.EQ.2) f123( I)  = XXc( I) 
             IF(MODE.EQ.3) f123( I)  = XXb( I) 
           ENDDO

      Return
      END

C----------------------------------------------------------------------
