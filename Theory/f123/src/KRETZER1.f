       SUBROUTINE KRETZER1()
C      PROGRAM KRETZER1
C-----------------------------------------------------------------------------
C      Program to COMPUTE DIS HQ PRODUCTION
C
C
C
C
C      11/05/99 FIO. MODIFY FRONT END
C      01/01/99 FIO. CORE CODE BY STEFAN KRETZER
C      S. Kretzer, I. Schienbein. Phys.Rev.D58:094035,1998  hep-ph/9805233
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      dimension acotlo(3), rterm(3), vterm(3)
      DIMENSION FHAD123Q(3)
      DIMENSION XMARRAY(6)
      PARAMETER (PI=3.141592653589793)
C-----------------------------------------------------------------------------
C--- ONLY COMMON BLOCK USED
      Common  / ActInt /  AERR, RERR, iActL, iActU
      common /fred/ xmc,xmb,HMASS,SINW2, XMW, xmz

C-----------------------------------------------------------------------------
C--- FOR NEUTRAL CURRENT
C      DATA SBIG,XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
c     >GLLEP,GRLEP,XMU,ISET,HMASS
C    >   /  98596.,0.1, 10.,1.6,1.6,0.5,0.50.5,0.5,,0.5,0.5,10,1,0.938/
C      DATA   XLEPOL /2.0/
C      DATA IPARTIN, IPARTOUT,SCALE  /  4, 4,  -1  /
C-----------------------------------------------------------------------------
C--- FOR CHARGED CURRENT
      DATA SBIG,XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,
     >   GLLEP,GRLEP,XMU,ISET,HMASS
     >   / 98596.,0.1, 10.,0.5,1.6,1.0,0.0,1.0,0.0,1.0,0.0,10.0,1,0.938/
      DATA   XLEPOL /1.0/
      DATA IPARTIN, IPARTOUT,SCALE  /  3, 4,  -1  /
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
       DATA XMARRAY  /  0.1, 0.1, 0.2, 1.6, 5.0, 170.0/
       DATA  AERR, RERR, iActL, iActU    /  0.0, 0.001,1,1/
       XMARRAY(4)=xmc
       XMARRAY(5)=xmb

C----------------------------------------------------------------------
C----------------------------------------------------------------------
C SETUP INTEGRATION PARAMETERS
C----------------------------------------------------------------------
C     RERR=1.E-8   !*** Setting this too small will yield unstable results ***
      WRITE(6,202)  AERR, RERR, iActL, iActU
202   FORMAT(' AERR= ',1PG14.7,' RERR= ',1PG14.7,
     >       ' IACTL= ',I3,' IACTU= ',I3)
      WRITE(6,*) ' ENTER:  AERR, RERR, iActL, iActU'
      READ (5,*)           AERR, RERR, iActL, iActU
C----------------------------------------------------------------------
C SETUP PARTON NUMBERS
C----------------------------------------------------------------------
      WRITE(6,203) IPARTIN, IPARTOUT
203   FORMAT(' IPARTIN= ',I3,' IPARTOUT= ',I3)
      WRITE(6,*) ' ENTER: IPARTIN, IPARTOUT '
      READ (5,*)          IPARTIN, IPARTOUT

       F1M=XMARRAY(ABS(IPARTIN))
       F2M=XMARRAY(ABS(IPARTOUT))

       IF(IPARTIN.EQ.IPARTOUT) THEN
         GLQ1=0.5
         GRQ1=0.5
         GLQ2=0.5
         GRQ2=0.5
         GLLEP=0.5
         GRLEP=0.5
         XLEPOL=2.0
       ENDIF

       WRITE(6,*)    GLQ1,GRQ1,GLQ2,GRQ2,GLLEP,GRLEP,SBIG,XLEPOL
       WRITE(6,*) '> GLQ1,GRQ1,GLQ2,GRQ2,GLLEP,GRLEP,SBIG,XLEPOL'
       READ( 5,*)    GLQ1,GRQ1,GLQ2,GRQ2,GLLEP,GRLEP,SBIG,XLEPOL

       E=(SBIG-HMASS**2)/(2.0*HMASS)

C----------------------------------------------------------------------
C SET THE LOOP
C----------------------------------------------------------------------
1     WRITE(6,102) XBJ,Q, E,F1M,F2M,ISET,HMASS,SCALE
102   FORMAT('   XBJ= ',1PG14.7,'     Q= ',1PG14.7,/,
     >       '     E= ',1PG14.7,'   F1M= ',1PG14.7,
     >       '   F2M= ',1PG14.7,'  ISET= ',I3, /,
     >       ' HMASS= ',1PG14.7,' SCALE= ',1PG14.7)
      WRITE(6,*) ' ENTER:  XBJ,Q, E,F1M,F2M,ISET,HMASS,SCALE  '
      READ (5,*)           XBJ,Q, E,F1M,F2M,ISET,HMASS,SCALE
C----------------------------------------------------------------------
C SET THE SCALES
C----------------------------------------------------------------------
* ... kinematics:
      xb = XBJ
      Q2 = Q*Q
      m1 = F1M
      m2 = F2M
      m1s = m1 * m1
      m2s = m2 * m2
C----------------------------------------------------
C   ***  COMPUTE MU SCALE
C----------------------------------------------------
       IF(SCALE.GT.0.0) THEN
           XMU = SCALE
       ELSEIF(SCALE.LT.0.0) THEN
           XMU = Q * Abs(scale)
       ENDIF
      Q2f = XMU*XMU
C      if ( Q2f .le. m2s ) Q2f = m2s
      Q2r = Q2f


C----------------------------------------------------
C   ***  call function
C----------------------------------------------------

      Call HADNLOQ123K(Isch,Ihad, XBJ,Q,F1M,F2M,GLQ1,GRQ1,GLQ2,GRQ2,XMU
     >                ,HMASS,IPARTIN,ISET,FHAD123Q)


C----------------------------------------------------
C   ***  PRINT OUT
C----------------------------------------------------

       WRITE(6 ,103) E,XBJ,Y,HMASS,ISET,
     >             SCALE,Q,Q*Q,W,XMU,RHO2
     >        ,(   FHAD123Q(N),N=1,3,1),   STOT


C----------------------------------------------------
103    FORMAT(/,
     >   '---------------------------------------------------',/,
     >   ' E    =',1PG14.7,' XBJ  =',1PG14.7,' Y   =',1PG14.7,/,
     >   ' HMASS=',1PG14.7,' JSET =',i4,10x ,' SCALE=',1PG14.7,/,
     >   ' Q    =',1PG14.7,' Q2   =',1PG14.7,' W    =',1PG14.7,/,
     >   ' FMU  =',1PG14.7,' RHO2 =',1PG14.7,/,
     >   '---------------------------------------------------',/,
     >   ' FHAD123Q   (123): ',4(1PG14.7,1X),/,
     >   '---------------------------------------------------'
     >   )
C----------------------------------------------------------------------
      GOTO 1
C----------------------------------------------------------------------

      end
