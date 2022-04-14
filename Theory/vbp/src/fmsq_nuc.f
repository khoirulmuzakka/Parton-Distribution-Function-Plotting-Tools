c ... 03.05.2007 (IJS/JYY)
c ... fmsq_nuc.f

      FUNCTION FMSQ_NUC(ZZ)
C                                                   -=-=- fmsq_nuc

C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: VbpAux
C Description: Auxilary modules
C ListOfFiles: fmsq_nuc fqq12 fqrhd
C===========================================================================

C #Header: /Net/d2a/wkt/1hep/3vbp/RCS/VbpAux.f,v 1.4 98/08/16 21:21:03 wkt Exp $ 
C #Log:	VbpAux.f,v $
c Revision 1.4  98/08/16  21:21:03  wkt
c Modifications induced by changes from EwkPac5 to EwkPac6.
c   SetEwk changed; Setup and Call of EW coupling constants changed.
c 
c Revision 1.3  98/08/12  14:10:37  wkt
c trivial /data/ statements added to avoid compiler warnings.
c 
c Revision 1.2  98/08/10  12:02:37  wkt
c Incompatible common block PRQCDC corrected
c 
c Revision 1.1  98/08/05  23:53:46  wkt
c Initial revision
c 

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      Logical lm,ll

      COMMON 
     1 / STCOM_NUC / RTS, Q, AMU, NBM, NTG, LSET1, LSET2, iOrd, iBsn
     1 / LOVBP_NUC / LM, LL, Ischm, Iscal
     1 / F0COM_NUC / XA, XB, FA0(-10:10), FB0(-10:10), IA, IB
     1 / E866xf/ ixfx

      DATA CF, TR / 1.33333333333333, 0.5 /
      
      Z = ZZ
      IF (ABS(1.-Z).LT.1E-15 .OR. Z.LT.1E-15) THEN
        PRINT *, 'Z out of range in FMSQ_NUC; Z =', Z
        call PAUSE('in fmsq_nuc.f')  
      ENDIF
cjfo   mods for E866 xf
      xffac=1.
      if(ixfx.eq.1)xffac=1./(xa+xb)
      
      TMP = 0.
C                                 Add mu dependent part: essentially Pff of A-P
      IF (LL) TMP = LOG(Q/AMU) * CF * (1. + Z**2)/(1.-Z)
C
C                                 Add DIS to MSbar part: essentially C2q Wilson
      IF (LM) TMP = TMP + CF / 2. *
     >              ( ( 1.+ Z**2 ) / ( 1. - Z ) * LOG( Z/(1.-Z) )
     >                + 3./2. / ( 1.- Z ) - 3. - 2.* Z )
cjfo      
      FMSQ_NUC = TMP*xffac
      RETURN
C                                            ****************************
      ENTRY FMSG_NUC(ZZ)

      Z = ZZ
      IF (ABS(1.-Z).LT.1E-15 .OR. Z.LT.1E-15) THEN
        PRINT *, 'Z out of range in FMSQ_NUC; Z =', Z
        call PAUSE('in fmsq_nuc.f')  
      ENDIF
cjfo   mods for E866 xf
      xffac=1.
      if(ixfx.eq.1)xffac=1./(xa+xb)
      
      TMP = 0.
C                                 Add mu dependent part: essentially Pfg of A-P
      IF (LL) TMP = LOG(AMU/Q) * TR * ( Z**2 + (1. - Z)**2 )
      
C                                 Add DIS to MSbar part: essentially C2g Wilson
      IF (LM) TMP = TMP + TR / 2.*
     >      ( ( Z**2 + (1.- Z )**2 ) * ( LOG( Z/(1.-Z) ) + 1. )
     >        - 6.* Z * ( 1.-Z ) )
C
cjfo
      FMSG_NUC = TMP*xffac
      RETURN
C                        ****************************
      END
