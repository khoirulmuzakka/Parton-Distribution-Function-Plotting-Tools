      FUNCTION ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adzint

C===========================================================================
C GroupName: Adzint
C Description: adaptive integration
C ListOfFiles: adzint adzcal adzspl intusz sglint totalz
C===========================================================================
C                                  Authors: Wu-Ki Tung and John C. Collins
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adzint.f,v 1.1 97/12/21 21:19:04 wkt Exp $
C #Log:	Adzint.f,v $
c Revision 1.1  97/12/21  21:19:04  wkt
c Initial revision
c 

C     FUNCTION   ADZINT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZSPL (F, I, IER)
C     SUBROUTINE ADZCAL (F,I)
C     SUBROUTINE SGLINT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOTALZ
C     FUNCTION   INTUSZ (X, FX)
C
C     COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), NUMINT,
C    > ICTA, ICTB, FA, FB, IB
C                        ------------------------

C     Adaptive integration routine which allows the integrand to be 
C     indeterminant at the lower and/or the upper ends of integration. 

C     Can self-adjust to any integrable singularity at the ends and compute 
C     the closest approximant, hence achieve the required accuracy efficiently
C     (provided the switch(s) IACTA (IACTB) are set to 2).
 
C     Input switches for end-treatment:
C        IACTA = 0 :   Use closed lower-end algorithm 
C                1 :   Open lower-end -- use open quadratic approximant
C                2 :   Open lower-end -- use adaptive singular approximant

C        IACTB = 0, 1, 2   (same as above, for the upper end)
 
C                Integral of F(X) from A to B, with error
C                less than ABS(AERR) + ABS(RERR*INTEGRAL)
C                Best estimate of error returned in ERREST.
CError code is IER: 0 :  o.k.
C                1 :  maximum calls to function reached before the 
C                     error criteria are met;
C                2 :  IACTA out of range, set to 1;
C                3 :  IACTB out of range, set to 1.
C                4 :  Error on Limits : B < A ; zero result returned.
C                5 :  Range of integration DX zero or close to roundoff
C                     returns DX * F(A+DX/2)
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
cf2py double precision y
cf2py y = f(y)
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZWRK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZWRK /
      DATA SMLL, Sml / 1E-20, 1E-12 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZINT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZINT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      DDX = B - A
      If (DDX .Le. 0D0) Then
        AdzInt = 0D0
        Ier = 4
        If (DDX .Lt. 0D0) 
     >     Print '(/A/)', 'B < A in AdzInt; check limits!!'
        Return
      ElseIf (DDX .Le. Sml) Then
        AdzInt = F(A + DDX/2) * DDX
        Ier = 5
        Return
      EndIf

      NUMINT = 3
      DX = DDX/ NUMINT
      DO 10  I = 1, NUMINT
          IF (I .EQ. 1)  THEN
             U(1) = A 
             IF (IACTA .EQ. 0) THEN
               FU(1) = F(U(1))
             ELSE 
C                                   For the indeterminant end point, use the
C                                   midpoint as a substitue for the endpoint.
               FA = F(A+DX/2.)
             ENDIF
          ELSE
              U(I) = V(I-1)
              FU(I) = FV(I-1)
          ENDIF

          IF (I .EQ. NUMINT) THEN
             V(I) = B
             IF (IACTB .EQ. 0) THEN
               FV(I) = F(V(I))
             ELSE
               IB = I
               FB = F(B-DX/2.)
             ENDIF
          ELSE
              V(I) = A + DX * I
              FV(I) = F(V(I))
          ENDIF
          CALL ADZCAL(F,I)
   10     CONTINUE
       CALL TOTALZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZSPL(F,I,IER)
   40         CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZINT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

