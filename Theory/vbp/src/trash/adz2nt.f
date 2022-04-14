      FUNCTION ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C                                                   -=-=- adz2nt
 
C===========================================================================
C GroupName: Adz2nt
C Description: second copy of adzint
C ListOfFiles: adz2nt adz2pl adz2al int2sz sgl2nt tot2lz
C=========================================================================== 
C #Header: /Net/cteq06/users/wkt/1hep/1utl/RCS/Adz2nt.f,v 1.1 97/12/21 21:19:00 wkt Exp $
C #Log:	Adz2nt.f,v $
c Revision 1.1  97/12/21  21:19:00  wkt
c Initial revision
c 

C List of GLOBAL Symbols

C     FUNCTION   ADZ2NT (F, A, B, AERR, RERR, ERREST, IER, IACTA, IACTB)
C     SUBROUTINE ADZ2PL (F, I, IER)
C     SUBROUTINE ADZ2AL (F,I)
C     SUBROUTINE SGL2NT (IACT, F1, F2, F3, DX, FINT, ESTER)
C     SUBROUTINE TOT2LZ
C     FUNCTION   INT2SZ (X, FX)
C
C     COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
C    > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
C    > ICTA, ICTB, NUMINT, IB
C                   ------------------------------------

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      EXTERNAL F
      PARAMETER (MAXINT = 1000)
C
C                   Work space:
      COMMON / ADZ2RK / U(MAXINT), V(MAXINT), FU(MAXINT), ERS, RES, 
     > FW(MAXINT), ERR(MAXINT), RESULT(MAXINT), FV(MAXINT), FA, FB,
     > ICTA, ICTB, NUMINT, IB

      SAVE / ADZ2RK /
      DATA SMLL / 1E-20 /
     
      IER = 0
      IF (AERR.LE.SMLL .AND. RERR.LE.SMLL)
     1 STOP 'Both Aerr and Rerr are zero in ADZ2NT!'
        
      IF (IACTA.LT.0 .OR. IACTA.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTA =', IACTA, ' IACTA set for regular open-end option.'
        IACTA = 1
        IER = 2
      ENDIF 
      IF (IACTB.LT.0 .OR. IACTB.GT.2) THEN
        PRINT '(A, I4/ A)', ' Illegal value of IACT in ADZ2NT call', 
     >  'IACTB =', IACTB, ' IACTB set for regular open-end option.'
        IACTB = 1
        IER = 3
      ENDIF
      ICTA = IACTA
      ICTB = IACTB
 
      NUMINT = 3
      DX = (B-A)/ NUMINT
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
          CALL ADZ2AL(F,I)
   10     CONTINUE
       CALL TOT2LZ
C                                                   Adaptive procedure:
   30     TARGET = ABS(AERR) + ABS(RERR * RES)
          IF (ERS .GT. TARGET)  THEN
              NUMOLD = NUMINT
              DO 40, I = 1, NUMINT
                  IF (ERR(I)*NUMOLD .GT. TARGET) CALL ADZ2PL(F,I,IER)
   40             CONTINUE
              IF (IER.EQ.0 .AND. NUMINT.NE.NUMOLD)  GOTO 30
              ENDIF
      ADZ2NT = RES
      ERREST = ERS
      RETURN
C                        ****************************
      END

