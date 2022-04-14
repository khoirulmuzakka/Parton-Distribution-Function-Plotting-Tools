      DOUBLE COMPLEX FUNCTION DiLog(X)
      DOUBLE COMPLEX X,XX
      DOUBLE PRECISION BER(10),B,P,PI
      DATA PI /3.14159265358979323846/
      DATA P  /1.6449340668482264365/
      DATA BER /0.166666666666666666667,  -0.033333333333333333333,
     >          0.0238095238095238095238, -0.033333333333333333333,
     >          0.075757575757575757576, -0.253113553113553113553,
     >          1.16666666666666666667,   -7.0921568627450980392,
     >          54.971177944862155388,   -529.12424242424242424/
C ***************************************************************************
      NFLAG=0
      XX=X
      LAND=1
      CALL COUNTRY(XX,NFLAG)

      IF (NFLAG.EQ.0) THEN
          XX=1.-X
          LAND=2
          CALL COUNTRY(XX,NFLAG)
      END IF

      IF (NFLAG.EQ.0) THEN
          XX=1./X
          LAND=3
          CALL COUNTRY(XX,NFLAG)
      END IF

      IF (NFLAG.EQ.0) THEN
          XX=1./(1.-X)
          LAND=4
          CALL COUNTRY(XX,NFLAG)
      END IF

      CALL SERIES(XX,DiLog)

      IF (LAND.EQ.1) RETURN

      IF (LAND.EQ.2) THEN
           DiLog=-DiLog+P
           B=ABS(1.-X)
           IF(B.GT.0.D0) DiLog=DiLog-LOG(X)*LOG(1.-X)
      END IF

      IF (LAND.EQ.3) THEN
          DiLog=-DiLog-P-0.5*LOG(-X)*LOG(-X)
      END IF

      IF (LAND.EQ.4) THEN
          DiLog=DiLog+2.*P-LOG(X)*LOG(1.-X)+0.5*(LOG(X-1.))**2
      END IF

      RETURN
      END
  
