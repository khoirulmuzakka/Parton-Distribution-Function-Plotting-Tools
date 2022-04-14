      FUNCTION TRNGLE (X,Y,Z, IRT)
C   "Triangle Function" of the three sides.
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

      DATA RDOFF / 1E-11 /

      IRT = 0

      AMX = MAX (X*X, Y*Y, Z*Z)

      TMP = X*X + Y*Y + Z*Z - 2.* (X*Y + Y*Z + Z*X)
      ATMP= ABS(TMP)

      IF     (ATMP .LT. AMX * RDOFF) THEN
        TMP = ATMP
        IRT = 1
      ELSEIF (TMP .LT. 0.) THEN
        PRINT '(A, 4(1PE12.3))', 'X,Y,Z, TMP =', X, Y, Z, TMP
        STOP 'Negative argument in TRNGLE function; check for errors!'
      ENDIF

      TRNGLE = SQRT (TMP)

      RETURN
C			****************************
      END
