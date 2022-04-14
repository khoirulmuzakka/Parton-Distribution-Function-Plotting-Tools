      SUBROUTINE KineLmt(Rs, p1m, p2m, yMx, ptMx)
C                                                   
C				                             Simple two-particle kinematics
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
C                                                 Max. y1 occurs at p_t = 0:
         CHYM = (Rs**2 + p1m**2 - p2m**2) / 2. / Rs / p1m
         SHYM = SQRT (CHYM ** 2 - 1D0)
         YMX = LOG (CHYM + SHYM)
C							                          Max. p_t = Abs(3-mom):
         ptMx2 = ((Rs**2 +p1m**2 -p2m**2) /2D0 /Rs) **2 - p1m**2
         ptMx = Sqrt(ptMx2) 

      RETURN
C			****************************
      END
C                                                          =-=-= VbpXsc
