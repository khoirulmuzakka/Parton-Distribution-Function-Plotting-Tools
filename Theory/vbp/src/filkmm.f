      Subroutine FILKMM 
C                                                   -=-=- filkmm
C                       Given the Mixing parameters for NGN generations, 
C                       calculate the KM matrix in the Cabbibo-KM-Maiani-
C                       -Wolfenstein-Chau-Keung...etc scheme.
C                       Used in SetEwk.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (D0=0D0, D1=1D0, D2=2D0, D3=3D0, D4=4D0, D10=1D1)
C      PARAMETER (NGN = 3, NANG= NGN*(NGN-1)/2, NPHS=(NGN-1)*(NGN-2)/2)
      PARAMETER (NGN = 3)
      PARAMETER (NANG= NGN*(NGN-1)/2)
      PARAMETER (NGNm1= (NGN-1))
      PARAMETER (NGNm2= (NGN-2))
      PARAMETER (NPHS= NGNm1*NGNm2/2)
 
      COMMON / KMATRX / VKM (NGN, NGN)
      COMMON / MIXPAR / CMX(NANG), DMX(NPHS)

C     (how to handle complex numbers due to phase factors efficiently,
C      bearing in mind that most applications do not need this info??)

C     For the moment, put the absolute values of the matrix elements in by hand

      DIMENSION AK (3, 3)

      DATA AK(1,1),AK(2,1),AK(3,1)  / 0.9755, 0.220, 0.007 /
      DATA AK(1,2),AK(2,2),AK(3,2)  / 0.220, 0.9744, 0.046 /
      DATA AK(1,3),AK(2,3),AK(3,3)  / 0.012, 0.045, 0.9993 /

      DO 5 I = 1, 3
      DO 6 J = 1, 3
        VKM (I,J) = AK (I,J)
    6 CONTINUE
    5 CONTINUE

      Return
C                        ****************************
      END

C                                                          =-=-= EwCpl0
