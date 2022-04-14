      SUBROUTINE UpC (A, La, UpA)
C                                                   -=-=- upc

C  5/29/94 WKT

C  This is a variation of the two old routines UC(A) and UpCase(A). Here
C  the converted value is return to the new variable UpA, rather than
C  the input dummy variable A, to avoid bombing the program with error
C  "access violation, reason mask=04" when the routine is given a CONSTANT
C  actual argument (rather than a character variable argument).

C  The inconvenience of this new version is that the calling program must
C  declare the length of UpA explicitly; and it better be .Ge. Len(A).

C  We trim any excess trailing characters in UpA and replace with spaces.
C  To be exact, the returned value should be considered as UpA(1:La).

      CHARACTER A*(*), UpA*(*), C*(1)
      INTEGER I, La, Ld

      La = Len(A)
      Lb = Len(UpA)

      If (Lb .Lt. La) Stop 'UpCase conversion length mismatch!'

      Ld = ICHAR('A')-ICHAR('a')

      DO 1 I = 1, Lb

        If (I .Le. La) Then
         c = A(I:I)
         IF ( LGE(C, 'a') .AND. LLE(C, 'z') ) THEN

           UpA (I:I) = CHAR(Ichar(c) + ld)
         Else
           UpA (I:I) = C

         ENDIF
        Else
         UpA (I:I) = ' '
        Endif

 1    CONTINUE
      
      RETURN
C               *************************
      END

C                                                          =-=-= Nurcpe
