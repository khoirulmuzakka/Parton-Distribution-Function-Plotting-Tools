       SUBROUTINE TESTDILOG()
C-----------------------------------------------------------------------------
C      SIMPLE PROGRAM TO TEST DiLog
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       DOUBLE COMPLEX X1,TEST,DiLog
       EXTERNAL DiLog

        DATA X1 /(-1.0,0.0)/

1      CONTINUE
       WRITE(6,*)          X1
       WRITE(6,*) ' ENTER: X1'
       READ (5,*)          X1


       TEST =DiLog(X1)
       WRITE(6,*) TEST
       GOTO 1

       END
