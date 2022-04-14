      SUBROUTINE WARNI (IWRN, NWRT, MSG, NMVAR, IVAB,
C                                                   -=-=- warni
     >                  IMIN, IMAX, IACT)

C     t++++++++++++++++++++++++++++++++     Routines to handle Warnings
C                                                  Integer version
      CHARACTER*(*) MSG, NMVAR

      Save Iw

      Data Nmax / 100 /

      IW = IWRN
      IV = IVAB
      
      IF  (IW .EQ. 0) THEN
         PRINT '(1X,A/1X, 2A,I10 /A,I4)', MSG, NMVAR, ' = ', IV,
     >         ' For all warning messages, check file unit #', NWRT
         IF (IACT .EQ. 1) THEN
         PRINT       '(A/2I10)', ' The limits are: ', IMIN, IMAX
         WRITE (NWRT,'(A/2I10)') ' The limits are: ', IMIN, IMAX
         ENDIF
      ENDIF

      If (Iw .LT. Nmax) Then
         WRITE (NWRT,'(1X,A/1X,2A, I10)') MSG, NMVAR, ' = ', IV
      Elseif (Iw .Eq. Nmax) Then
         Print '(/A/)', '!!! Severe Warning, Too many errors !!!'
         Print '(/A/)', '    !!! Check The Error File !!!'
         Write (Nwrt, '(//A//)')
     >     'Too many warnings, Message suppressed !!'
      Endif

      IWRN = IW + 1

      RETURN
C               *************************
      END

