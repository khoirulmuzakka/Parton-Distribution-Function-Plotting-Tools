      Subroutine SetVBP
C                                                   -=-=- setvbp

C These comments are included in the lead subprogram to survive forsplit.

C===========================================================================
C GroupName: SetVbp 
C Description: Setup VbpPac; initiate kinematic variables; interface common blks
C ListOfFiles: setvbp parvbp xsfact
C===========================================================================

C     Subroutine SetVBP (Ihd1, Ihd2, Jbsn, Iset, Jord, Ischm, Iscal)

C                Call SetVbp before the first time a VBP calculation is done.
C     Use ParVbp to change individual parameters in subsequent calculations.

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)

C                         force linking DatVbp and initiation of commons
      External DatVBP
C                   Set up EW coupling constants for Hard Xsec calculation
      Call StEwCpl2

      Return
C                    ********************
      End

