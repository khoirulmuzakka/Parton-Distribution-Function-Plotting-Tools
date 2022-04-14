      SUBROUTINE XSFACT (IXsec, Ibsn, FAC0, FAC01, FAC1, N1)
C                                                   -=-=- xsfact

C                 Returns the overall (kinematic-variable-independent) factors
C          factors for the various terms in the QCD-based parton model formula.
C
C          IXsec = 1  real vector boson production
C                  2  continuum lepton-pair production
C	
C           N1 is the # of different types of Next-to-Leading parton processes.
C                                     Ibsn is the ID for the Produced Particle.

C     VBP-Version                                       Vector-Boson Production
C                                                       -      -     -
C       N1 = 2  --- we have either the "Annihilation" or the "Compton" process.
C                                                       Ibsn is the VB code
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      PARAMETER (PI = 3.14159265359)
      CHARACTER MSG*65

      DIMENSION FAC1 (N1)
      COMMON / IOUNIT / NIN, NOUT, NWRT

      DATA IW1 / 0 /

      MSG =  'N1 .NE. 2 in XSFACT for VBP'
      IF (N1 .NE. 2) CALL WARNI (IW1, NWRT, MSG, 'N1', N1, 2, 2, 1)
C                                                ------------------------------
C                                                           Coupling parameters
      ALF = ALFEWK (IBSN)
C                                                ------------------------------
C                                 Section on Zeroth order process:  a + b --> P
c                                                    Overall kinematic factor:
      FK0 = PI                                                           
c                                                   pro-ind
C						    spin - color average - sum
      FA0 = 1.  / 2./ 2.  / 3.
C			    kinematic-variable-indep factor from matrix elments
C					(e.g. coupling constants .. etc)
      FM0 = ALF * 16.* PI
C
      FAC0 = FK0 * FA0 * FM0
C                                                  ----------------------------
C                               For LPP, multiply by appropriate overall factor

      IF (IXsec .EQ. 2) FAC0 = FAC0 * ALF / 3./ PI
C                                                  ----------------------------
C                                   We set FAC1=FAC0 since the formulas used in
C                                   the program proper are normalized this way.
      FAC1(1) = FAC0
      FAC1(2) = FAC0
C                        FAC01 is not used in this zero-mass subtraction scheme
      RETURN
C			****************************
      END

