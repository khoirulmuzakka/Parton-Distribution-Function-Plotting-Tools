      subroutine init(Fac0, Fac01, Fac1)
      implicit none
      double precision Fac0, Fac01, Fac1(2)
      integer Iscl, Ibsn, Iset1, Iset2, Jtgt, Jbem, IR, IXsc
      double precision AyXf
      integer :: N1 = 2
      double precision ThVbp
      integer LVbp
      double precision AlfQ, Qalf, aMcbt
      integer Isch, IorHrd, IorQcd, Lqcd, Iset
      COMMON
     >  /ThpVbp/ ThVbp(3), LVbp(3)
     >  /FitHrd/ AlfQ,Qalf,aMcbt(3),ISch,IorHrd,IorQcd,Lqcd(4),Iset

      Call SetVbp

c column 4 line 18 of fitinpnuc -- ThVbp(1) set in initYAML
      Iscl = Nint(ThVbp(1))
c in nCTEQ15 Ibsn 1 for all DY data
      Ibsn = 1
c the pdf sets have a different meaning here, depend on the user implementation of VBP_pdf_
c in examples/c++/example.cpp: 0 ... proton, 1 ... nucleus
      Iset1 = 0
      Iset2 = 1
c target specified in M04lis_nuc, always proton in this case
      Jtgt = 1 
c beam specified in M04lis_nuc, always nucleus in this case (A, Z) should be specified in the exp data files
      Jbem = 682 
c do the call
      Call SetVbpPar_nuc(Jtgt,Jbem,Ibsn,Iset1,Iset2,IorHrd,Isch,Iscl)
c set the vbp observable
c JSFN = 52 (which is the case for all DY data in nCTEQ15) implies AyXf = 0
      AyXf = 0
      Call ParVbp(1, 'Ixfx', AyXf, IR)

C                                   Kinematics independent factors; N1 is set as 2
c JSFN = 52 implies IXsc = 0
      IXsc = 2
      Call XsFact (IXsc, Ibsn, Fac0, Fac01, Fac1, N1)
      
c      Call lamcwz
      end

      subroutine calculate(Fac0, Fac01, Fac1, RS, X1, Q2, Rerr, doLO, Born, Onelop)
      implicit none
      double precision Fac0, Fac01, Fac1(2), RS, X1, Q2, Rerr, Born, Onelop
      double precision X2, TAU, Y, Q
      double precision yMx, ptMx
      double precision :: xMs = 2d0
      double precision xsec0, xsec(2), Scle
      integer Ir
      double precision Ord
      double precision dsigA1, dsigA2, dsigA1nlo, dsigA2nlo
      double precision VbpDy1_nuc
      logical*1 doLO
c........................................................................
c ... compute other variables from two independent variables
      X2 = Q2/(RS*RS*X1)   ! x2
      TAU = dsqrt(X1*X2)
      Y = 0.5d0 * dlog(X1/X2) ! rapidity y
      Q   = TAU * RS

c........................................................................
c ... kinematic variables
      Call SetVbpVar_nuc (RS, Q, Q, Rerr)
      Call KineLmt(RS*RS, Q, xMs, yMx, ptMx)

      Xsec0 = 0.38937966D6 ! X-section in nb - GeV**2
c        Elseif (IXsc .Eq. 2) Then
      Scle = 2.0 * TAU**2 ! For LPP, his leads to  Q**3 *d sig/dy d Q

      Xsec(1) = Xsec0 * SCLE * FAC0
      Xsec(2) = Xsec0 * SCLE * FAC1(1)

c........................................................................
c ... Calculation of cross section in leading order
      Ord = 1 ! Set hard-cross-section calculation to LO
      Call ParVbp_nuc(1, 'Iord', Ord, Ir)
      call vbp_mysetpdf(0)
      dsigA1 = VbpDy1_nuc(y) ! up to normalization which cancels in the ratio
      call vbp_mysetpdf(1)
      dsigA2 = VbpDy1_nuc(y) ! up to normalization which cancels in the ratio
      Born = dsigA1/dsigA2
      if (doLO.eqv..true.) then
        OneLop = Born
      else

c........................................................................
c ... NLO
          Ord = 2 ! Set hard-cross-section calculation to NLO
          Call ParVbp_nuc(1, 'Iord', Ord, Ir)
          call vbp_mysetpdf(0)
          dsigA1nlo = VbpDy1_nuc(y) ! up to normalization which cancels in the ratio
          call vbp_mysetpdf(1)
          dsigA2nlo = VbpDy1_nuc(y) ! up to normalization which cancels in the ratio
          Onelop = dsigA1nlo/dsigA2nlo
      endif
      end
