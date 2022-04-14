       SUBROUTINE Fgen123L(icharge,Mode, XBJ, Q, XMU, iTMC, A,F123L)
C-----------------------------------------------------------------------------
C      This is a front-end for Fgen123. "Fgen123L" simply adds on the "L" piece
C      Program to compute both CC and NC F123
C      25 April 2011: Call Fgen123L and compute "L" term
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123L(4)
      Dimension F123( 3)
      Integer iTMC,A
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123
      common /fred/ xmc,xmb,HMASS,SINW2, XMW, xmz

      call Fgen123(icharge,Mode,xbj,q,xmu,iTMC,A,F123)

c     copy arrays
      F123L(1)=F123(1)
      F123L(2)=F123(2)
      F123L(3)=F123(3)
C----------------------------------------------------------------------
C COMPUTE  FL
C----------------------------------------------------------------------
      rho=Sqrt(1.0d0+(2.0d0*hmass*xbj/Q)**2)  !*** Get Hmass from /fred/ common block
      FL=rho**2*F123(2)- 2.0d0*xbj*F123(1)
      F123L(4)=FL

      return
      end
