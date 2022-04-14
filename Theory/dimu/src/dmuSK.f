
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C-----------------------------------------------------------------------------
      SUBROUTINE signu2nlo(ilep,iorder,
     >   XMC,XBJ,Y,ENU,XMU,dsiglo,dsig)
C      
C      Front-end to Stefan Kretzer's NLO code
C      FIO: 07 March 2004
C      SK: 02 April 2004 supplied a switch so that code returns LO or NLO
C          cross section (iorder=1, 2)
C      
C-----------------------------------------------------------------------------
       Implicit Double Precision (A-H, O-Z)
       Dimension  flo(3), fnlo(3), f(3), born(3)
       Dimension  CKM(3,3)
       common /partonxi/partonxi
       common /counter/i
       common /PARAMDIMU/HMASS,amc,pi,xmw,gf,pbarn,ckm
       common /global/ q2,q2f,q2r,xlam,xi
       common /setup/ lep
       external gconv, qconv


C----------------------------------------------------
C   ***  INITIALIZE
C----------------------------------------------------
       dsiglo =  0.0d0
       dsig   =  0.0d0
C----------------------------------------------------
C   ***  CONVERT VARIABLES
C----------------------------------------------------
      lep=ilep
      q2r=xmu**2
      q2f=xmu**2
      xxmc=amc
      xxmc2=xxmc**2
      x=xbj
C----------------------------------------------------
C   ***  kinematics
C----------------------------------------------------

      q2 = 2.0d0 * hmass * enu * xbj * y
      xi = xbj*(1.d0+xxmc2/q2)   !***   RESCALING VARIABLE

      if (xi .ge. 1.d0) then
        !PRINT*,"[In DIMUON:] xi > 1 EXITING!",xxmc2,q2,xbj
        RETURN
      endif

      xlam = q2/(q2+xxmc2)        !***  SCALED MASS PARAMETER
      s = 2.d0*Enu*hmass           !***  ENERGY^2
      s = 2.d0*Enu*hmass + hmass**2       !*** Fred Patch:
      prop   = (1.d0+q2/xmw**2)**2   !***  BOSON PROPAGATOR CORRECTION
C----------------------------------------------------
C   ***  CALL ALPHAS OVER PI
C----------------------------------------------------

      alp = dimu_alpi(xmu)/2   !*** alpha-s/(2 pi)

C----------------------------------------------------
C   ***  CALL PDF'S
C----------------------------------------------------

      call DIMU_GETPDF (2,xi,Q2f,UVxi,DVxi,Ubxi,Dbxi,Sqxi,Sbxi,GLxi)
      dxi = dvxi + dbxi
      uxi = uvxi + ubxi

C----------------------------------------------------
C   ***  ISO-SPIN ROTATION
C----------------------------------------------------

      if ( iorder .eq. 2 ) then

      if ( lep .eq. 1) then
      partonxi = sqxi/xi*ckm(2,2)**2 + dxi/xi*ckm(2,1)**2
      elseif ( lep .eq. -1 ) then
      partonxi = sbxi/xi*ckm(2,2)**2 + dbxi/xi*ckm(2,1)**2
      endif
comment: 'parton' in function qconv has to be adjusted to 'partonxi'
      do 10 i = 1, 3
comment: Gottschalk eq.(116) q,g*C=F/2 -> factor 2
      f(i) = (alp*DINTEG(gconv,xi,1.D0,1.D-3)*(ckm(2,2)**2+ckm(2,1)**2)
     >     + alp*partonxi*xnocon(xi)
     >     + alp*dinteg(qconv,xi,1.d0,1.d-3)
     >                                        )*2.d0
  10  continue

      else
      do i = 1, 3
      f(i) = 0.d0
      enddo
      endif        ! if ( iorder .eq. 2 ) then


*+++++++++++++++++++++++++++++++++++++++++++++
      if ( lep .eq. 1) then
*        structure functions F1,2,3: addition of Born terms:
          born(1) = sqxi/xi*ckm(2,2)**2 + dxi/xi*ckm(2,1)**2
          born(2) = 2.d0*sqxi*ckm(2,2)**2 + 2.d0*dxi*ckm(2,1)**2
c         comment: F3(Gottschalk): +s(xi) in nu p collisions
          born(3) = 2.d0*sqxi/xi*ckm(2,2)**2 + 2.d0*dxi/xi*ckm(2,1)**2
      elseif (lep .eq. -1) then
*         structure functions F1,2,3: addition of Born terms:
          born(1) = sbxi/xi*ckm(2,2)**2 + dbxi/xi*ckm(2,1)**2
          born(2) = 2.d0*sbxi*ckm(2,2)**2 + 2.d0*dbxi*ckm(2,1)**2
c         comment: F3(Gottschalk): +s(xi) in nu p collisions
          born(3) = 2.d0*sbxi/xi*ckm(2,2)**2 + 2.d0*dbxi/xi*ckm(2,1)**2
      endif

*     conventional normalizition of F1 = .5*F1(Gottschalk)
      fnlo(1) = f(1)/2.d0 + born(1)
      fnlo(2) = xi*f(2)   + born(2)
      fnlo(3) = f(3)      + born(3)
*......................................................................

*+++++++++++++++++++++++++++++++++++++++++++++
*     cross section calculation (d2 sigma/dx dy) in pbarn:

         if (xi*s .le. (q2+xxmc2)) then
         write (*,*) ' y > 1 !!!'
         return
         endif

      dsig   = gf**2*Enu*hmass/pi/prop*pbarn
     >   * (y**2*xbj*fnlo(1)+(1.d0-y-xbj*y*hmass/(2*Enu))*fnlo(2)+(y-y**2/2.d0)*xbj*fnlo(3))
      dsiglo = gf**2*Enu*hmass/pi/prop*pbarn
     >   * (y**2*xbj*born(1)+(1.d0-y-xbj*y*hmass/(2*Enu))*born(2)+(y-y**2/2.d0)*xbj*born(3))    

c......................................................................
 11   continue

 999  return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function gconv(z)
*......................................................................
*     neutrino - gluon scattering
*     convolution g*Cg
*......................................................................
      implicit double precision (a-z)
      common /global/ q2,q2f,q2r,lam,xi
      external cg
      call DIMU_GETPDF (2, z, Q2f, UV, DV, Ub, Db, Sq, Sb, GL)
      gconv = gl/z**2*cg(xi/z)
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function cg(z)
*......................................................................
*     gluon initiated Wilson coefficient Cg
*......................................................................
      implicit double precision (a-z)
      integer i, counter
      common /global/ q2,q2f,q2r,lam,xi
      common /counter/ i
      external guniv
      if (i. eq. 1) then
      c1 = 3.d0-4.d0*(1.d0-lam)
      c2 = (1.d0-lam)*z/(1.d0-lam*z)-.5d0
      c3 = 2.d0
      c4 = -4.d0
      elseif (i .eq. 2) then
      c1 = 7.d0-18.d0*(1.d0-lam)+12.d0*(1.d0-lam)**2
      c2 = (1.d0-lam)/(1.d0-lam*z)-.5d0
      c3 = 6.d0*lam
      c4 = -12.d0*lam
      else
      c1 = -1.d0+2.d0*(1.d0-lam)
      c2 = .5d0
      c3 = -2.d0*(1.d0-z)
      c4 = 2.d0
      endif
      cg = guniv(z)+c1*z*(1.d0-z)+c2+(1.d0-lam)*z*Llam(z)*(c3+lam*z*c4)
* ... Kramer/Lampe/Spiesberger (incorrect)
*      cg = cg - .5d0*((1.d0-z)**2+z**2)*log(q2r/q2)
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function guniv(z)
*......................................................................
*     universal term of Cg (for all Fi, i=1,2,3):
*......................................................................
      implicit double precision (a-z)
      integer i, counter
      dimension vorz(3),ckm(3,3)
      common /PARAMDIMU/HMASS,amc,pi,xmw,gf,pbarn,ckm
      common /global/ q2,q2f,q2r,lam,xi
      common /counter/ i
      data vorz/ 1.d0, 1.d0, -1.d0/
      external Llam, pqg
      mc2 = amc**2
      guniv = (vorz(i)*Llam(z)+dlog((q2+mc2)/q2f)-1.d0)*pqg(z)
     >      + (2.d0*dlog(1.d0-z)-dlog(1.d0-lam*z)-dlog(z))*pqg(z)
C      write (*,*) z,guniv, Llam(z),pqg(z),lam
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function Llam(z)
*......................................................................
*     Gottschalk function eq.(59)
*......................................................................
      implicit double precision (a-z)
      common /global/ q2,q2f,q2r,lam,xi
      Llam = dlog((1.d0-lam*z)/((1.d0-lam)*z))
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function pqg(z)
*......................................................................
*     LO g->q splitting function:
*......................................................................
      implicit double precision (a-z)
      pqg = .5d0*(z**2+(1.d0-z)**2)
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function xnocon(z)
*......................................................................
*     neutrino - quark scattering
*     contributions without convolution from evaluating/integrating
*     '+'/delta distributions:
*......................................................................
      implicit double precision (a-z)
      integer i
      external func
      dimension ckm(3,3)
      common /global/ q2,q2f,q2r,lam,xi
      common /counter/i
      common /PARAMDIMU/HMASS,amc,pi,xmw,gf,pbarn,ckm

      mc2 = amc**2
      ka = 1.d0/lam*(1.d0-lam)*dlog(1.d0-lam)
      call table1(1.d0,a,b1,b2,b3)
*     contributions from hq1,2,3
      hnocon = -(4.d0+1.d0/2.d0/lam+(pi)**2/3.d0
     >       + ka*(1.d0+3.d0*lam)/2.d0/lam)+2.d0*(dlog(1.d0-z))**2
     >       + 2.d0*dinteg(func,0.d0,z,1.d-4)
     >       + a
     >       + b1*dlog((1.d0-z)/z)
comment:       b2 = 0.
     >       + b3/lam*(1.d0/lam*dlog(1.d0-lam)+1.d0)
*     contributions from Cq1,2,3 = prop.*Pqq(z)+... :
      xnocon = 2.d0*(1.d0+4.d0/3.d0*dlog((1.d0-z)/z))*dlog((q2+mc2)/q2f)
     >      + 4.d0/3.d0*hnocon
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function func(z)
      implicit double precision (a-z)
      common /global/ q2,q2f,q2r,lam,xi
      func = dlog((1.d0-lam*z)/(1.d0-z))
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine table1(z,a,b1,b2,b3)
*......................................................................
*     Gottschalk functions/parameters A;B1,2,3 ; eq.(43); tableI:
*......................................................................
      implicit double precision (a-z)
      integer i
      common /global/ q2,q2f,q2r,lam,xi
      common /counter/i
      a = 0.d0
      b3 = .5d0
      if (i .eq. 1) then
      b1 = 1.d0-4.d0*z+z**2
      b2 = z-z**2
      elseif (i .eq. 2) then
      a = .5d0/lam*(1.d0-lam)*dlog(1.d0-lam)
* ... typo in Gottschalk
      a=2.d0*a
      b1 = 2.d0-2.d0*z**2-2.d0/z
      b2 = 2.d0/z-1.d0-z
      else
      b1 = -1.d0-z**2
      b2  = 1.d0-z
      endif
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      function qconv(z)
*......................................................................
*     convolution q*Cq: remaining convolutions after
*     evaluation/integration of '+'/delta distributions
*......................................................................
      implicit double precision (a-z)
      integer i, lep, nf
      dimension ckm(3,3)
      common /global/ q2,q2f,q2r,lam,xi
      common /counter/i
      common /PARAMDIMU/HMASS,amc,pi,xmw,gf,pbarn,ckm
      common /partonxi/partonxi
      common /setup/ lep
      mc2 = amc**2

      call DIMU_GETPDF (2, z, Q2f, UV, DV, Ub, Db, Sq, Sb, GL)

      d = dv+db
      u = uv+ub


      if ( lep .eq. 1) then
      parton = sq/z*ckm(2,2)**2 + d/z*ckm(2,1)**2
      elseif (lep .eq. -1) then
      parton = sb/z*ckm(2,2)**2 + db/z*ckm(2,1)**2
      endif
      call table1(xi/z,a,b1,b2,b3)
      call table1(1.d0,a1,b11,b21,b31)
      xi2 = xi**2
      z2 = z**2
*     universal terms for i=1,2,3
      quniv = 4.d0/3.d0*((1.d0+xi2/z2)*parton-2.d0*partonxi)/(z-xi)
     >      * dlog((q2+mc2)/q2f) + 4.d0/3.d0*((2.d0*dlog(1.d0-xi/z)
     >      - dlog(1.d0-lam*xi/z))/(z-xi)*((1.d0+xi2/z2)*parton
     >      - 2.d0*xi/z*partonxi)-parton/(z-xi)*(1.d0+xi2/z2)
     >      * dlog(xi/z))
      qconv = quniv + 4.d0/3.d0 * (
     >        (parton*b1-partonxi*b11)/(z-xi)
     >      + parton/(z-lam*xi)*b2
     >      + parton/z*b3*(1.d0-xi/z)/(1.d0-lam*xi/z)**2 )
      return
      end
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
C23456789012345678901234567890123456789012345678901234567890123456789012
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      subroutine DIMU_GETPDF(ictq,x,q2f,xuv,xdv,xub,xdb,xsq,xsb,xgl)
*......................................................................
* ... GRV-like CTEQ routine by M.-H. Reno (received May 2002)
*     ... Stefan made chane so that SetCtq6(iset) is called whenever ictq
*         is changed
* ... and seperated strange from anti-strange
*......................................................................
      implicit real*8 (a-h,o-z)
      logical first
      save isav, first
      data first /.true./


      q = dsqrt(q2f)
      xub =   dimu_pdf(-1,x,q)       * x
      xdb =   dimu_pdf(-2,x,q)       * x
      xuv =   dimu_pdf( 1,x,q) * x -xub
      xdv =   dimu_pdf( 2,x,q) * x -xdb
      xsq =   dimu_pdf( 3,x,q)       * x
      xsb =   dimu_pdf(-3,x,q)       * x
      xgl =   dimu_pdf( 0,x,q)       * x

      return
      end


