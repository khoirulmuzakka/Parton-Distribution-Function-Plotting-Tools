
      subroutine dimu_calculate(xbj,y,ehad,eps,iproc,doLO,Born,Onelop)
c Note that the old emc factor has been removed to avoid confusion and
c errors. 
      implicit none
      integer iproc,imass,Iset
      double precision cmass,xbj,ehad,y,res1,eps
      double precision JR,Q,xmu,xi,xmfac,ffac,dsiglo,dsig
      double precision Born, Onelop,q2,prop
      double precision pbarn, gf,pi
      double precision HMASS,fredLO,xmw,amc,ckm
      dimension ckm(3,3)
      logical*1 doLO
      common /PARAMDIMU/HMASS,amc,pi,xmw,gf,pbarn,ckm

C    ---- Always calculate LO
        cmass=amc
        Q=SQRT(2.0*HMASS*Ehad*XBJ*Y)
        xmu=Q

        call signu2nlo(iproc, 1,
     >   cmass,XBJ,Y,Ehad,XMU,dsiglo,dsig)

        q2=Q**2
        xi=xbj*(1.0d0+cmass**2/q2)
        prop=(1.0d0+q2/xmw**2)

        ffac=xi/xbj*8.0d0*prop/pbarn  !*  convert fred's norm
        xmfac=1/(gf**2 * hmass * ehad/pi)  !*  Max's factor

C   ***** make exp adjustment
        res1 = dsiglo * ffac * xmfac
        Born = res1 * eps 

        if (doLO.eqv..true.) then 
          OneLop = Born
        else

C         NLO calculation >>---->>

C   ***** stefan's NLO code

        call signu2nlo(iproc, 2,
     >   cmass,XBJ,Y,Ehad,XMU,dsiglo,dsig)

        res1 = dsig   * ffac * xmfac

        OneLop = res1 * eps
C                                                                     <<----<<
        EndIf    

      Return
      End



      subroutine sig_charm (xbj,y,enu,iproc,doLO,sig)
      implicit none
      integer iproc
      double precision dsig, dsiglo, sig, Q, xmu
      double precision cmass,xbj,enu,y
      double precision pbarn, gf,pi
      double precision HMASS, xmw,amc,ckm
      dimension ckm(3,3)
      logical*1 doLO
      common /PARAMDIMU/HMASS,amc,pi,xmw,gf,pbarn,ckm
      
        cmass=amc
        Q=SQRT(2.0*HMASS*Enu*XBJ*Y)
        xmu=Q

        if (doLo) then 
            print *, "leading order cross section calculated"
            call signu2nlo(iproc, 1, cmass,XBJ,Y,Enu,XMU,dsiglo,dsig)
            sig = dsiglo
        else 
            call signu2nlo(iproc, 2, cmass,XBJ,Y,Enu,XMU,dsiglo,dsig)
            sig = dsig   
        end if
      return
      end






c
c
c
