      function f123Lext (ic, md, xin, qin, ischin, xmuin, useKfac,
     > iTMC, A, idatain, idNucin, res)

      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      double precision f123Lext
      double precision, dimension (4)::res
      logical usekfac
      integer ic,md,ischin,iTMC,A
      Dimension F123L(4)
      Dimension F123(3)
c Intergace common blocks      
      Common /PREC/ Aerrin, rerrin, iActlin, iactuin
      Common /INTPDF/ isetin, iflgin, ihadin
      Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &ftarget 

c end interface common blocks
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s
      common /fred/ xmc,xmb,HMASS,SINW2, XMW, xmz

c Copied from FgentTest don't know if I have to keep it
      Common  / ActInt /  AERR, RERR, iActL, iActU

      Data icharge, Isch, Iset, Iflg, Ihad, Mode, X, Q
     >   / 4, 1,1,0,1,1,0.1d0,10.d0 /

      HMASS=hmassin
      XMC=xmcin
      XMB=xmbin
      SINW2=sinw2in
      XMW=xmwin
      XMZ=xmzin

      AERR=aerrin
      RERR=RERRin  !*** Setting this too small will yield unstable results ***
      iactl=iActlin ! 2
      iactu=iactu ! 2

c Begin F.'s modif      
      icharge=ic
      Mode = md
      xbj = xin
      q = qin
      Isch = ischin
      xmu = xmuin

      if (.not.usekfac) then
        call Fgen123(icharge,Mode,xbj,q,xmu,iTMC,A,F123)
c       copy arrays
        res(1)=F123(1)
        res(2)=F123(2)
        res(3)=F123(3)
        rho=Sqrt(1.0d0+(2.0d0*hmass*xbj/Q)**2)  !*** Get Hmass from /fred/ common block
        FL=rho**2*F123(2)- 2.0d0*xbj*F123(1)
        res(4)=FL
      else 
        idata=idatain+1
        idNuc=idNucin+1
        call Fgen123LK(idata,idNuc,icharge,Mode, xbj, q, xmu, iTMC,
     &                  A, F123L)
        res(1)=F123L(1)
        res(2)=F123L(2)
        res(3)=F123L(3)
        res(4)=F123L(4)
      endif



      f123Lext=1

      return

      End
