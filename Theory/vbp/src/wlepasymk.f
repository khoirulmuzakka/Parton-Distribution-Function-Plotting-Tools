      subroutine WlepAsymK (Ird, Ist, rs, y, ptmin, asym)
C                                                   -=-=- wlepasymk

      Implicit Double Precision(A-H, O-Z)

c---------------------------------------------------------------------------
c  Lai 6/25/04 correct the incomplete implementation for read in from file
c              and clear up some unused/unnecessary stuff
c  Lai 4/2/97, 9/2/97, 9/16/97
c  add Ird switch for NLO and Resum K-factor
c  The K-factor parametrizations are obtained by RESBOS using CTEQ4M 
c  or read in from kfac.wasy file.
c
c  Ird = 1, LO
c
c  6/25/04 
c  Ird = 2, NLO with average k-factor obtained from CTEQ4M 
c  Ird = 3, Resum with average k-factor obtained from CTEQ4M 
c  Ird = 4, read in k-factor from kfac.wasy 
c---------------------------------------------------------------------------
c  calculates the LO w lepton asymmetry. See Barger and Phillips 
c  "Collider Physics" pgs 254-256. Note: theta-hat on pgs 255, 256 
c  is - (theta-hat) on page 254!  
c
c  J.F. Owens    June 30, 1994 ; modified, wkt 7/2/94
c
c  rs = sqrt(s)  y=lepton rapidity
c

      dimension xa(24),argp(24),argm(24),qa(8),qb(8)
      Common  /WasyK/ Ak(0:6)

c old coef.
c      AkNloMrsa(x)= 1.0681+ .0140426*x+ .0366726*x**2+ .00459567*x**3
c     >    - .00477456*x**4 - .000444555*x**5 + .000204021*x**6
c      AkNloCt4m(x)= 1.06761+ .0226861*x+ .038713*x**2+ .00215193*x**3
c     >    - .00671079*x**4 - .000378377*x**5 + .000240478*x**6
c      AKNLO(x)= 1.06786+ .0183644*x+ .0376928* x**2+ .0033738*x**3
c     >    - .005742675*x**4 - .000411466*x**5 + .0002222495*x**6
c      AkResMrsa(x)= 1.00729+ .00718262*x+ .0378679*x**2+ .00935662*x**3
c     >    - .00352178*x**4 - .000843583*x**5 + .0000245564*x**6
c      AkResCt4m(x)= 1.01084+ .010643*x+ .0410257*x**2+ .00963999*x**3
c     >    - .00615128*x**4 - .00104759*x**5 + .0000935093*x**6
c      AKRES(x)= 1.009065+ .00891281*x+ .0394468*x**2+ .009498305*x**3
c     >    - .00483653*x**4 - .0009455865*x**5 + .00005903285*x**6
c      Ak1NloCt3m(x)= 1.00596+ .011023*x+ .025569*x**2+ .005997*x**3
c     >    - .0017145*x**4 - .00052765*x**5
c      Ak1ResCt3m(x)= 0.92455- .00049911*x+ .022254*x**2+ .011930*x**3
c     >    + .00064675*x**4 - .0011885*x**5 - .00034656*x**6

      Ak1NloCt4m(x)= 1.01788+ .0067692*x+ .027310*x**2+ .0074929*x**3
     >    - .0019677*x**4 - .00067043*x**5
      Ak1ResCt4m(x)= 1.07949- .0041764*x+ .017066*x**2+ .016434*x**3
     >    + .00074073*x**4 - .0015393*x**5 - .00040555*x**6
      Akread(x)=Ak(0)+Ak(1)*x+Ak(2)*x**2+Ak(3)*x**3+Ak(4)*x**4
     >          +Ak(5)*x**5+Ak(6)*x**6

c  ckm factors (using just the Cabibbo angle)
      c2tc=.95
      s2tc=.05
c  w mass
      xmw=80.
      q = xmw
      s=rs**2
c
c  CDF uses pt(muon)>25 x=xmw/(2.*ptmin)
c
      x=xmw / (2.*ptmin)
      w=x+sqrt(x**2-1.)
      xi= log(w)
c
c  In lowest order, the pt cut translates into simple 
c  limits on the xa integral
c
      xamin=xmw/rs*exp(y-xi)
      xamax=xmw/rs*exp(y+xi)
      If(xamax.gt.1D0) xamax=1D0
c
c  6-pt gaussian quadrature routine -- supplied below
c
      call gq11(xamin,xamax,4,xa,argp,ansp)
      do 100 j=1,24
      xb=xmw**2/(s*xa(j))
c
c  my parton distribution calls
c  replace with yours
c
      call dist(Ist, xa(j), Q, qa)
      call dist(Ist, xb,    Q, qb)
      yh=y-.5* log(xa(j)/xb)
      st=1./cosh(yh)
      ct=tanh(yh)
      facp=(1.-ct)**2*st**2
      facm=(1.+ct)**2*st**2
c
c  notation is u,d,s,ubar,dbar,sbar,g,c 
c  for i=1,2,3,4,5,6,7,8
c  modify as needed
c
      argp(j)=c2tc*((qa(1)*qb(2)+qa(8)*qb(3))*facp
     2+(qa(5)*qb(4)+qa(6)*qb(8))*facm)
     3+s2tc*((qa(1)*qb(3)+qa(8)*qb(2))*facp
     4+(qa(6)*qb(4)+qa(5)*qb(8))*facm)
      argm(j)=c2tc*((qa(2)*qb(1)+qa(3)*qb(8))*facm
     2+(qa(4)*qb(5)+qa(8)*qb(6))*facp)
     3+s2tc*((qa(3)*qb(1)+qa(2)*qb(8))*facm
     4+(qa(4)*qb(6)+qa(8)*qb(5))*facp)
      argp(j)=argp(j)/xa(j)
      argm(j)=argm(j)/xa(j)
  100 continue
      call gq11(xamin,xamax,0,xa,argp,ansp)
      call gq11(xamin,xamax,0,xa,argm,ansm)
      If(Ird.eq.2) then
         ansp=ansp*Ak1NloCt4m(y)
         ansm=ansm*Ak1NloCt4m(-y)
      ElseIf(Ird.eq.3) then
         ansp=ansp*Ak1ResCt4m(y)
         ansm=ansm*Ak1ResCt4m(-y)
      ElseIf(Ird.eq.4) then
         ansp=ansp*Akread(y)
         ansm=ansm*Akread(-y)
      Endif
      asym=(ansp-ansm)/(ansp+ansm)
      return
      end
