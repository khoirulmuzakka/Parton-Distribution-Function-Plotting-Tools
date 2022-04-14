       SUBROUTINE Fgen123(icharge,Mode, XBJ, Q, XMU, iTMC, A, F123)
C-----------------------------------------------------------------------------
C      Program to compute both CC and NC F123
C      05/07/2007  Include Z-Z, and G-Z terms
C
C
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Dimension F123(3)
      Integer iTMC,A
      logical, save :: first = .true.
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123


      if(icharge.eq.0) then !*** Neutral Current  (only photon at this point)
          call Fnc123(icharge,Mode,xbj,q,xmu,F123)

      elseif(icharge.eq. 4) then !*** Neutral Current BOTH GAMMA & Z
       Call Fnc123(icharge,Mode, XBJ, Q,XMU, F123)

      elseif(icharge.eq.+1) then !*** Charged Current (W+)
       Call Fcc123(icharge,Mode, XBJ, Q,XMU, F123)

      elseif(icharge.eq.-1) then !*** Charged Current (W-)
       Call Fcc123(icharge,Mode, XBJ, Q,XMU, F123)

      else
         write(6,*) ' error: icharge =',icharge,' not implemented'
         stop
      endif

c     -------------------------------------
      hmass=0.938 !### proton mass
      if(iTMC>0) then
        call doTMC(iTMC,A,xbj,q,hmass,F123)
        if (first) then
           write(6,*) ' TMC corrections triggered with iTMC', iTMC
           first = .false.
        endif
      endif
c     -------------------------------------

      return
      end

      
C =================================================================
      SUBROUTINE doTMC(iTMC,A,x,q,hmass,F123)
C-----------------------------------------------------------------------------
C      Program to compute both TMC corretoins
C      02/13/2020  FIO: This low-level program pick up both CC and NC
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Integer iTMC, A
      Dimension F123(3), CORR(3)
      Double Precision H(0:2)
      Data H / -3.2874, 1.9274, -2.0701 /  ! From CJ15 Table 

       r=  Sqrt(1.0d0 + (4.0d0 * x**2 * hmass**2/q**2))
       xi= 2.0d0 * x /(1.0d0 + r)
       xmu=(hmass/q)**2

C      TMC Pre-factors: =============================
       facF1 = x/(xi*r)
       facF2 = x**2/(xi**2 * r**3)
       facF3 = x/(xi * r**2)

C      Fred's Pre-factors: [From the ACOT paper using light-cone coords]
       rho= r !*** Notation: 
       fredFacF1 = 1.0d0  !*** Factor of "1/2" is part of helicity definition 
       fredFacF2 = 1.0d0/rho**2 !*** Factor of "x" is part of helicity definition 
       fredFacF3 = 1.0d0/rho

C      These are the Higher-Twist factors: Approximation
C       ... leading order in M^2/Q^2 only
C      GUESS FOR F1: NEED TO VERIFY
       f1tmc = 1.0d0 + 1.0d0 * xmu * x * xi /r * (1-xi)**2
C      Code EQ.61
       f2tmc = 1.0d0 + 6.0d0 * xmu * x * xi /r * (1-xi)**2
C      Code EQ.62
       f3tmc = 1.0d0 -  xmu * x * xi /r *((1-xi)*Log(xi))


      select case (iTMC)
          case (1) !***  Only  Remove ACOT Pre-factor and replace w/ TMC factor
         CORR(1) =  facF1/fredFacF1
         CORR(2) =  facF2/fredFacF2
         CORR(3) =  facF3/fredFacF3

          case (2) !***  Apply TMC Corrections & Remove ACOT Pre-factor 
         CORR(1) =  f1tmc*facF1/fredFacF1
         CORR(2) =  f2tmc*facF2/fredFacF2
         CORR(3) =  f3tmc*facF3/fredFacF3
             
C        case (3) !*** Future; apply mu^2 TMC  & Remove ACOT Pre-factor 
C        case (4) !*** Future: Exact TMC integrals  & Remove ACOT Pre-factor 

         case (5) !*** Use CJ15 TMC from Alberto with A^(1/3) dependence
            cjHT=A**(1.d0/3.d0) * h(0) * x**h(1) * (1.0d0 + h(2)*x) ! CJ15 HT correction
            cj15=(1.0d0+cjHT/(q*q))
         CORR(1) =  cj15
         CORR(2) =  cj15
         CORR(3) =  cj15            
         
          case default
          print*, 'iTMC=',iTMC,' case is not implemented. Exiting !'
         CORR(1) =  1.0d0
         CORR(2) =  1.0d0
         CORR(3) =  1.0d0                
         call exit(-1)
      end select

C      Apply Correction to F123
       F123(1) = F123(1) * CORR(1)
       F123(2) = F123(2) * CORR(2)
       F123(3) = F123(3) * CORR(3)

      return
      end
      
      
