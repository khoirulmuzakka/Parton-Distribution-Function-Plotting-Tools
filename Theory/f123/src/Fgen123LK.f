       SUBROUTINE Fgen123LK(idata,idNuc,
     & icharge,Mode, XBJ, Q, XMU, iTMC, A, F123L)

C-----------------------------------------------------------------------------
C      This is a front-end for Fgen123L
C      On first call of idata point number, it computes and stores the K-factor
C      25 April 2011:
C
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      PARAMETER (ndata=10000,n123L=4,nMode=3,nNuclei = 2,
     &nTotal=ndata*n123L*nMode*nNuclei)  !*** Number of data points for K-factor array
      Dimension XKFACTOR(ndata,n123L,nMode, nNuclei)  !*** ndata points,  (F123L)=(1,2,3,4), Mode= (F,Fc,Fb)=(1,2,3)
      data XKFACTOR /nTotal*0.d0/
      save XKFACTOR

      integer iTMC,A
      Dimension F123L(4),F123LLO(4)
      Common /Ischeme/ Isch, Iset, Iflg, Ihad  !*** pass info out to Fnc123 and Fcc123
      common /fred/ xmc,xmb,HMASS,SINW2, XMW, xmz

C-----------------------------------------------------------------------------
c    if idata>ndata, increase k-factor table
C-----------------------------------------------------------------------------
      if(idata.gt.ndata)  then  !**** over-ride and use full calculation:
         write(6,*) ' Error: idata =',idata,' > ',ndata
         write(6,*) ' Increase ndata '
         stop
      endif

C-----------------------------------------------------------------------------
c    if idata=0 skip k-factor table and use full calculation
C-----------------------------------------------------------------------------
      if(idata.eq.0)  then  !**** over-ride and use full calculation:
         call Fgen123L(icharge,Mode,xbj,q,xmu,iTMC,A,F123L)
         return
      endif

C-----------------------------------------------------------------------------
c     FIRST TIME THROUGH: FILL K-FACTOR
C-----------------------------------------------------------------------------
      if(XKFACTOR(idata,1,mode,idnuc).eq.0.d0)  then  !***  FIRST TIME THROUGH: FILL K-FACTOR  ===
         call Fgen123L(icharge,Mode,xbj,q,xmu,iTMC,A,F123L)
         IschORIG=Isch
         Isch=5  !*** Massive LO Calculation
         Call Fgen123L(icharge,Mode,XBJ,Q,XMU,iTMC,A,F123Llo)
         Isch=IschORIG  !*** Reset Ischeme
C     Generate K-Factor
         do i=1,4
            if(F123Llo(i).eq.0.0d0) then
            XKFACTOR(Idata,i,MODE, idnuc)=1.0d0  !**** Default if denom is zero
            else
            XKFACTOR(Idata,i,MODE, idnuc)=F123l(i)/F123llo(i)
            endif
         enddo
c
C-----------------------------------------------------------------------------
         else  !***  NOT FIRST TIME THROUGH: USE K-FACTOR ======================
            IschORIG=Isch
            Isch=5  !*** Massive LO Calculation
            Call Fgen123L(icharge,Mode, XBJ, Q, XMU, iTMC, A,F123Llo)
            Isch=IschORIG  !*** Reset Ischeme
C     Use K-Factor
            do i=1,4
               F123L(i)=F123Llo(i) * XKFACTOR(Idata,i,mode, idnuc)
            enddo
         endif
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      return
      end
