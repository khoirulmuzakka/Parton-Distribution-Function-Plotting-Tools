      SUBROUTINE HADNLOQ123K(Isch,Ihad,XBJin,Qin,F1Min,F2Min,
     >  GLQin1,GRQin1,GLQin2,GRQin2,XMUin,HMASS,IPARTin,ISETin,FHAD123Q)
C-----------------------------------------------------------------------------
C      Computes Hadron Helicities
C      Calls XINTNLO to do convolution
C
C      HEAVY QUARK INITIATED CONTRIBUTIONS TO DEEP INELASTIC STRUCTURE FUNCTIONS.
C      S. Kretzer, I. Schienbein. Phys.Rev.D58:094035,1998  hep-ph/9805233
C-----------------------------------------------------------------------------
C-----------------------------------------------------------------------------
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      Double Precision nlo,nloq,nlog
      dimension acotlo(3), rterm(3), vterm(3)
      DIMENSION FHAD123Q(3)
      PARAMETER (PI=3.141592653589793)
C----------------------------------------------------------------------
      external subgconv, subqconv, gintegrand, qintegrand
      external PDF
C----------------------------------------------------------------------
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihadout
      Common  / ActInt /  AERR, RERR, iActL, iActU
      common /str/ str
      common /counter/ i
      common /delta1/ delta1
      common /sud/ sud
      common /subq/ subq
      common /acotlo/ acotlo
      common /vecax/ Sp, Sm, Rqp, Rqm
      integer Ftarget, Ftarget2
      Common /VARS/ hmassin, xmcin,xmbin,sinw2in, xmwin, xmzin,
     &Ftarget 
      
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget2

C----------------------------------------------------------------------
      DELTA(A,B,C) = SQRT(A**2 + B**2 + C**2 - 2.0*(A*B+B*C+C*A))
C----------------------------------------------------------------------
C----------------------------------------------------------------------
C     INITIALIZATION
C----------------------------------------------------------------------
      DO I=1,3,1
         FHAD123Q(I)=0.0D0
      ENDDO

* ... colour factor
      cf = 4.d0 / 3.d0

      iset=isetIN
      ipart=ipartIN
      IPARTONin=ipart
      IPARTONout=99  !*** Not used
      ihadout=ihad
      Ftarget2=Ftarget

      xb=XBJin
      Q2= Qin*Qin
      m1= F1Min
      m2= F2Min
      m1s=m1*m1
      m2s=m2*m2

      Q2f = XMUin*XMUin
      Q2r = Q2f
C----------------------------------------------------------------------
C     TEST
C     PUT TARGET MASS CORRECTION HERE : 5/7/07
C----------------------------------------------------------------------
        xbj=xb
        hmass=hmassin
        q=sqrt(q2)
        Discr = 1.0 + (4.0*HMASS**2*XBJ**2)/Q**2
        if (Ftarget.eq.1) then
          ETA   = 2.0/(1.0 + Sqrt(Discr))
        else 
          ETA = 1.0
        endif
C        write(6,*) ' eta = ',eta
        xb=xb*eta
C----------------------------------------------------------------------
C----------------------------------------------------------------------


      xi = xb * (1.d0+m2s/Q2)
      chi = xb/(2.d0*Q2) * ( (Q2-m1s+m2s) + delta(m1s,m2s,-Q2) )
      chit= xb/(2.d0*Q2) * ( (Q2-m1s+m2s) - delta(m1s,m2s,-Q2) )
      axb = ( 1.d0 + (m1+m2)**2/Q2 ) * xb
C----------------------------------------------------------------------
C     TEST
C----------------------------------------------------------------------
      IF(CHI.GE.1.0D0) RETURN


C----------------------------------------------------------------------
C*... vector and axial couplings:  General CURRENT PROCESS:
C      v1 = GLQin + GRQin
C      v2 = GLQin + GRQin
C      a1 = GLQin - GRQin
C      a2 = GLQin - GRQin
C----------------------------------------------------------------------
C----------------------------------------------------------------------
*... vector and axial couplings:  General CURRENT PROCESS:
      v1 = GLQin1 + GRQin1
      v2 = GLQin2 + GRQin2
      a1 = GLQin1 - GRQin1
      a2 = GLQin2 - GRQin2
C----------------------------------------------------------------------

      call couplings(a1,a2,v1,v2,Sp,Sm,Rqp,Rqm)

* ... LO amplitudes including kinematic mass effects:
      acotlo(1) = ( Sp*(Q2+m1s+m2s) - 2.d0*Sm*m1*m2 )
     >          / 2.d0 / delta(m1s,m2s,-Q2)
      acotlo(2) = Sp * xb * delta(m1s,m2s,-Q2)/Q2
      acotlo(3) = 2.d0 * Rqp

C      write(6,*) '(acotlo(i),i=1,3,1),xi,chi,chit,axb',
C     >   (acotlo(i),i=1,3,1),xi,chi,chit,axb


* ... real soft and virtual contributions from NLO quark graphs
      call svnlo(rterm,vterm,Sp,Sm,Rqp,Rqm)

C      str = Ctq3Pd(iset, ipart, chi, dsqrt(Q2f),irt)
       str = PDF(iset, Ihad, ipart, chi, dsqrt(Q2f),irt)

C       write(6,*) str,chi, ipart,iset, dsqrt(Q2f)

C----------------------------------------------------
* ... for structure function F_i
      do 10 i = 1, 3

      sud = fq(i)

* ... NLO contributions: gluon AND quark initiated
C      nlog = dinteg(gintegrand,axb,1.d0,1.d-3)
       nlog = 0.0
C **** NEED PROTECT THE UPPER ENDPOINT: iActU= 1 OR 2, NOT 0
C      nloq = dinteg(qintegrand,chi,1.d0,1.d-4)
       nloq= 0.0d0
       if(chi.le.1.d0) then
        nloq = AdzInt(qintegrand,chi,1.d0,
     >              AERR, RERR, ErrEst, iErr, iActL, iActU)
          else
          Write(6,*) " error: integration limits >1 in HADNLOQ123K"
       endif

* ... sum of quark/gluon convolutions and soft/virtual terms
      nlo = nlog + acotlo(i) * ( nloq
     >    + cf * (rterm(i)+vterm(i)) * str
* ... including the Sudakov log from the 1/(1-xi)+ distibution
     >    + cf * sud * str * dlog(1.d0-chi) )

      nloq = nlo - nlog

comment: F_i's normalized to +acotlo(i)*s(chi)+O(alpha_s) for all i
      alphas= alQCD(Sqrt(Q2r),iset)
      nloq  = nloq  * alphas/2.d0/PI

      FHAD123Q(i)=   nloq

C----------------------------------------------------
C  *** ADJUST NORMALIZATION CONVENTIONS TO MATCH OURS:
      IF(I.EQ.2)  FHAD123Q(i) =  FHAD123Q(i)/(1.0+M2S/Q2)
C----------------------------------------------------

10    continue

C----------------------------------------------------

C      write(6,*)  'FHAD123Q', (FHAD123Q(i),i=1,3,1)

      return
      end
