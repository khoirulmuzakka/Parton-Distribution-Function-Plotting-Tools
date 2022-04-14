      subroutine svnlo(rterm,vterm,Sp,Sm,Rqp,Rqm)
      implicit double precision (a-z)
      INTEGER iErr, iActL, iActU
      dimension rterm(3), vterm(3)
      external helpf
      Common  / ActInt /  AERR, RERR, iActL, iActU
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
      common /delta1/ delta1

      pi2 = dacos(-1.d0)*dacos(-1.d0)
      delta1 = delta(m1s,m2s,-Q2)
      sum  = m1s+m2s+Q2+delta1
      diff = m1s+m2s+Q2-delta1
      ln1  = dlog(sum/diff)

c ... arguments of dilogs as lower bounds for integral evaluation
      lb1 = 2.d0*delta1/sum
      lb2 = -diff/2.d0/delta1

C *** FRED: PATCH IN ADZINT  19 MAY 1999
C *** FRED: ensure integration limits a<b are satisfied: 02 may 2005
C     dilog1 = dinteg(helpf,lb1,0.d0,1.d-4)
C     dilog2 = dinteg(helpf,lb2,0.d0,1.d-4)
      dilog1 =-AdzInt(helpf,0.d0,lb1,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog2 = AdzInt(helpf,lb2,0.d0,AERR,RERR,ErrEst,iErr,iActL,iActU)

c ... soft real contributions
      rterm(1) = 2.d0 + (m1s+m2s+Q2)/delta1 * ( ln1 - dilog1 - dilog2
     >         - 0.5d0 * (dlog(diff/2.d0/delta1))**2 - pi2/6.d0 )
     >         + dlog(Q2*(chi-chit)**2/m2s/xb**2)
     >         * ( (m1s+m2s+Q2)/delta1*ln1 - 2.d0 )
c ... app is m1->0 limit of rterm
c      app = (dlog((Q2+m2s)**2/m2s/m1s)-2.d0)*dlog((Q2+m2s)**2/Q2/m2s)
c     >    +2.d0+dlog((Q2+m2s)**2/m1s/m2s)-pi2/3.d0-0.5d0
c     >    *(dlog((Q2+m2s)**2/m1s/m2s))**2
c      write(6,*) rterm(1), app
      rterm(2) = rterm(1)
      rterm(3) = rterm(1)

c ... virtual contributions
      lb3 = ( delta1 - (Q2-m1s+m2s) ) / 2.d0/delta1
      lb4 = ( delta1 +  Q2+m1s-m2s  ) / 2.d0/delta1
      lb5 = ( delta1 +  Q2-m1s+m2s  ) / 2.d0/delta1
      lb6 = ( delta1 - (Q2+m1s-m2s) ) / 2.d0/delta1

C23456789012345678901234567890123456789012345678901234567890123456789012
C *** FRED: PATCH IN ADZINT  19 MAY 1999
C *** FRED: ensure integration limits a<b are satisfied: 02 may 2005
C      dilog3 = dinteg(helpf,lb3,0.d0,1.d-4)
C      dilog4 = dinteg(helpf,lb4,0.d0,1.d-4)
C      dilog5 = dinteg(helpf,lb5,0.d0,1.d-4)
C      dilog6 = dinteg(helpf,lb6,0.d0,1.d-4)

      dilog3 =-AdzInt(helpf,0.d0,lb3,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog4 =-AdzInt(helpf,0.d0,lb4,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog5 =-AdzInt(helpf,0.d0,lb5,AERR,RERR,ErrEst,iErr,iActL,iActU)
      dilog6 =-AdzInt(helpf,0.d0,lb6,AERR,RERR,ErrEst,iErr,iActL,iActU)




      lns1 = ( dlog( ( delta1 - (Q2-m1s+m2s) ) / 2.d0/Q2 ) )**2
      lns2 = ( dlog( ( delta1 +  Q2+m1s-m2s  ) / 2.d0/Q2 ) )**2
      lns3 = ( dlog( ( delta1 +  Q2-m1s+m2s  ) / 2.d0/Q2 ) )**2
      lns4 = ( dlog( ( delta1 - (Q2+m1s-m2s) ) / 2.d0/Q2 ) )**2
      c0gl5 = 1/delta1 * ln1*( (Q2+m1s+m2s)*(dlog(Q2/delta1)+1.d0 )
     >      + delta1**2/2.d0/Q2 )
     >      + (Q2+m1s+m2s)/delta1 * ( 0.5d0 * ( lns1-lns2-lns3+lns4 )
     >      - dilog3 + dilog4 + dilog5 - dilog6 )
     >      - 0.5d0 * (m1s-m2s)/Q2 * dlog(m1s/m2s) + dlog(m1s*m2s/Q2/Q2)
     >      - 4.d0
c ... app is m1->0 limit of c0gl5
c      app = dlog((Q2+m2s)/m1s)*(0.5d0-dlog((Q2+m2s)/Q2))
c     >    + 0.5d0*(dlog((Q2+m2s)/m1s))**2+dlog((Q2+m2s)/m2s)
c     >    *(0.5d0+m2s/Q2-dlog((Q2+m2s)/Q2))+0.5d0
c     >    *(dlog((Q2+m2s)/m2s))**2+2.d0*dlog((Q2+m2s)/Q2)
c     >    +2.d0*dinteg(helpf,Q2/(Q2+m2s),0.d0,1.d-3)-4.d0
c      write(6,*) app, c0gl5

      c0gr5 = 2.d0*m1*m2/delta1 * ln1

      c0p1l5 = -1.d0/Q2*( (Q2-m1s+m2s)/delta1 * ln1
     >       +             dlog(m1s/m2s) )


      c0p1r5 = -1.d0/Q2*( (Q2+m1s-m2s)/delta1 * ln1
     >       -             dlog(m1s/m2s) )


      vterm(1) = c0gl5
     >         + ( Sm*(m1s+m2s+Q2) - 2.d0*Sp*m1*m2 )
     >         / ( Sp*(m1s+m2s+Q2) - 2.d0*Sm*m1*m2 ) * c0gr5

      vterm(2) = c0gl5 + (m1s * c0p1r5 + m2s * c0p1l5)/2.d0
     >         + Sm/Sp * ( c0gr5 + m1*m2/2.d0*(c0p1l5+c0p1r5) )

C ***  PATCH: FIO 5/19/99: IF Rqp=0, F3 SHOULD BE ZERO
      vterm(3) = 0.0
      IF(Rqp.NE.0.0)   vterm(3) = c0gl5 + Rqm/Rqp * c0gr5

      return
      end
