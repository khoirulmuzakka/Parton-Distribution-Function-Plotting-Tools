      function fq(i)
      implicit double precision (a-z)
      integer i
      dimension acotlo(3)
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget
      common /delta1/ delta1
      common /acotlo/ acotlo

      fq = 2.d0 /delta1
     >   *  ( (m1s+m2s+Q2)*dlog((m1s+m2s+Q2+delta1)/(m1s+m2s+Q2-delta1))
     >   - 2.d0*delta1  )

      return
      end
