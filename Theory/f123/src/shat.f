      function shat(xp)
      implicit double precision (a-z)
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget
      shat = Q2/xp * chi/xb - m1s*xp * xb/chi + m1s-Q2
      return
      end
