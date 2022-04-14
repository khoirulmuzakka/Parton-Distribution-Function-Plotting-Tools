      function that(xp,zp)
      implicit double precision (a-z)
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
      that = ( shat(xp)-m1s+Q2 )
     >     * ( zp - ( 1.d0 + 2.d0*m1s/(shat(xp)-m1s+Q2) ) ) + m1s
      return
      end
