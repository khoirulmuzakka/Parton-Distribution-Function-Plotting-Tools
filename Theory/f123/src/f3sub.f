      function f3sub(z,Q2,m1,m2)
      implicit double precision (a-z)
      integer k
      dimension fp(3)
      common /vecax/ Sp, Sm, Rqp, Rqm

      mi = m1
      mo = m2
      qp = Sp
      qm = Sm

      v = sqrt(1.d0-(mi+mo)**2/q2*z/(1.d0-z))
      vb = sqrt(1.d0-(mi-mo)**2/q2*z/(1.d0-z))

      fp(1)=v*vb*(mi**2-mo**2)/q2*2.d0*z*(1.d0-z)
      do 20 k = 2, 3
      Ll = dlog((1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)+v*vb)
     >  /     (1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)-v*vb))
       fp(k)=-(.5d0-z*(1.d0-z)+(mi**2-mo**2)/q2*z*(1.d0-2.d0*z)
     > -(mi**4-mo**4)/q2**2*z**2)*Ll
      mi = m2
      mo = m1
  20  continue

      f3sub = -( fp(1) + fp(2) - fp(3) )
      f3sub = Rqp * f3sub

      return
      end
