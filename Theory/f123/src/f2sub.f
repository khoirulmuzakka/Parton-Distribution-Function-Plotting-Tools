      function f2sub(z,Q2,m1,m2)
      implicit double precision (a-z)
      integer  k
      dimension fp(3)
      common /vecax/ Sp, Sm, Rqp, Rqm

      mi = m1
      mo = m2
c      qp = 2.d0
c      qm = 0.d0
      qp = Sp
      qm = Sm

      v = sqrt(1.d0-(mi+mo)**2/q2*z/(1.d0-z))
      vb = sqrt(1.d0-(mi-mo)**2/q2*z/(1.d0-z))

      fp(1)=qp*z*(v*vb*(-.5d0+4.d0*z*(1.d0-z)-((mi**2+mo**2)/q2-
     > (mi**2-mo**2)**2/q2**2)*z*(1.d0-z)))
      do 20 k = 2, 3
      Ll = dlog((1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)+v*vb)
     >  /     (1.d0+(mi**2-mo**2)/q2*z/(1.d0-z)-v*vb))
       fp(k)=qp*z*.25d0*(1.d0-2.d0*z*(1.d0-z)+mi**2/q2
     > *(1.d0+8.d0*z-18.d0*z**2)+mo**2/q2*(1.d0-4.d0*z+6.d0*z**2)-qm/qp
     > * 2.d0*mi*mo/q2
     > -(mi**4+mo**4)/q2**2*2.d0*z*(1.d0-3.d0*z)+mi**2*mo**2/q2**2
     > *4.d0*z
     > *(1.d0-5.d0*z)+(mi**6-mi**4*mo**2-mi**2*mo**4+mo**6)/q2**3
     > *2.d0*z**2)*Ll
      mi = m2
      mo = m1
  20  continue

      f2sub = fp(1) + fp(2) + fp(3)

      return
      end
