      function qintegrand(xhat)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      common /counter/ i
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
      common /str/ str
      common /sud/ sud
      common /subq/ subq
      common /delta1/ delta1

c ... partonic Mandelstam s=(p_quark+q)**2
      sh = shat(xhat)

c.......................................................................

C      sb =  Ctq3Pd(1,            3, chi/xhat, dsqrt(Q2f), irt)
       sb  = PDF(iSet, Ihad, IPARTONin, chi/xhat, dsqrt(Q2f), irt)

      qintegrand = 4.d0/3.d0 /xhat * 1.d0/(1.d0-xhat)
     >           * (  sb*hq1(xhat,sh) - xhat*str*sud  )

      if ( i .eq. 2 ) qintegrand = qintegrand
     >                           + 4.d0/3.d0/xhat
     >                           * sb * hq21(xhat,sh)

      if ( i .eq. 3 ) qintegrand = qintegrand
     >                           + 4.d0/3.d0/xhat
     >                           * sb * hq31(xhat,sh)

* .... alternatively subtract convolution part here:
c      qintegrand = qintegrand - subqconv(xhat)*4.d0/3.d0

      return
      end
