      function subqconv(x)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
      common /str/ str
C      sb =  Ctq3Pd(1,           3, chi/x, dsqrt(Q2f), irt)
       sb  = PDF(iSet, Ihad,IPARTONin, chi/x ,dsqrt(Q2f), irt)

      subqconv = 1.d0/x * ( (1.d0+x**2)/(1.d0-x)
     >         * ( dlog(Q2f/m1s) - 1.d0-2.d0*dlog(1.d0-x) )
     >         * ( sb - x*str ) )

      return
      end
