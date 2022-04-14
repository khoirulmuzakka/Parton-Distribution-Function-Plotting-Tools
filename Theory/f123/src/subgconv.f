      function subgconv(x)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
C      gl =  Ctq3Pd(1,    0, x, dsqrt(Q2f),irt)
       gl  = PDF(iSet, Ihad, 0, x, dsqrt(Q2f),irt)

      subgconv = 1.d0/x * gl * pqg(xi/x)
      return
      end
