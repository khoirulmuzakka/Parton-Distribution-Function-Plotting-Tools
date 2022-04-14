      function gintegrand(xhat)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s
      external PDF
      Common  / CFRED /  iSet, IPARTONin, IPARTONout,Ihad
      common /counter/ i
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 

C      gl =  xhat * Ctq3Pd (1,    0, xhat, dsqrt(Q2f), irt)
       gl  =  xhat * PDF(iSet, Ihad, 0, xhat, dsqrt(Q2f), irt)

      if ( i .ne. 2 ) gl = gl / xhat
      gintegrand = 1.d0/xhat * fsub(xb/xhat,Q2,i,m1,m2) * gl

      return
      end
