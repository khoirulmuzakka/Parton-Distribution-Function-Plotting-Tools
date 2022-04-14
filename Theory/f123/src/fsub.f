      function fsub(x,Q2,i,m1,m2)
      Implicit Double Precision (A-H, O-Z)
      Double Precision m1,m2,m1s,m2s

      if ( i .eq. 1) fsub = f1sub(x,Q2,m1,m2)
      if ( i .eq. 2) fsub = f2sub(x,Q2,m1,m2)
      if ( i .eq. 3) fsub = f3sub(x,Q2,m1,m2)

      fsub = 2.d0 * fsub

      return
      end
