      function delta(a,b,c)
      implicit double precision (a-z)
      delta = dsqrt( a**2 + b**2 + c**2 - 2.d0 * ( a*b + a*c + b*c ) )
      return
      end
