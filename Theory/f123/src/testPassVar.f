      function TestFunction(aa,bb,q)
      integer aa,bb
      double precision q
      double precision TestFunction
      integer ii,jj
      double precision whynot
      common /ABC/ ii,jj, whynot

      write(*,*)'debug',aa,bb,q, ii,jj,whynot

c      TestFunction(1) = aa*q
c      TestFunction(2) = bb*q
c      TestFunction(3) =0d0
c      TestFunction(4) = 0d0

      TestFunction = aa*bb*q
      q = aa

      end
