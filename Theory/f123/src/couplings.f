      subroutine couplings(a1,a2,v1,v2,Sp,Sm,Rqp,Rqm)
      implicit double precision (a-z)
      Sp  =  v1*v2 + a1*a2
      Sm  =  v1*v2 - a1*a2
      Rqp = (a1*v2 + a2*v1) / 2.d0
      Rqm = (a1*v2 - a2*v1) / 2.d0
      return
      end
