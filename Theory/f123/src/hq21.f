      function hq21(xhat,s)
      implicit double precision (a-z)
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
      common /powers/ Q4,Q6,Q8,m14,m16,m18,m24,m26,m28
      common /delta1/ delta1
      common /vecax/ Sp, Sm, Rqp, Rqm

      delta2 = delta(m1s,s,-Q2)
      delta12 = delta1  * delta1
      delta14 = delta12 * delta12
      delta22 = delta2  * delta2
      delta24 = delta22 * delta22
      spp = Q2 + m2s + m1s
      smm = Q2 - m2s - m1s
      spm = Q2 + m2s - m1s
      smp = Q2 - m2s + m1s
      spm2 = spm * spm
      smp2 = smp * smp
      diff = m1s+s+Q2-delta2
      sum  = m1s+s+Q2+delta2
      LN = dlog(diff/sum)

       int2 =
     > ((16.d0*(-2.d0*delta14*Sp*(m2s/(s-m2s)**2+s/(s-m2s)**2+
     > (LN*s*spp)/(delta2*(s-m2s)**2))+
     > 2.d0*m1*m2*Sm*((LN*s*((s-m2s)**3+2.d0*(s-m2s)**2*spm+
     > ((s-m2s)*(-delta12+3.d0*spm2))/2.d0))/
     > (delta2*(s-m2s))-((s-m2s)*(-2.d0*(m1s-m2s+2.d0*Q2)+(s-m2s)+
     > (-delta12+3.d0*smp2)/(2.d0*(s-m2s)))*((s-m2s)+spp))/(2.d0*s)+
     > 2.d0*(delta12-(2.d0*m1s-2.d0*m2s+Q2)*(s-m2s)+
     > (s-m2s)**2-3.d0*Q2*spp))+
     > Sp*(-2.d0*(m1s+m2s)*(s-m2s)**2-9.d0*m2s*spm2-
     > (2.d0*delta12*(delta12+2.d0*m2s*spm))/(s-m2s)+
     > 2.d0*(s-m2s)*(2.d0*delta12+(m1s-5.d0*m2s)*spm)+
     > (s*(-4.d0*delta12*spm+(s-m2s)*(delta12-3.d0*spm2)))/(s-m2s)+
     > ((s-m2s)*(-2.d0*(m1s-m2s+2.d0*Q2)+(s-m2s)+(-delta12+
     > 3.d0*smp2)/(2.d0*(s-m2s)))*spp*
     > ((s-m2s)+spp))/(2.d0*s)+delta12*(-m2s+2.d0*spp)+
     > (LN*s*(2.d0*delta12*(-3.d0*delta12+4.d0*m1s*smp)-
     > 4.d0*(s-m2s)**2*(m1s*(-m1s+m2s)+spm2)-(s-m2s)**3*spp+
     > (s-m2s)*(delta12*(8.d0*m1s-7.d0*spp)+18.d0*m1s*Q2*spp)))/
     > (delta2*(s-m2s)))))/delta24)

       int1 =
     > ((8.d0*(-(delta12*(m2s/(s-m2s)**2+s/(s-m2s)**2+
     > (LN*s*spp)/(delta2*(s-m2s)**2))*
     > (-2.d0*m1*m2*Sm+Sp*spp))+
     > 2.d0*m1*m2*Sm*((s*((s-m2s)+2.d0*spm))/(s-m2s)+
     > (delta12+(-3.d0*m1s+2.d0*m2s+Q2)*(s-m2s)+2.d0*m2s*spm)/(s-m2s)+
     > (LN*s*(3.d0*delta12+(s-m2s)**2-4.d0*m1s*smp+(s-m2s)*
     > (m1s+3.d0*spm)))/(delta2*(s-m2s))+
     > ((m2s+(s-m2s))*((s-m2s)+spp))/(2.d0*s))+
     > Sp*(4.d0*m14+2.d0*m1s*(s-m2s)-spm*(m2s+spm)-
     > ((delta12+2.d0*m2s*spm)*spp)/(s-m2s)-
     > (s*spm*((s-m2s)+2.d0*spp))/(s-m2s)+
     > (((s-m2s)+spp)*(delta12-4.d0*m1s*(s-m2s)+
     > (s-m2s)**2-2.d0*m2s*spp))/(4.d0*s)+
     > (LN*s*(-(s-m2s)**3-4.d0*(s-m2s)**2*spm+(s-m2s)*
     > (4.d0*m1s*m2s-7.d0*spm*spp)+
     > 2.d0*spp*(-delta12-2.d0*spm*spp)))/
     > (2.d0*delta2*(s-m2s)))))/delta22)

      hq21 = (s-m2s) / (8.d0*s) *
     >     ( delta2**2/2.d0/Sp/delta1           * int2
     >     - 2.d0*delta1/(Sp*spp-2.d0*Sm*m1*m2) * int1 )

      return
      end
