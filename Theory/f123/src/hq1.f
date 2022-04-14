      function hq1(xhat,s)
      implicit double precision (a-z)
      dimension acotlo(3)
      integer Ftarget
      common /kinematic/ xb, Q2, m1, m2, Q2f, xi, chi, chit, m1s, m2s,
     &Ftarget 
      common /powers/ Q4,Q6,Q8,m14,m16,m18,m24,m26,m28
      common /acotlo/ acotlo
      common /vecax/ Sp, Sm, Rqp, Rqm
      common /delta1/ delta1

      Q4 = Q2 * Q2
      Q6 = Q4 * Q2
      Q8 = Q6 * Q2
      m14  = m1s * m1s
      m16  = m14 * m1s
      m18  = m16 * m1s
      m24  = m2s * m2s
      m26  = m24 * m2s
      m28  = m26 * m2s
      delta2 = delta(m1s,s,-Q2)
      delta12 = delta1 * delta1
      delta22 = delta2 * delta2
      spp = Q2 + m2s + m1s
      smm = Q2 - m2s - m1s
      spm = Q2 + m2s - m1s
      smp = Q2 - m2s + m1s
      diff = m1s+s+Q2-delta2
      sum  = m1s+s+Q2+delta2
      LN = dlog(diff/sum)

      hq1 = (1.d0-xhat)/8.d0 * (s-m2s)/s *
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

      hq1 = hq1 / acotlo(1)

      return
      end
