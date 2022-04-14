#ifndef NCTEQPP_VBP_H
#define NCTEQPP_VBP_H

extern "C" {

  void vbp_init(double& Fac0, double& Fac01, double* Fac1);

  void vbp_calculate(double Fac0, double Fac01, double* Fac1, double RS, double X1, double Q2, double Rerr, bool& doLO, double& Born, double& Onelop);
  
}

#endif //NCTEQPP_VBP_H