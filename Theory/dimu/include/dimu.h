//
// Created by Florian Lyonnet on 9/1/17.
//

#ifndef NCTEQPP_DIMU_H
#define NCTEQPP_DIMU_H

extern "C" {

	void dimu_calculate(double x, double y, double ehad, double eps, int iproc, bool& doLO, double& Born, double& Onelop );
	void sig_charm(double x, double y, double enu, int iproc, bool& doLO, double& sig );
};

#endif //NCTEQPP_DIMU_H
