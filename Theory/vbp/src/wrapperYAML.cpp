#include "YAMLSettings.h"
#include <iostream>

using namespace std;

extern "C"{
	// /ThpVbp/ ThVbp(3), LVbp(3)
	extern struct THPVBP{
		double THVBP[3];
		int LVBP[3];
	} thpvbp_;


	void initYAML_(){
		string key("vbpNCTEQ15");
		SettingDefs settingDefs; 
		// If Iscl = 0 : keep Amu = Q in Vbp calculations
		//         = 1 : treat Amu as independent, and vary qscl in Amu = qscal * Q
		settingDefs[key].push_back(make_unique<ScalarSetting<double>>("IRscl", SettingRange<double>{0.0,1.0}, 0.0, thpvbp_.THVBP[0]));

		// load settings from a file
		string fname("vbp.yaml");
		YAML::Node inputSettings;
		try {
			cout << "Loading the settings from `" << fname << "`..." << endl;
			 inputSettings = YAML::LoadFile(fname);
			cout << "...done loading the settings from `" << fname << "`." << endl;
		}
		catch (...) {
			cout << "Error loading file `" << fname << "`. Exiting!" << endl;
			exit(-1) ;
		}
		// process the settings
		YAMLSettings yamlSettings(inputSettings, settingDefs);
	}
}

