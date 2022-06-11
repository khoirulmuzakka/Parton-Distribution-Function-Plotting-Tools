Naming convention of files:
For eta_dijet distribuiton, we have a typical file with name CMS5020LO.avgetaj12_ptavg1.dat
CMS5020LO means the LO results for CMS 5020 GeV measurement.
The .avgetaj12_ptavg1. is for one of the five ptavg fiducial regions.
Please find more information about the definition of each region in the plots like CMS5020_avgetaj12_ptavg1.pdf

The data and plot are with normalized or original value. From CMS5020_avgetaj12_ptavg1.pdf one could see the eta_dijet distribution at LO and NLO while the normalized version is inside CMS5020_avgetaj12norm_ptavg1.pdf. To quanitfy the dominant numberical uncertainties, more details of the original value can be find in CMS5020_avgetaj12_ptavg1_detail.pdf. 
The normalized data inside for example CMS5020LO.avgetaj12norm_ptavg1.dat satisfy the relation \sum_bin binValu*binWidth = 1.

Data structure of files:
Differential results for the 7 scales (mu_f,mur) \in {(1,1), (.5,.5), (2,2), (1,.5), (1,2), (.5,1), (2,1)}*p_T^{avg} are provided in the files with the format:
BinLowerEdge,BinCentral,BinUpperEdge,  (dsigma/dpt and numerical error) for the 7 scales in the order given above. The values are normalized to the binwidth defined by (BinUpperEdge-BinLowerEdge).
