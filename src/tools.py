import lhapdf
import numpy as np
from math import *


lhapdf.pathsAppend('./../LHAPDF-files/')

class PDFs : 
    def __init__(self, path, errorType, A=1, Z=1) : 
        '''
        args : 
            path  : ( string )folder name of the lhapdf file. 
            errorType : (string)
                        'sym'  - symmetric Hessian errors: delta = +/- 0.5*\sqrt{\sum_{i=1}^N (f_i^{+} - f_i^{-})^2}
                        'Asym' - symmetric Hessian errors: delta_+ = \sqrt(\sum_i(max(f_i^{+}-f_0, f_i^{-}-f_0)^2))
                                                        delta_- = \sqrt(\sum_i(max(f_0-f_i^{+}, f_0-f_i^{-})^2))
                        'mc2hessian' - symmetric Hessian errors but assuming that error PDFs for only one direction
                                    (plus or minus) are dumped as is the case when MC2Hessian is used to convert
                                    from MC replicas: delta = +/- \sqrt{\sum_{i=1}^N (f_i^{+} - f_0)^2}
                        'MCreplica' - for Monte Carlo replicas: delta = +/- \sqrt{1/Nrep*\sum_{i=1}^{Nrep} (f_i^{+} - f_0)^2}
        '''
        self.pdfset = lhapdf.getPDFSet(path).mkPDFs()
        self.errorType = errorType 
        self.A=A
        self.Z =Z
        self.N= A-Z
    
    def getPDFvalue (self, pdf, obsName, x, Q, fullNuc=False) : 
        '''
        Parton name can only be : "g", "t", "b", "c", "s", "u", d", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar", 
                                      "uv", "dv", "ssb", "smsb", "ubdb", "dboub", "kappa"
        here "kappa=ssb/ubdb"
        
        '''
        if (obsName == "g") : 
            return pdf.xfxQ(21,x,Q)
        elif (obsName == "t")  :
            return pdf.xfxQ(6,x,Q)
        elif (obsName == "tbar")  :
            return pdf.xfxQ(-6,x,Q)
        elif (obsName == "b")  :
            return pdf.xfxQ(5,x,Q)
        elif (obsName == "bbar")  :
            return pdf.xfxQ(-5,x,Q)
        elif (obsName == "c")  :
            return pdf.xfxQ(4,x,Q)
        elif (obsName == "cbar")  :
            return pdf.xfxQ(-4,x,Q)
        elif (obsName == "s")  :
            return pdf.xfxQ(3,x,Q)
        elif (obsName == "sbar")  :
            return pdf.xfxQ(-3,x,Q)
        elif (obsName =="ssb") : 
            return pdf.xfxQ(3,x,Q)+ pdf.xfxQ(-3,x,Q)
        elif (obsName =="smsb") : 
            return pdf.xfxQ(3,x,Q)- pdf.xfxQ(-3,x,Q)   
        elif (obsName == "u") : 
            if (fullNuc) :
                return (self.Z*pdf.xfxQ(2,x, Q) + self.N*pdf.xfxQ(1,x, Q))/self.A
            else :  
                return pdf.xfxQ(2,x, Q)
        elif (obsName =="ubar") :
            if (fullNuc) : 
                return (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A
            else : 
                return pdf.xfxQ(-2,x, Q)
        elif (obsName =="d") :
            if (fullNuc) : 
                return (self.Z*pdf.xfxQ(1,x, Q) + self.N*pdf.xfxQ(2,x, Q))/self.A
            else : 
                return pdf.xfxQ(1,x, Q)
        elif (obsName =="dbar") :
            if (fullNuc) : 
                return (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(-1,x, Q)
        elif (obsName == "uv") : 
            if (fullNuc) : 
                return (self.Z*pdf.xfxQ(2,x, Q) + self.N*pdf.xfxQ(1,x, Q))/self.A - (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A
            else : 
                return pdf.xfxQ(2,x, Q)-pdf.xfxQ(-2,x, Q)
        elif (obsName =="dv") : 
            if (fullNuc) : 
                return (self.Z*pdf.xfxQ(1,x, Q) + self.N*pdf.xfxQ(2,x, Q))/self.A - (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(1,x, Q) - pdf.xfxQ(-1,x, Q)
        elif (obsName== "ubdb") : 
            if (fullNuc) : 
                return (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A + (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(-2,x, Q)+ pdf.xfxQ(-1,x, Q)
        elif (obsName =="kappa") : 
            if (fullNuc) : 
                return (pdf.xfxQ(3,x,Q)+ pdf.xfxQ(-3,x,Q))/((self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A + (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A)
            else : 
                return (pdf.xfxQ(3,x,Q)+ pdf.xfxQ(-3,x,Q))/(pdf.xfxQ(-2,x, Q)+ pdf.xfxQ(-1,x, Q))
        else : 
            print ("Unknown flavor! Exiting...")
            exit
    

    def getPDFerrors(self, xList, Q, obs, fullNuc=False) : 
        central = []
        DeltaPlus = []
        DeltaMinus = []
        for x in xList : 
            cent = self.getPDFvalue(self.pdfset[0], obs, x, Q, fullNuc)
            central.append(cent)
            if (self.errorType == "Asym") : 
                dp = np.sqrt(np.sum(np.array([
                     max( self.getPDFvalue(pdfp, obs, x, Q, fullNuc)-cent, self.getPDFvalue(pdfm, obs, x, Q, fullNuc)-cent,0.  )**2 
                     for pdfp, pdfm in zip(self.pdfset[1::2], self.pdfset[2::2]) 
                     ] )))
                
                dm = np.sqrt(np.sum(np.array([ 
                     max( cent-self.getPDFvalue(pdfp, obs, x, Q, fullNuc), cent-self.getPDFvalue(pdfm, obs, x, Q, fullNuc),0. )**2
                                     for pdfp, pdfm in zip(self.pdfset[1::2], self.pdfset[2::2] )
                     ] )))

            elif (self.errorType == "sym") : 
                dp = 0.5*np.sqrt(np.sum(np.array([ 
                        ( self.getPDFvalue(pdfp, obs, x, Q, fullNuc) - self.getPDFvalue(pdfm, obs, x, Q, fullNuc) )**2
                                     for pdfp, pdfm in zip(self.pdfset[1::2], self.pdfset[2::2] )
                    ] )))
                dm = dp

            elif (self.errorType == "mc2hessian") : 
                dp = np.sqrt(np.sum(np.array([ 
                        ( self.getPDFvalue(pdf, obs, x, Q, fullNuc) - cent)**2
                        for pdf in self.pdfset[1::1] 
                    ] )))
                dm = dp

            elif (self.errorType == "MCreplica") : 
                Nrep = len(self.pdfset[1::1] )
                dp = np.sqrt(np.sum(np.array([ 
                        ( self.getPDFvalue(pdf, obs, x, Q, fullNuc) - cent)**2
                        for pdf in self.pdfset[1::1] 
                    ] )))/Nrep
                dm = dp
            else : 
                print("Unknown error type! Exiting ...")
                exit
            
            DeltaPlus.append(dp)
            DeltaMinus.append(dm)
            
        central = np.array(central)
        DeltaPlus = np.array(DeltaPlus)
        DeltaMinus = np.array(DeltaMinus)

        return central, DeltaPlus, DeltaMinus

















