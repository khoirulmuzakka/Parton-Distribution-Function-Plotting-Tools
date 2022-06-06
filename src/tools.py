import lhapdf
import numpy as np
from math import *
import sys


lhapdf.pathsAppend('./../LHAPDF-files/')

class PDFs : 
    def __init__(self, path, errorType, A=1, Z=1, fullNuc = False) : 
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
        self.Nmembers = len(self.pdfset[0::1])
        self.fullNuc = fullNuc
    
    def getPDFvalue (self, pdf, obsName, x, Q) : 
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
            if (self.fullNuc) :
                return (self.Z*pdf.xfxQ(2,x, Q) + self.N*pdf.xfxQ(1,x, Q))/self.A
            else :  
                return pdf.xfxQ(2,x, Q)
        elif (obsName =="ubar") :
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A
            else : 
                return pdf.xfxQ(-2,x, Q)
        elif (obsName =="d") :
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(1,x, Q) + self.N*pdf.xfxQ(2,x, Q))/self.A
            else : 
                return pdf.xfxQ(1,x, Q)
        elif (obsName =="dbar") :
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(-1,x, Q)
        elif (obsName == "uv") : 
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(2,x, Q) + self.N*pdf.xfxQ(1,x, Q))/self.A - (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A
            else : 
                return pdf.xfxQ(2,x, Q)-pdf.xfxQ(-2,x, Q)
        elif (obsName =="dv") : 
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(1,x, Q) + self.N*pdf.xfxQ(2,x, Q))/self.A - (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(1,x, Q) - pdf.xfxQ(-1,x, Q)
        elif (obsName== "ubdb") : 
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A + (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(-2,x, Q)+ pdf.xfxQ(-1,x, Q)
        elif (obsName== "dboub") : 
            if (self.fullNuc) : 
                denom = (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))
                if (denom==0.0) : 
                    return 0.0
                else : 
                    return (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q)) /denom
            else : 
                denom = pdf.xfxQ(-2,x, Q)
                if (denom==0.0) : 
                    return 0.0
                else : 
                    return pdf.xfxQ(-1,x, Q)/denom
        elif (obsName== "ubmdb") : 
            if (self.fullNuc) : 
                return (self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A - (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A
            else : 
                return pdf.xfxQ(-2,x, Q)- pdf.xfxQ(-1,x, Q)
        elif (obsName =="kappa") : 
            if (self.fullNuc) : 
                return (pdf.xfxQ(3,x,Q)+ pdf.xfxQ(-3,x,Q))/((self.Z*pdf.xfxQ(-2,x, Q) + self.N*pdf.xfxQ(-1,x, Q))/self.A + (self.Z*pdf.xfxQ(-1,x, Q) + self.N*pdf.xfxQ(-2,x, Q))/self.A)
            else : 
                return (pdf.xfxQ(3,x,Q)+ pdf.xfxQ(-3,x,Q))/(pdf.xfxQ(-2,x, Q)+ pdf.xfxQ(-1,x, Q))
        else : 
            print ("Unknown flavor! Exiting...")
            sys.exit()
    

    def getPDFerrors(self, xList, Q, obs) : 
        """
        Get PDFs and their uncertainties. 
        Arguments : 
              xList : list of x values
              Q     : Virtuality value
              obs   : flavor, possible values : ""g", "t", "b", "c", "s", "u", d", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar", 
                                      "uv", "dv", "ssb", "smsb", "ubdb", "dboub", "kappa"
              FullNuc : boolean flag to swicth on averaging procedure for u, ubar, d, dbar flavors. The full PDF is then evaluated as : 
                       f_i^A = (Z f^{p/A}+ N f^{n/A})/A
            
        return : central, deltaplus, deltaminus 
              central : central PDF (numpy array)
              deltaplus : Upper error (numpy array)
              deltaminus : lower errror (numpy array)
        """
        central = []
        DeltaPlus = []
        DeltaMinus = []
        for x in xList : 
            cent = self.getPDFvalue(self.pdfset[0], obs, x, Q)
            central.append(cent)
            if (self.errorType == "Asym") : 
                dp = np.sqrt(np.sum(np.array([
                     max( self.getPDFvalue(pdfp, obs, x, Q)-cent, self.getPDFvalue(pdfm, obs, x, Q)-cent,0.  )**2 
                     for pdfp, pdfm in zip(self.pdfset[1::2], self.pdfset[2::2]) 
                     ] )))
                
                dm = np.sqrt(np.sum(np.array([ 
                     max( cent-self.getPDFvalue(pdfp, obs, x, Q), cent-self.getPDFvalue(pdfm, obs, x, Q),0. )**2
                                     for pdfp, pdfm in zip(self.pdfset[1::2], self.pdfset[2::2] )
                     ] )))

            elif (self.errorType == "sym") : 
                dp = 0.5*np.sqrt(np.sum(np.array([ 
                        ( self.getPDFvalue(pdfp, obs, x, Q) - self.getPDFvalue(pdfm, obs, x, Q) )**2
                                     for pdfp, pdfm in zip(self.pdfset[1::2], self.pdfset[2::2] )
                    ] )))
                dm = dp

            elif (self.errorType == "mc2hessian") : 
                dp = np.sqrt(np.sum(np.array([ 
                        ( self.getPDFvalue(pdf, obs, x, Q) - cent)**2
                        for pdf in self.pdfset[1::1] 
                    ] )))
                dm = dp

            elif (self.errorType == "MCreplica") : 
                Nrep = len(self.pdfset[1::1] )
                dp = np.sqrt(np.sum(np.array([ 
                        ( self.getPDFvalue(pdf, obs, x, Q) - cent)**2
                        for pdf in self.pdfset[1::1] 
                    ] ))/(Nrep-1))
                dm = dp
            else : 
                print("Unknown error type! Exiting ...")
                sys.exit()
            
            DeltaPlus.append(dp)
            DeltaMinus.append(dm)
            
        central = np.array(central)
        DeltaPlus = np.array(DeltaPlus)
        DeltaMinus = np.array(DeltaMinus)

        return central, DeltaPlus, DeltaMinus

    

    def getPDF(self, x, Q, obs) : 
        """
        Function to get PDF value for a given x, Q values. 
        Arguments : 
            x : momentum fraction x
            Q : Factorization scale
            obs : flavor, possible values : ""g", "t", "b", "c", "s", "u", d", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar", 
                                      "uv", "dv", "ssb", "smsb", "ubdb", "dboub", "kappa"
            fullNuc = False : boolean flag to set whether nucleon averaging procedure (to get the full nuclear PDFs) should be performed. 
        return : 
            numpy array of size N, where N is the number of PDF members. The zeroth element is the central, while the rest are from error PDFs
        """
        N=len(self.pdfset[0::1] )
        a=[]
        for i in range(N) : 
            val = self.getPDFvalue(self.pdfset[i], obs, x, Q)
            a.append(val)
        return np.array(a)

    

















