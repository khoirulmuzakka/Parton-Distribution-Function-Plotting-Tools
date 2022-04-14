import numpy as np
import sys

class F123LO : 
    def __init__(self) : 
        self.flv = ["t", "b", "c", "s", "u", "d", "dbar", "ubar", "sbar", "cbar", "bbar", "tbar"]

    def F1 (self, x, Q, probe, pdf) :    
        pdfval ={}
        for fl in self.flv : 
            pdfval[fl] = pdf.getPDF(x, Q, fl)

        if (probe == "gamma") : 
            arr = ( 4.*(pdfval["u"]+ pdfval["ubar"] + pdfval["c"] + pdfval["cbar"]+ pdfval["t"]+ pdfval["tbar"]) +
                     pdfval["d"]+ pdfval["dbar"] + pdfval["s"] + pdfval["sbar"]+ pdfval["b"]+ pdfval["bbar"])/18
            cent = arr[0]
            pl = arr[1::2]
            mi = arr[2::2]
            err = 0.5*np.sqrt(np.sum(np.array([
                (plus-minus)**2 for plus, minus in zip(pl, mi)
            ])))
            return cent, err
        elif (probe == "W+") : 
            arr = pdfval["d"]+ pdfval["s"]+pdfval["b"] + pdfval["ubar"]+ pdfval["cbar"]+ pdfval["tbar"]
            cent = arr[0]
            pl = arr[1::2]
            mi = arr[2::2]
            err = 0.5*np.sqrt(np.sum(np.array([
                (plus-minus)**2 for plus, minus in zip(pl, mi)
            ])))
            return cent, err
        elif (probe == "W-") : 
            arr = pdfval["dbar"]+ pdfval["sbar"]+pdfval["bbar"] + pdfval["u"]+ pdfval["c"]+ pdfval["t"]
            cent = arr[0]
            pl = arr[1::2]
            mi = arr[2::2]
            err = 0.5*np.sqrt(np.sum(np.array([
                (plus-minus)**2 for plus, minus in zip(pl, mi)
            ])))
            return cent, err
        else : 
            print("Unknown probe! Probe must be ’gamma’, ’W+’, or ’W-’")
            sys.exit()
        
    def F2 (self, x, Q, probe, pdf) : 
        cent, err = self.F1(x, Q, probe, pdf)
        return 2.0*x*cent, 2.0*x*err
    



class DISCrossSectionLO : 
    pass

