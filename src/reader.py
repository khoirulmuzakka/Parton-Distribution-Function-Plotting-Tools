import numpy as np
import pandas as pd
import yaml



class Data : 
    def __init__(self, dtlist, norm, penalty) : 
        self.norm = norm 
        self.penalty = penalty 
        pt = dtlist[0]
        self.ID = pt["IDDataSet"]
        self.typeExp = pt["TypeExp"]#
        self.A1 = pt["A1"]
        self.A2 = pt["A2"]
        self.Z1 = pt["Z1"]
        self.Z2 = pt["A1"]
        self.kinvarname = list(pt["KinVarVals"].keys() )
        self.table = self._createDF(dtlist)
        self.chi2nopen = np.sum(self.table["Chi2"])
        self.chi2 = self.chi2nopen+ penalty
        self.N = len(dtlist)
        self.chi2dof = self.chi2/self.N
        #display(self.table)

    def _createDF (self, dt) : 
        mat = []
        for p in dt : 
            row =list(p["KinVarVals"].values())+ [p["Data"], p["Theo"],p["StatError"], p["TotError"]]+ list(p["CorrError"].values()) + [p["Chi2"]]
            mat.append(row.copy())
        return pd.DataFrame(mat, columns= self.kinvarname+ ["Data", "Theory", "StatError", "TotError"]+ list(p["CorrError"].keys()) + ["Chi2"] )



class DataSets : 
    def __init__(self, file) : 
        ostream = open(file)
        self.yout = yaml.safe_load(ostream)
        self.Ntot =0
        self._getListData()
        
    
    def _getListData (self) : 
        pts = self.yout["Chi2Fcn"]["Snapshots"][0]["perPointBreakdown"]
        self.Ntot = len(pts)
        data = []
        dt = []
        ID = 999999999 # int(pts[0]["IDDataSet"]) 
        self.IDlist = [] #[ID]
        penPerID = {}  #{ID: 0.0}
        normPerID = {} #{ID:1.0}
        for i, pt in enumerate(pts) : 
            if (int(pt["IDDataSet"]) != ID) : 
                ID = int(pt["IDDataSet"])
                self.IDlist.append(ID)
                penPerID[ID] = 0.0
                normPerID[ID] = 1.0
                if (len(dt) !=0) : 
                    data.append(dt.copy())
                dt = []
                
            dt.append(pt) 
            if i== len(pts)-1 : 
                data.append(dt.copy())
        normStats = {}
        pts = self.yout["Chi2Fcn"]["NormInfo"]
        for p in pts : 
            p_ids = list(p["IDs"])
            N_ids = len(p_ids)
            pen_perid = p["Penalty"]/N_ids 
            norm_ids = p["Value"]
            for id in p_ids : 
                penPerID[id] = pen_perid
                normPerID[id] = norm_ids

        self.datasets = {}
        assert(len(self.IDlist) == len(data))
        for id, d in zip(self.IDlist, data) : 
            self.datasets[id] = Data (d, normPerID[id], penPerID[id])

        chi2tot =0.0
        chi2wopen =0.0
        for id, d in self.datasets.items() : 
            chi2tot = chi2tot + d.chi2 
            chi2wopen = chi2wopen+ d.chi2nopen 

        print("Chi2 total : ", chi2tot)
        print("Chi2 w/o penalty : ", chi2wopen)
        self.chi2dof =  chi2tot/self.Ntot

    def filterByExp (self, exp) : 
        dlist = []
        N = 0.0
        chi2 =0.0
        for id, d in self.datasets.items() : 
            if d.typeExp == exp : 
                dlist.append(d)
                N = N+ d.N 
                chi2 = chi2 + d.chi2 

        return [dlist, N, chi2, chi2/N]
    
    def getChi2Data (self) : 
        IDl = []
        Nl = []
        chi2L = []
        chi2dofL = []
        for id, dt in self.datasets.items() : 
            IDl.append(id) 
            Nl.append(dt.N)
            chi2L.append (dt.chi2)
            chi2dofL.append(dt.chi2dof)

        return pd.DataFrame({"ID" : IDl, "N": Nl, "Chi2Tot" : chi2L, "Chi2/N" : chi2dofL})

