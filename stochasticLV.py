"""
Stochastic version of the classic LV species competition model simulated using the Gillespie
algorithm. 

Author: Ata Kalirad
Date: 2025-07-15
"""
import numpy as np

class SAA(object):
    
    def __init__(self, N1_0, N2_0, r1, r2, α12, α21, k):
        self.r1 = r1
        self.r2 = r2
        self.α12 = α12
        self.α21 = α21
        self.k = k
        self.set_inti_cond(N1_0, N2_0)
        self.gen_stoichiometry()
        self.init_propenisties()
        self.terminate = False


    def set_inti_cond(self, N1, N2):
        self.init_cond = {
                            "N1": [N1],
                            "N2": [N2],
                            "time": [0.0],
                        }
        
    def init_propenisties(self):
        self.prop =  {0: lambda pop: self.r1*pop['N1'][-1],
                      1: lambda pop: self.r1*pop['N1'][-1]*np.divide(pop['N1'][-1] + self.α12*pop['N2'][-1], self.k),
                      2: lambda pop: self.r2*pop['N2'][-1],
                      3: lambda pop: self.r2*pop['N2'][-1]*np.divide(pop['N2'][-1] + self.α21*pop['N1'][-1], self.k)}

    def get_propensities(self):
        return [i(self.init_cond) for i in self.prop.values()]
        

    def gen_stoichiometry(self):
        self.stoichiometry = {
                        0: {"N1": +1, "N2": 0},
                        1: {"N1": -1, "N2": 0},
                        2: {"N1": 0, "N2": +1},
                        3: {"N1": 0, "N2": -1},
                    }

    def execute_next_reaction(self):
        prop_list = self.get_propensities()
        non_zeros = [i for i,j in enumerate(prop_list) if j != 0]
        if len(non_zeros) > 0:
            times = [(1/prop_list[i])*np.log(1/np.random.uniform()) for i in non_zeros]
            nex_reaction = np.argmin(times)
            for i in self.init_cond.keys():
                if i != 'time':
                    self.init_cond[i].append(self.init_cond[i][-1] + self.stoichiometry[non_zeros[nex_reaction]][i])
                else:
                    self.init_cond[i].append(self.init_cond[i][-1] + times[nex_reaction])
        else:
            self.terminate = True

    def simulate(self, t_lim):
        while self.init_cond['time'][-1] < t_lim:
                self.execute_next_reaction()