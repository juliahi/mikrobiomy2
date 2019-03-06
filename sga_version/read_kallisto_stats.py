
""" Class for analysing kallisto stats.txt files """

import pandas



class Stats:
    def __init__(self, filename):
        names = ['kmers1', 'kmers2', 'contigs1', 'contigs2', 'contigs1+2',
                 'tu1', 'tu2', 'tu1+2', 'ti1', 'ti2', 'ti1+2']
        self.data = pandas.read_csv(filename, sep='\t', names=names)
        
        self.n = float(len(self.data))
        # print filename, self.n, len(self.get_mapped())
        
    def mapped(self):
        return self.data["ti1+2"] >= 1

    def nonunique(self):    
        return self.data["ti1+2"] > 1

    def unique(self):
        return self.data["ti1+2"] == 1
    
    def zerokmers(self):
        return (self.data["kmers1"] == 0) & (self.data["kmers2"] == 0)
        
    def conflicts(self):
        return (self.data["ti1+2"] == 0) & (self.data["tu1+2"] > 0)
    
    def get_mapped(self):
        return self.data[self.mapped()]
    
    def fraction_mapped(self):
        return sum(self.mapped())/self.n 
    
    def get_nonunique(self):
        return self.data[self.nonunique()]
    
    def fraction_nonunique(self):
        return sum(self.nonunique())/self.n

    def fraction_unique(self):
        return sum(self.unique())/self.n
    
    def get_zerokmers(self):
        return self.data[self.zerokmers()]
    
    def fraction_zerokmers(self):
        return sum(self.zerokmers())/self.n
        
    def get_conflicts(self):
        return self.data[self.conflicts()]

    def fraction_conflicts(self):
        return sum(self.conflicts())/self.n

    

    
