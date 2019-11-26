class Contig:
    def __init__(self,name,length):
        self.Name=name
        self.Length=length
        self.RDWindows=[]
        self.SampleNames=[]
    
    def addSample(self,name=""):
        self.RDWindows.append([0]*self.Length)
        self.SampleNames.append(name)
    """
    def __getitem__(self,i):
        return self.RDWindows[i]
    def __setitem__(self,i,v):
        self.RDWindows[i]=v
        return self.RDWindows[i]
    """

class Genome:
    def __init__(self,name=""):
        self.Name=name
        self.Contigs=[]
        self.Index={}
        self.SampleNames=[]
        self.SampleN=0
    def append(self,c):
        self.Contigs.append(c)
        self.Index[c.Name]=len(self.Contigs)-1
    def __getitem__(self,i):
        return self.Contigs[i]
    def __setitem__(self,i,v):
        self.Contigs[i]=v
        self.Index[v.Name]=i
        return self.Contigs[i]
    def get(self,ContigName):
        return self.Contigs[self.Index[ContigName]]
    def addSample(self,name):
        for c in self.Contigs:
            c.addSample(name)
        self.SampleNames.append(name)
        self.SampleN+=1
    def __len__(self):
        return len(self.Contigs)
