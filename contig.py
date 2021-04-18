from array import array
import multiprocessing.sharedctypes

class Contig:
    def __init__(self,name,length,RDWindowSize,Genome):
        self.Name=name
        self.Length=int(length/RDWindowSize)+(1 if length%RDWindowSize!=0 else 0)
        self.NLength=length
        self.RDWindowSize=RDWindowSize
        self.RDWindows=[]
        self.SampleNames=[]
        self.MixedRDRs=None
        self.RDWindowStandards=None
        self.ContigSampleReadCounts=[]
        self.ContigReadCount=None
        self.SampleReadCount=self.ContigSampleReadCounts
        self.Genome=Genome
        self.Ploidies=[]
        self.DRPs=[]
    
    def calcContigReadCount(self):
        self.ContigReadCount=0
        for rc in self.ContigSampleReadCounts:
            self.ContigReadCount+=rc
    
    def genVacant(self,TheGenome):#gen with no sample
        new=Contig(self.Name,self.NLength,self.RDWindowSize,TheGenome)
        return new       

    def addSample(self,name=""):
        self.RDWindows.append(array("f",[0]*self.Length))
        #self.RDWindows.append(multiprocessing.sharedctypes.RawArray("f",[0]*self.Length))
        self.SampleNames.append(name)
        self.ContigSampleReadCounts.append(0)
        self.DRPs.append([])
        self.Ploidies.append(2)
    
    def addSampleWithData(self,SampleName,SampleRDWindows,ContigSampleReadCount,SampleDRPs,Ploidy):
        self.RDWindows.append(SampleRDWindows)
        self.SampleNames.append(SampleName)
        self.ContigSampleReadCounts.append(ContigSampleReadCount)
        self.DRPs.append(SampleDRPs)
        self.Ploidies.append(Ploidy)
    
    def changeSampleName(self, Index, Name):
        self.SampleNames[Index]=Name
    
    def clearMemory(self):
        self.RDWindows=None
        self.MixedRDRs=None
        self.DRPs=None
        self.RDWindowsAcc=None
    
    def getNew(self,TheGenome):
        return Contig(self.Name,self.NLength,self.RDWindowSize,TheGenome)
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
        self.RefID=[]
        self.GenomeLength=None
    
    '''
    def duplicate(self):
        new=Genome()
        new.Name=self.Name
        new.Contigs=self.Contigs.copy()
        new.Index=self.Index.copy()
        new.SampleNames=self.SampleNames.copy()
        new.SampleN=self.SampleN
        new.RefID=self.RefID.copy()
        return new'''
    def genVacant(self):#gen with no sample
        new=Genome()
        new.Name=self.Name
        for c in self.Contigs:
            new.Contigs.append(c.genVacant(new))
        new.Index=self.Index.copy()
        new.RefID=self.RefID.copy()
        new.GenomeLength=self.GenomeLength
        return new

    def getRefID(self,ContigID):
        return RefID[ContigID]
    
    def getContigID(self,ReferenceID):
        ID=-1
        for i in range(len(self.RefID)):
            if self.RefID[i]==ReferenceID:
                ID=i
                break
        return ID
    def append(self,c):
        self.Contigs.append(c)
        self.Index[c.Name]=len(self.Contigs)-1
    def __getitem__(self,i):
        return self.Contigs[i]
    def __setitem__(self,i,v):
        self.Contigs[i]=v
        self.Index[v.Name]=i
        return self.Contigs[i]
    def hasContig(self,ContigName):
        for c in self.Contigs:
            if c.Name==ContigName:
                return True
        return False
    def get(self,ContigName):
        try:
            Result=self.Contigs[self.Index[ContigName]]
        except:
            return None
        return Result
    def addSample(self,name):
        for c in self.Contigs:
            c.addSample(name)
        self.SampleNames.append(name)
        self.SampleN+=1
    def changeSampleName(self, Index, Name):
        for c in self.Contigs:
            c.changeSampleName(Index, Name)
        self.SampleNames[Index]=Name
    def __len__(self):
        return len(self.Contigs)
    def calcContigReadCounts(self):
        for c in self.Contigs:
            c.calcContigReadCount()
    def setGenomeLength(self,RefLength):
        self.GenomeLength=RefLength
    def clearMemory(self):
        for c in self.Contigs:
            c.clearMemory()