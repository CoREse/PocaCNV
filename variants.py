class CNV():
    def __init__(self,BreakLeft=0,BreakRight=0,Phased=False,Alleles=[],Samples=[],Score=None):#Samples is like [(SampleID,(0,1)),(SampleID,(1,1))]
        self.BreakLeft=BreakLeft
        self.BreakRight=BreakRight
        self.SVLength=BreakRight-BreakLeft
        self.Phased=Phased
        self.Alleles=Alleles
        self.Samples=Samples
        self.Score=Score
        self.SVType="NONE"
        self.deductSVType()
    def deductSVType(self):
        if len(self.Alleles)==0:
            return
        if len(self.Alleles)>1:
            self.SVType="CNV"
        elif self.Alleles[0]=="<CN0>":
            self.SVType="DEL"
        else:
            self.SVType="DUP"