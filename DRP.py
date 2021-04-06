class DRP:
    def __init__(self,start,end,innerstart,innerend,SVT):#SVT:SupportedVariantType=0,loss,1,gain
        #self.SampleName=samplename
        self.Start=start
        self.End=end
        self.InnerStart=innerstart
        self.InnerEnd=innerend
        self.isize=end-start
        self.SupportedVariantType=SVT
    
    def __repr__(self):
        return "%s-%s(%s-%s)"%(self.Start,self.End,self.InnerStart,self.InnerEnd)