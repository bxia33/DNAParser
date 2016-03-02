'''
Created on Feb 28, 2016

@author: BX
'''
class target:
    def __init__(self, sequence, PAM, GC, gene, copys):
        self.seq = sequence
        self.PAM = PAM
        self.GC = GC
        self.gene = gene
        self.copys = copys
        self.score = 0
        self.caculator()
        return
    
    def __str__(self):
        return str(self.score)
    
    
    #I assume 3 scoring requirements are equally important
    # each have 100 points
    # if PAM not start with G or C, -100
    # each GC contents difference will lose 1 point
    # each repeat copy will lose 10 points
    def caculator(self):
        total = 300
        if self.PAM[0] == 'G' or self.PAM[0] == 'C':
            pass
        else:
            total -=100
        total -= abs(self.GC -45)
        
        total -= 10*self.copys
        self.score = total
        return

            
        
        
        
        
        