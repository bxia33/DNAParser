'''
Created on Feb 28, 2016

@author: BX
'''
import collections
class parser:
    def __init__(self, file_path, length, repeat, endSeq, limit_seq, up, bot,unique_start, unique_end):
        self.fpath = file_path
        self.DNA = ""                      #DNA sequence
        self.length = length                #query length which 23 in this case
        self.repeat = repeat                #limit repeat, which is 5
        self.endSeq = endSeq                # end sequence which is "GG"
        self.limit_seq = limit_seq          # no "ATG"
        self.up = up                        # assume 20% inclusive 
        self.bot = bot                      # assume 80% inclusive
        self.unique_start = unique_start    # unique sequence start at index of 13-1
        self.unique_end = unique_end        # unique sequence end at index of 20-1
        self.candidates = {}                                    # store candidate position and GC contents
        self.candidate_repeats =collections.defaultdict(int)     # store sequence repeat times 13-20 nucleotides
        self.DNAparser()
        return
    
    
    def DNAparser(self):
        with open(self.fpath) as f:
            for line in f:
                rows = line.strip()
                rows.replace(' ','')
                self.DNA+= rows
            f.close()
                    
        self.candidates = self.GCcalculater(self.DNA, self.length, self.up, self.bot, self.unique_start, self.unique_end, self.repeat)
        self.candidates = self.endwithSeq(self.DNA, self.candidates, self.length, self.endSeq)
#         self.candidates = self.isHomopolymers(self.DNA, self.candidates, self.length,self.repeat)
        self.candidates = self.containSeq(self.DNA, self.candidates, self.length,self.limit_seq)
        return 

    def GCcalculater(self, sequence, length, up, bot, start, end, R):
        result = {}
        size = len(sequence)
        i = 0
        bot = length*up
        up = length*bot
        
        repeat =1
        
        while i < size - length+1:
            seed = sequence[i+start-1:i+end]
            self.candidate_repeats[seed]+=1        
            if i-1 in result:
                count = result[i-1]
                if sequence[i+22] == 'G' or sequence[i+22]=='C':
                    count += 1
                if sequence[i-1] == 'G' or sequence[i-1] =='C':
                    count -=1
                result[i]=count
            else:
                count = 0
                for j in xrange(i, i +23):
                    if sequence[j] == 'G' or sequence[j]=='C':
                        count += 1
                result[i] = count       
            i+=1
            # try to delete homopolymers during the process
            if i-1>=0:
                if sequence[i]==sequence[i-1]:
                    repeat +=1
                    if repeat >=R:
                        for k in xrange(max(0,i-length+1),i-R+2):
                            if k in result:
                                del result[k]        
                else:
                    repeat =1
        # remove additional homopolymers in the last 23bp
        for i in range(size-length+1, size):
            if i-1>=0:
                if sequence[i]==sequence[i-1]:
                    repeat +=1
                    if repeat >=R:
                        for k in xrange(max(0,i-length+1),i-R+2):
                            if k in result:
                                del result[k]        
                else:
                    repeat =1            
        keys = result.keys()
        for key in keys:
            if bot<=result[key] <=up:
                pass
            else:
                if key in result:
                    del result[key]
        return result

    def endwithSeq(self, sequence, candidates, length, seq):
        end = 0 - len(seq)
        keys = candidates.keys()
        for key in keys:
            seed = sequence[key:key+length]
            if seed[end:length] == seq:
                pass
            else:
                del candidates[key]
        return candidates

    # currently using brute force to remove homopolymers, need a better algo.
    def isHomopolymers(self, DNA, candidates, length, R):
        keys = candidates.keys()
        for key in keys:
            sequence = DNA[key: key +length]
            repeat = 1
            for i in range(1,length):
                if sequence[i] == sequence[i-1]:
                    repeat +=1
                    if repeat >= R and key in candidates:
                        del candidates[key]
                else:
                    repeat =1
        return candidates

    def containSeq(self, DNA, candidates, lengh,seq):
        keys = candidates.keys()
        for key in keys:
            sequence = DNA[key: key +lengh]
            if sequence.find(seq) != -1:
                del candidates[key]
        return candidates






            
         
    


