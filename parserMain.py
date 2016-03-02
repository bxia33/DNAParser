'''
Created on Feb 28, 2016

@author: BX
'''
import DNAparser as dp
from AnnotationParser import anno_parser as ap
from Target import target
import sys

if __name__ == "__main__":
    f = "Sequence.fa"
    a = "Annotation.gff3"
    result = []
    DNA_parser = dp.parser(f, 23, 5, "GG" , "ATG",0.2, 0.8 , 13, 20)
    ap_parser = ap(a)
    for seed in DNA_parser.candidates:
        eightmer = DNA_parser.DNA[seed+13-1:seed+20]
        seed_start = seed
        seed_end = seed+23
        hit = ap_parser.hit(seed_start, seed_end)
        copys = DNA_parser.candidate_repeats[eightmer]
        
        dna_target = target(DNA_parser.DNA[seed:seed+20], \
                            DNA_parser.DNA[seed+20:seed+23],\
                            int(float(DNA_parser.candidates[seed])/23*100),\
                            hit, copys)
        result.append(dna_target)
    
    result.sort(key=lambda x: x.score, reverse=True)
    
    orig_stdout = sys.stdout
    f = file('output.txt', 'w')
    sys.stdout = f
    
    print "rank", '\t', "sequence", '\t','\t','\t','\t',"PAM",'\t',"GC content",'\t',"Gene",'\t',"Copy Number",'\t',"Score"
    for r in range(len(result)):
        result[r] 
        print r+1, '\t', '\t',result[r].seq,'\t', result[r].PAM, '\t',result[r].GC, '\t','\t','\t',result[r].gene,\
        '\t', result[r].copys, '\t','\t','\t','\t',result[r].score
    
    sys.stdout = orig_stdout
    f.close()

    
    
    
