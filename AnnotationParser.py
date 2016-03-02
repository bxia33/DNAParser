'''
Created on Feb 28, 2016

@author: BX
'''

# assume gff3 have 9 columns
class anno_parser:
    def __init__(self, filepath):        
        self.path = filepath
        self.annotations = []
        self.parseGFF3()
        return
    
    
    def parseGFF3(self):
        gff3format = ["seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"]
        with open(self.path) as f:
            for line in f:
                if line.startswith("#"): 
                    continue
                rows = line.strip().split("\t")
                
                # if the column number is not 9, skip that line
                if len(rows)!= len(gff3format):
                    continue
                else:
                #Normalize data
                    information = {
                        "seqid": None if rows[0] == "." else rows[0],
                        "source": None if rows[1] == "." else rows[1],
                        "type": None if rows[2] == "." else rows[2],
                        "start": None if rows[3] == "." else int(rows[3]),
                        "end": None if rows[4] == "." else int(rows[4]),
                        "score": None if rows[5] == "." else float(rows[5]),
                        "strand": None if rows[6] == "." else rows[6],
                        "phase": None if rows[7] == "." else rows[7],
                        "attributes": self.parseGFFAttributes(rows[8])
                    }
                    self.annotations.append(information)
        f.close()    
        return
    
    def parseGFFAttributes(self, attributeString):
        if attributeString == ".": 
            return {}
        attrs = {}
        for attribute in attributeString.split(";"):
            key, value = attribute.split("=")
            attrs[key] = value
        return attrs 
    
    # in description, it said whether hit a gene, not sure gene means the "type" in annotation file or whether it hit a
    # biological mean gene, because exon, CDS, mRNA should be considered either part of a gene or expression product of a gene.
    # here I assume "gene" means "type" in annotation file
    def hit(self, start, end):
        for dic in self.annotations:
            if dic["type"] =="gene":
                if start >= dic["start"] and end<=dic["end"]:
                    return True
        return False
                
                