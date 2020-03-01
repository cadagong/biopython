from Bio.Seq import Seq

#better and simpler method for measure GC content in DNA strand
#uses Bio module from biopython library
def measureGCcontent(seq):
    Ccount = seq.count("C")
    Gcount = seq.count("G")
    GCcount = Ccount + Gcount
    seqLength = len(seq)

    GCcontent = (GCcount/seqLength)*100
    print("GC content: {:.2f}%".format(GCcontent))   


if __name__ == '__main__':
    print('Please enter your DNA or RNA sequence:')
    mySeq = input()    
    mySeq = mySeq.upper()
    measureGCcontent(mySeq)

   
        
        
