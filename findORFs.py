import os
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from compareProteinProduct import *

#dictionary: key is codon, value is coresponding amino acid
codonTable = {
	'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
	'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
	'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
	'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
	'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
	'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
	'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
	'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
	'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
	'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
	'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
	'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
	'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
	'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
	'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
	'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

variantsArray = []


#get all possible protein sequences from input dna strand
def getORFs(seq, proteinList):    
    protein = ''
    isORF = False            
    for c in range(0, len(seq), 3):
        codon = seq[c:c+3]            
        if codonTable.get(codon) == 'M' and isORF == False:
            isORF = True
            protein += 'M'                                   
        elif isORF == True: 
            if codonTable.get(codon) == '_':                    
                proteinList.append(protein)
                protein = ''
                isORF = False
            elif codon in codonTable.keys():                                
                protein += codonTable.get(codon)                                              
            else:
                pass
        else:            
            pass

    

#generate the three possible reading frames for input dna strand
def generateReadingFrames(seq):
    frame1 = seq
    frame2 = seq[1:len(seq)]
    frame3 = seq[2:len(seq)]

    frameslist = [frame1, frame2, frame3]       
    return frameslist


#get complementary DNA strand
def getSequenceComplement(seq):
    seq = seq[len(seq)::-1] #reverse strand    
    compSeq = ''
    #get complement
    for n in seq:
        if n == 'G':
            compSeq += 'C'
        elif n == 'C':
            compSeq += 'G'
        elif n == 'A':
            compSeq += 'T'
        elif n == 'T':
            compSeq += 'A'    
    return compSeq

#method to delete duplicate proteins in a list of proteins
def deleteDupes(proteinlst):
    listNoDupes = []
    for p in proteinlst:
        if p not in listNoDupes:
            listNoDupes.append(p)
        else:
            pass
    return listNoDupes

#return longest protein in a list of proteins
def getLongestProtein(proteinList):
    longestProtein = proteinList[0]    
    for p in proteinList:                
        if len(p) > len(longestProtein):
            longestProtein = p
    return longestProtein

#write identified protein sequence to file
def writeProteinToFile(protein):
    print("\nEnter desired file name to store protein sequence in:")
    inputfile = input()
    inputfile = './sequence_files/protein/' + inputfile
    f = open(inputfile, "w")
    f.write(protein)
    f.close

#import DNA sequence
def importSequence():     
    inputfile = input("Enter text file containing DNA sequence:\n")    
    inputfile = "./sequence_files/genes/" + inputfile
    f = open(inputfile, "r")     
    for line in f:
        if line[0] == '>':
            variantsArray.append(line)                      
    f.seek(0)
    return f       

def identifyVariants(f):    
    return
    
if __name__ == '__main__':
    #import dna sequence
    f = importSequence()

    print("\nIdentified sequence variants:")
    index = 1
    for header in variantsArray:
        print(index, "-", header)
        index += 1

    print("Which variant would you like to translate?")
    print("Please enter the correspoding number: ", end='')
    variantNum = int(input())
    variant = variantsArray[variantNum-1]
    for line in f:
        if variant in line:
            break
         
    dnaSeq = ''    
    for line in f:  
        if line[0] is not '>':                  
            dnaSeq = dnaSeq + line 
        else:
            break    
               
    
    f.close()
    dnaSeq = dnaSeq.replace("\n", '')
    print("--------------------------------------------------------------")    
    
    #get comlementary dna sequence
    dnaSeqComp = getSequenceComplement(dnaSeq)         
    proteinList = []

    #generate the three possible reading for each strand (original and comp)
    directFrames = generateReadingFrames(dnaSeq)   
    for frame in directFrames:
        ORFs = getORFs(frame, proteinList)         

    complementFrames = generateReadingFrames(dnaSeqComp)    
    for frame in complementFrames:
        revORFs = getORFs(frame, proteinList)            

    #delete duplicates, print # of proteins and sequence of longest one
    proteinListNoDupes = deleteDupes(proteinList)
    print("Number of possible protein products found:", len(proteinListNoDupes), "\n")
    longestProtein = ''
    if len(proteinListNoDupes) > 0:
        longestProtein = getLongestProtein(proteinListNoDupes)
        print("Most likely protein product:")
        print(longestProtein)
        print("--------------------------------------------------------------")

        print("\nWould you like to save this sequence in a text file?")
        while 1:            
            answer = input()
            answer = answer.upper()
            if answer == "YES":
                writeProteinToFile(longestProtein)
                break
            elif answer == "NO":
                break
            else:
                print("Please enter 'Yes' or 'No'.\n")
        
    print("\nWould you like to compare your protein product with another?")
    while 1:        
        answer = input()
        answer = answer.upper()
        if answer == "NO":
            break
        elif answer == "YES":                
            print("Please enter the file path of the other protein:")
            otherProtFile = input()
            otherProtFile = "./sequence_files/protein/" + otherProtFile
            f = open(otherProtFile, "r")
            otherProt = f.read()
            compareProtein(longestProtein, otherProt)
            break
        else:
            print("Please enter 'Yes' or 'No'.\n")


        

	
