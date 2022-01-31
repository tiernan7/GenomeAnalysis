import numpy as np
import os
import timeit

start = timeit.default_timer()

baseDict = {'A' : 0,
            'C' : 1,
            'G' : 2,
            'T' : 3}

bases  = ['A', 'C', 'G', 'T']

dibases = ['AA', 'AC', 'AG', 'AT','CA', 'CC', 'CG', 'CT','GA', 'GC', 'GG', 'GT','TA', 'TC', 'TG', 'TT']

def countNucleotide(seq, base):
    count = 0
    for i in range(len(seq)):
        if seq[i] == base:
            count += 1
    return count

def countDinucleotide(seq, dibase):
    count = 0
    for i in range(0, len(seq) - 1):
        if seq[i:i+2] == dibase:
            count += 1
    return count



  
def dinucleotideToConditional(matrix):
    returnArray = np.zeros(16)
    for i in range(4):
        summ = 0
        for j in range(4):
            summ += matrix[4*i + j]
        for j in range(4):
            returnArray[4*i + j] = (matrix[4*i + j]/ summ)
    return returnArray

          
    
class Sequence:
    def __init__(self, name):
        self.nucleotideCounts = []
        self.nucleotideFrequencies = []
        self.dinucleotideCounts = []
        self.dinucleotideFrequencies = []
        self.conditionalFrequencies = []
            
        returnString = ""
        with open(name) as f:
            txt  = f.read().splitlines()  
        self.tagline = txt[0]
        for i in range(1, len(txt)):
            returnString += str.rstrip(txt[i])
        self.seq = returnString
        self.l = len(self.seq)
        self.fileName = name

        for i in range(len(bases)):
            self.nucleotideCounts += [countNucleotide(self.seq, bases[i])]
        for i in range(len(self.nucleotideCounts)):
            self.nucleotideFrequencies += [float(self.nucleotideCounts[i]) / float(sum(self.nucleotideCounts))]
    
        

        

        for i in range(0, len(dibases)):
            self.dinucleotideCounts += [countDinucleotide(self.seq, dibases[i])]
        for i in range(0, len(self.dinucleotideCounts)):
            self.dinucleotideFrequencies += [float(self.dinucleotideCounts[i])/ float(sum(self.dinucleotideCounts))]


                    
        self.conditionalFrequencies = dinucleotideToConditional(self.dinucleotideFrequencies)

         
    
    def getNucleotideHistogram(self):
        returnString = ""
        returnString += "*=" + str(self.l) + '\n'
        returnString += "A=" + str(self.nucleotideCounts[0]) + '\n'
        returnString += "C=" + str(self.nucleotideCounts[1]) + '\n'
        returnString += "G=" + str(self.nucleotideCounts[2]) + '\n'
        returnString += "T=" + str(self.nucleotideCounts[3]) + '\n'
        returnString += '\n'

        return returnString

    def summarizeSequence(self):

        returnString = ""

        returnString += "Fasta: " + self.fileName  + '\n'
        returnString += self.tagline + '\n'
        returnString += self.getNucleotideHistogram()

        returnString += "Nucleotide Frequencies:" + '\n'
        for char in bases:
            returnString += char + "=" + str("%.4f" % round(self.nucleotideFrequencies[baseDict[char]],4)) + '\n'
        returnString += '\n'

        
        returnString += "Dinucleotide Count Matrix:" + '\n'
        for i in range(0,16):
            if (i % 4 == 0):
                returnString += dibases[i][0] + "="
            returnString += str(self.dinucleotideCounts[i]) + " "
            if (i in [3,7,11,15]):
                returnString += '\n'
        returnString += '\n'

        returnString += "Dinucleotide Frequency Matrix:" + '\n'
        for i in range(0,16):
            if (i % 4 == 0):
                returnString += dibases[i][0] + "="
            returnString += str("%.4f" % round(self.dinucleotideFrequencies[i], 4)) + " "
            if (i in [3,7,11,15]):
                returnString += '\n'
        returnString += '\n'

        returnString += "Conditional Frequency Matrix:" + '\n'
        for i in range(0,16):
            if (i % 4 == 0):
                returnString += dibases[i][0] + "="
            returnString += str("%.4f" % round(self.conditionalFrequencies[i], 4)) + " "
            if (i in [3,7,11,15]):
                returnString += '\n'
        returnString += '\n'


        return returnString


def zerothOrderMarkov(seq):
    probArray = seq.nucleotideFrequencies
    return np.random.choice(bases, p = probArray)

        


def firstOrderMarkov(seq, prevBase):
    if prevBase == "":
        return np.random.choice(dibases, p = seq.dinucleotideFrequencies)
    else:
        row = baseDict[prevBase]
        probArray = seq.conditionalFrequencies[4*row : 4*row + 4]
        return np.random.choice(bases, p = probArray)



mouseGenome = Sequence("mm10_chr12_100Mb_110Mb.fa")


outfile1 = open("o0MM" + mouseGenome.fileName, "a")
outlines = []
outlines += (">" + "0MM_Simulated_" + mouseGenome.tagline[1:]) + '\n'
leng = 0
seqlen = mouseGenome.l
nextline = ""
while leng < seqlen:
    nextline += zerothOrderMarkov(mouseGenome)
    if (leng % 60 == 0 and leng != 0):
        outlines += [nextline + '\n']
        nextline = ""
    leng += 1
 
outlines += [nextline]
outfile1.writelines(outlines)
outfile1.close()


outfile2 = open("o1MM" + mouseGenome.fileName, "a")
outlines = []
outlines += (">" + "1MM_Simulated_" + mouseGenome.tagline[1:]) + '\n'
leng = 0
seqlen = mouseGenome.l
nextline = ""
nextline += firstOrderMarkov(mouseGenome, "")
leng += 2
last = nextline[-1]
while leng < seqlen:
    nextline += firstOrderMarkov(mouseGenome, last)
    if (leng % 60 == 0 and leng != 0):
        last = nextline[-1]
        outlines += [nextline + '\n']
        nextline = ""
    else:
        last = nextline[-1]
    leng += 1

outlines += [nextline]
outfile2.writelines(outlines)
outfile2.close()


simulated1 = Sequence("o0MM" + mouseGenome.fileName)
simulated2 = Sequence("o1MM" + mouseGenome.fileName)

print(simulated1.l)
print(simulated2.l)

outfile = open("TIERNAN_KENNEDY_HOMEWORK2.txt", "w")
outlines = []

outlines += mouseGenome.summarizeSequence()
outlines += simulated1.summarizeSequence()
outlines += simulated2.summarizeSequence()

outfile.writelines(outlines)
outfile.close()

stop = timeit.default_timer()

print("Runtime: " + str(stop - start) + "s")

