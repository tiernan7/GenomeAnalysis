from Bio import SeqIO
import re
import numpy as np
import math
import os
import timeit


#A function writes certain program checks to a log file
def log(text):
    with open("log","a") as logfile:
        logfile.write(text)


start = timeit.default_timer()
np.set_printoptions(suppress=True)

#  record = SeqIO.read("Armadillidium_vulgare_contig1.gbff","genbank")

                         
f = "Armadillidium_vulgare_contig1.gbff" #filename
o = "TIERNAN_KENNEDY_HOMEWORK3.txt" #output filename

#Makes it easier to develope when writing in append mode
if os.path.exists(o):
    os.remove(o)

#Using biopython to extract the sequence from the genbank file
record = SeqIO.read(f,"genbank")
seq = str(record.seq)


#Counts the number of occurences of a character in a string
def countNucleotide(seq, base):
    count = 0
    for i in range(len(seq)):
        if seq[i] == base:
            count += 1
    return count

log("The are are " + str(countNucleotide("TGATGGTATGNNTGATGK",'G')) + " Gs in TGATGGTATGNNTGATGK" + '\n'
)

#The four bases in DNA
bases = ['A', 'C', 'G', 'T']

#A dictionary for finding complements
baseDict = {'A' : 'T',
            'C' : 'G',
            'T' : 'A',
            'G' : 'C',
            'N' : 'N'}

#A hashing of the characters to ints
baseToNum = {'A' : 0,
             'C' : 1,
             'G' : 2,
             'T' : 3}

#Returns the reverse complement of a string
def reverseComplement(s):
    return "".join(np.flip([baseDict[char] for char in s]))

log("The reverse complement of T is: " + str(reverseComplement('T')) + '\n'
    )


#Prints the nucleotide counts and frequencies
def getNucleotideHistogram(sequence):
    nucleotideCounts = np.zeros(5)
    for base in bases:
        nucleotideCounts[baseToNum[base]] = countNucleotide(sequence, base)
    s = sum(nucleotideCounts)
    l = len(sequence)
    returnString = ""
    returnString += "Background Counts" + '\n'
    returnString += "A=" + str(int(nucleotideCounts[0])) + '\n'
    returnString += "C=" + str(int(nucleotideCounts[1])) + '\n'
    returnString += "G=" + str(int(nucleotideCounts[2])) + '\n'
    returnString += "T=" + str(int(nucleotideCounts[3])) + '\n'
    returnString += "N=" + str(int(l - s)) + '\n'
    returnString += '\n'
    returnString += "Background Frequencies" + '\n'
    returnString += "A=" + ("%.4f" % round(float((nucleotideCounts[0] + nucleotideCounts[3])) / (2*s), 4)) + '\n'
    returnString += "C=" + ("%.4f" % round(float((nucleotideCounts[1]+nucleotideCounts[2] )) / (2*s), 4)) + '\n'
    returnString += "G=" + ("%.4f" % round(float((nucleotideCounts[2]+nucleotideCounts[1])) / (2*s), 4)) + '\n'
    returnString += "T=" + ("%.4f" % round(float((nucleotideCounts[3] +nucleotideCounts[0]))/ (2*s), 4)) + '\n'

    return returnString

log("The nucleotide histogram for ATGCTCANNNNNNNKJ is: " + str(getNucleotideHistogram('ATGCTCANNNNNNNKJ')) + '\n'
)

#Opening the genbank file and reading it into a list of lines
file = open(f, "r")
allLines = file.readlines()
file.close()


#Removes all empty characters from a string
def removeEmptyCharacters(array):
    return [x for x in array if x != '']

log("The following text is on one line: " + str(removeEmptyCharacters("Line 1 \n Line 2")) + '\n'
)



#Takes in an array of words from a line of a genbank file
#Returns True if the any of the words contain an ambiguous coordinate signature
def badCoordinates(wordsArray):
    return (re.search("<+|>+", wordsArray[-1]) != None)
log("It is " + str(badCoordinates(["10002000>303"])) + " that 10002000>303 is a bad coordinate" + '\n'
)


#Takes in a CDS description containing a join statment
#Returns an ordered list of the bases around the start of the CDS
def brokenIntervals(intervals, comp): #Intervals is a list of each integer within a join statement in the orer of appearance, records whether or not the complemet
    if not(comp): #Handles the case where the join intervals are the on forward strand
        returnIterator = [] 
        start = int(intervals[0]) - 1 #-1 adjustment tp translate genbank coordinates to array coordinates
        returnIterator = np.arange(start- 10 , start + 1) #generates the sites of the ten bases before and upto the start
        ranges = [int(i) - 1 for i in intervals] #converts the entire list of intervals to array coordinates
        forwardLength = [] 
        ind = 1#distance from start site of current interval
        r = 0#index of which set of intervals to consider
        start = ranges[0]
        while len(forwardLength) < 10: #Generates the ten bases after the start
            if (start + ind > ranges[r+1]): #If you're pointing to an adress out of the current interval...
                r += 2
                start = ranges[r]#...start at the next site
                ind = 0#Reset the distance of the pointer from the current start site
            forwardLength += [start + ind]#Add the array coordinates of that site in the interval to the output array
            ind+=1
        returnIterator = np.append(returnIterator, forwardLength)#merge the two halves of the sequence
        return returnIterator

    else:#Handles the case where the join intervals are on the reverse strand     
        returnIterator = []
        a = np.flip([int(i) - 1 for i in intervals])#Reverses the list of coordinates to account for navigating the interval in reverse and subtraction to convert to ar#ray coordinates
        start = int(a[0])#start at the array coordinate of the last element in forward strand on the coding interval
        returnIterator = np.flip(np.arange(start , start + 11))#generates the coordinat running from 11 before the start site to the start site
        forwardLength = []
        ind = -1#The distcance from the start site
        r = 0#The index of the interval in question
        while len(forwardLength) < 10:#While you havent added enough stuff
            if (start + ind < a[r+1]):#check to see if youre going to leave the interval
                r += 2 
                start = a[r] #if so go to the next interval
                ind = 0
            forwardLength += [start + ind]#Add the array coordinates of the site to the sequence
            ind+= -1
        returnIterator = np.append(returnIterator, forwardLength)#Add the second half of the sequence to the first
        return returnIterator

log("The array sites for complement(join(1..10,11..30)) is " + str(brokenIntervals([1,10,11,30], True)) + '\n'
)

#Returns the column of an array
def column(arr, i):
    col = []
    for st in arr:
        col += [st[i]]
    return col
       

log("The column second column of [[1,2,3],[4,5,6],[7,8,9]] is "  + str(column([[1,2,3],[4,5,6],[7,8,9]],1)) + '\n'
)

#Takes a string text from a CDS description, but not multiple lines
def positionToChars(positionString):
    if (re.search('complement', positionString) == None):# if there are no complement tags
        if (re.search('join', positionString) == None): #and no join tags
            start = int(re.findall("\d+", positionString)[0]) - 1 #then start site is the first numeral in the position string minus 1 to convert to array coordinates
            return [seq[i] for i in np.arange(start - 10, start + 11)]#and we can return the ten elements on either side
        else: #otherwise if there are no complement tags but there is a join statement
            intervals = re.findall("\d+", positionString)# convert the string to a flat list of position numeral for the brokenIntervals function
            return [seq[i] for i in brokenIntervals(intervals, False)] #Use the broken intervals helper function to find the list of coordinates
    else:#otherwise if there are complement tags
        if (re.search('join', positionString) == None):#but no join statements
            start = int(re.findall("\d+", positionString)[-1]) - 1 #start at the last numeral in the CDS description, and subtract 1 for array coordinates
            return[baseDict[seq[i]] for i in np.flip( np.arange(start - 10, start + 11))]#return the reversed order of the ten elements around the start site

        else:# and if there are no join tags but the start is on the complementary strand
            intervals = re.findall("\d+", positionString)# get the list of numeraks in the CDS desription
            return [baseDict[seq[i]] for i in brokenIntervals(intervals, True)]#Use the broken intervals function to generates, which are reversed automatically by roken intervals
log("The characters for (1..10) is :" + str(positionToChars("(1..10)")) + '\n'
    )

#Takes in the text from any line containing CDS and extracts the portion of it containing the description and finds any overflow lines
positions = []
for i, lines in enumerate(allLines):
    cdsCheck = (re.search("CDS", lines) != None)#gets all lines containing "CDS"
    if cdsCheck and re.search('[0-9]',lines) != None:#but only takes those with numerals (To avoid the case when the 'CDS' protein sequence is accidentally parsed)
        words = re.split("\s", lines) #Split the line into words
        words = removeEmptyCharacters(words)#and remove empty characters
        if (badCoordinates(words)):#check to see if there are any ambiguous start or end sites
            pass #and if so ignor these
        else: #If the coorinates are okay
            delete = False
            if (re.search("complement|join", words[-1]) != None and words[-1][-1] != ")"):#checks to see if there are any descs that run onto two lines
                ind = 1
                delete = False
                while words[-1][-1] != ")": #A loop that continues until the end of the description has been found, start by checking first line
                    words[-1] += allLines[i+ind].strip()#adds the next line of sequence
                    if re.search(">|<", allLines[i+ind]) != None:#and marks it in case there was a variable end position that wasnt caught in the first check
                        delete = True 
                    i+=1#Goes to the next line
            positions += [words[-1]]#Adds the full description to the position array
            if delete:
                positions.pop(-1)#Deletes the position if it was marked as such



backgroundCounts = np.zeros(4)
backgroundFrequencies = np.zeros(4)

#counts the number of nucleotides of A,T,C and G in the genbank sequence
for base in bases:
    backgroundCounts[baseToNum[base]] = countNucleotide(seq,base)

for b  in bases:#Converts the frequencies of nucleotides to their ratios 
    backgroundFrequencies[baseToNum[b]] = float(backgroundCounts[baseToNum[b]] + backgroundCounts[baseToNum[baseDict[b]]]) / (2* sum(backgroundCounts)) #nucleotide counts are averaged with the sum of their complement to account for reverse stran


allSequences = [positionToChars(pos) for pos in positions]#Creating an array of windows of ten nucleotides on either side of each CDS feature

siteCounts = np.zeros((21,4))
siteFrequencies = np.zeros((21,4))
weightMatrix = np.zeros((21,4))


#counts the number of nucleotides at each position in each sequence
for i in range(len(siteCounts)):
    for base in bases:
        siteCounts[i][baseToNum[base]] = countNucleotide(column(allSequences, i),base )

#converts the counts from above to frequencies
for i in range(len(siteFrequencies)):
    s = sum(siteCounts[i])
    for j in range(4):
           siteFrequencies[i][j] = siteCounts[i][j] / s

#Perform a small correction and calculate weights for each base at each site
for i in range(len(siteFrequencies)):
    for base in bases:
        if siteFrequencies[i][baseToNum[base]] == 0:
            weightMatrix[i][baseToNum[base]] = -99.0
        else:
            weightMatrix[i][baseToNum[base]] = np.log2(siteFrequencies[i][baseToNum[base]] / backgroundFrequencies[baseToNum[base]])


#A function to score a sequence of 21 nucleotides against the weight matrix
def scoreSequence(s):
    assert len(s) == 21, "The sequence: " + str(s) + " does not have length 21"
    score = 0
    for i,char in enumerate(s):
        if char in bases:
            score += weightMatrix[i][baseToNum[char]] #Compute the dot product of the sequence against the weight matrix
        else: 
            score += 0
    return score

log("The score of NNNNNNNNNNNNNNNNNNNNA is " + str(scoreSequence("NNNNNNNNNNNNNNNNNNNNA")) + '\n'
    )

#A function that returns a histrogram of the scores of each CDS site as a string
def scoreHistogramCDS():
    returnString = ""
    returnString += "Score Histogram CDS:" + '\n'
    scores = []

    for s in allSequences:#For each sequence we collected from the CDS sites, we add its score, rounded down, to a list
        scores += [math.floor(scoreSequence(s))]

      
    returnString += "lt-50 : " + str(len([i for i in scores if i < -50])) + '\n'#for each bin in the histogram we...
    for i in range(-50,100):
        count = len([j for j in scores if j == i]) #count the number of each score in the list and..
        if count != 0:
            returnString += str(i) + " " + str(count) + '\n' #Add that row to the histogram if the bin is nonempty
    returnString += '\n'
      
    return returnString


#A function that computes outputs the histogram of scores for each sequence in the complex and also returns the list of high scorers
def scoreAll():
    returnString = ""
    returnString += "Score Histogram All" + '\n'
    scores =[]
    l = len(seq) 
    stretches = [seq[i:i+21] for i in range(len(seq) - 20)]


    for i, s in enumerate(stretches):
        scores += [[scoreSequence(s), i + 1 + 10, 0]]
        scores += [[scoreSequence(reverseComplement(s)), i + 1 + 10, 0]]


  

    
        

    returnString += "lt-50 : " + str( len([math.floor(i[0]) for i in scores if math.floor(i[0]) < -50])) + '\n'
    for i in range(-50,100):
        count = len([j[0] for j in scores if math.floor(j[0]) == i])
        if count != 0:
            returnString += str(i) + " " + str(count) + '\n'
    returnString += '\n'

    returnString += "Position List:"
    matches = np.array([sc for sc in scores if sc[0] >= 10])
    sortedMatches = matches[np.argsort(column(matches,1))]
    for m in sortedMatches:
        returnString += str(m[1]) + " " + str(m[2]) + " " + str("%.4f" % round(m[0], 4)) + '\n'
    returnString += '\n'

    return returnString


with  open(o, 'a') as file:
    file.write(getNucleotideHistogram(seq))
    file.write('\n')

    file.write("Count Matrix:" + '\n')
    for i in range(len(siteCounts)):
        file.write(str(i-10) + " ")
        for j in range(4):
            file.write(str(round(siteCounts[i][j])) + " ")
        file.write('\n')
    file.write('\n')

    file.write("Frequency Matrix:" + '\n')
    for i in range(len(siteFrequencies)):
        file.write(str(i-10) + " ")
        for j in range(4):
            file.write(str("%.4f" % round(siteFrequencies[i][j],4)) + " ")
        file.write('\n')
    file.write('\n')

    file.write("Weight Matrix:" + '\n')
    for i in range(len(weightMatrix)):
        file.write(str(i-10) + " ")
        for j in range(4):
            file.write(str("%.4f" % round(weightMatrix[i][j],4)) + " ")
        file.write('\n')
    file.write('\n')

    file.write("Max score: " + str("%.10f" % round(sum([max(i) for i in weightMatrix]),10)) + '\n')
    file.write('\n')

    file.write(scoreHistogramCDS() + '\n')
    file.write('\n')
    
    file.write(scoreAll())
    file.write('\n')
    file.write("Program:" + '\n')

stop = timeit.default_timer()
print("Runtime: " + str(stop - start) + "s")

    
        


