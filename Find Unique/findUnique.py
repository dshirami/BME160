#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

import sys
class FastAreader :

    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):

        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
                if not line: # we are at EOF
                    return header, sequence
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith ('>'):
                    yield header,sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else :
                    sequence += ''.join(line.rstrip().split()).upper()

        yield header,sequence

########################################################################
# Main
# Here is the main program
#
########################################################################
def main(inCL=None):
    '''
    findUniques takes in a filename from the commandline and prints
    a report of the unique genes in the tRNA mitochondria
    '''
    myReader = FastAreader()
    tRNAs = [] #this will be where I store a list of ala 22 [head][seq]
    for head, seq in myReader.readFasta():#get head and seq from file
        tRNA = []
        tRNA.append(head.replace(' ','')) #get head without white spaces
        stripSeq = seq.replace('.','-')
        stripSeq = stripSeq.replace('_','-')
        tRNA.append(stripSeq.replace('-','')) #get seq without "._-"
        tRNAs.append(tRNA) #add it to where I keep them all
    tRNAs.sort() #sort by object sort/ alhabetically

    subStrings = [] #list of each tRNA substring set
    for seqs in tRNAs: #for all 22 seqs
    #create all possible substrings
        newSet = []
        for pos in range(0, len(seqs[1])):
            for count in range(pos+1, len(seqs[1])+1):
                newSet.append(seqs[1][pos:count])
        subStrings.append(set(newSet))#add it to the subStrings list
        #all these lists will be ordered to the same position
        #that is how i keep track what is where

    allOtherSets = [] #set of all other substrings
    for seqs in subStrings:
    #create union of all other tRNA sets
        unionSets = set()
        for sets in range(0, len(subStrings)):
            #iterate through everything except the set we're on
            if seqs != subStrings[sets]:
                unionSets = unionSets.union(subStrings[sets])
        allOtherSets.append(unionSets)
        #add it to the list of all other substrings

    for pos in range(0,len(subStrings)):
        subStrings[pos] = subStrings[pos].difference(allOtherSets[pos])
        #using allOtherSets we take the difference from subStrings to get
        #the unique substrings in each

    for sets in range(0, len(subStrings)): #for each tRNA substring
        #find the extensions of the unique substrings
        extensions = [] #keep track of what to take out
        for look in subStrings[sets]: #nested for loop to compare
            for comp in subStrings[sets]:
                if look in comp and look is not comp:
                    #is a part of but not exactlty
                    extensions.append(comp)
        extensions = set(extensions) #make it a set so no repeats
        subStrings[sets] = subStrings[sets].difference(extensions)
        #take diffference

    for sets in range(0, len(subStrings)):
        subStrings[sets] = list(subStrings[sets])
        #make the subStrings a list so we can sort and iterate through

    '''
    so to sort the substrings I made each element in substrings
    to be a list for example:

    subStrings[0] = [[where it belongs,unique substring],[]]

    each element would now be a list with the 0 position being how many
    dots/where its found and 1 position being the unique substrings
    '''
    for pos, sets in enumerate(subStrings):
        for num,findPos in enumerate(sets):
            posUniq = []
            numDots = tRNAs[pos][1].find(findPos)
            #using find() to see where it belongs
            posUniq.append(numDots)
            #add the two elements in this order to a temp list
            posUniq.append(findPos)
            subStrings[pos][num] = posUniq #reassign each element in subStrings
        subStrings[pos] = sorted(subStrings[pos], key=lambda x: x[0])
        #sort it based on the number/position of the substrings

    for pos, tRna in enumerate(tRNAs):
        #print the report
        print(tRna[0])
        print(tRna[1])
        for uniq in subStrings[pos]:
            print('{}{}'.format('.'*uniq[0], uniq[1]))

    pass
if __name__ == "__main__":
    main()

