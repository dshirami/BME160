#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

class NucParams:
    rnaCodonTable = {
    # RNA codon table
    # U
    'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
    'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
    'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
    'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
    # C
    'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
    'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
    'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
    'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
    # A
    'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
    'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
    'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
    'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
    # G
    'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
    'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
    'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
    'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }
    dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

    def __init__ (self, inString=''):
        '''initialize object. MOstly dictionaries for other functions'''
        self.myInstring = inString.upper()
        self.nucComp = {}
        self.codonComp = {}
        pass

    def addSequence (self, inSeq):
        '''Add new sequence to existing sequence. This case just concats the
        strings but a list would be better but I am assuming that the input
        file will be correct and a mutiple of 3'''
        #self.myInstring.append(inSeq.upper())
        self.myInstring += inSeq
        self.__init__(self.myInstring)
        pass

    def aaComposition(self):
        '''create a dictionary that has the count of each 1 letter abreviation
        in the pretein'''
        self.aaComp = {
            'A': 0,  'G': 0,  'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0,  'F': 0, 'L': 0, 'R': 0, 'Y': 0, '-' : 0
            } #create a dict with all 20 and the stop '-' Amino acids
        codon = ''
        aaSeq = ''
        for codons in self.myInstring: #create a string of just the 1 let abrv.
            codon += codons
            if len(codon) == 3:
                if codon in NucParams.dnaCodonTable:
                    aaSeq += NucParams.dnaCodonTable[codon]
                elif codon in NucParams.rnaCodonTable:
                    aaSeq += NucParams.rnaCodonTable[codon]
                codon = ''
        for aa in aaSeq: #go through the 1 let abvr and add to count in the dict
            for aaKey in self.aaComp:
                if aa == aaKey:
                    self.aaComp[aa] += 1
        return self.aaComp

    def nucComposition(self):
        '''create a dict that has the count of ACGT or ACGU depending on DNA or
        RNA '''
        validNuc = 'ACGTU' #set what is a valid nucleotide
        for nuc in self.myInstring: #create or add to count of the nucleotide
            if nuc in validNuc: #check if its a valid nuc
                if nuc in self.nucComp: #try to add one if it exists
                    self.nucComp[nuc] += 1
                else: #create if not
                    self.nucComp[nuc] = 1
        return self.nucComp

    def codonComposition(self):
        '''create a dict that counts how many of the 3 letter codon is contained
        in the genome sequence'''
        tempCod = ''
        RNAonly = self.myInstring.replace('T','U') #change seq to be for RNA
        for code in RNAonly: #create the dictionary of the codons
            tempCod += code
            if len(tempCod) == 3: #every 3 bases add or create to the dictionary
                if tempCod in NucParams.rnaCodonTable:
                    if tempCod in self.codonComp:
                        self.codonComp[tempCod] += 1
                    else: #create if can't
                        self.codonComp[tempCod] = 1
                tempCod = ''
        return self.codonComp

    def nucCount(self):
        return sum(self.nucComp)

import sys
class FastAreader :
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen (self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta (self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''

        with self.doOpen() as fileH:

            header = ''
            sequence = ''

            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>') :
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

class ProteinParam :
# These tables are for calculating:
#     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
#     absorbance at 280 nm (aa2abs280)
#     pKa of positively charged Amino Acids (aa2chargePos)
#     pKa of negatively charged Amino acids (aa2chargeNeg)
#     and the constants aaNterm and aaCterm for pKa of the respective termini
#  Feel free to move these to appropriate methods as you like

# As written, these are accessed as class attributes, for example:
# ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093,  'G': 75.067,  'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225,  'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
        }

    mwH2O = 18.015
    aa2abs280= {'Y':1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R':12.4, 'H':6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__ (self, protein):
        '''
           initiate protein object
           called from main
        '''
        self.myProtein = protein
        return

    def aaCount (self):
        '''
           Counts how many of the characters are actually proteins in the string
        '''
        aaCounter = 0 #initiate counter
        for aa in self.myProtein.upper(): #iterate through protein string
            if aa in ProteinParam.aa2mw: #if its in the protein dictionary increase counter by 1
                aaCounter += 1
        return aaCounter

    def pI (self):
        '''
           Calculates the theoretical isometric point by seeing which pH has a charge closest to 0
        '''
        isoP = 0 #object to store the pI value
        lowestCharge = 100000 #object to compare the charge closest to 0
        pH = 0 #initiate pH object
        while(pH <= 14): #iterate through every pH every .01 to 14
            charge = self._charge_(pH)
            if(charge < 0): #absolute value to check for closeness to 0
                charge *= -1
            if (charge < lowestCharge): #if the charge is closer to 0/ less than the last one closest,
                lowestCharge = charge #make lowest charge the current charge
                isoP = pH #make isoP the current pH
            pH += .01
        return isoP

    def aaComposition (self) :
        '''
           Returns a dictionary that has the values of how many of each of the proteins are in in the string
        '''
        aaComp = {
            'A': 0,  'G': 0,  'M': 0, 'S': 0, 'C': 0,
            'H': 0, 'N': 0, 'T': 0, 'D': 0, 'I': 0,
            'P': 0, 'V': 0, 'E': 0, 'K': 0, 'Q': 0,
            'W': 0,  'F': 0, 'L': 0, 'R': 0, 'Y': 0
            } #initiate dictionary that hold the values for counting
        for aaKey in aaComp: #iterate through dictionary
            count = 0 #reset count for every aa
            for proKey in self.myProtein: #iterate through the protein string
                if aaKey == proKey: #if the object in aaComp is equal to the object in the protein string
                    count += 1 #increase count by 1
            aaComp[aaKey] = count #set the dictionary object value to the count value
        return aaComp # return the dictionary

    def _charge_ (self, pH):
        '''
           method which has the math behind finding the net charge using pH and aa charges
        '''
        posNet = 0 #initiate pos and neg objects
        negNet = 0
        aaComp = self.aaComposition() # get the composition of the protein stand
        for posCharge in ProteinParam.aa2chargePos: #find positive charge
            posNet += aaComp.get(posCharge)*(10**ProteinParam.aa2chargePos[posCharge] / (10**ProteinParam.aa2chargePos[posCharge] + 10**pH))
        for negCharge in ProteinParam.aa2chargeNeg:#find negative charge
            negNet += aaComp.get(negCharge)*(10**pH / (10**ProteinParam.aa2chargeNeg[negCharge] + 10**pH))
        posNet += 10**ProteinParam.aaNterm / (10**ProteinParam.aaNterm + 10**pH)
        negNet += 10**pH / (10**ProteinParam.aaCterm + 10**pH)
        return posNet-negNet

    def molarExtinction (self):
        '''
           finds the extinction coefficient which indicates how much light a protein absorbs at a certain wavelength
        '''
        aaComp = self.aaComposition() #composition of the protein strand
        extinc = aaComp.get('Y')*1490 + aaComp.get('W')*5500 + aaComp.get('C')*125 #equation to find the molar extinction
        return extinc

    def massExtinction (self):
        '''
           calculate mass extinction using molar Extinction
        '''
        myMW =  self.molecularWeight() #get molar weight to divide the molar extinction by
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight (self):
        '''
           find molecular weight using composition and the molecular weights
        '''
        molWeight = 0
        for aa in self.myProtein: #iterate through protein strand
            molWeight += ProteinParam.aa2mw.get(aa) #add the molecular weight of the aa to the total weight
        molWeight -= ((self.aaCount()-1)*ProteinParam.mwH2O)
        return molWeight

class OrfFinder:
    '''
    This class finds Open Reading frames in a sequence
    when findORF is called. This can give All or the longest ORF
    in the sequence
    '''

    def findTopFrame(position):
        '''function to find what frame the codon is in for top strand'''
        frame = (position%3)+1
        return frame

    def appendTopORFs(ORFList, ORF):
        '''function to append found ORFs from top stand to a ORFList
        made it so that in ORFList - [[frame1],[f2],[f3],[f-3],f[-2],[f-1]]'''
        ORF = ORF + [ORF[1]-ORF[0]+1] + [OrfFinder.findTopFrame(ORF[1])]
        #I added 2 more arguments to the ORF object
        #[left position, right pos, length, frame]
        if OrfFinder.findTopFrame(ORF[1]) == 1:
            if ORF not in ORFList[0]:#so u dont repeat ORFs
                ORFList[0].append(ORF)
        if OrfFinder.findTopFrame(ORF[1]) == 2:
            if ORF not in ORFList[1]:
                ORFList[1].append(ORF)
        if OrfFinder.findTopFrame(ORF[1]) == 3:
            if ORF not in ORFList[2]:
                ORFList[2].append(ORF)

    def findBotFrame(length, position):
        '''function to find what frame the codon is in for bot strand'''
        frame = -((length-position)%3+1)
        return frame

    def appendBotORFs(ORFList, length, ORF):
        '''function to append found ORFs from bot stand to a ORFList
        made it so that in ORFList - [[frame1],[f2],[f3],[f-3],f[-2],[f-1]]'''
        ORF = ORF + [ORF[1]-ORF[0]+1] + [OrfFinder.findBotFrame(length, ORF[1])]
        #I added 2 more arguments to the ORF object
        #[left position, right pos, length, frame]
        if OrfFinder.findBotFrame(length, ORF[1]) == -1:
        #made it so that in ORFList - [[frame1],[f2],[f3],[f-3],f[-2],[f-1]]
            if ORF not in ORFList[-1]:
                ORFList[-1].append(ORF)
        if OrfFinder.findBotFrame(length, ORF[1]) == -2:
            if ORF not in ORFList[-2]:
                ORFList[-2].append(ORF)
        if OrfFinder.findBotFrame(length, ORF[1]) == -3:
            if ORF not in ORFList[-2]:
                ORFList[-3].append(ORF)
        return ORFList

    def __init__(self, seq, long):
        '''
        initiate seq object and if longest(true/false) based on commandline
        input
        '''
        self.mySeq = seq
        self.longest = long

    def findORF(self):
        '''
        This is the main function to find the ORFs in a Seqeuence
        '''
        #initialize ORFList ordered as [[frame1],[f2],[f3],[f-3],[f-2],[f-1]]
        ORFList = [[],[],[],[],[],[]]
        for pos in range(0,len(self.mySeq)):
            #iterate through every position in the sequence
            codon = self.mySeq[pos:pos+3]
            #set codon to be the 3 letter codon from where youre checking
            ORF = [1, len(self.mySeq)] #set ORF to a default(entire thing)
            if codon == 'ATG':#look for starts
                ORF[0] = pos+1#set start position
                #check frame and create a list of codons based on frame
                if OrfFinder.findTopFrame(pos) == 1:
                    newFrame = [self.mySeq[p:p+3] for p in
                            range(pos,len(self.mySeq),3)]
                if OrfFinder.findTopFrame(pos) == 2:
                    newFrame = [self.mySeq[0]] + [self.mySeq[p:p+3] for p in
                        range(pos,len(self.mySeq),3)]
                if OrfFinder.findTopFrame(pos) == 3:
                    newFrame = [self.mySeq[0:2]] + [self.mySeq[p:p+3] for p in
                        range(pos,len(self.mySeq),3)]
                #in the new framed list look for a stop codon
                for startCheck, findStart in enumerate(newFrame, start=1):
                    if findStart=='TAA' or findStart=='TAG' or findStart=='TGA':
                        #when u find a stop set right position to ORF
                        if OrfFinder.findTopFrame(pos) == 1:
                            ORF[1] = startCheck*3+pos
                        if OrfFinder.findTopFrame(pos) == 2:
                            ORF[1] = startCheck*3+pos-3
                        if OrfFinder.findTopFrame(pos) == 3:
                            ORF[1] = startCheck*3+pos-3
                        #add it to ORFList and stop checking for a stop
                        OrfFinder.appendTopORFs(ORFList, ORF)
                        break
                #check if the start doesnt have a stop and is a dangling case
                newFrame = [self.mySeq[p:p+3] for p in range(pos,len(self.mySeq),3)]
                #make a frame and check if a stop is in the list
                if 'TAA' not in newFrame and 'TAG' not in newFrame and 'TGA' not in newFrame:
                    OrfFinder.appendTopORFs(ORFList, ORF)
            #check if theres a stop but no start (hanging case at the front)
            if codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
                #check stops to see if they have a start before them
                if OrfFinder.findTopFrame(pos) == 1:
                    newFrame = [self.mySeq[p:p+3] for p in
                            range(0,pos+3,3)]
                if OrfFinder.findTopFrame(pos) == 2:
                    newFrame = [self.mySeq[0]] + [self.mySeq[p:p+3] for p in
                        range(1,pos+3,3)]
                if OrfFinder.findTopFrame(pos) == 3:
                    newFrame = [self.mySeq[0:2]] + [self.mySeq[p:p+3] for p in
                        range(3,pos+3,3)]
                if 'ATG' not in newFrame:
                    #if it doesnt have a start add ORF from front-stop
                    ORF[0] = 1
                    ORF[1] = pos+3
                    OrfFinder.appendTopORFs(ORFList, ORF)
        #Bot strand
            #reset ORF so you can use it again
            ORF = [1,len(self.mySeq)]
            #look for stops but BACKWARDS because its the bot strand
            if codon == 'TTA' or codon == 'CTA' or codon == 'TCA':
                #if u find a stop make it the start position
                ORF[0] = pos+1
                #make a new list in the right frame
                newFrame = [self.mySeq[p:p+3] for p in range(pos,len(self.mySeq),3)]
                #look through the list and keep track how far into the list you
                #went
                for stopCheck, findStop in enumerate(newFrame[1:], start=2):
                    if findStop=='TTA' or findStop=='CTA' or findStop=='TCA':
                        #if u find a stop in the list its NOT a ORF so exit the
                        #loop/ its a new ORF so look at the next time around
                        break
                    if findStop == 'CAT':
                        #if u find a start, set the right postion
                        ORF[1] = stopCheck*3+pos
                        OrfFinder.appendBotORFs(ORFList, len(self.mySeq), ORF)
                        #keep looking until you hit another STOP because that
                        #gets all of the ORF in this ORF not just the shortest
                #hanging case
                if 'CAT' not in newFrame:
                    #if there isnt a STOP then add [stop-end] to ORFList
                    ORF = ORF + [ORF[1]-ORF[0]+1] + [OrfFinder.findBotFrame(len(self.mySeq), ORF[0]+2)]
                    ORFList[ORF[3]].append(ORF)
            #reset ORF so you can use it again
            ORF = [1, len(self.mySeq)]
            #check for hanging ORF at the begining
            if codon == 'CAT': #check every start codon
                ORF[1] = pos+3
                #make a frame based on the start codon
                newFrame = [self.mySeq[p:p+3] for p in range(0,pos+2,3)]
                #if there isnt a STOP BEFORE the start codon add it to ORFList
                #as an ORF from [1-start]
                if 'TTA' not in newFrame and 'CTA' not in newFrame and 'TCA' not in newFrame:
                    ORF = ORF + [ORF[1]-ORF[0]+1] + [OrfFinder.findBotFrame(len(self.mySeq), ORF[0]+2)]
                    ORFList[ORF[3]].append(ORF)
        '''
        The next 2 for loops are for the case that nothing was found in a
        specific frame for the sequence. It adds an ORF of the length of the
        entire strand to ORFList at that Frame.
        '''
        for frame in range(0,3):
            if not ORFList[frame]:
                ORFList[frame].append([1,len(self.mySeq), len(self.mySeq), frame+1])
        for frame in range(-3,0):
            if not ORFList[frame]:
                ORFList[frame].append([1, len(self.mySeq), len(self.mySeq), frame])
        #now sorting and returning the List of ORFs
        newORFList = []
        '''if -lG is in the Command Line then I want a list with just the longest
        ORF. I did this by utilizing how I added ORFs to ORFList. I can just
        compare the ORFList[1] and ORFList[2](top strand) and
        ORFList[0] and ORFList[2](bot strand) and take the highest/longest ORF
        from each ORF (ORFList[2] within ORFList[0 or 1]) and make a new list
        with those
        '''
        if self.longest:
            #for top frames
            for listORF in ORFList[0:3]:
                lis = []
                #sort by left position and then by length of the position
                l = sorted(listORF, key=lambda c:[c[1],c[2]], reverse=True)
                for ORF in l:#iterate through newly sorted list
                    #check if the ORF has been added already
                    if ORF[1] not in lis:
                        newORFList.append(ORF)
                    #once you add ORF, mark it so you dont add a smaller ORF
                    lis.append(ORF[1])
            #for bot frames
            for listORF in ORFList[3:]:
                lis = []
                #sort by left position and then by length of the position
                l = sorted(listORF, key=lambda c:[c[0],c[2]], reverse=True)
                for ORF in l:
                    #check if the ORF has been added already
                    if ORF[0] not in lis:
                        newORFList.append(ORF)
                    #once you add ORF, mark it so you dont add a smaller ORF
                    lis.append(ORF[0])
        #if the use wants every ORF make a new list from ORFList
        if self.longest == False:
            for orf in ORFList:
                for o in orf:
                    newORFList.append(o)
        #sort list of ORFs based on Length and then left postion
        s = sorted(newORFList, key=lambda c:[-c[2],c[0]])
        return s
