#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

import math
import sys
import os
class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''

    def __init__(self, inOpts=None) :
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(description = 'Translates Quality lines to desired PHRED mapping',
                                             epilog = '',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s  <input >output [options] option1[default]'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('outFile', action = 'store', help='output file name')
        self.parser.add_argument('-P33in', '--PHRED33input', action = 'store', nargs='?', const=True, default=False)
        self.parser.add_argument('-P64in', '--PHRED64input', action = 'store', nargs='?', const=True, default=False)
        self.parser.add_argument('-P64Bin', '--PHRED64Boffin', action = 'store', nargs='?', const=True, default=False)
        self.parser.add_argument('-P64SOLin', '--PHRED64SOLin', action = 'store', nargs='?', const=True, default=False)
        self.parser.add_argument('-P33out', '--PHRED33output', action = 'store', nargs='?', const=True, default=True)
        self.parser.add_argument('-P64out', '--PHRED64output', action = 'store', nargs='?', const=True, default=False)
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
########################################################################

class FastAreader :
    '''
    Define objects to read FastA or FastQ files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq,header2, quality)
    '''
    def __init__ (self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname
        #Check if file is empty
        if os.stat(self.fname).st_size == 0:
            sys.stderr.write('Error: Empty File \n')
            sys.exit(1)#exit if it is

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
            header2 = ''
            quality = ''
            # skip to first fasta header
            line = fileH.readline()
            while (
                    not line.startswith('>') or
                    not line.startswith('@') or
                    not line.startswith('+')
                    ):
                if (
                    line.startswith('>') or
                    line.startswith('@') or
                    line.startswith('+')
                    ):
                    break
                line = fileH.readline()
            #if file is FastA format
            if line.startswith('>'):
                while not line.startswith('>'):
                    line = fileH.readline()
                header = line.rstrip()

                for line in fileH:
                    if line.startswith('>'):
                        yield header,sequence, header2, quality
                        header = line.rstrip()
                        sequence = ''
                    else :
                        sequence += ''.join(line.rstrip().split())
            #if file is FastQ format
            elif line.startswith('@') or line.startswith('+'):
                while line:
                    #cant look for @'s because @ is in ASCII table
                    #so just add line by line
                    header = line.rstrip()
                    line = fileH.readline()
                    sequence += ''.join(line.rstrip().split()).upper()
                    line = fileH.readline()
                    header2 = line.rstrip()
                    line = fileH.readline()
                    quality += ''.join(line.rstrip().split())
                    yield header,sequence, header2, quality
                    header,sequence,header2,quality = '','','',''
                    line = fileH.readline()
        if header.startswith('>'):
            yield header,sequence, header2, quality

class FastQTranslate:
    '''
    Class to find Q scores (Q0-Q40)
    and to convert the new Q scores into the desired output
    '''
    def __init__(self, listSeqs):
        '''initialize list of sequences to use'''
        justQs = []
        #get quality lines to self.myQualitySec
        if len(listSeqs[0]) == 2:
            for seq in listSeqs:
                justQs.append(seq[1])
        if len(listSeqs[0]) == 4:
            for seq in listSeqs:
                justQs.append(seq[3])
        self.myQualitySec = justQs

    def getQScore(self, isP33, isP64, isP64B):
        '''Get q scores for P33 P64 and P64B'''
        listQScores = []
        if isP33: #set to the offset of PHRED
            qRange = 33
        if isP64:
            qRange = 64
        if isP64B:
            qRange = 64
        for PHRED in self.myQualitySec:
            qScores = []
            #empty list to store Q scores
            for pQual in PHRED:
                #subtract to get new Q score from Q0-Q40
               qScores.append(ord(pQual)-qRange)
            listQScores.append(qScores)
        return listQScores

    def getP64SOLQScore(self):
        '''
        get P64 solexa q score
        solexa's Q score has a different equation so
        we need to find the P(error) and use that to find
        Q score of PHRED
        Qscore(SOLEXA) -> Qscore(PHRED)
        using Q = -10*log10(P(e))
        '''
        listQScores = []
        #for every quality line
        for PHRED in self.myQualitySec:
            qScores = []
            for pQual in PHRED:
                solexaQual = (ord(pQual)-64) #solexa Q score
                #find PHRED Q score
                qScores.append(round(10*math.log(10**(solexaQual/10.0)+1,10)))
            listQScores.append(qScores)
        return listQScores

    def toP33(self, listQScores):
        '''
        Changing Q score to PHRED33 mapping
        just need to add 33 to every Q score
        '''
        list33PHRED = []
        for scores in listQScores:
            P33 = []
            for Q in scores:
                P33.append(chr(Q+33))
            list33PHRED.append(P33)
        return list33PHRED

    def toP64(self, listQScores):
        '''
        Changing Q score to PHRED64 mapping
        just need to add 64 to every Q score
        '''
        list64PHRED = []
        for scores in listQScores:
            P64 = []
            for Q in scores:
                P64.append(chr(Q+64))
            list64PHRED.append(''.join(P64))
        return list64PHRED

    def fromP64B(self, listQScores, P33Out, P64Out):
        '''
        Changing Q scores to desired output from
        PHRED65 with B offset
        add the desired off set value
        '''
        PHREDList = []
        if P33Out:
            qScore = 33
        if P64Out:
            qScore = 64
        for scores in listQScores:
            newQScores = []
            for Q in scores:
                if Q<=2: #because the lowest Q=2
                    newQScores.append(chr(qScore))
                else: #otherwise just find the PHRED char
                    newQScores.append(chr(Q+qScore))
            PHREDList.append(''.join(newQScores))
        return PHREDList

    def fromP64SOL(self, listQScores, P33Out, P64Out):
        '''
        Changing Q scores to desired output from
        PHRED64 Solexa

        since we found the appropriate Q scores to convert
        to PHRED mapping we can just add the PHRED offsets to
        the Q scores
        '''
        PHREDList = []
        if P33Out:
            qScore = 33
        if P64Out:
            qScore = 64
        for scores in listQScores:
            newQScores = []
            for Q in scores:
                newQScores.append(chr(Q+qScore))
            PHREDList.append(newQScores)
        return PHREDList

def main(inCL=None):
    '''Main function to Translate Q lines'''
    #command line parsing
    if inCL is None:
        myCommandLine = CommandLine()
    else:
        myCommandLine = CommandLine(inCL)

    print(myCommandLine.args)
    # myCommandLine.args.inFile has input file name
    # myCommandLine.args.outFile has output file name
    # myCommandLine.args.PHRED33input true if entered
    # myCommandLine.args.PHRED64input true if entered
    # myCommandLine.args.PHRED64Boffin true if entered
    # myCommandLine.args.PHRED64SOLin true if entered
    # myCommandLine.args.PHRED33output true if desired
    # myCommandLine.args.PHRED64output true if desired
    myReader = FastAreader(myCommandLine.args.inFile)

    # for the printing part later
    what2WhatP = ['0','33']
    if myCommandLine.args.PHRED64output:
        myCommandLine.args.PHRED33output = False
        what2WhatP[1] = '64'
    if myCommandLine.args.PHRED33input:
        what2WhatP[0] = '33'
    if myCommandLine.args.PHRED64input:
        what2WhatP[0] = '64'
    if myCommandLine.args.PHRED64Boffin:
        what2WhatP[0] = '64B'
    if myCommandLine.args.PHRED64SOLin:
        what2WhatP[0] = '64 SOLEXA'

    #Read fastA or fastQ file and create a list of inputs
    fastQ = False
    sequences = []
    for head, seq, H2, Q in myReader.readFasta():
        addQ = []
        #if its FastQ
        if head[0] == '@' and seq and Q:
            if H2[0] == '+':
                fastQ = True
                addQ.append(head[1:])
                stripSeq = seq.replace(' ','')
                stripSeq = stripSeq.replace('*', 'N')
                addQ.append(stripSeq.replace('.','N'))
                addQ.append(H2[1:])
                addQ.append(Q)
                sequences.append(addQ)
            #check for missing 2nd header
            if H2[0] != '+':
                sys.stderr.write('\'+\' not in repeat title line.\n')
        #if its fastA
        if head[0] == '>':
            addQ.append(head[1:])
            addQ.append(seq)
            sequences.append(addQ)
        #check for missing sequence
        if seq == '':
            sys.stderr.write('\'{}\' has no sequence.\n'.format(head))
        #check for missing ID markers here
        if head[0] != '@' and head[0] !='>':
            if H2:
                sys.stderr.write('\'@\' not in ID - {}\n.'.format(head))
            if not H2:
                sys.stderr.write('\'>\' not in ID - {}\n'.format(head))

    '''
    This block of code is to see if the sequence and qline match up
    and if the IDs match up
    '''
    difSeqQual = 0
    difIDs = 0
    #keep track of nonvalid sequences
    unUsableSeqs = []
    for seqs in sequences:
        if len(seqs) == 4:
            if len(seqs[1]) != len(seqs[3]):
                difSeqQual += 1
                if seqs not in unUsableSeqs:
                    unUsableSeqs.append(seqs)
                sys.stderr.write(
                                '{} input Sequence and'.format(difSeqQual) +
                                ' Quality Lines do not match in length.\n'
                                )
                tempSeq = []
            if seqs[2] and seqs[0] != seqs[2]:
                difIDs += 1
                if seqs not in unUsableSeqs:
                    unUsableSeqs.append(seqs)
                sys.stderr.write(
                                '{} input Sequence and '.format(difIDs)+
                                'Quality IDs do not match.\n'
                                )
                tempSeq = []

    #remove nonvalid Sequences
    for qDontWork in unUsableSeqs:
        sequences.remove(qDontWork)

    #if there are no sequences exit program
    if sequences == []:
        sys.stderr.write('Error: No valid Sequences.\n')
        sys.exit(1)

    #now we can find New PHREDs
    myPHRED = FastQTranslate(sequences)#class to find Q scores
    if (
        myCommandLine.args.PHRED33input or
        myCommandLine.args.PHRED64input or
        myCommandLine.args.PHRED64Boffin
        ): #qscores for PHRED33, 64 or 64B
        qScore = myPHRED.getQScore(
                                    myCommandLine.args.PHRED33input,
                                    myCommandLine.args.PHRED64input,
                                    myCommandLine.args.PHRED64Boffin,
                                    )
    #Qscore for finding PHRED64 SOLEXA
    if myCommandLine.args.PHRED64SOLin:
        qScore = myPHRED.getP64SOLQScore()

    #add list of qscore to the sequence
    for num,scores in enumerate(qScore):
        sequences[num].append(qScore[num])

    #now get the qscores for desired output
    newPHRED = []
    if myCommandLine.args.PHRED33input: #PHRED33 -> P33 or P64
        if myCommandLine.args.PHRED33output:
            newPHRED = myPHRED.toP33(qScore)
        if myCommandLine.args.PHRED64output:
            newPHRED = myPHRED.toP64(qScore)

    if myCommandLine.args.PHRED64input: #PHRED64 -> P33 or P64
        if myCommandLine.args.PHRED33output:
            newPHRED = myPHRED.toP33(qScore)
        if myCommandLine.args.PHRED64output:
            newPHRED = myPHRED.toP64(qScore)

    #PHRED64 B off set -> P33 to P64
    if myCommandLine.args.PHRED64Boffin:
        if fastQ:
            terminalBs = []
            for pos in range(0,len(sequences)):
                temp = list(sequences[pos][1])
                #this is to replace all the terminal complementary B bases on
                #the sequence
                for Qual in range(0, len(sequences[pos][3])):
                    if sequences[pos][3][Qual:]=='B'*(len(sequences[pos][3])-Qual):
                        temp[Qual:] = 'N'*(len(sequences[pos][3])-Qual)
                        sequences[pos][1] = ''.join(temp)
                        break
        newPHRED = myPHRED.fromP64B(
                                    qScore,
                                    myCommandLine.args.PHRED33output,
                                    myCommandLine.args.PHRED64output
                                    )

    #PHRED 64 SOLEXA -> P33 or P64
    if myCommandLine.args.PHRED64SOLin:
        newPHRED = myPHRED.fromP64SOL(
                                    qScore,
                                    myCommandLine.args.PHRED33output,
                                    myCommandLine.args.PHRED64output
                                    )
    #Print to results to outfile
    file = open(myCommandLine.args.outFile, "w")
    file.write(
            'Quality Mappings of Sequences translated from ' +
            'PHRED{} to PHRED{}\n'.format(what2WhatP[0],what2WhatP[1])
            )
    for pos,translatedQmap in enumerate(sequences):
        if fastQ:
            file.write('\n')
            file.write('Sequence ID     - ')
            for char in translatedQmap[0]:
                file.write(str(char))
            file.write('\n')
            file.write('Sequence        - ')
            for char in translatedQmap[1]:
                file.write(str(char))
            file.write('\n')
            file.write('Repeat ID       - ')
            for char in translatedQmap[2]:
                file.write(str(char))
            file.write('\n')
            file.write('Original Q Map  - ')
            for char in translatedQmap[3]:
                file.write(str(char))
            file.write('\n')
            file.write('New Quality Map - ')
            for newQ in range(0,len(newPHRED[pos])):
                file.write(newPHRED[pos][newQ])
        if not fastQ:
            file.write('\n')
            file.write('\n')
            file.write('Sequence ID     - ')
            for char in translatedQmap[0]:
                file.write(str(char))
            file.write('\n')
            file.write('Original Q Map  - ')
            for char in translatedQmap[1]:
                file.write(str(char))
            file.write('\n')
            file.write('New Quality Map - ')
            for newQ in range(0,len(newPHRED[pos])):
                file.write(newPHRED[pos][newQ])
        file.write('\n')
    file.close()



if __name__ == "__main__":
    main()
