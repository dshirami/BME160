#!/usr/bin/env python3
# Name: Dylan Shiramizu (dsirami)
# Group Members: None

'''
Read a seqname line of a FASTQ and parse through each field and displays each on a new line

Example:
input: @EAS139:136:FC706VJ:2:2104:15343:197393
output:
Instrument = EAS139
Run ID = 136
Flow Cell ID = FC706VJ
Flow Cell Lane = 2
Tile Number = 2104
X-coord = 15343
Y-coord = 197393

'''

class FastqString (str):
    '''
    Classes to parse through the input.
    '''
    def parse(self):
        '''
        get rid of the @ and split the string at the ":"s and return a list
        Assuming the input is always right there is no need to check if it contains 7 different fields
        '''
        listFASTQ = self[1:]
        listFASTQ = listFASTQ.split(":")
        return listFASTQ
        pass

def main():
    ''' Turn the string input into a list using .parse() and print each field.'''
    FASTQ = input("FASTQ Seqname: ")
    fastSeq = FastqString(FASTQ)
    fastSeq = fastSeq.parse()
    print("Instrument = {0}" .format(fastSeq[0]))
    print("Run ID = {0}" .format(fastSeq[1]))
    print("Flow Cell ID = {0}" .format(fastSeq[2]))
    print("Flow Cell Lane = {0}" .format(fastSeq[3]))
    print("Tile Number = {0}" .format(fastSeq[4]))
    print("X-Coord = {0}" .format(fastSeq[5]))
    print("Y-Coord = {0}" .format(fastSeq[6]))
    pass

main()
