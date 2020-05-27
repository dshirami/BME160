#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

'''
Read a DNA string from user input and return a collapsed substring of embedded Ns to: {count}.

Example:
 input: AaNNNNNNGTC
output: AA{6}GTC

Any lower case letters are converted to uppercase
'''

class DNAstring (str):
    def length (self):
        return (length(self))

    def purify(self):
        ''' Return an upcased version of the string, collapsing a single run of Ns.'''
        UpperDNA = self.upper() # make the DNA string to be all uppercase
        """Check to see if there are Ns in DNA string. But since we are assuming there can only be 1 block of Ns at MOST,
           An if statement will do.
        """
        if "N" in UpperDNA:
            nCount = UpperDNA.count("N") # count the # of "N"s in the DNA string
            indexN = UpperDNA.index("N") # look for where the N block starts
            DNAList = list(UpperDNA) #convert string to list so we can edit the DNA string
            del DNAList[indexN:nCount+indexN] # delete the Ns
            DNAList.insert(indexN, "{{{0}}}" .format(nCount)) #add {"number of Ns"} where the N block was
            UpperDNA = ''.join(DNAList) #convert list back to string so we can print it accordingly
        return UpperDNA


def main():
    ''' Get user DNA data and clean it up.'''
    data = input('DNA data?')
    thisDNA = DNAstring (data)
    pureData = thisDNA.purify()
    print (pureData)

main()
