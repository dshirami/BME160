#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

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

# Please do not modify any of the following.  This will produce a standard output that can be parsed

import sys
def main():
    inString = input('protein sequence?')
    while inString :
        myParamMaker = ProteinParam(inString)
        myAAnumber = myParamMaker.aaCount()
        print ("Number of Amino Acids: {aaNum}".format(aaNum = myAAnumber))
        print ("Molecular Weight: {:.1f}".format(myParamMaker.molecularWeight()))
        print ("molar Extinction coefficient: {:.2f}".format(myParamMaker.molarExtinction()))
        print ("mass Extinction coefficient: {:.2f}".format(myParamMaker.massExtinction()))
        print ("Theoretical pI: {:.2f}".format(myParamMaker.pI()))
        print ("Amino acid composition:")
        myAAcomposition = myParamMaker.aaComposition()
        keys = list(myAAcomposition.keys())
        keys.sort()
        if myAAnumber == 0 : myAAnumber = 1  # handles the case where no AA are present
        for key in keys :
            print ("\t{} = {:.2%}".format(key, myAAcomposition[key]/myAAnumber))

        inString = input('protein sequence?')

if __name__ == "__main__":
    main()
