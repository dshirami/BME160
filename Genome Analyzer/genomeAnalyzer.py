#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

from sequenceAnalysis import NucParams, FastAreader
def main ():
    myReader = FastAreader() # make sure to change this to use stdin
    myNuc = NucParams()
    for head, seq in myReader.readFasta() :
        myNuc.addSequence(seq)

    #Finding sequence length
    aaComp = myNuc.aaComposition()
    nucBases = 0
    for base in aaComp:
        nucBases += (aaComp[base]*3)
    print('Sequence Length = {:0.2f} Mb'.format(nucBases/1000000))

    print()

    #Finding GC content
    nucComp = myNuc.nucComposition()
    totalGC = nucComp['G'] + nucComp['C']
    totalAAs = sum(nucComp.values())
    print('GC content = {:0.1f}%'.format((totalGC/totalAAs)*100))

    print()

    # sort codons in alpha order, by Amino Acid
    nucs = {key:val for key, val in
            sorted(myNuc.rnaCodonTable.items(),
                key=lambda item:(item[1],item[0]))}

    # calculate relative codon usage for each codon and print
    val = 0
    thisCodonComp = myNuc.codonComposition()
    for nucI in nucs:
        aa = nucs[nucI]
        total = aaComp[aa]
        used = thisCodonComp[nucI]
        if total>0:
            val =  used/total
        print ('{:s} : {:s} {:5.1f} ({:6d})'.format(nucI,
            nucs[nucI], val*100, thisCodonComp[nucI]))

if __name__ == "__main__":
    main()
