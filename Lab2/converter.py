#!/usr/bin/env python3
# Name: Dylan Shiramizu (dshirami)
# Group Members: None

'''
Convert between sequence information and amino acid representation
Examples:
input: ATG
output: ATG = MET
input: UAG
output: UAG = ---
input: E
output: E = GLU
input: Asp
output: ASP = D

also works with lowercase inputs
'''
short_AA = {
            'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'
            }

long_AA = {value:key for key,value in short_AA.items()}

rnaCodonTable = {
# Second Base
# U             C             A             G
#U
'UUU': 'Phe', 'UCU': 'Ser', 'UAU': 'Tyr', 'UGU': 'Cys',
'UUC': 'Phe', 'UCC': 'Ser', 'UAC': 'Tyr', 'UGC': 'Cys',
'UUA': 'Leu', 'UCA': 'Ser', 'UAA': '---', 'UGA': '---',
'UUG': 'Leu', 'UCG': 'Ser', 'UAG': '---', 'UGG': 'Trp',
#C
'CUU': 'Leu', 'CCU': 'Pro', 'CAU': 'His', 'CGU': 'Arg',
'CUC': 'Leu', 'CCC': 'Pro', 'CAC': 'His', 'CGC': 'Arg',
'CUA': 'Leu', 'CCA': 'Pro', 'CAA': 'Gln', 'CGA': 'Arg',
'CUG': 'Leu', 'CCG': 'Pro', 'CAG': 'Gln', 'CGG': 'Arg',
#A
'AUU': 'Ile', 'ACU': 'Thr', 'AAU': 'Asn', 'AGU': 'Ser',
'AUC': 'Ile', 'ACC': 'Thr', 'AAC': 'Asn', 'AGC': 'Ser',
'AUA': 'Ile', 'ACA': 'Thr', 'AAA': 'Lys', 'AGA': 'Arg',
'AUG': 'Met', 'ACG': 'Thr', 'AAG': 'Lys', 'AGG': 'Arg',
#G
'GUU': 'Val', 'GCU': 'Ala', 'GAU': 'Asp', 'GGU': 'Gly',
'GUC': 'Val', 'GCC': 'Ala', 'GAC': 'Asp', 'GGC': 'Gly',
'GUA': 'Val', 'GCA': 'Ala', 'GAA': 'Glu', 'GGA': 'Gly',
'GUG': 'Val', 'GCG': 'Ala', 'GAG': 'Glu', 'GGG': 'Gly'
}
dnaCodonTable = {key.replace('U','T'):value for key, value in rnaCodonTable.items()}

def main():
    '''
    Take the sequence or representation uppercase all of them and run them in:
    shortlook
    shortlook but with keys and values swapped
    longlook
    longlook but with keys and values swapped
    '''
    lookUp = input("Enter codon: ")
    lookUp = lookUp.upper() #make the input all upper
    finalCodon = 'unknown'
    shortLook = short_AA.get(lookUp) #look thorugh shortlook for the key
    if shortLook != None: #if the value exists make finalCodon the value
        finalCodon = shortLook
    rnaLook = rnaCodonTable.get(lookUp) #look through backwards shortlook for the key
    if rnaLook != None:#if the value exists make finalCodon the value
        finalCodon = rnaLook
    longLook = long_AA.get(lookUp) #look through longlook for the key
    if longLook != None:#if the value exists make finalCodon the value
        finalCodon = longLook
    dnaLook = dnaCodonTable.get(lookUp) #look through backwards longlook for the key
    if dnaLook != None:#if the value exists make finalCodon the value
        finalCodon = dnaLook
    print("{0} = {1}" .format(lookUp, finalCodon))
    pass

main()
