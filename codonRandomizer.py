'''
Summary

This script takes an amino acid sequence and returns a randomized DNA sequence that codes for that amino acid sequence.
The amino acid sequence can be given either directly as an argument or loaded from a text file. 
The likelyhood of different codons being used can be changed by modifying the codonWeight variable.

---------------------------

Example
    The following terminal command will create 20 randomized DNA sequences of the amino acid sequence QQSSVA.
    python3 codonRandomizer.py 20 QQSSVA

    If the amino acid sequence was instead stored in a text file called aminoAcidSequence.txt, the following command is equivalent.
    python3 codonrandomizer.py 20 aminoAcidSequence.txt
    
---------------------------

Arguments

[1] : Interger
    Number of DNA sequences to be generated

[2] : String or text file location
    The amino acid sequence to be turned into a DNA sequence.
    If the amino acid sequence is in a text file

---------------------------

Output

By default the DNA sequence will be saved to output.txt, which can be changed by modifying the outfile variable.
If multiple DNA sequences are being generated they will be separated by dashed lines.


'''
import os
import sys
import random

outfile = 'output.txt'

codonDict = {
        'A': ['GCT', 'GCC', 'GCA', 'GCG'],                  #Alanine (Ala)
        'R': ['CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],    #Arginine (Arg)
        'N': ['AAT', 'AAC'],                                #Asparagine (Asn)
        'D': ['GAT', 'GAC'],                                #Aspartic acid (Asp)
        'C': ['TGT', 'TGC'],                                #Cysteine (Cys)
        'Q': ['CAA', 'CAG'],                                #Glutamine (Gln)
        'E': ['GAA', 'GAG'],                                #Glutamic acid (Glu)
        'G': ['GGT', 'GGC', 'GGA', 'GGG'],                  #Glycine (Gly)
        'H': ['CAT', 'CAC'],                                #Histidine (His)
        'I': ['ATT', 'ATC', 'ATA'],                         #Isoleucine (Ile)
        'L': ['CTT', 'CTC', 'CTA', 'CTG', 'TTA', 'TTG'],    #Leucine (Leu)
        'K': ['AAA', 'AAG'],                                #Lysine (Lys)
        'M': ['ATG'],                                       #Methionine (Met)
        'F': ['TTT', 'TTC'],                                #Phenylalanine (Phe)
        'P': ['CCT', 'CCC', 'CCA', 'CCG'],                  #Proline (Pro)
        'S': ['TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'],    #Serine (Ser)
        'T': ['ACT', 'ACC', 'ACA', 'ACG'],                  #Threonine (Thr)
        'W': ['TGG'],                                       #Tryptophan (Trp)
        'Y': ['TAT', 'TAC'],                                #Tyrosine (Tyr)
        'V': ['GTT', 'GTC', 'GTA', 'GTC'],                  #Valine (Val)
        's': ['TAA', 'TGA', 'TAG']                          #Stop
        }
codonWeights = {
        #Empty list means giving even probabilities across all possible codons
        #Weights are normalized, so [1, 2, 0] is the same as [3, 6, 0]
        #Note that weights must be integers
        #example: this would give equal probabilities for GCT and GCC codons when producing Alanine
        #'A': [1, 1, 0, 0]
        'A': [1, 1, 1, 0],
        'R': [],
        'N': [],
        'D': [],
        'C': [],
        'Q': [],
        'E': [],
        'G': [1, 1, 1, 0],
        'H': [],
        'I': [1, 1, 0],
        'L': [0, 0, 0, 0, 1, 1],
        'K': [],
        'M': [],
        'F': [],
        'P': [1, 1, 1, 0],
        'S': [1, 1, 1, 0, 0, 0],
        'T': [],
        'W': [],
        'Y': [],
        'V': [1, 0, 1, 0],
        's': []
        }


#Check if number of codons in codonDict and CodonWeights are equal
for aa in codonDict:
    if not(len(codonWeights[aa]) == 0 or len(codonWeights[aa]) == len(codonDict[aa])):
        print('Error: number of weights and number of codons do not match for ', aa)
        exit()

n = sys.argv[1]
try: n = int(n)
except Exception as error:
    print('Error: First argument must be a number (int or float)')
    print(error)
    exit()

#verify n
if n <= 0:
    print('Error: number of DNA sequences to be generated must be greater than 0, you gave ' + str(n))
    exit()
#Build amino acid sequence
if os.path.isfile(sys.argv[2]):
    aaSequence = ''
    with open(sys.argv[2], 'r') as f:
        for line in f:
            if line[-1:] in ['\n', ' ', '\t']: aaSequence += line[:-1]
            else: aaSequence += line
else: aaSequence = sys.argv[2]

#Replace [] in codon weights with [1,1..]
for acid in codonWeights:
    if codonWeights[acid] == []:
        codonWeights[acid] = [1]*len(codonDict[acid])

#Verify that input amino acid sequence is viable
for acid in aaSequence:
    if acid == '\n':
        continue
    if acid not in codonDict:
        print('Error: found an amino acid that is not in the codon dictionary: ' + str(acid))
        exit()

#generate codon sequence and save it
out = open(outfile, 'w+')
for i in range(n):
    codonSequence = ''
    for acid in aaSequence:
        choiceList = []
        for weightIndex, weight in enumerate(codonWeights[acid]):
            choiceList += [codonDict[acid][weightIndex]]*weight
        codonSequence += random.choice(choiceList)
    out.write(codonSequence + '\n')
    out.write('-------------------------- \n')
out.close()
print('Successfully saved codon sequences at ' + str(outfile))
