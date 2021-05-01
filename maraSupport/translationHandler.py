
translationTable = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'
}

def __makeReverseTranslationTable__():
    reverseTanslationTable = {}
    for codon, aminoAcid in translationTable.items():
        if not aminoAcid in reverseTanslationTable:
            reverseTanslationTable[aminoAcid] = []
        reverseTanslationTable[aminoAcid].append(codon)
    return reverseTanslationTable

reverseTranslationTable = __makeReverseTranslationTable__()


def characterizeMismatches(string1:str, string2:str):
    if not len(string1) == len(string2):
        raise ValueError("Can only characterize mismatches on strings of identical length")
    matchingMatrix = []
    mismatchCounter = 0
    for i in range(len(string1)):
        if string1[i] == string2[i]:
            matchingMatrix.append(True)
        else:
            mismatchCounter += 1
            matchingMatrix.append(False)
    return mismatchCounter, matchingMatrix


def findPotentialSubstitutions(referenceCodon:str, alternateAminoAcid:str):
    potentialAlternativeCodons = reverseTranslationTable[alternateAminoAcid.upper()]
    alternativeCharacterizations = {1:[], 2:[], 3:[]}
    for potentialAlternativeCodon in potentialAlternativeCodons:
        if potentialAlternativeCodon == referenceCodon:
            continue
        mismatches, matchingMatrix = characterizeMismatches(referenceCodon, potentialAlternativeCodon)
        alternativeCharacterizations[mismatches].append((potentialAlternativeCodon, matchingMatrix))
    return alternativeCharacterizations


def findLikeliestSubstitution(referenceCodon:str, alternateAminoAcid:str):
    alternativeCharacterizations = findPotentialSubstitutions(referenceCodon,alternateAminoAcid)
    mostLikelyMismatchCount = 0
    for potentialMismatchCount in [1, 2, 3]:
        if alternativeCharacterizations[potentialMismatchCount]:
            mostLikelyMismatchCount = potentialMismatchCount
            break
    if mostLikelyMismatchCount == 0:
        raise ValueError("Tried to find an alternative codon for %s, but no alternative appears to exist.")
    potentialAlternativeCodons = []
    matchingMatrices = []
    for potentialAlternativeCodon, matchingMatrix in alternativeCharacterizations[mostLikelyMismatchCount]:
        potentialAlternativeCodons.append(potentialAlternativeCodon)
        matchingMatrices.append(matchingMatrix)
    potentialMismatchSiteMatrix = []
    for i in range(len(matchingMatrices[0])):
        potentialForMismatch = False
        for matchingMatrix in matchingMatrices:
            if not matchingMatrix[i]:
                potentialForMismatch = True
                break
        potentialMismatchSiteMatrix.append(potentialForMismatch)
    return potentialAlternativeCodons, potentialMismatchSiteMatrix, alternativeCharacterizations


if __name__ == "__main__":
    mostLikelyAltCodon, mostLikelyMismatchMatrix, allAlternativeCodons = findLikeliestSubstitution("GAT", "G")
    print("something")