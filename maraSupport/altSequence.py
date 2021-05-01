
try:
    import codonGetter
    import translationHandler
except ImportError:
    from . import codonGetter
    from . import translationHandler

class RefAltPair:

    def __init__(self, ref, alt):
        self.ref = ref
        self.alt = alt


class PotentialAlternativeSequence:

    def __init__(self, start:int, end:int, aminoAcids:RefAltPair, codons:RefAltPair):
        self.start = start
        self.end = end
        self.aminoAcids = aminoAcids
        self.codons = codons

    @property
    def potentialMismatchMatrix(self):
        if not len(self.codons.ref) == len(self.codons.alt):
            raise ValueError("Got unmatching codon lengths for comparison: %s and %s" %(self.codons.ref, self.codons.alt))
        mismatchMatrix = []
        for baseSet in zip(self.codons.ref, self.codons.alt):
            if baseSet[0] == baseSet[1]:
                mismatchMatrix.append(False)
            else:
                mismatchMatrix.append(True)
        return mismatchMatrix

    @property
    def potentialMismatchSites(self):
        potentialMismatchSites = []
        currentPosition = self.start  #working with a base-zero coordinate here
        for positionMismatched in self.potentialMismatchMatrix:
            if positionMismatched:
                potentialMismatchSites.append(currentPosition)
            currentPosition += 1
        return potentialMismatchSites



def getPotentialAlternatives(gene:str, ref:str, codon:int, alt:str):
    codonStart = codonGetter.getGenomicPositionFromCodon(gene, codon)
    codonEnd = codonStart + 3
    refCodon = codonGetter.getCodonRefSeq(gene, codon)
    refAminoAcidFromSequence = translationHandler.translationTable[refCodon]
    if not refAminoAcidFromSequence == ref.upper():
        print("WARNING: Mutation %s:%s%s%s appears to claim %s as a reference amino acid, but the genome shows %s" %(
            gene,
            ref,
            codon,
            alt,
            ref,
            refAminoAcidFromSequence
            )
        )
    mostLikelyAltCodon, mostLikelyMismatchMatrix, allAlternativeCodons = translationHandler.findLikeliestSubstitution(refCodon, alt)
    mostLikelyAltReturn = PotentialAlternativeSequence(codonStart, codonEnd, RefAltPair(ref, alt), RefAltPair(refCodon, mostLikelyAltCodon), mostLikelyMismatchMatrix)
    otherPotentialAltReturn = {}
    for mismatchCount in allAlternativeCodons:
        alternativesToAdd = []
        if not allAlternativeCodons[mismatchCount]:
            continue
        for altCodon, matchMatrix in allAlternativeCodons[mismatchCount]:
            if altCodon in mostLikelyAltCodon:
                continue
            else:
                alternativesToAdd.append((altCodon, matchMatrix))
        if alternativesToAdd:
            otherPotentialAltReturn[mismatchCount] = alternativesToAdd.copy()
    return mostLikelyAltCodon, mostLikelyMismatchMatrix, otherPotentialAltReturn


if __name__ == "__main__":
    likelyCodon, likelyMismatch, otherPotentials = getPotentialAlternatives("S", "D", 614, "G")
    print("something")