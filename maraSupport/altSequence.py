
import json

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

    def __str__(self):
        return "%s/%s" %(self.ref, self.alt)

    @property
    def dictionary(self):
        returnDict = {
            "ref": self.ref,
            "alt": self.alt
        }
        return returnDict


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

    @property
    def dictionary(self):
        returnDict = {
            "codon start": self.start,
            "codon end": self.end,
            "amino acids": self.aminoAcids.dictionary,
            "codons": self.codons.dictionary,
            "potential mismatch sites": self.potentialMismatchSites,
            "potential mismatch matrix": self.potentialMismatchMatrix
        }
        return returnDict

    @property
    def json(self):
        return json.dumps(self.dictionary, indent=2)

    def __str__(self):
        mismatchString = ""
        for site in self.potentialMismatchMatrix:
            if not site:
                mismatchString += "_"
            else:
                mismatchString += "*"
        return "%s-%s-%s %s" %(self.start + 1, self.codons.alt, self.end, mismatchString)



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
    mostLikelyAltCodons, mostLikelyMismatchMatrix, allAlternativeCodons = translationHandler.findLikeliestSubstitution(refCodon, alt)
    mostLikelyAltsReturn = []
    for mostLikelyAltCodon in mostLikelyAltCodons:
        mostLikelyAltsReturn.append(PotentialAlternativeSequence(codonStart, codonEnd, RefAltPair(ref, alt), RefAltPair(refCodon, mostLikelyAltCodon)))
    otherPotentialAltReturn = {}
    for mismatchCount in allAlternativeCodons:
        alternativesToAdd = []
        if not allAlternativeCodons[mismatchCount]:
            continue
        for altCodon, matchMatrix in allAlternativeCodons[mismatchCount]:
            if altCodon in mostLikelyAltCodons:
                continue
            else:
                alternativesToAdd.append(PotentialAlternativeSequence(codonStart, codonEnd, RefAltPair(ref, alt), RefAltPair(refCodon, altCodon)))
        if alternativesToAdd:
            otherPotentialAltReturn[mismatchCount] = alternativesToAdd.copy()
    return mostLikelyAltsReturn, otherPotentialAltReturn


def parseSubstitutionVariant(mutation:str):
    gene, variation = mutation.split(":")
    ref = variation[0]
    alt = variation[-1]
    codonNumber = int(variation[1:-1])
    return gene, ref, alt, codonNumber


def getAltCodons(mutation:str):
    gene, ref, alt, codonNumber = parseSubstitutionVariant(mutation)
    mostLikelyAlts, otherPotentialAlts = getPotentialAlternatives(gene, ref, codonNumber, alt)
    return mostLikelyAlts, otherPotentialAlts


def validateMutation(mutation:str):
    if not ":" in mutation:
        raise ValueError("Mutation must be formatted as 'gene:change")
    gene, change = mutation.split(":")
    ref = change[0]
    codon = change[1:-1]
    alt = change[-1]
    if not ref.isalpha():
        raise ValueError("Reference amino acid must be a letter")
    try:
        codon = int(codon)
    except:
        raise ValueError("Codon number must be an integer")
    if not alt.isalpha():
        raise ValueError("Alternative amino acid must be a letter")
    return True


if __name__ == "__main__":
    mostLikelys, otherPotentials = getPotentialAlternatives("S", "D", 614, "G")
    print(mostLikelys[0])
    print("something")