import os
from Bio import SeqIO
import io
try:
    import referenceGenomeHandler
except ImportError:
    from . import referenceGenomeHandler

referenceGenomePath = os.path.join(os.path.split(__file__)[0], "reference", "NC_045512.2.gb")
aliasTable = referenceGenomeHandler.getAliasTable()
genomeCoordinateTable = referenceGenomeHandler.genomeCoordinateTable


class RefAltPair:

    def __init__(self, ref, alt):
        self.ref = ref
        self.alt = alt


class PotentialAlternativeSequence:

    def __init__(self, start:int, end:int, aminoAcids:RefAltPair, codons:RefAltPair, potentialMismatchMatrix:list):
        self.start = start
        self.end = end
        self.aminoAcids = aminoAcids
        self.codons = codons
        self.potentialMismatchMatrix = potentialMismatchMatrix



def getGenomeSequence(gbkPath:str=referenceGenomePath):
    if not os.path.isfile(gbkPath):
        raise FileNotFoundError("Unable to find reference file %s" %gbkPath)
    file = open(gbkPath, 'r')
    gbkStream = io.StringIO()
    gbkStream.write(file.read())
    file.close()
    gbkStream.seek(0)
    gbHandle = SeqIO.parse(gbkStream, 'genbank')
    recordList = [record for record in gbHandle]
    if not recordList:
        raise RuntimeError("Reference genome %s returned no records" % gbkPath)
    gbRecord = recordList[0]
    return str(gbRecord.seq)


def getGenomicPositionFromCodon(gene:str, codon:int):
    try:
        geneID = aliasTable[gene.upper()]
    except KeyError:
        raise KeyError("Unable to find a gene matching alias %s" %gene)
    geneStart, geneStop = genomeCoordinateTable[geneID]
    positionWithinGene = (codon - 1) * 3
    codonStart = geneStart + positionWithinGene
    return codonStart


def getCodonRefSeq(gene:str, codon:int, genomeSequence:str):
    codonStart = getGenomicPositionFromCodon(gene, codon)
    codonEnd = codonStart + 3
    return genomeSequence[codonStart:codonEnd]






if __name__ == "__main__":
    genomicSequence = getGenomeSequence()
    refCodon = getCodonRefSeq("S", 614, genomicSequence)
    print("something")