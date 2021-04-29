import os
from Bio import SeqIO, Seq
import io
import csv

referenceGenomePath = os.path.join(os.path.split(__file__)[0], "reference", "NC_045512.2.gb")
genomeLocusTablePath = os.path.join(os.path.split(__file__)[0], "reference", "covid genome map table.csv")
geneAliasTablePath = os.path.join(os.path.split(__file__)[0], "reference", "aliasList.csv")



class ReferenceGenome:

    def __init__(self, gbkPath:str=referenceGenomePath):
        self.gbkPath = gbkPath
        if not os.path.isfile(gbkPath):
            raise FileNotFoundError("Unable to find reference file %s" %gbkPath)
        file = open(gbkPath, 'r')
        self.__gbkStream = io.StringIO()
        self.__gbkStream.write(file.read())
        file.close()
        self.__gbkStream.seek(0)
        self.__gbHandle = SeqIO.parse(self.__gbkStream, 'genbank')
        recordList = [record for record in self.__gbHandle]
        if not recordList:
            raise RuntimeError("Reference genome %s returned no records" % self.gbkPath)
        self.gbRecord = recordList[0]
        if len(recordList) > 1:
            print("WARNING: %s appears to have multiple records" % gbkPath)
        self.makeAttributes()

    def makeAttributes(self):
        self.geneTable = self.getGeneTable()
        self.geneTableByLocation = self.getGeneTableByLocation()
        self.sequenceTable = self.getSequenceTableByGene()
        self.translations = self.getTranslationTableByGene()

    def getGeneTable(self):
        def getGeneList():
            genes = {}
            for feature in self.gbRecord.features:
                if feature.type == "gene":
                    genes[feature.qualifiers["gene"][0]] = feature
                elif feature.type == "mat_peptide":
                    genes[feature.qualifiers["product"][0]] = feature
                elif "UTR" in feature.type:
                    genes[feature.type] = feature
            return genes
        geneList = getGeneList()
        return geneList

    def getGeneTableByLocation(self):
        geneTableByLocation = {}
        for geneID, geneData in self.geneTable.items():
            geneStart = geneData.location.nofuzzy_start
            geneEnd = geneData.location.nofuzzy_end
            geneLocationTuple = (geneStart, geneEnd)
            geneTableByLocation[geneLocationTuple] = geneID
        return geneTableByLocation

    def getSequenceTableByGene(self):
        sequenceTable = {}
        for coordinates, geneID in self.geneTableByLocation.items():
            start, end = coordinates
            sequence = self.gbRecord.seq[start:end]
            sequenceTable[geneID] = str(sequence)
        return sequenceTable

    def getTranslationTableByGene(self):
        translationTable = {}
        for geneID, sequence in self.sequenceTable.items():
            translationTable[geneID] = Seq.translate(sequence)
        return translationTable


def getAliasTable(aliasTable:str=geneAliasTablePath):
    if not os.path.isfile(aliasTable):
        raise FileNotFoundError("Unable to find alias table at %s" %aliasTable)
    fileHandle = open(aliasTable, 'r')
    csvHandle = csv.reader(fileHandle)
    aliasTable = {}
    for line in csvHandle:
        if line[0].startswith("#"):
            continue
        filteredLine = [item.strip().upper() for item in line if item]
        if not filteredLine:
            continue
        identifier = filteredLine[0]
        aliases = filteredLine[1:]
        aliasTable[identifier] = identifier
        for alias in aliases:
            aliasTable[alias] = identifier
    return aliasTable


def getAdditionalLociFromList(locusTable:str=genomeLocusTablePath):
    if not os.path.isfile(locusTable):
        raise FileNotFoundError("Unable to find locus table at %s" %locusTable)
    fileHandle = open(locusTable, 'r')
    csvHandle = csv.reader(fileHandle)
    locusTable = {}
    for line in csvHandle:
        if line[0].startswith("#"):
            continue
        filteredLine = [item.strip().upper() for item in line if item]
        if not filteredLine:
            continue
        gene, start, end = filteredLine[:3]
        gene = gene.upper()
        start = int(start) - 1
        end = int(end)
        locusTable[gene] = (start, end)
    return locusTable


def makeGenomeCoordinateTable(referenceGenome:ReferenceGenome, aliasTable:dict=None, extendedLocusTable:dict=None):
    if aliasTable is None:
        aliasTable = getAliasTable()
    if extendedLocusTable is None:
        extendedLocusTable = getAdditionalLociFromList()
    genomeCoordinateTable = {}
    for coordinates, gene in referenceGenome.geneTableByLocation.items():
        identifier = aliasTable[gene.upper()]
        genomeCoordinateTable[identifier] = coordinates
    for gene, coordinates in extendedLocusTable.items():
        identifier = aliasTable[gene]
        if identifier in aliasTable:
            continue
        else:
            genomeCoordinateTable[identifier] = coordinates
    return genomeCoordinateTable

genomeCoordinateTable = makeGenomeCoordinateTable(ReferenceGenome())

if __name__ == "__main__":
    reference = ReferenceGenome()
    genomeCoordinateTable = makeGenomeCoordinateTable(reference)
    print("something")


