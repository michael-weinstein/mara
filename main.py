
import maraSupport
import sys
import typing


def getMutationFromArguments():
    args = sys.argv
    if not len(args) == 3:
        raise ValueError("This program takes only 1 argument: the mutation being analyzed")
    mutation = args[2]


def analyzeMutation(mutation:str):
    if not maraSupport.altSequence.validateMutation(mutation):
        raise ValueError("Unable to validate the mutation as readable")
    mostLikelyMutations, otherPossibleMutations = maraSupport.altSequence.getAltCodons(mutation)


def printMutationResults(codonList:typing.List[maraSupport.altSequence.PotentialAlternativeSequence]):
    for possibleMutation in codonList:
        print(possibleMutation.json)


def printResults(mostLikelyCodons:list, otherPossibleCodons:list):
    print("MOST LIKELY CHANGES")
    printMutationResults(mostLikelyCodons)
    print()
    print("OTHER POSSIBLE CHANGES")
    printMutationResults(otherPossibleCodons)


def main():
    mutation = getMutationFromArguments()
    mostLikelyMutations, otherPossibleMutations = analyzeMutation(mutation)
    printResults(mostLikelyMutations, otherPossibleMutations)
