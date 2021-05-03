
import maraSupport
import sys
import typing
import json


def getMutationFromArguments():
    args = sys.argv
    if not len(args) == 2:
        raise ValueError("This program takes only 1 argument: the mutation being analyzed")
    mutation = args[1]
    return mutation


def analyzeMutation(mutation:str):
    if not maraSupport.altSequence.validateMutation(mutation):
        raise ValueError("Unable to validate the mutation as readable")
    mostLikelyMutations, otherPossibleMutations = maraSupport.altSequence.getAltCodons(mutation)
    return mostLikelyMutations, otherPossibleMutations


def getResultDict(mostLikelyCodons:typing.List[maraSupport.altSequence.PotentialAlternativeSequence], otherPossibleCodons:typing.Dict[int, typing.List[maraSupport.altSequence.PotentialAlternativeSequence]]):
    resultDict = {
        "best": [],
        "other": {}
    }
    for potentialCodon in mostLikelyCodons:
        resultDict["best"].append(potentialCodon.dictionary)
    for changeCount, potentialCodonList in otherPossibleCodons.items():
        resultDict["other"][changeCount] = [potentialCodon.dictionary for potentialCodon in potentialCodonList]
    return resultDict


def printResultDict(resultDict):
    print(json.dumps(resultDict, indent=2))


def main():
    mutation = getMutationFromArguments()
    mostLikelyMutations, otherPossibleMutations = analyzeMutation(mutation)
    resultDict = getResultDict(mostLikelyMutations, otherPossibleMutations)
    printResultDict(resultDict)


if __name__ == "__main__":
    main()