import sys, os
import argparse
from BLUE_object import *
from combTools import *

import default_files

from LHC_object import mergeMap_default, renameMap_default

default_file = default_files.default_file_CMS

def main():

    parser = argparse.ArgumentParser(description='specify options')

    parser.add_argument('-f',action='store',type=str, default=default_file, help='input file')
    parser.add_argument('--excludeMeas',action='store', help='provide list of measurements to be excluded. Example: --exclude \'meas 1, meas 2\'')
    parser.add_argument('--excludeSyst',action='store', help='provide list of systematics to be excluded. Example: --exclude \'syst 1, syst 2\'')
    parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat', default=0)
    parser.add_argument('--toysIndividualSyst',action='store_true', help='also run toys for each individual (relevant) systematic')
    parser.add_argument('--scanCorrAll',action='store_true', help='scan all correlations with simple assumptions')
    parser.add_argument('--excludeMeasOneByOne',action='store_true', help='exclude measurements one-by-one')
    parser.add_argument('--excludeSystOneByOne',action='store_true', help='exclude uncertainties one-by-one')
    parser.add_argument('--noSigns',action='store_true', help='remove correlation signs')
    parser.add_argument('--deriveImpactSigns',action='store_true', help='derive signs of the impacts')
    parser.add_argument('--subCombinations',action='store_true', help='perform sub-combinations')
    parser.add_argument('--onlyWeightsAbove',action='store',type=float, help='re-perform combination with ony weights above given value')

    args = parser.parse_args()

    if args.f == default_file:
        print(('\n WARNING: using default file "{}"'.format(default_file)))

    excludeMeas = []
    if not args.excludeMeas is None:
        excludeMeas = args.excludeMeas.split(',')
        excludeMeas = [removeUselessCharachters(e) for e in excludeMeas]

    excludeSyst = []
    if not args.excludeSyst is None:
        excludeSyst = args.excludeSyst.split(',')
        excludeSyst = [removeUselessCharachters(e) for e in excludeSyst]

    base_obj = BLUE_object(args.f,excludeMeas=excludeMeas,excludeSyst=excludeSyst)
    for old, new in renameMap_default['CMS'].items():
        base_obj.renameSyst(old,new)

    base_obj.printFullCorrTable()
    base_obj.printPullWeightsTable()

    if args.noSigns:
        base_obj.removeSigns()

    if args.deriveImpactSigns:
        print('\nestimating signs of impacts, this will take a short while...\n')
        base_obj.deriveSignedImpacts()

        
    base_obj.printSummaryTable(CMS_grid=True)
    
    clone = base_obj.clone()
    for merged, original_l in list(mergeMap_default['CMS'].items()):
        clone.mergeSyst(merged,original_l)
    clone.printSummaryTable()


    base_obj.printResults()
    clone.printImpactsSorted()
    base_obj.printStats()
    drawWeights(base_obj,path='plots')
    base_obj.printSummaryCorrTableCMS()

    if args.excludeMeasOneByOne:
        excludeMeasOneByOne(base_obj)

    if args.excludeSystOneByOne:
        excludeSystOneByOne(base_obj)

    if args.scanCorrAll:
        makeCorrelationScans(base_obj)
        plotScanSummary(base_obj)

    if args.nToys > 0:
        base_obj.prepareForToys('inputs/MCstat_CMS.txt')
        base_obj.throwToys(args.nToys)
        getToyResults(base_obj)
        if args.toysIndividualSyst:
            for syst in base_obj.systForToys:
                getToyResults(base_obj,[syst])

    if args.subCombinations:
        l_7TeV = [m for m in base_obj.usedMeas if 'CMS11' in m]
        l_8TeV = [m for m in base_obj.usedMeas if m not in l_7TeV]
        obsDict = {'7TeV':l_7TeV, '8TeV':l_8TeV}
        base_obj.doSubCombination(obsDict=obsDict,printResults=True)

    if not args.onlyWeightsAbove is None:
        base_obj.doCombinationWeightsAbove(wmin=args.onlyWeightsAbove,printout=True)


if __name__ == "__main__":
    main()
