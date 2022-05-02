import sys, os
import argparse
from BLUE_object import *
from combTools import *

default_file = 'original_inputs/CMS_negCorr_CMSGrid_V5_Dec.txt'

def main():

    parser = argparse.ArgumentParser(description='specify options')

    parser.add_argument('-f',action='store',type=str, default=default_file, help='input file')
    parser.add_argument('--excludeMeas',action='store', help='provide list of measurements to be excluded. Example: --exclude \'meas 1, meas 2\'')
    parser.add_argument('--excludeSyst',action='store', help='provide list of systematics to be excluded. Example: --exclude \'syst 1, syst 2\'')
    parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat')
    parser.add_argument('--toysIndividualSyst',action='store_true', help='also run toys for each individual (relevant) systematic')
    parser.add_argument('--scanCorrAll',action='store_true', help='scan all correlations with simple assumptions')
    parser.add_argument('--excludeMeasOneByOne',action='store_true', help='exclude measurements one-by-one')
    parser.add_argument('--excludeSystOneByOne',action='store_true', help='exclude uncertainties one-by-one')
    parser.add_argument('--noSigns',action='store_true', help='remove correlation signs')
    parser.add_argument('--deriveImpactSigns',action='store_true', help='derive signs of the impacts')

    args = parser.parse_args()

    if args.f == default_file:
        print('\n WARNING: using default file "{}"'.format(default_file))

    excludeMeas = []
    if not args.excludeMeas is None:
        excludeMeas = args.excludeMeas.split(',')
        excludeMeas = [removeUselessCharachters(e) for e in excludeMeas]

    excludeSyst = []
    if not args.excludeSyst is None:
        excludeSyst = args.excludeSyst.split(',')
        excludeSyst = [removeUselessCharachters(e) for e in excludeSyst]

    base_obj = BLUE_object(args.f,excludeMeas,excludeSyst)
    if args.noSigns:
        base_obj.removeSigns()

    if args.deriveImpactSigns:
        print '\nestimating signs of impacts, this will take a short while...\n'
        base_obj.deriveSignedImpacts()

    base_obj.printResults()
    base_obj.printImpactsSorted()
    drawWeights(base_obj,path='plots')

    if args.excludeMeasOneByOne:
        excludeMeasOneByOne(base_obj)

    if args.excludeSystOneByOne:
        excludeSystOneByOne(base_obj)

    if args.scanCorrAll:
        makeCorrelationScans(base_obj)
        plotScanSummary(base_obj)

    if args.nToys > 0:
        base_obj.prepareForToys('MCstat_CMS.txt')
        base_obj.throwToys(args.nToys)
        getToyResults(base_obj)
        if args.toysIndividualSyst:
            for syst in base_obj.systForToys:
                getToyResults(base_obj,[syst])

    return

if __name__ == "__main__":
    main()
