
from BLUE_object import BLUE_object
from LHC_object import LHC_object
from LHC_tools import makeAllCorrelationScansLHC, flipAmbiguousSigns
from combTools import getToyResults, getToyResultsLHCobj, excludeMeasOneByOne, makeCorrelationScans, plotScanSummary
import argparse, sys, copy

import default_files

f_ATLAS = default_files.default_file_ATLAS
f_CMS = default_files.default_file_CMS

def makeLHC_MCstat_file(obj_ATLAS,obj_CMS):
    systForToys = copy.deepcopy(obj_ATLAS.systForToys)
    for syst in obj_CMS.systForToys:
        if not syst in systForToys:
            systForToys.append(syst)
    f = open('MCstat_LHC.txt','w')
    for syst in systForToys:
        f.write('\t{}'.format(syst))
    f.write('\n')
    for meas in list(obj_ATLAS.MCstat_d.keys()):
        f.write(meas)
        for syst in systForToys:
            if syst in list(obj_ATLAS.MCstat_d[meas].keys()):
                f.write(' {}'.format(obj_ATLAS.MCstat_d[meas][syst]))
            else:
                f.write(' 0')
        f.write('\n')
    for meas in list(obj_CMS.MCstat_d.keys()):
        f.write(meas)
        for syst in systForToys:
            if syst in list(obj_CMS.MCstat_d[meas].keys()):
                f.write(' {}'.format(obj_CMS.MCstat_d[meas][syst]))
            else:
                f.write(' 0')
        f.write('\n')

    return

def main():

    parser = argparse.ArgumentParser(description='specify options')

    parser.add_argument('--scanAllCorr',action='store_true', help='scan all correlations with simple assumptions for both methods')
    parser.add_argument('--flipSigns',action='store_true', help='flip all ambiguous signs in LHC correlations')
    parser.add_argument('--unblind',action='store_true', help='do not blind the LHC combination')
    parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat', default=0)
    parser.add_argument('--subCombinations',action='store_true', help='perform sub-combinations')
    parser.add_argument('--onlyWeightsAbove',action='store',type=float, help='re-perform combination with ony weights above given value')
    parser.add_argument('--excludeMeasOneByOne',action='store_true', help='exclude measurements one-by-one')

    args = parser.parse_args()

    obj_ATLAS = BLUE_object(f_ATLAS, ATLAS=True)
    # obj_ATLAS.printResults()
    # obj_ATLAS.printImpactsSorted()

    obj_CMS = BLUE_object(f_CMS,ATLAS=False)
    # obj_CMS.printResults()`
    # obj_CMS.printImpactsSorted()

    LHC_full_unblind = LHC_object(obj_ATLAS, obj_CMS, blind=False, separateCombinations=False)
    LHC_sep_unblind = LHC_object(obj_ATLAS, obj_CMS, blind=False, separateCombinations=True)

    LHC_full_unblind.printCorrTables(draw=True)
    LHC_full_unblind.makeSummaryPlot(blind=not args.unblind)
    
    print('CMS combination \n')
    LHC_full_unblind.obj_d['CMS'].simplePrint()
    LHC_full_unblind.obj_d['CMS'].printImpactsSorted()


    if not args.unblind:
        LHC_full = LHC_object(obj_ATLAS, obj_CMS, blind=True, separateCombinations=False)
        LHC_full.printResults()
        LHC_full.printImpactsSorted()

        LHC_sep = LHC_object(obj_ATLAS, obj_CMS, blind=True, separateCombinations=True)
        LHC_sep.printResults()
        LHC_sep.printImpactsSorted()
    
    else:
        LHC_full_unblind.printResults()
        LHC_full_unblind.printImpactsSorted()
        LHC_sep_unblind.printResults()
        LHC_sep_unblind.printImpactsSorted()

    print()
    print('|separate / full -1 |')
    print(abs(LHC_sep_unblind.getBlueObject().results.mt/LHC_full_unblind.getBlueObject().results.mt -1)*100, '%\n')
    print()
    print('|full - separate|')
    print(abs(LHC_sep_unblind.getBlueObject().results.mt - LHC_full_unblind.getBlueObject().results.mt), 'GeV\n')


    LHC_full_unblind.BLUE_obj.printWeights(prefix='LHC_A')
    LHC_sep_unblind.BLUE_obj.printWeights(prefix='LHC_B')

    if args.unblind:
        LHC_full_unblind.BLUE_obj.printPulls(prefix='LHC_A')
        LHC_sep_unblind.BLUE_obj.printPulls(prefix='LHC_B')


    if args.scanAllCorr:
        makeAllCorrelationScansLHC(LHC_full_unblind,LHC_sep_unblind,blind=not args.unblind)
        makeCorrelationScans(LHC_full_unblind.BLUE_obj,blind=not args.unblind)
        plotScanSummary(LHC_full_unblind.BLUE_obj,blind=not args.unblind)

    if args.flipSigns:
        flipAmbiguousSigns(LHC_full_unblind,LHC_sep_unblind)

    if args.nToys > 0:

        # full combination
        LHC_full_unblind.CMS_obj.prepareForToys('MCstat_CMS_forLHC.txt')
        LHC_full_unblind.ATLAS_obj.prepareForToys('MCstat_ATLAS.txt')
        makeLHC_MCstat_file(LHC_full_unblind.ATLAS_obj,LHC_full_unblind.CMS_obj)        

        LHC_obj_full = LHC_full_unblind.getBlueObject()
        LHC_obj_full.prepareForToys('MCstat_LHC.txt')
        LHC_obj_full.throwToys(args.nToys)
        getToyResults(LHC_obj_full,plotToys=False,blind=not args.unblind)

        # separate combiantions
        LHC_sep_unblind.CMS_obj.prepareForToys('MCstat_CMS_forLHC.txt')
        LHC_sep_unblind.ATLAS_obj.prepareForToys('MCstat_ATLAS.txt')
        LHC_sep_unblind.CMS_obj.throwToys(args.nToys)
        LHC_sep_unblind.ATLAS_obj.throwToys(args.nToys)
        getToyResultsLHCobj(LHC_sep_unblind,blind=not args.unblind)
        
    if args.subCombinations:
        CMS = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if 'CMS' in meas]
        ATLAS = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if not meas in CMS]
        obsDict = {'ATLAS':ATLAS, 'CMS':CMS}
        LHC_full_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)

        CMS = [meas for meas in LHC_sep_unblind.BLUE_obj.usedMeas if 'CMS' in meas]
        ATLAS = [meas for meas in LHC_sep_unblind.BLUE_obj.usedMeas if not meas in CMS]
        obsDict = {'ATLAS':ATLAS, 'CMS':CMS}
        LHC_sep_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)

        from combTools import measToTex
        ll = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if measToTex(meas)=='$ll$']
        lj = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if measToTex(meas)=='$lj$']
        aj = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if measToTex(meas)=='$aj$']
        other = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if not meas in ll and not meas in lj and not meas in aj]
        obsDict = {'ll':ll, 'lj':lj, 'aj':aj, 'other':other}
        if args.unblind:
            LHC_full_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)
        else:
            LHC_full.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)

    if not args.onlyWeightsAbove is None:

        r = LHC_full_unblind.BLUE_obj.doCombinationWeightsAbove(wmin=args.onlyWeightsAbove,printout = args.unblind)

        print ('only measurements with weight above {}%'.format(args.onlyWeightsAbove*100))
        print ('n. used meas = {}'.format(len(r.weights.keys())))
        print ('mt - mt_orig = {:.2f} GeV'.format(r.mt-LHC_sep_unblind.BLUE_obj.results.mt))
        print ('tot uncert = {:.2f} GeV\n'.format(r.tot))

    if args.excludeMeasOneByOne:
        excludeMeasOneByOne(LHC_full_unblind.BLUE_obj,blind=not args.unblind)


    return


if __name__ == "__main__":
    main()


