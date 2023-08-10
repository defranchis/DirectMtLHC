
from BLUE_object import BLUE_object
from LHC_object import LHC_object
from LHC_tools import makeAllCorrelationScansLHC, flipAmbiguousSigns
from combTools import getToyResults, getToyResultsLHCobj, excludeMeasOneByOne, makeCorrelationScans, plotScanSummary, drawWeights, measToTex, measToROOT
import argparse, sys, copy, os
import numpy as np

import default_files

f_ATLAS = default_files.default_file_ATLAS
f_CMS = default_files.default_file_CMS
PU_hack = True

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
    parser.add_argument('--scanbJES',action='store_true', help='scan b JES correlations with both methods')
    parser.add_argument('--flipSigns',action='store_true', help='flip all ambiguous signs in LHC correlations')
    parser.add_argument('--blind', dest='unblind', action='store_false', help='blind the LHC combination')
    parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat', default=0)
    parser.add_argument('--subCombinations',action='store_true', help='perform sub-combinations')
    parser.add_argument('--onlyWeightsAbove',action='store',type=float, help='re-perform combination with ony weights above given value')
    parser.add_argument('--excludeMeasOneByOne',action='store_true', help='exclude measurements one-by-one')

    args = parser.parse_args()

    obj_ATLAS = BLUE_object(f_ATLAS, ATLAS=True, PU_hack=PU_hack)
    # obj_ATLAS.printResults()
    # obj_ATLAS.printImpactsSorted()

    obj_CMS = BLUE_object(f_CMS,ATLAS=False,PU_hack=PU_hack)
    # obj_CMS.printResults()`
    # obj_CMS.printImpactsSorted()

    LHC_full_unblind = LHC_object(obj_ATLAS, obj_CMS, blind=False, separateCombinations=False, PU_hack=PU_hack)
    LHC_sep_unblind = LHC_object(obj_ATLAS, obj_CMS, blind=False, separateCombinations=True, PU_hack=PU_hack)

    LHC_full_unblind.printCorrTables(draw=True)
    LHC_full_unblind.BLUE_obj.printPullWeightsTable(blind=not args.unblind)
    if args.unblind:
        LHC_full_unblind.printPullWeightsComparisonTable()
    LHC_full_unblind.printSummaryTableLHC()
    LHC_full_unblind.printAllImpactsSorted()
    # LHC_full_unblind.makeSummaryPlot(blind=not args.unblind)
    
    LHC_sep_unblind.BLUE_obj.printPullWeightsTable(blind=not args.unblind)

    print('CMS combination \n')
    # LHC_full_unblind.obj_d['CMS'].simplePrint()
    # LHC_full_unblind.obj_d['CMS'].printImpactsSorted()


    if not args.unblind:
        LHC_full = LHC_object(obj_ATLAS, obj_CMS, blind=True, separateCombinations=False, PU_hack=PU_hack)
        LHC_full.printResults()
        LHC_full.printImpactsSorted()

        LHC_sep = LHC_object(obj_ATLAS, obj_CMS, blind=True, separateCombinations=True, PU_hack=PU_hack)
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


    drawWeights(LHC_full_unblind.BLUE_obj)
    produceSummaryTable(LHC_full_unblind,blind=not args.unblind)


    if args.scanAllCorr or args.scanbJES:
        makeAllCorrelationScansLHC(LHC_full_unblind,LHC_sep_unblind,blind=not args.unblind, only_bJES = args.scanbJES)
        plotScanSummary(LHC_full_unblind.BLUE_obj,blind=not args.unblind,syst_list=list(LHC_full_unblind.corrMap.keys()))

    if args.flipSigns:
        flipAmbiguousSigns(LHC_full_unblind,LHC_sep_unblind)

    if args.nToys > 0:

        # full combination
        LHC_full_unblind.CMS_obj.prepareForToys('inputs/MCstat_CMS_forLHC.txt')
        LHC_full_unblind.ATLAS_obj.prepareForToys('MCstat_ATLAS.txt')
        makeLHC_MCstat_file(LHC_full_unblind.ATLAS_obj,LHC_full_unblind.CMS_obj)        

        LHC_obj_full = LHC_full_unblind.getBlueObject()
        LHC_obj_full.prepareForToys('MCstat_LHC.txt')
        LHC_obj_full.throwToys(args.nToys)
        getToyResults(LHC_obj_full,plotToys=False,blind=not args.unblind)

        # separate combiantions
        LHC_sep_unblind.CMS_obj.prepareForToys('inputs/MCstat_CMS_forLHC.txt')
        LHC_sep_unblind.ATLAS_obj.prepareForToys('MCstat_ATLAS.txt')
        LHC_sep_unblind.CMS_obj.throwToys(args.nToys)
        LHC_sep_unblind.ATLAS_obj.throwToys(args.nToys)
        getToyResultsLHCobj(LHC_sep_unblind,blind=not args.unblind)
        
    if args.subCombinations:
        CMS = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if 'CMS' in meas]
        ATLAS = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if not meas in CMS]
        obsDict = {'ATLAS':ATLAS, 'CMS':CMS}
        LHC_full_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True,jsonForPlot=True)

        CMS = [meas for meas in LHC_sep_unblind.BLUE_obj.usedMeas if 'CMS' in meas]
        ATLAS = [meas for meas in LHC_sep_unblind.BLUE_obj.usedMeas if not meas in CMS]
        obsDict = {'ATLAS':ATLAS, 'CMS':CMS}
        LHC_sep_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)

        ll = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if measToTex(meas)=='$dil$']
        lj = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if measToTex(meas)=='$lj$']
        aj = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if measToTex(meas)=='$aj$']
        other = [meas for meas in LHC_full_unblind.BLUE_obj.usedMeas if not meas in ll and not meas in lj and not meas in aj]
        obsDict = {'ll':ll, 'lj':lj, 'aj':aj, 'other':other}
        if args.unblind:
            LHC_full_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)
        else:
            res, unc = LHC_full_unblind.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=False)
            LHC_full.BLUE_obj.doSubCombination(obsDict=obsDict,printResults=True)
            print('\ndifferences:')
            for ch in obsDict.keys():
                print('{} - full = {:.2f} GeV'.format(ch,res[ch]-LHC_full_unblind.BLUE_obj.results.mt))
            


    if not args.onlyWeightsAbove is None:

        r = LHC_full_unblind.BLUE_obj.doCombinationWeightsAbove(wmin=args.onlyWeightsAbove,printout = args.unblind)

        print ('only measurements with weight above {}%'.format(args.onlyWeightsAbove*100))
        print ('n. used meas = {}'.format(len(r.weights.keys())))
        print ('mt - mt_orig = {:.2f} GeV'.format(r.mt-LHC_sep_unblind.BLUE_obj.results.mt))
        print ('tot uncert = {:.2f} GeV\n'.format(r.tot))

    if args.excludeMeasOneByOne:
        excludeMeasOneByOne(LHC_full_unblind.BLUE_obj,blind=not args.unblind)


    print(LHC_full_unblind.BLUE_obj.chi2)
    print(LHC_full_unblind.BLUE_obj.prob)
    print()
    print(LHC_sep_unblind.BLUE_obj.chi2)
    print(LHC_sep_unblind.BLUE_obj.prob)


    return

def getTotUncFromDict(ud,no_stat=False):
    unc = np.array([ud[syst] for syst in list(ud.keys())])
    if not no_stat: return (np.sum(unc**2))**.5
    else: return (np.sum(unc**2) - ud['Stat']**2)**.5

def produceSummaryTable(LHC_obj,blind=True):

    obj = LHC_obj.BLUE_obj
    ll = [meas for meas in obj.usedMeas if measToTex(meas)=='$dil$']
    lj = [meas for meas in obj.usedMeas if measToTex(meas)=='$lj$']
    aj = [meas for meas in obj.usedMeas if measToTex(meas)=='$aj$']
    other = [meas for meas in obj.usedMeas if not meas in ll and not meas in lj and not meas in aj]
    obsDict = {'dil':ll, 'lj':lj, 'aj':aj, 'other':other}
    res, unc = obj.doSubCombination(obsDict=obsDict,printResults=False)

    of = open('summary_table_LHC.txt','w')
    of.write('input name, \t mt, \t stat, \t syst, \t tot\n\n')
    
    all_dict = dict()

    for meas in obj.usedMeas:
        mt, stat, syst, tot = obj.value[meas], obj.uncert[meas]['Stat'], getTotUncFromDict(obj.uncert[meas],no_stat=True), getTotUncFromDict(obj.uncert[meas])
        of.write('{},\t {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format(measToROOT(meas).replace('#',''),mt,stat,syst,tot))
        all_dict[measToROOT(meas).replace('#','')] = [mt,stat,syst,tot]

    of.write('\n')
    of.write('{},\t {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format('ATLAS_comb',LHC_obj.ATLAS_obj.results.mt,LHC_obj.ATLAS_obj.results.stat,LHC_obj.ATLAS_obj.results.syst,LHC_obj.ATLAS_obj.results.tot))
    of.write('{},\t {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format('CMS_comb',LHC_obj.CMS_obj.results.mt,LHC_obj.CMS_obj.results.stat,LHC_obj.CMS_obj.results.syst,LHC_obj.CMS_obj.results.tot))
    all_dict['ATLAS_comb'] = [LHC_obj.ATLAS_obj.results.mt,LHC_obj.ATLAS_obj.results.stat,LHC_obj.ATLAS_obj.results.syst,LHC_obj.ATLAS_obj.results.tot]
    all_dict['CMS_comb'] = [LHC_obj.CMS_obj.results.mt,LHC_obj.CMS_obj.results.stat,LHC_obj.CMS_obj.results.syst,LHC_obj.CMS_obj.results.tot]
    of.write('\n')

    of.write('\n')
    for ch in obsDict.keys():
        mt = res[ch] if not blind else 172 + (res[ch]-obj.results.mt)
        stat = unc[ch]['Stat']
        tot = getTotUncFromDict(unc[ch])
        of.write('{},\t {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format(ch+' comb',mt,stat,(tot**2-stat**2)**.5,tot))
        all_dict[ch+' comb'] = [mt,stat,(tot**2-stat**2)**.5,tot]
    of.write('\n')
    mt_comb = obj.results.mt if not blind else 172
    of.write('{},\t {:.2f}, {:.2f}, {:.2f}, {:.2f}\n'.format('full comb',mt_comb,obj.results.stat,obj.results.syst,obj.results.tot))
    all_dict['full comb'] = [mt_comb,obj.results.stat,obj.results.syst,obj.results.tot]
    import json
    with open('summary_LHC.json','w') as j:
        j.write(json.dumps(all_dict))

    os.system('python3 plotter/summaryPlot.py')
        
    return


if __name__ == "__main__":
    main()


