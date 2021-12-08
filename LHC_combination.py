
from BLUE_object import BLUE_object
from LHC_object import LHC_object
from LHC_tools import makeAllCorrelationScansLHC, flipAmbiguousSigns
import argparse, sys, copy

f_ATLAS = 'original_inputs/ATLAS_signed_2021_10_05.txt'
f_CMS = 'original_inputs/CMS_negCorr_CMSGrid_V5_Oct18.txt'

def main():

    parser = argparse.ArgumentParser(description='specify options')
    parser.add_argument('--scanAllCorr',action='store_true', help='scan all correlations with simple assumptions for both methods')
    parser.add_argument('--flipSigns',action='store_true', help='flip all ambiguous signs in LHC correlations')
    parser.add_argument('--unblind',action='store_true', help='do not blind the LHC combination')
    args = parser.parse_args()

    obj_ATLAS = BLUE_object(f_ATLAS, ATLAS=True)
    obj_ATLAS.printResults()
    # obj_ATLAS.printImpactsSorted()

    obj_CMS = BLUE_object(f_CMS,ATLAS=False)
    obj_CMS.printResults()
    # obj_CMS.printImpactsSorted()

    excludeSyst = ['PU'] # temporarily hard-coded

    if not args.unblind:
        LHC_full = LHC_object(obj_ATLAS, obj_CMS, blind=True, separateCombinations=False, excludeSyst=excludeSyst)
        LHC_full.getBlueObject().printResults()

        LHC_sep = LHC_object(obj_ATLAS, obj_CMS, blind=True, separateCombinations=True, excludeSyst=excludeSyst)
        LHC_sep.getBlueObject().printResults()
    
    LHC_full_unblind = LHC_object(obj_ATLAS, obj_CMS, blind=False, separateCombinations=False, excludeSyst=excludeSyst)
    LHC_sep_unblind = LHC_object(obj_ATLAS, obj_CMS, blind=False, separateCombinations=True, excludeSyst=excludeSyst)

    if args.unblind:
        LHC_full_unblind.getBlueObject().printResults()
        LHC_sep_unblind.getBlueObject().printResults()

    print
    print 'separate / full'
    print (LHC_sep_unblind.getBlueObject().results.mt/LHC_full_unblind.getBlueObject().results.mt -1)*100, '%\n'

    if args.scanAllCorr:
        makeAllCorrelationScansLHC(LHC_full_unblind,LHC_sep_unblind,blind=not args.unblind)
    if args.flipSigns:
        flipAmbiguousSigns(LHC_full_unblind,LHC_sep_unblind)

    return


if __name__ == "__main__":
    main()


