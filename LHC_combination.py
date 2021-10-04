
from BLUE_object import BLUE_object
from LHC_object import LHC_object
from combTools import makeCorrelationScansLHC
import argparse

f_ATLAS = 'original_inputs/ATLASnew.txt'
f_CMS = 'original_inputs/CMS_negCorr_V4_v2_CMSGrid.txt'


def main():

    parser = argparse.ArgumentParser(description='specify options')
    parser.add_argument('--scanCorr',action='store_true', help='scan all correlations with simple assumptions')
    parser.add_argument('--scanBrutal',action='store_true', help='scan all correlations wildly')
    args = parser.parse_args()

    obj_ATLAS = BLUE_object(f_ATLAS,ATLAS=True)
    obj_ATLAS.printResults()
    obj_ATLAS.printImpactsSorted()

    obj_CMS = BLUE_object(f_CMS,ATLAS=False,excludeMeas=['CMS12_SVX'])
    obj_CMS.printResults()
    obj_CMS.printImpactsSorted()

    LHC = LHC_object(obj_ATLAS,obj_CMS)
    # LHC = LHC_object(obj_ATLAS,obj_CMS,separateCombinations=False)
    obj_LHC = LHC.getBlueObject()
    obj_LHC.printResults()
    obj_LHC.printImpactsSorted()

    if args.scanCorr or args.scanBrutal:
        makeCorrelationScansLHC(obj_LHC,args.scanBrutal)

    return


if __name__ == "__main__":
    main()


