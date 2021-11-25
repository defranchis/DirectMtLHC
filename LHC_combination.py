
from BLUE_object import BLUE_object
from LHC_object import LHC_object
from combTools import makeCorrelationScansLHC
import argparse

f_ATLAS = 'original_inputs/ATLAS_signed_2021_10_05.txt'
# f_ATLAS = 'original_inputs/ATLASnew.txt'
f_CMS = 'original_inputs/CMS_negCorr_CMSGrid_V5_Oct18.txt'
# f_CMS = 'original_inputs/CMS_negCorr_V4_v2_CMSGrid.txt'

def main():

    # parser = argparse.ArgumentParser(description='specify options')
    # parser.add_argument('--scanCorr',action='store_true', help='scan all correlations with simple assumptions')
    # parser.add_argument('--scanBrutal',action='store_true', help='scan all correlations wildly')
    # args = parser.parse_args()

    obj_ATLAS = BLUE_object(f_ATLAS, ATLAS=True)
    # # obj_ATLAS = BLUE_object(f_ATLAS,ATLAS=True)
    obj_ATLAS.printResults()
    # obj_ATLAS.printImpactsSorted()
    # print obj_ATLAS.matrix['JESFLV']

    obj_CMS = BLUE_object(f_CMS,ATLAS=False)
    obj_CMS.printResults()
    # obj_CMS.printImpactsSorted()

    LHC_full = LHC_object(obj_ATLAS,obj_CMS,blind=True, separateCombinations = False)
    # LHC = LHC_object(obj_ATLAS,obj_CMS,separateCombinations=False)
    LHC_full.getBlueObject().printResults()
    LHC_full.getBlueObject().printImpactsSorted()
    # obj_LHC.printResults()
    # obj_LHC.printImpactsSorted()

    LHC_sep = LHC_object(obj_ATLAS,obj_CMS,blind=True, separateCombinations = True)
    # LHC = LHC_object(obj_ATLAS,obj_CMS,separateCombinations=False)
    LHC_sep.getBlueObject().printResults()
    LHC_sep.getBlueObject().printImpactsSorted()


    # if args.scanCorr or args.scanBrutal:
    #     makeCorrelationScansLHC(obj_LHC,args.scanBrutal)

    return


if __name__ == "__main__":
    main()


