
from BLUE_object import BLUE_object
from LHC_object import LHC_object
from combTools import makeCorrelationScansLHC

f_ATLAS = 'original_inputs/ATLASnew.txt'
f_CMS = 'original_inputs/CMS_negCorr_V4_v2_CMSGrid.txt'


def main():

    obj_ATLAS = BLUE_object(f_ATLAS,ATLAS=True)
    obj_ATLAS.printResults()

    obj_CMS = BLUE_object(f_CMS,ATLAS=False,excludeMeas=['CMS12_SVX'])
    obj_CMS.printResults()

    LHC = LHC_object(obj_ATLAS,obj_CMS)
    obj_LHC = LHC.getBlueObject()
    obj_LHC.printResults()

    makeCorrelationScansLHC(obj_LHC)

    return


if __name__ == "__main__":
    main()


