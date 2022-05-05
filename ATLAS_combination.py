from BLUE_object import BLUE_object
from combTools import *
import argparse

infile = 'original_inputs/ATLASallinputs_2022_05_05.txt'

def makeATLAS_MCstat_file(obj):
    f = open('MCstat_ATLAS.txt','w')
    for syst in obj.usedSyst:
        f.write('\t{}'.format(syst))
    f.write('\n')
    for meas in obj.usedMeas:
        fmeas = open('../top-mass-combination/inputs/ATLASforLHC/ATLAS_{}.txt'.format(meas),'r')
        lines = fmeas.read().splitlines()
        f.write(meas)
        for syst in obj.usedSyst:
            for line in lines:
                if line.split()[0] == syst:
                    break
            f.write(' {}'.format(line.split()[-1]))
        f.write('\n')
    return

def main():

    parser = argparse.ArgumentParser(description='specify options')

    parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat')
    parser.add_argument('--toysIndividualSyst',action='store_true', help='also run toys for each individual (relevant) systematic')

    args = parser.parse_args()

    base_obj = BLUE_object(infile,ATLAS=True)
    base_obj.printResults()

    if args.nToys > 0:
        makeATLAS_MCstat_file(base_obj)
        base_obj.prepareForToys('MCstat_ATLAS.txt')
        base_obj.throwToys(args.nToys)
        getToyResults(base_obj,plotToys=False)
        if args.toysIndividualSyst:
            for syst in base_obj.systForToys:
                getToyResults(base_obj,[syst])

    return


if __name__ == "__main__":
    main()


