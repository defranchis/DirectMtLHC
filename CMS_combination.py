import sys, os
import copy
import numpy as np
import ROOT
from ROOT import TH1F, TCanvas, TGraph
import argparse
from BLUE_object import *

ROOT.gROOT.SetBatch(True)
np.random.seed(1)

def main():

    parser = argparse.ArgumentParser(description='specify options')

    parser.add_argument('-f',action='store',type=str, required=True, help='<Required> input file')
    parser.add_argument('--excludeMeas',action='store', help='provide list of measurements to be excluded. Example: --exclude \'meas 1, meas 2\'')
    parser.add_argument('--excludeSyst',action='store', help='provide list of systematics to be excluded. Example: --exclude \'syst 1, syst 2\'')
    parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat')
    parser.add_argument('--toysIndividualSyst',action='store_true', help='also run toys for each individual (relevant) systematic')
    parser.add_argument('--scanCorrAll',action='store_true', help='scan all correlations with simple assumptions')
    parser.add_argument('--excludeMeasOneByOne',action='store_true', help='exclude measurements one by one')

    args = parser.parse_args()

    if (not args.nToys is None) and args.noBLUE:
        print '\nERROR: cannot use nToys option with noBLUE option'
        print 'exiting...\n'
        sys.exit()

    excludeMeas = []
    if not args.excludeMeas is None:
        excludeMeas = args.excludeMeas.split(',')
        excludeMeas = [removeUselessCharachters(e) for e in excludeMeas]

    excludeSyst = []
    if not args.excludeSyst is None:
        excludeSyst = args.excludeSyst.split(',')

    obj = BLUE_object(args.f,excludeMeas,excludeSyst)
    obj.printResults()


if __name__ == "__main__":
    main()
