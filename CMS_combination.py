import sys, os
import copy
import numpy as np
import ROOT
from ROOT import TH1F, TCanvas, TGraph
import argparse
from BLUE_object import *

ROOT.gROOT.SetBatch(True)

toys_dir = 'toys_workdir'
scan_dir = 'scan_workdir'

def getToyResults(base_obj,l=[]):
    l_mt, l_tot, l_stat, l_syst, d_weights, d_syst = base_obj.getToyResults(l)
    print np.array(l_mt).mean(), np.array(l_mt).std()
    #then do the plotting, etc
    return

def excludeMeasOneByOne(base_obj):
    for meas in base_obj.usedMeas:
        obj = base_obj.clone()
        obj.addExcludeMeas([meas])
        print 'excluded:', obj.excludeMeas
        obj.simplePrint()
        print
    return


def decreaseAllStrongCorrelations(base_obj):
    step = 0.01
    red_corrs = list(np.arange(0.7,1+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()

    for i,red_corr in enumerate(red_corrs):
        obj = base_obj.clone()
        obj.reduceCorrelations(red_corr,'ALL')
        h_mt.SetPoint(i,red_corr,obj.results.mt)
        h_tot.SetPoint(i,red_corr,obj.results.tot)
        h_stat.SetPoint(i,red_corr,obj.results.stat)
        h_syst.SetPoint(i,red_corr,obj.results.syst)
        
    h_mt.SetTitle('correlation scan (reduced); correlation; mt [GeV]')
    h_tot.SetTitle('correlation scan (reduced); correlation; tot uncert [GeV]')
    h_stat.SetTitle('correlation scan (reduced); correlation; stat uncert [GeV]')
    h_syst.SetTitle('correlation scan (reduced); correlation; syst uncert [GeV]')

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    c = TCanvas()
    h_mt.Draw('apl')
    c.SaveAs('{}/mt_red.png'.format(scan_dir))
    c.Clear()

    h_tot.Draw('apl')
    c.SaveAs('{}/tot_red.png'.format(scan_dir))
    c.Clear()

    h_stat.Draw('apl')
    c.SaveAs('{}/stat_red.png'.format(scan_dir))
    c.Clear()

    h_syst.Draw('apl')
    c.SaveAs('{}/syst_red.png'.format(scan_dir))
    c.Clear()

    return

def increaseAllWeakCorrelations(base_obj):
    step = 0.01
    incr_corrs = list(np.arange(-.3,0.3+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()
    for i,incr_corr in enumerate(incr_corrs):
        obj = base_obj.clone()
        obj.increaseCorrelations(incr_corr,'ALL')
        h_mt.SetPoint(i,incr_corr,obj.results.mt)
        h_tot.SetPoint(i,incr_corr,obj.results.tot)
        h_stat.SetPoint(i,incr_corr,obj.results.stat)
        h_syst.SetPoint(i,incr_corr,obj.results.syst)
        
    h_mt.SetTitle('correlation scan (increased); correlation; mt [GeV]')
    h_tot.SetTitle('correlation scan (increased); correlation; tot uncert [GeV]')
    h_stat.SetTitle('correlation scan (increased); correlation; stat uncert [GeV]')
    h_syst.SetTitle('correlation scan (increased); correlation; syst uncert [GeV]')

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    c = TCanvas()
    h_mt.Draw('apl')
    c.SaveAs('{}/mt_incr.png'.format(scan_dir))
    c.Clear()

    h_tot.Draw('apl')
    c.SaveAs('{}/tot_incr.png'.format(scan_dir))
    c.Clear()

    h_stat.Draw('apl')
    c.SaveAs('{}/stat_incr.png'.format(scan_dir))
    c.Clear()

    h_syst.Draw('apl')
    c.SaveAs('{}/syst_incr.png'.format(scan_dir))
    c.Clear()

def makeCorrelationScan(base_obj,syst_scan):
    if not os.path.exists(scan_dir+'/syst'):
        os.makedirs(scan_dir+'/syst')
    step = 0.01
    red_corrs = list(np.arange(0.7,1+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()
    print '\n scanning systematics {} \n'.format(syst_scan)
    for i,red_corr in enumerate(red_corrs):
        obj = base_obj.clone()
        obj.reduceCorrelations(red_corr,syst_scan)
        h_mt.SetPoint(i,red_corr,obj.results.mt)
        h_tot.SetPoint(i,red_corr,obj.results.tot)
        h_stat.SetPoint(i,red_corr,obj.results.stat)
        h_syst.SetPoint(i,red_corr,obj.results.syst)
        
    h_mt.SetTitle('mt scan - {}; correlation; mt [GeV]'.format(syst_scan))
    h_tot.SetTitle('tot uncert scan - {}; correlation; tot uncert [GeV]'.format(syst_scan))
    h_stat.SetTitle('stat uncert scan - {}; correlation; stat uncert [GeV]'.format(syst_scan))
    h_syst.SetTitle('syst uncert scan - {}; correlation; syst uncert [GeV]'.format(syst_scan))

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    c = TCanvas()
    h_mt.Draw('apl')
    c.SaveAs('{}/syst/mt_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/mt_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    h_tot.Draw('apl')
    c.SaveAs('{}/syst/tot_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/tot_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    h_stat.Draw('apl')
    c.SaveAs('{}/syst/stat_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/stat_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    h_syst.Draw('apl')
    c.SaveAs('{}/syst/syst_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/syst_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    return

def makeCorrelationScans(base_obj):

    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)

    decreaseAllStrongCorrelations(base_obj)
    increaseAllWeakCorrelations(base_obj)

    f = open('{}/slides.tex'.format(scan_dir),'w')

    for syst in base_obj.usedSyst:
        if syst == 'Stat': continue
        makeCorrelationScan(base_obj,syst)
        f.write('\\begin{{frame}}{{correlation scan for systematics {}}}\n'.format(syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/mt_{}.pdf}}\n'.format(scan_dir,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/tot_{}.pdf}}\\\\\n'.format(scan_dir,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/stat_{}.pdf}}\n'.format(scan_dir,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/syst_{}.pdf}}\n'.format(scan_dir,syst))
        f.write('\\end{frame}\n\n')

    return

    

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

    excludeMeas = []
    if not args.excludeMeas is None:
        excludeMeas = args.excludeMeas.split(',')
        excludeMeas = [removeUselessCharachters(e) for e in excludeMeas]

    excludeSyst = []
    if not args.excludeSyst is None:
        excludeSyst = args.excludeSyst.split(',')

    base_obj = BLUE_object(args.f,excludeMeas,excludeSyst)
    base_obj.printResults()

    if args.excludeMeasOneByOne:
        excludeMeasOneByOne(base_obj)

    if args.scanCorrAll:
        makeCorrelationScans(base_obj)

    if args.nToys > 0:
        base_obj.prepareForToys('MCstat.txt')
        base_obj.throwToys(args.nToys)
        getToyResults(base_obj)
        if args.toysIndividualSyst:
            for syst in base_obj.systForToys:
                getToyResults(base_obj,[syst])

if __name__ == "__main__":
    main()
