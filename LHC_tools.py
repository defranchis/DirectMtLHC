import os, sys, copy
import numpy as np
import ROOT
from ROOT import TCanvas, TGraph, TLegend
from BLUE_object import BLUE_object
from LHC_object import LHC_object

ROOT.gROOT.SetBatch(True)
np.random.seed(1)
LHC_dir = 'LHC_scans'

def makeAllCorrelationScansLHC(full,sep,blind=False):

    if not os.path.exists(LHC_dir):
        os.makedirs(LHC_dir)

    if full.blind or sep.blind:
        if not blind:
            print 'ERROR: unblinded scan requested but objects are blinded'
            sys.exit()
        else:
            print '*******'
            print 'WARNING: one or more input objects are blind - mt scan will be meaningless'
            print '*******'

    mt_smear = 0
    rand_offset_small = 0
    if blind:
        mt_smear = np.random.normal(0,1) # 1 GeV smear
        rand_offset_small = np.random.normal(0,0.1)

    LHC_full = full.clone()
    LHC_sep = sep.clone()
    orig_corrMap = copy.deepcopy(LHC_full.corrMap)
    if LHC_sep.corrMap != orig_corrMap:
        print 'ERROR: correlation maps should be the same for both methods'
        sys.exit()
    for syst in LHC_full.LHCsyst:
        if syst == 'Stat': continue
        makeCorrelationScan(LHC_full,LHC_sep,orig_corrMap,syst,mt_smear,rand_offset_small)
    return

def makeCorrelationScan(LHC_full,LHC_sep,corrMap,syst,mt_smear=0,rand_offset_small=0):
    print 'scanning', syst

    obj_d = {'full' : LHC_full, 'separate' : LHC_sep}

    orig_corrMap = copy.deepcopy(corrMap)
    if not syst in orig_corrMap:
        orig_corrMap[syst] = 0.

    corr = orig_corrMap[syst]
    step = 0.05
    halfrange = min(0.25, 1-abs(corr))
    corrs = list(np.arange(corr-halfrange,corr+halfrange+step/2,step))
    corrs = [round(c,3) for c in corrs]

    scan_d = dict()
    for meth, obj in obj_d.items():
        result_l = scanOne(obj,orig_corrMap,syst,corrs)
        scan_d[meth] = result_l

    plotScanResults(corrs,scan_d,syst,'tot')
    plotScanResults(corrs,scan_d,syst,'mt',mt_smear,rand_offset_small)

    return

def scanOne(obj,corrMap,syst,corrs):
    orig_corrMap = copy.deepcopy(corrMap)
    result_l = []
    for corr in corrs:
        cm = copy.deepcopy(orig_corrMap)
        cm[syst] = corr
        tmp_obj = obj.clone()
        tmp_obj.setNewLHCcorrMap(cm)
        result_l.append(tmp_obj.getBlueObject().results)
    return result_l

def plotScanResults(corrs,scan_d,syst,variable,mt_smear=0,rand_offset_small=0):
    c = TCanvas()
    leg = TLegend(.15,.15,.4,.3)
    leg.SetBorderSize(0)
    gd = dict()
    for meth, scan in scan_d.items():
        g = TGraph()
        g.SetName(meth)
        for i, corr in enumerate(corrs):
            if variable == 'tot':
                g.SetPoint(i,corr,scan[i].tot)
            elif variable == 'mt':
                g.SetPoint(i,corr,scan[i].mt+mt_smear)
            else:
                print 'ERROR: variable {} not supported'.format(variable)
                sys.exit()
        if meth == 'full':
            g.SetLineColor(ROOT.kBlue)
        else:
            g.SetLineColor(ROOT.kRed)
        leg.AddEntry(g,meth,'l')
        gd[meth] = g
    for i, meth in enumerate(gd.keys()):
        if i==0:
            if variable == 'tot':
                gd[meth].SetMinimum(0.27)
                gd[meth].SetMaximum(0.31)
                gd[meth].SetTitle('{}; correlation; total uncertainty [GeV]'.format(syst))
            else:
                gd[meth].SetMinimum(round(172.5+mt_smear+rand_offset_small-0.1,1))
                gd[meth].SetMaximum(round(172.5+mt_smear+rand_offset_small+0.1,1))
                gd[meth].SetTitle('{}; correlation; combined mt [GeV]'.format(syst))
            gd[meth].Draw('al')
        else:
            gd[meth].Draw('l same')
    leg.Draw('same')
    c.SaveAs('{}/scan_{}_{}.png'.format(LHC_dir,variable,syst))
    return
