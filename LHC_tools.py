import os, sys, copy
import numpy as np
import ROOT
from ROOT import TCanvas, TGraph, TLegend, TLatex, TFile
from BLUE_object import BLUE_object
from LHC_object import LHC_object

ROOT.gROOT.SetBatch(True)
np.random.seed(1)
from combTools import scan_dir_LHC

def makeAllCorrelationScansLHC(full,sep,blind=False):

    os.makedirs(scan_dir_LHC+'/syst',exist_ok=True)

    LHC_full = full.clone()
    LHC_sep = sep.clone()
    orig_corrMap = copy.deepcopy(LHC_full.corrMap)
    if LHC_sep.corrMap != orig_corrMap:
        print('ERROR: correlation maps should be the same for both methods')
        sys.exit()
    for syst in LHC_full.LHCsyst:
        if syst == 'Stat': continue
        makeCorrelationScan(LHC_full,LHC_sep,orig_corrMap,syst,blind)
    return

def makeCorrelationScan(LHC_full,LHC_sep,corrMap,syst,blind=True):
    print('scanning', syst)

    obj_d = {'full' : LHC_full, 'separate' : LHC_sep}

    orig_corrMap = copy.deepcopy(corrMap)
    if not syst in orig_corrMap:
        orig_corrMap[syst] = 0.

    corr = orig_corrMap[syst]
    step = 0.005 if syst == 'JESFLV' else 0.05
    halfrange = min(0.25, 1-abs(corr))
    corrs = list(np.arange(corr-halfrange,corr+halfrange+step/2,step))
    corrs = [round(c,3) for c in corrs]

    scan_d = dict()
    for meth, obj in list(obj_d.items()):
        result_l = scanOne(obj,orig_corrMap,syst,corrs)
        if blind:
            for i, corr in enumerate(corrs):
                result_l[i].mt -= obj.BLUE_obj.results.mt
        scan_d[meth] = result_l

    of = '{}/syst/{}.root'.format(scan_dir_LHC,syst)
    f = TFile(of,'recreate')

    g_tot = plotScanResults(corrs,scan_d,syst,'tot')
    g_mt = plotScanResults(corrs,scan_d,syst,'mt',blind)

    g_tot.Write()
    g_mt.Write()

    f.Close()

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

def plotScanResults(corrs,scan_d,syst,variable,blind=True):
    c = TCanvas()
    c.SetLeftMargin(0.15)
    leg = TLegend(.2,.15,.45,.3)
    leg.SetBorderSize(0)
    gd = dict()
    for meth, scan in list(scan_d.items()):
        g = TGraph()
        g.SetName(meth)
        g.GetYaxis().SetTitleOffset(1.7)
        for i, corr in enumerate(corrs):
            if variable == 'tot':
                g.SetPoint(i,corr,scan[i].tot)
            elif variable == 'mt':
                if blind:
                    g.SetPoint(i,corr,scan[i].mt)
            else:
                print('ERROR: variable {} not supported'.format(variable))
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
                gd[meth].SetMinimum(0.28)
                gd[meth].SetMaximum(0.32)
                gd[meth].SetTitle('{}; correlation; total uncertainty [GeV]'.format(syst))
            else:
                offset = 172.5 if not blind else 0
                gd[meth].SetMinimum(offset-0.04)
                gd[meth].SetMaximum(offset+0.04)
                if not blind:
                    gd[meth].SetTitle('{}; correlation; combined mt [GeV]'.format(syst))
                else:
                    gd[meth].SetTitle('{}; correlation; combined mt - nominal mt [GeV]'.format(syst))
            gd[meth].Draw('al')
        else:
            gd[meth].Draw('l same')
    leg.Draw('same')
    latexLabel = TLatex()
    from systNameDict import systNameDict
    latexLabel.SetTextSize(0.045)
    latexLabel.DrawLatexNDC(.15,.92,'scan for: {}'.format(systNameDict[syst]))
    c.SaveAs('{}/scan_{}_{}.png'.format(scan_dir_LHC,variable,syst))

    gd['full'].SetName(variable)

    return gd['full']

def flipAmbiguousSigns(full,sep):

    if full.blind or sep.blind:
        print('*******')
        print('WARNING: one or more input objects are blind - mt results will be meaningless')
        print('*******')

    LHC_full = full.clone()
    LHC_sep = sep.clone()
    orig_corrMap = copy.deepcopy(LHC_full.corrMap)

    if LHC_sep.corrMap != orig_corrMap:
        print('ERROR: correlation maps should be the same for both methods')
        sys.exit()
    if LHC_full.noSignsOnImpacts != LHC_sep.noSignsOnImpacts:
        print('ERROR: list of ambiguous signs different in the two combinations')
        sys.exit()
    
    ambiguous_l = LHC_full.noSignsOnImpacts['ATLAS']
    for amb in LHC_full.noSignsOnImpacts['CMS']:
        if not amb in ambiguous_l:
            ambiguous_l.append(amb)

    flipAllSignsLHC(LHC_full,LHC_sep,ambiguous_l,orig_corrMap)
    for syst in ambiguous_l:
        flipSignLHC(LHC_full,LHC_sep,syst,orig_corrMap)
    print()

    return

def flipAllSignsLHC(LHC_full,LHC_sep,ambiguous_l,orig_corrMap):

    corrMap = copy.deepcopy(orig_corrMap)
    full = LHC_full.clone()
    sep = LHC_sep.clone()

    obj_d = {'full': full, 'separate': sep}
    for syst in ambiguous_l:
        if syst in list(corrMap.keys()):
            corrMap[syst] *= -1
    print()
    print()
    print('*sign flip for all ambiguous systematics')

    for meth, obj in list(obj_d.items()):
        print()
        print('-> method =', meth)
        print('uncertainty original signs =', obj.getBlueObject().results.tot, 'GeV')
        mt_orig = obj.getBlueObject().results.mt
        obj.setNewLHCcorrMap(corrMap)
        print('uncertainty flipped signs =', obj.getBlueObject().results.tot, 'GeV')
        print('mt(flip) - mt(original) =', round(obj.getBlueObject().results.mt-mt_orig,3) , 'GeV')

def flipSignLHC(LHC_full,LHC_sep,syst,orig_corrMap):

    if not syst in list(orig_corrMap.keys()):
        print('\nWARNING: systematics {} (for sign flip) not in correlation map. Nothing done'.format(syst))
        return

    corrMap = copy.deepcopy(orig_corrMap)
    full = LHC_full.clone()
    sep = LHC_sep.clone()

    obj_d = {'full': full, 'separate': sep}
    corrMap[syst] *= -1
    
    print()
    print()
    print('*sign flip for systematics {}*'.format(syst))

    
    for meth, obj in list(obj_d.items()):
        print()
        print('-> method =', meth)
        print('uncertainty original signs =', obj.getBlueObject().results.tot, 'GeV')
        mt_orig = obj.getBlueObject().results.mt
        obj.setNewLHCcorrMap(corrMap)
        print('uncertainty flipped signs =', obj.getBlueObject().results.tot, 'GeV')
        print('mt(flip) - mt(original) =', round(obj.getBlueObject().results.mt-mt_orig,3) , 'GeV')

    return
