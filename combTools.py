import os, sys, copy
import numpy as np
import ROOT
from ROOT import TH1F, TCanvas, TGraph, TMath, TLegend
from BLUE_object import *

ROOT.gROOT.SetBatch(True)

toys_dir = 'toys_workdir'
scan_dir = 'scan_workdir'
scan_dir_CMS = 'scan_workdir'
scan_dir_LHC = 'LHC_scan_workdir'


def getToyResults(base_obj,l=[],plotToys=True,blind=False):
    l_mt, l_tot, l_stat, l_syst, d_weights, d_syst = base_obj.getToyResults(l)
    
    print('\nvar\tmean\t\trms\tnom')
    if not blind:
        print('mt\t{}\t\t{}\t{}'.format(round(np.array(l_mt).mean(),3),round(np.array(l_mt).std(),3),base_obj.results.mt))
    else:
        print('mt\t{}\t\t{}\t{}'.format(round(np.array(l_mt).mean()/base_obj.results.mt,5),round(np.array(l_mt).std(),3),base_obj.results.mt/base_obj.results.mt))
    print('tot\t{}\t\t{}\t{}'.format(round(np.array(l_tot).mean(),3),round(np.array(l_tot).std(),3),base_obj.results.tot))
    print('stat\t{}\t\t{}\t{}'.format(round(np.array(l_stat).mean(),3),round(np.array(l_stat).std(),3),base_obj.results.stat))
    print('syst\t{}\t\t{}\t{}'.format(round(np.array(l_syst).mean(),3),round(np.array(l_syst).std(),3),base_obj.results.syst))
    print('\n-> weights\n')
    print('meas\t\tmean\trms\tnom')
    for meas in list(d_weights.keys()):
        print('{}\t{}\t{}\t{}'.format(meas,round(np.array(d_weights[meas]).mean(),2),round(np.array(d_weights[meas]).std(),2),base_obj.results.weights[meas]))
    print()
    
    if len(l)==0 and plotToys and not blind:
        plotToyResults(l_mt, l_tot, l_stat, l_syst, d_weights, d_syst, base_obj)
    return

def getToyResultsLHCobj(LHC_base_obj,blind=False):
    l_mt, l_tot, l_stat, l_syst = LHC_base_obj.getToyResults()
    
    print('\nvar\tmean\t\trms\tnom')
    if not blind:
        print('mt\t{}\t\t{}\t{}'.format(round(np.array(l_mt).mean(),3),round(np.array(l_mt).std(),3),LHC_base_obj.getBlueObject().results.mt))
    else:
        print('mt\t{}\t\t{}\t{}'.format(round(np.array(l_mt).mean()/LHC_base_obj.getBlueObject().results.mt,5),round(np.array(l_mt).std(),3),LHC_base_obj.getBlueObject().results.mt/LHC_base_obj.getBlueObject().results.mt))
    print('tot\t{}\t\t{}\t{}'.format(round(np.array(l_tot).mean(),3),round(np.array(l_tot).std(),3),LHC_base_obj.getBlueObject().results.tot))
    print('stat\t{}\t\t{}\t{}'.format(round(np.array(l_stat).mean(),3),round(np.array(l_stat).std(),3),LHC_base_obj.getBlueObject().results.stat))
    print('syst\t{}\t\t{}\t{}'.format(round(np.array(l_syst).mean(),3),round(np.array(l_syst).std(),3),LHC_base_obj.getBlueObject().results.syst))
    print()
    
    return


def plotToyResults(l_mt, l_tot, l_stat, l_syst, d_weights, d_syst, base_obj):

    if not os.path.exists(toys_dir):
        os.makedirs(toys_dir)

    f = open('{}/slides.tex'.format(toys_dir),'w')

    h_mt = TH1F('h_mt','h_mt',70,np.array(l_mt).mean()-3*np.array(l_mt).std(),np.array(l_mt).mean()+3*np.array(l_mt).std())
    h_mt.SetTitle('toys: mt; mt [GeV]; a.u.')
    for t in l_mt:
        h_mt.Fill(t)
    c = TCanvas()
    l = ROOT.TLine(base_obj.results.mt,0,base_obj.results.mt,0.03)
    l.SetLineWidth(2)
    l.SetLineColor(ROOT.kGreen)
    h_mt.DrawNormalized()
    l.Draw('same')
    c.SaveAs('{}/toy_mt.png'.format(toys_dir))
    c.SaveAs('{}/toy_mt.pdf'.format(toys_dir))
    c.Clear()
    f.write('\\begin{{frame}}{{toys for variable mt}}\n')
    f.write('\\centering\\includegraphics[width=.75\\textwidth]{{{}/toy_mt.pdf}}\n'.format(toys_dir))
    f.write('\\end{frame}\n\n')

    h_tot = TH1F('h_tot','h_tot',70,np.array(l_tot).mean()-3*np.array(l_tot).std(),np.array(l_tot).mean()+3*np.array(l_tot).std())
    h_tot.SetTitle('toys: tot; tot [GeV]; a.u.')
    for t in l_tot:
        h_tot.Fill(t)
    c = TCanvas()
    l = ROOT.TLine(base_obj.results.tot,0,base_obj.results.tot,0.03)
    l.SetLineWidth(2)
    l.SetLineColor(ROOT.kGreen)
    h_tot.DrawNormalized()
    l.Draw('same')
    c.SaveAs('{}/toy_tot.png'.format(toys_dir))
    c.SaveAs('{}/toy_tot.pdf'.format(toys_dir))
    c.Clear()
    f.write('\\begin{{frame}}{{toys for variable tot}}\n')
    f.write('\\centering\\includegraphics[width=.75\\textwidth]{{{}/toy_tot.pdf}}\n'.format(toys_dir))
    f.write('\\end{frame}\n\n')

    h_stat = TH1F('h_stat','h_stat',70,np.array(l_stat).mean()-3*np.array(l_stat).std(),np.array(l_stat).mean()+3*np.array(l_stat).std())
    h_stat.SetTitle('toys: stat; stat [GeV]; a.u.')
    for t in l_stat:
        h_stat.Fill(t)
    c = TCanvas()
    l = ROOT.TLine(base_obj.results.stat,0,base_obj.results.stat,0.03)
    l.SetLineWidth(2)
    l.SetLineColor(ROOT.kGreen)
    h_stat.DrawNormalized()
    l.Draw('same')
    c.SaveAs('{}/toy_stat.png'.format(toys_dir))
    c.SaveAs('{}/toy_stat.pdf'.format(toys_dir))
    c.Clear()
    f.write('\\begin{{frame}}{{toys for variable stat}}\n')
    f.write('\\centering\\includegraphics[width=.75\\textwidth]{{{}/toy_stat.pdf}}\n'.format(toys_dir))
    f.write('\\end{frame}\n\n')

    h_syst = TH1F('h_syst','h_syst',70,np.array(l_syst).mean()-3*np.array(l_syst).std(),np.array(l_syst).mean()+3*np.array(l_syst).std())
    h_syst.SetTitle('toys: syst; syst [GeV]; a.u.')
    for t in l_syst:
        h_syst.Fill(t)
    c = TCanvas()
    l = ROOT.TLine(base_obj.results.syst,0,base_obj.results.syst,0.03)
    l.SetLineWidth(2)
    l.SetLineColor(ROOT.kGreen)
    h_syst.DrawNormalized()
    l.Draw('same')
    c.SaveAs('{}/toy_syst.png'.format(toys_dir))
    c.SaveAs('{}/toy_syst.pdf'.format(toys_dir))
    c.Clear()
    f.write('\\begin{{frame}}{{toys for variable syst}}\n')
    f.write('\\centering\\includegraphics[width=.75\\textwidth]{{{}/toy_syst.pdf}}\n'.format(toys_dir))
    f.write('\\end{frame}\n\n')

    if not os.path.exists('{}/syst'.format(toys_dir)):
        os.makedirs('{}/syst'.format(toys_dir))

    for syst in base_obj.usedSyst:
        if base_obj.results.impacts[syst] == 0: continue
        mean = np.array(d_syst[syst]).mean()
        rms = np.array(d_syst[syst]).std()
        h = TH1F(syst,syst,70,mean/base_obj.results.impacts[syst]-3*rms/mean,mean/base_obj.results.impacts[syst]+3*rms/mean)
        h_abs = TH1F(syst+'_abs',syst+'_abs',70,mean-3*rms,mean+3*rms)
        h.SetTitle('{}: normalised to nominal impact; {} / nominal; a.u.'.format(syst,syst))
        h_abs.SetTitle('{}: toys; impact of {} [GeV]; a.u.'.format(syst,syst))
        for j in range(0,len(d_syst[syst])):
            h.Fill(d_syst[syst][j]/base_obj.results.impacts[syst])
            h_abs.Fill(d_syst[syst][j])
        l = ROOT.TLine(1,0,1,0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/syst/{}_rel.png'.format(toys_dir,syst))
        c.SaveAs('{}/syst/{}_rel.pdf'.format(toys_dir,syst))
        l = ROOT.TLine(base_obj.results.impacts[syst],0,base_obj.results.impacts[syst],0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h_abs.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/syst/{}_abs.png'.format(toys_dir,syst))
        c.SaveAs('{}/syst/{}_abs.pdf'.format(toys_dir,syst))
        f.write('\\begin{{frame}}{{toys for systematics {}: absolute and relative variation}}\n'.format(syst))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/syst/{}_abs.pdf}}\n'.format(toys_dir,syst))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/syst/{}_rel.pdf}}\n'.format(toys_dir,syst))
        f.write('\\end{frame}\n\n')

    if not os.path.exists('{}/weights'.format(toys_dir)):
        os.makedirs('{}/weights'.format(toys_dir))


    for meas in base_obj.usedMeas:
        mean = np.array(d_weights[meas]).mean()
        rms = np.array(d_weights[meas]).std()
        h = TH1F(meas,meas,70,mean/base_obj.results.weights[meas]-3*rms/mean,mean/base_obj.results.weights[meas]+3*rms/mean)
        h_abs = TH1F(meas+'_abs',meas+'_abs',70,mean-3*rms,mean+3*rms)
        h.SetTitle('{}: normalised to nominal weight; weight / nominal; a.u.'.format(meas))
        h_abs.SetTitle('{}: weights; weight of measurement {}; a.u.'.format(meas,meas))
        for j in range(0,len(d_weights[meas])):
            h.Fill(d_weights[meas][j]/base_obj.results.weights[meas])
            h_abs.Fill(d_weights[meas][j])

        l = ROOT.TLine(1,0,1,0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/weights/weight_{}_rel.png'.format(toys_dir,meas))
        c.SaveAs('{}/weights/weight_{}_rel.pdf'.format(toys_dir,meas))
        l = ROOT.TLine(base_obj.results.weights[meas],0,base_obj.results.weights[meas],0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h_abs.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/weights/weight_{}_abs.png'.format(toys_dir,meas))
        c.SaveAs('{}/weights/weight_{}_abs.pdf'.format(toys_dir,meas))

        f.write('\\begin{{frame}}{{weight of measurement {}: absolute and relative variation}}\n'.format(meas.replace('_',' ')))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/weights/weight_{}_abs.pdf}}\n'.format(toys_dir,meas))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/weights/weight_{}_rel.pdf}}\n'.format(toys_dir,meas))
        f.write('\\end{frame}\n\n')


    return

def excludeMeasOneByOne(base_obj,blind=False):
    h = TH1F('h','excluding measurements one-by-one; excluded measurement; combined mt [GeV]',len(base_obj.usedMeas), -0.5, -.5+len(base_obj.usedMeas))
    if blind:
        h.GetYaxis().SetTitle('combined mt - nominal mt [GeV]')
    for i, meas in enumerate(base_obj.usedMeas):
        obj = base_obj.clone()
        obj.addExcludeMeas([meas])
        print('excluded:', obj.excludeMeas)
        if not blind:
            obj.simplePrint()
        print()
        central = obj.results.mt
        if blind:
            central -= base_obj.results.mt
        h.SetBinContent(i+1,central)
        h.SetBinError(i+1,obj.results.tot)
        h.GetXaxis().SetBinLabel(i+1,measToROOT(meas))

    if not os.path.exists('plots'):
        os.makedirs('plots')

    ROOT.gStyle.SetOptStat(0000)
    central = base_obj.results.mt if not blind else 0
    l=ROOT.TLine(0,central,len(base_obj.usedMeas)-1,central)
    l.SetLineColor(ROOT.kGreen+1)
    l.SetLineWidth(2)
    c = TCanvas()
    h.Draw('E1')
    l.Draw('same')
    if base_obj.LHC: suffix = 'LHC'
    elif base_obj.ATLAS: suffix = 'ATLAS'
    elif base_obj.CMS: suffix = 'CMS'
    c.SetBottomMargin(0.13)
    h.GetXaxis().SetTitleOffset(1.6)
    if not blind:
        h.GetYaxis().SetTitleOffset(1.3)
    else:
        h.GetYaxis().SetTitleOffset(1.1)
    c.SetBottomMargin(0.13)
    c.SaveAs('plots/excludeOneByOne_{}.png'.format(suffix))
    c.SaveAs('plots/excludeOneByOne_{}.pdf'.format(suffix))

    
    return

def excludeSystOneByOne(base_obj):
    for syst in base_obj.usedSyst:
        if syst == 'Stat': continue
        obj = base_obj.clone()
        obj.addExcludeSyst([syst])
        print('excluded:', obj.excludeSyst)
        obj.simplePrint()
        print()
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
    scan_dir = scan_dir_LHC if base_obj.LHC else scan_dir_CMS
    if not os.path.exists(scan_dir+'/syst'):
        os.makedirs(scan_dir+'/syst')
    step = 0.01
    red_corrs = list(np.arange(0.7,1+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()
    print('\n scanning systematics {} \n'.format(syst_scan))
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

    rf = ROOT.TFile('{}/syst/{}.root'.format(scan_dir,syst_scan),'recreate')
    h_mt.SetName('mt')
    h_tot.SetName('tot')
    h_stat.SetName('stat')
    h_syst.SetName('syst')
    h_mt.Write()
    h_tot.Write()
    h_stat.Write()
    h_syst.Write()
    rf.Close()

    return

def plotScanSummary(base_obj,blind=False,syst_list=[]):

    scan_dir = scan_dir_LHC if base_obj.LHC else scan_dir_CMS
    
    mt_maxdiff = {}
    tot_maxdiff = {}
    mt_m = 0
    mt_M = 0
    tot_m = 0
    tot_M = 0

    if syst_list == []:
        syst_list = base_obj.usedSyst

    for i, syst in enumerate(syst_list):
        if syst == 'Stat': continue
        rf = ROOT.TFile('{}/syst/{}.root'.format(scan_dir,syst),'read')
        mt = rf.Get('mt')
        tot = rf.Get('tot')
        mt_diff = TMath.MaxElement(mt.GetN(),mt.GetY()) - TMath.MinElement(mt.GetN(),mt.GetY())
        tot_diff = TMath.MaxElement(tot.GetN(),tot.GetY()) - TMath.MinElement(tot.GetN(),tot.GetY())
        if not blind:
            pass_mt = mt_diff/base_obj.results.mt > 0.00004
        else:
            pass_mt = mt_diff > 0.003
        if pass_mt:
            mt_maxdiff[syst] = mt_diff
            if len(list(mt_maxdiff.keys()))==1:
                mt_m = TMath.MinElement(mt.GetN(),mt.GetY())
                mt_M = TMath.MaxElement(mt.GetN(),mt.GetY())
            else:
                if TMath.MinElement(mt.GetN(),mt.GetY()) < mt_m:
                    mt_m = TMath.MinElement(mt.GetN(),mt.GetY())
                if TMath.MaxElement(mt.GetN(),mt.GetY()) > mt_M:
                    mt_M = TMath.MaxElement(mt.GetN(),mt.GetY())
        if not base_obj.LHC:
            pass_tot = tot_diff/base_obj.results.tot > 0.006
        else:
            pass_tot = tot_diff/base_obj.results.tot > 0.006
        if pass_tot:
            tot_maxdiff[syst] = tot_diff
            if len(list(tot_maxdiff.keys()))==1:
                tot_m = TMath.MinElement(tot.GetN(),tot.GetY())
                tot_M = TMath.MaxElement(tot.GetN(),tot.GetY())
            else:
                if TMath.MinElement(tot.GetN(),tot.GetY()) < tot_m:
                    tot_m = TMath.MinElement(tot.GetN(),tot.GetY())
                if TMath.MaxElement(tot.GetN(),tot.GetY()) > tot_M:
                    tot_M = TMath.MaxElement(tot.GetN(),tot.GetY())

    c = TCanvas()
    leg = TLegend(.5,.6,.87,.87) if not base_obj.LHC else TLegend(.35,.7,.62,.87)
    leg.SetBorderSize(0)
    from systNameDict import systNameDict

    list_mt = mt_maxdiff.keys()
    list_tot = tot_maxdiff.keys()
    if base_obj.LHC:
        list_mt = list_tot = [s for s in list_mt if s in list_tot]

    for i,syst in enumerate(list_mt):
        rf = ROOT.TFile('{}/syst/{}.root'.format(scan_dir,syst),'read')
        mt = rf.Get('mt')
        mt.SetMarkerColor(i+1)
        mt.SetLineColor(i+1)
        mt.SetLineWidth(2)
        leg.AddEntry(mt,systNameDict[syst].replace('\\',''),'pl')
        if i==0:
            mt.SetTitle('mt dependence on correlations')
            mt.GetYaxis().SetTitleOffset(1.8)
            mt.GetYaxis().SetRangeUser(mt_m-.003,mt_M+.003)
            mt.Draw('apl')
        else:
            mt.Draw('pl same')
    leg.Draw('same')
    c.SetLeftMargin(0.15)
    if base_obj.LHC:
        c.SaveAs(scan_dir+'/scan_mt_summary_LHC.png')
        c.SaveAs(scan_dir+'/scan_mt_summary_LHC.pdf')
    else:
        c.SaveAs(scan_dir+'/scan_mt_summary.png')
        c.SaveAs(scan_dir+'/scan_mt_summary.pdf')

    c.Clear()
    leg.Clear()

    for i,syst in enumerate(list_tot):
        rf = ROOT.TFile('{}/syst/{}.root'.format(scan_dir,syst),'read')
        tot = rf.Get('tot')
        tot.SetMarkerColor(i+1)
        tot.SetLineColor(i+1)
        tot.SetLineWidth(2)
        leg.AddEntry(tot,systNameDict[syst],'pl')
        if i==0:
            tot.SetTitle('total uncertainty - dependence on correlations')
            tot.GetYaxis().SetRangeUser(tot_m-.001,tot_M+.005)
            tot.GetYaxis().SetTitleOffset(1.8)
            tot.Draw('apl')
        else:
            tot.Draw('pl same')
    leg.Draw('same')
    if base_obj.LHC:
        c.SaveAs(scan_dir+'/scan_tot_summary_LHC.png')
        c.SaveAs(scan_dir+'/scan_tot_summary_LHC.pdf')
    else:
        c.SaveAs(scan_dir+'/scan_tot_summary.png')
        c.SaveAs(scan_dir+'/scan_tot_summary.pdf')
        
    leg = TLegend(.4,.6,.8,.85)
    leg.SetBorderSize(0)

    if base_obj.LHC:
        c.Clear()
        c.Divide(1,2)
        p1 = c.cd(1)
        p1.SetBottomMargin(0.01)
        p1.SetRightMargin(0.015)
        p1.SetLeftMargin(0.1)
        for i,syst in enumerate(list_mt):
            rf = ROOT.TFile('{}/syst/{}.root'.format(scan_dir,syst),'read')
            mt = rf.Get('mt')
            mt.SetMarkerColor(i+1)
            mt.SetLineColor(i+1)
            mt.SetLineWidth(2)
            leg.AddEntry(mt,systNameDict[syst],'l')
            if i==0:
                mt.GetYaxis().SetRangeUser(mt_m-.008,mt_M+.005)
                mt.GetYaxis().SetTitleSize(mt.GetYaxis().GetTitleSize()*1.7)
                mt.GetYaxis().SetTitleOffset(.8)
                mt.GetYaxis().SetTitle('m_{t} - m_{t}^{central}')
                mt.GetYaxis().SetLabelSize(tot.GetYaxis().GetLabelSize()*1.5)
                mt.GetXaxis().SetLabelSize(0)
                mt.Draw('apl')
            else:
                mt.Draw('pl same')
        leg.Draw('same')

        p2 = c.cd(2)
        p2.SetBottomMargin(0.14)
        p2.SetTopMargin(0.05)
        p2.SetRightMargin(0.015)
        p2.SetLeftMargin(0.1)
        for i,syst in enumerate(list_tot):
            rf = ROOT.TFile('{}/syst/{}.root'.format(scan_dir,syst),'read')
            tot = rf.Get('tot')
            tot.SetMarkerColor(i+1)
            tot.SetLineColor(i+1)
            tot.SetLineWidth(2)
            if i==0:
                tot.GetYaxis().SetRangeUser(tot_m-.001,tot_M+.005)
                tot.GetYaxis().SetTitleSize(tot.GetYaxis().GetTitleSize()*1.7)
                tot.GetYaxis().SetTitleOffset(.8)
                tot.GetYaxis().SetTitle('total uncertainty')
                tot.GetXaxis().SetLabelSize(tot.GetXaxis().GetLabelSize()*1.5)
                tot.GetYaxis().SetLabelSize(tot.GetYaxis().GetLabelSize()*1.5)
                tot.GetXaxis().SetTitleSize(tot.GetXaxis().GetTitleSize()*1.7)
                tot.Draw('apl')
            else:
                tot.Draw('pl same')

        c.SaveAs(scan_dir+'/scan_full_summary.png')
        c.SaveAs(scan_dir+'/scan_full_summary.pdf')

    return

def makeCorrelationScans(base_obj,blind=False):

    scan_dir = scan_dir_LHC if base_obj.LHC else scan_dir_CMS

    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)

    if not blind:
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

def makeCorrelationScanLHC(base_obj,syst_scan,brutal=False):

    step = 0.01
    orig_corr = round(base_obj.matrix[syst_scan]['CMS_comb']['ATLAS_comb'],1)
    if not brutal:
        corrs = list(np.arange(max(0,orig_corr -.5),min(1,orig_corr+.5)+step/2,step))
    else:
        corrs = list(np.arange(-1.*min(1,orig_corr+.5),min(1,orig_corr+.5)+step/2,step))

    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()

    print('\n scanning systematics {} \n'.format(syst_scan))

    for i,corr in enumerate(corrs):
        obj = base_obj.clone()
        obj.setCorrelationLHC(corr,syst_scan)
        if not blind:
            h_mt.SetPoint(i,corr,obj.results.mt)
        else:
            h_mt.SetPoint(i,corr,obj.results.mt-base_obj.results.mt)
        h_tot.SetPoint(i,corr,obj.results.tot)
        h_stat.SetPoint(i,corr,obj.results.stat)
        h_syst.SetPoint(i,corr,obj.results.syst)
        
    h_mt.SetTitle('mt scan - {}; correlation; mt [GeV]'.format(syst_scan))
    h_tot.SetTitle('tot uncert scan - {}; correlation; tot uncert [GeV]'.format(syst_scan))
    h_stat.SetTitle('stat uncert scan - {}; correlation; stat uncert [GeV]'.format(syst_scan))
    h_syst.SetTitle('syst uncert scan - {}; correlation; syst uncert [GeV]'.format(syst_scan))

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    href_mt = TGraph()
    href_tot = TGraph()
    href_stat = TGraph()
    href_syst = TGraph()

    if not blind:
        href_mt.SetPoint(0,orig_corr,base_obj.results.mt)
    else:
        href_mt.SetPoint(0,orig_corr,0)
    href_tot.SetPoint(0,orig_corr,base_obj.results.tot)
    href_stat.SetPoint(0,orig_corr,base_obj.results.stat)
    href_syst.SetPoint(0,orig_corr,base_obj.results.syst)

    href_mt.SetMarkerStyle(8)
    href_tot.SetMarkerStyle(8)
    href_stat.SetMarkerStyle(8)
    href_syst.SetMarkerStyle(8)

    href_mt.SetMarkerColor(ROOT.kRed)
    href_tot.SetMarkerColor(ROOT.kRed)
    href_stat.SetMarkerColor(ROOT.kRed)
    href_syst.SetMarkerColor(ROOT.kRed)


    c = TCanvas()
    h_mt.Draw('apl')
    href_mt.Draw('p same')
    c.SaveAs('{}/mt_{}.pdf'.format(scan_dir_LHC,syst_scan))
    c.SaveAs('{}/mt_{}.png'.format(scan_dir_LHC,syst_scan))
    c.Clear()

    h_tot.Draw('apl')
    href_tot.Draw('p same')
    c.SaveAs('{}/tot_{}.pdf'.format(scan_dir_LHC,syst_scan))
    c.SaveAs('{}/tot_{}.png'.format(scan_dir_LHC,syst_scan))
    c.Clear()

    h_stat.Draw('apl')
    href_stat.Draw('p same')
    c.SaveAs('{}/stat_{}.pdf'.format(scan_dir_LHC,syst_scan))
    c.SaveAs('{}/stat_{}.png'.format(scan_dir_LHC,syst_scan))
    c.Clear()

    h_syst.Draw('apl')
    href_syst.Draw('p same')
    c.SaveAs('{}/syst_{}.pdf'.format(scan_dir_LHC,syst_scan))
    c.SaveAs('{}/syst_{}.png'.format(scan_dir_LHC,syst_scan))
    c.Clear()

    return


def makeCorrelationScansLHC(base_obj,brutal=False):
    if not os.path.exists(scan_dir_LHC):
        os.makedirs(scan_dir_LHC)

    f = open('{}/slides.tex'.format(scan_dir_LHC),'w')

    for syst in list(base_obj.matrix.keys()):
        makeCorrelationScanLHC(base_obj,syst,brutal)
        f.write('\\begin{{frame}}{{correlation scan for systematics {}}}\n'.format(syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/mt_{}.pdf}}\n'.format(scan_dir_LHC,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/tot_{}.pdf}}\\\\\n'.format(scan_dir_LHC,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/stat_{}.pdf}}\n'.format(scan_dir_LHC,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst_{}.pdf}}\n'.format(scan_dir_LHC,syst))
        f.write('\\end{frame}\n\n')

    return

def drawWeights(base_obj,path='plots'):
    h = TH1F('h','LHC combination ; ; measurement weight',len(base_obj.usedMeas), -0.5, -.5+len(base_obj.usedMeas))
    for i, meas in enumerate(base_obj.usedMeas):
        h.Fill(i,base_obj.results.weights[meas])
        if base_obj.CMS:
            h.GetXaxis().SetBinLabel(i+1,nameForCMSPlots(meas))
        else:
            h.GetXaxis().SetBinLabel(i+1,measToROOT(meas))

    if not os.path.exists(path):
        os.makedirs(path)

    ROOT.gStyle.SetOptStat(0000)
    c = TCanvas()
    h.SetFillColor(ROOT.kYellow)
    h.Draw('hist')
    c.SaveAs('{}/weights_{}.png'.format(path,base_obj.experiment()))
    c.SaveAs('{}/weights_{}.pdf'.format(path,base_obj.experiment()))


    return

def nameForCMSPlots(name):
    if 'CMS11' in name:
        name += ' 7TeV'
    elif 'CMS12' in name:
        name += ' 8TeV'
    return name.replace('CMS11_','').replace('CMS12_','')


def measToTex(meas):
    if 'dil' in meas or 'MT2' in meas: return '$ll$'
    elif 'lj' in meas: return '$lj$'
    elif 'allhad' in meas or 'aj' in meas: return '$aj$'
    elif 'sto' in meas: return '$t$'
    elif 'SVX' in meas: return '$vtx$'
    elif 'jps' in meas: return '$J/\\psi$'
    elif meas == 'ATLAS_comb' or meas == 'CMS_comb': return meas.replace('_comb','')
    else: return 'ERROR'

def measToROOT(meas):
    name = 'CMS ' if 'CMS' in meas else 'ATLAS '
    if 'dil' in meas or 'MT2' in meas: name+='ll '
    elif 'lj' in meas: name+='lj '
    elif 'allhad' in meas or 'aj' in meas: name+='aj '
    elif 'sto' in meas: name+='t '
    elif 'SVX' in meas: name+='vtx '
    elif 'jps' in meas: name+='J/#Psi '
    else: return 'ERROR'
    if 'CMS' in name:
        if '11' in meas:
            return name+'7 TeV'
        else:
            return name+'8 TeV'
    else:
        if '7' in meas:
            return name+'7 TeV'
        else:
            return name+'8 TeV'


def removeUselessCharachters(name):
    name = name.replace(' ','_')
    if name.endswith('_'):
        name = name[:-1]
    if name.endswith('_'):
        name = removeUselessCharachters(name)
    if name.startswith('_'):
        name = name[1:]
    if name.startswith('_'):
        name = removeUselessCharachters(name)
    return name

def isSymmetricMatrix(matrix):
    if len(matrix[0]) != len(matrix[:][0]) : return False
    for i in range(0,len(matrix[0])):
        for j in range(0,len(matrix[0])):
            if matrix[i][j] != matrix[j][i]: return False
    return True

def isPositiveDefinite(m):
    if not isInvertible(m):
        return False
    w,v = np.linalg.eig(m)
    return (w > 0).all()

def isNonNegativeDefinite(m):
    if not isInvertible(m):
        return False
    w,v = np.linalg.eig(m)
    return (w >= 0).all()

def isInvertible(m):
    if np.isfinite(np.linalg.cond(m)):
        return True
    else:
        return False
