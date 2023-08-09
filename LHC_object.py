from BLUE_object import BLUE_object, sqrts
from combTools import measToTex, measToROOT
import copy
import sys,os
import numpy as np
import itertools
from operator import itemgetter
import systNameDict as snd

import ROOT as rt
from ROOT import TH2D, TCanvas, gStyle, TGraphErrors, TLatex

from array import *

corrMap_default = {'JES3': 0.5, 'JESFLV': 0.85, 'RAD': 0.5, 'MCGEN': 0.5, 'BKMC': .85, 'PDF': .85 , 'BTAG': 0.5, 'UE': .85, 'CR': .85, 'JESflavresLHC': .85, 'HADR':-.5}
mergeMap_default = {'ATLAS':{}, 'CMS': {'RAD': ['Q','JPS']}}
renameMap_default = {'ATLAS':{'bJES':'JESFLV', 'JESflavres':'JESflavresLHC', 'JESflavcomp':'JESlight', 'Other':'Extra'} ,
                     'CMS': {'JES5':'JESFLV', 'JES4':'JESflavresLHC', 'BFRAG':'HADR', 'JES6': 'JESlight'}}

noSignsOnImpacts = {'ATLAS':['BKMC', 'BTAG', 'PDF'], 'CMS': []}

mergeImpacts_default = {}
# mergeImpacts_default = {'JESlight':['JES6','JESflavresLHC','JESflavcomp'],'HADR':['HADR','SLEPB']}

tab_dir = 'tables'
plot_dir = 'plots'

class LHC_object:

    def __init__(self,ATLAS_obj, CMS_obj, excludeSyst = [], separateCombinations = False, corrMap = None, mergeMap = None, renameMap = None, blind = False, 
                 merge_and_rename=True,mergeImpacts=mergeImpacts_default,PU_hack=False):
        self.ATLAS_obj = ATLAS_obj.clone()
        self.CMS_obj = CMS_obj.clone()
        self.obj_d = {'ATLAS':self.ATLAS_obj, 'CMS':self.CMS_obj}

        self.excludeSyst = copy.deepcopy(excludeSyst)
        if len(self.excludeSyst) > 0:
            for obj in list(self.obj_d.values()):
                obj.addExcludeSyst(self.excludeSyst)
            
        self.experiments = list(self.obj_d.keys())
        self.blind = blind
        if blind:
            for obj in list(self.obj_d.values()):
                obj.makeBlind()
            
        self.separateCombinations = separateCombinations
        self.LHC_file = 'LHC_input_sep.txt' if self.separateCombinations else 'LHC_input_full.txt'
        self.merge_and_rename = merge_and_rename
        self.mergeImpacts = mergeImpacts
        self.PU_hack = PU_hack

        # self.removeZeroImpacts()
        # for exp in self.experiments:
        #     print '->', exp
        #     self.obj_d[exp].simplePrint()

        self.noSignsOnImpacts = copy.deepcopy(noSignsOnImpacts)

        if corrMap is None:  self.corrMap = copy.deepcopy(corrMap_default)
        else: self.corrMap = copy.deepcopy(corrMap)

        if self.PU_hack and not 'PU' in list(self.corrMap.keys()):
            print('\nWARNING: adding PU correlation of 0.85\n')
            self.corrMap['PU'] = 0.85

        if mergeMap is None: self.mergeMap = copy.deepcopy(mergeMap_default)
        else: self.mergeMap = copy.deepcopy(mergeMap)
        if renameMap is None: self.renameMap = copy.deepcopy(renameMap_default)
        else: self.renameMap = copy.deepcopy(renameMap)

        if self.merge_and_rename:
            self.renameAllSyst()
            self.mergeAllSyst()

        if len(self.mergeImpacts) > 0:
            self.ATLAS_obj.setMergeImpacts(self.mergeImpacts)
            self.CMS_obj.setMergeImpacts(self.mergeImpacts)

        self.commonSyst = self.getCommonSyst()
        self.ATLAS_only, self.CMS_only, self.LHCsyst = self.prepareLHCcombination()

        self.redefineNegativeSignsLHC()

        if self.separateCombinations:
            for obj in list(self.obj_d.values()):
                if len(list(obj.results.signedImpacts.keys())) == 0:
                    obj.deriveSignedImpacts()
        else:
            self.LHCmeas = self.getLHCmeas()
            self.LHCmatrix = self.createLHCmatrix()

        for exp, obj in list(self.obj_d.items()):
            for syst in self.noSignsOnImpacts[exp]:
                if syst in obj.uncertWithSign:
                    print('ERROR: something wrong in signs of syst {} in {}'.format(syst,exp))
                    sys.exit()
            for syst in obj.uncertWithSign:
                if syst in self.noSignsOnImpacts[exp]:
                    print('ERROR: something wrong in signs of syst {} in {}'.format(syst,exp))
                    sys.exit()

        self.writeBLUEinputCMS(separateCombinations=self.separateCombinations)
        self.BLUE_obj = BLUE_object(self.LHC_file,LHC=True,blind=self.blind,mergeImpacts=self.mergeImpacts,PU_hack=self.PU_hack)

        self.BLUE_obj.renamed = {**self.obj_d['ATLAS'].renamed,**self.obj_d['CMS'].renamed}
        for renamed in self.BLUE_obj.renamed.keys():
            if list(self.BLUE_obj.renamed.keys()).count(renamed) > 1:
                print('ERROR: ambiguity')
                sys.exit()

        # self.BLUE_obj.checkAllSystMatrices()
        if not self.separateCombinations:
            self.sortUsedMeas()


    def clone(self):
        return copy.deepcopy(self)

    def update(self):
        self.writeBLUEinputCMS(separateCombinations=self.separateCombinations)
        self.BLUE_obj = BLUE_object(self.LHC_file,LHC=True,blind=self.blind,mergeImpacts=self.mergeImpacts,PU_hack=self.PU_hack)


    def renameAllSyst(self):
        print()
        for exp in self.experiments:
            for old, new in list(self.renameMap[exp].items()):
                print('{}: renaming syst {} to {}'.format(exp,old,new))
                self.obj_d[exp].renameSyst(old,new)
        return

    def mergeAllSyst(self):
        print()
        for exp in self.experiments:
            for merged, original_l in list(self.mergeMap[exp].items()):
                print('{}: creating new syst {} from sources {}'.format(exp,merged,original_l))
                signs_propagated = self.obj_d[exp].mergeSyst(merged,original_l)
                self.noSignsOnImpacts[exp].append(merged)
        return


    def removeZeroImpacts(self):
        for obj in list(self.obj_d.values()):
            for syst in obj.usedSyst:
                if obj.results.impacts[syst] == 0:
                    obj.addExcludeSyst([syst])
        return

    def getCommonSyst(self):
        l = []
        for syst in self.ATLAS_obj.usedSyst:
            if syst in self.CMS_obj.usedSyst:
                l.append(syst)
        return l

    def prepareLHCcombination(self):
        ATLAS_only = list(set(self.ATLAS_obj.usedSyst) - set(self.commonSyst))
        CMS_only = list(set(self.CMS_obj.usedSyst) - set(self.commonSyst))
        LHCsyst = ATLAS_only + CMS_only + self.commonSyst
        
        if self.merge_and_rename:
            print()
            print('-> unique to ATLAS')
            print(ATLAS_only)
            print()
            print('-> unique to CMS')
            print(CMS_only)
            print()

        return ATLAS_only, CMS_only, LHCsyst

    def redefineNegativeSignsLHC(self):
        for syst, corr in self.corrMap.items():
            if corr < 0:
                print('WARNING: redifining sign for syst {}'.format(syst))
                self.CMS_obj.flipAllSignsSyst(syst)
                self.corrMap[syst] *= -1
        return

    def getLHCmeas(self):
        allMeas = []
        for exp, obj in list(self.obj_d.items()):
            for meas in obj.usedMeas:
                allMeas.append('{}_{}'.format(meas,exp))
        return allMeas

    def createLHCmatrix(self):
        LHCmatrix = dict()
        for syst in self.LHCsyst:
            LHCmatrix[syst] = self.getLHCsystDict(syst)
        return LHCmatrix

    def getLHCsystDict(self,syst):
        m_d = dict()
        for meas1 in self.LHCmeas:
            m_dd = dict()
            for meas2 in self.LHCmeas:
                if meas1 == meas2:
                    corr = 1.
                else:
                    if '_ATLAS' in meas1 and '_ATLAS' in meas2:
                        if syst in list(self.obj_d['ATLAS'].matrix.keys()):
                            corr = self.obj_d['ATLAS'].matrix[syst][meas1.replace('_ATLAS','')][meas2.replace('_ATLAS','')]
                            if syst in self.noSignsOnImpacts['ATLAS']:
                                pass
                            elif syst in self.ATLAS_obj.uncertWithSign:
                                corr = abs(corr)

                        else:
                            corr = 0.
                    elif '_CMS' in meas1 and '_CMS' in meas2:
                        if syst in list(self.obj_d['CMS'].matrix.keys()):
                            corr = self.obj_d['CMS'].matrix[syst][meas1.replace('_CMS','')][meas2.replace('_CMS','')]
                            if not syst in self.noSignsOnImpacts['CMS']:
                                corr = abs(corr)
                        else:
                            corr = 0.
                    else:
                        if syst in list(self.corrMap.keys()):
                            corr = self.corrMap[syst]
                        else:
                            corr = 0.
                    # if abs(corr) > 0.85:
                    #     corr = 0.85 * (corr/abs(corr))
                m_dd[meas2] = corr
            m_d[meas1] = m_dd

        if syst == 'PU' and self.PU_hack:
            for meas1 in self.LHCmeas:
                for meas2 in self.LHCmeas:
                    if sqrts(meas1) != sqrts(meas2):
                        m_d[meas1][meas2] = 0.0

        return m_d

    def writeSeparateCombinationsInput(self,f):

        f.write('\./combine<<!\n\'LHC Top Mass combination\'\n-1 -1 0 internal & minuit debug level, dependency flag\n')
        f.write('1 2 {}  # of observables, measurements, error classes\n'.format(len(self.LHCsyst)))
        f.write('\'Mtop\'     name of observable\n\n')

        syst_l = copy.deepcopy(self.LHCsyst)
        syst_l.remove('Stat')
        syst_l = ['Stat'] + sorted(syst_l)

        for syst in syst_l:
            f.write('\'{}\'\t'.format(syst))
        f.write('\n')


        for exp in self.experiments:
            f.write('\'{}\' \'Mtop\' {}'.format(exp+' comb',self.obj_d[exp].results.mt))
            for syst in syst_l:
                if syst in list(self.obj_d[exp].results.impacts.keys()):
                    factor = 1.
                    if syst != 'Stat' and self.obj_d[exp].results.signedImpacts[syst] < 0:
                        factor = -1.
                    f.write(' {}'.format(factor*self.obj_d[exp].results.impacts[syst]))
                else:
                    f.write(' 0.0')
            f.write('\n')

        f.write('\n')

        used = []
        for syst in list(self.corrMap.keys()):
            if not syst in syst_l and not syst in self.excludeSyst:
                print('ERROR: systematics {} (in correlation map) not found in LHC grid'.format(syst))
                sys.exit()
        for syst in syst_l:
            for i,meas1 in enumerate(self.experiments):
                for j,meas2 in enumerate(self.experiments):
                    corr = 0.0
                    if syst in list(self.corrMap.keys()):
                        corr = self.corrMap[syst]
                        if not syst in used:
                            used.append(syst)
                    if i==j:
                        corr = 1.0
                    f.write('{} '.format(corr))
                if i == 0:
                    f.write('\'{}\''.format(syst))
                f.write('\n')
            f.write('\n')

        for syst in list(self.corrMap.keys()):
            if not syst in used and not syst in self.excludeSyst:
                print('ERROR: systematic {} (in correlation map) not used'.format(syst))
                sys.exit()

        f.write('!\n')
        return
        
    def writeSimultaneousCombinationInput(self,f):
        f.write('\./combine<<!\n\'LHC Top Mass combination\'\n-1 -1 0 internal & minuit debug level, dependency flag\n')
        f.write('1 {} {}  # of observables, measurements, error classes\n'.format(len(self.ATLAS_obj.usedMeas)+len(self.CMS_obj.usedMeas),len(self.LHCsyst)))
        f.write('\'Mtop\'     name of observable\n\n')

        syst_l = copy.deepcopy(self.LHCsyst)
        syst_l.remove('Stat')
        syst_l = ['Stat'] + sorted(syst_l)

        for syst in syst_l:
            f.write('\'{}\' '.format(syst))
        f.write('\n')

        for meas in self.LHCmeas:
            if '_ATLAS' in meas:
                obj = self.obj_d['ATLAS']
                meas_exp = meas.replace('_ATLAS','')
            else:
                obj = self.obj_d['CMS']
                meas_exp = meas.replace('_CMS','')

            f.write('\'{}\' \'Mtop\' {}'.format(meas_exp.replace('_',' '), obj.value[meas_exp]))

            for syst in syst_l:
                if syst in list(obj.matrix.keys()) or syst == 'Stat':
                    f.write(' {}'.format(obj.uncert[meas_exp][syst]))
                else:
                    f.write(' 0.0')
            f.write('\n')
        f.write('\n')


        for syst in syst_l:
            for i,meas1 in enumerate(self.LHCmeas):
                for meas2 in self.LHCmeas:
                    f.write('{} '.format(self.LHCmatrix[syst][meas1][meas2]))
                if i == 0:
                    f.write('\'{}\''.format(syst))
                f.write('\n')
            f.write('\n')
        f.write('!\n')

        return


    def writeBLUEinputCMS(self,tmprun=False,separateCombinations=True):

        if tmprun:
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)
            f = open('{}/LHC_combination_tmp.txt'.format(tmp_dir),'w')
        else:
            f = open(self.LHC_file,'w')

        if separateCombinations:
            self.writeSeparateCombinationsInput(f)
        else:
            self.writeSimultaneousCombinationInput(f)

        return

    def getBlueObject(self):
        return self.BLUE_obj

    def printResults(self):
        self.BLUE_obj.printResults()

    def setNewLHCcorrMap(self,corrMap):
        self.corrMap = corrMap
        if not self.separateCombinations:
            self.LHCmatrix = self.createLHCmatrix()
        self.update()

    def updateCorr(self,syst,corr):
        self.corrMap[syst] = corr
        if not self.separateCombinations:
            self.LHCmatrix = self.createLHCmatrix()
        self.update()

    def simplePrint(self):
        self.BLUE_obj.simplePrint()

    def printImpactsSorted(self):
        self.BLUE_obj.printImpactsSorted()

    def getToyResults(self):
        if not self.ATLAS_obj.toysThrown:
            print('ERROR: throw ATLAS toys first')
            sys.exit()
        if not self.CMS_obj.toysThrown:
            print('ERROR: throw CMS toys first')
            sys.exit()
        l_mt = []
        l_tot = []
        l_stat = []
        l_syst = []
        if self.ATLAS_obj.nToys != self.CMS_obj.nToys:
            print('ERROR: number of toys must be the same for ATLAS and CMS')
            sys.exit()
        nToys = self.ATLAS_obj.nToys
        for nToy in range(0,nToys):
            if nToy % 10 == 0:
                print('toy n.', nToy)
            toy_uncert_ATLAS = self.ATLAS_obj.getToyUncert(nToy,self.ATLAS_obj.systForToys)
            toy_uncert_CMS = self.CMS_obj.getToyUncert(nToy,self.CMS_obj.systForToys)
            tmpobj_ATLAS = self.ATLAS_obj.clone()
            tmpobj_CMS = self.CMS_obj.clone()
            tmpobj_ATLAS.setNewUncertainties(toy_uncert_ATLAS)
            tmpobj_CMS.setNewUncertainties(toy_uncert_CMS)
            tmpobj_LHC = LHC_object(tmpobj_ATLAS, tmpobj_CMS, blind=False, separateCombinations=True, merge_and_rename = False)
            l_mt.append(tmpobj_LHC.getBlueObject().results.mt)
            l_tot.append(tmpobj_LHC.getBlueObject().results.tot)
            l_stat.append(tmpobj_LHC.getBlueObject().results.stat)
            l_syst.append(tmpobj_LHC.getBlueObject().results.syst)
        print()
        
        return l_mt, l_tot, l_stat, l_syst

    def printCorrTables(self,draw=False):

        if self.separateCombinations:
            print ('WARNING: LHC correlation tables are trivial')
            sys.exit()

        if not os.path.exists(tab_dir):
            os.makedirs(tab_dir)

        usedMeas = self.usedMeas_sorted
    
        # for syst in list(self.corrMap.keys()):
        for syst in self.BLUE_obj.usedSyst:
            self.printCorrTable(usedMeas,syst,tab_dir)

        corr = self.BLUE_obj.printFullCorrTable()

        if draw:

            stops = array('d', [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
            red   = array('d')
            green = array('d')
            blue  = array('d')
            colors = [[103,0,31],
                      [178,24,43],
                      [214,96,77],
                      [244,165,130],
                      [253,219,199],
                      [247,247,247],
                      [209,229,240],
                      [146,197,222],
                      [67,147,195],
                      [33,102,172],
                      [5,48,97]]
            for color in colors:
                red.append(color[0]/255.)
                green.append(color[1]/255.)
                blue.append(color[2]/255.)
            rt.TColor.CreateGradientColorTable(11, stops, red[::-1], green[::-1], blue[::-1], 30)

            gStyle.SetPaintTextFormat("0.2f")
            h = TH2D('h','h',len(usedMeas),-.5,-.5+len(usedMeas),len(usedMeas),-.5,-.5+len(usedMeas))
            h.GetXaxis().SetTickLength(0)
            h.GetYaxis().SetTickLength(0)
            for i in range(0,len(usedMeas)):
                h.GetXaxis().SetBinLabel(i+1,measToROOT(usedMeas[i]))
                h.GetYaxis().SetBinLabel(i+1,measToROOT(usedMeas[len(usedMeas)-i-1]))
                for j in range(0,len(usedMeas)):
                    h.SetBinContent(i+1,j+1,corr[i][len(usedMeas)-j-1])
            h.GetZaxis().SetRangeUser(-1,1)
            c = TCanvas()
            c.SetLeftMargin(0.15)
            c.SetRightMargin(0.12)
            c.SetBottomMargin(0.1)
            c.SetTopMargin(0.1)
            h.Draw('COLZ TEXT')

            LHCLabel1 = TLatex()
            LHCLabel1.SetTextSize(0.05)
            LHCLabel1.SetTextColor(1)
            LHCLabel1.SetTextFont(42)
            LHCLabel1.DrawLatex(.01, 15, "#font[72]{ATLAS+CMS} Internal")
            

            c.SaveAs('{}/LHC_corr_plot.png'.format(plot_dir))
            c.SaveAs('{}/LHC_corr_plot.pdf'.format(plot_dir))

        return


    def getCovariance(self,syst,usedMeas=[]):
        if len(usedMeas) == 0:
            usedMeas = self.BLUE_obj.usedMeas
        u = np.array([self.BLUE_obj.p_uncert[meas][syst] for meas in usedMeas])
        corr = np.diag(np.ones(len(usedMeas)))
        for i,meas1 in enumerate(usedMeas):
            for j,meas2 in enumerate(usedMeas):
                corr[i][j] = self.BLUE_obj.p_matrix[syst][meas1][meas2]
        m = np.matmul(np.diag(u),np.matmul(corr,np.diag(u)))
        return m

    def printCorrTable(self,usedMeas,syst,tab_dir):

        m = self.BLUE_obj.p_matrix[syst]
        o = open('{}/LHC_corr_{}.tex'.format(tab_dir,syst),'w')

        f_start = open('templates/LHC_syst_start.tex')
        o.write(f_start.read().replace('#SYST#',syst).replace('#SYSTNAME#',snd.systNameDict[syst]))

        for i, meas in enumerate(usedMeas):
            if i==0:
                o.write('\t& {} '.format(measToTex(meas)))
            else:
                o.write('& {} '.format(measToTex(meas)))
        o.write('\\\\\n')
        for i, meas1 in enumerate(usedMeas):
            if not i%3:
                o.write('\\hline\n')
            o.write(measToTex(meas1)+' ')
            for meas2 in usedMeas:
                o.write('& {:.2f} '.format(m[meas1][meas2]))
            o.write('\\\\\n')

        f_end = open('templates/end.tex')
        o.write(f_end.read())
            
        return

    def sortUsedMeas(self):
        
        if self.separateCombinations:
            print('\nfunction sortUsedMeas is meant to be used for full combination only. Exiting...\n')
            sys.exit()
        
        usedMeas = [m for m in self.BLUE_obj.usedMeas if 'dil7' in m and not 'CMS' in m]
        usedMeas += [m for m in self.BLUE_obj.usedMeas if '7' in m and not 'CMS' in m and not m in usedMeas]
        usedMeas += [m for m in self.BLUE_obj.usedMeas if '8' in m and not 'CMS' in m]
        usedMeas += [m for m in self.BLUE_obj.usedMeas if 'CMS' in m]

        if sorted(usedMeas) != sorted(self.BLUE_obj.usedMeas):
            print('ERROR in printCorrTable: re-ordering of input measurements did not work')
            sys.exit()

        self.usedMeas_sorted = usedMeas

        return

    def getTotalUncertainty(self,meas):
        d = self.BLUE_obj.uncert[meas]
        tot = 0
        for syst in d.keys():
            tot += d[syst]**2
        return tot**.5

    def getUpDownRangeSyst(self,syst):
        if syst in self.ATLAS_only + self.CMS_only or syst == 'JES1':
            return None, None
        if not syst in self.corrMap.keys():
            return -.25, .25
        if self.corrMap[syst] == 0.5:
            return .25, .75
        if self.corrMap[syst] == 0.85:
            return .5, 1
        else:
            return None, None

    def getMassTotCorrPointSyst(self,syst,corr):
        clone = self.clone()
        clone.updateCorr(syst,corr)
        return np.array([clone.BLUE_obj.results.mt, clone.BLUE_obj.results.tot])

    def getDeltaScan(self, up, down, syst):
        diff = (abs(self.getMassTotCorrPointSyst(syst,up) - self.getMassTotCorrPointSyst(syst,down))/2*1000).round()
        return list(diff)

    def printSummaryTableLHC(self):

        unc_list = self.BLUE_obj.getUncertaintyListForTable()

        all_unc = []
        for l in unc_list:
            all_unc.extend(l)
        for syst in self.BLUE_obj.usedSyst:
            if not syst in all_unc and syst !='Stat':
                print('ERROR: uncertainty {} missing from categorised list. Please add it'.format(syst))
                sys.exit()

        o = open('{}/summary_table_LHC.tex'.format(tab_dir),'w')
        
        f_start = open('templates/summary_LHC.tex')
        o.write(f_start.read())

        for i,l in enumerate(unc_list):
            if i==0: o.write('\\hline')
            for syst in l:
                if not syst in self.BLUE_obj.usedSyst: continue
                o.write('\n')
                o.write(snd.systNameDict[syst])
                o.write(' & {} &'.format(self.corrMap[syst] if syst in self.corrMap.keys() else 0 if syst not in self.ATLAS_only + self.CMS_only else '\\NA'))
                if syst in self.ATLAS_only + self.CMS_only or syst in ['JES1','METH','Extra']:
                    o.write(' \\NA & \\NA & \\NA ')
                else:
                    down, up = self.getUpDownRangeSyst(syst)
                    o.write(' $[{:+}, {:+}]$'.format(down,up))
                    deltaM, deltaTot = self.getDeltaScan(up,down,syst)
                    o.write(' & {:.0f}'.format(deltaM) if deltaM>0 else ' & $<$1') 
                    o.write(' & {:.0f} '.format(deltaTot) if deltaTot>0 else ' & $<$1 ') 
                    
                o.write('\\\\')
            o.write(' [\\cmsTabSkip]')

        o.write('\n\\end{scotch}')

        return



    def makeSummaryPlot(self,blind=True):

        tge_stat = TGraphErrors()
        tge_tot = TGraphErrors()

        y_positions = dict()
        tot_points = len(self.usedMeas_sorted) + len(self.experiments) + 1
        offset = 0
        for i,meas in enumerate(self.usedMeas_sorted):
            mass, stat, tot = self.BLUE_obj.value[meas], self.BLUE_obj.uncert[meas]['Stat'], self.getTotalUncertainty(meas)
            if 'CMS' in meas and offset ==0:
                offset = 2
            y_pos = tot_points - i - offset
            y_positions[meas] = y_pos
            tge_tot.SetPoint(i,mass,y_pos)
            tge_tot.SetPointError(i,tot,0)
            tge_stat.SetPoint(i,mass,y_pos)
            tge_stat.SetPointError(i,stat,0)
            

        c = TCanvas()
        y_min = 168
        h = c.DrawFrame(y_min,-3,177,tot_points+2)
        h.GetYaxis().SetLabelSize(0)
        h.GetYaxis().SetTickLength(0)
        h.GetXaxis().SetTitle('m_{t} [GeV]')
        c.Update()
        tge_tot.SetMarkerStyle(8)
        tge_tot.Draw('p same')
        tge_stat.Draw('p same')

        latexLabel1 = TLatex()
        latexLabel1.SetTextSize(0.035)

        offset_tex = .2

        for meas in self.usedMeas_sorted:
            latexLabel1.DrawLatex(y_min+offset_tex/2.,y_positions[meas]-offset_tex,measToROOT(meas))

        tge_stat_comb = TGraphErrors()
        tge_tot_comb = TGraphErrors()
        
        y_ATLAS = tot_points-len(self.ATLAS_obj.usedMeas)
        tge_stat_comb.SetPoint(0,self.ATLAS_obj.results.mt,y_ATLAS)
        tge_stat_comb.SetPointError(0,self.ATLAS_obj.results.stat,0)
        tge_tot_comb.SetPoint(0,self.ATLAS_obj.results.mt,y_ATLAS)
        tge_tot_comb.SetPointError(0,self.ATLAS_obj.results.tot,0)

        y_CMS = y_ATLAS - len(self.CMS_obj.usedMeas) - offset
        tge_stat_comb.SetPoint(1,self.CMS_obj.results.mt,y_CMS)
        tge_stat_comb.SetPointError(1,self.CMS_obj.results.stat,0)
        tge_tot_comb.SetPoint(1,self.CMS_obj.results.mt,y_CMS)
        tge_tot_comb.SetPointError(1,self.CMS_obj.results.tot,0)


        tge_stat_comb_LHC = TGraphErrors()
        tge_tot_comb_LHC = TGraphErrors()

        y_LHC = -1
        if blind:
            mt = 172.5
        else:
            mt = self.BLUE_obj.results.mt
        tge_stat_comb_LHC.SetPoint(2,mt,y_LHC)
        tge_stat_comb_LHC.SetPointError(2,self.BLUE_obj.results.stat,0)
        tge_tot_comb_LHC.SetPoint(2,mt,y_LHC)
        tge_tot_comb_LHC.SetPointError(2,self.BLUE_obj.results.tot,0)


        tge_stat_comb.SetMarkerColor(rt.kRed)
        tge_stat_comb.SetLineColor(rt.kRed)
        tge_tot_comb.SetMarkerColor(rt.kRed)
        tge_tot_comb.SetLineColor(rt.kRed)
        tge_tot_comb.SetMarkerStyle(8)
        tge_stat_comb.Draw('p same')
        tge_tot_comb.Draw('p same')

        tge_stat_comb_LHC.SetMarkerColor(rt.kBlue)
        tge_stat_comb_LHC.SetLineColor(rt.kBlue)
        tge_tot_comb_LHC.SetMarkerColor(rt.kBlue)
        tge_tot_comb_LHC.SetLineColor(rt.kBlue)
        tge_tot_comb_LHC.SetMarkerStyle(8)
        tge_stat_comb_LHC.Draw('p same')
        tge_tot_comb_LHC.Draw('p same')

        latexLabel1.DrawLatex(y_min+offset_tex/2.,y_ATLAS-offset_tex,'#color[2]{ATLAS combination}')
        latexLabel1.DrawLatex(y_min+offset_tex/2.,y_CMS-offset_tex,'#color[2]{CMS combination}')
        if blind:
            latexLabel1.DrawLatex(y_min+offset_tex/2.,y_LHC-offset_tex,'#color[4]{LHC combination (dummy)}')
        else:
            latexLabel1.DrawLatex(y_min+offset_tex/2.,y_LHC-offset_tex,'#color[4]{LHC combination}')


        if not os.path.exists(plot_dir):
            os.makedirs(plot_dir)

        c.SaveAs('{}/summary_plot.png'.format(plot_dir))
        c.SaveAs('{}/summary_plot.pdf'.format(plot_dir))
        
                
        return


    def printPullWeightsComparisonTable(self):

        if not os.path.exists(tab_dir):
            os.makedirs(tab_dir)

        o = open('{}/LHC_weight_comparison.tex'.format(tab_dir),'w')
        O = open('{}/LHC_BLUE_weights_pulls.tex'.format(tab_dir),'w')

        f_start = open('templates/LHC_start.tex')
        o.write(f_start.read().replace('Correlation matrix','LHC pulls and weights, and weights comparison with ATLAS and CMS').replace('tab:corr','tab:pulls_weights_comparison').replace(': LHC combination',''))
        F_start = open('templates/LHC_start_paper.tex')
        O.write(F_start.read())

        for i, meas in enumerate(self.BLUE_obj.usedMeas):
            if i==0:
                o.write('\t& {} '.format(measToTex(meas)))
            else:
                o.write('& {} '.format(measToTex(meas)))
        o.write('\\\\\n\\hline\nLHC pulls ')
        O.write('Pull ')
        for meas in self.BLUE_obj.usedMeas:
            o.write('& {:.2f} '.format(self.BLUE_obj.results.pulls[meas]))
            O.write('& ${:+.2f}$ '.format(self.BLUE_obj.results.pulls[meas]))
        o.write('\\\\\nLHC weights ')
        O.write('\\\\\nWeight ')
        for meas in self.BLUE_obj.usedMeas:
            o.write('& {:.2f} '.format(self.BLUE_obj.results.weights[meas]))
            O.write('& ${:+.2f}$ '.format(self.BLUE_obj.results.weights[meas]))
        o.write('\\\\\n\\hline\nATLAS weights/2 ')
        for meas in self.BLUE_obj.usedMeas:
            if not meas in self.ATLAS_obj.usedMeas:
                o.write('& -- ')
            else:
                o.write('& {:.2f} '.format(self.ATLAS_obj.results.weights[meas]/2))
        o.write('\\\\\nCMS weights/2 ')
        for meas in self.BLUE_obj.usedMeas:
            if not meas in self.CMS_obj.usedMeas:
                o.write('& -- ')
            else:
                o.write('& {:.2f} '.format(self.CMS_obj.results.weights[meas]/2))
        o.write('\\\\\n')            

        f_end = open('templates/end.tex')
        o.write(f_end.read())
        O.write(' \\\\\n\end{scotch}}')

        
        return

    def printAllImpactsSorted(self):

        if not os.path.exists(tab_dir):
            os.makedirs(tab_dir)

        o = open('{}/impacts_all.tex'.format(tab_dir),'w')
        f_start = open('templates/impacts_all_start.tex')
        o.write(f_start.read())

        for k, v in sorted(list(self.BLUE_obj.results.mergedImpacts.items()), key=itemgetter(1), reverse = True):
            if k == 'Stat': continue
            #TODO: make sure it works when systematics are merged
            o.write('\n')
            a = self.ATLAS_obj.results.mergedImpacts[k] if k in self.ATLAS_obj.results.mergedImpacts.keys() else 999
            c = self.CMS_obj.results.mergedImpacts[k] if k in self.CMS_obj.results.mergedImpacts.keys() else 999
            o.write('{:>25}\t&\t{:.2f}\t&\t{:.2f}\t&\t{:.2f} \\\\'.format(snd.systNameDict[k],v,a,c).replace('0.00','\makebox[0pt][r]{$<$}0.01').replace('999.00','\\NA'))


        o.write(' [\\cmsTabSkip]\n')
        o.write('{:>25}\t&\t{:.2f}\t&\t{:.2f}\t&\t{:.2f} \\\\\n'.format('Total systematics',self.BLUE_obj.results.syst,self.ATLAS_obj.results.syst,self.CMS_obj.results.syst))
        o.write('{:>25}\t&\t{:.2f}\t&\t{:.2f}\t&\t{:.2f} \\\\'.format('Statistical',self.BLUE_obj.results.stat,self.ATLAS_obj.results.stat,self.CMS_obj.results.stat))
        o.write(' [\\cmsTabSkip]\n')
        o.write('{:>25}\t&\t{:.2f}\t&\t{:.2f}\t&\t{:.2f} \\\\\n'.format('Total',self.BLUE_obj.results.tot,self.ATLAS_obj.results.tot,self.CMS_obj.results.tot))

        o.write('\n\\end{scotch}')
        o.close()

        return
