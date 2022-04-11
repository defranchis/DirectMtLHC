from BLUE_object import BLUE_object
import copy
import sys,os
import numpy as np

# corrMap_default = {'JES3': 0.5, 'JESFLV': 0.5, 'RAD': 0.5, 'MCGEN': 0.5, 'BKMC': 1.0, 'PDF': 1.0 , 'BTAG': 0.5, 'UE': 1.0, 'PU': 1.0, 'CR': 1.0}
corrMap_default = {'JES3': 0.5, 'JESFLV': 0.5, 'RAD': 0.5, 'MCGEN': 0.5, 'BKMC': .85, 'PDF': .85 , 'BTAG': 0.5, 'UE': .85, 'CR': .85}
mergeMap_default = {'ATLAS':{}, 'CMS': {'RAD': ['Q','JPS'],'HADR':['SLEPB','BFRAG']}}
renameMap_default = {'ATLAS':{} ,'CMS': {'CMSFL1':'JESFLV'}}

noSignsOnImpacts = {'ATLAS':['JESFLV', 'BKMC', 'BTAG', 'PDF'], 'CMS': []}

class LHC_object:

    def __init__(self,ATLAS_obj, CMS_obj, excludeSyst = [], separateCombinations = False, corrMap = None, mergeMap = None, renameMap = None, blind = False, merge_and_rename=True,mergeImpacts={}):
        self.ATLAS_obj = ATLAS_obj.clone()
        self.CMS_obj = CMS_obj.clone()
        self.obj_d = {'ATLAS':self.ATLAS_obj, 'CMS':self.CMS_obj}

        self.excludeSyst = copy.deepcopy(excludeSyst)
        if len(self.excludeSyst) > 0:
            for obj in self.obj_d.values():
                obj.addExcludeSyst(self.excludeSyst)
            
        self.experiments = self.obj_d.keys()
        self.blind = blind
        if blind:
            for obj in self.obj_d.values():
                obj.makeBlind()
            
        self.separateCombinations = separateCombinations
        self.merge_and_rename = merge_and_rename
        self.mergeImpacts = mergeImpacts

        # self.removeZeroImpacts()
        # for exp in self.experiments:
        #     print '->', exp
        #     self.obj_d[exp].simplePrint()

        self.noSignsOnImpacts = copy.deepcopy(noSignsOnImpacts)

        if corrMap is None:  self.corrMap = copy.deepcopy(corrMap_default)
        else: self.corrMap = copy.deepcopy(corrMap)
        if mergeMap is None: self.mergeMap = copy.deepcopy(mergeMap_default)
        else: self.mergeMap = copy.deepcopy(mergeMap)
        if renameMap is None: self.renameMap = copy.deepcopy(renameMap_default)
        else: self.renameMap = copy.deepcopy(renameMap)

        if self.merge_and_rename:
            self.renameAllSyst()
            self.mergeAllSyst()

        self.commonSyst = self.getCommonSyst()
        self.LHCsyst = self.prepareLHCcombination()

        if self.separateCombinations:
            for obj in self.obj_d.values():
                if len(obj.results.signedImpacts.keys()) == 0:
                    obj.deriveSignedImpacts()
        else:
            self.LHCmeas = self.getLHCmeas()
            self.LHCmatrix = self.createLHCmatrix()

        for exp, obj in self.obj_d.items():
            for syst in self.noSignsOnImpacts[exp]:
                if syst in obj.uncertWithSign:
                    print 'ERROR: something wrong in signs of syst {} in {}'.format(syst,exp)
                    sys.exit()
            for syst in obj.uncertWithSign:
                if syst in self.noSignsOnImpacts[exp]:
                    print 'ERROR: something wrong in signs of syst {} in {}'.format(syst,exp)
                    sys.exit()

        self.writeBLUEinputCMS(separateCombinations=self.separateCombinations)
        self.LHC_obj = BLUE_object('LHC_input.txt',LHC=True,blind=self.blind,mergeImpacts=self.mergeImpacts)
        # self.LHC_obj.checkAllSystMatrices()

    def clone(self):
        return copy.deepcopy(self)

    def update(self):
        self.writeBLUEinputCMS(separateCombinations=self.separateCombinations)
        self.LHC_obj = BLUE_object('LHC_input.txt',LHC=True,blind=self.blind,mergeImpacts=self.mergeImpacts)


    def renameAllSyst(self):
        print
        for exp in self.experiments:
            for old, new in self.renameMap[exp].items():
                print '{}: renaming syst {} to {}'.format(exp,old,new)
                self.obj_d[exp].renameSyst(old,new)
        return

    def mergeAllSyst(self):
        print
        for exp in self.experiments:
            for merged, original_l in self.mergeMap[exp].items():
                print '{}: creating new syst {} from sources {}'.format(exp,merged,original_l)
                signs_propagated = self.obj_d[exp].mergeSyst(merged,original_l)
                self.noSignsOnImpacts[exp].append(merged)
        return


    def removeZeroImpacts(self):
        for obj in self.obj_d.values():
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
            print
            print '-> unique to ATLAS'
            print ATLAS_only
            print
            print '-> unique to CMS'
            print CMS_only
            print

        return LHCsyst

    def getLHCmeas(self):
        allMeas = []
        for exp, obj in self.obj_d.items():
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
                        if syst in self.obj_d['ATLAS'].matrix.keys():
                            corr = self.obj_d['ATLAS'].matrix[syst][meas1.replace('_ATLAS','')][meas2.replace('_ATLAS','')]
                            if syst in self.noSignsOnImpacts['ATLAS']:
                                pass
                            elif syst in self.ATLAS_obj.uncertWithSign:
                                corr = abs(corr)

                        else:
                            corr = 0.
                    elif '_CMS' in meas1 and '_CMS' in meas2:
                        if syst in self.obj_d['CMS'].matrix.keys():
                            corr = self.obj_d['CMS'].matrix[syst][meas1.replace('_CMS','')][meas2.replace('_CMS','')]
                            if not syst in self.noSignsOnImpacts['CMS']:
                                corr = abs(corr)
                        else:
                            corr = 0.
                    else:
                        if syst in self.corrMap.keys():
                            corr = self.corrMap[syst]
                        else:
                            corr = 0.
                    # if abs(corr) > 0.85:
                    #     corr = 0.85 * (corr/abs(corr))
                m_dd[meas2] = corr
            m_d[meas1] = m_dd

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
                if syst in self.obj_d[exp].results.impacts.keys():
                    factor = 1.
                    if syst != 'Stat' and self.obj_d[exp].results.signedImpacts[syst] < 0:
                        factor = -1.
                    f.write(' {}'.format(factor*self.obj_d[exp].results.impacts[syst]))
                else:
                    f.write(' 0.0')
            f.write('\n')

        f.write('\n')

        used = []
        for syst in self.corrMap.keys():
            if not syst in syst_l and not syst in self.excludeSyst:
                print 'ERROR: systematics {} (in correlation map) not found in LHC grid'.format(syst)
                sys.exit()
        for syst in syst_l:
            for i,meas1 in enumerate(self.experiments):
                for j,meas2 in enumerate(self.experiments):
                    corr = 0.0
                    if syst in self.corrMap.keys():
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

        for syst in self.corrMap.keys():
            if not syst in used and not syst in self.excludeSyst:
                print 'ERROR: systematic {} (in correlation map) not used'.format(syst)
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
                if syst in obj.matrix.keys() or syst == 'Stat':
                    f.write(' {}'.format(obj.uncert[meas_exp][syst]))
                else:
                    f.write(' 0.0')
            f.write('\n')
        f.write('\n')


        for i,meas1 in enumerate(self.LHCmeas):
            for j,meas2 in enumerate(self.LHCmeas):
                corr = 0.0
                if i==j:
                    corr = 1.0
                f.write('{} '.format(corr))
            if i == 0:
                f.write('\'Stat\'')
            f.write('\n')
        f.write('\n')


        for syst in syst_l:
            if syst == 'Stat': 
                continue
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
            f = open('LHC_input.txt','w')

        if separateCombinations:
            self.writeSeparateCombinationsInput(f)
        else:
            self.writeSimultaneousCombinationInput(f)

        return

    def getBlueObject(self):
        return self.LHC_obj

    def printResults(self):
        self.LHC_obj.printResults()

    def setNewLHCcorrMap(self,corrMap):
        self.corrMap = corrMap
        if not self.separateCombinations:
            self.LHCmatrix = self.createLHCmatrix()
        self.update()

    def printResults(self):
        self.LHC_obj.printResults()

    def simplePrint(self):
        self.LHC_obj.simplePrint()

    def printImpactsSorted(self):
        self.LHC_obj.printImpactsSorted()

    def getToyResults(self):
        if not self.ATLAS_obj.toysThrown:
            print 'ERROR: throw ATLAS toys first'
            sys.exit()
        if not self.CMS_obj.toysThrown:
            print 'ERROR: throw CMS toys first'
            sys.exit()
        l_mt = []
        l_tot = []
        l_stat = []
        l_syst = []
        if self.ATLAS_obj.nToys != self.CMS_obj.nToys:
            print 'ERROR: number of toys must be the same for ATLAS and CMS'
            sys.exit()
        nToys = self.ATLAS_obj.nToys
        for nToy in range(0,nToys):
            if nToy % 10 == 0:
                print 'toy n.', nToy
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
        print

        return l_mt, l_tot, l_stat, l_syst
