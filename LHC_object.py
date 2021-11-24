from BLUE_object import BLUE_object
import copy
import sys,os

corrMap_default = {'JES3': 0.5, 'JESFLV': 0.5, 'RAD': 0.5, 'MCGEN': 0.5, 'BKMC': 1.0, 'PDF': 1.0 , 'BTAG': 0.5, 'UE': 1.0, 'PU': 1.0, 'CR': 1.0}
mergeMap_default = {'ATLAS':{}, 'CMS': {'RAD': ['Q','JPS']}}
renameMap_default = {'ATLAS':{} ,'CMS': {'CMSFL1':'JESFLV'}}

noSignsOnImpacts = {'ATLAS':['JESFLV', 'BKMC', 'BTAG', 'PDF'], 'CMS': []}

class LHC_object:

    def __init__(self,ATLAS_obj, CMS_obj, separateCombinations = True, corrMap = None, mergeMap = None, renameMap = None, blind = False):
        self.ATLAS_obj = ATLAS_obj.clone()
        self.CMS_obj = CMS_obj.clone()
        self.obj_d = {'ATLAS':self.ATLAS_obj, 'CMS':self.CMS_obj}
        self.experiments = self.obj_d.keys()
        self.blind = blind
        if blind:
            for obj in self.obj_d.values():
                obj.makeBlind()
        self.separateCombinations = separateCombinations
        # self.removeZeroImpacts()
        for exp in self.experiments:
            print '->', exp
            self.obj_d[exp].simplePrint()

        self.noSignsOnImpacts = noSignsOnImpacts

        if corrMap is None:  self.corrMap = corrMap_default
        else: self.corrMap = corrMap
        if mergeMap is None: self.mergeMap = mergeMap_default
        else: self.mergeMap = mergeMap
        if renameMap is None: self.renameMap = renameMap_default
        else: self.renameMap = renameMap

        self.renameAllSyst()
        self.mergeAllSyst()

        if self.separateCombinations:
            for obj in self.obj_d.values():
                if len(obj.results.signedImpacts.keys()) == 0:
                    obj.deriveSignedImpacts()

        self.commonSyst = self.getCommonSyst()
        self.LHCsyst = self.prepareLHCcombination()
        self.writeBLUEinputCMS(separateCombinations=self.separateCombinations) # to be changed (merging happens before)
        self.LHC_obj = BLUE_object('LHC_input.txt',LHC=True,blind=self.blind)

    def renameAllSyst(self):
        for exp in self.experiments:
            for old, new in self.renameMap[exp].items():
                self.obj_d[exp].renameSyst(old,new)
        return

    def mergeAllSyst(self):
        for exp in self.experiments:
            for merged, original_l in self.mergeMap[exp].items():
                signs_propagated = self.obj_d[exp].mergeSyst(merged,original_l) # function to be implemented
                if not signs_propagated:
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

        print
        print '-> unique to ATLAS'
        print ATLAS_only
        print
        print '-> unique to CMS'
        print CMS_only
        print

        return LHCsyst

    def writeSeparateCombinationsInput(self,f):

        f.write('\./combine<<!\n\'LHC Top Mass combination\'\n-1 -1 0 internal & minuit debug level, dependency flag\n')
        f.write('1 2 {}  # of observables, measurements, error classes\n'.format(len(self.LHCsyst)))
        f.write('\'Mtop\'     name of observable\n\n')

        syst_l = self.LHCsyst
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
            if not syst in syst_l:
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
            if not syst in used:
                print 'ERROR: systematic {} (in correlation map) not used'.format(syst)
                sys.exit()

        f.write('!\n')
        return
        
    def writeSimultaneousCombinationInput(self,f):
        # to implement
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

