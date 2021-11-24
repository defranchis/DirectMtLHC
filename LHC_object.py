from BLUE_object import BLUE_object
import copy
import sys,os

corrMap_default = {'JES3': 0.5, 'JESFLV': 0.5, 'RAD': 0.5, 'MCGEN': 0.5, 'BKMC': 1.0, 'PDF': 1.0 , 'BTAG': 0.5, 'UE': 1.0, 'PU': 1.0, 'CR': 1.0}
mergeMap_default = {'HADR': ['HADR','LHCHAD'], 'JESFLV': ['CMSFL1','JES4','JES5','JES6'], 'RAD': ['RAD','Q','JPS']}

class LHC_object:

    def __init__(self,ATLAS_obj, CMS_obj, separateCombinations = True, corrMap = None, mergeMap = None, blind = False):
        self.ATLAS_obj = ATLAS_obj.clone()
        self.CMS_obj = CMS_obj.clone()
        self.blind = blind
        if blind:
            self.ATLAS_obj.makeBlind()
            self.CMS_obj.makeBlind()
        self.separateCombinations = separateCombinations
        if separateCombinations:
            if len(self.ATLAS_obj.results.signedImpacts.keys()) == 0:
                self.ATLAS_obj.deriveSignedImpacts()
            if len(self.CMS_obj.results.signedImpacts.keys()) == 0:
                self.CMS_obj.deriveSignedImpacts()
        # self.removeZeroImpacts()
        self.ATLAS_obj.simplePrint()
        self.CMS_obj.simplePrint()

        if corrMap is None:  self.corrMap = corrMap_default
        else: self.corrMap = corrMap
        if mergeMap is None: self.mergeMap = mergeMap_default
        else: self.mergeMap = mergeMap


        self.commonSyst = self.getCommonSyst()
        self.LHCmap = self.prepareLHCcombination()
        self.writeBLUEinputCMS(separateCombinations=self.separateCombinations)
        self.LHC_obj = BLUE_object('LHC_input.txt',LHC=True,signsOnImpacts=True,blind=self.blind)

    def removeZeroImpacts(self):
        for syst in self.ATLAS_obj.usedSyst:
            if self.ATLAS_obj.results.impacts[syst] == 0:
                self.ATLAS_obj.addExcludeSyst([syst])
        for syst in self.CMS_obj.usedSyst:
            if self.CMS_obj.results.impacts[syst] == 0:
                self.CMS_obj.addExcludeSyst([syst])

    def getCommonSyst(self):
        l = []
        for syst in self.ATLAS_obj.usedSyst:
            if syst in self.CMS_obj.usedSyst:
                l.append(syst)
        return l

    def prepareLHCcombination(self):
        ATLAS_only = list(set(self.ATLAS_obj.usedSyst) - set(self.commonSyst))
        CMS_only = list(set(self.CMS_obj.usedSyst) - set(self.commonSyst))
        map_d = dict()
        for syst in self.commonSyst:
            map_d[syst] = [syst]
        for syst in self.mergeMap.keys():
            if syst in self.commonSyst:
                print 'ERROR: systematics {} (in merge map) is common systematic'.format(syst)
                sys.exit()
            map_d[syst] = self.mergeMap[syst]
        for syst in map_d.keys():
            for orig_syst in map_d[syst]:
                if orig_syst in ATLAS_only:
                    ATLAS_only.remove(orig_syst)
                if orig_syst in CMS_only:
                    CMS_only.remove(orig_syst)
        print
        print '-> unique to ATLAS'
        print ATLAS_only
        print
        print '-> unique to CMS'
        print CMS_only
        print
        for syst in ATLAS_only:
            map_d[syst] = [syst]
        for syst in CMS_only:
            map_d[syst] = [syst]
        
        return map_d

    def writeBLUEinputCMS(self,tmprun=False,separateCombinations=True):

        if tmprun:
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)
            f = open('{}/LHC_combination_tmp.txt'.format(tmp_dir),'w')
        else:
            f = open('LHC_input.txt','w')
        f.write('\./combine<<!\n\'LHC Top Mass combination\'\n-1 -1 0 internal & minuit debug level, dependency flag\n')
        f.write('1 2 {}  # of observables, measurements, error classes\n'.format(len(self.LHCmap.keys())))
        f.write('\'Mtop\'     name of observable\n\n')

        exp = dict()
        exp['ATLAS'] = self.ATLAS_obj
        exp['CMS'] = self.CMS_obj
        syst_l = self.LHCmap.keys()
        syst_l.remove('Stat')
        syst_l = ['Stat'] + sorted(syst_l)

        for syst in syst_l:
            f.write('\'{}\'\t'.format(syst))
        f.write('\n')


        for meas in exp.keys():
            f.write('\'{}\' \'Mtop\' {}'.format(meas+' comb',exp[meas].results.mt))
            for syst in syst_l:
                if syst in exp[meas].results.impacts.keys():
                    #tochange!
                    # f.write(' {}'.format(exp[meas].results.impacts[syst]))
                    f.write(' {}'.format(exp[meas].results.signedImpacts[syst]))
                else:
                    comb = 0
                    for orig_syst in self.LHCmap[syst]:
                        if orig_syst in exp[meas].results.impacts.keys():
                            comb += exp[meas].results.impacts[orig_syst]**2
                    comb = comb**.5
                    f.write(' {}'.format(comb))
            f.write('\n')

        f.write('\n')

        used = []
        for syst in self.corrMap.keys():
            if not syst in syst_l:
                print 'ERROR: systematics {} (in correlation map) not found in LHC grid'.format(syst)
                sys.exit()
        for syst in syst_l:
            for i,meas1 in enumerate(exp.keys()):
                for j,meas2 in enumerate(exp.keys()):
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

    def getBlueObject(self):
        return self.LHC_obj

    def printResults(self):
        self.LHC_obj.printResults()

