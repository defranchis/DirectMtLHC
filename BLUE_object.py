import sys, os, copy
import numpy as np
from operator import itemgetter

np.random.seed(1)

tmp_dir = 'tmp_workdir'

# tocheck = ['MCGEN', 'CR', 'METH', 'RAD', 'AtlFastFull', 'TRIG', 'JES3', 'UE', 'JES8', 'HADR']


def isSymmetricMatrix(matrix):
    if len(matrix[0]) != len(matrix[:][0]) : return False
    for i in range(0,len(matrix[0])):
        for j in range(0,len(matrix[0])):
            if matrix[i][j] != matrix[j][i]: return False
    return True

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

def isPositiveDefinite(m):
    if not isInvertible(m):
        return False
    w,v = np.linalg.eig(m)
    return (w > 0).all()

def isInvertible(m):
    if np.isfinite(np.linalg.cond(m)):
        return True
    else:
        return False

class result_object:

    def __init__(self):
        self.mt = -999
        self.tot = -999
        self.stat = -999
        self.syst = -999
        self.impacts = dict()
        self.weights = dict()
        self.signedImpacts = dict()

class BLUE_object:

    def __init__(self,inputfile = None, excludeMeas = [], excludeSyst = [], ATLAS = False, LHC = False, signsOnImpacts = False):

        if inputfile is None:
            print 'ERROR: please provide input file'
            sys.exit()
        if ATLAS and LHC:
            print 'ERROR! either ATLAS or LHC combination'
            sys.exit()
            
        self.inputfile = inputfile
        self.ATLAS = ATLAS
        self.LHC = LHC
        self.CMS = not (self.ATLAS or self.LHC)
        self.signsOnImpacts = self.CMS or signsOnImpacts
        self.lines = open(inputfile,'r').read().splitlines()
        self.excludeMeas = excludeMeas
        self.excludeSyst = excludeSyst

        if self.ATLAS:
            self.systnames, self.measurements, self.value, self.uncert, self.matrix = self.getAllATLASInfo()
        else:
            self.systnames = self.getSystNames()
            self.measurements = self.getMeasurements()
            self.value, self.uncert = self.getAllResults()
            self.matrix = self.getFullCorrelationMatrix()

        if self.signsOnImpacts:
            self.p_matrix, self.p_uncert = self.propagateNegativeCorrelations()
            self.checkMatrices()
        else:
            self.p_matrix = self.matrix
            self.p_uncert = self.uncert

        self.usedMeas, self.usedSyst = self.getUsedMeasSyst()
        self.nMeas_orig = len(self.measurements)
        self.writeBLUEinputCMS()
        self.log = self.runCombination()
        self.results = self.readResults()
        self.toysInitialised = False
        self.toysThrown = False


    def update(self):
        if self.signsOnImpacts:
            self.p_matrix, self.p_uncert = self.propagateNegativeCorrelations()
            self.checkMatrices()
        else:
            self.p_matrix = self.matrix
            self.p_uncert = self.uncert            
        self.usedMeas, self.usedSyst = self.getUsedMeasSyst()
        self.log = self.runCombination()
        self.results = self.readResults()

    def clone(self):
        return copy.deepcopy(self)

    def getUsedMeasSyst(self):

        for excl in self.excludeMeas:
            if not excl in self.measurements:
                print 'ERROR! measurement {} (in exclude list) not found in input file'.format(excl)
                sys.exit()
        for excl in self.excludeSyst:
            if not excl in self.systnames:
                print 'ERROR! systematics {} (in exclude list) not found in input file'.format(excl)
                sys.exit()

        m = [meas for meas in self.measurements if not meas in self.excludeMeas]
        s = [syst for syst in self.systnames if not syst in self.excludeSyst]

        return m,s

    def getSystNames(self):
        for line in self.lines:
            if 'Stat' in line:
                systnames = line.split()
                for i in range(0,len(systnames)):
                    systnames[i] = systnames[i].replace('\'','')
                return systnames
        return ['ERROR!']
        sys.exit()

    def getMeasurements(self):
        measurements = []
        for line in self.lines:
            if not 'Mtop' in line: continue
            measurement = line.split('\'')[1].replace(' ','_')
            if not 'Mtop' in measurement:
                measurements.append(removeUselessCharachters(measurement))
        return measurements

    def getMeasurementResult(self,measurement):
        uncertainties = dict()
        for line in self.lines:
            if measurement.replace('_',' ') in line:
                uncert = line.replace(measurement.replace('_',' '),'').replace('Mtop','').replace('\'','').split()
                central = float(uncert[0])
                uncert = uncert[1:]
                if len(uncert) != len(self.systnames):
                    print 'ERROR! {} uncertainties provided with {} systnames'.format(len(uncert),len(self.systnames))
                    sys.exit()
                for i, systname in enumerate(self.systnames):
                    uncertainties[systname] = float(uncert[i])
                break
        return [central, uncertainties]

    def getAllResults(self):
        all_uncertainties = dict()
        all_central = dict()
        for measurement in self.measurements:
            central, uncertainties = self.getMeasurementResult(measurement)
            all_central[measurement] = central
            all_uncertainties[measurement] = uncertainties
        return all_central, all_uncertainties

    def getFullCorrelationMatrix(self):
        matrix_dict = dict()
        for systname in self.systnames:
            if systname == 'Stat': continue
            matrix_dict[systname] = self.getCorrelationMatrixSyst(systname)

        if not self.checkAllMatrixDict(matrix_dict):
            print 'ERROR! matrix dictionary is unphysical'
            sys.exit()
        return matrix_dict

    def getCorrelationMatrixSyst(self,systname):
        found = False
        for i, line in enumerate(self.lines):
            if '\'{}\''.format(systname) in line and '1.0' in line:
                matrix = self.lines[i:i+len(self.measurements)]
                found = True
                break
        if not found:
            print 'ERROR! matrix for syst:', systname, 'not found'
            sys.exit()
        if len(matrix)!= len(self.measurements):
            print 'logic ERROR', systname
            sys.exit()
        if matrix[0].split()[-1].replace('\'','') != systname:
            print 'ERROR! Wrong format in correlation matrix', systname
            sys.exit()
        for i, m_line in enumerate(matrix):
            if not '1.0' in m_line:
                print 'ERROR! Line missing in correlation matrix for syst:', systname
                sys.exit()
            matrix[i] = m_line.split()
            if i==0 : matrix[i] = matrix[i][0:-1]


        if not isSymmetricMatrix(matrix):
            print 'ERROR! matrix not symmetric for syst:', systname
            sys.exit()
        all_corr_dict = dict()
        for i, meas_i in enumerate(self.measurements):
            corr_dict = dict()
            for j, meas_j in enumerate(self.measurements):
                corr_dict[meas_j] = float(matrix[i][j])
            all_corr_dict[meas_i] = corr_dict

        return all_corr_dict


    def checkAllMatrixDict(self,matrix):
        for systname in self.systnames:
            if systname == 'Stat': continue
            if not self.checkMatrixDict(matrix,systname):
                return False
        return True

    def checkMatrixDict(self,matrix,systname):
        matrix = matrix[systname]
        for meas_i in self.measurements:
            if matrix[meas_i][meas_i]!=1: return False
            for meas_j in self.measurements:
                if matrix[meas_i][meas_j] != matrix[meas_j][meas_i]: return False
        return True

    def writeBLUEinputCMS(self,tmprun=False):

        if tmprun:
            if not os.path.exists(tmp_dir):
                os.makedirs(tmp_dir)
            f = open('{}/CMS_combination_tmp.txt'.format(tmp_dir),'w')
        else:
            outdir = 'signed_files'
            if not os.path.exists(outdir):
                os.makedirs(outdir)
            f = open('{}/{}_signed.txt'.format(outdir,(self.inputfile).split('/')[-1].replace('.txt','')),'w')
        f.write('\./combine<<!\n\'Top Mass combination\'\n-1 -1 0 internal & minuit debug level, dependency flag\n')
        f.write('1 {} {}  # of observables, measurements, error classes\n'.format(len(self.usedMeas),len(self.usedSyst)))
        f.write('\'Mtop\'     name of observable\n\n')
        for syst in self.usedSyst:
            f.write('\'{}\'\t'.format(syst))
        f.write('\n')


        for meas in self.usedMeas:
            f.write('\'{}\' \'Mtop\' {}'.format(meas,self.value[meas]))
            for syst in self.usedSyst:
                f.write(' {}'.format(self.p_uncert[meas][syst]))
            f.write('\n')

        f.write('\n')

        for i,meas1 in enumerate(self.usedMeas):
            for j,meas2 in enumerate(self.usedMeas):
                corr = 0.0
                if i==j:
                    corr = 1.0
                f.write('{} '.format(corr))
            if i == 0:
                f.write('\'Stat\'')
            f.write('\n')
        f.write('\n')


        for syst in self.usedSyst:
            if syst == 'Stat': 
                continue
            for i,meas1 in enumerate(self.usedMeas):
                for meas2 in self.usedMeas:
                    f.write('{} '.format(self.p_matrix[syst][meas1][meas2]))
                if i == 0:
                    f.write('\'{}\''.format(syst))
                f.write('\n')
            f.write('\n')
        f.write('!\n')

        return

    def propagateNegativeCorrelations(self):
        matrix = copy.deepcopy(self.matrix)
        uncert = copy.deepcopy(self.uncert)
        for syst in self.systnames:
            if syst == 'Stat': 
                continue
            for meas1 in self.measurements:
                for meas2 in self.measurements:
                    if matrix[syst][meas1][meas2] == 0:
                        continue
                    u1 = uncert[meas1][syst]
                    u2 = uncert[meas2][syst]
                    if u1 != 0 and u2 != 0:
                        matrix[syst][meas1][meas2] *= u1*u2/abs(u1*u2)
                    elif u1 != 0:
                        matrix[syst][meas1][meas2] *= u1/abs(u1)
                    elif u2 != 0:
                        matrix[syst][meas1][meas2] *= u2/abs(u2)

        for meas in self.measurements:
            for syst in self.systnames:
                if uncert[meas][syst] < 0:
                    uncert[meas][syst] *= -1.
        return matrix, uncert

    def checkMatrices(self):
        if not (self.checkFullMatrix(self.matrix, self.uncert) == self.checkFullMatrix(self.p_matrix,self.p_uncert)).all():
            print 'ERROR! something went wrong in sign propagation'
            sys.exit()
        return

    def checkFullMatrix(self,matrix,uncert):

        m_stat = self.getStatMatrix(matrix,uncert)
        m_tot = m_stat

        for syst in self.systnames:
            if syst == 'Stat':
                continue
            m_syst = self.getSystMatrix(syst,matrix[syst],uncert)
            m_tot += m_syst

        if not isPositiveDefinite(m_tot):
            print 'ERROR: full matrix is not positive definite\n'
            sys.exit()

        return m_tot

    def getSystMatrix(self,syst,matrix,uncert):
        m = np.zeros((len(self.measurements),len(self.measurements)))
        for i,meas1 in enumerate(self.measurements):
            for j,meas2 in enumerate(self.measurements):
                m[i,j] = matrix[meas1][meas2]*uncert[meas1][syst]*uncert[meas2][syst]
        return m

    def getStatMatrix(self,matrix,uncert):
        m =  np.zeros((len(self.measurements),len(self.measurements)))
        for i,meas in enumerate(self.measurements):
            m[i,i] = uncert[meas]['Stat']*uncert[meas]['Stat']
        return m

    def runCombination(self):
        self.writeBLUEinputCMS(tmprun=True)
        os.system('source {}/CMS_combination_tmp.txt > {}/log.log'.format(tmp_dir,tmp_dir))
        result = (open('{}/log.log'.format(tmp_dir),'r')).read()
        if not 'End of derived quantities.' in result.splitlines()[-1]:
            print 'ERROR: combination failed'
            sys.exit()
        return result
        
    def readResults(self):
        resobj = result_object()
        result = (self.log).splitlines()
        res_l = []
        res_syst_l = []
        w_index_l = []
        for j,r in enumerate(result):
            if r.replace(' ','').startswith('Mtop') and '+-' in r:
                res_l.append(r)
            if 'Errors:' in r:
                res_syst_l.append(result[j+1])
            if 'Relative weights of measurements' in r:
                w_index_l.append(j+3)
        res = res_l[-1]
        res_syst = res_syst_l[-1]

        res = res.replace('+-','')
        res = res.replace('[','')
        res = res.replace(']','')
        res = res.replace('|','')
        blah, mt, tot, stat, syst, blahblah = res.split()

        resobj.mt = float(mt)
        resobj.tot = float(tot)
        resobj.stat = float(stat)
        resobj.syst = float(syst)

        res_syst = res_syst.split()
        res_syst = [float(r) for r in res_syst if r!='Mtop']        
        
        for i,syst in enumerate(self.usedSyst):
            resobj.impacts[syst] = res_syst[i]

        w_index = w_index_l[-1]
        for z in range(w_index,w_index+len(self.usedMeas)):
            m, w = result[z].split()
            resobj.weights[m] = float(w)

        return resobj

    def printResults(self):
        print '\n-> combination results\n'
        if len(self.excludeMeas)>0:
            print 'excluded measurements:', self.excludeMeas
            print
        if len(self.excludeSyst)>0:
            print 'excluded systematics:', self.excludeSyst
            print
        self.simplePrint()
        print
        if len(self.results.signedImpacts.keys())>0:
            print self.results.signedImpacts
        else:
            print self.results.impacts
        print
        print self.results.weights
        print

    def simplePrint(self):
        print 'mt\t\ttot\tstat\tsyst\t[GeV]'
        print '{}\t\t{}\t{}\t{}'.format(self.results.mt,self.results.tot,self.results.stat,self.results.syst)
    
    def printOutput(self):
        print self.log

    
    def reduceCorrelations(self,red_corr,syst_scan='ALL'):
        if syst_scan != 'ALL' and not syst_scan in self.systnames:
            print 'ERROR: systematics {} (to scan) not in input file '
            sys.exit()
        for syst in self.systnames:
            if syst == 'Stat': continue
            if syst_scan != 'ALL' and syst != syst_scan: continue
            for m1 in self.measurements:
                for m2 in self.measurements:
                    if m1 == m2: continue
                    if self.matrix[syst][m1][m2] == 1:
                        self.matrix[syst][m1][m2] = round(red_corr,3)
        self.update()
        return

    def setCorrelationLHC(self,corr,syst):
        for m1 in self.measurements:
            for m2 in self.measurements:
                if m1 == m2: continue
                self.matrix[syst][m1][m2] = round(corr,3)
        self.update()
        return

    def increaseCorrelations(self,incr_corr,syst_scan='ALL'):
        if syst_scan != 'ALL' and not syst_scan in self.systnames:
            print 'ERROR: systematics {} (to scan) not in input file '
            sys.exit()
        for syst in self.systnames:
            if syst == 'Stat': continue
            if syst_scan != 'ALL' and syst != syst_scan: continue
            for m1 in self.measurements:
                for m2 in self.measurements:
                    if m1 == m2: continue
                    if self.matrix[syst][m1][m2] == 0:
                        self.matrix[syst][m1][m2] = round(incr_corr,3)
        self.update()
        return
    def addExcludeMeas(self,l):
        for meas in l:
            if not meas in self.measurements:
                print 'ERROR! measurement {} (to be added to exclude list) not found in input file'.format(meas)
                sys.exit()
            self.excludeMeas.append(meas)
        self.update()
        return

    def addExcludeSyst(self,l):
        for syst in l:
            if not syst in self.systnames:
                print 'ERROR! systematics {} (to be added to exclude list) not found in input file'.format(syst)
                sys.exit()
            self.excludeSyst.append(syst)
        self.update()
        return

    def prepareForToys(self,fn):
        f = open(fn,'r')
        toy_lines = f.read().splitlines()
        self.systForToys = toy_lines[0].split()
        for s in self.systForToys:
            if not s in self.systnames:
                print 'ERROR: systematic {} in file MCstat.txt not found in input file'.format(s)
                sys.exit()
            if not s in self.usedSyst:
                print 'logic Error: systematic {} (used for toys) excluded from combination'.format(s)
                sys.exit()

        self.MCstat_d = {}
        for i, line in enumerate(toy_lines):
            if i==0: continue
            l = line.split()
            thisMeas = l[0]
            if not thisMeas in self.measurements:
                print 'ERROR: measurement {} in file MCstat.txt not found in input file'.format(thisMeas)
                sys.exit()
            if not thisMeas in self.usedMeas:
                print 'WARNING: measurement {} (used for toys) excluded from combination'.format(thisMeas)

            MCstat_dd = {}
            for j,ll in enumerate(l):
                if ll == thisMeas: continue
                MCstat_dd[self.systForToys[j-1]] = float(ll)
            self.MCstat_d[thisMeas] = MCstat_dd

        self.toysInitialised = True
        # self.printToysInfo()
        return
        
    def printToysInfo(self):            
        for meas in self.MCstat_d.keys():
            print '\n-> ', meas, '\n'
            print 'syst\tnom\tMCstat'
            for syst in self.systForToys:
                print syst,'\t', self.uncert[meas][syst],'\t', self.MCstat_d[meas][syst]
        print'\n'
        return

    def throwToys(self,nToys):
        if not self.toysInitialised:
            print 'ERROR: toys must be initialised first\n'
            sys.exit()
        self.nToys = nToys
        if self.toysThrown:
            print 'WARNING: throwing new toys with nToys = {}\n'.format(nToys)
        else:
            print 'INFO: throwing toys with nToys = {}\n'.format(nToys)
        self.toy_d = {}
        for meas in self.MCstat_d.keys():
            toy_dd = {}
            for syst in self.systForToys:
                t = np.random.normal(self.uncert[meas][syst],self.MCstat_d[meas][syst],nToys)
                toy_dd[syst] = list(t)
            self.toy_d[meas] = toy_dd
        self.toysThrown = True
        return

    def getToyResults(self,l=[]):
        if not self.toysThrown:
            print 'ERROR: throw toys first'
            sys.exit()
        systForToys = self.systForToys
        print '**********'
        if len(l) > 0:
            systForToys = l
            for syst in systForToys:
                if not syst in self.systForToys:
                    print 'ERROR: systematics {} (requested for toys) not toy input file'.format(syst)
                    sys.exit()
            print 'toys restricted to systematics: {}'.format(systForToys)
        else:
            print 'toys for all available systematics: {}'.format(systForToys)
        print '**********'
        l_mt = []
        l_tot = []
        l_stat = []
        l_syst = []
        d_weights = {}
        d_syst = {}

        for meas in self.usedMeas:
            d_weights[meas] = []
        for syst in self.usedSyst:
            d_syst[syst] = []
        for nToy in range(0,self.nToys):
            toy_uncert = self.getToyUncert(nToy,systForToys)
            tmpobj = self.clone()
            tmpobj.setNewUncertainties(toy_uncert)
            l_mt.append(tmpobj.results.mt)
            l_tot.append(tmpobj.results.tot)
            l_stat.append(tmpobj.results.stat)
            l_syst.append(tmpobj.results.syst)
            for syst in tmpobj.results.impacts.keys():
                d_syst[syst].append(tmpobj.results.impacts[syst])
            for meas in tmpobj.results.weights.keys():
                d_weights[meas].append(tmpobj.results.weights[meas])
        return l_mt, l_tot, l_stat, l_syst, d_weights, d_syst

    def getToyUncert(self,nToy,systForToys):
        toy_uncert = copy.deepcopy(self.uncert)
        for meas in self.MCstat_d.keys():
            for syst in systForToys:
                toy_uncert[meas][syst] = self.toy_d[meas][syst][nToy]
        return toy_uncert
        
    def setNewUncertainties(self,uncert):
        self.uncert = uncert
        self.update()
        return
        
    def removeSigns(self):
        for meas in self.measurements:
            for syst in self.systnames:
                self.uncert[meas][syst] = abs(self.uncert[meas][syst])
        self.update()
        return

    def getAllATLASInfo(self):
        l = self.lines
        for i, line in enumerate(l):
            if i==0: continue
            if l[i-1] == '# Systematics': break
        systnames = line.split(',')
        start = l.index('# Measurements')
        stop = l.index('# Correlations')
        measurements = []
        all_uncertainties = dict()
        all_central = dict()
        for line in l[start+1:stop]:
            measurement = line.split(' ')[0]
            measurements.append(measurement)
            all_values = line.split(' ')[1].split(',')
            all_central[measurement] = float(all_values[0])
            uncertainties = dict()
            for i, syst in enumerate(systnames):
                uncertainties[syst]=float(all_values[i+1])
            all_uncertainties[measurement] = uncertainties

        start = l.index('# Correlations')
        tempdict = dict()
        for line in l[start+1:]:
            name = line.split(' ')[0]
            corrs = line.split(' ')[1].split(',')
            tempdict[name] = corrs

        all_corr_dict = dict()
        for s, syst in enumerate(systnames):
            if syst=='Stat': continue
            i_corr_dict = dict()
            for i, meas_i in enumerate(measurements):
                j_corr_dict = dict()
                for j, meas_j in enumerate(measurements):
                    if meas_i==meas_j:
                        j_corr_dict[meas_j] = 1.0
                    else:
                        corr = float(tempdict['{}_{}'.format(meas_i,meas_j)][s])
                        if self.signsOnImpacts:
                            corr = abs(corr)
                        j_corr_dict[meas_j] = corr
                i_corr_dict[meas_i] = j_corr_dict
            all_corr_dict[syst]=i_corr_dict

        goodMatrix = True
        for syst in systnames:
            if syst=='Stat': continue
            for i in measurements:
                if all_corr_dict[syst][i][i]!=1:
                    goodMatrix=False
                for j in measurements:
                    if all_corr_dict[syst][i][j]!=all_corr_dict[syst][j][i]:
                        goodMatrix=False

        if not goodMatrix:
            print 'ERROR! matrix dictionary is unphysical'
            sys.exit()

        return systnames, measurements,all_central,all_uncertainties,all_corr_dict

    def printImpactsSorted(self):
        for k, v in sorted(self.results.impacts.items(), key=itemgetter(1), reverse = True):
            print '{}\t{}'.format(k,round(v,2))
        print

    def deriveSignedImpact(self,syst):

        obj_up = self.clone()
        for meas in obj_up.usedMeas:
            obj_up.value[meas] += obj_up.uncert[meas][syst]
        obj_up.update()
        up = obj_up.results.mt - self.results.mt

        obj_down = self.clone()
        for meas in obj_down.usedMeas:
            obj_down.value[meas] -= obj_down.uncert[meas][syst]
        obj_down.update()
        down = obj_down.results.mt - self.results.mt

        return up, down



    def deriveSignedImpacts(self):
        for syst in self.results.impacts.keys():
            if syst == 'Stat' or self.results.impacts[syst] == 0:
                self.results.signedImpacts[syst] = self.results.impacts[syst]
                continue
            up, down = self.deriveSignedImpact(syst)
            # if syst in tocheck:
            #     print syst, round(up,3), round(down,3), round(self.results.impacts[syst],3)
            if up*down > 0:
                print '\nWARNING: sign of impact of {} not defined\n'.format(syst)
            elif up*down == 0:
                self.results.signedImpacts[syst] = 0
            else:
                self.results.signedImpacts[syst] = (up/abs(up)) * self.results.impacts[syst]
        return

