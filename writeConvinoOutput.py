import numpy as np
import copy
        

def writeConfig(outdir,systnames,measurements,uncert,merged=[]):
    f = open(outdir+'/mt_config.txt','w')
    f.write('[global]\n\n[end global]\n\n[inputs]\n\n')
    f.write('\tnFiles = {} \n\n'.format(len(measurements)))
    for i, meas in enumerate(measurements):
        f.write('\tfile{} = {}.txt\n'.format(i,meas))
    f.write('\n[end inputs]\n\n[observables]\n\n\tmt_comb = ')
    for meas in measurements:
        f.write('mt_{} '.format(meas))
        if meas!= measurements[-1]: 
            f.write('+ ')
        else: f.write('\n\n[end observables]\n\n')
    f.write('[correlations]\n\n\t#!FILE  =  extra_correlations.txt\n\n[end correlations]\n\n[uncertainty impacts]\n\n')
    for syst in systnames:
        if syst == 'Stat': continue
        if syst in merged: continue
        f.write('\t{} = '.format(syst))
        first = True
        for meas in measurements:
            if uncert[meas][syst] == 0:
                continue
            if first:
                f.write('{}_{} '.format(syst,meas))
                first = False
            else:
                f.write('+ {}_{} '.format(syst,meas))

        f.write('\n')
    f.write('\n\n[end uncertainty impacts]\n')
    f.close()
    return

def writeConfigMerged(outdir,systnames,measurements,uncert,merged=[]):
    f = open(outdir+'/mt_config.txt','w')
    f.write('[global]\n\n[end global]\n\n[inputs]\n\n')
    f.write('\tnFiles = 1 \n')
    f.write('\tfile0 = all_inputs.txt\n')
    f.write('\n[end inputs]\n\n[observables]\n\n\tmt_comb = ')
    for meas in measurements:
        f.write('mt_{} '.format(meas))
        if meas!= measurements[-1]: 
            f.write('+ ')
        else: f.write('\n\n[end observables]\n\n')
    f.write('[correlations]\n\n\t#!FILE  =  extra_correlations.txt\n\n[end correlations]\n\n[uncertainty impacts]\n\n')
    for syst in systnames:
        if syst == 'Stat': continue
        if syst in merged: continue
        f.write('\t{} = '.format(syst))
        first = True
        for meas in measurements:
            if uncert[meas][syst] == 0:
                continue
            if first:
                f.write('{}_{} '.format(syst,meas))
                first = False
            else:
                f.write('+ {}_{} '.format(syst,meas))
        f.write('\n')
    f.write('\n')
    for syst in merged:
        nonzero = False
        for meas in measurements:
            if uncert[meas][syst] != 0:
                nonzero = True
                break
        if nonzero:
            f.write('\t{} = {}\n'.format(syst,syst))
    f.write('\n')
    f.write('\n\n[end uncertainty impacts]\n')
    f.close()
    return

def writeFileMerged(outdir,systnames,measurements,value,uncert,merged):
    f = open('{}/all_inputs.txt'.format(outdir),'w')
    f.write('\n\n[hessian]\n\n[end hessian]\n\n[correlation matrix]\n\n[end correlation matrix]\n\n[not fitted]\n\n\t')
    all_syst = []
    for syst in systnames:
        if 'stat' in syst.lower():
            continue
        for meas in measurements:
            if uncert[meas][syst] == 0:
                continue
            syststr = syst
            if not syst in merged:
                syststr += '_{}'.format(meas)
            if not syststr in all_syst:
                all_syst.append(syststr)
                f.write('{} '.format(syststr))
    f.write('stat ')
    for meas in measurements:
        f.write('\n\tmt_{} '.format(meas))
        for s in all_syst:
            if s in merged:
                f.write('{} '.format(uncert[meas][s]))
            else:
                this_syst = s.split('_')[0]
                this_meas = s.replace(this_syst+'_','')
                if this_meas == meas:
                    f.write('{} '.format(uncert[meas][this_syst]))
                else:
                    f.write('0 ')
        f.write('{}'.format(uncert[meas]['Stat']))
    f.write('\n\n[end not fitted]\n\n[systematics]\n\n[end systematics]\n\n[estimates]\n\n')
    f.write('\tn_estimates = {}\n\n'.format(len(measurements)))
    for i, meas in enumerate(measurements):
        f.write('\tname_{} = mt_{}\n'.format(i,meas))
        f.write('\tvalue_{} = {}\n\n'.format(i,value[meas]))
    f.write('[end estimates]\n')
    f.close()

    return

def writeMeasurementFile(outdir,systnames,measurement,value,uncert,merged=[]):
    f = open('{}/{}.txt'.format(outdir,measurement),'w')
    f.write('\n\n[hessian]\n\n[end hessian]\n\n[correlation matrix]\n\n[end correlation matrix]\n\n[not fitted]\n\n\t')
    systnames.remove('Stat')
    systnames.append('Stat')
    for syst in systnames:
        if uncert[measurement][syst] == 0:
            continue
        if syst == 'Stat':
            f.write('stat ')
        else:
            if syst in merged:
                f.write('{} '.format(syst))
            else:
                f.write('{}_{} '.format(syst,measurement))
    f.write('\n\tmt_{} '.format(measurement))
    for syst in systnames:
        if uncert[measurement][syst] == 0:
            continue
        f.write('{} '.format(uncert[measurement][syst]))
    f.write('\n\n[end not fitted]\n\n[systematics]\n\n[end systematics]\n\n[estimates]\n\n')
    f.write('\tn_estimates = 1\n\n')
    f.write('\tname_0 = mt_{}\n'.format(measurement))
    f.write('\tvalue_0 = {}\n'.format(value[measurement]))
    f.write('\n[end estimates]\n')
    f.close()
    return

def writeAllFiles(outdir,systnames,measurements,value,uncert,merged=[]):
    for measurement in measurements:
        writeMeasurementFile(outdir,systnames,measurement,value,uncert,merged)
    return

def writeCorrelations(outdir,systnames,measurements,matrix,uncert,merged=[]):
    f = open(outdir+'/extra_correlations.txt','w')
    for syst in systnames:
        if syst == 'Stat': continue
        if syst in merged: continue
        for i, meas_i in enumerate(measurements):
            if uncert[meas_i][syst] == 0:
                continue
            for j, meas_j in enumerate(measurements):
                if not j>i: continue
                if uncert[meas_j][syst] == 0:
                    continue
                if matrix[syst][meas_i][meas_j] == 0:
                    continue
                f.write('{}_{} = ({}) {}_{}\n'.format(syst,meas_i,matrix[syst][meas_i][meas_j],syst,meas_j))
        f.write('\n')
    f.close()
    return

def propagateNegativeCorrelations(matrix,systnames,measurements,uncert):
    for syst in systnames:
        if syst == 'Stat': 
            continue
        for meas1 in measurements:
            for meas2 in measurements:
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

    for meas in measurements:
        for syst in systnames:
            if uncert[meas][syst] < 0:
                uncert[meas][syst] *= -1.
    return matrix, uncert

def mergeCorrelations(systnames,measurements,uncert,matrix):

    orig_m, orig_u = propagateNegativeCorrelations(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))

    ready_s = []
    eligible_s = []

    for syst in systnames:
        if syst == 'Stat': continue
        ready = True
        eligible = True
        for meas1 in measurements:
            if uncert[meas1][syst] == 0:
                continue
            for meas2 in measurements:
                if matrix[syst][meas1][meas2] != 1:
                    if uncert[meas1][syst] != 0 and uncert[meas2][syst] != 0:
                        ready = False
                if abs(matrix[syst][meas1][meas2]) != 1:
                    if uncert[meas1][syst] != 0 and uncert[meas2][syst] != 0:
                        eligible = False
        if ready:
            ready_s.append(syst)
        elif eligible:
            eligible_s.append(syst)

    ref = measurements[0]
    for syst in eligible_s:
        signs = matrix[syst][ref]
        for meas in measurements:
            if uncert[meas][syst] == 0:
                continue            
            uncert[meas][syst] *= signs[meas]
        for m1 in measurements:
            if uncert[m1][syst] == 0:
                continue        
            for m2 in measurements:
                matrix[syst][m1][m2] = 1.0

    ready_s.extend(eligible_s)

    test_m, test_u = propagateNegativeCorrelations(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))

    if test_u != orig_u:
        print 'WARNING: something inconsistent in input uncertainties. Please check'
        for meas in measurements:
            for syst in systnames:
                if test_u[meas][syst] != orig_u[meas][syst]:
                    print meas, syst, test_u[meas][syst], orig_u[meas][syst]

    if test_m != orig_m:
        for syst in systnames:
            if syst == 'Stat': continue
            if not syst in eligible_s: continue
            if test_m[syst] != orig_m[syst]:
                print 'ERROR: correlations for systematic {} inconsistent'.format(syst)

    return ready_s, uncert, matrix


def isInvertible(m):
    if np.isfinite(np.linalg.cond(m)):
        return True
    else:
        return False

def isPositiveDefinite(m):
    if not isInvertible(m):
        return False
    w,v = np.linalg.eig(m)
    return (w > 0).all()


def getSystMatrix(syst,matrix,measurements,uncert):
    m = np.zeros((len(measurements),len(measurements)))
    for i,meas1 in enumerate(measurements):
        for j,meas2 in enumerate(measurements):
            m[i,j] = matrix[meas1][meas2]*uncert[meas1][syst]*uncert[meas2][syst]
    return m

def getStatMatrix(matrix,measurements,uncert):
    m =  np.zeros((len(measurements),len(measurements)))
    for i,meas in enumerate(measurements):
        m[i,i] = uncert[meas]['Stat']*uncert[meas]['Stat']
    return m

def checkFullMatrix(matrix,systnames,measurements,uncert):

    m_stat = getStatMatrix(matrix,measurements,uncert)
    m_tot = m_stat

    for syst in systnames:
        if syst == 'Stat':
            continue
        m_syst = getSystMatrix(syst,matrix[syst],measurements,uncert)
        m_tot += m_syst

    if isInvertible(m_tot):
        print '\ngood! full matrix is invertible'
    else:
        print '\nERROR! full matrix is not invertible\n'
        return
    if isPositiveDefinite(m_tot):
        print 'good! full matrix is positive definite\n'
    else:
        print 'ERROR: full matrix is not positive definite\n'

    return m_tot


def checkExternalCorrelations(outdir):
    f = open(outdir+'/extra_correlations.txt','r')
    lines = f.read().splitlines()

    paras = []
    
    for l in lines:
        if not '(' in l: continue
        corr = float(l.split('(')[-1].split(')')[0])
        if corr == 0: continue
        s1 = l.split(' ')[0]
        s2 = l.split(' ')[-1]
        if not s1 in paras:
            paras.append(s1)
        if not s2 in paras:
            paras.append(s2)
            
    m = np.zeros((len(paras),len(paras)))    

    for i, p in enumerate(paras):
        m[i,i] = 1

    for l in lines:
        if not '(' in l: continue
        corr = float(l.split('(')[-1].split(')')[0])
        if corr == 0: continue
        i = paras.index(l.split(' ')[0])
        j = paras.index(l.split(' ')[-1])
        if abs(corr) == 1:
            corr *= .999
        m[i,j] = corr
        m[j,i] = corr
        
    
    print '\n-> final checks on parameter correlation matrix\n'
    print 'determinant =', np.linalg.det(m)
    print 'is matrix invertible:', isInvertible(m)
    print 'is matrix positive definite:', isPositiveDefinite(m)
    print 

    if not isInvertible(m) or not isPositiveDefinite(m):
        print 'matrix is not healthy, performing further checks...\n'
        return detailedMatrixCheck(m,paras)

    return


def checkMatrixSingleSource(m,source,paras):
    for i,par in enumerate(paras):
        if not source in par:
            for j, dummy in enumerate(paras):
                if j!=i:
                    m[i,j] = 0
                    m[j,i] = 0
    det = np.linalg.det(m)
    if det <= 0:
        print '{}\t{}\t{}\t{}'.format(source, isInvertible(m), isPositiveDefinite(m), det)

    return isPositiveDefinite(m)


def checkSourceWithoutMeasurement(m,source,meas,paras):
    for i,par in enumerate(paras):
        if not source in par:
            for j, dummy in enumerate(paras):
                if j!=i:
                    m[i,j] = 0
                    m[j,i] = 0
        else:
            if meas in par:
                for j, dummy in enumerate(paras):
                    if j!=i:
                        m[i,j] = 0
                        m[j,i] = 0

    det = np.linalg.det(m)
    if det > 0:
        print '{}\t{}\t{}\t{}\t{}'.format(source, meas, isInvertible(m), isPositiveDefinite(m), det)

    return isPositiveDefinite(m)

    return

def detailedMatrixCheck(m,paras):
    sources = []
    measurements = []
    for par in paras:
        source = par.split('_')[0]
        meas = par.replace(source+'_','')
        if not source in sources:
            sources.append(source)
        if not meas in measurements:
            measurements.append(meas)
        
    print '\n-> problematic sources:'
    print '\nsource\tinvertible\tpos. def\tdet\n'
    
    failed = []
    for source in sources:
        success = checkMatrixSingleSource(copy.deepcopy(m),source,paras)
        if not success:
            failed.append(source)
    print

    print '\n-> trying excluding some inputs:'
    print '\nsource\texcluded\tinvertible\tpos. def\tdet\n'
    
    for source in failed:
        for meas in measurements:
            checkSourceWithoutMeasurement(copy.deepcopy(m),source,meas,paras)
    print
    return failed

def getSingleCorrelation(matrix,syst,measurements,trick=True):
    systmatrix = matrix[syst]
    m = np.zeros((len(measurements),len(measurements)))
    for i, m1 in enumerate(measurements):
        for j,m2 in enumerate(measurements):
            m[i,j] = systmatrix[m1][m2]
            if i!=j and abs(m[i,j]) == 1 and trick:
                m[i,j] *= .999
    return m

def getSingleCovariance(matrix,syst,measurements,uncert):
    systmatrix = matrix[syst]
    m = np.zeros((len(measurements),len(measurements)))
    for i, m1 in enumerate(measurements):
        for j,m2 in enumerate(measurements):
            m[i,j] = systmatrix[m1][m2]*uncert[m1][syst]*uncert[m2][syst]
            if i!=j and abs(systmatrix[m1][m2]) == 1:
                m[i,j] *= .999
    return m

def printSingleCorrelation(matrix,syst,measurements,trick=True):
    print '->',syst
    print measurements
    print
    m = getSingleCorrelation(matrix,syst,measurements,trick)
    print 'det = ',np.linalg.det(m)
    print m
    print


def printCombinedCorrelation(matrix,sources,measurements,uncert):
    m = np.zeros((len(measurements),len(measurements)))
    for source in sources:
        m += getSingleCovariance(matrix,source,measurements,uncert)
    for i in range(0,len(measurements)):
        for j in range(0,len(measurements)):
            if i==j: continue
            m[i,j] = m[i,j]/((m[i,i]*m[j,j])**.5)
    for i in range(0,len(measurements)):
        m[i,i] = 1
    print 'det = ', np.linalg.det(m)
    print isPositiveDefinite(m)
    print m
    return

def printCombinedCorrelations(matrix,tomerge,measurements,uncert):
    for key in tomerge.keys():
        print key
        printCombinedCorrelation(matrix,tomerge[key],measurements,uncert)
        print
    return
