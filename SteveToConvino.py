import sys, os

def getSystNames(lines):
    for line in lines:
        if 'Stat' in line:
            systnames = line.split()
            for i in range(0,len(systnames)):
                systnames[i] = systnames[i].replace('\'','')
            return systnames
    return ['ERROR!']

def removeUselessCharachters(name):
    if name.endswith('_'):
        name = name[:-1]
    if name.endswith('_'):
        name = removeUselessCharachters(name)
    return name

def getMeasurements(lines):
    measurements = []
    for line in lines:
        if not 'Mtop' in line: continue
        measurement = line.split('\'')[1].replace(' ','_')
        if not 'Mtop' in measurement: measurements.append(removeUselessCharachters(measurement))
    return measurements

def getMeasurementResult(lines,measurement,systnames):
    uncertainties = dict()
    for line in lines:
        if measurement.replace('_',' ') in line:
            uncert = line.replace(measurement.replace('_',' '),'').replace('Mtop','').replace('\'','').split()
            central = float(uncert[0])
            uncert = uncert[1:]
            if len(uncert) != len(systnames):
                print 'ERROR!'
                sys.exit()
            for i, systname in enumerate(systnames):
                uncertainties[systname] = float(uncert[i])
            break
    return [central, uncertainties]

def getAllResults(lines,measurements,systnames):
    all_uncertainties = dict()
    all_central = dict()
    for measurement in measurements:
        central, uncertainties = getMeasurementResult(lines,measurement,systnames)
        all_central[measurement] = central
        all_uncertainties[measurement] = uncertainties
    return all_central, all_uncertainties

def isSymmetricMatrix(matrix):
    if len(matrix[0]) != len(matrix[:][0]) : return False
    for i in range(0,len(matrix[0])):
        for j in range(0,len(matrix[0])):
            if matrix[i][j] != matrix[j][i]: return False
    return True

def getCorrelationMatrixSyst(lines,measurements,systname):
    for i, line in enumerate(lines):
        if systname in line and '1.0' in line:
            matrix = lines[i:i+len(measurements)]
            break
    if len(matrix)!= len(measurements):
        print 'ERROR!', systname
        sys.exit()
    if matrix[0].split()[-1].replace('\'','') != systname:
        print 'ERROR!', systname
        sys.exit()
    for i, m_line in enumerate(matrix):
        matrix[i] = m_line.split()
        if i==0 : matrix[i] = matrix[i][0:-1]
    
    if not isSymmetricMatrix(matrix):
        print 'ERROR! matrix not symmetric for syst:', systname
        sys.exit()
    all_corr_dict = dict()
    for i, meas_i in enumerate(measurements):
        corr_dict = dict()
        for j, meas_j in enumerate(measurements):
            corr_dict[meas_j] = float(matrix[i][j])
        all_corr_dict[meas_i] = corr_dict
        
    return all_corr_dict

def checkMatrixDict(matrix,systname,measurements):
    matrix = matrix[systname]
    for meas_i in measurements:
        if matrix[meas_i][meas_i]!=1: return False
        for meas_j in measurements:
            if matrix[meas_i][meas_j] != matrix[meas_j][meas_i]: return False
    return True

def checkAllMatrixDict(matrix,systnames,measurements):
    for systname in systnames:
        if systname == 'Stat': continue
        if not checkMatrixDict(matrix,systname,measurements):
            return False
    return True

def getFullCorrelationMatrix(lines,measurements,systnames):
    matrix_dict = dict()
    for systname in systnames:
        if systname == 'Stat': continue
        matrix_dict[systname] = getCorrelationMatrixSyst(lines,measurements,systname)

    if not checkAllMatrixDict(matrix_dict,systnames,measurements):
        print 'ERROR! matrix dictionary is unphysical'
        sys.exit()
    return matrix_dict

def writeConfig(outdir,systnames,measurements):
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
        f.write('\t{} = '.format(syst))
        for meas in measurements:
            f.write('{}_{} '.format(syst,meas))
            if meas != measurements[-1]:
                f.write('+ ')
            else:
                f.write('\n')
    f.write('\n\n[end uncertainty impacts]\n')
    f.close()
    return

def writeMeasurementFile(outdir,systnames,measurement,value,uncert):
    f = open('{}/{}.txt'.format(outdir,measurement),'w')
    f.write('\n\n[hessian]\n\n[end hessian]\n\n[correlation matrix]\n\n[end correlation matrix]\n\n[not fitted]\n\n\t')
    systnames.remove('Stat')
    systnames.append('Stat')
    for syst in systnames:
        if syst == 'Stat':
            f.write('stat ')
        else:
            f.write('{}_{} '.format(syst,measurement))
    f.write('\n\tmt_{} '.format(measurement))
    for syst in systnames:
        f.write('{} '.format(uncert[measurement][syst]))
    f.write('\n\n[end not fitted]\n\n[systematics]\n\n[end systematics]\n\n[estimates]\n\n')
    f.write('\tn_estimates = 1\n\n')
    f.write('\tname_0 = mt_{}\n'.format(measurement))
    f.write('\tvalue_0 = {}\n'.format(value[measurement]))
    f.write('\n[end estimates]\n')
    f.close()
    return

def writeAllFiles(outdir,systnames,measurements,value,uncert):
    for measurement in measurements:
        writeMeasurementFile(outdir,systnames,measurement,value,uncert)
    return

def writeCorrelations(outdir,systnames,measurements,matrix):
    f = open(outdir+'/extra_correlations.txt','w')
    for syst in systnames:
        if syst == 'Stat': continue
        for i, meas_i in enumerate(measurements):
            for j, meas_j in enumerate(measurements):
                if not j>i: continue
                f.write('{}_{} = ({}) {}_{}\n'.format(syst,meas_i,matrix[syst][meas_i][meas_j],syst,meas_j))
        f.write('\n')
    f.close()
    return

# infilename = 'allCMS_legacy_result.txt'
infilename = sys.argv[1]
lines = open(infilename,'r').read().splitlines()
systnames = getSystNames(lines)
measurements = getMeasurements(lines) 
value, uncert = getAllResults(lines,measurements,systnames)
matrix = getFullCorrelationMatrix(lines,measurements,systnames)

outdir = 'inputsConvinoCMS/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

writeConfig(outdir,systnames,measurements)
writeAllFiles(outdir,systnames,measurements,value,uncert)
writeCorrelations(outdir,systnames,measurements,matrix)
