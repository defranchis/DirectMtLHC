import sys, os
import itertools
from writeConvinoOutput import *
import copy

simpleInput = True
permutations = False
merge_syst = True

def removeUselessCharachters(name,useless):
    if useless in name:
        name = name.replace(useless,'')
    if useless in name:
        name = removeUselessCharachters(name,useless)
    return name


def getBLUEoutdir(infilename):
    f = open(infilename,'r')
    l = f.read().splitlines()
    for line in l:
        if 'writeInputs' in line: break
    line = line.split('\'')[1].split('\'')[0]
    return [line.replace(line.split('/')[-1],''),line.split('/')[-1]]


def getAllInfo(infilename):
    f = open(infilename,'r')
    l = f.read().splitlines()
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
        all_central[measurement] = all_values[0]
        uncertainties = dict()
        for i, syst in enumerate(systnames):
            uncertainties[syst]=abs(float(all_values[i+1]))
        all_uncertainties[measurement] = uncertainties

    start = l.index('# Correlations')
    tempdict = dict()
    for line in l[start+1:]:
        name = line.split(' ')[0]
        corrs = line.split(' ')[1].split(',')
        tempdict[name] = corrs

    # everything is read in
    # now creating matrix in standard format


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
                    j_corr_dict[meas_j] = float(tempdict['{}_{}'.format(meas_i,meas_j)][s])
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

    return [systnames, measurements,all_central,all_uncertainties,all_corr_dict]


infilename = sys.argv[1]
if simpleInput:
    systnames, measurements, value, uncert, matrix = getAllInfo(infilename)
else:
    reldir = infilename.replace(infilename.split('/')[-1],'')
    blue_outdir, allinputs = getBLUEoutdir(infilename)
    systnames, measurements, value, uncert, matrix = getAllInfo(reldir+blue_outdir+allinputs)

outdir = 'inputsConvinoATLAS/'
if permutations:
    outdir = 'inputsConvinoATLAS_debug/'
if not os.path.exists(outdir):
    os.makedirs(outdir)

m_orig = checkFullMatrix(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))

if not permutations:
    
    if merge_syst:
        merged, uncert, matrix = mergeCorrelations(systnames,measurements,copy.deepcopy(uncert),copy.deepcopy(matrix))
        m_merge = checkFullMatrix(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))
        if not (m_orig==m_merge).all():
            print 'ERROR: something wrong in merging'
            print 'please disable merging and check consistency of input correlations'
            print
            sys.exit()
    else:
        merged = []

    writeConfigMerged(outdir,systnames,measurements,uncert,merged)
    writeFileMerged(outdir,systnames,measurements,value,uncert,merged)
    writeCorrelations(outdir,systnames,measurements,matrix,uncert,merged)
    failed = checkExternalCorrelations(outdir)
    if failed is not None:
        print '\nprinting problematic matrices...\n'
        for syst in failed:
            printSingleCorrelation(matrix,syst,measurements)
    
    # tomerge = {
    #     'JES7':['JES8pilepT','JES8pilerho','JESpilemu','JESpilenvtx'],
    #     'AtlJetUnc':['JEff','JVF'],
    #     'MET':['METRes','METScale','METCellOut','METPile'],
    # }
    # printCombinedCorrelations(matrix,tomerge,measurements,uncert)        

    # for syst in systnames:
    #     if syst.lower() == 'stat': continue
    #     printSingleCorrelation(matrix,syst,measurements,trick=False)

else:
    of = open('run_all_debug.sh','w')
    of.write('mkdir -p logs_debug\n\n')

    ol = open('log_list.txt','w')

    od_orig = outdir+'input'
    for r in range(2,len(measurements)):
        comb_list = itertools.combinations(measurements,r)
        for comb in comb_list:
            comb = list(comb)
            od = od_orig
            for c in comb:
                od += '_{}'.format(c)
            if not os.path.exists(od):
                os.makedirs(od)

            m_orig = checkFullMatrix(copy.deepcopy(matrix),systnames,comb,copy.deepcopy(uncert))
            if merge_syst:
                merged, uncert, matrix =  mergeCorrelations(systnames,comb,copy.deepcopy(uncert),copy.deepcopy(matrix))
                m_merge = checkFullMatrix(copy.deepcopy(matrix),systnames,comb,copy.deepcopy(uncert))
                if not (m_orig==m_merge).all():
                    print 'ERROR: something wrong in merging'
                    sys.exit()
            else:
                merged = []
                
            writeConfig(od,systnames,comb,uncert,merged)
            writeAllFiles(od,systnames,comb,value,uncert,merged)
            writeCorrelations(od,systnames,comb,matrix,uncert,merged)
            of.write('nohup convino --prefix ATLAS_{0} {1}/mt_config.txt --neyman &> logs_debug/log_{0}.log & \n'.format(od.replace('{}/input_'.format(outdir.replace('/','')),''),od))
            ol.write('log_{}.log\n'.format(od.replace('{}/input_'.format(outdir.replace('/','')),'')))


