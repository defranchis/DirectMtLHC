import sys, os
from writeConvinoOutput import *

measurements = ['dil7','ljets7','allhad7','dil8','ljets8','allhad8']
fullInfoD = dict()

def readInputFile(measurement):
    infilename = 'inputs/ATLAS{}TeV.txt'.format(measurement)
    infile = open(infilename,'r')
    lines = infile.read().splitlines()
    central = lines[0].split()[-1]
    stat = lines[1].split()[-1]
    print central, stat, measurement
    lines = lines[2:]
    systD = dict()
    MCstatD = dict()
    for line in lines:
        systname, impact, MCstatImpact = line.split()
        systD[systname] = impact
        MCstatD[systname] = MCstatImpact
    fullInfoD[measurement] = [central, stat, systD, MCstatD]


def writeConfig(outdir,measurements):
    f = open(outdir+'/mt_config.txt','w')
    f.write('[global]\n\n[end global]\n\n[inputs]\n\n')
    f.write('\tnFiles = 1 \n\n')
    f.write('\tfile0 = allinputs.txt\n')
    f.write('\n[end inputs]\n\n[observables]\n\n\tmt_comb = ')
    for meas in measurements:
        f.write('mt_{} '.format(meas))
        if meas!= measurements[-1]: 
            f.write('+ ')
        else: f.write('\n\n[end observables]\n\n')
    f.write('[correlations]\n\n[end correlations]\n\n\n[uncertainty impacts]\n\n[end uncertainty impacts]\n')
    f.close()
    return


def getFullSystNames(measurements,fullInfoD):
    fullsystnames = fullInfoD[measurements[0]][2].keys()
    for meas in measurements:
        systnames = fullInfoD[meas][2].keys()
        for systname in systnames:
            if not systname in fullsystnames:
                fullsystnames.append(systname)
    fullsystnames.append('stat')
    return fullsystnames

def completeSystList(systD,MCstatD,fullsystnames):
    systnames = systD.keys()
    for syst in fullsystnames:
        if not syst in systnames:
            systD[syst] = 0
            MCstatD[syst] = 0
    return [systD, MCstatD]


def writeMeasurementFile(outdir,measurements,fullInfoD):
    f = open('{}/allinputs.txt'.format(outdir),'w')
    f.write('\n\n[hessian]\n\n[end hessian]\n\n[correlation matrix]\n\n[end correlation matrix]\n\n[not fitted]\n\n\t')
    
    systnames = getFullSystNames(measurements,fullInfoD)
    for syst in systnames:
        f.write('{} '.format(syst))

    for measurement in measurements:
        central, stat, systD, MCstatD = fullInfoD[measurement]
        systD['stat'] = stat
        systD, MCstatD = completeSystList(systD,MCstatD,systnames)
        f.write('\n\n\tmt_{} '.format(measurement))
        for syst in systnames:
            f.write('{} '.format(systD[syst]))
    f.write('\n\n[end not fitted]\n\n[systematics]\n\n[end systematics]\n\n[estimates]\n\n')
    f.write('\tn_estimates = {}\n\n'.format(len(measurements)))
    for i, meas in enumerate(measurements):
        f.write('\tname_{} = mt_{}\n'.format(i,meas))
        f.write('\tvalue_{} = {}\n'.format(i,fullInfoD[meas][0]))
    f.write('\n[end estimates]\n')
    f.close()

    return


for measurement in measurements:
    readInputFile(measurement)

outdir = 'inputsConvino_ATLASbreakdown'
if not os.path.exists(outdir):
    os.makedirs(outdir)


writeConfig(outdir,measurements)
writeMeasurementFile(outdir,measurements,fullInfoD)
