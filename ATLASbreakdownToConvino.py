import sys, os
import itertools

allmeasurements = ['dil7','ljets7','allhad7','dil8','ljets8','allhad8']
fullInfoD = dict()

permutations = True

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

    return

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


def getFullSystNames():
    fsn = fullInfoD[allmeasurements[0]][2].keys()
    for meas in allmeasurements:
        systnames = fullInfoD[meas][2].keys()
        for systname in systnames:
            if not systname in fsn:
                fsn.append(systname)
    fsn.append('stat')
    return fsn

def completeSystList(systD,MCstatD):
    systnames = systD.keys()
    for syst in fullsystnames:
        if not syst in systnames:
            systD[syst] = 0
            MCstatD[syst] = 0
    return [systD, MCstatD]


def writeMeasurementFile(outdir,measurements):
    f = open('{}/allinputs.txt'.format(outdir),'w')
    f.write('\n\n[hessian]\n\n[end hessian]\n\n[correlation matrix]\n\n[end correlation matrix]\n\n[not fitted]\n\n\t')
    
    for syst in fullsystnames:
        f.write('{} '.format(syst))

    for measurement in measurements:
        central, stat, systD, MCstatD = fullInfoD[measurement]
        systD['stat'] = stat
        systD, MCstatD = completeSystList(systD,MCstatD)
        f.write('\n\n\tmt_{} '.format(measurement))
        for syst in fullsystnames:
            f.write('{} '.format(systD[syst]))
    f.write('\n\n[end not fitted]\n\n[systematics]\n\n[end systematics]\n\n[estimates]\n\n')
    f.write('\tn_estimates = {}\n\n'.format(len(measurements)))
    for i, meas in enumerate(measurements):
        f.write('\tname_{} = mt_{}\n'.format(i,meas))
        f.write('\tvalue_{} = {}\n'.format(i,fullInfoD[meas][0]))
    f.write('\n[end estimates]\n')
    f.close()

    return



outdir = 'inputsConvino_ATLASbreakdown/'
if permutations: 
    outdir = 'inputsConvino_ATLASbreakdown_debug/'

if not os.path.exists(outdir):
    os.makedirs(outdir)

for measurement in allmeasurements:
    readInputFile(measurement)
fullsystnames = getFullSystNames()

if not permutations:
    writeConfig(outdir,allmeasurements)
    writeMeasurementFile(outdir,allmeasurements)

else:
    of = open('run_all_debug_breakdown.sh','w')
    of.write('mkdir -p logs_debug_breakdown\n\n')

    ol = open('log_list_breakdown.txt','w')

    outdir += 'input'
    for r in range(2,len(allmeasurements)):
        comb_list = itertools.combinations(allmeasurements,r)
        for comb in comb_list:
            comb = list(comb)
            od = outdir
            for c in comb:
                od += '_{}'.format(c)
            if not os.path.exists(od):
                os.makedirs(od)
            writeConfig(od,comb)
            writeMeasurementFile(od,comb)

            of.write('nohup convino --prefix ATLAS_{0} {1}/mt_config.txt -d --neyman &> logs_debug_breakdown/log_{0}.log &\nsleep 20\n\n'.format(od.replace('inputsConvino_ATLASbreakdown_debug/input_',''),od))
            ol.write('log_{}.log\n'.format(od.replace('inputsConvino_ATLASbreakdown_debug/input_','')))

            


