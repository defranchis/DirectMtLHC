import sys, os
import numpy as np
import copy

from writeConvinoOutput import isInvertible, isPositiveDefinite

indir = '../../convino_install/from_git/Convino'

outdir = 'inputsConvinoLHC'
if not os.path.exists(outdir):
    os.makedirs(outdir)

corr_d = {'MCGEN': 0.5, 'CR': 1, 'METH': 0, 'RAD': 0.5, 'JES1': 0, 'BKMC': 1, 'JER': 0, 'TRIG': 0, 'PDF': 1, 'BTAG': 0.5, 'PU': 1, 'JES2': 0, 'JES3': 0.5, 'LES': 0, 'JES6': 0, 'JES4': 1, 'JES5': 0.5, 'UE': 1, 'BKDT': 0, 'HADR': 0, 'MET': 0}

def readPostFitInfo(experiment):
    f = open('{}/{}_result.txt'.format(indir,experiment),'r')
    l_orig = f.read().splitlines()
    f.close()

    l = l_orig[l_orig.index('[full covariance matrix]') +1 : l_orig.index('[end full covariance matrix]')]
    syst_list = []
    cov = np.zeros((len(l),len(l)))
    for i,line in enumerate(l):
        entries = line.split()
        syst_list.append(entries[0])
        for j in range(0,len(entries)-1):
            cov[i,j] = float(entries[j+1])        

    for line in l_orig:
        if 'mt_comb: ' in line:
            break
    foo, mt, up, down = line.split()

    return syst_list, cov, float(mt)

def writeConvinoInputFromCovariance(experiment, hess, systs, mt):
    f = open('{}/allinputs_{}.txt'.format(outdir,experiment),'w')
    f.write('\n\n[hessian]\n')
    for i,s1 in enumerate(systs):
        f.write('{}_{} '.format(s1,experiment))
        for j,s2 in enumerate(systs):
            if j>i: continue
            f.write('{} '.format(hess[i,j]))
        f.write('\n')
            
    f.write('\n[end hessian]\n\n[correlation matrix]\n\n[end correlation matrix]\n\n[not fitted]\n\n\t')
    f.write('\n\n[end not fitted]\n\n[systematics]\n\n[end systematics]\n\n[estimates]\n\n')
    f.write('\tn_estimates = 1\n\n')
    f.write('\tname_0 = mt_{}\n'.format(experiment))
    f.write('\tvalue_0 = {}\n\n'.format(mt))
    f.write('[end estimates]\n')
    f.close()
    return

def invertMatrix(m):
    if not isInvertible(m):
        print 'Error: matrix is not invertible'
        sys.exit()
    if not isPositiveDefinite(m):
        print 'Error: matrix is not positive definite'
        sys.exit()
    return np.linalg.inv(m)

def getMergedSystList(systlist):
    merged = []
    for syst in systlist:
        if not syst.split('_')[0] in merged:
            merged.append(syst.split('_')[0])
    return merged

def makeMergedMap(mergedlist,systlist):
    merge = {}
    actually_merged = []
    for m in mergedlist:
        l = []
        for s in systlist:
            if s.split('_')[0] == m:
                l.append(s)
        merge[m] = l
        if len(l)>1:
            actually_merged.append(m)
    return merge, actually_merged

def getMergedCovariance(cov_m,mergemap,systlist,mergedlist):
    merged_m = np.zeros((len(mergedlist),len(mergedlist)))
    for i,ms1 in enumerate(mergedlist):
        for j,ms2 in enumerate(mergedlist):
            l1 = mergemap[ms1]
            l2 = mergemap[ms2]
            i1 = [systlist.index(l) for l in l1]
            i2 = [systlist.index(l) for l in l2]
            entry = 0
            for ii1 in i1:
                for ii2 in i2:
                    entry += cov_m[ii1,ii2]
            merged_m[i,j] = entry
    return merged_m

def getCorrelationFromCovariance(cov):
    corr = np.zeros((cov.shape[0],cov.shape[0]))
    for i in range(0,cov.shape[0]):
        for j in range(0,cov.shape[0]):
            if i==j:
                corr[i,j] = 1
            else:
                corr[i,j]=cov[i,j]/((cov[i,i]*cov[j,j])**.5)
    return corr

def printLargeMergedCorrelations(corr,syst,actual):
    th = 0.1
    l_large = []
    print '\nnon-negligible correlations after merging (> {}) \n'.format(th)
    for i in range(0,corr.shape[0]-1):
        for j in range(i+1,corr.shape[0]-1):
            if corr[i,j] > th:
                if syst[i] in actual or syst[j] in actual:
                    print syst[i],syst[j],round(corr[i,j],2)
                    l_large.append([syst[i],syst[j]])
    if len(l_large) == 0:
        print 'None'
    return l_large

def printOriginalCorrelations(large_corr,corr_m,mergemap,systlist):
    print '\nprinting original correlations\n'
    for s1,s2 in large_corr:
        l1 = mergemap[s1]
        l2 = mergemap[s2]
        for ll1 in l1:
            for ll2 in l2:
                print ll1, ll2, corr_m[systlist.index(ll1),systlist.index(ll2)]
        print '\n'

    return

def normalizeInputsForConvino(cov):
    norm = copy.deepcopy(cov)
    unc = []
    for i in range(0,cov.shape[0]-1):
        unc.append(cov[i,i]**.5)
    unc.append(1) #do not normalize mt
    for i in range(0,cov.shape[0]):
        for j in range(0,cov.shape[0]):
            norm[i,j] /= unc[i]*unc[j]
    return norm

def renameMergedListCMS(mergedlist):
    l = []
    for m in mergedlist:
        if m == 'mt': 
            l.append(m)
            continue
        r = m.upper()
        r = r.replace('LHC','')
        if r == 'HAD': r = 'HADR'
        l.append(r)
    return l

def mergeInputsForLHC(experiment):

    systlist, cov_m, mt_comb  = readPostFitInfo(experiment)
    corr_m = getCorrelationFromCovariance(cov_m)
    mergedlist = getMergedSystList(systlist)
    mergemap, actually_merged = makeMergedMap(mergedlist,systlist)
    cov_merged = getMergedCovariance(cov_m,mergemap,systlist,mergedlist)

    if round(cov_m.sum(),5) != round(cov_merged.sum(),5):
        print 'ERROR: something wrong in merging'
        sys.exit()

    corr_merged = getCorrelationFromCovariance(cov_merged)
    large_corr = printLargeMergedCorrelations(corr_merged,mergedlist,actually_merged)
    if len(large_corr) > 0:
        printOriginalCorrelations(large_corr,corr_m,mergemap,systlist)

    cov_norm = normalizeInputsForConvino(cov_merged)    
    hessian = invertMatrix(cov_norm)
    if experiment == 'CMS':
        mergedlist = renameMergedListCMS(mergedlist)
    writeConvinoInputFromCovariance(experiment,hessian, mergedlist, mt_comb)

    return mergedlist

def writeConvinoConfig(experiments, mergedlist_d):
    f = open('{}/mt_config.txt'.format(outdir),'w')

    f.write('[global]\n\n[end global]\n\n[inputs]\n\n')
    f.write('\tnFiles = {} \n'.format(len(experiments)))
    for i, exp in enumerate(experiments):
        f.write('\tfile{} = allinputs_{}.txt\n'.format(i,exp))
    f.write('\n[end inputs]\n\n[observables]\n\n\tmt_comb = ')
    for exp in experiments:
        f.write('mt_{} '.format(exp))
        if exp!= experiments[-1]: 
            f.write('+ ')
        else: f.write('\n\n[end observables]\n\n')
    f.write('[correlations]\n\n\t#!FILE  =  extra_correlations.txt\n\n[end correlations]\n\n[uncertainty impacts]\n\n')

    common = getCommonSyst(mergedlist_d)
    
    for c in common:
        f.write('\t{} = '.format(c))
        first = True
        for exp in experiments:
            if first:
                f.write('{}_{} '.format(c,exp))
                first = False
            else:
                f.write('+ {}_{} '.format(c,exp))
        f.write('\n')
    f.write('\n')

    for exp in experiments:
        l = mergedlist_d[exp]
        for ll in l:
            if ll == 'mt':
                continue
            if not ll in common:
                f.write('\t{} = {}_{}\n'.format(ll,ll,exp))
        f.write('\n')
    f.write('\n[end uncertainty impacts]\n')

    f.close()
    return

def getCommonSyst(mergedlist_d):
    common = []
    experiments = mergedlist_d.keys()
    onelist = mergedlist_d[experiments[0]]
    for s in onelist:
        if s == 'mt': continue
        iscommon = True
        for i in range(1,len(experiments)):
            if not s in mergedlist_d[experiments[i]]:
                iscommon = False
                break
        if iscommon:
            common.append(s)
    return common

def writeExternalCorrelations(common):
    f = open(outdir+'/extra_correlations.txt','w')
    for c in common:
        if corr_d[c] != 0:
            f.write('{}_ATLAS = ({}) {}_CMS\n'.format(c,corr_d[c],c))
    f.close()
    return

def main():
    experiments = ['ATLAS','CMS']
    mergedlist_d = dict()
    for exp in experiments:
        print '\n\n-> processing {} inputs \n'.format(exp)
        ml = mergeInputsForLHC(exp)
        mergedlist_d[exp] = ml

    writeConvinoConfig(experiments,mergedlist_d)   
    writeExternalCorrelations(getCommonSyst(mergedlist_d))

    return

if __name__ == "__main__":
    main()
