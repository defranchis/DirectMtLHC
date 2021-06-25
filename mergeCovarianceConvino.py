import sys, os
import numpy as np

from writeConvinoOutput import isInvertible, isPositiveDefinite

indir = '../../convino_install/from_git/Convino'

outdir = 'inputsConvinoLHC'
if not os.path.exists(outdir):
    os.makedirs(outdir)

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
    for m in mergedlist:
        l = []
        for s in systlist:
            if s.split('_')[0] == m:
                l.append(s)
        merge[m] = l
    return merge        

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

def mergeInputsForLHC(experiment):
    systlist, cov_m, mt_comb  = readPostFitInfo(experiment)
    mergedlist = getMergedSystList(systlist)
    mergemap = makeMergedMap(mergedlist,systlist)
    cov_merged = getMergedCovariance(cov_m,mergemap,systlist,mergedlist)

    if round(cov_m.sum(),5) != round(cov_merged.sum(),5):
        print 'ERROR: something wrong in merging'
        sys.exit()

    hessian = invertMatrix(cov_merged)
    writeConvinoInputFromCovariance(experiment,hessian, mergedlist, mt_comb)
    return


def main():
    for experiment in ['ATLAS','CMS']:
        mergeInputsForLHC(experiment)
    return

if __name__ == "__main__":
    main()
