import sys, os
from writeConvinoOutput import *
import copy
import numpy as np
import ROOT
from ROOT import TH1F, TCanvas, TGraph
import argparse

ROOT.gROOT.SetBatch(True)

parser = argparse.ArgumentParser(description='specify options')

parser.add_argument('-f',action='store',type=str, required=True, help='<Required> input file')
parser.add_argument('--noConvino',action='store_true', help='disable Convino output')
parser.add_argument('--noBLUE',action='store_true', help='disable BLUE output')
parser.add_argument('--noMergeSyst',action='store_true', help='do not merge fully correlated sources')
parser.add_argument('--exclude',action='store', help='provide list of measurements to be excluded. Example: --exclude \'meas 1, meas 2\'')
parser.add_argument('--nToys',action='store',type=int, help='number of toys for MC stat')
parser.add_argument('--toysIndividualSyst',action='store_true', help='also run toys for each individual (relevant) systematic')
parser.add_argument('--scanCorrAll',action='store_true', help='scan all correlations with simple assumptions')
parser.add_argument('--excludeOneByOne',action='store_true', help='exclude measurements one by one')

args = parser.parse_args()

toys_dir = 'toys_workdir'
scan_dir = 'scan_workdir'
tmp_dir = 'tmp_workdir'
np.random.seed(1)


if (not args.nToys is None) and args.noBLUE:
    print '\nERROR: cannot use nToys option with noBLUE option'
    print 'exiting...\n'
    sys.exit()

def getSystNames(lines):
    for line in lines:
        if 'Stat' in line:
            systnames = line.split()
            for i in range(0,len(systnames)):
                systnames[i] = systnames[i].replace('\'','')
            return systnames
    return ['ERROR!']

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
                print 'ERROR! {} uncertainties provided with {} systnames'.format(len(uncert),len(systnames))
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
    found = False
    for i, line in enumerate(lines):
        if '\'{}\''.format(systname) in line and '1.0' in line:
            matrix = lines[i:i+len(measurements)]
            found = True
            break
    if not found:
        print 'ERROR! matrix for syst:', systname, 'not found'
        sys.exit()
    if len(matrix)!= len(measurements):
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

def newWriteOutputSteve(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,toys=False,scan=False,tmprun=False):
    lines = copy.deepcopy(lines)
    if tmprun and (toys or scan):
        print 'ERROR: either tmp run or toys or scan'
        sys.exit()
    if toys and scan:
        print 'ERROR: either toys or scan'
        sys.exit()
    if toys:
        f = open('{}/CMS_Steve_negCorr_toy.txt'.format(toys_dir),'w')
    elif scan:
        f = open('{}/CMS_Steve_negCorr_scan.txt'.format(scan_dir),'w')
    elif tmprun:
        if not os.path.exists(tmp_dir):
            os.makedirs(tmp_dir)
        f = open('{}/CMS_Steve_negCorr_tmp.txt'.format(tmp_dir),'w')
    else:
        outdir = 'signed_files' #fromhere
        if not os.path.exists(outdir):
            os.makedirs(outdir)
        f = open('{}/{}_signed.txt'.format(outdir,(args.f).split('/')[-1].replace('.txt','')),'w')
    for i in range (0,len(lines)):
        if '# of observables' in lines[i]:
            lines[i] = lines[i].replace(' {} '.format(nMeas_orig),' {} '.format(len(measurements)))
        if (toys or scan or tmprun) and '\./combine' in lines[i]:
            lines[i] = lines[i].replace('\./combine','\../../code_Steve/combine')
        
        f.write('{}\n'.format(lines[i]))
        if 'Stat' in lines[i]: break

    f.write('\n')

    for meas in measurements:
        f.write('\'{}\' \'Mtop\' {}'.format(meas,value[meas]))
        for syst in systnames:
            f.write(' {}'.format(abs(uncert[meas][syst])))
        f.write('\n')

    f.write('\n')

    for i,meas1 in enumerate(measurements):
        for j,meas2 in enumerate(measurements):
            corr = 0.0
            if i==j:
                corr = 1.0
            f.write('{} '.format(corr))
        if i == 0:
            f.write('\'Stat\'')
        f.write('\n')
    f.write('\n')


    for syst in systnames:
        if syst == 'Stat': 
            continue
        for i,meas1 in enumerate(measurements):
            for meas2 in measurements:
                f.write('{} '.format(matrix[syst][meas1][meas2]))
            if i == 0:
                f.write('\'{}\''.format(syst))
            f.write('\n')
        f.write('\n')
    f.write('!\n')

    return

def writeOutputSteve(lines,systnames,measurements,matrix,exclude,nMeas_orig):

    f = open('CMS_Steve_negCorr.txt','w')
    for i in range (0,len(lines)):
        if '# of observables' in lines[i]:
            lines[i] = lines[i].replace(' {} '.format(nMeas_orig),' {} '.format(len(measurements)))
        f.write('{}\n'.format(lines[i]))
        if 'Stat' in lines[i]: break

    f.write('\n')
    start = False

    for i in range(0,len(lines)):
        if 'Stat' in lines[i]:
            start = True
        if start and 'Mtop' in lines[i]:
            excluded = False
            for e in exclude:
                if e.replace('_',' ') in lines[i] or e in lines[i]:
                    excluded = True
                    break
            if not excluded:
                f.write('{}\n'.format(lines[i].replace('-','')))

    f.write('\n\n')

    for i,meas1 in enumerate(measurements):
        for j,meas2 in enumerate(measurements):
            corr = 0.0
            if i==j:
                corr = 1.0
            f.write('{} '.format(corr))
        if i == 0:
            f.write('\'Stat\'')
        f.write('\n')
    f.write('\n')


    for syst in systnames:
        if syst == 'Stat': 
            continue
        for i,meas1 in enumerate(measurements):
            for meas2 in measurements:
                f.write('{} '.format(matrix[syst][meas1][meas2]))
            if i == 0:
                f.write('\'{}\''.format(syst))
            f.write('\n')
        f.write('\n')
    f.write('!\n')
    return
    
def excludeMeasurements(measurements, exclude):
    exclude = [removeUselessCharachters(e) for e in exclude]
    for e in exclude:
        if not e in measurements:
            print 'ERROR: measurement "{}" not found (to be excluded)\n'.format(e)
            sys.exit()
        measurements.remove(e)

    return measurements


def throwToys(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,MCstat_d,systForToys):

    if not os.path.exists(toys_dir):
        os.makedirs(toys_dir)

    t_mt = []
    t_tot = []
    t_stat = []
    t_syst = []

    variables = ['mt','tot','stat','syst']
    nominals = []
    toys = []
    for i, v in enumerate(variables):
        toys.append([])

    if len(systForToys)>1 or not args.toysIndividualSyst:
        toy_d = {}
        print 'throwing new toys...\n'
        for meas in MCstat_d.keys():
            toy_dd = {}
            for syst in systForToys:
                t = np.random.normal(uncert[meas][syst],MCstat_d[meas][syst],args.nToys)
                toy_dd[syst] = list(t)
            toy_d[meas] = toy_dd
    else:
        toy_d = toy_dg


    toy_syst_d = {}
    nom_syst_d = {}

    toy_weight_d = {}
    nom_weight_d = {}

    for syst in systnames:
        toy_syst_d[syst] = []

    for meas in measurements:
        toy_weight_d[meas] = []

    for i in range(-1,args.nToys):
        toy_uncert = copy.deepcopy(uncert)
        if i>=0:
            for meas in MCstat_d.keys():
                for syst in systForToys:
                    toy_uncert[meas][syst] = toy_d[meas][syst][i]

        #don't forget to propagate negative correlations again!
        p_matrix, p_uncert = propagateNegativeCorrelations(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(toy_uncert))
        newWriteOutputSteve(lines,systnames,measurements,p_matrix,exclude,nMeas_orig,value,toy_uncert,True)
        os.system('source {}/CMS_Steve_negCorr_toy.txt > {}/log.log'.format(toys_dir,toys_dir))
        result = (open('{}/log.log'.format(toys_dir),'r')).read().splitlines()
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
        res = res.split()
    
        res_syst = res_syst.split()
        res_syst = [float(r) for r in res_syst if r!='Mtop']        

        w_index = w_index_l[-1]
        

        if i<0:
            for z in range(1,5):
                nominals.append(float(res[z]))
            for z,syst in enumerate(systnames):
                nom_syst_d[syst] = res_syst[z]
            for z in range(w_index,w_index+len(measurements)):
                m, w = result[z].split()
                nom_weight_d[m] = float(w)

        else:
            for z, v in enumerate(variables):
                toys[z].append(float(res[z+1]))
            for z, syst in enumerate(systnames):
                toy_syst_d[syst].append(res_syst[z])
            for z in range(w_index,w_index+len(measurements)):
                m, w = result[z].split()
                toy_weight_d[m].append(float(w))

    if args.toysIndividualSyst:
        print systForToys
    print 'var, mean, rms, nominal\n'
    for z, v in enumerate(variables):
        if z==0:
            print v,'\t', round(np.array(toys[z]).mean(),3), '  ' , round(np.array(toys[z]).std(),3), '\t', round(nominals[z],3)
        else:
            print v,'\t', round(np.array(toys[z]).mean(),3), '\t  ', round(np.array(toys[z]).std(),3), '\t', round(nominals[z],3)
    print

    if args.toysIndividualSyst:
        return toy_d

    f = open('{}/slides.tex'.format(toys_dir),'w')
    
    for z,v in enumerate(variables):
        h = TH1F(v,v,70,np.array(toys[z]).mean()-3*np.array(toys[z]).std(),np.array(toys[z]).mean()+3*np.array(toys[z]).std())
        h.SetTitle('{}; {} [GeV]; a.u.'.format(v,v))
        for j, t in enumerate(toys[z]):
            h.Fill(t)
        c = TCanvas()
        l = ROOT.TLine(nominals[z],0,nominals[z],0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        h.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/{}.png'.format(toys_dir,v))
        c.SaveAs('{}/{}.pdf'.format(toys_dir,v))
        f.write('\\begin{{frame}}{{toys for variable {}}}\n'.format(v))
        f.write('\\centering\\includegraphics[width=.75\\textwidth]{{{}/{}.pdf}}\n'.format(toys_dir,v))
        f.write('\\end{frame}\n\n')

    if not os.path.exists('{}/syst'.format(toys_dir)):
        os.makedirs('{}/syst'.format(toys_dir))

    for syst in systnames:
        if nom_syst_d[syst] == 0: continue
        mean = np.array(toy_syst_d[syst]).mean()
        rms = np.array(toy_syst_d[syst]).std()
        h = TH1F(syst,syst,70,mean/nom_syst_d[syst]-3*rms/mean,mean/nom_syst_d[syst]+3*rms/mean)
        h_abs = TH1F(syst+'_abs',syst+'_abs',70,mean-3*rms,mean+3*rms)
        h.SetTitle('{}: normalised to nominal impact; {} / nominal; a.u.'.format(syst,syst))
        h_abs.SetTitle('{}: toys; impact of {} [GeV]; a.u.'.format(syst,syst))
        for j in range(0,args.nToys):
            h.Fill(toy_syst_d[syst][j]/nom_syst_d[syst])
            h_abs.Fill(toy_syst_d[syst][j])
        l = ROOT.TLine(1,0,1,0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/syst/{}_rel.png'.format(toys_dir,syst))
        c.SaveAs('{}/syst/{}_rel.pdf'.format(toys_dir,syst))
        l = ROOT.TLine(nom_syst_d[syst],0,nom_syst_d[syst],0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h_abs.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/syst/{}_abs.png'.format(toys_dir,syst))
        c.SaveAs('{}/syst/{}_abs.pdf'.format(toys_dir,syst))
        f.write('\\begin{{frame}}{{toys for systematics {}: absolute and relative variation}}\n'.format(syst))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/syst/{}_abs.pdf}}\n'.format(toys_dir,syst))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/syst/{}_rel.pdf}}\n'.format(toys_dir,syst))
        f.write('\\end{frame}\n\n')


    if not os.path.exists('{}/weights'.format(toys_dir)):
        os.makedirs('{}/weights'.format(toys_dir))


    for meas in measurements:
        mean = np.array(toy_weight_d[meas]).mean()
        rms = np.array(toy_weight_d[meas]).std()
        h = TH1F(meas,meas,70,mean/nom_weight_d[meas]-3*rms/mean,mean/nom_weight_d[meas]+3*rms/mean)
        h_abs = TH1F(meas+'_abs',meas+'_abs',70,mean-3*rms,mean+3*rms)
        h.SetTitle('{}: normalised to nominal weight; weight / nominal; a.u.'.format(meas))
        h_abs.SetTitle('{}: weights; weight of measurement {}; a.u.'.format(meas,meas))
        for j in range(0,args.nToys):
            h.Fill(toy_weight_d[meas][j]/nom_weight_d[meas])
            h_abs.Fill(toy_weight_d[meas][j])

        l = ROOT.TLine(1,0,1,0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/weights/weight_{}_rel.png'.format(toys_dir,meas))
        c.SaveAs('{}/weights/weight_{}_rel.pdf'.format(toys_dir,meas))
        l = ROOT.TLine(nom_weight_d[meas],0,nom_weight_d[meas],0.03)
        l.SetLineWidth(2)
        l.SetLineColor(ROOT.kGreen)
        c.Clear()
        h_abs.DrawNormalized()
        l.Draw('same')
        c.SaveAs('{}/weights/weight_{}_abs.png'.format(toys_dir,meas))
        c.SaveAs('{}/weights/weight_{}_abs.pdf'.format(toys_dir,meas))

        f.write('\\begin{{frame}}{{weight of measurement {}: absolute and relative variation}}\n'.format(meas.replace('_',' ')))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/weights/weight_{}_abs.pdf}}\n'.format(toys_dir,meas))
        f.write('\\includegraphics[width=.495\\textwidth]{{{}/weights/weight_{}_rel.pdf}}\n'.format(toys_dir,meas))
        f.write('\\end{frame}\n\n')
        
    return

def getResultReducedCorrelation(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,red_corr,syst_scan='ALL'):
    m_red = copy.deepcopy(matrix)
    for syst in systnames:
        if syst == 'Stat': continue
        if syst_scan != 'ALL' and syst != syst_scan: continue
        for m1 in measurements:
            for m2 in measurements:
                if m1 == m2: continue
                if m_red[syst][m1][m2] == 1:
                    m_red[syst][m1][m2] = round(red_corr,3)
    p_matrix, p_uncert = propagateNegativeCorrelations(copy.deepcopy(m_red),systnames,measurements,copy.deepcopy(uncert))
    newWriteOutputSteve(lines,systnames,measurements,p_matrix,exclude,nMeas_orig,value,uncert,False,True)
    os.system('source {}/CMS_Steve_negCorr_scan.txt > {}/log.log'.format(scan_dir,scan_dir))
    result = (open('{}/log.log'.format(scan_dir),'r')).read().splitlines()
    res_l = []
    for j,r in enumerate(result):
        if r.replace(' ','').startswith('Mtop') and '+-' in r:
            res_l.append(r)
    res = res_l[-1]
    res = res.replace('+-','')
    res = res.replace('[','')
    res = res.replace(']','')
    res = res.replace('|','')
    blah, mt, tot, stat, syst, blahblah = res.split()

    return float(mt), float(tot), float(stat), float(syst)


def getResultIncreasedCorrelation(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,incr_corr,syst_scan='ALL'):
    m_incr = copy.deepcopy(matrix)
    for syst in systnames:
        if syst == 'Stat': continue
        if syst_scan != 'ALL' and syst != syst_scan: continue
        for m1 in measurements:
            for m2 in measurements:
                if m_incr[syst][m1][m2] == 0:
                    m_incr[syst][m1][m2] = round(incr_corr,3)
    p_matrix, p_uncert = propagateNegativeCorrelations(copy.deepcopy(m_incr),systnames,measurements,copy.deepcopy(uncert))
    newWriteOutputSteve(lines,systnames,measurements,p_matrix,exclude,nMeas_orig,value,uncert,False,True)
    os.system('source {}/CMS_Steve_negCorr_scan.txt > {}/log.log'.format(scan_dir,scan_dir))
    result = (open('{}/log.log'.format(scan_dir),'r')).read().splitlines()
    res_l = []
    for j,r in enumerate(result):
        if r.replace(' ','').startswith('Mtop') and '+-' in r:
            res_l.append(r)
    res = res_l[-1]
    res = res.replace('+-','')
    res = res.replace('[','')
    res = res.replace(']','')
    res = res.replace('|','')
    blah, mt, tot, stat, syst, blahblah = res.split()

    return float(mt), float(tot), float(stat), float(syst)

def makeCorrelationScan(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,syst_scan):
    if not os.path.exists(scan_dir+'/syst'):
        os.makedirs(scan_dir+'/syst')
    step = 0.01
    red_corrs = list(np.arange(0.7,1+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()
    print '\n scanning systematics {} \n'.format(syst_scan)
    for i,red_corr in enumerate(red_corrs):
        mt, tot, stat, syst = getResultReducedCorrelation(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,red_corr,syst_scan)
        print red_corr, mt, tot, stat, syst
        h_mt.SetPoint(i,red_corr,mt)
        h_tot.SetPoint(i,red_corr,tot)
        h_stat.SetPoint(i,red_corr,stat)
        h_syst.SetPoint(i,red_corr,syst)
        
    h_mt.SetTitle('mt scan - {}; correlation; mt [GeV]'.format(syst_scan))
    h_tot.SetTitle('tot uncert scan - {}; correlation; tot uncert [GeV]'.format(syst_scan))
    h_stat.SetTitle('stat uncert scan - {}; correlation; stat uncert [GeV]'.format(syst_scan))
    h_syst.SetTitle('syst uncert scan - {}; correlation; syst uncert [GeV]'.format(syst_scan))

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    c = TCanvas()
    h_mt.Draw('apl')
    c.SaveAs('{}/syst/mt_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/mt_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    h_tot.Draw('apl')
    c.SaveAs('{}/syst/tot_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/tot_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    h_stat.Draw('apl')
    c.SaveAs('{}/syst/stat_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/stat_{}.png'.format(scan_dir,syst_scan))
    c.Clear()

    h_syst.Draw('apl')
    c.SaveAs('{}/syst/syst_{}.pdf'.format(scan_dir,syst_scan))
    c.SaveAs('{}/syst/syst_{}.png'.format(scan_dir,syst_scan))
    c.Clear()
    

def makeCorrelationScans(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert):

    if not os.path.exists(scan_dir):
        os.makedirs(scan_dir)

    step = 0.01
    red_corrs = list(np.arange(0.7,1+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()
    for i,red_corr in enumerate(red_corrs):
        mt, tot, stat, syst = getResultReducedCorrelation(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,red_corr)
        print red_corr, mt, tot, stat, syst
        h_mt.SetPoint(i,red_corr,mt)
        h_tot.SetPoint(i,red_corr,tot)
        h_stat.SetPoint(i,red_corr,stat)
        h_syst.SetPoint(i,red_corr,syst)
        
    h_mt.SetTitle('correlation scan (reduced); correlation; mt [GeV]')
    h_tot.SetTitle('correlation scan (reduced); correlation; tot uncert [GeV]')
    h_stat.SetTitle('correlation scan (reduced); correlation; stat uncert [GeV]')
    h_syst.SetTitle('correlation scan (reduced); correlation; syst uncert [GeV]')

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    c = TCanvas()
    h_mt.Draw('apl')
    c.SaveAs('{}/mt_red.png'.format(scan_dir))
    c.Clear()

    h_tot.Draw('apl')
    c.SaveAs('{}/tot_red.png'.format(scan_dir))
    c.Clear()

    h_stat.Draw('apl')
    c.SaveAs('{}/stat_red.png'.format(scan_dir))
    c.Clear()

    h_syst.Draw('apl')
    c.SaveAs('{}/syst_red.png'.format(scan_dir))
    c.Clear()


###########


    step = 0.01
    incr_corrs = list(np.arange(-.3,0.3+step/2,step))
    h_mt = TGraph()
    h_tot = TGraph()
    h_stat = TGraph()
    h_syst = TGraph()
    for i,incr_corr in enumerate(incr_corrs):
        mt, tot, stat, syst = getResultIncreasedCorrelation(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,incr_corr)
        print incr_corr, mt, tot, stat, syst
        h_mt.SetPoint(i,incr_corr,mt)
        h_tot.SetPoint(i,incr_corr,tot)
        h_stat.SetPoint(i,incr_corr,stat)
        h_syst.SetPoint(i,incr_corr,syst)
        
    h_mt.SetTitle('correlation scan (increased); correlation; mt [GeV]')
    h_tot.SetTitle('correlation scan (increased); correlation; tot uncert [GeV]')
    h_stat.SetTitle('correlation scan (increased); correlation; stat uncert [GeV]')
    h_syst.SetTitle('correlation scan (increased); correlation; syst uncert [GeV]')

    h_mt.SetMarkerStyle(8)
    h_tot.SetMarkerStyle(8)
    h_stat.SetMarkerStyle(8)
    h_syst.SetMarkerStyle(8)

    c = TCanvas()
    h_mt.Draw('apl')
    c.SaveAs('{}/mt_incr.png'.format(scan_dir))
    c.Clear()

    h_tot.Draw('apl')
    c.SaveAs('{}/tot_incr.png'.format(scan_dir))
    c.Clear()

    h_stat.Draw('apl')
    c.SaveAs('{}/stat_incr.png'.format(scan_dir))
    c.Clear()

    h_syst.Draw('apl')
    c.SaveAs('{}/syst_incr.png'.format(scan_dir))
    c.Clear()

    f = open('{}/slides.tex'.format(scan_dir),'w')

    for syst in systnames:
        if syst == 'Stat': continue
        makeCorrelationScan(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,syst)

        f.write('\\begin{{frame}}{{correlation scan for systematics {}}}\n'.format(syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/mt_{}.pdf}}\n'.format(scan_dir,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/tot_{}.pdf}}\\\\\n'.format(scan_dir,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/stat_{}.pdf}}\n'.format(scan_dir,syst))
        f.write('\\centering\\includegraphics[width=.49\\textwidth]{{{}/syst/syst_{}.pdf}}\n'.format(scan_dir,syst))
        f.write('\\end{frame}\n\n')


    return

def getCombinationResult(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert):
    p_matrix, p_uncert = propagateNegativeCorrelations(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))
    newWriteOutputSteve(lines,systnames,measurements,p_matrix,exclude,nMeas_orig,value,uncert,False,False,True)
    os.system('source {}/CMS_Steve_negCorr_tmp.txt > {}/log.log'.format(tmp_dir,tmp_dir))
    result = (open('{}/log.log'.format(tmp_dir),'r')).read().splitlines()
    res_l = []
    for j,r in enumerate(result):
        if r.replace(' ','').startswith('Mtop') and '+-' in r:
            res_l.append(r)
    res = res_l[-1]
    res = res.replace('+-','')
    res = res.replace('[','')
    res = res.replace(']','')
    res = res.replace('|','')
    blah, mt, tot, stat, syst, blahblah = res.split()
    return float(mt), float(tot), float(stat), float(syst)

def excludeOneByOne(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert):
    for meas in measurements:
        if meas in exclude: continue
        exclude.append(meas)
        print 'excluded:', exclude
        m_red = copy.deepcopy(measurements)
        m_red.remove(meas)
        mt, tot, stat, syst = getCombinationResult(lines,systnames,m_red,matrix,exclude,nMeas_orig,value,uncert)
        exclude.remove(meas)
        print 'mt\ttot\tstat\tsyst'
        print mt, tot, stat, syst
        print
    return

infilename = args.f
lines = open(infilename,'r').read().splitlines()
systnames = getSystNames(lines)
measurements = getMeasurements(lines) 
value, uncert = getAllResults(lines,measurements,systnames)
matrix = getFullCorrelationMatrix(lines,measurements,systnames)

nMeas_orig = len(measurements)

if not args.exclude is None:
    exclude = args.exclude.split(',')
    exclude = [removeUselessCharachters(e) for e in exclude]
    measurements = excludeMeasurements(measurements,exclude)
else:
    exclude = []

m_orig = checkFullMatrix(matrix,systnames,measurements,uncert)

if not args.noBLUE:
    p_matrix, p_uncert = propagateNegativeCorrelations(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))
    m_prop = checkFullMatrix(p_matrix,systnames,measurements,p_uncert)
    if not (m_orig==m_prop).all():
        print 'ERROR: something wrong in propagaiton of negative uncertainties'
        sys.exit()
    # writeOutputSteve(lines,systnames,measurements,p_matrix,exclude,nMeas_orig)
    newWriteOutputSteve(lines,systnames,measurements,p_matrix,exclude,nMeas_orig,value,uncert)

if not args.noConvino:
    outdir = 'inputsConvinoCMS/'
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    if not args.noMergeSyst:
        merged, uncert, matrix = mergeCorrelations(systnames,measurements,copy.deepcopy(uncert),copy.deepcopy(matrix))
        m_merge = checkFullMatrix(copy.deepcopy(matrix),systnames,measurements,copy.deepcopy(uncert))
        if not (m_orig==m_merge).all():
            print 'ERROR: something wrong in merging'
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

if args.nToys > 0:
    f = open('MCstat.txt','r')
    toy_lines = f.read().splitlines()
    systForToys = toy_lines[0].split()
    for s in systForToys:
        if not s in systnames:
            print 'ERROR: systematic {} in file MCstat.txt not found in input file'.format(s)
            sys.exit()

    MCstat_d = {}
    for i, line in enumerate(toy_lines):
        if i==0: continue
        l = line.split()
        thisMeas = l[0]
        if not thisMeas in measurements:
            print 'ERROR: measurement {} in file MCstat.txt not found in input file'.format(thisMeas)
            sys.exit()
        MCstat_dd = {}
        for j,ll in enumerate(l):
            if ll == thisMeas: continue
            MCstat_dd[systForToys[j-1]] = float(ll)
        MCstat_d[thisMeas] = MCstat_dd

    for meas in MCstat_d.keys():
        print '\n', meas, '\n'
        print 'syst\tnom\tMCstat'
        for syst in systForToys:
            print syst,'\t', uncert[meas][syst],'\t', MCstat_d[meas][syst]
    print'\n'
    
    toy_dg = throwToys(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,MCstat_d,systForToys)
    if args.toysIndividualSyst:
        for st in systForToys:
            throwToys(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert,MCstat_d,[st])


if args.scanCorrAll:
    makeCorrelationScans(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert)

if args.excludeOneByOne:
    excludeOneByOne(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert)

print 'central result: mt, tot, stat, syst'
print getCombinationResult(lines,systnames,measurements,matrix,exclude,nMeas_orig,value,uncert)
print
