systNameDict = {
    'Stat':'Statistical',
    'MCGEN':'ME generator',
    'JES1':'LHC JES 1',
    'BTAG':'b tagging',
    'HADR':'LHC HAD',
    'JESFLV':'LHC b-JES',
    'TOPPT':'Top quark \pt',
    'PDF':'PDF',
    'JES2':'LHC JES 2',
    'RAD':'LHC RAD',
    'CR':'Color reconnection',
    'JER':'Jet energy resolution',
    'BKMC':'Background (simulation)',
    'LES':'Lepton energy scale',
    'METH':'Method',
    'UE':'Underlying event',
    'PU':'Pileup',
    'JES3':'LHC JES 3',
    'JES7':'LHC JES 7', # to remove
    'BKDT':'Background (data)',
    'AtlJetUnc':'ATLAS jet reconstruction',
    'MET':'Missing transverse energy',
    'Extra':'Extra', # to check
    'AtlFastFull':'ATLAS fast simulation',
    'TRIG':'Trigger efficiency',
    'JES8':'LHC JES 8', # to remove
    'CMSJES1':'CMS JES 1',
    'BFRAG':'b quark fragmentation',
    'SLEPB':'Semileptonic B-hadron BR',
    'Q':'ME scale',
    'JPS':'ME/PS matching',
    'CMSFL1':'LHC JES 4',
    #new
    'bJES':'ATLAS-b',
    'JESflavres':'ATLAS-res',
    'JESflavcomp':'ATLAS-flavcomp',
    'JES4':'CMS-g',
    'JES5':'CMS-b',
    'JES6':'CMS-l',
    'JESlight': 'LHC l-JES',
    'JESflavresLHC':'LHC g-JES',
}

syst_exp = ['JES1', 'JES2', 'JES3', 'JESFLV', 'JESflavresLHC', 'JESlight', 'CMSJES1', 'JER', 'LES', 'BTAG', 'MET', 'PU', 'TRIG']
syst_mod = ['MCGEN', 'RAD', 'HADR', 'SLEPB', 'CR', 'UE', 'PDF', 'TOPPT']
syst_bkg = ['BKDT', 'BKMC']
syst_oth = ['AtlFastFull', 'AtlJetUnc', 'METH', 'Extra']
