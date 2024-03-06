
# script from Mark Owen

from math import sqrt
from ROOT import TGraphAsymmErrors, TBox, TLine
from ROOT import kGray

Comb = 'Comb'
ATLAS = 'ATLAS'
CMS = 'CMS'

marker_size = 0.8
marker_size_stat = 0.1
line_width = 1
line_width2 = 2

class topMeas:
    def setType(self, name):
        if 'ATLAS' in name:
            self._type = ATLAS
        elif 'CMS' in name:
            self._type = CMS
        elif 'comb' in name:
            self._type = Comb


    def setDisplayName(self, name):
        self._displayName = ''
        if 'ATLAS' in name:
            self._displayName += 'ATLAS'
        elif 'CMS' in name:
            self._displayName += 'CMS'
        elif 'comb' in name:
            self._displayName += 'LHC combined'
        
        if 'full ' in name:
            self._displayName += ''
        elif 'dil ' in name:
            self._displayName += ', dilepton'
        elif 'lj ' in name:
            self._displayName += ', lepton+jets'
        elif 'aj ' in name:
            self._displayName += ', all-jets'
        elif ' vtx ' in name:
            self._displayName += ', secondary vertex'
        elif 'Psi ' in name:
            self._displayName += ', J/\\psi'
        elif ' t ' in name:
            self._displayName += ', single top'
        elif 'other' in name:
            self._displayName += ', other'
        elif 'run 1' in name:
            self._displayName += ', combined'
        else:
            print('Could not make nice name for ', name)
        
        bits = name.split('TeV')
        if len(bits) > 1:
            energy = bits[0].split(' ')[-2]
            self._displayName +=  ' ' + energy + ' TeV'
        

    def __init__(self, name, val, stat, syst):
        self._name = name
        self._val = val
        self._stat = stat
        self._syst = syst
        self._tot = sqrt(self._stat*self._stat + self._syst*self._syst)
        self.setType(name)
        self.setDisplayName(name)

    def __init__(self, name, fromJson):
        self._name = name
        self._val = fromJson[0]
        self._stat = fromJson[1]
        self._syst = fromJson[2]
        self._tot = fromJson[3]
        self.setType(name)
        self.setDisplayName(name)

    
        
    def getVal(self):
        return self._val
    def getStatErr(self):
        return self._stat
    def getSystErr(self):
        return self._syst
    def getTotalErr(self):
        return self._tot
    def getType(self):
        return self._type
    def getDisplayName(self, noExp=False):
        if noExp:
            try:
                return ' ' + self._displayName.split(',')[1]
            except IndexError:
                return '  combined'
        return self._displayName
    def getName(self):
        return self._name
    
    def getDisplayValErr(self):
        displayString = "{:3.2f}".format(self.getVal())
        displayString += ' #pm '
        displayString += "{:1.2f}".format(self.getTotalErr())
        displayString += ' (#pm '
        displayString += "{:1.2f}".format(self.getStatErr())
        displayString += ' #pm '
        displayString += "{:1.2f}".format(self.getSystErr())
        displayString += ')'
        return displayString

class measHolder:

    def __init__(self, jsondata):
        self._measurements = []
        self._ATLASmeasurements = []
        self._CMSmeasurements = []
        self._COMBmeasurements = []
        self._current_y = 0.90
        self._yPos = {}
        self._needTitle = {}
        self._ygap = 0.034
        for i in jsondata:
            self._measurements.append(topMeas(i,jsondata[i]))
            thisMeas = self._measurements[-1]
            if thisMeas.getType()==ATLAS:
                self._ATLASmeasurements.append(thisMeas)
            elif thisMeas.getType()==CMS:
                self._CMSmeasurements.append(thisMeas)
            elif thisMeas.getType()==Comb:
                self._COMBmeasurements.append(thisMeas)

    def getYpos(self, measName):
        return self._yPos[measName]

    def makeBoxGraph(self, ymax):
        for m in self._COMBmeasurements:
            if m.getName() == "full comb":
                statbox = TBox( m.getVal() - m.getStatErr(), 0.0, m.getVal() + m.getStatErr(), ymax)
                statbox.SetFillColor(kGray+1)
                statbox.SetLineColor(kGray+1)
                totbox = TBox( m.getVal() - m.getTotalErr(), 0.0, m.getVal() + m.getTotalErr(), ymax)
                totbox.SetFillColor(kGray)
                totbox.SetLineColor(kGray)
                line = TLine( m.getVal(), 0.0, m.getVal(), ymax)
                line.SetLineStyle(3)
                line.SetLineWidth(2)
                line.SetLineColor(1)
                return line, statbox, totbox
        return None

    def makeGraph(self, whichtype, doIndent=False):
        if whichtype==ATLAS:
            g, g_stat = self._makeGraph(self._ATLASmeasurements, doIndent)
            return g, g_stat
        if whichtype==CMS:
            g, g_stat = self._makeGraph(self._CMSmeasurements, doIndent)
            return g, g_stat
        if whichtype==Comb:
            g, g_stat = self._makeGraph(self._COMBmeasurements, doIndent)
            return g, g_stat

        return None

   

    def _makeGraph(self, measurements, doIndent=False):
        g = TGraphAsymmErrors(len(measurements))
        g_stat = TGraphAsymmErrors(len(measurements))
        space = 0.0
        if doIndent:
            self._current_y = self._current_y - self._ygap
        for im, m in enumerate(measurements):
            #y = first_y - 0.0317 * im - (space-0.18)
            y = self._current_y - self._ygap
            self._yPos[m.getName()] = y
            self._needTitle[m.getName()] = (doIndent and im==0)
            g.SetPoint(im,m.getVal(),y)
            g.SetPointError(im, m.getTotalErr(), m.getTotalErr(), 0.0, 0.0)
            g_stat.SetPoint(im, m.getVal(), y)
            g_stat.SetPointError(im, m.getStatErr(), m.getStatErr(), 0.0, 0.0)
            self._current_y = y
        g.SetMarkerSize(marker_size)
        g.SetLineWidth(line_width2)
        g.SetLineStyle(1)
        g_stat.SetLineWidth(line_width)
        g_stat.SetMarkerSize(marker_size_stat)
        return g, g_stat
