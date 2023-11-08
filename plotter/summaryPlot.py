
# script from Mark Owen

from ROOT import gROOT, TCanvas, TLatex, gPad, gStyle, TH1F, TGraphAsymmErrors, TLegend
from ROOT import kRed, kBlue, kBlack, kGray
import json
import itertools
from topMeas import measHolder, ATLAS, CMS, Comb, marker_size, line_width2, line_width

gROOT.SetMacroPath(".:$(HOME)/atlasrootstyle:");
#gROOT.SetDynamicPath(".:$(ROOTSYS)/lib:$(HOME)/atlasrootstyle");

gROOT.LoadMacro("AtlasStyle.C")
gROOT.ProcessLine('SetAtlasStyle()')


gROOT.SetBatch(True)

# Define a bunch of constants
defaultFont = 42
gStyle.SetPadTickX(0)
gStyle.SetOptStat(0) 
gStyle.SetTextFont(defaultFont)  
gStyle.SetMarkerStyle(20) 
  
color_theory_scale = kGray+1
color_theory_scale_transparency = 0.5
color_theory_transparency = 0.6
color_theory = kGray

end_error_size = 6

markerstyle_ATLAS = 21
markerstyle_CMS = 20
markerstyle_comb = 23

color_ATLAS = kBlue
color_CMS = kRed
color_comb = kBlack

xmin = 162.3
xmax = 185.0

text_size = 0.026

# open input file and load measurements
try:
    infile = open('summary_LHC.json','r')
except IOError:
    print('Could not open input json')
    exit(1)
inputdata = json.load(infile)
measurements = measHolder(inputdata)

# create canvas
c = TCanvas("c_topmass", "c_topmass",850,850)
c.SetBorderSize(0)
gPad.SetLeftMargin(0.05)
gPad.SetRightMargin(0.05)
gPad.SetTopMargin(0.01)
gPad.SetBottomMargin(0.09)

# Draw axis
haxis = TH1F('haxis',';m_{t} [GeV]',1,xmin,xmax)
haxis.SetLineWidth(0)
haxis.Draw()
haxis.GetYaxis().SetNdivisions(0)
haxis.GetXaxis().CenterTitle(1)
haxis.GetXaxis().SetTitleSize(0.037)
haxis.GetXaxis().SetTitleOffset(1.1)
haxis.GetXaxis().SetLabelSize(0.037)

# Draw LHC result as shaded box
line, statbox, totbox = measurements.makeBoxGraph(0.87)
totbox.Draw()
statbox.Draw()
line.Draw()

# Draw ATLAS measurements
ag, ag_stat = measurements.makeGraph(ATLAS, True)
gStyle.SetEndErrorSize(end_error_size)
ag.SetMarkerStyle(markerstyle_ATLAS); 
ag.SetMarkerColor(color_ATLAS)
ag.SetLineColor(color_ATLAS)
ag.Draw('p')
ag_stat.SetMarkerStyle(markerstyle_ATLAS)
ag_stat.SetLineColor(color_ATLAS)
ag_stat.SetMarkerColor(color_ATLAS)
ag_stat.Draw("P")
ag.Draw("P")

# Draw CMS measurements
cg, cg_stat = measurements.makeGraph(CMS, True)
cg.SetMarkerStyle(markerstyle_CMS); 
cg.SetMarkerColor(color_CMS)
cg.SetLineColor(color_CMS)
cg.Draw('p')
cg_stat.SetMarkerStyle(markerstyle_CMS)
cg_stat.SetLineColor(color_CMS)
cg_stat.SetMarkerColor(color_CMS)
cg_stat.Draw("P")
cg.Draw("P")

# Draw combinations
cog, cog_stat = measurements.makeGraph(Comb, True)
cog.SetMarkerStyle(markerstyle_comb)
cog.SetLineColor(color_comb)
cog.SetMarkerColor(color_comb)
cog_stat.SetMarkerStyle(markerstyle_comb)
cog_stat.SetMarkerColor(color_comb)
cog_stat.Draw('P')
cog.Draw('P')


# ATLAS / CMS labels
y_labels = 1.01;
LHCLabel1 = TLatex()
LHCLabel1.SetTextSize(text_size*1.1)
LHCLabel1.SetTextColor(1)
LHCLabel1.SetTextFont(defaultFont)
LHCLabel1.DrawLatex(xmin + (xmax-xmin)*0.02, y_labels, "#font[62]{ATLAS+CMS}")
LHCLabel1.SetTextFont(defaultFont);
# LHCLabel1.SetTextSize(0.9*text_size);
# LHCLabel1.DrawLatex(xmin + (xmax-xmin)*0.02, y_labels-0.03,"LHC#font[52]{#scale[1.2]{top}}WG")
LHCLabel1.SetTextSize(text_size)
LHCLabel1.DrawLatex(xmin + (xmax-xmin)*0.8, y_labels, "#sqrt{s}=7,8 TeV")

# Legend for boxes
boxLeg = TLegend(xmin + (xmax-xmin)*0.02,  y_labels-0.11, xmin + (xmax-xmin)*0.4, y_labels-0.04,"","NB")
#boxLeg = TLegend(0.05, 0.2, 0.8, 0.9)
boxLeg.SetFillStyle(0)
boxLeg.SetLineWidth(0)
boxLeg.SetShadowColor(0)
boxLeg.SetTextSize(text_size)
boxLeg.SetTextFont(defaultFont)
boxLeg.AddEntry(line, "ATLAS+CMS combined", "l")
boxLeg.AddEntry(statbox, "stat uncertainty","f")
boxLeg.AddEntry(totbox, "total uncertainty","f")
boxLeg.Draw()


# redraw axis
gPad.RedrawAxis()

# Header for m(top) values
latexMeasLabel = TLatex()
latexMeasLabel.SetTextSize(text_size*0.8); 
latexMeasLabel.DrawLatex( xmin + (xmax-xmin)*0.7, measurements.getYpos(measurements._measurements[0].getName()) + 0.02, "m_{t} #pm total (#pm stat #pm syst) [GeV]")

# Measurement labels
for m in itertools.chain(measurements._measurements):
    latexMeasLabel.SetTextFont(defaultFont)
    latexMeasLabel.SetTextSize(text_size)
    print( m.getName(), measurements.getYpos(m.getName()))
    if measurements._needTitle[m.getName()]:
        title = m.getType()
        if title==Comb:
            title = "ATLAS+CMS"
            LHCLabel1.SetTextSize(0.85*text_size);
            LHCLabel1.DrawLatex(xmin + (xmax-xmin)*0.19,  measurements.getYpos(m.getName()) - 0.01 + measurements._ygap,"LHC#font[52]{#scale[1.2]{top}}WG")
        latexMeasLabel.SetTextFont(62)
        latexMeasLabel.DrawLatex(xmin + (xmax-xmin)*0.02, measurements.getYpos(m.getName()) - 0.01 + measurements._ygap, title)
    latexMeasLabel.SetTextFont(defaultFont)
    if m.getName() == "full comb" or 'run 1' in m.getName():
        latexMeasLabel.SetTextFont(62)
        latexMeasLabel.SetTextSize(1.05*text_size)
    latexMeasLabel.DrawLatex(xmin + (xmax-xmin)*0.02, measurements.getYpos(m.getName()) - 0.01, m.getDisplayName(True))
    latexMeasLabel.SetTextSize(text_size*0.8);
    if m.getName() == "full comb" or 'run 1' in m.getName() or m.getName().endswith('_comb'):
        latexMeasLabel.SetTextFont(62)
        latexMeasLabel.SetTextSize(latexMeasLabel.GetTextSize()*1.05); 
    latexMeasLabel.DrawLatex(xmin + (xmax-xmin)*0.7, measurements.getYpos(m.getName()) - 0.01, m.getDisplayValErr())    

# Label for stat syst bars
err_size = 2.0
err_leg_y = 0.95
err_leg = TGraphAsymmErrors(1)
err_leg.SetPoint(0, 180.0, err_leg_y)
err_leg.SetMarkerStyle(markerstyle_comb)
err_leg.SetLineColor(kGray+2)
err_leg.SetMarkerColor(kGray+2)
err_leg.SetMarkerSize(marker_size)
err_leg.SetLineWidth(line_width2)
err_stat_leg = TGraphAsymmErrors(err_leg)
err_stat_leg.SetLineWidth(line_width)
err_leg.SetPointError(0, err_size, err_size, 0.0, 0.0)
err_stat_leg.SetPointError(0, err_size/3.0, err_size/3.0, 0.0, 0.0)
err_leg.Draw("p")
err_stat_leg.Draw("p")
latexMeasLabel.SetTextFont(defaultFont)
latexMeasLabel.SetTextSize(text_size)
latexMeasLabel.SetTextColor(kGray+2)
latexMeasLabel.DrawLatex(err_leg.GetPointX(0)-err_size*1.25, err_leg_y+0.02, "total")
latexMeasLabel.DrawLatex(err_leg.GetPointX(0)-err_size*0.55, err_leg_y-0.032, "stat")

gPad.Update()

c.Print('plots/mtopSummary.pdf')
