import ROOT, copy, os, sys
import MuonPerformance.MuonAnalyser.CMS_lumi as CMS_lumi
import MuonPerformance.MuonAnalyser.tdrstyle as tdrstyle
from MuonPerformance.MuonAnalyser.histoHelper import *
ROOT.gROOT.SetBatch(True)
tdrstyle.setTDRStyle()

def setMarkerStyle(h,color,style):
    h.SetMarkerColor(color)
    h.SetMarkerStyle(style)
    h.SetMarkerSize(1.5)
    h.SetLineColor(color)
    h.SetLineWidth(2)

def drawSampleName(samplename):
    tex2 = ROOT.TLatex()
    tex2.SetNDC()
    tex2.SetTextFont(62)
    tex2.SetTextSize(0.03)
    #tex2.DrawLatex(0.2, 0.2, samplename)
    #tex2.DrawLatex(0.2, 0.9, samplename)
    tex2.DrawLatex(0.63, 0.2, samplename)
    #tex2.DrawLatex(0.65, 0.88, samplename)

def draw(h_init, y_name, hlists, name, text):
    #Plot style
    setMarkerStyle(hlists[0], 1, 32) 
    setMarkerStyle(hlists[1], 1, 22) 
    setMarkerStyle(hlists[2], 3, 32) 
    setMarkerStyle(hlists[3], 3, 22) 
    setMarkerStyle(hlists[4], 4, 32) 
    setMarkerStyle(hlists[5], 4, 22) 
    setMarkerStyle(hlists[6], 6, 32) 
    setMarkerStyle(hlists[7], 6, 22) 

    #Set canvas
    #canv = makeCanvas(plotvar+name, False)
    #setMargins(canv, False)
    canv = ROOT.TCanvas()
    #canv.SetGrid()
    h_init.GetYaxis().SetTitle(y_name)
    h_init.Draw()
    drawSampleName(text)

    #Legend and drawing
    leg = ROOT.TLegend(0.17,0.16,0.61,0.462)
    #leg = ROOT.TLegend(0.55,0.16,0.9,0.45)
    for h in hlists:
        h.Draw("e1same")
        leg.AddEntry(h,h.GetTitle(),"pl")
    hlists[0].Draw("e1same")
    leg.SetTextFont(62)
    leg.SetTextSize(0.025)
    leg.SetBorderSize(0)
    #leg.SetFillStyle(0)
    leg.Draw()

    #CMS_lumi setting
    iPos = 0
    iPeriod = 0
    if( iPos==0 ): CMS_lumi.relPosX = 0.12
    CMS_lumi.extraText = "Preliminary"
    CMS_lumi.lumi_sqrtS = ""
    CMS_lumi.CMS_lumi(canv, iPeriod, iPos)

    #canv.SetLogx()
    canv.Modified()
    canv.Update()
    canv.SaveAs("eff_%s.png"%(name))
    canv.SaveAs("eff_%s.pdf"%(name))

datadir = './'
fname = "sliceTest"

#Set extra text
#text = "Slice Test 2018C"
#text = "#splitline{Slice Test 2018C}{All Muon, #Delta BX = 0, #phi Fiducial}"
#text = "#splitline{#splitline{Slice Test 2018C}{Muon p_{T} > 20, #Delta BX = 0}}{}"
text = "#splitline{#splitline{Slice Test 2018C}{Muon p_{T} > 20, #Delta BX = 0,}}{#phi Fiducial}"
#plotvar = "pt"
plotvar = "inPhi"

f = ROOT.TFile(datadir+fname+".root")
t = f.Get("SliceTestEfficiencyAnalysis")
hl_eff = []; hl_pos = []; hl_pos_matched = [];
c = ROOT.TCanvas()
for ch in [27,28,29,30]:
    for lay in [1,2]:
        h = t.Get("%s ch %d lay %d"%(plotvar,ch,lay))
        h2 = t.Get("%s_matched ch %d lay %d"%(plotvar,ch,lay))
        #h.Rebin(2)
        #h2.Rebin(2)
        eff = ROOT.TEfficiency(h2,h)
        mean = 0
        if h.Integral() !=0: mean = h2.Integral()/float(h.Integral())*100
        eff.SetTitle("Chamber %d Layer %d (%.2f%%)"%(ch,lay,mean))
        eff.SetLineWidth(2)
        hl_eff.append(eff)


if plotvar == "pt": h_init = ROOT.TH1F("","pT Efficiency;Muon p_{T} [GeV];",20,0,200)
if plotvar == "inStrip": h_init = ROOT.TH1F("","Strip Efficiency;Strip;",24,0,384)
if plotvar == "inRoll": h_init = ROOT.TH1F("","iEta Efficiency;iEta;",8,0.5,8.5)
if plotvar == "inPhi": h_init = ROOT.TH1F("","Phi Efficiency;#phi;",36,-1.85,-1.1)
h_init.GetYaxis().SetTitleOffset(1)
h_init.GetXaxis().SetTitleOffset(1.1)
h_init.GetYaxis().SetTitleOffset(1.5)
h_init.GetXaxis().SetTitleSize(0.05)
h_init.GetXaxis().SetLabelSize(0.038)
h_init.GetYaxis().SetTitleSize(0.05)
h_init.GetYaxis().SetLabelSize(0.038)

h_init.SetMaximum(1)
h_init.SetMinimum(0)
draw(h_init, "Efficiency ", hl_eff, plotvar, text)

