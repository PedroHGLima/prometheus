
__all__ = ["PlotQuantities", "GetStatistics","QuadrantConfig" ]

from Gaugi.monet.PlotFunctions import *
from Gaugi.monet.TAxisFunctions import *
from copy import copy
from itertools import product
import ROOT



def AddTopLabels(can,legend, legOpt = 'p', quantity_text = '', etlist = None
                     , etalist = None, etidx = None, etaidx = None, legTextSize=10
                     , runLabel = '', extraText1 = None, legendY1=.66, legendY2=.97
                     , maxLegLength = 19, logger=None):
    text_lines = []
    text_lines += [GetAtlasInternalText()]
    text_lines.append( GetSqrtsText(13) )
    _etlist = copy(etlist)
    _etalist = copy(etalist)
    if runLabel: text_lines.append( runLabel )
    if extraText1: text_lines.append( extraText1 )
    DrawText(can,text_lines,.40,.68,.70,.93,totalentries=4)
    if legend:
        MakeLegend( can,.67,legendY1,.93,legendY2,textsize=legTextSize
                  , names=legend, option = legOpt, squarebox=False
                  , totalentries=0, maxlength=maxLegLength )
    try:
        extraText = []
        if _etlist and etidx is not None:
            # add infinity in case of last et value too large
            if _etlist[-1]>9999:  _etlist[-1]='#infty'
            binEt = (str(_etlist[etidx]) + ' < E_{T} [GeV] < ' + str(_etlist[etidx+1]) if etidx+1 < len(_etlist) else
                                     'E_{T} > ' + str(_etlist[etidx]) + ' GeV')
            extraText.append(binEt)
        if quantity_text: 
            if not isinstance(quantity_text,(tuple,list)): quantity_text = [quantity_text]
            extraText += quantity_text
        if _etalist and etaidx is not None:
            binEta = (str(_etalist[etaidx]) + ' < #eta < ' + str(_etalist[etaidx+1]) if etaidx+1 < len(_etalist) else
                                        str(_etalist[etaidx]) + ' < #eta < 2.47')
            extraText.append(binEta)
        DrawText(can,extraText,.14,.68,.35,.93,totalentries=4)
    except NameError as e:
        if logger:
          logger.warning("Couldn't print test due to error: %s", e)
        pass


class QuadrantConfig(object):
  def __init__(self, name_a, expression_a, name_b, expression_b):
    self._name_a=name_a; self._name_b=name_b
    self._expression_a=expression_a; self._expression_b=expression_b
  
  def name_a(self):
    return self._name_a

  def expression_a(self):
    return self._expression_a

  def name_b(self):
    return self._name_b

  def expression_b(self):
    return self._expression_b




#
# Plot quadrant 
#
def PlotQuantities( sg, basepath, key, outname, legends, drawopt='hist', divide='B', etidx=None, etaidx=None, xlabel='', runLabel='',addbinlines=False, etBins=[], etaBins=[]):


  import ROOT
  ROOT.gROOT.SetBatch(ROOT.kTRUE)
  ROOT.gErrorIgnoreLevel=ROOT.kWarning
  ROOT.TH1.AddDirectory(ROOT.kFALSE)

  from Gaugi.monet.utils import sumHists


  if (etidx is not None) and (etaidx is not None):
    hists = [
              sg.histogram(basepath+'/passed_passed/'+key),
              sg.histogram(basepath+'/passed_rejected/'+key),
              sg.histogram(basepath+'/rejected_passed/'+key),
              sg.histogram(basepath+'/rejected_rejected/'+key)
            ]
  else:
    passed_passed = []; passed_rejected = []; rejected_passed = []; rejected_rejected = []
    for etBinIdx, etaBinIdx in product(range(len(etBins)-1),range(len(etaBins)-1)):
      binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx) 
      passed_passed.append( sg.histogram(basepath+'/'+binning_name+'/passed_passed/'+key) ) 
      passed_rejected.append( sg.histogram(basepath+'/'+binning_name+'/passed_rejected/'+key) )
      rejected_passed.append( sg.histogram(basepath+'/'+binning_name+'/rejected_passed/'+key) )
      rejected_rejected.append( sg.histogram(basepath+'/'+binning_name+'/rejected_rejected/'+key) )

    hists = [
              sumHists(passed_passed),
              sumHists(passed_rejected),
              sumHists(rejected_passed),
              sumHists(rejected_rejected),
            ]

  ref_hist = sumHists( hists )
  #ref_hist = sumHists( [hists[1],hists[2]])

  from ROOT import kBlack,kRed,kGreen,kGray,kMagenta,kBlue
  outcan = RatioCanvas( outname, outname, 500, 500)
  pad_top = outcan.GetPrimitive('pad_top')
  pad_bot = outcan.GetPrimitive('pad_bot')

  pad_top.SetLogy()
  collect=[]
  divs=[]
  #outcan.GetPrimitive('pad_bot').SetLogy()
  these_colors = [kBlack,kRed+1, kBlue+2,kGray+1]
  #these_colors = [kBlack,kGray+1]
  these_transcolors=[]
  for c in these_colors:
    these_transcolors.append(ROOT.TColor.GetColorTransparent(c, .5))

  divs = []
  for idx, hist in enumerate(hists):
    
    hist.SetMarkerStyle(20)
    hist.SetMarkerSize(0.55)
    hist.SetLineColor(these_colors[idx])
    hist.SetMarkerColor(these_colors[idx])
    hist.SetFillColor(these_transcolors[idx])
    #hist.SetFillColor(ROOT.TColor.GetColorTransparent( these_colors[idx], .5)  )
    AddHistogram( pad_top, hist, 'histE2 L same', False, None, None)
    #AddHistogram( pad_top, hist, 'hist L same', False, None, None)


    div = hist.Clone(); div.Divide(div,ref_hist,1.,1.,'b'); div.Scale(100.); collect.append(div)
    div.SetMarkerSize(0.5)
    div.SetMarkerColor(these_colors[idx])
    # TODO: Check why error bar still here. Force error bar equal zero
    for ibin in range(div.GetNbinsX()):  div.SetBinError(ibin,0.0)
    divs.append( div )
    # add left axis
    if idx == 1 or idx == 2: AddHistogram( pad_bot, div , 'p', False, None, None)
    #if idx == 2: AddHistogram( pad_bot, div , 'p', False, None, None)


  AddTopLabels(outcan, legends, runLabel=runLabel, legOpt='p',legTextSize=12,
               etlist=etBins,etalist=etaBins,etidx=etidx,etaidx=etaidx)

  SetAxisLabels(outcan,xlabel,'Count','Disagreement [%]')
  #SetAxisLabels(outcan,xlabel,'Count','(Red or Blue)/Total [%]')
  FormatCanvasAxes(outcan, XLabelSize=18, YLabelSize=18, XTitleOffset=0.87, YTitleOffset=1.5, YTitleSize=16)
  #FormatCanvasAxes(outcan, XLabelSize=18, YLabelSize=18, XTitleOffset=0.87, YTitleOffset=1.5)
  AutoFixAxes(pad_top,ignoreErrors=False)
  AutoFixAxes(pad_bot,ignoreErrors=False)
  FixYaxisRanges(pad_bot, ignoreErrors=True, yminc=-eps )
  
  if addbinlines:
    AddBinLines(pad_top,hists[0],useHistMax=True,horizotalLine=0.)
    
  #AddRightAxisObj(pad_bot, [divs[1]], drawopt="p,same", equate=[0., max([d.GetBinContent(h.GetMaximumBin()) for d,h in zip([divs[1]], [hists[1]])])]
  #               , drawAxis=True, axisColor=(ROOT.kRed+1), ignorezeros=False, ignoreErrors=True, label = "Ringer Rejected [%]")

  outcan.SaveAs( outname+'.C' ) 
  outname_png = outname+'.png'
  outname = outname+'.pdf'
  outcan.SaveAs( outname ) 
  outcan.SaveAs( outname_png ) 
  return outname




def GetStatistics( sg ,basepath, key, etidx=None, etaidx=None, etBins=[], etaBins=[] ):
  
  # get all quadrant histograms
  from Gaugi.monet.utils import sumHists
  if (etidx is not None) and (etaidx is not None):
    hists = [
              sg.histogram(basepath+'/passed_passed/'+key),
              sg.histogram(basepath+'/passed_rejected/'+key),
              sg.histogram(basepath+'/rejected_passed/'+key),
              sg.histogram(basepath+'/rejected_rejected/'+key)
            ]
  else:
    from itertools import product
    passed_passed = []; passed_rejected = []; rejected_passed = []; rejected_rejected = []
    for etBinIdx, etaBinIdx in product(range(len(etBins)-1),range(len(etaBins)-1)):
      binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx) 
      passed_passed.append( sg.histogram(basepath+'/'+binning_name+'/passed_passed/'+key) ) 
      passed_rejected.append( sg.histogram(basepath+'/'+binning_name+'/passed_rejected/'+key) )
      rejected_passed.append( sg.histogram(basepath+'/'+binning_name+'/rejected_passed/'+key) )
      rejected_rejected.append( sg.histogram(basepath+'/'+binning_name+'/rejected_rejected/'+key) )

    hists = [
              sumHists(passed_passed),
              sumHists(passed_rejected),
              sumHists(rejected_passed),
              sumHists(rejected_rejected),
            ]

  # NOTE: Follow the statistics definitions for each case.
  # passed is 1 and rejected is zero
  # expression A is i and expression B is j
  # Contigency table:
  #      | hi=0  hi=1
  # hj=0 |  a     c
  # hj=1 |  b     d
  a = hists[3].GetEntries()
  d = hists[0].GetEntries()
  b = hists[2].GetEntries()
  c = hists[1].GetEntries()
  m = a+b+c+d
  
  Qij=0;Pij=0;Kp=0;dis_ij=0

  try:
    # Q statistics
    Qij = (a*d-b*c) / float(a*d+b*c)

    # correlation coef
    Pij = (a*d-b*c)/np.sqrt(( (a+b)*(a+c)*(c+d)*(b+d) ))

    # kappa-statistics
    Q1 = (a+d)/float(m)
    Q2 = ( (a+b)*(a+c)+(c+d)*(b+d) ) / float(m*m)
    Kp = (Q1-Q2)/ float(1-Q2)
    
    dis_ij = (b+c)/float(m)
  except:
    pass
  return  {'Qij':Qij,'Pij':Pij,'Kp':Kp, 'dis_ij':dis_ij}





