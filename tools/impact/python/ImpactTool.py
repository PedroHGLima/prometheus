__all__ = ['ImpactTool']


# core includes
from Gaugi import StatusCode
from Gaugi import Algorithm
from Gaugi import mkdir_p
from Gaugi import GeV
from Gaugi import progressbar
from Gaugi.tex.TexAPI       import *
from Gaugi.tex.BeamerAPI    import *
from Gaugi.monet.AtlasStyle import SetAtlasStyle
from prometheus import Dataframe as DataframeEnum


# tool includes
from EventSelectionTool import RetrieveBinningIdx
from ProfileTools.analysis.constants import *
from ROOT import TH1F
from functools import reduce
from itertools import product
import os, gc, time, math
import numpy as np

from ImpactTools.drawers import *

#
# Analysis tool
#
class ImpactTool( Algorithm ):

  #
  # Constructor
  #
  def __init__(self, name, dataframe, selection_list_labels=None, **kw):

    
    Algorithm.__init__(self, name)

    # declare all properties with default values  
    self.declareProperty( "Basepath", "Event/ImpactTool", "Impact base path histogram." )
    self.declareProperty( "EtBinningValues" , [], "Et bin selection for data selection." )
    self.declareProperty( "EtaBinningValues", [], "Et bin selection for data selection." )


    # Set all properties values from the contructor args
    for key, value in kw.items():
      self.setProperty( key, value )

    if selection_list_labels is None:
      # default selection names definition
      self.__selections = [
                      'ringer',
                      'no_ringer',
                    ]
    else:
      self.__selections = selection_list_labels

    self.__selectionFeatures = list()
   

  #
  # Add selection configuration
  #
  def add_selection( self, name_a, expression_a, name_b, expression_b):

    self.__selectionFeatures.append( SelectionConfig(name_a, expression_a, name_b, expression_b) )


  #
  # Set et binning values 
  #
  def setEtBinningValues( self, etbins ):
    self.setProperty( "EtBinningValues", etbins )
 
  #
  # Set eta binning values
  #
  def setEtaBinningValues( self, etabins ):
    self.setProperty( "EtaBinningValues", etabins )


  #
  # Initialize method
  #
  def initialize(self):
    
    Algorithm.initialize(self)
    sg = self.getStoreGateSvc()

    basepath = self.getProperty("Basepath")
    etBins = self.getProperty( "EtBinningValues" )
    etaBins = self.getProperty( "EtaBinningValues" )

    etabins = default_etabins

    for feat in self.__selectionFeatures:
      # hold selection name
      selection_name = feat.name_a()+'_VS_'+feat.name_b()

      ### loopover ets...
      for etBinIdx in range(len(etBins)-1):
        ### loop over etas...
        for etaBinIdx in range(len(etaBins)-1):
          ### loop over selections...
          for selection in self.__selections:
            # hold binning name
            binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx)
              
            dirname = basepath+'/'+selection_name+'/'+binning_name+'/'+selection
            sg.mkdir( dirname )
            sg.addHistogram(TH1F('et', ('%s;%s;Count')%(basicInfoQuantities['et'],basicInfoQuantities['et']),
            basicInfoNBins['et'],basicInfoLowerEdges['et'],basicInfoHighEdges['et']) )
            sg.addHistogram(TH1F('eta' , ('%s;%s;Count')%(basicInfoQuantities['eta'],basicInfoQuantities['eta']),
                            len(etabins)-1, np.array(etabins)) )
            sg.addHistogram(TH1F('phi', ('%s;%s;Count')%(basicInfoQuantities['phi'],basicInfoQuantities['phi']),
                            20, -3.2, 3.2) )
            sg.addHistogram(TH1F('avgmu', ('%s;%s;Count')%(basicInfoQuantities['avgmu'],basicInfoQuantities['avgmu']),
                            16,0,80) )
            sg.addHistogram(TH1F('nvtx', ('%s;%s;Count')%(basicInfoQuantities['nvtx'],basicInfoQuantities['nvtx']),
                            len(nvtx_bins)-1,np.array(nvtx_bins)) )

            # we want measure this
            for itype in ['off', 'hlt']:
              for key in ['f1', 'f3', 'weta2', 'wtots1', 'rphi', 'deltaEta1', 'deltaPhiRescaled2', 'eratio', 'reta', 'rhad']:
                sg.addHistogram(TH1F('%s_%s' %(itype, key), 
                  ('%s;%s;Count')%(electronQuantities[key],electronQuantities[key]),
                  standardQuantitiesNBins[key],
                  standardQuantitiesLowerEdges[key],
                  standardQuantitiesHighEdges[key]))

            for fc_key in ['f1', 'f3', 'weta2', 'reta', 'eratio']:
              sg.addHistogram(TH1F('fc_%s' %(fc_key), 
                  ('%s;%s;Count')%(electronQuantities[fc_key],electronQuantities[fc_key]),
                  standardQuantitiesNBins[fc_key],
                  standardQuantitiesLowerEdges[fc_key],
                  standardQuantitiesHighEdges[fc_key]))
            # loop over selections

    # loop over pairs
    self.init_lock()
    return StatusCode.SUCCESS



  #
  # Execute method
  #
  def execute(self, context):
    

    basepath = self.getProperty("Basepath")
    etBins = self.getProperty( "EtBinningValues" )
    etaBins = self.getProperty( "EtaBinningValues" )


    # Retrieve container
    if self._dataframe is DataframeEnum.Electron_v1:
      off_elCont    = context.getHandler( "ElectronContainer" )
      hlt_elCont    = context.getHandler( "HLT__ElectronContainer" )
      fastcalo_Cont = context.getHandler( "HLT__TrigEMClusterContainer")
    elif self._dataframe is DataframeEnum.Photon_v1:
      elCont    = context.getHandler( "PhotonContainer" )
    else:
      elCont    = context.getHandler( "ElectronContainer" )
    
    evt = context.getHandler( "EventInfoContainer" )
    eta = math.fabs(off_elCont.eta())
    et = off_elCont.et()/GeV
    
    # get tracks containers
    # Retrieve container
    if self._dataframe is DataframeEnum.Electron_v1:
      off_track = off_elCont.trackParticle()
      track = False
    elif self._dataframe is DataframeEnum.Photon_v1:
      track = False

    dec = context.getHandler( "MenuContainer" )

    evt = context.getHandler("EventInfoContainer")
    pw = evt.MCPileupWeight()
    sg = self.getStoreGateSvc()

    if et < etBins[0]:
      return StatusCode.SUCCESS
    if eta > 2.47:
      return StatusCode.SUCCESS

    etBinIdx, etaBinIdx = RetrieveBinningIdx(et,eta,etBins, etaBins, logger=self._logger )
    binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx)

    for feat in self.__selectionFeatures:
      # create a list of histogram paths
      dir_list       = []
      selection_name = feat.name_a()+'_VS_'+feat.name_b()
      # check the decision for this expression
      passed_a       = bool(dec.accept( feat.expression_a() ))
      if passed_a: # check if passed for expression_a
        # if passed then append the path in dir_list
        dir_list.append(basepath+'/'+selection_name+'/'+binning_name+'/ringer')
      # now do the same for expression_b
      passed_b       = bool(dec.accept( feat.expression_b() ))
      if passed_b:
        dir_list.append(basepath+'/'+selection_name+'/'+binning_name+'/no_ringer')
      
      # now loop over the dir_list to fill all histograms
      for ipath in dir_list:
        # getting the path to fill.
        dirname  = ipath

        pw=1
        # Fill basic infos
        sg.histogram(dirname+'/et').Fill(et,pw)
        sg.histogram(dirname+'/eta').Fill(off_elCont.eta(),pw)
        sg.histogram(dirname+'/phi').Fill(off_elCont.phi(),pw)
        sg.histogram(dirname+'/avgmu').Fill(evt.avgmu(),pw)
        sg.histogram(dirname+'/nvtx').Fill(evt.nvtx(),pw)
        # Fill shower shapes
        # off
        sg.histogram(dirname+'/off_f1').Fill(off_elCont.f1(),pw)
        sg.histogram(dirname+'/off_f3').Fill(off_elCont.f3(),pw)
        sg.histogram(dirname+'/off_weta2').Fill(off_elCont.weta2(),pw)
        sg.histogram(dirname+'/off_wtots1').Fill(off_elCont.wtots1(),pw)
        sg.histogram(dirname+'/off_reta').Fill(off_elCont.reta(),pw)
        sg.histogram(dirname+'/off_rhad').Fill(off_elCont.rhad(),pw)
        sg.histogram(dirname+'/off_rphi').Fill(off_elCont.rphi(),pw)
        sg.histogram(dirname+'/off_eratio').Fill(off_elCont.eratio(),pw)
        sg.histogram(dirname+'/off_deltaEta1').Fill(off_elCont.deltaEta1(),pw)
        sg.histogram(dirname+'/off_deltaPhiRescaled2').Fill(off_elCont.deltaPhiRescaled2(),pw)
        # hlt
        sg.histogram(dirname+'/hlt_f1').Fill(hlt_elCont.f1(),pw)
        sg.histogram(dirname+'/hlt_f3').Fill(hlt_elCont.f3(),pw)
        sg.histogram(dirname+'/hlt_weta2').Fill(hlt_elCont.weta2(),pw)
        sg.histogram(dirname+'/hlt_wtots1').Fill(hlt_elCont.wtots1(),pw)
        sg.histogram(dirname+'/hlt_reta').Fill(hlt_elCont.reta(),pw)
        sg.histogram(dirname+'/hlt_rhad').Fill(hlt_elCont.rhad(),pw)
        sg.histogram(dirname+'/hlt_rphi').Fill(hlt_elCont.rphi(),pw)
        sg.histogram(dirname+'/hlt_eratio').Fill(hlt_elCont.eratio(),pw)
        sg.histogram(dirname+'/hlt_deltaEta1').Fill(hlt_elCont.deltaEta1(),pw)
        sg.histogram(dirname+'/hlt_deltaPhiRescaled2').Fill(hlt_elCont.deltaPhiRescaled2(),pw)

        # fastcalo
        sg.histogram(dirname+'/fc_f1').Fill(fastcalo_Cont.f1(),pw)
        sg.histogram(dirname+'/fc_f3').Fill(fastcalo_Cont.f3(),pw)
        sg.histogram(dirname+'/fc_weta2').Fill(fastcalo_Cont.weta2(),pw)
        #sg.histogram(dirname+'/fc_wtots').Fill(fastcalo_Cont.wtots(),pw)
        sg.histogram(dirname+'/fc_reta').Fill(fastcalo_Cont.reta(),pw)
        sg.histogram(dirname+'/fc_eratio').Fill(fastcalo_Cont.eratio(),pw)

        # Fill track variables
        if track:
          # off
          sg.histogram(dirname+'/off_trackd0pvunbiased').Fill(off_track.d0(),pw)
          sg.histogram(dirname+'/off_d0significance').Fill(off_track.d0significance(),pw)
          sg.histogram(dirname+'/off_eProbabilityHT').Fill(off_track.eProbabilityHT(),pw)
          sg.histogram(dirname+'/off_TRT_PID').Fill(off_track.trans_TRT_PID(),pw)
          sg.histogram(dirname+'/off_DeltaPOverP').Fill(off_track.DeltaPOverP(),pw)
          # hlt
          sg.histogram(dirname+'/hlt_trackd0pvunbiased').Fill(hlt_track.d0(),pw)
          sg.histogram(dirname+'/hlt_d0significance').Fill(hlt_track.d0significance(),pw)
          sg.histogram(dirname+'/hlt_eProbabilityHT').Fill(hlt_track.eProbabilityHT(),pw)
          sg.histogram(dirname+'/hlt_TRT_PID').Fill(hlt_track.trans_TRT_PID(),pw)
          sg.histogram(dirname+'/hlt_DeltaPOverP').Fill(hlt_track.DeltaPOverP(),pw)
     

    return StatusCode.SUCCESS

  
  #
  # Finalize method
  #
  def finalize(self):
    self.fina_lock()
    return StatusCode.SUCCESS


  #
  # Standalone plot method
  #
  def plot(self, dirnames, pdfoutputs, pdftitles, runLabel='' ,doPDF=True):
    


    SetAtlasStyle()
    beamer_plots = {}
    global tobject_collector

    basepath = self.getProperty("Basepath")
    etBins = self.getProperty( "EtBinningValues" )
    etaBins = self.getProperty( "EtaBinningValues" )




    for idx, feat in enumerate(self.__selectionFeatures):
      
      dirname = os.getcwd()+'/'+dirnames[idx]
      mkdir_p(dirname)
      # hold selection name
      selection_name = feat.name_a()+'_Vs_'+feat.name_b()
      # For beamer... 
      if not selection_name in beamer_plots.keys():
        beamer_plots[selection_name]={}
        beamer_plots[selection_name]['integrated']={}

      ### Plot binning plots  
      if (len(etBins) * len(etaBins)) > 1:
        for etBinIdx, etaBinIdx in progressbar(product(range(len(etBins)-1),range(len(etaBins)-1)),
                                               (len(etBins)-1)*(len(etaBins)-1),
                                               prefix = "Plotting... ", logger=self._logger):
          # hold binning name
          binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx)
          # for beamer...
          if not binning_name in beamer_plots[selection_name].keys():
            beamer_plots[selection_name][binning_name]={}
          
          ### loop over standard quantities
          for key in standardQuantitiesNBins.keys(): 
            outname = dirname+'/'+selection_name.replace('_Vs_','_')+'_'+ key + '_' + binning_name
            out = PlotQuantities(basepath+'/'+selection_name+'/'+binning_name, key, 
                outname,etidx=etBinIdx,etaidx=etaBinIdx,xlabel=electronQuantities[key],divide='b',runLabel=runLabel)
            beamer_plots[selection_name][binning_name][key] = out
            #del tobject_collector[:]
        
          ### loop over info quantities
          for key in basicInfoQuantities.keys():
            outname = dirname+'/'+selection_name.replace('_Vs_','_')+'_'+ key + '_' + binning_name
            out = PlotQuantities(basepath+'/'+selection_name+'/'+binning_name, key, 
                outname, etidx=etBinIdx,etaidx=etaBinIdx,xlabel=basicInfoQuantities[key],divide='b', runLabel=runLabel)
            beamer_plots[selection_name][binning_name][key] = out
            #del tobject_collector[:]

          
          beamer_plots[selection_name][binning_name]['statistics'] = GetStatistics(basepath+'/'+selection_name+'/'+binning_name, \
                                                                                        'avgmu',etidx=etBinIdx,etaidx=etaBinIdx)
      

      #### Plot integrated histograms
      ### loop over standard quantities
      for key in standardQuantitiesNBins.keys(): 
        outname = dirname+'/'+selection_name.replace('_Vs_','_')+'_'+ key
        out = PlotQuantities(basepath+'/'+selection_name, key, 
              outname,xlabel=electronQuantities[key],divide='b',runLabel=runLabel,
              addbinlines=True)
        beamer_plots[selection_name]['integrated'][key] = out
        tobject_collector = []
        gc.collect()
      ### loop over info quantities
      for key in basicInfoQuantities.keys():
        outname = dirname+'/'+selection_name.replace('_Vs_','_')+'_'+ key + '_' + binning_name
        out = PlotQuantities(basepath+'/'+selection_name, key, 
            outname,xlabel=basicInfoQuantities[key],divide='b', runLabel=runLabel,
            addbinlines=True)
        beamer_plots[selection_name]['integrated'][key] = out
        tobject_collector = []
        gc.collect()
      
      beamer_plots[selection_name]['integrated']['statistics'] = GetStatistics(basepath+'/'+selection_name, 'avgmu')




