
__all__ = ['ImpactTool']


# core
from Gaugi import GeV
from Gaugi import Algorithm
from Gaugi import StatusCode
from Gaugi.messenger.macros import *
from prometheus import Dataframe as DataframeEnum


# External
from ROOT import TH1F, TH2F, TProfile, TProfile2D
from ProfileTools import zee_etbins, jpsiee_etbins, default_etabins, nvtx_bins, high_nvtx_bins
from TrigEgammaEmulationTool import Chain, Group, TDT
import numpy as np

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
# Efficiency tool
#
class ImpactTool( Algorithm ):
  
  __triggerLevels = ['L1Calo','L2Calo','L2','EFCalo','HLT']

  #
  # Constructor
  #
  def __init__(self, name, dojpsiee=False, **kw):

    Algorithm.__init__(self, name)
    self.__groups = list()

    # declare all props here
    self.declareProperty( "Basepath", "Event/ImpactTool", "Histograms base path for the efficiency tool"      )
    self.declareProperty( "DoJpisee", dojpsiee                 , "Use the J/psiee et bins in the eff et histograms." )
    self.declareProperty( "EtBinningValues" , [], "Et bin selection for data selection." )
    self.declareProperty( "EtaBinningValues", [], "Et bin selection for data selection." )
 
    # Set property values using the constructor args
    for key, value in kw.items():
      self.setProperty(key, value)

  #
  # Add trigger group to the monitoring list
  #
  def addGroup( self, group ): 

    from Gaugi import ToolSvc
    emulator = ToolSvc.retrieve( "Emulator" )
    if not emulator.isValid( group.chain().name() ):
      emulator+=group.chain()
    self.__groups.append( group )  


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
    
    basepath = self.getProperty( "Basepath" )
    doJpsiee = self.getProperty( "DoJpisee" )
    etBins = self.getProperty( "EtBinningValues" )
    etaBins = self.getProperty( "EtaBinningValues" )

    sg = self.getStoreGateSvc()
    
    #et_bins  = zee_etbins
    eta_bins = default_etabins
    nvtx_bins.extend(high_nvtx_bins)
    #eta_bins = [0,0.6,0.8,1.15,1.37,1.52,1.81,2.01,2.37,2.47]
    et_bins = jpsiee_etbins if doJpsiee else [4.,7.,10.,15.,20.,25.,30.,35.,40.,45.,50.,60.,80.,150.] 
    
    etabins = default_etabins

    for group in self.__groups:
      # Get the chain object
      chain = group.chain()
      ### loopover ets...
      for etBinIdx in range(len(etBins)-1):
        ### loop over etas...
        for etaBinIdx in range(len(etaBins)-1):
          binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx)
        
          m_dirname = basepath+'/'+chain.name()+ '/' + binning_name
          sg.mkdir( m_dirname )

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

    self.init_lock()
    return StatusCode.SUCCESS 

    

  def execute(self, context):
  

    basepath = self.getProperty( "Basepath" )
    etBins = self.getProperty( "EtBinningValues" )
    etaBins = self.getProperty( "EtaBinningValues" )

    # Retrieve Electron container
    if self._dataframe is DataframeEnum.Electron_v1:
      off_elCont    = context.getHandler( "ElectronContainer" )
      hlt_elCont    = context.getHandler( "HLT__ElectronContainer" )
      fastcalo_Cont = context.getHandler( "HLT__TrigEMClusterContainer")
    elif self._dataframe is DataframeEnum.Photon_v1:
      elCont    = context.getHandler( "PhotonContainer" )
    else:
      elCont    = context.getHandler( "ElectronContainer" )
      
    dec = context.getHandler( "MenuContainer" )

    for group in self.__groups:
      chain = group.chain()
      accept = dec.accept( chain.name() )
      if accept.getCutResult("HLT"):
        for off_el, hlt_el, fc_cluster in zip(off_elCont, hlt_elCont, fastcalo_Cont):

          eta = math.fabs(off_el.eta())
          et = off_el.et()/GeV
          # get the bin idx
          etBinIdx, etaBinIdx = RetrieveBinningIdx(et,eta,etBins, etaBins, logger=self._logger )
          binning_name = ('et%d_eta%d') % (etBinIdx,etaBinIdx)

          if off_el.et()  < (group.etthr()- 5)*GeV:  continue 
          if abs(off_el.eta())>2.47: continue
          if (etaBinIdx == -1) or (etBinIdx == -1): continue

          dirname = basepath+'/'+chain.name()+'/' + binning_name
          
          # fill histograms
          self.fillHistograms(dirname, off_el, hlt_el, fc_cluster, group.etthr(), group.pidname())  

    return StatusCode.SUCCESS 

  #
  # Fill histograms
  #
  def fillHistograms( self, dirname, off_obj, hlt_obj, fastcalo_obj, etthr, pidword ):
  
    sg = self.getStoreGateSvc()
    
    pid = off_obj.accept(pidword) if pidword else True

    eta = off_obj.caloCluster().etaBE2()
    et = off_obj.et()/GeV
    phi = off_obj.phi()
    evt = self.getContext().getHandler("EventInfoContainer")

    avgmu = evt.avgmu()
    nvtx = evt.nvtx()
    pw = evt.MCPileupWeight()

    if pid: 
      if et > etthr+1.0:
        # Fill basic infos
        sg.histogram(dirname+'/et').Fill(et,pw)
        sg.histogram(dirname+'/eta').Fill(off_obj.eta(),pw)
        sg.histogram(dirname+'/phi').Fill(off_obj.phi(),pw)
        sg.histogram(dirname+'/avgmu').Fill(evt.avgmu(),pw)
        sg.histogram(dirname+'/nvtx').Fill(evt.nvtx(),pw)
        # Fill shower shapes
        # off
        sg.histogram(dirname+'/off_f1').Fill(off_obj.f1(),pw)
        sg.histogram(dirname+'/off_f3').Fill(off_obj.f3(),pw)
        sg.histogram(dirname+'/off_weta2').Fill(off_obj.weta2(),pw)
        sg.histogram(dirname+'/off_wtots1').Fill(off_obj.wtots1(),pw)
        sg.histogram(dirname+'/off_reta').Fill(off_obj.reta(),pw)
        sg.histogram(dirname+'/off_rhad').Fill(off_obj.rhad(),pw)
        sg.histogram(dirname+'/off_rphi').Fill(off_obj.rphi(),pw)
        sg.histogram(dirname+'/off_eratio').Fill(off_obj.eratio(),pw)
        sg.histogram(dirname+'/off_deltaEta1').Fill(off_obj.deltaEta1(),pw)
        sg.histogram(dirname+'/off_deltaPhiRescaled2').Fill(off_obj.deltaPhiRescaled2(),pw)
        # hlt
        sg.histogram(dirname+'/hlt_f1').Fill(hlt_obj.f1(),pw)
        sg.histogram(dirname+'/hlt_f3').Fill(hlt_obj.f3(),pw)
        sg.histogram(dirname+'/hlt_weta2').Fill(hlt_obj.weta2(),pw)
        sg.histogram(dirname+'/hlt_wtots1').Fill(hlt_obj.wtots1(),pw)
        sg.histogram(dirname+'/hlt_reta').Fill(hlt_obj.reta(),pw)
        sg.histogram(dirname+'/hlt_rhad').Fill(hlt_obj.rhad(),pw)
        sg.histogram(dirname+'/hlt_rphi').Fill(hlt_obj.rphi(),pw)
        sg.histogram(dirname+'/hlt_eratio').Fill(hlt_obj.eratio(),pw)
        sg.histogram(dirname+'/hlt_deltaEta1').Fill(hlt_obj.deltaEta1(),pw)
        sg.histogram(dirname+'/hlt_deltaPhiRescaled2').Fill(hlt_obj.deltaPhiRescaled2(),pw)

        # fastcalo
        sg.histogram(dirname+'/fc_f1').Fill(fastcalo_obj.f1(),pw)
        sg.histogram(dirname+'/fc_f3').Fill(fastcalo_obj.f3(),pw)
        sg.histogram(dirname+'/fc_weta2').Fill(fastcalo_obj.weta2(),pw)
        #sg.histogram(dirname+'/fc_wtots').Fill(fastcalo_Cont.wtots(),pw)
        sg.histogram(dirname+'/fc_reta').Fill(fastcalo_obj.reta(),pw)
        sg.histogram(dirname+'/fc_eratio').Fill(fastcalo_obj.eratio(),pw)
        #sg.histogram( dirname+'/eta' ).Fill(eta, pw)
        #sg.histogram( dirname+'/phi' ).Fill(phi, pw)
        #sg.histogram( dirname+'/mu' ).Fill(avgmu, pw)
        #sg.histogram( dirname+'/nvtx' ).Fill(nvtx, pw)
        #sg.histogram( dirname+'/etVsEta' ).Fill(et,eta, pw)

  #
  # Finalize method
  #
  def finalize(self):
    
    basepath = self.getProperty( "Basepath" )
    sg = self.getStoreGateSvc()

    self.fina_lock()
    return StatusCode.SUCCESS 