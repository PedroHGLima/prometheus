
__all__ = ['EfficiencyTool']


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


#
# Efficiency tool
#
class EfficiencyTool( Algorithm ):
  
  __triggerLevels = ['L1Calo','L2Calo','L2','EFCalo','HLT']

  #
  # Constructor
  #
  def __init__(self, name, dojpsiee=False, **kw):

    Algorithm.__init__(self, name)
    self.__groups = list()

    # declare all props here
    self.declareProperty( "Basepath", "Event/EfficiencyTool", "Histograms base path for the efficiency tool"      )
    self.declareProperty( "DoJpisee", dojpsiee                 , "Use the J/psiee et bins in the eff et histograms." )
 
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
  # Initialize method
  #
  def initialize(self):
    
    basepath = self.getProperty( "Basepath" )
    doJpsiee = self.getProperty( "DoJpisee" )

    sg = self.getStoreGateSvc()
    
    #et_bins  = zee_etbins
    eta_bins = default_etabins
    nvtx_bins.extend(high_nvtx_bins)
    #eta_bins = [0,0.6,0.8,1.15,1.37,1.52,1.81,2.01,2.37,2.47]
    et_bins = jpsiee_etbins if doJpsiee else [4.,7.,10.,15.,20.,25.,30.,35.,40.,45.,50.,60.,80.,150.,300.] 
    
    for group in self.__groups:
      # Get the chain object
      chain = group.chain()
      for dirname in ( self.__triggerLevels if (type(chain) is Chain or type(chain) is TDT) else ['Selector'] ):

        sg.mkdir( basepath+'/'+chain.name()+'/Efficiency/'+dirname )
        sg.addHistogram(TH1F('et','E_{T} distribution;E_{T};Count', len(et_bins)-1, np.array(et_bins)))
        sg.addHistogram(TH1F('eta','#eta distribution;#eta;Count', len(eta_bins)-1, np.array(eta_bins)))
        sg.addHistogram(TH1F("phi", "#phi distribution; #phi ; Count", 20, -3.2, 3.2));
        sg.addHistogram(TH1F('mu' ,'<#mu> distribution;<#mu>;Count', 20, 0, 100))
        sg.addHistogram(TH1F('nvtx' ,'N_{vtx} distribution;N_{vtx};Count', len(nvtx_bins)-1, np.array(nvtx_bins)))
        sg.addHistogram(TH1F('match_et','E_{T} matched distribution;E_{T};Count', len(et_bins)-1, np.array(et_bins)))
        sg.addHistogram(TH1F('match_eta','#eta matched distribution;#eta;Count', len(eta_bins)-1, np.array(eta_bins)))
        sg.addHistogram(TH1F("match_phi", "#phi matched distribution; #phi ; Count", 20, -3.2, 3.2));
        sg.addHistogram(TH1F('match_mu' ,'<#mu> matched distribution;<#mu>;Count', 20, 0, 100))
        sg.addHistogram(TH1F('match_nvtx' ,'N_{vtx} matched distribution;N_{vtx};Count', len(nvtx_bins)-1, np.array(nvtx_bins)))
        sg.addHistogram(TProfile("eff_et", "#epsilon(E_{T}); E_{T} ; Efficiency" , len(et_bins)-1, np.array(et_bins)))
        sg.addHistogram(TProfile("eff_eta", "#epsilon(#eta); #eta ; Efficiency"  , len(eta_bins)-1,np.array(eta_bins)))
        sg.addHistogram(TProfile("eff_phi", "#epsilon(#phi); #phi ; Efficiency", 20, -3.2, 3.2));
        sg.addHistogram(TProfile("eff_mu", "#epsilon(<#mu>); <#mu> ; Efficiency", 20, 0, 100));	
        sg.addHistogram(TProfile("eff_nvtx", "#epsilon(N_{vtx}); N_{vtx} ; Efficiency", len(nvtx_bins)-1, np.array(nvtx_bins)));	
        sg.addHistogram( TH2F('match_etVsEta', "Passed;E_{T};#eta;Count", len(et_bins)-1, np.array(et_bins), len(eta_bins)-1, np.array(eta_bins)) )
        sg.addHistogram( TH2F('etVsEta' , "Total;E_{T};#eta;Count", len(et_bins)-1,  np.array(et_bins), len(eta_bins)-1, np.array(eta_bins)) )
        sg.addHistogram( TProfile2D('eff_etVsEta' , "Total;E_{T};#eta;Count", len(et_bins)-1,  np.array(et_bins), len(eta_bins)-1, np.array(eta_bins)) )

        # for boosted analysis
        deltaR_bins = np.arange(0, 2.5, step=0.1)
        sg.addHistogram(TH1F('deltaR','#\Delta R distribution;#\Delta R;Count', len(deltaR_bins)-1, deltaR_bins))
        sg.addHistogram(TH1F('match_deltaR','#\Delta R matched distribution;#\Delta R;Count', len(deltaR_bins)-1, deltaR_bins))
        sg.addHistogram(TProfile('eff_deltaR', "#\epsilon(#\Delta R); #\Delta R ; Efficiency" , len(deltaR_bins)-1, deltaR_bins))

    self.init_lock()
    return StatusCode.SUCCESS 

    

  def execute(self, context):
  

    basepath = self.getProperty( "Basepath" )
    # Retrieve Electron container
    if self._dataframe is DataframeEnum.Electron_v1:
      elCont    = context.getHandler( "ElectronContainer" )
    elif self._dataframe is DataframeEnum.Photon_v1:
      elCont    = context.getHandler( "PhotonContainer" )
    else:
      elCont    = context.getHandler( "ElectronContainer" )
      
    dec = context.getHandler( "MenuContainer" )


    for group in self.__groups:
      chain = group.chain()
      if type(chain) is Chain or type(chain) is TDT:
        for el in elCont:
          if el.et()  < (group.etthr()- 5)*GeV:  continue 
          if abs(el.eta())>2.47: continue
          dirname = basepath+'/'+chain.name()+'/Efficiency'
          accept = dec.accept( chain.name() )
          self.fillEfficiency(dirname+'/'+'L1Calo', el, group.etthr(), group.pidname(), accept.getCutResult("L1Calo") )
          self.fillEfficiency(dirname+'/'+'L2Calo', el, group.etthr(), group.pidname(), accept.getCutResult("L2Calo") )
          self.fillEfficiency(dirname+'/'+'L2'    , el, group.etthr(), group.pidname(), accept.getCutResult("L2")     )
          self.fillEfficiency(dirname+'/'+'EFCalo', el, group.etthr(), group.pidname(), accept.getCutResult("EFCalo") )
          self.fillEfficiency(dirname+'/'+'HLT'   , el, group.etthr(), group.pidname(), accept.getCutResult("HLT")    )  

      else:
        MSG_FATAL( self, "Chain type not reconized." )

    return StatusCode.SUCCESS 



  #
  # Fill efficiency histograms
  #
  def fillEfficiency( self, dirname, el, etthr, pidword, isPassed ):
  
    sg = self.getStoreGateSvc()
    
    pid = el.accept(pidword) if pidword else True

    eta = el.caloCluster().etaBE2()
    et = el.et()/GeV
    phi = el.phi()
    evt = self.getContext().getHandler("EventInfoContainer")
    avgmu = evt.avgmu()
    nvtx = evt.nvtx()
    pw = evt.MCPileupWeight()

    deltaR = el.deltaR()
    if pid: 
      sg.histogram( dirname+'/et' ).Fill(et, pw)
      if et > etthr+1.0:
        sg.histogram( dirname+'/eta' ).Fill(eta, pw)
        sg.histogram( dirname+'/phi' ).Fill(phi, pw)
        sg.histogram( dirname+'/mu' ).Fill(avgmu, pw)
        sg.histogram( dirname+'/nvtx' ).Fill(nvtx, pw)
        sg.histogram( dirname+'/etVsEta' ).Fill(et,eta, pw)

        sg.histogram( dirname+'/deltaR' ).Fill(deltaR, pw)

      if isPassed:
        sg.histogram( dirname+'/match_et' ).Fill(et, pw)
        sg.histogram( dirname+'/eff_et' ).Fill(et,1, pw)
        
        if et > etthr+1.0:
          sg.histogram( dirname+'/match_eta' ).Fill(eta, pw)
          sg.histogram( dirname+'/match_phi' ).Fill(phi, pw)
          sg.histogram( dirname+'/match_mu' ).Fill(avgmu, pw)
          sg.histogram( dirname+'/match_nvtx' ).Fill(nvtx, pw)
          
          sg.histogram( dirname+'/match_deltaR' ).Fill(deltaR, pw)
          
          sg.histogram( dirname+'/match_etVsEta' ).Fill(et,eta, pw)
          sg.histogram( dirname+'/eff_eta' ).Fill(eta,1, pw)
          sg.histogram( dirname+'/eff_phi' ).Fill(phi,1, pw)
          sg.histogram( dirname+'/eff_mu' ).Fill(avgmu,1, pw)
          sg.histogram( dirname+'/eff_nvtx' ).Fill(nvtx,1, pw)
          sg.histogram( dirname+'/eff_etVsEta' ).Fill(et,eta,1, pw)

          sg.histogram( dirname+'/eff_deltaR' ).Fill(deltaR, 1, pw)

      else:
        sg.histogram( dirname+'/eff_et' ).Fill(et,0)
        if et > etthr+1.0:
          sg.histogram( dirname+'/eff_eta' ).Fill(eta,0, pw)
          sg.histogram( dirname+'/eff_phi' ).Fill(phi,0, pw)
          sg.histogram( dirname+'/eff_mu' ).Fill(avgmu,0, pw)
          sg.histogram( dirname+'/eff_nvtx' ).Fill(nvtx,0, pw)
          sg.histogram( dirname+'/eff_etVsEta' ).Fill(et,eta,0, pw)
          # for boosted 
          sg.histogram( dirname+'/eff_deltaR' ).Fill(deltaR,0, pw)



  #
  # Finalize method
  #
  def finalize(self):
    
    basepath = self.getProperty( "Basepath" )
    sg = self.getStoreGateSvc()
    for group in self.__groups:
      chain = group.chain()
      MSG_INFO( self, '{:-^78}'.format((' %s ')%(chain.name())))
      if type(chain) is Chain or type(chain) is TDT:
        for trigLevel in self.__triggerLevels:
          dirname = basepath+'/'+chain.name()+'/Efficiency/'+trigLevel
          total  = sg.histogram( dirname+'/eta' ).GetEntries()
          passed = sg.histogram( dirname+'/match_eta' ).GetEntries()
          eff = passed/float(total) * 100. if total>0 else 0
          eff=('%1.2f')%(eff); passed=('%d')%(passed); total=('%d')%(total)
          stroutput = '| {0:<30} | {1:<5} ({2:<5}, {3:<5}) |'.format(trigLevel,eff,passed,total)
          MSG_INFO( self, stroutput)
      else:
        pass
        
      MSG_INFO(self, '{:-^78}'.format((' %s ')%('-')))
    self.fina_lock()
    return StatusCode.SUCCESS 




