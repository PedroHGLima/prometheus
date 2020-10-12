
__all__ = ["Collector"]


from prometheus import Dataframe as DataframeEnum
from Gaugi import StatusCode, NotSet, retrieve_kw, progressbar
from Gaugi import csvStr2List, expandFolders, save, load
from Gaugi.messenger.macros import *
from Gaugi.messenger  import Logger
from Gaugi.constants import GeV
from Gaugi import EnumStringification
from Gaugi import Algorithm
import numpy as np


class Collector( Algorithm ):

  def __init__(self, name, **kw):
    
    Algorithm.__init__(self, name)
    self._event = {}
    self._event_label = []
    self._save_these_bins = list()
    self._extra_features = list()
    
    self.declareProperty( "OutputFile", 'sample.pic', "The output file name"       )
    self.declareProperty( "DoTrack"   , False        , "Dump all track variables." )


    for key, value in kw.items():
      self.setProperty(key, value)


  def setEtBinningValues( self, etbins ):
    self._etbins = etbins


  def setEtaBinningValues( self, etabins ):
    self._etabins = etabins


  def AddFeature( self, key ):
    self._extra_features.append( key )


  def initialize(self):
    Algorithm.initialize(self)
    for etBinIdx in range(len(self._etbins)-1):
      for etaBinIdx in range(len(self._etabins)-1):
        self._event[ 'et%d_eta%d' % (etBinIdx,etaBinIdx) ] = None

    doTrack = self.getProperty("DoTrack")

    self._event_label.append( 'avgmu' )

    # Add fast calo ringsE
    self._event_label.extend( [ 'L2Calo_ring_%d'%r for r in range(100) ] )

    self._event_label.extend( [ 'L2Calo_et',
                                'L2Calo_eta',
                                'L2Calo_phi',
                                'L2Calo_reta',
                                'L2Calo_ehad1', # new
                                'L2Calo_eratio',
                                'L2Calo_f1', # new
                                'L2Calo_f3', # new
                                'L2Calo_weta2', # new
                                'L2Calo_wstot', # new
                                'L2Calo_e2tsts1', # new
                                ] )


    if doTrack:
      self._event_label.extend( ['L2_hasTrack',
                                 'L2_pt',
                                 'L2_eta',
                                 'L2_phi',
                                 'L2_trkClusDeta',
                                 'L2_trkClusDphi',
                                 'L2_etOverPt'] )


    self._event_label.extend( [
                                # Offline variables
                                'et',
                                'eta',
                                'phi',
                                'eratio',
                                'reta',
                                'rphi',
                                'f1',
                                'f3',
                                'rhad',
                                'rhad1',
                                'wtots1',
                                'weta1',
                                'weta2',
                                'e277',
                                'deltaE',
                                'deltaR', # for boosted 
                                'eeMass', # for boosted
                                ] )



    if self._dataframe is DataframeEnum.Electron_v1:
      self._event_label.extend( [
                                # Offline variables
                                'el_lhtight',
                                'el_lhmedium',
                                'el_lhloose',
                                'el_lhvloose',
                                ] )
    elif self._dataframe is DataframeEnum.Photon_v1:
      self._event_label.extend( [
                                # Offline variables
                                'ph_tight',
                                'ph_medium',
                                'ph_loose',
                                ] )
    else:
      self._event_label.extend( [
                                # Offline variables
                                'el_lhtight',
                                'el_lhmedium',
                                'el_lhloose',
                                'el_lhvloose',
                                ] )


    self._event_label.extend( self._extra_features )

    return StatusCode.SUCCESS


  def fill( self, key , event ):

    if self._event[key]:
      self._event[key].append( event )
    else:
      self._event[key] = [event]


  #
  # execute 
  #
  def execute(self, context):

    doTrack = self.getProperty("DoTrack")
   

    if self._dataframe is DataframeEnum.Electron_v1:
      elCont    = context.getHandler( "ElectronContainer" )
    
    elif self._dataframe is DataframeEnum.Photon_v1:
      elCont    = context.getHandler( "PhotonContainer" )
    
    eventInfo = context.getHandler( "EventInfoContainer" )
    fc        = context.getHandler( "HLT__FastCaloContainer" )
    
   
    # For some reason, all trk object are the closes one. Maybe some bug into the ntuple.
    # But, here, we are interest to get the closest track object w.r.t the cluster. So,
    # You can use any trk position inside of the container.
    trkCont   = context.getHandler( "HLT__FastElectronContainer" )
    
    
    hasTrack = True if trkCont.size()>0 else False



    from PileupCorrectionTools.utilities import RetrieveBinningIdx
    etBinIdx, etaBinIdx = RetrieveBinningIdx( fc.et()/1000., abs(fc.eta()), self._etbins, self._etabins, logger=self._logger )
    if etBinIdx < 0 or etaBinIdx < 0:
      return StatusCode.SUCCESS


    key = ('et%d_eta%d') % (etBinIdx, etaBinIdx)

    event_row = list()
    # event info
    event_row.append( eventInfo.avgmu() )

    # fast calo features
    event_row.extend( fc.ringsE()   )
    event_row.append( fc.et()       )
    event_row.append( fc.eta()      )
    event_row.append( fc.phi()      )
    event_row.append( fc.reta()     )
    event_row.append( fc.ehad1()    )
    event_row.append( fc.eratio()   )
    event_row.append( fc.f1()       )
    event_row.append( fc.f3()       )
    event_row.append( fc.weta2()    )
    event_row.append( fc.wstot()    )
    event_row.append( fc.e2tsts1()  )



    # fast electron features
    if doTrack:

      if hasTrack:
        event_row.append( hasTrack)
        event_row.append( trkCont.pt() )
        event_row.append( trkCont.eta() )
        event_row.append( trkCont.phi() )
        event_row.append( trkCont.trkClusDeta() )
        event_row.append( trkCont.trkClusDphi() )
        event_row.append( trkCont.etOverPt() )
      else:
        event_row.extend( [False, -1, -1, -1, -1, -1, -1] )




    from EventAtlas import EgammaParameters
    
    #if fc.wstot() < 0:    
    #  wstot = elCont.showerShapeValue(EgammaParameters.wtots1)
    #  e2tsts1 = fc.e2tsts1() 
    #  eratio = elCont.showerShapeValue(EgammaParameters.Eratio)
    #  weta2 = elCont.showerShapeValue(EgammaParameters.weta2)
    #  print( 'reta = %1.2f, ehad1 = %1.2f, eratio = %1.2f (%1.2f), f1 = %1.2f, f3 = %1.2f, weta2 = %1.2f (%1.2f), wstot = %1.2f (%1.2f), e2tsts1= %1.2f' % 
    #           (fc.reta(),fc.ehad1(),fc.eratio(),eratio,fc.f1(),fc.f3(),fc.weta2(),weta2,fc.wstot(),wstot, e2tsts1) )
      
    # Offline Shower shapes
    event_row.append( elCont.et() )
    event_row.append( elCont.eta() )
    event_row.append( elCont.phi() )
    event_row.append( elCont.showerShapeValue( EgammaParameters.Eratio ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.Reta ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.Rphi ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.f1 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.f3 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.Rhad ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.Rhad1 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.wtots1 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.weta1 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.weta2 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.e277 ) )
    event_row.append( elCont.showerShapeValue( EgammaParameters.DeltaE ) )
    event_row.append( elCont.deltaR() )
    event_row.append( elCont.eeMass() )
    
    if self._dataframe is DataframeEnum.Electron_v1:
      event_row.append( elCont.accept( "el_lhtight"  ) )
      event_row.append( elCont.accept( "el_lhmedium" ) )
      event_row.append( elCont.accept( "el_lhloose"  ) )
      event_row.append( elCont.accept( "el_lhvloose" ) )


    elif self._dataframe is DataframeEnum.Photon_v1:
      event_row.append( elCont.accept( "ph_tight"  ) )
      event_row.append( elCont.accept( "ph_medium" ) )
      event_row.append( elCont.accept( "ph_loose"  ) )
      
    dec = context.getHandler("MenuContainer")

    for feature in self._extra_features:
      passed = dec.accept(feature).getCutResult('Pass')
      event_row.append( passed )


    self.fill(key , event_row)


    return StatusCode.SUCCESS


  def finalize( self ):

    from Gaugi import save, mkdir_p

    outputname = self.getProperty("OutputFile")

    for etBinIdx in range(len(self._etbins)-1):
      for etaBinIdx in range(len(self._etabins)-1):

        key =  'et%d_eta%d' % (etBinIdx,etaBinIdx)
        mkdir_p( outputname )
        if self._event[key] is None:
          continue

        d = {
            "features"  : self._event_label,
            "etBins"    : self._etbins,
            "etaBins"   : self._etabins,
            "etBinIdx"  : etBinIdx,
            "etaBinIdx" : etaBinIdx
            }

        d[ 'pattern_'+key ] = np.array( self._event[key] )
        MSG_INFO( self, 'Saving %s with : (%d, %d)', key, d['pattern_'+key].shape[0], d['pattern_'+key].shape[1] )
        save( d, outputname+'/'+outputname+"_"+key , protocol = 'savez_compressed')
    return StatusCode.SUCCESS





