
__all__ = ['EventSimulator']

from Gaugi.messenger import Logger, LoggingLevel
from Gaugi.messenger.macros import *
from prometheus.enumerations import Dataframe as DataframeEnum
from Gaugi import StatusCode
from Gaugi.gtypes import NotSet
from Gaugi import TEventLoop

# Import all root classes
import ROOT


# The main framework base class for SIM analysis.
# This class is responsible to build all containers object and
# manager the storegate and histogram services for all classes.
class EventSimulator( TEventLoop ):
    
  def __init__(self, name, **kw):
    # Retrieve all information needed
    TEventLoop.__init__(self, name, **kw)
    # Loading libraries
    ROOT.gSystem.Load('libprometheus')
    
  # Initialize all services
  def initialize( self ):

    MSG_INFO( self, 'Initializing EventSimulator...')
    
    if super(EventSimulator,self).initialize().isFailure():
      MSG_FATAL( self, "Impossible to initialize the TEventLoop services.")

   
    # RingerPhysVal hold the address of required branches
    if self._dataframe is DataframeEnum.Lorenzett_v1:
      #self._t.SetBranchStatus("*", False)
      from ROOT import edm
      from EventSimulator import CaloCellCollection, CaloRings, CaloCluster, EventInfo
      self._event = edm.Lorenzett_v1()
      self._t.GetEntry(0)
    else:
      return StatusCode.FATAL

    MSG_INFO( self, "Creating containers...")
    
    # Initialize the base of this container. 
    # Do not change this key names!
    # NOTE: Do not change this order.
    # we must retrieve the cells first and the reco other features
    # event dataframe containers
    self._containersSvc['EventInfoContainer']           = EventInfo()
    self._containersSvc['Truth__CaloCellsContainer']    = CaloCellCollection() 
    self._containersSvc['Truth__CaloRingsContainer']    = CaloRings()
    self._containersSvc['Truth__CaloClusterContainer']  = CaloCluster()
    self._containersSvc['CaloCellsContainer']           = CaloCellCollection() 
    self._containersSvc['CaloRingsContainer']           = CaloRings()
    self._containersSvc['CaloClusterContainer']         = CaloCluster()


                           

    # configure all EDMs needed
    for key, edm in self._containersSvc.items():

      # enable truth property by the container key name. 
      # Using hlt flag to avoid reimplementation
      edm.is_hlt = False if 'Truth__' in key else True


      # attach the EDM pointer into the context list
      self.getContext().setHandler(key,edm)
      # add properties
      edm.dataframe = self._dataframe
      edm.tree  = self._t
      edm.level = self._level
      edm.event = self._event
      edm.setContext(self.getContext())
      # If initializations is failed, we must remove this from the container 
      # service
      if(edm.initialize().isFailure()):
        MSG_WARNING( self, 'Impossible to create the EDM: %s',key)

    self.getContext().initialize()

    return StatusCode.SUCCESS


  def getEntry( self, entry ):
    self._t.GetEntry( entry )



