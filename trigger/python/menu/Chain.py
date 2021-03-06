
__all__ = ['Chain']

from Gaugi import Algorithm
from Gaugi import StatusCode
from Gaugi import GeV
from Gaugi.messenger.macros import *
from EventAtlas import Accept
from TrigEgammaEmulationTool.menu import TriggerInfo


#
# Chain definition
#
class Chain( Algorithm ):


  #
  # Constructor
  #
  def __init__(self, name, L1Item, chain):

    Algorithm.__init__(self, name)
    self.__trigInfo = TriggerInfo(chain)
    self.__l1item = L1Item
    self.__trigger = chain



  #
  # Initialize method
  #
  def initialize(self):

    self.__signature = self.__trigInfo.signature()

    # Configure the L2Calo hypo step
    from TrigEgammaEmulationTool.TrigEgammaL1CaloHypoTool import configure
    self.__l1caloItem = configure( self.__l1item )


    if self.__trigInfo.signature() == 'electron':

      # Configure L2Calo step
      if self.__trigInfo.ringer():
        if self.__trigInfo.etthr() > 15: #GeV
          from TrigEgammaEmulationTool import installElectronRingerZeeFromVersion as configure
        else:
          from TrigEgammaEmulationTool import installElectronRingerJpsieeFromVersion as configure
        version = self.__trigInfo.ringerVersion()
        names = configure(version)
        # define like tight, medium, loose and vloose
        self.__l2caloItem = names[self.__trigInfo.pidnameIdx()]


        version = self.__trigInfo.ringerVersion("fast_el")
        if version is not None:
          names = configure( version, step="fast_el" )
          self.__l2Item = names[self.__trigInfo.pidnameIdx()]
        else:
          from TrigEgammaEmulationTool.TrigEgammaL2ElectronHypoTool import configure
          self.__l2Item = configure( self.__trigger )

      else:
        # Configure the L2Calo hypo step
        from TrigEgammaEmulationTool.TrigEgammaL2CaloHypoTool import configure
        self.__l2caloItem = configure( self.__trigger )

        # Configure the L2 hypo step
        from TrigEgammaEmulationTool.TrigEgammaL2ElectronHypoTool import configure
        self.__l2Item = configure( self.__trigger )



      # Configure the HLT hypo step
      from TrigEgammaEmulationTool.TrigEgammaElectronHypoTool import configure
      self.__hltItem = configure( self.__trigger )

      if self.__trigInfo.isolated():
        self.__applyIsolation=True
        from TrigEgammaEmulationTool.TrigEgammaElectronIsolationHypoTool import configure
        self.__hltIsoItem = configure( self.__trigger )
      else:
        self.__applyIsolation=False


    elif self.__trigInfo.signature() == 'photon':

      # Configure L2Calo step
      if self.__trigInfo.ringer():
        from TrigEgammaEmulationTool import installPhotonRingerZradFromVersion as configure
        version = self.__trigInfo.ringerVersion()
        names = configure(version)
        # define like tight, medium, loose and vloose
        self.__l2caloItem = names[self.__trigInfo.pidnameIdx()]

      else:
        # Configure the L2Calo hypo step
        from TrigEgammaEmulationTool.TrigEgammaL2CaloHypoTool import configure
        self.__l2caloItem = configure( self.__trigger )


      # Configure the EFCalo hypo step
      from TrigEgammaEmulationTool.TrigEgammaL2PhotonHypoTool import configure
      self.__l2Item = configure( self.__trigger )

      # Configure the HLT hypo step
      from TrigEgammaEmulationTool.TrigEgammaPhotonHypoTool import configure
      self.__hltItem = configure( self.__trigger )

      if self.__trigInfo.isolated():
        self.__applyIsolation=True
        from TrigEgammaEmulationTool.TrigEgammaPhotonIsolationHypoTool import configure
        self.__hltIsoItem = configure( self.__trigger )
      else:
        self.__applyIsolation=False

    else:
      MSG_FATAL(self, "Signature not found")


    # configure et cuts
    self.__l2caloEtCut = (self.__trigInfo.etthr() - 3 ) * GeV
    self.__efcaloEtCut = (self.__trigInfo.etthr()) * GeV
    self.__hltEtCut = (self.__trigInfo.etthr()) * GeV

    from Gaugi import ToolSvc
    emulator = ToolSvc.retrieve("Emulator")
    if emulator.retrieve( self.__l1caloItem ).initialize().isFailure():
      MSG_FATAL( self, "It's not possible to initialize the tool with name %s", self.__l1caloItem )
    if emulator.retrieve( self.__l2caloItem ).initialize().isFailure():
      MSG_FATAL( self, "It's not possible to initialize the tool with name %s", self.__l2caloItem )
    if emulator.retrieve( self.__l2Item ).initialize().isFailure():
      MSG_FATAL( self, "It's not possible to initialize the tool with name %s", self.__l2Item )
    if emulator.retrieve( self.__hltItem ).initialize().isFailure():
      MSG_FATAL( self, "It's not possible to initialize the tool with name %s", self.__hltItem )


    if self.__applyIsolation:
      if emulator.retrieve( self.__hltIsoItem ).initialize().isFailure():
        MSG_FATAL( self, "It's not possible to initialize the tool with name %s", self.__hltItem )


    # Print chain info steps
    MSG_INFO( self, "")
    MSG_INFO( self, "+ Chain with name   : %s", self.name() )
    MSG_INFO( self, "|--> L1Calo       : %s", self.__l1caloItem)
    MSG_INFO( self, "|--> L2Calo EtCut : %d", self.__l2caloEtCut)
    MSG_INFO( self, "|--> L2Calo       : %s", self.__l2caloItem)
    MSG_INFO( self, "|--> L2           : %s", self.__l2Item)
    MSG_INFO( self, "|--> EFCalo EtCut : %d", self.__efcaloEtCut)
    MSG_INFO( self, "|--> HLT EtCut    : %d", self.__hltEtCut)
    MSG_INFO( self, "|--> HLT          : %s", self.__hltItem)
    if self.__applyIsolation:
      MSG_INFO( self, "|--> HLTIso       : %s", self.__hltIsoItem)


    self.init_lock()
    return StatusCode.SUCCESS



  #
  # Finalize method
  #
  def finalize(self):
    self.fina_lock()
    return StatusCode.SUCCESS


  #
  # Accept method
  #
  def accept( self, context ):

    dec = context.getHandler( "MenuContainer" )

    accept = self.getAcceptInfo()

    passedL1 = bool(dec.accept( self.__l1caloItem ))


    # Is passed by L1?
    if not passedL1:
      return accept

    accept.setCutResult( 'L1Calo' , True )

    # Is passed by L2Calo et cut? AND hypo cut
    em = context.getHandler("HLT__TrigEMClusterContainer")

    if  not ( ( em.et() > self.__l2caloEtCut ) and bool(dec.accept( self.__l2caloItem )) ):
      return accept

    accept.setCutResult( 'L2Calo' , True )



    if self.__trigInfo.signature() == 'electron':
      cont = context.getHandler("HLT__TrigElectronContainer")
    else:
      cont = context.getHandler("HLT__PhotonContainer")

    # Is passed by L2 electron/photon, treat events with container with size equal zero
    passedL2 = bool(dec.accept( self.__l2Item )) if cont.size() > 0 else False

    if not passedL2:
      return accept


    accept.setCutResult( 'L2' , True )

    # Is passed by EF calo et cut
    clCont = context.getHandler( "HLT__CaloClusterContainer" )
    current = clCont.getPos()
    passedEFCalo = False
    for cl in clCont:
      if cl.et() > self.__efcaloEtCut:
        passedEFCalo=True
        break

    if not passedEFCalo:
      return accept


    accept.setCutResult( 'EFCalo' , True )



    # Is passed by HLT electron/photon et cut
    passedHLT_etcut = False

    if self.__signature == 'electron':
      cont = context.getHandler("HLT__ElectronContainer")
      for el in cont:
        if el.et() > self.__hltEtCut:
          passedHLT_etcut = True; break

    elif self.__signature == 'photon':
      cont = context.getHandler("HLT__PhotonContainer")
      for ph in cont:
        if ph.et() > self.__hltEtCut:
          passedHLT_etcut = True; break

    else:
      MSG_FATAL( self, "signature not reconized to emulate the HLT et cut step" )


    if not passedHLT_etcut:
      return accept


    # check the HLT decision
    passedHLT = bool( dec.accept( self.__hltItem ) )

    if passedHLT and self.__applyIsolation:
      # Apply the isolation cut and overwrite the HLT previus decision
      passedHLT = bool( dec.accept( self.__hltIsoItem ) )


    accept.setCutResult( 'HLT', passedHLT )
    accept.setCutResult( 'Pass', passedHLT )

    return accept



  def getAcceptInfo(self):

    accept = Accept( self.name() )
    accept.addCut( 'L1Calo' )
    accept.addCut( 'L2Calo' )
    accept.addCut( 'L2'     )
    accept.addCut( 'EFCalo' )
    accept.addCut( 'HLT'    )
    accept.addCut( 'Pass'   )
    return accept










