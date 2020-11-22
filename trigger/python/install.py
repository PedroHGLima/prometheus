

__all__ =  [
            # Old T2Calo selector
            "installTrigEgammaL2CaloSelectors",
            "installTrigEgammaL2ElectronSelectors",

            # Zee for electrons signatures
            "installElectronL2CaloRingerSelector_v6",
            "installElectronL2CaloRingerSelector_v8",
            "installElectronL2CaloRingerSelector_v10",
            
            #jpsiee for electron signatures
            'installLowEnergyElectronL2CaloRingerSelector_v1',

            #zrad for photon signatures
            'installPhotonL2CaloRingerSelector_v1',


            # install helper
            'installElectronRingerZeeFromVersion',
            'installElectronRingerJpsieeFromVersion',
            'installPhotonRingerZradFromVersion',
           ]






def installElectronRingerZeeFromVersion( key , useOnnx=False, step="fast_calo"):
  
  versions =  {
                "fast_calo" : {
                  # Zee
                  "v6"                 : installElectronL2CaloRingerSelector_v6(useOnnx),
                  "v8"                 : installElectronL2CaloRingerSelector_v8(useOnnx),
                  "v10"                : installElectronL2CaloRingerSelector_v10(useOnnx),
                  "v11"                : installElectronL2CaloRingerSelector_v11(useOnnx),
                },
                
                "fast_el" : {
                  #"v1_el"                 : installElectronL2RingerSelector_v1_el(useOnnx),
                  }

             }

  return versions[step][key]




def installPhotonRingerZradFromVersion( key , useOnnx=True, step="fast_calo"):
  
  versions =  {
               "fast_calo" : {
                'v1'                 : installPhotonL2CaloRingerSelector_v1(useOnnx)
                },
               "fast_ph" : {

                }
              }
  return versions[step][key]






def installElectronRingerJpsieeFromVersion( key , useOnnx=True, step="fast_calo" ):
  
  versions =  {
                "fast_calo" : {
                  # Jpsiee
                  'v1'                 : installLowEnergyElectronL2CaloRingerSelector_v1(useOnnx),
                },
                "fast_el" : {
                
                }
              }

  return versions[step][key]




#
# Helper to avoid to much repetition code into this file
#
def attach( hypos ):
  from Gaugi import ToolSvc
  emulator = ToolSvc.retrieve( "Emulator" )
  names = []
  for hypo in hypos:
    names.append( hypo.name() )
    if not emulator.isValid( hypo.name() ):
      emulator+=hypo
  return names






def installTrigEgammaL2CaloSelectors():

  from TrigEgammaEmulationTool import TrigEgammaL2CaloSelectorTool
  hypos = [
        TrigEgammaL2CaloSelectorTool("T0HLTElectronT2CaloTight"   , OperationPoint ='lhtight'  ) ,
        TrigEgammaL2CaloSelectorTool("T0HLTElectronT2CaloMedium"  , OperationPoint ='lhmedium' ) ,
        TrigEgammaL2CaloSelectorTool("T0HLTElectronT2CaloLoose"   , OperationPoint ='lhloose'  ) ,
        TrigEgammaL2CaloSelectorTool("T0HLTElectronT2CaloVLoose"  , OperationPoint ='lhvloose' ) ,
      ]

  return attach(hypos)




def installTrigEgammaL2ElectronSelectors():

  from TrigEgammaEmulationTool import TrigEgammaL2ElectronSelectorTool
  hypos = [
        TrigEgammaL2ElectronSelectorTool("T0HLTElectronL2") 
      ]
  return attach(hypos)





###########################################################
################## Official 2017 tuning ###################
###########################################################
def installElectronL2CaloRingerSelector_v6( useOnnx=False ):

  from TrigEgammaEmulationTool import RingerSelectorTool
  import os
  calibpath = os.environ['PRT_PATH'] + '/trigger/data/zee/TrigL2_20170505_v6'
 
  def getPatterns( context ):
    def norm1( data ):
      return (data/abs(sum(data))).reshape((1,100))
    fc = context.getHandler("HLT__FastCaloContainer")
    rings = norm1( fc.ringsE() )
    return [rings]

 
  hypos = [
      RingerSelectorTool("T0HLTElectronRingerTight_v6"    ,getPatterns,ConfigFile = calibpath+'/ElectronRingerTightTriggerConfig.conf'    ,UseOnnx=useOnnx), 
      RingerSelectorTool("T0HLTElectronRingerMedium_v6"   ,getPatterns,ConfigFile = calibpath+'/ElectronRingerMediumTriggerConfig.conf'   ,UseOnnx=useOnnx), 
      RingerSelectorTool("T0HLTElectronRingerLoose_v6"    ,getPatterns,ConfigFile = calibpath+'/ElectronRingerLooseTriggerConfig.conf'    ,UseOnnx=useOnnx), 
      RingerSelectorTool("T0HLTElectronRingerVeryLoose_v6",getPatterns,ConfigFile = calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf',UseOnnx=useOnnx), 
    ]
  return attach(hypos)



###########################################################
################## Official 2018 tuning ###################
###########################################################
def installElectronL2CaloRingerSelector_v8( useOnnx=False ):

  from TrigEgammaEmulationTool import RingerSelectorTool
  import os
  calibpath = os.environ['PRT_PATH'] + '/trigger/data/zee/TrigL2_20180125_v8'

  def getPatterns( context ):
    def norm1( data ):
      return (data/abs(sum(data))).reshape((1,100))
    fc = context.getHandler("HLT__FastCaloContainer")
    rings = norm1( fc.ringsE() )
    return [rings]


  hypos = [
      RingerSelectorTool( "T0HLTElectronRingerTight_v8"    ,getPatterns, ConfigFile = calibpath+'/ElectronRingerTightTriggerConfig.conf'     , UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerMedium_v8"   ,getPatterns, ConfigFile = calibpath+'/ElectronRingerMediumTriggerConfig.conf'    , UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerLoose_v8"    ,getPatterns, ConfigFile = calibpath+'/ElectronRingerLooseTriggerConfig.conf'     , UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v8",getPatterns, ConfigFile = calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf' , UseOnnx=useOnnx), 
    ]

  return attach(hypos)



  
###########################################################
################## Testing 2020 tuning  ###################
###########################################################
def installElectronL2CaloRingerSelector_v10( useOnnx=False ):

  from TrigEgammaEmulationTool import RingerSelectorTool
  import os
  calibpath = os.environ['PRT_PATH'] + '/trigger/data/zee/TrigL2_20200715_v10'

  def getPatterns( context ):
    def norm1( data ):
      return (data/abs(sum(data))).reshape((1,100))
    fc = context.getHandler("HLT__FastCaloContainer")
    rings = norm1( fc.ringsE() )
    return [rings]


  hypos = [
      RingerSelectorTool( "T0HLTElectronRingerTight_v10"    , getPatterns, ConfigFile = calibpath+'/ElectronRingerTightTriggerConfig.conf',     
          UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerMedium_v10"   , getPatterns, ConfigFile = calibpath+'/ElectronRingerMediumTriggerConfig.conf',    
          UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerLoose_v10"    , getPatterns, ConfigFile = calibpath+'/ElectronRingerLooseTriggerConfig.conf',     
          UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v10", getPatterns, ConfigFile = calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf', 
          UseOnnx=useOnnx), 
    ]

  return attach(hypos)





###########################################################
################## Testing 2020 tuning  ###################
###########################################################
def installElectronL2CaloRingerSelector_v11( useOnnx=False ):

  # Using shower shapes + rings here

  from TrigEgammaEmulationTool import RingerSelectorTool
  import os
  calibpath = os.environ['PRT_PATH'] + '/trigger/data/zee/TrigL2_20200715_v11'


  def getPatterns( context ):
    def norm1( data ):
      return (data/abs(sum(data))).reshape((1,100))
    fc = context.getHandler("HLT__FastCaloContainer")
    rings = norm1( fc.ringsE() )
    reta = fc.reta()
    eratio = fc.eratio()
    f1 = fc.f1()/0.6
    f3 = fc.f3()/0.04
    weta2 =fc.weta2()/0.02
    wstot = fc.wstot()
    if eratio>10.0:
      eratio = 0.0
    elif eratio>1.0:
      eratio=1.0
    if wstot<-99:
      wstot=0.0

    return [rings, [reta,eratio,f1,f3,weta2,wstot]]


  hypos = [
      RingerSelectorTool( "T0HLTElectronRingerTight_v11"    , getPatterns, ConfigFile = calibpath+'/ElectronRingerTightTriggerConfig.conf'     , UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerMedium_v11"   , getPatterns, ConfigFile = calibpath+'/ElectronRingerMediumTriggerConfig.conf'    , UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerLoose_v11"    , getPatterns, ConfigFile = calibpath+'/ElectronRingerLooseTriggerConfig.conf'     , UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v11", getPatterns, ConfigFile = calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf' , UseOnnx=useOnnx), 
    ]

  return attach(hypos)



###########################################################
###################  ZRad v1 tuning   #####################
###########################################################
def installPhotonL2CaloRingerSelector_v1( useOnnx=True ):
  '''
  This tuning is the very medium tuning which was adjusted to operate in the knee of the ROC curve given the best balance between PD and FR.
  '''
  from TrigEgammaEmulationTool import RingerSelectorTool
  import os
  calibpath = os.environ['PRT_PATH'] + '/trigger/data/zrad/TrigL2_20200909_v1'

  def getPatterns( context ):
    def norm1( data ):
      return (data/abs(sum(data))).reshape((1,100))
    fc = context.getHandler("HLT__FastCaloContainer")
    rings = norm1( fc.ringsE() )
    return [rings]


  hypos = [
              RingerSelectorTool( "T0HLTPhotonRingerTight_v1" ,getPatterns  , ConfigFile = calibpath+'/PhotonRingerTightTriggerConfig.conf' , UseOnnx=useOnnx), 
              RingerSelectorTool( "T0HLTPhotonRingerMedium_v1",getPatterns  , ConfigFile = calibpath+'/PhotonRingerMediumTriggerConfig.conf', UseOnnx=useOnnx), 
              RingerSelectorTool( "T0HLTPhotonRingerLoose_v1" ,getPatterns  , ConfigFile = calibpath+'/PhotonRingerLooseTriggerConfig.conf' , UseOnnx=useOnnx), 
    ]

  return attach(hypos)





###########################################################
################### jpsiee v1 tuning  #####################
###########################################################
def installLowEnergyElectronL2CaloRingerSelector_v1( useOnnx=False ):

  from TrigEgammaEmulationTool import RingerSelectorTool
  import os
  calibpath = os.environ['PRT_PATH'] + '/trigger/data/jpsiee/TrigL2_20200805_v1'


  def getPatterns( context ):
    def norm1( data ):
      return (data/abs(sum(data))).reshape((1,100))
    fc = context.getHandler("HLT__FastCaloContainer")
    rings = norm1( fc.ringsE() )
    return [rings]


  hypos = [
      RingerSelectorTool( "T0HLTLowEnergyElectronRingerTight_v1"    ,getPatterns, ConfigFile=calibpath+'/ElectronRingerTightTriggerConfig.conf'    ,UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTLowEnergyElectronRingerMedium_v1"   ,getPatterns, ConfigFile=calibpath+'/ElectronRingerMediumTriggerConfig.conf'   ,UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTLowEnergyElectronRingerLoose_v1"    ,getPatterns, ConfigFile=calibpath+'/ElectronRingerLooseTriggerConfig.conf'    ,UseOnnx=useOnnx), 
      RingerSelectorTool( "T0HLTLowEnergyElectronRingerVeryLoose_v1",getPatterns, ConfigFile=calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf',UseOnnx=useOnnx), 
    ]

  return attach(hypos)






