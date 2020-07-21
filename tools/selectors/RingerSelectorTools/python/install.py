

__all__ =  [
            #"installElectronL2CaloRingerSelector_v5", 
            "installElectronL2CaloRingerSelector_v6",
            "installElectronL2CaloRingerSelector_v8",
            "installElectronL2CaloRingerSelector_v10",
           ]
import os



# same as ringer v6 but use the output after the tansig TF function in 
# the last neuron
#def installElectronL2CaloRingerSelector_v5( toolname = "Emulator" ):
#
#  from RingerSelectorTools import RingerSelectorTool
#  # do not change this paths...
#  #calibpath = 'RingerSelectorTools/TrigL2_20170505_v6'
#  calibpath = os.environ['PRT_PATH'] + '/tools/Selectors/RingerSelectorTools/data/TrigL2_20170505_v6'
#
#  selectors = [
#      RingerSelectorTool( "T0HLTElectronRingerTight_v5", 
#                          calibpath+'/TrigL2CaloRingerElectronTightConstants.json', 
#                          calibpath+'/TrigL2CaloRingerElectronTightThresholds.json',
#      RingerSelectorTool( "T0HLTElectronRingerMedium_v5", 
#                          calibpath+'/TrigL2CaloRingerElectronMediumConstants.json', 
#                          calibpath+'/TrigL2CaloRingerElectronMediumThresholds.json', 
#                          remove_last_activation=False ), 
#      RingerSelectorTool( "T0HLTElectronRingerLoose_v5", 
#                          calibpath+'/TrigL2CaloRingerElectronLooseConstants.json', 
#                          calibpath+'/TrigL2CaloRingerElectronLooseThresholds.json', 
#                          remove_last_activation=False ), 
#      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v5", 
#                          calibpath+'/TrigL2CaloRingerElectronVeryLooseConstants.json', 
#                          calibpath+'/TrigL2CaloRingerElectronVeryLooseThresholds.json', 
#                          remove_last_activation=False ), 
#
#    ]
#
#  from Gaugi import ToolSvc as toolSvc
#  tool = toolSvc.retrieve( toolname )
#  if tool:
#    for sel in selectors:
#      tool+=sel
#  else:
#    raise RuntimeError( "%s not found into the ToolSvc." % toolname )



###########################################################
################## Official 2017 tuning ###################
###########################################################
def installElectronL2CaloRingerSelector_v6( toolname = "Emulator" ):

  from RingerSelectorTools import RingerSelectorTool
  # do not change this paths...
  #calibpath = 'RingerSelectorTools/TrigL2_20180125_v8'
  calibpath = os.environ['PRT_PATH'] + '/tools/Selectors/RingerSelectorTools/data/TrigL2_20170505_v6'

  selectors = [
      RingerSelectorTool( "T0HLTElectronRingerTight_v6",
                          calibpath+'/ElectronRingerTightTriggerConfig.conf'), 
      RingerSelectorTool( "T0HLTElectronRingerMedium_v6", 
                          calibpath+'/ElectronRingerMediumTriggerConfig.conf'), 
      RingerSelectorTool( "T0HLTElectronRingerLoose_v6", 
                          calibpath+'/ElectronRingerLooseTriggerConfig.conf'), 
      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v6", 
                          calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf'), 

    ]

  from Gaugi import ToolSvc as toolSvc
  tool = toolSvc.retrieve( toolname )
  if tool:
    for sel in selectors:
      tool+=sel
  else:
    raise RuntimeError( "%s not found into the ToolSvc." % toolname )




###########################################################
################## Official 2018 tuning ###################
###########################################################
def installElectronL2CaloRingerSelector_v8( toolname = "Emulator" ):

  from RingerSelectorTools import RingerSelectorTool
  # do not change this paths...
  #calibpath = 'RingerSelectorTools/TrigL2_20180125_v8'
  calibpath = os.environ['PRT_PATH'] + '/tools/Selectors/RingerSelectorTools/data/TrigL2_20180125_v8'

  selectors = [
      RingerSelectorTool( "T0HLTElectronRingerTight_v8",
                          calibpath+'/ElectronRingerTightTriggerConfig.conf'), 
      RingerSelectorTool( "T0HLTElectronRingerMedium_v8", 
                          calibpath+'/ElectronRingerMediumTriggerConfig.conf'), 
      RingerSelectorTool( "T0HLTElectronRingerLoose_v8", 
                          calibpath+'/ElectronRingerLooseTriggerConfig.conf'), 
      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v8", 
                          calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf'), 

    ]

  from Gaugi import ToolSvc as toolSvc
  tool = toolSvc.retrieve( toolname )
  if tool:
    for sel in selectors:
      tool+=sel
  else:
    raise RuntimeError( "%s not found into the ToolSvc." % toolname )



  
###########################################################
################## Testing 2020 tuning  ###################
###########################################################
def installElectronL2CaloRingerSelector_v10( toolname = "Emulator" ):

  from RingerSelectorTools import RingerSelectorTool
  # do not change this paths...
  #calibpath = 'RingerSelectorTools/TrigL2_20180125_v8'
  calibpath = os.environ['PRT_PATH'] + '/tools/Selectors/RingerSelectorTools/data/TrigL2_20200715_v10'

  
  def norm1_and_reshape( data ):
      return (data/abs(sum(data))).reshape((1,100, 1))


  selectors = [
      RingerSelectorTool( "T0HLTElectronRingerTight_v10",
                          calibpath+'/ElectronRingerTightTriggerConfig.conf',
                          norm1_and_reshape), 
      RingerSelectorTool( "T0HLTElectronRingerMedium_v10", 
                          calibpath+'/ElectronRingerMediumTriggerConfig.conf', 
                          norm1_and_reshape), 
      RingerSelectorTool( "T0HLTElectronRingerLoose_v10", 
                          calibpath+'/ElectronRingerLooseTriggerConfig.conf', 
                          norm1_and_reshape), 
      RingerSelectorTool( "T0HLTElectronRingerVeryLoose_v10", 
                          calibpath+'/ElectronRingerVeryLooseTriggerConfig.conf', 
                          norm1_and_reshape), 

    ]

  from Gaugi import ToolSvc as toolSvc
  tool = toolSvc.retrieve( toolname )
  if tool:
    for sel in selectors:
      tool+=sel
  else:
    raise RuntimeError( "%s not found into the ToolSvc." % toolname )





