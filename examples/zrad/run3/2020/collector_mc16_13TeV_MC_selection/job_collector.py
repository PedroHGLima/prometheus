from prometheus import EventATLAS
from prometheus.enumerations import Dataframe as DataframeEnum
from Gaugi.messenger import LoggingLevel, Logger
from Gaugi import ToolSvc, ToolMgr
import argparse
mainLogger = Logger.getModuleLogger("job")
parser = argparse.ArgumentParser(description = '', add_help = False)
parser = argparse.ArgumentParser()


parser.add_argument('-i','--inputFiles', action='store',
    dest='inputFiles', required = True, nargs='+',
    help = "The input files that will be used to generate the plots")

parser.add_argument('-o','--outputFile', action='store',
    dest='outputFile', required = False, default = None,
    help = "The output store name.")

parser.add_argument('-n','--nov', action='store',
    dest='nov', required = False, default = -1, type=int,
    help = "Number of events.")

parser.add_argument('--Zee', action='store_true',
    dest='doZee', required = False,
    help = "Do Zee collection.")

parser.add_argument('--Jpsi', action='store_true',
    dest='doJpsi', required = False,
    help = "Do Jpsi collection.")

parser.add_argument('--Zrad', action='store_true',
    dest='doZrad', required = False,
    help = "Do Zrad collection.")

parser.add_argument('--egam7', action='store_true',
    dest='doEgam7', required = False,
    help = "The colelcted sample came from EGAM7 skemma.")


parser.add_argument('--fakes', action='store_true',
    dest='doFakes', required = False,
    help = "The colelcted sample came from EGAM7 skemma.")

import sys,os
if len(sys.argv)==1:
  parser.print_help()
  sys.exit(1)

args = parser.parse_args()


if args.doZee or args.doJpsi:
  signature = 'electron'
elif args.doZrad or args.doFakes:
  signature = 'photon'
else:
  signature = 'electron'

acc = EventATLAS(  "EventATLASLoop",
                  inputFiles = args.inputFiles,
                  treePath= '*/HLT/PhysVal/Egamma/fakes' if args.doFakes else '*/HLT/PhysVal/Egamma/photons',
                  dataframe = DataframeEnum.Electron_v1 if signature == 'electron' else DataframeEnum.Photon_v1,
                  outputFile = args.outputFile,
                  # outputFile = 'dummy.root',
                  level = LoggingLevel.INFO,
                )



from EventSelectionTool import EventSelection, SelectionType, EtCutType

evt = EventSelection('EventSelection')
# evt.setCutValue( SelectionType.SelectionOnlineWithRings )

# Do not change this!
if args.doEgam7:
  #pidname = '!VeryLooseLLH_DataDriven_Rel21_Run2_2018'
  pidname = '!el_lhvloose'
elif args.doZee or args.doJpsi:
  #pidname = 'MediumLLH_DataDriven_Rel21_Run2_2018'
  pidname = 'el_lhmedium'
elif args.doZrad:
  pidname = 'ph_medium'
elif args.doFakes:
  pidname = '!ph_loose'
else:
  pidname = 'el_lhmedium'

if args.doZee:
    evt.setCutValue( SelectionType.SelectionPID, pidname )
    evt.setCutValue( EtCutType.L2CaloAbove , 15)
    ToolSvc += evt
elif args.doJpsi:
    evt.setCutValue( SelectionType.SelectionPID, pidname )
    evt.setCutValue( EtCutType.L2CaloAbove, 4 )
    evt.setCutValue( EtCutType.L2CaloBelow, 15 )
    evt.setCutValue( EtCutType.OfflineAbove, 2 )
    ToolSvc += evt
elif args.doZrad:
    evt.setCutValue( SelectionType.SelectionPID, pidname )
    evt.setCutValue( EtCutType.L2CaloAbove , 15)
    ToolSvc += evt
elif args.doFakes:
    evt.setCutValue( SelectionType.SelectionPID, pidname )
    evt.setCutValue( EtCutType.L2CaloAbove , 15)
    ToolSvc += evt
else:
    evt.setCutValue( SelectionType.SelectionPID, pidname )
    evt.setCutValue( EtCutType.L2CaloAbove , 15)
    ToolSvc += evt

if args.doZee:
#    from RingerSelectorTools import installElectronL2CaloRingerSelector_v6
#    installElectronL2CaloRingerSelector_v6()
    from RingerSelectorTools import installElectronL2CaloRingerSelector_v8
    installElectronL2CaloRingerSelector_v8()

from TrigEgammaEmulationTool import installTrigEgammaL2CaloSelectors
installTrigEgammaL2CaloSelectors()


from CollectorTool import Collector
alg = Collector( 'Collector' , OutputFile = args.outputFile.replace('.root',''), 
                )

etabins = [0.0, 0.8, 1.37, 1.54, 2.37, 2.50]
if args.doZee:
    etbins  = [15.0, 20.0, 30.0, 40.0, 50.0, 1000000.0]
elif args.doJpsi:
    etbins = [0.0, 7.0, 10.0, 15.0]
else:
    etbins  = [15.0, 20.0, 30.0, 40.0, 50.0, 1000000.0]
alg.setEtBinningValues( etbins   )
alg.setEtaBinningValues( etabins )
alg.doTrigger  = True


alg.AddFeature( "T0HLTElectronT2CaloTight"        )
alg.AddFeature( "T0HLTElectronT2CaloMedium"       )
alg.AddFeature( "T0HLTElectronT2CaloLoose"        )
alg.AddFeature( "T0HLTElectronT2CaloVLoose"       )
if args.doZee:
#    alg.AddFeature( "T0HLTElectronRingerTight_v6"     )
#    alg.AddFeature( "T0HLTElectronRingerMedium_v6"    )
#    alg.AddFeature( "T0HLTElectronRingerLoose_v6"     )
#    alg.AddFeature( "T0HLTElectronRingerVeryLoose_v6" )
    alg.AddFeature( "T0HLTElectronRingerTight_v8"     )
    alg.AddFeature( "T0HLTElectronRingerMedium_v8"    )
    alg.AddFeature( "T0HLTElectronRingerLoose_v8"     )
    alg.AddFeature( "T0HLTElectronRingerVeryLoose_v8" )
#alg.AddFeature( "HLT__isLHTight"                  )
#alg.AddFeature( "HLT__isLHMedium"                 )
#alg.AddFeature( "HLT__isLHLoose"                  )
#alg.AddFeature( "HLT__isLHVLoose"                 )
ToolSvc += alg

acc.run(args.nov)
