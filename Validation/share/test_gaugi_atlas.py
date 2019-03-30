
from Gaugi import EventATLASLoop
from Gaugi.enumerations import Dataframe as DataframeEnum
from Gaugi.messenger import LoggingLevel
from Gaugi import ToolSvc, ToolMgr
import argparse

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

import sys,os
if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()


ToolMgr += EventATLASLoop(  "EventATLASLoop",
                            inputFiles = args.inputFiles, 
                            treePath = '*/HLT/Physval/Egamma/probes', 
                            nov = args.nov,
                            dataframe = DataframeEnum.PhysVal_v2, 
                            outputFile = args.outputFile,
                            level = LoggingLevel.DEBUG
                          )


from Gaugi import Algorithm
ToolSvc += Algorithm( "AlgTest" )

from Gaugi import job
job.initialize()
job.execute()
job.finalize()

sys.exit(0)
