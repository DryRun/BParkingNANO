import os
import sys
import glob
import numpy as np
import math
import datetime

redo_splitting = True
subjob_filelist_paths = []
if redo_splitting:
	input_files = sorted(glob.glob("/home/dryu/store/BParkingMC/BdToKstarJpsi_ToKPiMuMu_probefilter_part*/*root"))
	files_per_job = 20
	nsubjobs = int(math.ceil(1. * len(input_files) / files_per_job))
	subjob_filelists = np.array_split(np.array(input_files), nsubjobs)
	for isubjob, subjob_files in enumerate(subjob_filelists):
		subjob_filelist_path = "subjobfiles_{}.txt".format(isubjob)
		with open(subjob_filelist_path, "w") as subjob_filelist:
			for subjob_file in subjob_files:
				subjob_filelist.write("file:{}\n".format(subjob_file))
		subjob_filelist_paths.append(os.path.abspath(subjob_filelist_path))
else:
	subjob_filelist_paths = glob.glob("/user_data/dryu/BFrag/CMSSW_10_2_15_skim/src/PhysicsTools/BParkingNano/test/Bu2PiJPsi2PiMuMu/subjobfiles*txt")

ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
submitdir = os.path.expandvars("$HOME/BFrag/CMSSW_10_2_15_skim/src/PhysicsTools/BParkingNano/test/Bu2PiJPsi2PiMuMu/condor/job{}".format(ts))
os.system("mkdir -pv {}".format(submitdir))

with open("{}/run_Bu2PiJpsiNano.sh".format(submitdir), "w") as run_script:
	run_script.write("#!/bin/bash\n")
	run_script.write("source /cvmfs/cms.cern.ch/cmsset_default.sh\n")
	run_script.write("FILELISTS=({})\n".format(" ".join(subjob_filelist_paths)))
	run_script.write("tar -xzf CMSSW_10_2_15_skim.tar.gz\n")
	run_script.write("cd CMSSW_10_2_15_skim/src\n")
	run_script.write("eval `scram runtime -sh`\n")
	run_script.write("cd $_CONDOR_SCRATCH_DIR\n")
	run_script.write("echo DEBUG\n")
	run_script.write("ls -lrth\n")
	run_script.write("cmsRun run_nano_BuKMuMu_cfg.py inputFiles_load=${FILELISTS[$1]} isMC=True tag=\"subjob$1\"\n")
files_to_transfer = ["{}/run_Bu2PiJpsiNano.sh".format(submitdir), "$HOME/BFrag/CMSSW_10_2_15_skim/src/PhysicsTools/BParkingNano/test/run_nano_BuKMuMu_cfg.py", "$HOME/BFrag/CMSSW_10_2_15_skim.tar.gz"]
csub_command = "csub run_Bu2PiJpsiNano.sh -F {} -n {} -d {}".format(
	",".join(files_to_transfer),
	len(subjob_filelist_paths),
	submitdir
	)
print(csub_command)
os.system(csub_command)
