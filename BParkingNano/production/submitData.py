from CRABClient.UserUtilities import config, ClientException, getUsernameFromCRIC
from PhysicsTools.BParkingNano.skim_version import skim_version
#from input_crab_data import dataset_files
import yaml
import datetime
from fnmatch import fnmatch
from argparse import ArgumentParser

production_tag = datetime.date.today().strftime('%Y%b%d')

if __name__ == '__main__':

  from CRABAPI.RawCommand import crabCommand
  from CRABClient.ClientExceptions import ClientException
  from httplib import HTTPException
  from multiprocessing import Process

  def submit(config):
      try:
          crabCommand('submit', config=config)
      except HTTPException as hte:
          print "Failed submitting task: %s" % (hte.headers)
          print hte
      except ClientException as cle:
          print "Failed submitting task: %s" % (cle)


  parser = ArgumentParser()
  parser.add_argument('-y', '--yaml', default = 'samples_ffr_data.yaml', help = 'File with dataset descriptions')
  parser.add_argument('-f', '--filter', default='*', help = 'filter samples, POSIX regular expressions allowed')
  args = parser.parse_args()

  with open(args.yaml) as f:
    doc = yaml.load(f) # Parse YAML file
    common = doc['common'] if 'common' in doc else {'data' : {}, 'mc' : {}}
    
    # loop over samples
    for sample, info in doc['samples'].iteritems():
      print("\n\n*** Sample {} ***".format(sample))
      # Given we have repeated datasets check for different parts
      parts = info['parts'] if 'parts' in info else [None]
      for part in parts:
        name = sample % part if part is not None else sample
        
        # filter names according to what we need
        if not fnmatch(name, args.filter): continue
        print 'submitting', name

        isMC = info['isMC']

        this_config = config()
        this_config.section_('General')
        this_config.General.transferOutputs = True
        this_config.General.transferLogs = True
        this_config.General.workArea = 'BParkingNANO_%s' % production_tag

        this_config.section_('Data')
        this_config.Data.publication = False
        #this_config.Data.outLFNDirBase = '/store/group/cmst3/group/bpark/%s' % (this_config.General.workArea)
        this_config.Data.outLFNDirBase = '/store/user/{}/BParkingNANO/{}/'.format(getUsernameFromCRIC(), skim_version)

        this_config.Data.inputDBS = 'global'

        this_config.section_('JobType')
        this_config.JobType.pluginName = 'Analysis'
        this_config.JobType.psetName = '../test/run_nano_FFR_AllJpsiMuMu_cfg.py'
        #this_config.JobType.maxJobRuntimeMin = 2750
        this_config.JobType.allowUndistributedCMSSW = True

        this_config.section_('User')
        this_config.section_('Site')
        this_config.Site.storageSite = 'T3_US_Brown'

        this_config.Data.inputDataset = info['dataset'] % part \
                                   if part is not None else \
                                      info['dataset']

        this_config.General.requestName = name
        common_branch = 'mc' if isMC else 'data'
        #this_config.Data.splitting = 'FileBased' if isMC else 'LumiBased'
        this_config.Data.splitting = 'FileBased' if isMC else 'Automatic'
        if not isMC:
            this_config.Data.lumiMask = info.get(
                'lumimask', 
                common[common_branch].get('lumimask', None)
            )
        else:
            this_config.Data.lumiMask = ''

        this_config.Data.unitsPerJob = info.get(
            'splitting',
            common[common_branch].get('splitting', None)
        )

        if "totalUnits" in info:
          this_config.Data.totalUnits = info.get("totalUnits")

        globaltag = info.get(
            'globaltag',
            common[common_branch].get('globaltag', None)
        )
        
        this_config.JobType.pyCfgParams = [
            'isMC=%s' % isMC, 'reportEvery=1000',
            'tag=%s' % production_tag,
            'globalTag=%s' % globaltag,
        ]
        
        this_config.JobType.outputFiles = ['_'.join(['BParkNANO', 'mc' if isMC else 'data', production_tag])+'.root']
        
        print this_config
        p = Process(target=submit, args=(this_config,))
        p.start()
        p.join()
        #submit(this_config)
      print("*** Done with Sample {} ***\n\n".format(sample))

