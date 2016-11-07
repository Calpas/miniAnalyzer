from WMCore.Configuration import Configuration
config = Configuration()

config.section_('General')
config.General.requestName     = 'SN'
config.General.transferOutputs = True

config.section_('JobType')
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = '../python/ConfFile_cfg.py'
config.JobType.outputFiles = ['output.root']
config.JobType.maxMemoryMB = 2500

config.section_('Data')
config.Data.inputDataset = 'DSET'
config.Data.inputDBS = 'global'
config.Data.unitsPerJob = 20000 # evt per job
#config.Data.totalUnits = 100000000 # tot evt to run one. Comment to run on all evt
config.Data.totalUnits = 100000000 # tot evt to run one. Comment to run on all evt. Max jos==10000!!!
config.Data.splitting = 'EventAwareLumiBased'
config.Data.outLFNDirBase = '/store/user/calpas/hto2taus/ntuple/v1'

config.section_('User')

config.section_('Site')
config.Site.storageSite = 'T2_PT_NCG_Lisbon'
