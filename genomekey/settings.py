import os, sys

"""
#from cosmos.config import settings as cosmos_settings


#def _get_drmaa_native_specification(jobAttempt):
#    task = jobAttempt.task
#    drm  = cosmos_settings['DRM']

    cpu_req  = task.cpu_requirement
    mem_req  = task.memory_requirement
    time_req = task.time_requirement
    queue    = task.workflow.default_queue
    
    if drm == 'LSF':           # for Orchestra Runs
        if time_req <= 12*60: queue = 'rodeo_unlimited'
        else:                 queue = 'long'
                
#        return '-R rusage[mem={0}] span[hosts=1] -n {1} -W 0:{2} -q {3}'.format(mem_req, cpu_req, time_req, queue)

    elif drm == 'GE':
        return '-l spock_mem={mem_req}M,num_proc={cpu_req}'.format(mem_req=mem_req, cpu_req=cpu_req)

    else:
        raise Exception('DRM not supported')


svr  = cosmos_settings['server_name']

# get path to current genomekey installation
genomekey_path = os.path.abspath(os.path.dirname(sys.argv[0]))

if svr == 'orchestra':
    ref_path   = '/groups/cbi/WGA/reference'
    tools_path = '/groups/cbi/WGA/tools'

elif svr == 'aws':
    ref_path   = '/WGA/reference'
    tools_path = '/WGA/tools'

elif svr == 'gce':
    ref_path   = '/pseq/WGA/ref'  # In shared disk
    tools_path = '/tools/'        # In boot   disk
else:
    raise Exception('Unknown server_name {0} in Cosmos configuration: must be one of [orchestra, aws, gce]'.format(svr))


opj = os.path.join
"""
ref_path   = '/groups/cbi/WGA/reference'
tools_path = '/groups/cbi/WGA/tools'

settings = {
    'date'                  : '$(date "+%T %D")',
    'java'                  : opj(tools_path, 'java -d64 -XX:ParallelGCThreads=2 -XX:+UseParallelOldGC -XX:+AggressiveOpts'),
    'scratch'               : '/tmp',

    'bamUtil'               : opj(tools_path, 'bamUtil.v1.0.11'),
    'bwa'                   : opj(tools_path, 'bwa.v0.7.7'),              
    'gatk'                  : opj(tools_path, 'gatk.v3.1.1.jar'),
    'picard_dir'            : opj(tools_path, 'picard.v1.109'),  
    'samtools'              : opj(tools_path, 'samtools.v0.1.19'),

    'annovarext'            : opj(tools_path, 'annovarext'),

    'reference_fasta'       : opj(ref_path,   'human_g1k_v37.fasta'),
    'dbsnp_vcf'             : opj(ref_path,   'dbsnp_137.b37.excluding_sites_after_129.vcf'),
    'hapmap_vcf'            : opj(ref_path,   'hapmap_3.3.b37.vcf'),
    'mills_vcf'             : opj(ref_path,   'Mills_and_1000G_gold_standard.indels.b37.vcf'),
    '1ksnp_vcf'             : opj(ref_path,   '1000G_phase1.snps.high_confidence.b37.vcf'),
    '1komni_vcf'            : opj(ref_path,   '1000G_omni2.5.b37.vcf'),
    '1kindel_vcf'           : opj(ref_path,   '1000G_phase1.indels.b37.vcf'),
    # keep in home directory for the moment, should ultimately be part of the AMI
    #'empty_sam'             : opj(ref_path,   'empty.sam'),
    'empty_sam'             : opj(genomekey_path, '../test/empty.sam'),

    'gatk_realigntarget'    : '',
    'gatk_indelrealign'     : '',

    'drmaa_spec'            : _get_drmaa_native_specification
}
"""
if   svr == 'aws':
    settings['scratch'] = '/mnt'
elif svr == 'gce':
    settings['scratch'] = '/mnt'

"""
