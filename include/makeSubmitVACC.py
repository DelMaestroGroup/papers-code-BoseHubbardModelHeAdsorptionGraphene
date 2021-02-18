#!/usr/bin/env python
# Filename: makeSubmitVacc.py
import os
class submit:
    inputFile = None
    fileName = None

    def __init__(self, inputFile=None, fileName='submit.pbs'):
        """Fix this description""" #FIXME
        self.inputFile = input('input file: ') if inputFile is None else inputFile
        self.fileName = fileName

    def pbs_script(self, nodes=1, ppn=1, walltime='30:00:00', jobName='JOBNAME', email=None, emailOptions='bea', queueType=None):
        """Fix this description""" #FIXME
        print("FIXME: newer version of this implemented on some julia notebook")

        # Open the pbs file and write its header
        fileName = self.fileName
        pbsFile = open(fileName,'w')

        header = '''#!/bin/bash
#PBS -S /bin/bash
#PBS -l nodes={0}:ppn={1}
#PBS -l walltime={2}
#PBS -N {3}
#PBS -V
#PBS -j oe
#PBS -o out/${{PBS_JOBNAME}}-${{PBS_ARRAYID}}-${{PBS_JOBID}}'''.format(nodes,ppn,walltime,jobName)

        if email:
            header += '''
#PBS -M {0}
#PBS -m bea'''.format(email)

        if queueType:
            header += '''
#PBS -q {0}'''.format(queueType)
    
        header += '''
cd $PBS_O_WORKDIR

echo "Starting PBS script submit-pimc.pbs at:`date`" 
echo "  host:       ${PBS_O_HOST}"
echo "  node:       `cat ${PBS_NODEFILE}`"
echo "  jobid:      ${PBS_JOBID}"
echo "  array job:  ${PBS_ARRAYID}" 

case ${PBS_ARRAYID} in
'''
        
        pbsFile.write(header)
        f = open(self.inputFile, 'r')
        idx = 0
        for line in f:
            pbsFile.write('{0})\n{1};;\n'.format(idx,line))
            idx += 1
        pbsFile.write('esac\necho \"Finished run at:`date`\"')
        pbsFile.close()
        f.close()

        print("Copy {} to the VACC.".format(fileName))
        print('\nSubmit jobs with: qsub -t 0-{0} {1}\n'.format(idx-1,os.path.basename(fileName)))
        self.totalNumCommands = idx



# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
if __name__ == "__main__":
    import argparse
    # parse the command line options
    #args = docopt(__doc__)
    #FIXME add argument parsing
    print("FIXME: this is not implemented properly")
    aaa = submit('./test.sh')
    aaa.pbs_script(email='test')
