import os
import sys
import pwd


def write_bash_script(command, output_sample_name, dag_dir):
    bash_name = '{}/{}.sh'.format(dag_dir+'inputs', output_sample_name)
    bashScript = '#!/bin/bash\n invalue=$(<$INPUT)\n echo "$invalue"\n'
    bashScript += 'echo $invalue\n'
    bashScript += 'value=${invalue:5}\n echo $value\n'
    copy_command = 'cp $value input_file.root'.replace(
        '/hdfs', '')
    bashScript += retryLogic(copy_command)
    bashScript += command
    bashScript += '\n rm input_file.root'
    bashScript += '\n'
    with open(bash_name, 'w') as file:
        file.write(bashScript)
    os.system('chmod +x {}'.format(bash_name))
    return bash_name


def retryLogic(command):
    return '''\nn=0
until [ $n -ge 5 ]
do
\techo "attempting copy for the ${{n}} time"
\t{} && break
\tn=$[$n+1]
done\n
'''.format(command)


def default_farmout(jobName, input_name, output_dir, bash_name, submit_dir, dag_dir, filesperjob):
    farmoutString = 'farmoutAnalysisJobs --infer-cmssw-path --fwklite --input-file-list=%s' % (
        input_name)
    farmoutString += ' --submit-dir=%s --output-dag-file=%s --output-dir=%s' % (
        submit_dir, dag_dir, output_dir)
    farmoutString += ' --input-files-per-job=%i %s %s ' % (
        filesperjob, jobName, bash_name)
    farmoutString += '--use-hdfs --memory-requirements=3000 --vsize-limit=8000'
    return farmoutString


def format_extra(extra_inputs, nfs_sample_dir):
    command = ''
    if extra_inputs != None and len(extra_inputs) > 0:
        command += ' --extra-inputs={}'.format(extra_inputs[0])
        for extra in extra_inputs:
            os.system('cp {} {}'.format(extra, nfs_sample_dir))
            if extra == extra_inputs[0]:
                continue
            command += ',{}'.format(extra)
    return command


def need_big_boi(bigboi):
    command = ''
    if bigboi:
        command += ' --vsize-limit=10000'
    return command


def submit_command(command, jobName, input_sample_dir, output_sample_name=None, dryrun=False, extra_inputs=None, bigboi=False):
    '''
    Submit a command to run on Wisconsin Condor

    This function is used to make condor submission easier. Provided a command and an input directory,
    submit_command will create all of the condor configs needed to run the command on all files in the
    input directory using condor. Additional parameters are needed for simple things like the job name
    and the name of the directory where you want the files to appear.

    USAGE: 
        from condor_handler import submit_command
        my_command = 'python ${{CMSSW_BASE}}/bin/${{SCRAM_ARCH}}/condor_classify.py -t {} -o $OUTPUT -f $value '.format(args.channel)
        inputs = 'path/to/inputs/DYJetsLongName'
        extra_files = ['file1.hdf5', 'file2.root']
        submit_command(command=my_command, jobName='myJobName', input_sample_dir=inputs, output_sample_name='DYJets', extra_inputs=extra_files)

    Parameters
    ----------
    command : str
        Command to be run on condor. No additional formatting is done after being provided
    jobName : str
        Name for this condor job
    input_sample_dir : str
        Path to the input files
    output_sample_name : str
        Name for the output samples. Can be used to shorten names like 
        DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8_v6_ext1-v2 -> DYJets
        Defaults to the input sample name if none are provided.
    dryrun : bool
        Create the farmoutAnalysisJob command, but don't submit
    extra_inputs : list of strings
        Extra input files to be copied to /nfs_scratch and then submitted with the job
    bigboi : bool
        Bump vsize up from 3GB to 10GB

    Returns
    -------
    None
        Nothing is returned, but the jobs are submitted to condor

    '''
    input_sample_name = os.path.basename(input_sample_dir)
    print "input_sample_name: ", input_sample_name
    if output_sample_name == None:
        output_sample_name = input_sample_name
    print "output_sample_name:", output_sample_name
    if input_sample_name == '':
        print "input_sample_name not defined, check for trailing '/' on input_sample_dir path"
        return
    else:
        nfs_sample_dir = '/nfs_scratch/{}/{}/{}'.format(
            pwd.getpwuid(os.getuid())[0], jobName, output_sample_name)

    # create submit dir
    submit_dir = '%s/submit' % (nfs_sample_dir)
    if os.path.exists(submit_dir):
        print('Submission directory exists for {} {}.'.format(
            jobName, output_sample_name))

    # create dag dir
    dag_dir = '{}/dags/dag'.format(nfs_sample_dir)
    os.system('mkdir -p {}'.format(os.path.dirname(dag_dir)))
    os.system('mkdir -p {}'.format(dag_dir+'inputs'))

    # output dir
    output_dir = 'gsiftp://cms-lvs-gridftp.hep.wisc.edu:2811//hdfs/store/user/{}/{}/{}/'.format(
        pwd.getpwuid(os.getuid())[0], jobName, output_sample_name)

    # create file list
    filelist = ['{}/{}'.format(input_sample_dir, x)
                for x in os.listdir(input_sample_dir)]
    filesperjob = 1
    input_name = '{}/{}.txt'.format(dag_dir+'inputs', output_sample_name)
    with open(input_name, 'w') as file:
        for f in filelist:
            file.write('%s\n' % f.replace('/hdfs', '', 1))

    # create bash script
    bash_name = write_bash_script(command, output_sample_name, dag_dir)

    # create farmout command
    farmoutString = default_farmout(
        jobName, input_name, output_dir, bash_name, submit_dir, dag_dir, filesperjob)
    farmoutString += format_extra(extra_inputs, nfs_sample_dir)
    farmoutString += need_big_boi(bigboi)

    if not dryrun:
        print('Submitting {}'.format(output_sample_name))
        os.system(farmoutString)
    else:
        print farmoutString

    return


if __name__ == "__main__":
    print 'This script cannot be used by itself. Please import it into your own script where you can provide an input command'
