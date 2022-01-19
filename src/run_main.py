#
# run_main
#
# Template wrapper script that runs a command-line program within a Docker container.
#
import os, subprocess, sys
from datetime import datetime
# sys.path.append('global_utils/src/')
sys.path.append('../../global_utils/src/')
import module_utils


def runOtherPre( input_dir, output_dir, run_json ):
    """ This function is used to run any other commands BEFORE the main program has run.
    run_json has most of what you might need to run other commands, and has the following structure:

    run_json = {'module': module_name, 'local_input_dir': <LOCAL_INPUT_DIR>, 'local_output_dir': <LOCAL_OUT_DIR>, \
		'remote_input_dir': remote_input_directory, 'remote_output_dir': remote_output_directory, \
                'program_arguments': program_arguments, 'run_arguments': run_arguments_json, \
                'module_instance_json': module_instance_json}

    LOCAL_INPUT_DIR has any downloaded files. LOCAL_OUT_DIR has any output data or log files that will be uploaded.
    
    If you are not running any other commands or post-processing, then leave this function blank.

    Note that this function may change run_json - e.g., converting file formats.
    """
    # change BAM to SAM in program arguments -> BWA MEM outputs SAM but we will convert to BAM in post-processing
    pargs = run_json['program_arguments']
    pargs = pargs.replace('.bam ', '.sam ')
    run_json['program_arguments'] = pargs
    return run_json


def runOtherPost( input_dir, output_dir, run_json ):
    """ This function is used to run any other commands AFTER the main program has run.
    run_json has most of what you might need to run other commands, and has the structure shown above.

    If you are not running any other commands or pre-processing, then leave this function blank.
    """
    # run samtools sort
    samtools_sort_cmd = ['samtools sort -o {} {}'.format(run_json['local_output_file'], str(run_json['local_output_file']).replace('.bam','.sam'))]
    subprocess.call(samtools_sort_cmd, shell=True)
    
    # run samtools index
    samtools_index_cmd = ['samtools index {}'.format(run_json['local_output_file'])]
    subprocess.call(samtools_index_cmd, shell=True)
    return


def runMain():
    # time the duration of module container run
    run_start = datetime.now()
    print('Container running...')
    
    # initialize program run
    run_json = module_utils.initProgram()
    
    # do any pre-processing (specific to module)
    run_json = runOtherPre( run_json['local_input_dir'], run_json['local_output_dir'], run_json )
    
    # run main program
    module_utils.runProgram( run_json['program_arguments'], str(run_json['local_output_file']).replace('.bam','.sam') )
    
    # do any post-processing
    runOtherPost( run_json['local_input_dir'], run_json['local_output_dir'], run_json )
    
    # create run log that includes program run duration
    run_end = datetime.now()
    run_json['module_run_duration'] = str(run_end - run_start)
    module_utils.logRun( run_json, run_json['local_output_dir'] )
    
    # upload output data files
    module_utils.uploadOutput( run_json['local_output_dir'], run_json['remote_output_dir'] )
    print('DONE!')
    
    return


if __name__ == '__main__':
    runMain()
