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
    def writeOtherBAM(run_json, bam_stats_file, samflags_filter_string, read_type, index_reads=False):
        bam_outfile = run_json['local_output_file'].replace('.bam','.{}.bam'.format(read_type)))
        other_bam_cmd = ['samtools view -h -b -o {} {} {}'.format(bam_outfile, samflags_filter_string, run_json['local_output_file'])]
        subprocess.call(other_bam_cmd, shell=True)
        run_json['postrun_commands'].append(other_bam_cmd[0])
        
        # write read statistics
        p1 = subprocess.Popen(['samtools view {}'.format(bam_outfile)], shell=True, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(['wc -l'], shell=True, stdin=p1.stdout, stdout=subprocess.PIPE)
        num_reads = str(p2.communicate()[0].decode('utf-8')).lstrip(' \t').rstrip(' \t\n')
        bam_stats_file.write('{},{}\n'.format(read_type, num_reads))
        
        if index_reads == True:
            samtools_index_cmd = ['samtools index {}'.format(bam_outfile)]
            subprocess.call(samtools_index_cmd, shell=True)
            run_json['postrun_commands'].append(samtools_index_cmd[0])
        
        return run_json

    
    run_json['postrun_commands'] = []
    
    # run samtools sort
    samtools_sort_cmd = ['samtools sort -o {} {}'.format(run_json['local_output_file'], str(run_json['local_output_file']).replace('.bam','.sam'))]
    subprocess.call(samtools_sort_cmd, shell=True)
    run_json['postrun_commands'].append(samtools_sort_cmd[0])
    
    # run samtools index
    samtools_index_cmd = ['samtools index {}'.format(run_json['local_output_file'])]
    subprocess.call(samtools_index_cmd, shell=True)
    run_json['postrun_commands'].append(samtools_index_cmd[0])

    # can remove SAM file
    rm_sam_cmd = ['rm {}'.format(str(run_json['local_output_file']).replace('.bam','.sam'))]
    subprocess.call(rm_sam_cmd, shell=True)
    run_json['postrun_commands'].append(rm_sam_cmd[0])

    # output basic read mapping statistics (samtools flagstat just counts lines without considering multimapping / supplemental are the same read)
    bam_stats_file = open(str(run_json['local_output_file']).replace('.bam','.alignment_stats.csv'),'w')
    bam_stats_file.write('read_type,count\n')
    
    p1 = subprocess.Popen(['samtools view -b -F 256 -F 2048 {}'.format(str(run_json['local_output_file']))], shell=True, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(['wc -l'], shell=True, stdin=p1.stdout, stdout=subprocess.PIPE)
    num_total_reads = str(p2.communicate()[0].decode('utf-8')).lstrip(' \t').rstrip(' \t\n')
    bam_stats_file.write('total,{}\n'.format(num_total_reads))
    
    # create any other BAM files specified in options
    # choices: mapped,unmapped, paired,properly_paired,paired_and_mapped,secondary,supplementary,duplicates,singleton,chimeric,R1,R2
    if 'options' in run_json['module_instance_json'] and run_json['module_instance_json']['options'] != '':
        other_bams = run_json['module_instance_json']['options'].split(',')
        if 'mapped' in other_bams:
            # mapped BAM file - may include secondary (multimap) and supplementary alignments - index also
            run_json = writeOtherBAM(run_json, bam_stats_file, '-F 0x4', 'mapped', True)
            # create uniquely mapped BAM file - no unmapped, secondary, or supplementary alignments - index also
            # see - https://www.biostars.org/p/138116/            
            run_json = writeOtherBAM(run_json, bam_stats_file, '-F 904', 'uniquely_mapped', True)
        
        if 'paired' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1', 'paired')            
            
        if 'properly_paired' in other_bams or 'properly-paired' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1 -f 0x2 -F 0x4', 'properly_paired')
        
        if 'mapped_and_paired' in other_bams or 'paired_and_mapped' in other_bams \
           or 'mapped-and-paired' in other_bams or 'paired-and-mapped' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1 -F 0x4', 'mapped_and_paired')                        
            
        if 'secondary' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x100', 'secondary')            
            
        if 'supplementary' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x800', 'supplementary')
            
        if 'duplicates' in other_bams or 'duplicate' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x400', 'duplicates')
            
        if 'singletons' in other_bams or 'singleton' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1 -f 0x8 -F 0x4', 'singletons')
            
        if 'chimeric' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1 -F 0x4 -F 0x8', 'chimeric')
            
        if 'R1' in other_bams or 'read1' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1 -f 0x40', 'R1')
            
        if 'R2' in other_bams or 'read2' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x1 -f 0x80', 'R2')
        
        if 'unmapped' in other_bams:
            run_json = writeOtherBAM(run_json, bam_stats_file, '-f 0x4', 'unmapped')
    
    bam_stats_file.close()
    return run_json


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
    run_json = runOtherPost( run_json['local_input_dir'], run_json['local_output_dir'], run_json )
    
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
