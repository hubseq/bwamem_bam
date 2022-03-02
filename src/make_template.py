import json
MODULE = 'bwamem_bam'

mi_template_json = {'module_version': '00.00.00', 'program_name': 'bwa', 'program_subname': 'mem', 'program_version': '0.7.17', 'compute': {'environment': 'aws', 'language': 'Python', 'language_version': '3.7', 'vcpus': 2, 'memory': 6000}, 'program_arguments': '-S -t 2', 'program_input': [{'input_type': 'file', 'input_file_type': 'FASTQ', 'input_position': -1, 'input_prefix': '-i'}, {'input_type': 'file', 'input_file_type': 'FASTQ.GZ', 'input_position': -1, 'input_prefix': ''}], 'program_output': [{'output_type': 'file', 'output_file_type': 'BAM', 'output_position': 0, 'output_prefix': '-o'}], 'alternate_inputs': [{'input_type': 'file', 'input_file_type': 'BED', 'input_position': 0, 'input_prefix': '-L'}, {'input_type': 'folder', 'input_file_type': 'FASTA', 'input_position': -2, 'input_prefix': ''}], 'alternate_outputs': [{'output_type': 'file', 'output_file_type': 'BAM.BAI', 'output_position': -100, 'output_prefix': ''}], 'options': 'mapped', "defaults": {"alternate_inputs": ["s3://npipublicinternal/test/genomes/hg38/bwaindex/hg38.fasta"], "output_file": "<SAMPLE_ID>.bwamem_bam.bam"}}
with open(MODULE+'.template.json','w') as fout:
    json.dump(mi_template_json, fout)

io_dryrun_json = {'input': ['s3://npipublicinternal/test/dnaseq_targeted/run_test1/fastq/dnaseq_test_R1.fastq.gz', 's3://npipublicinternal/test/dnaseq_targeted/run_test1/fastq/dnaseq_test_R2.fastq.gz'], 'output': ['s3://npipublicinternal/test/dnaseq_targeted/run_test1/bwamem_bam/dnaseq_test.bwamem.bam'], 'alternate_inputs': ['s3://npipublicinternal/test/genomes/hg38/chr9/chr9.fasta'], 'alternate_outputs': ['s3://npipublicinternal/test/dnaseq_targeted/run_test1/bwamem_bam/dnaseq_test.bwamem.bam.bai'], 'program_arguments': '', 'sample_id': MODULE+'_test', 'dryrun': ''}
io_json = {'input': ['s3://npipublicinternal/test/dnaseq_targeted/run_test1/fastq/dnaseq_test_R1.fastq.gz', 's3://npipublicinternal/test/dnaseq_targeted/run_test1/fastq/dnaseq_test_R2.fastq.gz'], 'output': ['s3://npipublicinternal/test/dnaseq_targeted/run_test1/bwamem_bam/dnaseq_test.bwamem.bam'], 'alternate_inputs': ['s3://npipublicinternal/test/genomes/hg38/chr9/chr9.fasta'], 'alternate_outputs': ['s3://npipublicinternal/test/dnaseq_targeted/run_test1/bwamem_bam/dnaseq_test.bwamem.bam.bai'], 'program_arguments': '', 'sample_id': MODULE+'_test'}
io_dryrun_local_json = {'input': ['/Users/jerry/icloud/Documents/ngspipelines/bwamem/test/dnaseq_test_R1.fastq.gz', '/Users/jerry/icloud/Documents/ngspipelines/bwamem/test/dnaseq_test_R2.fastq.gz'], 'output': ['/Users/jerry/icloud/Documents/ngspipelines/bwamem/testout/dnaseq_test.bam'], 'alternate_inputs': ['/Users/jerry/icloud/Documents/ngspipelines/genomes/hg38/chr9/chr9.fasta'], 'alternate_outputs': ['/Users/jerry/icloud/Documents/ngspipelines/bwamem/testout/dnaseq_test.bam.bai'], 'program_arguments': '', 'sample_id': MODULE+'_test', 'dryrun': ''}
with open(MODULE+'.dryrun_test.io.json','w') as fout:
    json.dump(io_dryrun_json, fout)
with open(MODULE+'.test.io.json','w') as fout:
    json.dump(io_json, fout)
with open(MODULE+'.dryrun_local_test.io.json','w') as fout:
    json.dump(io_dryrun_local_json, fout)

# job info test JSONs                                                                                                        
job_json = {"container_overrides": {"command": ["--module_name", MODULE, "--run_arguments", "s3://npipublicinternal/test/modules/"+MODULE+"/job/"+MODULE+".test.job.json", "--working_dir", "/home/"]}, "jobqueue": "batch_scratch_queue", "jobname": "job_"+MODULE+"_test"}
with open(MODULE+'.test.job.json','w') as fout:
    json.dump(io_json, fout)
