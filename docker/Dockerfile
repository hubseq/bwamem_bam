# Use Python image base from Dockerhub
FROM python:3.7

# Metadata
LABEL container.base.image = "python:3.7"
LABEL software.name = "bwamem_bam"

# Install system packages and AWS CLI
# RUN printf "deb http://archive.debian.org/debian/ jessie main\ndeb-src http://archive.debian.org/debian/ jessie main\ndeb http://security.debian.org jessie/updates main\ndeb-src http://security.debian.org jessie/updates main" > /etc/apt/sources.list
RUN apt-get -y update && \
    apt-get -y install python3-pip bwa
RUN pip install awscli boto3

WORKDIR /

# install BWA (manual)
# RUN mkdir /bwa-0.7.17
# COPY misc/bwa-precompiled-0-7-17.tar.gz /bwa-0.7.17
# RUN tar xvf /bwa-0.7.17/bwa-precompiled-0-7-17.tar.gz
# COPY misc/bwa-0.7.17.tar.bz2 /
# RUN tar xjf bwa-0.7.17.tar.bz2 && cd bwa-0.7.17/ && make

# install samtools
RUN mkdir /samtools
ENV PATH="/samtools/:$PATH"
COPY samtools/bin/* /samtools/

# copy relevant scripts to root dir
RUN mkdir -p global_utils/src
COPY src/run_main.py /
COPY src/global_utils/src/module_utils.py /global_utils/src/
COPY src/global_utils/src/file_utils.py /global_utils/src/
COPY src/global_utils/src/global_keys.py /global_utils/src/
COPY src/global_utils/src/aws_s3_utils.py /global_utils/src/
ENV PATH="/usr/local/bin/:$PATH"

ENTRYPOINT ["python","/run_main.py"] 
