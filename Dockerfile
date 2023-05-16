FROM ubuntu:latest
# installing without interactive dialogue
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt install --yes  --no-install-recommends \
# install text editors
nano \
vim \
# install python3
python3 \
python3-pip \
# install perl
perl \
#install java
openjdk-17-jre-headless \
#install git and cmake for running setup.py
git-all \
cmake \
# install prerequisites for samtools 
wget \
make \
gcc \
zlib1g-dev \
libbz2-dev \
liblzma-dev \
libcurl4-gnutls-dev \
libssl-dev \
libncurses5-dev \
bzip2 \
#install ncbi-entrez-direct
ncbi-entrez-direct \
#install prerequisites for sratoolkit
libxml2-dev \
#install rsem for bam_to_bigWig
rsem \
#install cwl runner
cwltool \
tabix


#cpan prerequisites for sratoolkit
RUN cpan install -T XML::LibXML && cpan install -T URI

RUN wget https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2\
  && tar xvjf samtools-1.17.tar.bz2 \
  && cd samtools-1.17 \
  && make 
ENV PATH="/samtools-1.17:${PATH}"

# install sratoolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.0.2/sratoolkit.3.0.2-ubuntu64.tar.gz \
  && tar -xf sratoolkit.3.0.2-ubuntu64.tar.gz 
ENV PATH="/sratoolkit.3.0.2-ubuntu64/bin:${PATH}"

# install wigToBigWig
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
RUN chmod 777 wigToBigWig
RUN mv wigToBigWig /usr/local/bin/

# install bwa
RUN wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2 \
  && tar xvjf bwa-0.7.17.tar.bz2 \
  && cd bwa-0.7.17 \
  && make CC='gcc -fcommon' 
ENV PATH="/bwa-0.7.17:${PATH}"
 
# install pysam
RUN pip3 install pysam

# echo container env path
RUN echo $PATH  

RUN rm *.gz \
   && rm *.bz2

# create symlink to python3 and alias python to python3 
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN echo 'alias python='python3'' >> ~/.bashrc
RUN /bin/bash -c "source ~/.bashrc" 

# create one folder for storaging repo files
RUN mkdir /opt/RNA_repo
COPY . /opt/RNA_repo/
WORKDIR /opt/RNA_repo
RUN python3 setup.py install

RUN mv picard /usr/local/bin/ \
   && mv trimmomatic /usr/local/bin/ \
   && chmod +x /usr/local/bin/picard \
   && chmod +x /usr/local/bin/trimmomatic

# add gatk to path
ENV PATH="/opt/RNA_repo/rnannot/lib/gatk-4.4.0.0:${PATH}"
RUN chmod +x /opt/RNA_repo/rnannot/lib/gatk-4.4.0.0/gatk

# add FastQC to path 
ENV PATH="/opt/RNA_repo/rnannot/lib/FastQC:${PATH}"
RUN chmod +x /opt/RNA_repo/rnannot/lib/FastQC/fastqc
# create one folder for binding to local and let it as working directory
RUN mkdir /opt/output
WORKDIR /opt/output






 
