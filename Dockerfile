FROM ubuntu:latest
# installing without interactive dialogue
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install --yes \
# install text editors
nano \
vim \
# install python3
python3 \
python3-pip \
# install perl
perl \
#install java
default-jre-headless \
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
rsem

#cpan prerequisites for sratoolkit
RUN cpan install -T XML::LibXML && cpan install -T URI

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.16.1/samtools-1.16.1.tar.bz2 \
  && bunzip2 -q samtools-1.16.1.tar.bz2 \
  && tar -xf samtools-1.16.1.tar \
  && cd samtools-1.16.1 \
  && make
ENV PATH="/samtools-1.16.1:${PATH}"

# install sratoolkit
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-ubuntu64.tar.gz \
  && tar -xf sratoolkit.2.10.0-ubuntu64.tar.gz
ENV PATH="/sratoolkit.2.10.0-ubuntu64/bin:${PATH}"

# install wigToBigWig
RUN wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/wigToBigWig
RUN chmod 777 wigToBigWig
RUN mv wigToBigWig /usr/local/bin/
 
# install pysam
RUN pip3 install pysam

# echo container env path
RUN echo $PATH  

# create symlink to python3 and alias python to python3 
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN echo 'alias python='python3'' >> ~/.bashrc
RUN /bin/bash -c "source ~/.bashrc" 

# create one folder for storaging repo files
RUN mkdir /opt/RNA_repo
COPY . /opt/RNA_repo/
WORKDIR /opt/RNA_repo
RUN python3 setup.py install

# create one folder for binding to local and let it as working directory
RUN mkdir /opt/output
WORKDIR /opt/output






 
