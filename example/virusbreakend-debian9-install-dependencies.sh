# Use conda if you can.
# "conda create -n virusbreakend samtools bwa kraken2 repeatmasker"
# is much easier and doesn't require root permissions
#
# Example script of how to install virusbreakend dependencies into /opt/tools/*
# on a clean Debian 9 installation
sudo apt-get -y install dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 perl-modules make
sudo apt-get -y install openjdk-8-jdk r-base libcurl4-gnutls-dev libxml2-dev libssl-dev cpanminus
sudo apt-get -y install time libbz2-dev liblzma-dev libz-dev libssl-dev libcurl4-openssl-dev
sudo apt-get -y install python3-h5py
sudo cpanm install Text::Soundex
sudo apt-get -y install build-essential git automake
## GRIDSS
gridss_version=2.10.0
sudo mkdir -p /opt/tools/gridss/${gridss_version}
wget https://github.com/PapenfussLab/gridss/releases/download/v${gridss_version}/gridss-${gridss_version}.tar.gz
sudo tar zxvf gridss-${gridss_version}.tar.gz  -C /opt/tools/gridss/${gridss_version} 
## Kraken2
kraken_version=2.0.9
sudo mkdir -p /opt/tools/kraken2/${kraken_version}
wget https://github.com/DerrickWood/kraken2/archive/v${kraken_version}-beta.tar.gz
tar zxvf v${kraken_version}-beta.tar.gz
cd kraken2-${kraken_version}-beta
sudo ./install_kraken2.sh /opt/tools/kraken2/${kraken_version}
cd -
## RepeatMasker
# Tandem Repeat Finder
sudo mkdir -p /opt/tools/trf/4.0.9
wget https://github.com/Benson-Genomics-Lab/TRF/releases/download/v4.09.1/trf409.linux64
sudo mv trf409.linux64 /opt/tools/trf/4.0.9/trf
sudo chmod a+x /opt/tools/trf/4.0.9/trf
# RMBlast
sudo mkdir -p /opt/tools/rmblast/2.10.0
wget http://www.repeatmasker.org/rmblast-2.10.0+-x64-linux.tar.gz
tar zxvf rmblast-2.10.0+-x64-linux.tar.gz
sudo mv rmblast-2.10.0/bin/* /opt/tools/rmblast/2.10.0/
# HMMER
sudo mkdir -p /opt/tools/hmmer/3.3
wget http://eddylab.org/software/hmmer/hmmer-3.3.tar.gz
tar zxvf hmmer-3.3.tar.gz
cd hmmer-3.3
./configure -prefix=/opt/tools/hmmer/3.3 && make
sudo make install
cd -
# RepeatMasker
sudo mkdir -p /opt/tools/repeatmasker
wget http://www.repeatmasker.org/RepeatMasker-4.1.1.tar.gz
sudo tar zxvf RepeatMasker-4.1.1.tar.gz -C /opt/tools/repeatmasker/ 
sudo mv /opt/tools/repeatmasker/RepeatMasker /opt/tools/repeatmasker/4.1.1
cd /opt/tools/repeatmasker/4.1.1
sudo perl configure -default_search_engine rmblast -rmblast_dir /opt/tools/rmblast/2.10.0 -trf_prgm /opt/tools/trf/4.0.9/trf -hmmer_dir /opt/tools/hmmer/3.3
cd -
# create the human library cache
echo "AAAAAAAAAAAAAAAAAAAAA" > rmcache.fa
sudo /opt/tools/repeatmasker/4.1.1/RepeatMasker -species human rmcache.fa
sudo rm -r rmcache.fa*
# NCBI SRA toolkit
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.8/setup-apt.sh
chmod +x setup-apt.sh
sudo ./setup-apt.sh

#source /etc/profile.d/sra-tools.sh
# vdb-config --interactive # set up for GCE downloads (but not charges)