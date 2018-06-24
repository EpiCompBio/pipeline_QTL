##################################################
# Dockerfile for
# https://github.com/EpiCompBio/pipeline_QTL
##################################################


############
# Base image
############

FROM continuumio/miniconda3
# It runs on Debian GNU/Linux 8; use e.g. uname -a ; cat /etc/issue.net
# https://hub.docker.com/r/continuumio/miniconda/
# Or simply run:
# docker run --rm -ti continuumio/miniconda3
# docker run --rm -ti ubuntu


#########
# Contact
#########
MAINTAINER Antonio Berlanga-Taylor <a.berlanga@imperial.ac.uk>


#########################
# Update/install packages
#########################

# Install system dependencies
# If running on Debian and anaconda/miniconda image, use apt-get:
RUN apt-get update && apt-get upgrade -qy apt-utils

RUN apt-get install -qy gcc \
    g++ \
    tzdata \
    wget \
    bzip2 \
    unzip \
    sudo \
    bash \
    fixincludes

# Get plotting libraries:
RUN apt-get install -qy \
            inkscape \
            graphviz

#########################
# Install conda
#########################

# Miniconda:
#RUN cd /usr/bin \
#    && wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh \
#    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /usr/local/miniconda

#RUN export PATH="/usr/local/miniconda/bin:$PATH"

# Add conda channels, last to be added take priority
# Don't mix conda-forge and/or bioconda with defaults channel in R as packages
# will conflict with other and fail
# channels
RUN conda config --add channels bioconda ; \
    conda config --add channels conda-forge ; \
    conda config --remove channels defaults ; \
    conda config --remove channels r

# Update conda:
RUN conda update -y conda

#########################
# Install dependencies
#########################

##########
# Install all packages needed
# Major packages:
RUN conda install python=3.5 ; \
    conda install -y r ; \
    conda install -y git

# Install python packages:
RUN pip install --upgrade pip numpy ; \
    pip install cython ; \
    pip install pysam ; \
    pip install pandas ; \
    pip install future ruffus ; \
    conda install -y sphinx ; \
    pip install sphinxcontrib-bibtex

# Install CGAT tools:
#RUN wget --no-check-certificate https://raw.githubusercontent.com/CGATOxford/cgat/master/requires.txt ; \
#    cat requires.txt | grep -v "#" | xargs -n 1 pip install ; \
#    conda install -y alignlib-lite ; \
#    conda install -y bedtools ; \
#    conda install -y pybedtools ; \
#    conda install -y ucsc-wigtobigwig ; \
#    pip install git+git://github.com/AntonioJBT/CGATPipeline_core.git ; \
#    pip install cgat
 
# Install project specific packages:
RUN conda install -y r-docopt=0.4.5 r-data.table=1.10.4 r-ggplot2=2.2.1 ; \
    conda install -y r-matrixeqtl=2.1.1 -c bioconda ; \
    conda install -y plink ; \
    R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("cowplot", repos = "http://cran.us.r-project.org") ; library("cowplot")' ; \
    R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("ggthemes", repos = "http://cran.us.r-project.org") ; library("ggthemes")' ; \
    R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("readxl", repos = "http://cran.us.r-project.org") ; library("readxl")' ; \
#R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("bigpca", repos = "http://cran.us.r-project.org") ; library("bigpca")' ; \

# Install rpy2 with conda as pip version causes conflicts:
#RUN conda install -y rpy2
##########

##########
# Install packages that aren't on pip or conda:

# FlashPCA2 and dependencies:

# Eigen: http://eigen.tuxfamily.org/index.php?title=Main_Page
# Eigen has a conda recipe though: conda install -c menpo eigen
RUN cd ; \
    mkdir flashpca2 ; \
    wget http://bitbucket.org/eigen/eigen/get/3.3.4.tar.gz ; \
    tar xvfz 3.3.4.tar.gz \
    mv eigen-eigen-5a0156e40feb eigen ; \
    cd eigen ;

# TO DO: Get full paths to specify for flashpca

# Spectra: https://github.com/yixuan/spectra/
RUN cd ; \
    cd flashpca2 ; \
    wget https://github.com/yixuan/spectra/archive/v0.5.0.tar.gz ; \
    tar xvfz v0.5.0.tar.gz ; \
    mv spectra-0.5.0 spectra ; \
    cd spectra ;

# Boost has a conda recipe:
RUN conda install -y boost

# FlashPCA2: https://github.com/gabraham/flashpca/releases
RUN cd ; \
    cd flashpca2 ; \
    wget https://github.com/gabraham/flashpca/releases/download/v2.0/flashpca_2.0.tar.gz ; \
    tar xvfz flashpca_2.0.tar.gz ; \
    cd flashpca_2.0/ ; \
    make all EIGEN_INC=/XXXX/flashpca2/eigen \
    BOOST_INC=/XXXX/flashpca2/boost-1.59.0/include \
    BOOST_LIB=/XXXX/flashpca2/boost-1.59.0/lib \
    SPECTRA_INC=/XXXX/flashpca2/spectra

# In a Mac e.g.: make all EIGEN_INC=/XXXX/flashpca2/eigen-eigen-5a0156e40feb \
#                         BOOST_INC=/XXXX/anaconda/envs/r_test/include/boost \
#                         BOOST_LIB=/XXXX/anaconda/envs/r_test/lib \
#                         SPECTRA_INC=/XXXX/flashpca2/spectra-0.5.0/include



##########


##############################
# Install package of interest
##############################

RUN pip install git+git://github.com/EpiCompBio/pipeline_QTL.git

############################
# Default action to start in
############################
# Only one CMD is read (if several only the last one is executed)
#ENTRYPOINT ['/xxx']
#CMD echo "Hello world"
#CMD project_quickstart.py
#CMD ["/bin/bash"]
CMD ["/bin/bash"]

# To build run as:
#docker build --no-cache=true -t antoniojbt/pipe_tests_alpine .

# To run e.g.:
# docker run --rm -ti antoniojbt/pipe_tests

# If mounting a volume do e.g.:
# docker run -v /host/directory:/container/directory --rm -ti antoniojbt/pipe_tests
# docker run -v ~/Documents/github.dir/docker_tests.dir:/home/ --rm -ti antoniojbt/pipe_tests_alpine

# Create a shared folder between docker container and host
#VOLUME ["/shared/data"]
