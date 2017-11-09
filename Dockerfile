##################################################
# Dockerfile for project_quickstart 
# https://github.com/AntonioJBT/project_quickstart
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
# Don't use conda-forge with R as any packages will conflict with other
# channels
RUN /bin/bash -c 'conda config --add channels r \
    && conda config --add channels bioconda \
    && conda config --add channels defaults'

# There are dependency issues, specially rpy2, with clashes between conda-forge
# and default r (r-base), see here for icu 54 vs 56:
#https://github.com/conda-forge/conda-forge.github.io/issues/234
# but updates, new recipes, etc. changed things. 
# Solution was to update/re-install icu, which switches it to conda-forge icu58:
#RUN conda install -y icu

# Update conda:
RUN conda update -y conda

# Create a python 3.5 environment:
#RUN conda create -y -n py35 python=3.5

# Environments can't easily be sourced from a Dockerfile
# Mixing RUN with source errors, use instead this form of RUN:
# See:
# https://stackoverflow.com/questions/20635472/using-the-run-instruction-in-a-dockerfile-with-source-does-not-work
# https://docs.docker.com/engine/reference/builder/#run
# https://stackoverflow.com/questions/20635472/using-the-run-instruction-in-a-dockerfile-with-source-does-not-work/45087082#45087082
# Install everything in the virtual environment:
# Several cgat dependecies fail with pip so run with conda instead:
#RUN /bin/bash -c 'source activate py35 ; \

# Install all packages needed:
RUN conda install -y git ; \
    pip install --upgrade pip cython numpy ; \
    pip install pysam ; \
    pip install pandas ; \
    pip install future ruffus ; \
    conda install -y sphinx ; \
    pip install sphinxcontrib-bibtex ; \
    conda install -y r-base=3.3 r-MatrixEQTL=2.1.1 r-data.table=1.10.4 r-ggplot2=2.2.1 ; \
    # Run docopt second as some conflict with r versions if installed at the
    # same time:
    conda install -y r-docopt ; \
    #R --vanilla -e 'source("https://bioconductor.org/biocLite.R") ; install.packages("svglite", repos = "http://cran.us.r-project.org") ; library("svglite")' ; \
    wget --no-check-certificate https://raw.githubusercontent.com/CGATOxford/cgat/master/requires.txt ; \
    cat requires.txt | grep -v "#" | xargs -n 1 pip install ; \
    conda install -y alignlib-lite ; \
    conda install -y bedtools ; \
    conda install -y pybedtools ; \
    conda install -y icu -c conda-forge ; \
    pip install git+git://github.com/AntonioJBT/CGATPipeline_core.git ; \


# This fails with continuumio/miniconda3 in this file and manually:
#    pip install cgat
# This also fails:
#    pip install git+git://github.com/CGATOxford/cgat.git


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
