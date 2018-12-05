################## BASE IMAGE ######################
FROM biocontainers/biocontainers:latest

################## METADATA ######################
LABEL base_image="biocontainers:latest"
LABEL version="1"
LABEL software="spectra-cluster-py"
LABEL software.version="v0.1"
LABEL about.summary="spectra-cluster-py Python API to process .clustering result files"
LABEL about.home="https://github.com/spectra-cluster/spectra-cluster-py"
LABEL about.documentation="https://github.com/spectra-cluster/spectra-cluster-py"
LABEL about.license_file="https://github.com/spectra-cluster/spectra-cluster-py/LICENSE"
LABEL about.license="SPDX:Apache-2.0"
LABEL about.tags="Proteomics"

################## MAINTAINER ######################
MAINTAINER Johannes Griss <jgriss@ebi.ac.uk>

USER root

RUN apt-get update && apt-get install -y python3 python3-pip python3-docopt python3-numpy python3-lxml python3-pandas python3-networkx && pip3 install setuptools

RUN pip3 install git+https://github.com/hollenstein/maspy

RUN wget -O /tmp/spectra_cluster_py.zip https://github.com/spectra-cluster/spectra-cluster-py/zipball/master && \
    unzip /tmp/spectra_cluster_py.zip -d /tmp/spectra_cluster_py && \
    cd /tmp/spectra_cluster_py/spec* && \
    python3 ./setup.py install && \
    cd / && \
    rm -rf /tmp/spectra_cluster_py*

USER biodocker
