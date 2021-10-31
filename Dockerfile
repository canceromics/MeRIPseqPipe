FROM nfcore/base:1.9
LABEL authors="Kaiyu Zhu, Xiaoqiong Bao" \
    description="Docker image containing all software requirements for the MeRIPseqPipe pipeline"

# Install the conda environment
COPY environment.yml /
RUN conda env create -f /environment.yml && conda clean -a
# install subread
RUN conda create -n multiqc -c conda-forge -c bioconda python=3.7.8 multiqc=1.7 && conda clean -a

RUN conda env export --name meripseqpipe-1.0dev > meripseqpipe-1.0dev.yml
ENV PATH /mspc:$PATH
ENV PATH /opt/conda/bin:$PATH
ENV PATH /opt/conda/envs/multiqc/bin/:$PATH
ENV PATH /opt/conda/envs/meripseqpipe-1.0dev/bin:$PATH


# install MATK
RUN wget https://github.com/kingzhuky/MATK_backup/releases/download/v0.1dev/MATK-1.0.jar

# install QNB
RUN wget https://cran.r-project.org/src/contrib/Archive/QNB/QNB_1.1.11.tar.gz && \ 
    R CMD INSTALL QNB_1.1.11.tar.gz && \
    rm QNB_1.1.11.tar.gz

# install MeTDiff
RUN git clone https://github.com/compgenomics/MeTDiff.git && \
    R CMD build MeTDiff/ && \
    R CMD INSTALL MeTDiff_1.0.tar.gz && \
    rm -rf MeTDiff*

# install MeTPeak
RUN git clone https://github.com/compgenomics/MeTPeak.git && \
    R CMD build MeTPeak/ && \
    R CMD INSTALL MeTPeak_1.0.0.tar.gz && \
    rm -rf MeTPeak*

# install MSPC
RUN conda install -y unzip 
RUN wget -O mspc.zip "https://github.com/Genometric/MSPC/releases/download/v5.4.0/linux-x64.zip" && \
    unzip mspc.zip -d mspc && \
    chmod 775 mspc/mspc && \ 
    rm mspc.zip
