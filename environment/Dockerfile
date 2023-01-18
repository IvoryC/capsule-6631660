# hash:sha256:bebfa1ae8bcb22173080d767db848e9fd36478ff26b3a5c0fa9d6c45015935cd
FROM registry.codeocean.com/codeocean/r-base:3.6.3-ubuntu18.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        python3-pip \
        python3-setuptools \
        python3-wheel \
    && rm -rf /var/lib/apt/lists/*

#RUN pip3 install -U --no-cache-dir \
#    matplotlib \
#    numpy \
#    pandas \
#    scikit-learn \
#    scipy \
#    seaborn

RUN Rscript -e 'remotes::install_version("ggplot", "0.4.2")' \
    && Rscript -e 'remotes::install_version("ggvenn", "0.1.9")' \
    && Rscript -e 'remotes::install_version("glmnet", "4.1-4")' \
    && Rscript -e 'remotes::install_version("tidyverse", "1.3.1")' \
    && Rscript -e 'remotes::install_version("torch", "0.7.2")'

RUN Rscript -e 'options(warn=2); install.packages("torch")'
RUN Rscript -e 'options(warn=2); install.packages("tidyverse")'
RUN Rscript -e 'options(warn=2); install.packages("glmnet")'
RUN Rscript -e 'options(warn=2); install.packages("gbm")'
RUN Rscript -e 'options(warn=2); install.packages("doMC")'
RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/BiocManager/BiocManager_1.30.16.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.2.3.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2); install.packages("minqa")'
RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/lme4/lme4_1.1-27.1.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/plyr/plyr_1.8.6.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2); install.packages(c("e1071", "ModelMetrics", "pROC", "reshape2", "ipred", "timeDate", "PRROC", "cowplot"))'
RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/gower/gower_0.2.2.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2);  devtools::install_github("tidymodels/recipes")'


RUN Rscript -e 'options(warn=2); install.packages("https://cran.r-project.org/src/contrib/Archive/caret/caret_6.0-90.tar.gz", repos=NULL, type="source")'
RUN Rscript -e 'options(warn=2); install.packages("MLmetrics")'


RUN Rscript -e 'options(warn=2); BiocManager::install(c( \
        "decontam", \
        "philentropy", \
        "biomformat",\
        "limma",\
        "edgeR",\
       "vegan",\
        "snm"\
    ) )'



 #   scikit-learn \
 #   scipy \

#RUN pip3 install -U --no-cache-dir \
#    matplotlib \
#    numpy \
#    pandas \
#    seaborn

RUN pip3 install numpy
RUN pip3 install pandas
#RUN apt-get install libpython3-dev

#RUN apt-get install -y --no-install-recommends python3-dev
RUN pip3 install pillow==6.2.1
RUN pip3 install seaborn
RUN pip3 install scikit-learn
RUN pip3 install setuptools
# RUN pip3 install cmake
# RUN pip3 install xgboost==1.5.0


RUN Rscript -e 'options(warn=2); devtools::install_github("Shenhav-and-Korem-labs/SCRuB")'
RUN Rscript -e 'options(warn=2); devtools::install_github("donaldtmcknight/microDecon")'

COPY postInstall /
RUN /postInstall
