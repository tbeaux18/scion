## Emacs, make this -*- mode: sh; -*-

FROM ubuntu1804:base
MAINTAINER Timothy Baker <tbaker8@luc.edu>

RUN chmod 777 /tmp/

ENV DEBIAN_FRONTEND=noninteractive
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV R_BASE_VERSION 3.5.3


LABEL org.label-schema.license='GPL-2.0' \
      org.label-schema.vcs-url='https://github.com/rocker-org/r-base' \
      org.label-schema.vendor='Rocker Project' \
      maintainer='Dirk Eddelbuettel <edd@debian.org>'

RUN sed -Ei 's/^# deb-src /deb-src /' /etc/apt/sources.list \
  && apt-get update \
  && echo 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/' >> /etc/apt/sources.list \
  && gpg --keyserver hkp://keyserver.ubuntu.com:80 --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 \
  && gpg -a --export E298A3A825C0D65DFD57CBB651716619E084DAB9 | apt-key add - \
  && apt-get update \
  && apt-get -y install \
      libxml2-dev \
      libssh2-1-dev \
      libcurl4-openssl-dev \
      libmysqlclient-dev \
      mbuffer \
      parallel \
  && apt-get clean autoclean \
  && apt-get autoremove -y \
  && apt-get -y build-dep libcurl4-gnutls-dev \
  && apt-get -y install \
      libcurl4-gnutls-dev \
  && rm -rf /var/lib/apt/lists/*


## Configure default locale, see https://github.com/rocker-org/rocker/issues/19
RUN echo 'en_US.UTF-8 UTF-8' >> /etc/locale.gen \
	&& locale-gen en_US.utf8 \
	&& /usr/sbin/update-locale LANG=en_US.UTF-8


## Now install R and littler, and create a link for littler in /usr/local/bin
RUN apt-get update \
	&& apt-get install -y \
		littler \
    r-cran-littler \
		r-base=${R_BASE_VERSION}-* \
		r-base-dev=${R_BASE_VERSION}-* \
		r-recommended=${R_BASE_VERSION}-* \
	&& ln -s /usr/lib/R/site-library/littler/examples/install.r /usr/local/bin/install.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/install2.r /usr/local/bin/install2.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/installGithub.r /usr/local/bin/installGithub.r \
	&& ln -s /usr/lib/R/site-library/littler/examples/testInstalled.r /usr/local/bin/testInstalled.r \
	&& install.r docopt \
	&& rm -rf /tmp/downloaded_packages/ /tmp/*.rds \
	&& rm -rf /var/lib/apt/lists/*

# will move this up after ensuring cairo installation works
RUN apt-get update \
  && apt-get -y install \
  libcairo2-dev \
  libxt-dev

# installing R packages for zUMI to run
RUN \
  R -e "install.packages(c( \
          'yaml','shiny','shinythemes','shinyBS','ggplot2','mclust','dplyr', \
          'cowplot','Matrix','BiocManager','devtools','stringdist','data.table', \
          'plotly', 'crosstalk', 'inflection', 'DT', 'stringr', 'htmltools', 'parallel'))" \
  && R -e "BiocManager::install(c( \
            'GenomicRanges','GenomicFeatures', \
            'GenomicAlignments','AnnotationDbi','Rsubread', \
            'Glimma', 'edgeR', 'DESeq2'))" \
  && R -e "devtools::install_github('hadley/multidplyr')" \
  && R -e "devtools::install_github('VPetukhov/ggrastr')"
