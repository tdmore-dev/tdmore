# Use the Rocker R images
FROM rocker/tidyverse:latest

WORKDIR /app/tdmore

## Install package dependencies
# for RxODE
RUN apt-get update && apt-get install -y libudunits2-dev

## Install phantomjs / shinytest
RUN apt-get update && apt-get install -y phantomjs lbzip2
RUN R -e 'devtools::install_version("shinytest")'
RUN R -e 'webdriver::install_phantomjs()'
# Add local phantomjs to path
ENV PATH=/root/bin:$PATH
ENV QT_QPA_PLATFORM=minimal

## Install R package dependencies manually, to enable caching
RUN R -e 'devtools::install_version("deSolve")'
RUN R -e 'devtools::install_version("gridExtra")'
RUN R -e 'devtools::install_version("mvnfast")'
RUN R -e 'devtools::install_version("numDeriv")'
RUN R -e 'devtools::install_version("RcppEigen")'
RUN R -e 'devtools::install_version("zoo")'
RUN R -e 'devtools::install_version("dparser")'
RUN R -e 'devtools::install_version("farver")'
RUN R -e 'devtools::install_version("ggforce")'
RUN R -e 'devtools::install_version("inline")'
RUN R -e 'devtools::install_version("lotri")'
RUN R -e 'devtools::install_version("n1qn1")'
RUN R -e 'devtools::install_version("polyclip")'
RUN R -e 'devtools::install_version("PreciseSums")'
RUN R -e 'devtools::install_version("R.methodsS3")'
RUN R -e 'devtools::install_version("R.oo")'
RUN R -e 'devtools::install_version("R.utils")'
RUN R -e 'devtools::install_version("tweenr")'
RUN R -e 'devtools::install_version("units")'
RUN R -e 'devtools::install_version("classInt")'
RUN R -e 'devtools::install_version("e1071")'
RUN R -e 'devtools::install_version("fastGHQuad")'
RUN R -e 'devtools::install_version("flextable")'
RUN R -e 'devtools::install_version("gdtools")'
RUN R -e 'devtools::install_version("huxtable")'
RUN R -e 'devtools::install_version("lbfgs")'
RUN R -e 'devtools::install_version("minqa")'
RUN R -e 'devtools::install_version("officer")'
RUN R -e 'devtools::install_version("StanHeaders")'
RUN R -e 'devtools::install_version("vpc")'
RUN R -e 'devtools::install_version("zip")'
RUN R -e 'devtools::install_version("lbfgsb3c")'
RUN R -e 'devtools::install_version("mnormt")'
RUN R -e 'devtools::install_version("optextras")'
RUN R -e 'devtools::install_version("optimr")'
RUN R -e 'devtools::install_version("Rcgmin")'
RUN R -e 'devtools::install_version("Rvmmin")'
RUN R -e 'devtools::install_version("setRNG")'

## Installing development dependencies
RUN R -e 'devtools::install_version("bookdown")'
RUN R -e 'devtools::install_version("diffobj")'
RUN R -e 'devtools::install_version("fontBitstreamVera")'
RUN R -e 'devtools::install_version("fontLiberation")'
RUN R -e 'devtools::install_version("fontquiver")'
RUN R -e 'devtools::install_version("freetypeharfbuzz")'
RUN R -e 'devtools::install_version("miniUI")'
RUN R -e 'devtools::install_version("svglite")'
RUN R -e 'devtools::install_version("vdiffr")'

## Install package dependencies from DESCRIPTION
COPY DESCRIPTION .
RUN R -e 'devtools::install_deps()'
RUN R -e 'devtools::install_dev_deps()'

## Install full package
COPY . .
RUN R -e 'devtools::install(dependencies=TRUE)'

## BUGFIX: Manually update devtools
## Due to issue r-lib/devtools#2129
## Can be removed once devtools CRAN version is updated
RUN R -e 'if( packageVersion("devtools") != "2.2.1") stop("Devtools package version: ", packageVersion("devtools"), "; bugfix in Dockerfile not needed anymore, please remove!")'
RUN R -e 'remotes::install_github("r-lib/devtools")'

