# Use the Rocker R images
FROM rocker/tidyverse:latest

WORKDIR /app/tdmore

## Install compile-time package dependencies
RUN apt-get update && apt-get install -y libudunits2-dev lbzip2 libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev libmpfr-dev libgmp-dev cmake

## Initialize renv, see https://environments.rstudio.com/docker#r-packages
RUN R -e 'install.packages("renv")'
COPY renv.lock .
RUN R -e 'renv::restore(library=.libPaths()[1])'
## this restores the renv.lock to the system library
## because it is the system library, no .Rprofile or renv/ directory is created

## Install phantomjs
RUN apt-get update && apt-get install -y lbzip2
RUN R -e 'shinytest::installDependencies()'
ENV OPENSSL_CONF /etc/ssl/

## Ensure all required packages were installed through renv
COPY DESCRIPTION .
RUN R -e 'stopifnot( all(remotes::local_package_deps(dependencies=TRUE) %in% rownames( installed.packages() )) )'

## Install full package; dependencies were installed earlier
COPY . .
RUN R -e 'devtools::install(dependencies=FALSE)'
