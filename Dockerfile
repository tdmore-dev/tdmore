# Use the Rocker R images
FROM rocker/tidyverse:latest

WORKDIR /app/tdmore

## Install compile-time package dependencies
RUN apt-get update && apt-get install -y libudunits2-dev

## Initialize renv, see https://environments.rstudio.com/docker#r-packages
RUN R -e 'install.packages("renv")'
COPY renv.lock .
RUN R -e 'renv::restore(library=.libPaths()[1])'
## this installs an .Rprofile and an renv directory
## We can deactivate renv again, as we restored to the library directory anyway

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
