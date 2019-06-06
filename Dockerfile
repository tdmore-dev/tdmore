# Use the Rocker R images
FROM rocker/tidyverse:latest

# Copy the current package
WORKDIR /app/tdmore
COPY . /app/tdmore

# Install dependencies required for RxODE
RUN apt-get update
RUN apt-get install -y libudunits2-dev


# Install package
RUN R -e 'devtools::install(dependencies=TRUE)'
