# Use the Rocker R images
FROM rocker/tidyverse:latest

# Copy the current package
WORKDIR /app/tdmore
COPY . /app/tdmore

# Install package
RUN R -e 'devtools::install(dependencies=TRUE)'
