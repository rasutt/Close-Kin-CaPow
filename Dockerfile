FROM rocker/shiny

# Install Packages
RUN R -e 'install.packages(c("shiny"))'
RUN R -e 'install.packages(c("Matrix", "TMB"), type = "source")'

# Copy App
COPY ./app/ /srv/shiny-server/

# Compile objective function to avoid app startup time out
RUN R -e 'TMB::compile("/srv/shiny-server/TMB_files/UnifiedNLL.cpp")'

# run app
CMD ["/usr/bin/shiny-server"]