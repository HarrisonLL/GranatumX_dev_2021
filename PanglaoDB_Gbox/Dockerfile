FROM granatumx/gbox-r-sdk:1.0.0
# Change to FROM granatumx/gbox-py-sdk:1.0.0 for Python

COPY . .

# install required dependencies here
RUN Rscript install_packages.R

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz  package.yaml 

