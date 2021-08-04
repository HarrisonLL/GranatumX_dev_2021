FROM granatumx/gbox-py-sdk:1.0.0

COPY . .

# uncomment the following line if you need to install dependencies
RUN pip install --upgrade pip
# then install required dependencies here
RUN pip install synapseclient
RUN pip install psutil
RUN pip install sympy

RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz  package.yaml
