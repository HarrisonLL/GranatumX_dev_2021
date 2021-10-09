FROM granatumx/gbox-py-sdk:1.0.0

COPY . .


RUN pip install pandas==1.1.0
RUN pip install scanpy==1.8.1
RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
