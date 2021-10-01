FROM granatumx/gbox-py-sdk:1.0.0

COPY . .


RUN pip install scanpy==1.4.6
RUN ./GBOXtranslateVERinYAMLS.sh
RUN tar zcvf /gbox.tgz package.yaml yamls/*.yaml
