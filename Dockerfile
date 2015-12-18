FROM python:2.7
MAINTAINER dan.leehr@duke.edu

RUN useradd -m svr
ADD . /SVR_models
RUN chown -R svr /SVR_models
WORKDIR /SVR_models
RUN pip install -r requirements.txt
USER svr
