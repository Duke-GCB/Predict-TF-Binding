FROM python:3.7
LABEL maintainer="dan.leehr@duke.edu"

# Add a non-root user
RUN useradd -m svr
USER svr

ADD requirements.txt /opt/predict-tf-binding/requirements.txt
WORKDIR /opt/predict-tf-binding

# Install requirements as root
USER root
RUN pip install -r requirements.txt

# switch back to svr
USER svr
ADD . /opt/predict-tf-binding/
ENV PATH /opt/predict-tf-binding:$PATH
CMD predict_tf_binding.py
