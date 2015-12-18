FROM python:2.7
MAINTAINER dan.leehr@duke.edu

# install LIBSVM.
# Owned by root and placed in /opt
ENV LIBSVM_VER 321
RUN curl -SL https://github.com/cjlin1/libsvm/archive/v${LIBSVM_VER}.tar.gz | tar -xzC /opt # makes /opt/libsvm-321
WORKDIR /opt/libsvm-${LIBSVM_VER}
RUN make

# Switch to non-root user
RUN useradd -m svr
USER svr

ENV PATH $PATH:/opt/libsvm-${LIBSVM_VER}
ADD . /SVR_models
WORKDIR /SVR_models

# Install requirements as root
USER root
RUN pip install -r requirements.txt

# switch back to svr
USER svr
