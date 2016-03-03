#!/bin/bash

docker run -it -v `pwd`:`pwd` dukegcb/svr_models ./predict_genome.py $@
