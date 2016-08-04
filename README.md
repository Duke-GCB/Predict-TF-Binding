Predict-TF-Binding
==================

Python scripts for modeling and predicting TF-DNA binding using libsvm.

## Dependencies

1. Python with numpy
2. [libsvm](https://www.csie.ntu.edu.tw/~cjlin/libsvm/) providing binaries `svm-train` and `svm-predict`

## Training a model

To train a model (producing a `.model` file), use the **SVR_model_maker.py** script. See [SVR_modelmaker_README.md](./SVR_modelmaker_README.md) for details

`./SVR_model_maker.py -i scores.txt -c 1 -p 0.1 --searchstrings GCGG`

Where `scores.txt` is a tab-delimited file containing sequences and scores in the following format:

```
Name	ID	Sequence	Orientation1	Orientation2	Best-orientation	Replicate_intensity_difference
Bound2	00001	TAAGCCGGGCTTGCGGGCGCGCACACGTGGAACGCA	8.16296285181	8.58549476756	8.58549476756	-0.422531915752
Bound2	00002	TGGGCCGGGGACCTGGGCGCAGCCTCCCTCGCCGCA	6.37445679936	6.19073367654	6.37445679936	0.183723122823
Bound2	00003	TGTGCAGACGCCCGGCGCGCCTCCCGCTTAATCTGA	8.27945472383	8.19456431445	8.27945472383	0.0848904093827
Bound2	00004	TGTCCCTGCGTGCAGAGCGCGGTGAGAGTGGGTGGA	7.39775421771	7.79538055766	7.79538055766	-0.397626339952
```

## Generating predictions

Predictions can be generated on a FASTA sequence from a model file using the `predict_tf_binding.py` script:

    ./predict_tf_binding.py -g hg19.fa -m E2F1_SVR.model -c GCGG -w 36 -k 1 2 3 -o E2F1-GCGG.txt`

Predictions are generated in BED format, and may be processed from there. For full usage, see `predict_tf_binding.py -h`

