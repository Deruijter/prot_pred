#!/bin/bash

> forest.out
for filename in ../feature_sets/*.csv; do
	python 3_machine_learn_data_sets.py -d $(basename $filename)
done
