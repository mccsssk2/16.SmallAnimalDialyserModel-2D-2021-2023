#!/bin/bash
#
#
dia:
	gcc -o dialyser main.c -lm
	make run
	
run:
	./dialyser

viz:
	python3 viz.py
	
clean:
	rm -rf dialyser *.dat *.png
