#!/bin/bash

for i in $(ls *.jl); do
	#julia /home/ale/jocco/jocco.jl $i
  julia /home/ale/Dropbox/PhD/bocco/jocco.jl $i
done
