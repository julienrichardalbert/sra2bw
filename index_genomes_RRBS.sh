#!/bin/bash
for x in ./*; do echo $x; cd $x; bismark_genome_preparation . --parallel 20; cd ..; done
