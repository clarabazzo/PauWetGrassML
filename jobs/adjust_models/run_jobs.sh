#!/bin/bash
# Read a string with spaces using for loop
declare -a arr=(17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 113 114 115 116 117 118 119 120 121 122 123 124 125 126 127 128 209 210 211 212 213 214 215 216 217 218 219 220 221 222 223 224 305 306 307 308 309 310 311 312 313 314 315 316 317 318 319 320)
for i in "${arr[@]}"
do
sbatch job$i.slurm
done
