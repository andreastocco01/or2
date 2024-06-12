#!/bin/bash

#greedy_0-1
# python ../tools/newexperiment.py \
#         --output greedy.txt \
#         --executable ../main \
#         --timelimit 120 \
#         --ninstances 20 \
#         --nnodes 1500 \
#         --parallel 8 \
#         --costortime cost \
#        "[0,1]"

#vns_200-201-202-203-204-205-206
python ../tools/newexperiment.py \
        --output vns.txt \
        --executable ../main \
        --timelimit 120 \
        --ninstances 20 \
        --nnodes 1500 \
        --parallel 8 \
        --costortime cost \
        "[200,201,202,203,204,205,206]"

#fixten_3-4-5-6-7
# python ../tools/newexperiment.py \
#         --output fixten.txt \
#         --executable ../main \
#         --timelimit 120 \
#         --ninstances 20 \
#         --nnodes 1500 \
#         --parallel 8 \
#         --costortime cost \
#         "[3,4,5,6,7]"

#sinten_8-9-10-11-12-13-14-15
# python ../tools/newexperiment.py \
#         --output sinten.txt \
#         --executable ../main \
#         --timelimit 120 \
#         --ninstances 20 \
#         --nnodes 1500 \
#         --parallel 8 \
#         --costortime cost \
#         "[8,9,10,11,12,13,14,15]"

#b&c_17-18-19-20-21-22-23-24
# python ../tools/newexperiment.py \
#         --output bec.txt \
#         --executable ../main \
#         --timelimit 120 \
#         --ninstances 20 \
#         --nnodes 300 \
#         --parallel 1 \
#         --costortime time \
#         "[17,18,19,20,21,22,23,24]"

#diving_25-26-27-28-29-30
# python ../tools/newexperiment.py \
#         --output diving.txt \
#         --executable ../main \
#         --timelimit 120 \
#         --ninstances 20 \
#         --nnodes 1000 \
#         --parallel 1 \
#         --costortime cost \
#         "[25,26,27,28,29,30]"
