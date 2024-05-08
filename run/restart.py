#########################################################################################
#
#   Python Restart script for the submission of dissociation after electron impact
#   Written by O Bramley                                                 08.05.24
#   Maininly moving stuff from original run script to here
# 
#   Written to make submission of jobs to the HPC easier and quicker
#     
#
#########################################################################################
import json 

if __name__=="__main__":
    with open('inputs.json') as f:
        inputs=json.load(f)