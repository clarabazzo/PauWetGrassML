#!/bin/bash

wd = '~/Clara/'

#--- make run list based on dashboard_adjust
dash_run = read.csv(paste0(wd,'jobs/adjust_models/dashboard_adjust.csv'))

#--- list make_runs
dash_run = dash_run[dash_run$make_run==1,]
run_list = paste(dash_run$job, collapse = ' ')

bash_temp = c(
'#!/bin/bash',
'# Read a string with spaces using for loop',
'declare -a arr=(<list>)',
'for i in "${arr[@]}"',
'do',
'sbatch job$i.slurm',
'done')

bash_temp = gsub('<list>',run_list, bash_temp)

write(bash_temp, 'run_jobs.sh')
message('done')

