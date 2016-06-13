#!/bin/bash

# join_list takes the delimiter as the first argument and prints
# the list passed as second argument, separated with the delimiter
function join_list {
 local d=$1; shift; echo -n "$1"; shift; printf "%s" "${@/#/$d}";
}

# file with sbatch settings
SFILE="$1"; shift
# name of job is first parameter
JOBNAME="$1"
# add list of MATLAB commands here
COMMAND="addpath('mps','routines','util','workspaces');"
COMMAND+="launch_simulation("`join_list \',\' "'$@'"`");"
COMMAND+="exit;"

echo "$COMMAND"
sbatch -J "$JOBNAME" "$SFILE" "$COMMAND"

