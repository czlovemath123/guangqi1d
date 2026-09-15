#!/bin/bash

# check input
if [ $# -eq 0 ]
then
    echo "Need Input File"
    exit 1
fi

# check file whether exist
if [ ! -f "$1" ]
then
    echo "File '$1' not found!"
    exit 2
fi


file="$1"
temp="${1##*/}_tmp.txt"

count=0

tail -n +2 "$file" | sed 's/\r$//' > "$temp" 

while IFS= read -r line
do
    ((count++))
    python3 modify_parameters.py $line
    if [ $? -ne 0 ]; then
        echo "modify_parameters.py failed"
        exit 1
    else
        mkdir -p out
        mpiexec -np 4 ./guangqi </dev/null
        if [ $? -ne 0 ]; then
            echo "guangqi failed"
            rm "$temp"
            exit 1
        else
            echo "guangqi finished"
            python3 rename.py $line
            echo "Task $count finished"
        fi
    fi
done < "$temp"

rm "$temp"
