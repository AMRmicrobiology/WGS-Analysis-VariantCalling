#!/bin/bash

download_folder="data/in/"

if read -e -p "Please provide the full path to your Acc_List.txt file: " filepath && [ -s "$filepath" ]; then 
    while IFS= read -r accession
    do
        echo "Processing accession: '$accession'"
        if [ -n "$accession" ]; then
            fastq-dump --split-files --gzip "$accession" -O "$download_folder"
        else
            echo "Skipping empty line."
        fi
    done < "$filepath"
    echo "Downloaded files are stored in : $download_folder"
else
    echo "ERROR: the specified file is empty or does not exist"
fi
