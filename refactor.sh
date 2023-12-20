#!/bin/bash

find . -type f \( -name "*.m" -o -name "*.R" -o -name "*.sh" -o -name "*.py" \) -exec sh -c '
    for file do
        awk -v file="$file" '{
            gsub("/LEADS/", "/leads6mm/");
            gsub("data_f7p1", "data");
            gsub("script_f7p1", "code/leads_processing");
            gsub("freesurfer_processing", "freesurfer/7.1.0");
            gsub("/screening/", "/data/petonly/screening/");
            gsub("/screening/", "/data/petonly/m36/");
        }1' "$file" > "$file.tmp" && mv "$file.tmp" "$file"
    done
' sh {} +
