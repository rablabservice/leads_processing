#!/bin/bash

# Default values
processed_dir="/mnt/coredata/processing/leads/data/processed"
subjects=()
scan_type="all"
specified_subjects=false

# Help function
print_help() {
    echo "Usage: cat_log.sh [options]"
    echo
    echo "Options:"
    echo "  -d, --processed_dir DIR    Path to the processed directory (default: /mnt/coredata/processing/leads/data/processed)"
    echo "  -s, --subjects SUBJ1 SUBJ2 Specify subjects to process (default: all directories in processed_dir)"
    echo "  -t, --scan_type TYPE       Type of scans to process (must be 'all', 'pet', or 'mri'; default: 'all')"
    echo "  -h, --help                 Display this help message"
    echo
    echo "Examples:"
    echo "  ./cat_log.sh -d /path/to/processed_dir -s subject1 subject2 -t pet"
}

# Parse arguments
if [ "$#" -eq 0 ]; then
    print_help
    exit 0
fi

while [[ "$#" -gt 0 ]]; do
    case $1 in
        -d|--processed_dir) processed_dir="$2"; shift 2 ;;
        -s|--subjects) shift; specified_subjects=true; while [[ "$1" != -* && "$1" != "" ]]; do subjects+=("$1"); shift; done ;;
        -t|--scan_type) scan_type="$2"; shift 2 ;;
        -h|--help) print_help; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; print_help; exit 1 ;;
    esac
done

# Default subjects to all directories in processed_dir if not specified
if [ ${#subjects[@]} -eq 0 ]; then
    for dir in "$processed_dir"/*/; do
        subjects+=("$(basename "$dir")")
    done
fi

# Function to get scan directories based on scan_type
get_scan_dirs() {
    local subject_dir="$1"
    local scan_dirs=()
    case $scan_type in
        "mri")
            scan_dirs=($(find "$subject_dir" -mindepth 1 -maxdepth 1 -type d -name "MRI-T1*"))
            ;;
        "pet")
            scan_dirs=($(find "$subject_dir" -mindepth 1 -maxdepth 1 -type d ! -name "MRI-T1*"))
            ;;
        "all")
            scan_dirs=($(find "$subject_dir" -mindepth 1 -maxdepth 1 -type d))
            ;;
        *)
            echo "Invalid scan_type: $scan_type"; exit 1 ;;
    esac
    echo "${scan_dirs[@]}"
}

# Iterate over subjects
for subject in "${subjects[@]}"; do
    subject_dir="$processed_dir/$subject"
    scan_dirs=($(get_scan_dirs "$subject_dir"))
    scan_dirs=($(printf '%s\n' "${scan_dirs[@]}" | sort))  # Sort scan directories alphabetically
    log_files_found=false

    for scan_dir in "${scan_dirs[@]}"; do
        log_dir="$scan_dir/log"
        if [ -d "$log_dir" ]; then
            logf=$(find "$log_dir" -type f -name "*.log" -printf "%T@ %p\n" | sort -n | tail -1 | cut -d' ' -f2-)

            if [ -z "$logf" ]; then
                continue
            fi

            log_files_found=true

            # Print hyphenated line 4 characters longer than the length of the log file's basename
            echo
            echo "////////////////////////////////////////////////////////////////////////////////////////"
            echo $(basename "$logf")

            # Print the contents of the log file
            cat "$logf"
        fi
    done

    if $specified_subjects && ! $log_files_found; then
        echo "No log files found for subject $subject"
    fi
done
