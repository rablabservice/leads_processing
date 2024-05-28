#!/bin/bash

# Default settings
cut_columns=(1 2)
max_line=88
file=""

# Function to display help
usage() {
    echo "Usage: $0 [options] file"
    echo ""
    echo "Options:"
    echo "  -c, --cut col1 col2 ...    Columns to omit (default: 1 2)"
    echo "  -h, --help                 Display this help message"
    echo ""
    echo "Arguments:"
    echo "  file                       The CSV file to process"
    exit 0
}

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -c|--cut) shift; cut_columns=($@); shift ${#cut_columns[@]};;
        -h|--help) usage;;
        *) file=$1; shift;;
    esac
done

# Validate input
if [[ -z "$file" ]]; then
    usage
fi

# Determine column widths
declare -a col_widths
num_columns=$(head -1 "$file" | awk -F',' '{print NF}')

# Initialize col_widths array with zero
for ((i=1; i<=num_columns; i++)); do
    col_widths[$i]=0
done

# Update col_widths with max length for each column
awk -F',' -v cut_cols="${cut_columns[*]}" '
BEGIN {
    split(cut_cols, cuts, " ")
    for (i in cuts) cut_cols[cuts[i]] = 1
}
{
    for (i=1; i<=NF; i++) {
        if (!(i in cut_cols)) {
            if (length($i) > col_widths[i]) {
                col_widths[i] = length($i)
            }
        }
    }
}
END {
    for (i=1; i<=length(col_widths); i++) {
        if (!(i in cut_cols)) {
            print col_widths[i]
        }
    }
}' "$file" > col_widths.tmp

# Read column widths into array
readarray -t col_widths < col_widths.tmp

# Function to apply alternating colors
apply_color() {
    local color_index=$1
    local text=$2
    if (( color_index % 2 == 0 )); then
        # Light gray
        printf "\e[47m\e[30m%s\e[0m" "$text"
    else
        # White
        printf "\e[107m\e[30m%s\e[0m" "$text"
    fi
}

# Print header with column widths and colors
awk -F',' -v cut_cols="${cut_columns[*]}" -v col_widths="${col_widths[*]}" '
BEGIN {
    split(cut_cols, cuts, " ")
    split(col_widths, widths, " ")
    color_index = 0
}
NR == 1 {
    for (i=1; i<=NF; i++) {
        if (!(i in cut_cols)) {
            printf "%s", (i > 1 ? " " : "")
            apply_color(color_index, sprintf("%-" widths[i-1] "s", $i))
            color_index++
        }
    }
    printf "\n"
}
function apply_color(color_index, text) {
    if (color_index % 2 == 0) {
        # Light gray
        printf "\033[47m\033[30m%s\033[0m", text
    } else {
        # White
        printf "\033[107m\033[30m%s\033[0m", text
    }
}
' "$file"

# Print rows with formatted columns
awk -F',' -v cut_cols="${cut_columns[*]}" -v col_widths="${col_widths[*]}" -v max_line=$max_line '
BEGIN {
    split(cut_cols, cuts, " ")
    split(col_widths, widths, " ")
}
NR > 1 {
    color_index = 0
    line = ""
    for (i=1; i<=NF; i++) {
        if (!(i in cuts)) {
            text = sprintf("%-" widths[i-1] "s", $i)
            part = (color_index % 2 == 0 ? "\033[47m\033[30m" text "\033[0m" : "\033[107m\033[30m" text "\033[0m")
            color_index++
            if (length(line part) > max_line) {
                print line
                line = part " "
            } else {
                line = line part " "
            }
        }
    }
    print line
}
' "$file"

# Clean up
rm col_widths.tmp
