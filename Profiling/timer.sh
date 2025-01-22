#!/bin/bash

# Output file
output_file="cpu_times_new_solver_3.txt"

# Clear the output file if it exists
# > "$output_file"

# Run the command 5 times and extract the line with `Process 0 - CPU Time:`
for i in {1..10}
do
  ./build/SWE-MPI-Runner -x 500 -y 500 | grep "Process 0 - CPU Time:" >> "$output_file"
done