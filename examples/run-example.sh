#!/bin/bash

# Check if an argument is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <component-name>"
    exit 1
fi

# Build the component
cabal build $1

# Check if the build was successful
if [ $? -eq 0 ]; then
    # Execute the component
    cabal exec $1
else
    echo "Failed to build $1"
    exit 1
fi

