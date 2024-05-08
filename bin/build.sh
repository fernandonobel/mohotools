#!/bin/bash

VERSION=$(yq '.version' matlab_project.yaml)
DISTRIBUTION_NAME="sdist/mohotools-v$VERSION.zip"

if [ -f "$DISTRIBUTION_NAME" ]; then
    echo "Error: $DISTRIBUTION_NAME already exists."
    exit
fi

zip -r $DISTRIBUTION_NAME "+mohotools" examples CHANGELOG.md LICENSE matlab_project.yaml README.md 
