#!/bin/bash
OUTPUT_PATH=.
INPUT_PATH=..
MODULE_NAME=FenicsCode

epydoc --html -o $OUTPUT_PATH/html $1 $INPUT_PATH/$MODULE_NAME