#!/bin/bash

# check out argument numbers
if [ "$1" == "--help" -o "$1" == "-h" ]
then
echo "Hi there!"
echo "This script is designed to help you play with hammurabi X."
echo "author: Jiaxin Wang, emal: jiaxin.wang@situ.edu.cn"
echo "To execute hammurabi X you can either do:"
echo "run.sh [XML parameter file path]"
echo "or do"
echo "hamx [XML parameter file path]"
echo "an XML template file can be found in the templates directory."
exit 0
fi

if [ $# != 1 ]
then
echo "wrong $# input(s)!"
echo "hammurabi X requires the path to the XML parameter file."
echo "try ./run.sh -h for more details."
exit 1
else
if [[ -f "$1" ]]
then
hamx $1
else
echo "could not find $1"
fi
fi
