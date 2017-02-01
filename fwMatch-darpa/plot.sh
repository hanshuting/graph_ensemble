#!/bin/bash
# TODO: Error check

matlab -nodesktop -nosplash -r "plotExpt('${1}', ${2}); exit;"

