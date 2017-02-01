#!/bin/bash
# TODO: Error check

matlab -nodesktop -nosplash -r "runExpt('${1}', ${2}); plotExpt('${1}', ${2}); exit;"

