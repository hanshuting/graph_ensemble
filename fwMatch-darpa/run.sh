#!/bin/bash
# TODO: Error check

matlab -nodisplay -nodesktop -nosplash -r "runExpt('${1}', ${2});exit"

