#!/bin/bash

## # Initialization
## svn2git svn+ssh://clio.psy.vu.nl/home/r.pool/svn/BioCrates --verbose --username r.pool --exclude Run
## git remote add origin git@github.com:rpool/BioCrates.git
## git push origin master

# Sync
svn2git --rebase

exit
