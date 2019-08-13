#!/usr/bin/env bash

git lfs track lib/data/genetic_map_b37_withX.txt.gz
git add lib/data/genetic_map_b37_withX.txt.gz
git add .gitattributes
git commit -m 'lfs'

git remote add origin_gh_norm https://github.com/scimerc/genimpute.git
git push origin_gh_norm

