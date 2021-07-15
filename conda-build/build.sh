#!/bin/sh
set -e

mkdir -p              "${PREFIX}/bin"
mv ./data             "${PREFIX}/"
mv ./bin/phibers.py   "${PREFIX}/bin/phibers"
mv ./bin/sr_a.joblib  "${PREFIX}/bin"
