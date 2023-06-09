#!/bin/bash
spowtd simulate rise ekolongouma.sqlite3 curves_pars.yml --observations -o ekolongouma_curves_observations.yml && spowtd simulate recession ekolongouma.sqlite3 curves_pars.yml --observations >> ekolongouma_curves_observations.yml

if [ $? -ne 0 ]
then
  exit 1
fi
