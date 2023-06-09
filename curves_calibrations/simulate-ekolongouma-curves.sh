#!/bin/bash
spowtd simulate rise Ekolongouma1.5-5.sqlite3 /data/leuven/324/vsc32460/AC/spowtd/curves_calibrations/curves_parameters.yml --observations -o /data/leuven/324/vsc32460/AC/spowtd/curves_calibrations/Ekolongouma1.5-5/Ekolongouma1.5-5_curves_observations.yml && spowtd simulate recession Ekolongouma1.5-5.sqlite3 /data/leuven/324/vsc32460/AC/spowtd/curves_calibrations/curves_parameters.yml --observations >> /data/leuven/324/vsc32460/AC/spowtd/curves_calibrations/Ekolongouma1.5-5/Ekolongouma1.5-5_curves_observations.yml

if [ $? -ne 0 ]
then
  exit 1
fi
