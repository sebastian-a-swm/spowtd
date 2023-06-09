#!/bin/bash
spowtd simulate rise D_Ekolongouma001_55_GS_smooth_shift.sqlite3 curves_parameters_poros_shift.yml --observations -o curves_observations_poros_shift.yml && spowtd simulate recession D_Ekolongouma001_55_GS_smooth_shift.sqlite3 curves_parameters_poros_shift.yml --observations >> curves_observations_poros_shift.yml

if [ $? -ne 0 ]
then
  exit 1
fi
