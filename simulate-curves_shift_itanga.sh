#!/bin/bash
spowtd simulate rise D_Itanga001_55_GS_smooth_shift.sqlite3 curves_parameters_poros_shift_itanga.yml --observations -o curves_observations_poros_shift_itanga.yml && spowtd simulate recession D_Itanga001_55_GS_smooth_shift.sqlite3 curves_parameters_poros_shift_itanga.yml --observations >> curves_observations_poros_shift_itanga.yml

if [ $? -ne 0 ]
then
  exit 1
fi
