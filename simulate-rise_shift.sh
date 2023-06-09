#!/bin/bash
spowtd simulate rise Ekolongouma001_55_GS_smooth_shift.sqlite3 rise_parameters_poros_shift.yml -o rise_observations_poros_shift.yml --observations

if [ $? -ne 0 ]
then
  exit 1
fi
