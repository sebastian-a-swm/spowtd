#!/bin/bash
spowtd simulate rise Q_Itanga001_55_GS_smooth_shift.sqlite3 rise_parameters_poros_shift_itanga.yml -o rise_observations_poros_shift_itanga.yml --observations

if [ $? -ne 0 ]
then
  exit 1
fi
