#!/bin/bash
spowtd simulate rise E_Itanga001_55_GS_smooth.sqlite3 spline_rise_parameters_itanga.yml -o spline_rise_observations_itanga.yml --observations

if [ $? -ne 0 ]
then
  exit 1
fi
