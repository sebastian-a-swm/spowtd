#!/bin/bash
spowtd simulate rise E_Ekolongouma001_55_GS_smooth.sqlite3 spline_rise_parameters.yml -o spline_rise_observations.yml --observations

if [ $? -ne 0 ]
then
  exit 1
fi
