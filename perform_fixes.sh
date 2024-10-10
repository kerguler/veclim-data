#!/bin/bash

source veclim-data-server/.env
cp veclim-fixes/cbook__init__.py ${DIR_PYTHON}/matplotlib/cbook/__init__.py
cp veclim-fixes/geoaxes.py ${DIR_PYTHON}/cartopy/mpl/geoaxes.py