#!/bin/bash

# Runs with scattering
mkdir scat
rm scat/*
$OPENSN_INVERSE -i forward_ref.lua --lua c=0.5
$OPENSN_INVERSE -i constant.lua --lua c=0.5
$OPENSN_INVERSE -i constant_perturb.lua --lua c=0.5
$OPENSN_INVERSE -i adjoint_run.lua --lua c=0.5

# Runs without scattering
mkdir abs
rm abs/*
$OPENSN_INVERSE -i forward_ref.lua --lua c=0.0
$OPENSN_INVERSE -i constant.lua --lua c=0.0
$OPENSN_INVERSE -i constant_perturb.lua --lua c=0.0
$OPENSN_INVERSE -i adjoint_run.lua --lua c=0.0

python ndse_1dtest.py scat
python ndse_1dtest.py abs