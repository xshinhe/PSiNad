#!/bin/env python3

python ../../../kidsqmmm.py -i QMMM.in.GROMACS_MNDO -mm=gromacs -qm=mndo \
	-gro=real.gro,model-H.gro -top=real.top,model-H.top -c=real.crd

