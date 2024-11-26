#!/bin/env python3

python ../../../kidsqmmm.py -i QMMM.in.OPENMM_MNDO -mm=openmm -qm=mndo -c=state.xml \
   	-l=layer.info -xmlff=real.xml,model-H.xml -pdb=real.pdb,model-H.pdb


