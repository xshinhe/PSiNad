#!/bin/bash

############### CUSTOMIZE ###############
boxL=10
#########################################

QM_xyzfile=$1
flag=$2

if [ "$QM_xyzfile" == "" ]; then
	echo "ERROR: you should provide a XYZ file for QM part"
	echo "       like:"
	echo "       bash this_script.sh my_qm_sys.xyz"
	exit
fi

if [ "$flag" == "" ]; then
	flag=default
elif [ "$flag" == "ignore" ]; then
	echo "WARNING: You will ignore some check during automatically generating of files"
	echo "       And you should be responsible for all results"
	exit
fi

QM_dir=`echo $QM_xyzfile | sed 's/.xyz/_pro/g'`

if [ ! -d $QM_dir ]; then
	echo "# CREATE QM DIR: ${QM_dir}"
	mkdir $QM_dir
	echo "# COPYING QM XYZ: ${QM_xyzfile} as QM.xyz"
	cp $QM_xyzfile $QM_dir/QM.xyz
fi

cd ${QM_dir}
echo "# ENTERING QM DIR: $(pwd)"

natom_qm=$(head -n 1 QM.xyz)
echo "HERE ARE $natom_qm QM ATOMS"

if [ ! -f QM.mol2 ]; then
	echo "# TRYING TO CONVERT AS MOL2 FORMAT"
	obabel_exe=`which obabel`
	if [ "$obabel_exe" == "" ]; then 
		echo "ERROR: you should have a version of obabel"
		echo "       please install it by apt or something else"
		exit
	fi
	# build QM mol 
	obabel QM.xyz -O QM.mol2
fi 

if [ ! -f QM_ch.mol2 -o ! -f QM.frcmod ]; then
	echo "# TRYING TO CONVERT AS MOL2 FORMAT WITH BCC CHARGES BY GAFF"
	if [ ! -f "checked_mol2" -a ! "$flag" == "ignore" ]; then
		touch checked_mol2
		echo "WARNING: STOP HERE. You have to check mol2. Continue by repeating this script "
		exit
	fi
	antechamber_exe=`which antechamber`
	if [ "$antechamber_exe" == "" ]; then 
		echo "ERROR: you should have a version of anteamber"
		exit
	fi
	antechamber -i QM.mol2 -fi mol2 -o QM_ch.mol2 -fo mol2 -c bcc -at gaff
	parmchk -a Y -i QM_ch.mol2 -f mol2 -o QM.frcmod
fi 

if [ ! -f QM_solv.top ]; then
	echo "# TRYING TO BUILD TOPOLOGY BY TLEAP USING YOUR CUSTOMIZED CHAGRES"
	if [ ! -f "checked_ch_mol2" -a ! "$flag" == "ignore" ]; then
		touch checked_ch_mol2
		echo "WARNING: STOP HERE. You have to customize charges. Continue by repeating this script "
		exit
	fi
    echo """
# load parameters
source leaprc.water.tip3p
source leaprc.ff14SB
source leaprc.gaff

#
# load solute molecule
loadamberparams QM.frcmod
QM = loadmol2 QM_ch.mol2
list
check QM
saveoff QM QM.lib

# solvate solute in water (TIP3P)
# 14 is the minimal distance (in Angstroms) between
# the solute and the edge of the box
solvatebox QM TIP3PBOX ${boxL}

# save structure
saveamberparm QM QM_solv.top QM_solv.crd
quit
""" > tleap.inp
	echo "# TRYING TO BUILD TOPOLOGY BY TLEAP USING YOUR CUSTOMIZED CHAGRES"
	tleap -f tleap.inp # see details
    ambpdb -p QM_solv.top < QM_solv.crd > QM_solv.pdb
fi 


natom_total=`cat QM_solv.pdb | grep ATOM | wc -l`
natom_mm=$(($natom_total-$natom_qm))
echo "# TOTAL ATOMS = ${natom_total} AFTER SOLVATION"
if [ "${natom_total}" == "${natom_qm}" ]; then
	echo "ERROR: SOLVATION FALIED"
	rm QM_solv.top
	exit
fi

if [ ! -f min.info ]; then
	echo "# TRYING TO DO ENERGY MINIMIZATION"
	echo "# THE QUAMTUM PART WILL BE FROZEN DURING OPTIMIZATION"

	echo """
Minimize            ! line for comment
&cntrl              ! start the control sequence
   imin=1,          ! do minimization
   ntx=1,           ! read initial coordinates from filename
   irest=0,         ! that’s not a restart
   maxcyc=4000,     ! in total 4000 opt. cycles
   ncyc=2000,       ! last 2000 being Newton-Raphson steps
   ntpr=100,        ! Print energy every 100 steps
   ntwx=0,          ! No mdcrd file will be written
   cut=$(($boxL-2)),        ! cutoff for vdw/electrost. interactions (no more than boxL)
   ntb=1,           ! use periodic boundaries
   ntxo=1,          ! write output coordinates in ascii format
   ibelly=1,        ! use belly to frozen qm atoms
   bellymask=\"@$(($natom_qm+1))-$natom_total\", ! only mm atoms are movable
/                   ! end of control sequence
""" > min.in
	# minimalization
	sander -O -i min.in -o min.out -p QM_solv.top -c QM_solv.crd -r min.rst -inf min.info
	process_minout.perl min.out
	mkdir summ_min
	mv summary* summ_min
fi 

#exit

if [ -f min.rst -a ! -f heat.info ]; then
	echo "# TRYING TO HEAT"
	echo "# THE QUAMTUM PART WILL BE FROZEN DURING HEATING"

	echo """
Heat                     ! line for comment
&cntrl                   ! start the control sequence
   imin=0,               ! no minimization
   ntx=1,                ! read initial coordinates from file
   irest=0,              ! that’s not a restart
   nstlim=10000,         ! in total 10000 MD steps
   dt=0.002,             ! with a time step of 2 fs
   ntf=2,                ! ommit bond interactions involving H atoms
   ntc=2,                ! switch on SHAKE for hydrogens
   tempi=0.0,            ! initial temperature
   temp0=300.0,          ! target temperature
   ntpr=100,             ! print energy every 100 steps
   ntwx=100,             ! write coord. every 100 steps to mdcrd
   cut=$(($boxL-2)),        ! cutoff for vdw/electrost. interactions (no more than boxL)
   ntb=1,                ! use periodic boundaries
   !ntp=0,                ! No mdcrd file will be written
   ntt=5,                ! Use NOSE-HOOVER-CHAIN dynamics
   gamma_ln=2.0,         ! collision frequency for Langevin dynamics
   nmropt=1,             ! read NMR restraints/weights
   ig=-1,                ! random seed
   ntxo=1,               ! write output coordinates in ascii format
   ibelly=1,        ! use belly to frozen qm atoms
   bellymask=\"@$(($natom_qm+1))-$natom_total\", ! only mm atoms are movable
/                        ! end of control sequence
&wt type='TEMP0', istep1=0, istep2=8000,      ! bring the system from 0 to 300K in steps 0 to
value1=0.0, value2=300.0 /                    ! 9000
&wt type='TEMP0', istep1=8001, istep2=10000,  ! stay at 300 K for the next 1000 steps
value1=300.0, value2=300.0 /                  !
&wt type='END' /                              ! end of temperature input
""" > heat.in
	# heat
	sander -O -i heat.in -o heat.out -p QM_solv.top -c min.rst -r heat.rst -x heat.mdcrd -inf heat.info
	process_mdout.perl heat.out
	mkdir summ_heat
	mv summary* summ_heat
	
	echo """
parm QM_solv.top
trajin heat.rst
autoimage
trajout heat.rst7
go
quit
""" > boxize.in
	cpptraj -i boxize.in
	ambpdb -p QM_solv.top < heat.rst7 > heat.pdb
fi 

if [ -f heat.rst -a ! -f npt.info ]; then
	echo "# TRYING TO DO PRESSURE EQUILIBRIUM"
	echo "# THE QUAMTUM PART WILL BE FROZEN DURING PRESSURE EQUILIBRIUM"

	echo """
NPT SAMPLING
&cntrl
   imin=0,
   ntx=5,
   irest=1,
   nstlim=30000,
   dt=0.002,
   ntf=2,
   ntc=2,
   temp0=300.0,
   ntpr=100,
   ntwx=100,
   cut=$(($boxL-2)),        ! cutoff for vdw/electrost. interactions (no more than boxL)
   ntb=2, ! npt
   ntp=1, ! isobaric
   barostat=2, ! (1) berendsen (2) MonteCarlo
   ntt=5,                ! Use NOSE-HOOVER-CHAIN dynamics
   gamma_ln=2.0,
   ig=-1,
   ibelly=1,        ! use belly to frozen qm atoms
   bellymask=\"@$(($natom_qm+1))-$natom_total\", ! only mm atoms are movable
/                   ! end of control sequence
""" > npt.in
	# NPT ENSEMBLES
	sander -O -i npt.in -o npt.out -p QM_solv.top -c heat.rst -r npt.rst -x npt.mdcrd -inf npt.info
	process_mdout.perl npt.out
	mkdir summ_npt
	mv summary* summ_npt

	echo """
parm QM_solv.top
trajin npt.rst
autoimage
trajout npt.rst7
go
quit
""" > boxize.in
	cpptraj -i boxize.in
	ambpdb -p QM_solv.top < npt.rst7 > npt.pdb
fi


