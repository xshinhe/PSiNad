#!/usr/bin/perl -w

my ($BASEDIR,$input,@PDB,@RES,@high_layer,@medium_layer,@medium_layer1,@low_layer,$type,@hmlinks,@hllinks,@mlinks,@llinks);
my ($resid,@UNIQUE_RES,@UNIQUE_RES_NUM,$solvent,$solvent_name,@solvent_name,@solvent,$dist,@link_number,$AMBERHOME,$BABELHOME);
my (@nonstd_resid,$nonstd_resid,@std_resid,@ff,@predef_resid,$predef_resid,@model_resid,%PTE,$solshell,@REASSIGNED_RES);
my ($dirname,$iGeom,$nGeom,$rand,$rst,$nAt,@at,@vel,@xyz,%HL,%ML,%HML,%rHML,@HML,@LL,$au2ang,$amu2e,$au2eV,$nproc,%crd_2_pdb);
my ($mass,@cm,@xyz0,$mmmd,$MMswitch,$top,$eqrst,$numAt,$eqpdb,$amb2au,$gau2au,$switch,@resid_size,$heat,@equivalence,$print);
my (%aminoacids,%N_aminoacids,%C_aminoacids,@ALL_RESIDUES,$mode,$GAUSS_EXE,$GAUSS_DIR,$sampling,$numModes,$geomdir,$temp);
my ($MPIrun,$MPI,$new_resid,@alphabet,@sol_med_pdb2real,%pdb_2_crd,$constrain,$term_resid_counter,$inact_modes,$COBRAMM_PATH);

$AMBERHOME=$ENV{'AMBERHOME'};
$GAUSS_DIR=$ENV{'GAUSSIAN_DIR'}.'/'.$ENV{'GAUSSIAN_EXE'};
$GAUSS_EXE=$ENV{'GAUSSIAN_DIR'}.'/'.$ENV{'GAUSSIAN_EXE'}.'/'.$ENV{'GAUSSIAN_EXE'};
$COBRAMM_PATH=$ENV{'COBRAMM_PATH'};
$MPIHOME=$ENV{'MPIHOME'};

#Explicit exports
#$AMBERHOME="/usr/local/amber18"
#$MPIHOME="/usr/local/amber18/AmberTools/src/mpich2-1.4/instlld/bin"
#$GAUSS_EXE="/usr/local/g16/g16"
#$GAUSS_DIR="/usr/local/g16"
#$COBRAMM_PATH="/usr/local/cobramm-git"

$nproc = 1;
$_ = `echo "\$(( \$(lscpu | awk '/^Socket\\(s\\)/{ print \$2 }') * \$(lscpu | awk '/^Core\\(s\\) per socket/{ print \$4 }') ))"`;
if ($_ == 1 || $_ == 2){
        $nproc = 1;
	$MPIHOME = '';
} else {
        $nproc = $_ - 1;
}

if (defined $MPIHOME && $MPIHOME ne ''){
	$MPIrun = "$MPIHOME/mpirun -np $nproc";
	$MPI = '.MPI';
} else {
	$MPIrun = '';
	$MPI = '';
	$nproc = 1;
}

$BASEDIR    = &getcwd();
$input      = "$BASEDIR/cobram.parm";
use strict;
use warnings;
no warnings 'recursion';
use POSIX qw/floor/;
use POSIX qw/ceil/;
use Cwd;
use Env;
use feature 'switch';
use File::Copy;
use Array::Utils qw':all';
use Term::ANSIColor;
@high_layer = ();
@medium_layer = ();
@low_layer = ();
@solvent_name = ('WAT');
$solvent_name = 'WAT';
$solvent = 0.0;
$solshell = 'CoMH';
$eqrst = 'equilibration.rst';
$au2ang   = 0.529177249;
$amb2au   = 0.000935;
$gau2au   = 2.41889E-17;
$amu2e    = 1822.88853;
$au2eV    = 27.211324570273;
$MMswitch = 0;
$switch   = 0;
$heat = 0;
$print = 'yes';
$constrain = read_input("High and Medium layer MD constraints");
if ($constrain eq 'harmonic' || $constrain eq '0'){
	$constrain = 0;
} elsif ($constrain eq 'ibelly'){
	$constrain = 1;
} else {
	print color("red"), "Unrecognized value for constraints ", color("reset"), $constrain, "\n";
        print "Should be either 'harmonic' or 'ibelly'";
	exit;
}
$term_resid_counter = 0;
@alphabet = qw/A B C D E F G H I J K L M N O P Q R S T U V W X Y Z/;
%PTE = (H  => [1 ,  1.007825], He => [2 ,  4.002603], Li => [3 ,   7.016004], Be => [4 ,  9.012182], B  => [5 , 11.009305], C  => [6 , 12.000000], N  => [7 , 14.003074], O  => [8 , 15.994914], F  => [9 , 18.998403], Ne => [10, 19.992440], Na => [11, 22.989769], Mg => [12, 23.985041], Al => [13, 26.981538], Si => [14, 27.976926], P  => [15, 30.973761], S  => [16, 31.972071], Cl => [17, 34.968852], Ar => [18, 39.962383], K => [19, 38.963706], Ca => [20,  39.962590], Sc => [21 ,44.9559], Ti => [22, 47.9479], V => [23, 50.9439], Cr => [24, 51.94051], Mn => [25, 54.938050], Fe => [26, 55.93494], Co => [27, 58.9332], Ni => [28, 57.93534], Cu => [29, 62.92960], Zn => [30, 63.92914], Ga => [31, 70.92470], Ge => [32, 73.92117], As => [33, 74.92159], Se => [34, 79.91652], Br => [35, 78.91833], Kr => [36, 83.91150], Ru => [37, 84.91178], Sr => [38, 87.90561], Y => [39, 88.90584], Zr => [40, 89.90470], Nb => [41, 92.90637], Mo => [42, 97.90540], Tc => [43, 97.907216], Ru => [44, 101.90435], Rh => [45, 102.90550], Pd => [46, 105.90348], Ag => [47, 106.90509], Cd => [48, 113.90335], In => [49, 114.90387], Sn => [50, 119.90219], Sb => [51, 120.90381], Te => [52, 129.90622], I => [53, 126.90446], Xe => [54, 131.9042], W => [74, 183.95093], Re => [75, 186.95575], Os => [76, 191.96147], Ir => [77, 192.96292], Pt => [78, 194.96477], Au => [79, 196.96655], Hg => [80, 201.97062], Tl => [81, 204.97441], Pb => [82, 207.97663], Bi => [83, 208.98038],);
%aminoacids =   (ALA => 10, ARG => 24, ASH => 13, ASN => 14, ASP => 12, CYM => 10, CYS => 11, CYX => 10, GLH => 16, GLN => 17, GLU => 15, GLY => 7, HID => 17, HIE => 17, HIP => 18, ILE => 19, LEU => 19, LYN => 21, LYS => 22, MET => 17, PHE => 20, PRO => 14, SER => 11, THR => 14, TRP => 24, TYR => 21, VAL => 16);
%N_aminoacids = (ALA => 12, ARG => 26, ASH =>  0, ASN => 16, ASP => 14, CYM =>  0, CYS => 13, CYX => 12, GLH =>  0, GLN => 19, GLU => 17, GLY => 9, HID => 19, HIE => 19, HIP => 20, ILE => 21, LEU => 21, LYN =>  0, LYS => 24, MET => 19, PHE => 22, PRO => 16, SER => 13, THR => 16, TRP => 26, TYR => 23, VAL => 18);
%C_aminoacids = (ALA => 11, ARG => 25, ASH =>  0, ASN => 15, ASP => 13, CYM =>  0, CYS => 12, CYX => 11, GLH =>  0, GLN => 18, GLU => 16, GLY => 8, HID => 18, HIE => 18, HIP => 19, ILE => 20, LEU => 20, LYN =>  0, LYS => 23, MET => 18, PHE => 21, PRO => 15, SER => 12, THR => 15, TRP => 25, TYR => 22, VAL => 17);

#####################################
############### Main ################
#####################################

check_files();
if ($switch == 1 || $switch == 3){
	my $ff = read_input("force fields");
	@ff = split(/\s+/,$ff);
	read_equivalent();
        $predef_resid = read_input("pre-defined non-standard residues");
	$nonstd_resid = read_input("non-standard residues");
        if ($predef_resid ne '0' || $nonstd_resid ne '0'){
                prepare_nonstd();
        }
	check_pdb(); 
	read_coord();
	get_high_layer(); 
	get_medium_layer(); 
	get_low_layer();
	sol_specs();
	get_links();
	prepare_real_layers();
	prepare_model();
	prepare_real();
	check_atom_types();
	real_vs_libs();
        if ($switch == 3){
	  edit_velocity();
        }
	if ($solvent_name ne "WAT" && $solvent_name ne "CL3"){
		print color("yellow"), "Rattle implemented only for TIP3P water (WAT) and Chloroform (CL3) residues.\n", color("reset");
	} elsif (($solvent_name eq "WAT" || $solvent_name eq "CL3") && intersect(@solvent,@medium_layer)) {
		rattle_h2o_cl3();
	}
} elsif ($switch == 2) {
        my $nGeoms = read_input("number of geometries");
	$sampling = read_input("sampling method");
	if ($sampling eq '0' || $sampling eq "Wigner"){
		$sampling = "Wigner";
		print "Initial conditions will be generated by Wigner sampling.\n"; 
		$numModes=get_NM();
		print "Number of normal modes of the High(+Medium) Layer = $numModes\n";
	} elsif ($sampling eq 'ZPE') {
		print "Initial conditions will be generated by zero point energy (ZPE) sampling as available in Gaussian.\n";
	} else {
		print color("red"), "Please specify whether Wigner or ZPE sampling should be performed!\n", color("reset");	
		exit;
	}
	#not necessary when force constants and masses are read from .fchk
	#if ($sampling eq "Wigner" && -e "$BASEDIR/geometry.chk" && !-e "$BASEDIR/geometry.log"){
	#		generate_log();
	#}
	#if ($sampling eq "Wigner"){
	#		generate_molf();
	#}
	for $iGeom (1..$nGeoms){
	        print "\n###########################\nProceeding geometry $iGeom ...\n\n";
	        $dirname = sprintf("%s%s%03d", $BASEDIR, '/geom_', $iGeom);
	        mkdir $dirname;
	        chdir $dirname;
	        gen_init_cond($iGeom);
                if ($mode eq 'with_environment'){
	                replace_MM_HL($iGeom);
	                $eqpdb = sprintf("%s%03d%s", 'geom_', $iGeom, '.pdb');
	 	        if ($heat == 1){
	 	       		run_MM_heat();
	 	        }
	                run_MM_constr();
	                create_velocity();
                        prepare_switch3();
		} elsif ($mode eq 'gas-phase' && $sampling eq "Wigner") {
			read_wigner();
		} elsif ($mode eq 'gas-phase' && $sampling eq "ZPE") {
                        read_gaussian();
		}
	        chdir $BASEDIR;
	}
	if ($sampling eq "Wigner"){
		analyze_sampling();
	}
}
keep_files();

#####################################
########### Subroutines #############
#####################################

sub check_files{
if (!-e "$BASEDIR/cobram.parm"){
	print color("red"), "cobram.parm is missing!\n", color("reset");
        exit;
}
my $coord = read_input("coordinate file");
my $nGeoms = read_input("number of geometries");
my $rst = read_input("MM geometry/velocity file");
if ( length $coord && $coord ne 0 && -e "$BASEDIR/$coord" && ! -e "$BASEDIR/velocity.dat"){
	print "Input files will be generated for running optimizations with Cobram.\n";
	$switch = 1;
} elsif ($nGeoms != 0 && ! -e "$BASEDIR/velocity.dat" && ! -e "$BASEDIR/geom_*.pdb"){
	print "Initial conditions for running dynamics with Cobram will be generated.\n";
	$switch = 2;
        if (-e "$BASEDIR/$rst" || (-e "$BASEDIR/real_layers.xyz" && (`grep -c " L" $BASEDIR/real_layers.xyz` != 0 ))){
                $mode = "with_environment";
		$_ = read_input("MM input file");
		if ($_ eq '0'){
                	print color("red"), "Input file for dynamics is missing!", color("reset");
			exit;
		}
		$_ = read_input("MM topology file");
		if ($_ eq '0'){
		        print color("red"), "MM topology file is missing!", color("reset");
		        exit;
		}
                if (-e "$BASEDIR/$rst"){
                        length_of_rst($rst);
                } else {
                        $heat = 1;
                        if (!-e "$BASEDIR/real_layers.xyz"){
                                print color("red"), "MD heating and equilibration is required but real_layers.xyz is missing!\n", color("reset");
                                exit;
                        }
                        real_layers_2_crd();
                }
                if ($heat == 0 && !-e "$BASEDIR/crd_2_HM.map"){
                        print color("red"), "MD equilibration is required but crd_2_HM.map is missing! Mapping of the sampled HL/ML geometry on the MD .rst not possible\n", color("reset");
                        exit;
                }
        } else {
		$mode = "gas-phase";
	}
        $sampling = read_input("sampling method");
	if ($sampling eq "Wigner"  && !-e "$BASEDIR/geometry.fchk" && !-e "$BASEDIR/geometry.chk" && !-e "$BASEDIR/geometry.chk.gz"){
		print color("red"), "geometry.chk from a previous frequency calculation is missing!\n", color("reset");
		exit;
	} elsif ($sampling eq "ZPE" && !-e "$BASEDIR/geometry.chk" && !-e "$BASEDIR/geometry.chk.gz" && !-e "$BASEDIR/geometry.log"){
                print color("red"), "geometry.log or geometry.chk from a previous frequency calculation are missing.\n", color("reset");
                exit;
        } 
} elsif (-e "$BASEDIR/velocity.dat" && -e "$BASEDIR/$coord"){
	print "Input files for running dynamics with Cobram will be generated.\n";
	$switch = 3;
	open PDB, "<$BASEDIR/$coord" or die;
	my $numAt = 0;
	my $velAt = 0;
	while (<PDB>){
		if (/^\s*ATOM\s+/){
			$numAt++;
		}
	}
	close(PDB);
	open VEL, "<$BASEDIR/velocity.dat";
	while (<VEL>){
        	$velAt++;
        }
	close(VEL);
	if ($numAt != $velAt){
		print color("red"), "The velocity.dat file is incomplete ($velAt, should be $numAt). You are probably trying to repeat the preparation step. Copy a fresh velocity.dat file from the parent directory.\n", color("reset");
		exit;
	}
}
if ($switch == 0){
	print color("red"), "Cannot decide what to do! Possibly files missing, examine cobram.parm!\n", color("reset");
	exit; 
}
}

#######################################

sub length_of_rst{
my $count = 0;
my $rst = $_[0];
open IN, "<".$rst or die;
$_ = <IN>;
$_ = <IN>;
$_ =~ s/^\s*//;
$_ =~ s/\s*$//;
my ($numAt) = split(/\s+/,$_);
while (<IN>){
	$count++;
}
if ($count > int($numAt/2)+5){
	print "An $rst file with velocities is found. Only MD equilibration will be performed.\n";
} else {
	print color("yellow"), "An $rst file was found but without velocities. MD heating and equilibration will be performed.\n",color("reset");
	$heat = 1;
	if (!-e "$BASEDIR/real_layers.xyz"){
                print color("red"), "MD heating and equilibration is required but real_layers.xyz from the previous frequency calculation is missing!\n", color("reset");
                exit;
        }
	real_layers_2_crd();
}
close(IN);
}

#######################################

sub real_layers_2_crd{
$rst = "real.crd";
open XYZ, "<real_layers.xyz" or die;
my $count=0;
while(<XYZ>){
	$count++;
}
close(XYZ);
open XYZ, "<real_layers.xyz" or die;
open CRD, ">real.crd" or die;
print CRD "\n$count\n";
$count=0;
while(<XYZ>){
	$count++;
	$_ =~ s/^\s*//;
	$_ =~ s/\s*$//;
	@_ = split(/\s+/,$_);
	printf CRD "%12.7f%12.7f%12.7f",$_[2],$_[3],$_[4];
	if ($count % 2 == 0){
		$count=0;
		print CRD "\n";
	}
}
close(XYZ);
close(CRD);
}

#######################################

sub read_input{
my $key = $_[0];
open (INPX, "<".$input) or die;
while (<INPX>){
	if (/^$key/){
		$_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
                @_ = split(/=/,$_);
		if ($_[1]){
                	$_[1] =~ s/^\s*//;
			close(INPX);
                	return $_[1];
		} else {
			close(INPX);
			return 0;
		}
	}
}
return 0;
}

#####################################

sub sol_specs{
$solvent = read_input("automatic inclusion of residues to High and Medium layer");
if ($solvent eq "0"){
        print "No solvent is selected\n";
} elsif ($solvent =~ /([\w{3}\s+]+)\s+(\d+[\.\d+]*)\s+(\w+)/) {
        @solvent_name = split(/\s+/,$1);
	$solvent_name = $solvent_name[0];
        if (scalar(@solvent_name) > 1 && $solvent_name[1] eq 'ALL'){
		@solvent_name = @ALL_RESIDUES;
	}
        $solvent = $2;
        $solshell = $3;
	print "Residues @solvent_name within $solvent from $solshell will be included.\n";
        include_solvent($solvent);
} elsif ($solvent =~ /([\w{3}\s+]+)\s+(\d+[\.\d+]*)/) {
        @solvent_name = split(/\s+/,$1);
	$solvent_name = $solvent_name[0];
        if (scalar(@solvent_name) > 1 && $solvent_name[1] eq 'ALL'){
                @solvent_name = @ALL_RESIDUES;
        }
        $solvent = $2;
	print "Residues @solvent_name within $solvent from $solshell will be included.\n";
        include_solvent($solvent);
} elsif ($solvent =~ /([\w{3}\s+]+)\s+(\w+)/){
        @solvent_name = split(/\s+/,$1);
	$solvent_name = $solvent_name[0];
	if ($2 eq 'H' || $2 eq 'M'){
        	$solshell = $2;
	} else {
		push(@solvent_name,$2);
		$solshell = 'M';
	}
        if (scalar(@solvent_name) > 1 && $solvent_name[1] eq 'ALL'){
                @solvent_name = @ALL_RESIDUES;
        }
	print "H-bonded residues @solvent_name will be included to $solshell layer.\n";
        get_H_bonds();
} elsif ($solvent =~ /([\w{3}\s+]+)/){
	@solvent_name = split(/\s+/,$1);
	$solvent_name = $solvent_name[0];
        $solshell = 'M';
        if (scalar(@solvent_name) > 1 && $solvent_name[1] eq 'ALL'){
                @solvent_name = @ALL_RESIDUES;
        }
	print "H-bonded residues @solvent_name will be included to M layer\n";
        get_H_bonds();
}
}

#####################################

sub read_coord{
my ($old_type,$col1,$col2,$col3,$col4,$col5,$col6,$col7,$col8,$col9,$col10,$col11,$index);
my $coord = read_input("coordinate file");
open (PDB, "<".$coord) or die "Input file $coord is missing!";
print "Reading $coord ... ";
$resid = 0;
my $kilo = 0;
my $nonstd_count = 0;
my $nonstd_num = 0;
$old_type = '';
while (<PDB>){
	if ($_ =~ /^\s*ATOM\s+/){
		@_ = split(//,$_);
		$col2 = join('',@_[6..10]);	#atom number
		$col2 =~ s/^\s*//;
                $col2 =~ s/\s*$//;
		$col3 = join('',@_[12..15]);	#atom name
		$col3 =~ s/^\s*//;
                $col3 =~ s/\s*$//;
		$col4 = $_[16];
		$col5 = join('',@_[17..19]);	#residue name
		$col5 =~ s/^\s*//;
                $col5 =~ s/\s*$//;
		$col6 = $_[21];
		$col7 = join('',@_[22..25]);	#residue number
		$col7 =~ s/^\s*//;
                $col7 =~ s/\s*$//;
		$col8 = $_[26];
		$col9 = join('',@_[30..37]);	#x-coord
		$col9 =~ s/^\s*//;
                $col9 =~ s/\s*$//;
		$col10 = join('',@_[38..45]);	#y-coord
		$col10 =~ s/^\s*//;
                $col10 =~ s/\s*$//;
                $col11 = join('',@_[46..53]);	#z-coord
		$col11 =~ s/^\s*//;
                $col11 =~ s/\s*$//;
		if ($col6 eq '@'){
                        $kilo = 10000;
                }
		if ($col2 == 1 || $old_type eq 'TER'){
			$old_type = $col7;
		}
                push(@{$PDB[0]},$col5);
                push(@{$PDB[1]},$col9);
                push(@{$PDB[2]},$col10);
                push(@{$PDB[3]},$col11);
                push(@ALL_RESIDUES,$col5);
		$col3 =~ s/\d+//g;
		$col3 =~ s/\'//g;
		my $col3_bkp = $col3;
		if ($col3 ne 'EPW'){
		  until (grep {$_ eq $col3} keys %PTE){
	                $col3 =~ s/.$//g;
			if ($col3 eq ''){
				print  color("red"), "Error: can't assign atom $col3_bkp; Note that atom names in the .pdb (and .lib files) file should be defined by first capital letter (e.g Na, Mg) to differentiate them from atom types (i.e. NA)\n", color("reset");
				exit;
			}
		  }
		} 
		push(@{$PDB[4]},$col3);
		push(@{$PDB[5]},$col7+$kilo);
		if ($col7 ne $old_type){		# $resid increases by one for each new residue
			$resid++;
		}
	        push(@{$RES[$resid]},$col2);
		$old_type = $col7;
        }
        if ($_ =~ /^\s*TER/){
                $resid++;
		$old_type = 'TER';
        }
}
close (PDB);
my %dummy = ();
@ALL_RESIDUES = grep { ! $dummy{ $_ }++ } @ALL_RESIDUES;
#print "Resid size: @resid_size\n";
print "done!\n";
}

#####################################

sub get_high_layer{
my @range;
$_ = read_input("high layer atoms");
if ($_ eq '0'){
	return;
}
@_ = split(/\s+/,$_);
for my $i (0..$#_){
        @range = split(/-/,$_[$i]);
	if (!$range[1]){
		push(@high_layer,$pdb_2_crd{$range[0]});	
	} else {
        	for my $j ($range[0]..$range[1]){
                	push(@high_layer,$pdb_2_crd{$j});
		}
        }
	undef @range;
}
@high_layer = sort {$a <=> $b} @high_layer;
}

######################################

sub get_medium_layer{
$_ = read_input("medium layer atoms");
if ($_ eq '0'){
        return;
}
@_ = split(/\s+/,$_);
for my $i (0..$#_){
        my @range = split(/-/,$_[$i]);
	if (!$range[1]){
                push(@medium_layer,$pdb_2_crd{$range[0]});
        } else {
        	for my $j ($range[0]..$range[1]){
	                push(@medium_layer,$pdb_2_crd{$j});
        	}
	}
}
@medium_layer = sort {$a <=> $b} @medium_layer;
@medium_layer1 = @medium_layer;
}

#####################################

sub get_low_layer{
for my $j (1..$#{$PDB[0]}+1){
        if ((grep {$_ == $j} @medium_layer) || (grep {$_ == $j} @high_layer)){
        } else {
                push(@low_layer,$j);
        }
}
}

#####################################

sub get_H_bonds{
my (@ilayer);
my @heteroatoms = ("N", "O", "F", "P", "S", "Cl", "Br");
my $dist;
print "Looking for H-bonds ... ";
if ($solshell eq 'H'){
	@ilayer = (@high_layer);
} elsif ($solshell eq 'M'){
	@ilayer = (@high_layer,@medium_layer);
}
foreach my $i (@ilayer){
	if (grep {$_ eq ${$PDB[4]}[$i-1]} @heteroatoms){
		foreach my $j (@low_layer){
			if (${$PDB[4]}[$j-1] eq "H" && (grep {$_ eq ${$PDB[0]}[$j-1]} @solvent_name)){
				$dist = calculate_dist($i,$j);
				if ($dist <= 2.5){
					foreach my $k (@low_layer){
						if ($k == $j){
						} else {
							$dist = calculate_dist($k,$j);
							if (($dist <= 1.20) && (grep {$_ eq ${$PDB[4]}[$k-1]} @heteroatoms)){
								#print "Heteroatom: $i(${$PDB[4]}[$i-1])\n";
								push(@solvent,complete_solvent($j));
							}
						}
					}
				}
			}
		}
	}
}
print "done for heteroatoms in the high/medium layer (Het(H/M)---H-Het(M/L))\n\n";
print "Looking for H-bonds ... ";
foreach my $i (@high_layer){
        if (${$PDB[4]}[$i-1] eq "H"){
		foreach my $k (@high_layer){
			if ($k == $i){
			} else {
				$dist = calculate_dist($i,$k);
				if (($dist <= 1.20) && (grep {$_ eq ${$PDB[4]}[$k-1]} @heteroatoms)){
					foreach my $j (@low_layer){
                 			        if ((grep {$_ eq ${$PDB[4]}[$j-1]} @heteroatoms) && (grep {$_ eq ${$PDB[0]}[$j-1]} @solvent_name)){
                                			$dist = calculate_dist($i,$j);
                                			if ($dist <= 2.5){
								print "Hydrogen: $i(${$PDB[4]}[$i-1]) ... $k(${$PDB[4]}[$k-1])\n";
								push(@solvent,complete_solvent($j));
                                			}
                        			}
                			}
				}
			}
                }
        }
}
if ($solshell eq 'M'){
	foreach my $i (@medium_layer){
	        if (${$PDB[4]}[$i-1] eq "H"){
	                foreach my $k (@medium_layer){
	                        if ($k == $i){
	                        } else {
	                                $dist = calculate_dist($i,$k);
	                                if (($dist <= 1.20) && (grep {$_ eq ${$PDB[4]}[$k-1]} @heteroatoms)){
	                                        foreach my $j (@low_layer){
	                                                if ((grep {$_ eq ${$PDB[4]}[$j-1]} @heteroatoms) && (grep {$_ eq ${$PDB[0]}[$j-1]} @solvent_name)){
	                                                        $dist = calculate_dist($i,$j);
	                                                        if ($dist <= 2.5){
									print "Hydrogen: $i(${$PDB[4]}[$i-1]) ... $k(${$PDB[4]}[$k-1])\n";
									push(@solvent,complete_solvent($j));
	                                                        }
	                                                }
	                                        }
	                                }
	                        }
	                }
	        }
	}
}
print "done for H's attached to heteroatoms in the high/medium layer (Het(H/M)-H(H/M)---Het(L)\n\n";
if ($solshell eq 'H'){
	move_to_high_layer(@solvent);
} elsif ($solshell eq 'M') {
	move_to_medium_layer(@solvent);
}
print "High layer atoms: ";
for my $i (@high_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";

my @tmp_medium_layer = @medium_layer;
@medium_layer = ();
foreach my $i (@tmp_medium_layer){
        if (grep {$_ == $i} @high_layer){
        } else {
                push(@medium_layer,$i)
        }
}


print "Medium layer atoms: ";
for my $i (@medium_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";
#print "Low layer atoms: ";
##for my $i (@low_layer){
##        print "$i(${$PDB[0]}[$i-1]) ";
##}
##print "\n\n";
}

#####################################

sub include_solvent {
my $r = $_[0];
if (($solshell eq 'CoMH') || ($solshell eq 'CoMHM')){
	my @R = center_of_mass($solshell);
	print "Center of mass: @R\n";
	push(@{$PDB[1]},$R[0]);
	push(@{$PDB[2]},$R[1]);
	push(@{$PDB[3]},$R[2]);
	foreach my $j (@low_layer){
		if (grep {$_ eq ${$PDB[0]}[$j-1]} @solvent_name){
			$dist = calculate_dist($#{$PDB[1]}+1,$j);
			if ($dist <= $r){
				#print "Include low layer atom: $j(${$PDB[4]}[$j-1])\n";
				push(@solvent,complete_solvent($j));
			}
		}
	}
} elsif ($solshell eq 'DtH'){
	foreach my $i (@high_layer){
	        foreach my $j (@low_layer){
			if (grep {$_ eq ${$PDB[0]}[$j-1]} @solvent_name){
			        $dist = calculate_dist($i,$j);
			        if ($dist <= $r){
				        #print "High layer atom: $i(${$PDB[4]}[$i-1]) ... Low layer atom: $j(${$PDB[4]}[$j-1])\n";
				        push(@solvent,complete_solvent($j));
		                }
			}
	        }
	}
} elsif ($solshell eq 'DtHM'){
        foreach my $i (@high_layer,@medium_layer){
                foreach my $j (@low_layer){
			if (grep {$_ eq ${$PDB[0]}[$j-1]} @solvent_name){
		                $dist = calculate_dist($i,$j);
		                if ($dist <= $r){
		                        #print "High/Medium layer atom: $i(${$PDB[4]}[$i-1]) ... Low layer atom: $j(${$PDB[4]}[$j-1])\n";
		                        push(@solvent,complete_solvent($j));
		                }
			}
                }
        }
}

for my $i (0..$#medium_layer){
  if (${$PDB[0]}[$medium_layer[$i]] eq 'WAT'){
    print "Atom $medium_layer[$i] is WAT\n";
    push(@solvent,$medium_layer[$i]);
  }
}

move_to_medium_layer(@solvent);
print "High layer atoms: ";
for my $i (@high_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";

my @tmp_medium_layer = @medium_layer;
@medium_layer = ();
foreach my $i (@tmp_medium_layer){
        if (grep {$_ == $i} @high_layer){
        } else {
		push(@medium_layer,$i)
	}
}

print "Medium layer atoms: ";
for my $i (@medium_layer){
        print "$i(${$PDB[0]}[$i-1]) ";
}
print "\n\n";
}

#####################################

sub center_of_mass{
my @R = 0;
my $M = 0;
if ($_[0] eq 'CoMH'){
	@_ = @high_layer;
} elsif ($_[0] eq 'CoMHM'){
	@_ = (@high_layer,@medium_layer);
}
foreach my $i (@_){
	$M += ${$PTE{${$PDB[4]}[$i-1]}}[1];
}
foreach my $i (@_){
	$R[0] += ${$PTE{${$PDB[4]}[$i-1]}}[1]*${$PDB[1]}[$i-1]/$M;
	$R[1] += ${$PTE{${$PDB[4]}[$i-1]}}[1]*${$PDB[2]}[$i-1]/$M;
	$R[2] += ${$PTE{${$PDB[4]}[$i-1]}}[1]*${$PDB[3]}[$i-1]/$M;
}
return @R;
}

#####################################

sub calculate_dist{
my $a = $_[0];
my $b = $_[1];
$dist = sqrt((${$PDB[1]}[$a-1] - ${$PDB[1]}[$b-1])**2 + (${$PDB[2]}[$a-1] - ${$PDB[2]}[$b-1])**2 + (${$PDB[3]}[$a-1] - ${$PDB[3]}[$b-1])**2);
return $dist;
}

#####################################

sub complete_solvent{
my $atom = $_[0];
my @molecule = ();
#print "Solvent atom: $atom(${$PDB[4]}[$atom-1])\n";
for my $i (0..$#RES){
	if (grep {$_ == $atom} @{$RES[$i]}){
		for my $j (@{$RES[$i]}){
			push(@molecule,$j);
		}
		last;
	}
}
#print "Solvent molecule: ";
#foreach my $i (@molecule){
#	print "$i(${$PDB[4]}[$i-1]) ";
#}
#print "\n";
return @molecule;
}

#####################################

sub move_to_high_layer{
my %dummy = ();
my @solvent = grep { ! $dummy{ $_ }++ } @_; # removes duplicated solvent molecules
#printf "Number of molecules that will be added to the high layer:%4d\n\n", ($#solvent+1)/($#{$RES[-1]}+1);
foreach my $i (@solvent){
        push(@high_layer,$i);
        my $count = 0;
        foreach my $j (@low_layer){
                if ($j == $i){
                        last;
                }
                $count++;
        }
        splice(@low_layer,$count,1);
}
@high_layer = sort {$a <=> $b} @high_layer;
@low_layer = sort {$a <=> $b} @low_layer;
}

#####################################

sub move_to_medium_layer{
my %dummy = ();
my @solvent = grep { ! $dummy{ $_ }++ } @_; # removes duplicated solvent molecules
#printf "Number of molecules that will be added to the medium layer:%4d\n\n", ($#solvent+1)/($#{$RES[-1]}+1);
foreach my $i (@solvent){
        if (!grep {$_ == $i} @medium_layer){
  	  push(@medium_layer,$i);
        }
	my $count = 0;
	foreach my $j (@low_layer){
		if ($j == $i){
			last;
		}
		$count++;
	}
	splice(@low_layer,$count,1);
}
my @medium_layer = grep { ! $dummy{ $_ }++ } @_; # removes duplicated solvent molecules
@medium_layer = sort {$a <=> $b} @medium_layer;
@low_layer = sort {$a <=> $b} @low_layer;
}

#####################################

sub get_links{
my $link_atoms = read_input("links");
if ($link_atoms eq '0'){
        return;
}
@_ = split(/\s+/,$link_atoms);
for my $i (0..$#_){
        my @b_atoms = split(/-/,$_[$i]);
        if (grep {$_ == $b_atoms[0]} @high_layer){
		if (grep {$_ == $b_atoms[1]} @medium_layer){
	                push(@hmlinks,$pdb_2_crd{$b_atoms[0]});
                	push(@mlinks,$pdb_2_crd{$b_atoms[1]});
		} elsif (grep {$_ == $b_atoms[1]} @low_layer){
			push(@hllinks,$pdb_2_crd{$b_atoms[0]});
                        push(@llinks,$pdb_2_crd{$b_atoms[1]});
		} else {
			print color("red"), "Link atom $pdb_2_crd{$b_atoms[0]} was found in the high layer but link atom $pdb_2_crd{$b_atoms[1]} was found neither in the medium layer nor in the low layer!\n", color("reset"); 
			exit;
		}
        } elsif (grep {$_ == $b_atoms[0]} @medium_layer){
		if (grep {$_ == $b_atoms[1]} @high_layer){
                        push(@hmlinks,$pdb_2_crd{$b_atoms[1]});
                        push(@mlinks,$pdb_2_crd{$b_atoms[0]});
		} else {
			print color("red"), "Link atom $pdb_2_crd{$b_atoms[0]} was found in the medium layer but link atom $pdb_2_crd{$b_atoms[0]} was not found in the high layer!\n", color("reset"); 
			exit;
		}
        } elsif (grep {$_ == $b_atoms[0]} @low_layer){
                if (grep {$_ == $b_atoms[1]} @high_layer){
                        push(@hllinks,$pdb_2_crd{$b_atoms[1]});
                        push(@llinks,$pdb_2_crd{$b_atoms[0]});
                } else {
                        print color("red"), "Link atom $pdb_2_crd{$b_atoms[0]} was found in the low layer but link atom $pdb_2_crd{$b_atoms[0]} was not found in the high layer!\n", color("reset"); 
			exit;
                }
	} else {
                print color("red"), "Link atom $pdb_2_crd{$b_atoms[0]} wasn't found in any layers!\n", color("reset"); 
		exit;
        }
}
if (@mlinks){
	print "Link atoms:    H        M\n";
	foreach my $i (0..$#mlinks){
		print "              $hmlinks[$i](${$PDB[4]}[$hmlinks[$i]-1])    $mlinks[$i](${$PDB[4]}[$mlinks[$i]-1])\n";
	}
}
if (@llinks){
	print "Link atoms:    H        L\n";
	foreach my $i (0..$#llinks){
	        print "              $hllinks[$i](${$PDB[4]}[$hllinks[$i]-1])    $llinks[$i](${$PDB[4]}[$llinks[$i]-1])\n";
	}
}
print "\n";
}

########################################

sub prepare_real_layers{
my (@dist,$r,@H_pos);
my $default;
$default->{'C'} = 1.089;
$default->{'N'} = 1.008;
$default->{'O'} = 0.947;

print "Generating the real_layers.xyz and performing H/M/L layer sorting ... ";
open (MAPH, ">crd_2_H.map") or die;
print MAPH "PDB atoms moved to the High layer (High Layer -> .crd):\n";
open (MAPM, ">crd_2_M.map") or die;
print MAPM "PDB atoms moved to the Medium layer (Medium Layer -> .crd):\n";
open (OUT, ">real_layers.xyz") or die;
open (XYZ, ">model-H.xyz") or die;
open (XYZ1, ">HighMediumLayer.xyz") or die;
open (PDB1, ">model-H1.pdb") or die;
printf XYZ "%6d\n\n",scalar(@high_layer)+scalar(@hmlinks)+scalar(@hllinks);
printf XYZ1 "%6d\n\n",scalar(@high_layer)+scalar(@medium_layer);
my $count = 0;
my $high_count = 0;
for my $i (1..$#{$PDB[0]}+1){
	if (grep {$_ == $i} @high_layer){
		 $count++;
                 $high_count++;
		 printf MAPH "%10d  %10d\n", $count, $crd_2_pdb{$i};
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  H\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
		 printf XYZ "%s    %12.6f  %12.6f  %12.6f\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
                 printf XYZ1 "%s    %12.6f  %12.6f  %12.6f\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
                 printf PDB1 "ATOM  %5d %3s  LIG     1    %8.3f%8.3f%8.3f\n", $high_count, ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
	} elsif ((grep {$_ == $i} @medium_layer) && (grep {$_ == $i} @mlinks)){
		 $count++;
		 printf MAPM "%10d  %10d\n", $count, $crd_2_pdb{$i};
                 my $k = 0;
                 while ($k <= $#mlinks){
 	               if ($i == $mlinks[$k]){
	                       last;
                       }
                       $k++;
                 }
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  M H  %10d\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1], $hmlinks[$k];
                 printf XYZ1 "%s    %12.6f  %12.6f  %12.6f\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
		 push(@dist,${$PDB[1]}[$hmlinks[$k]-1] - ${$PDB[1]}[$i-1]);
		 push(@dist,${$PDB[2]}[$hmlinks[$k]-1] - ${$PDB[2]}[$i-1]);
                 push(@dist,${$PDB[3]}[$hmlinks[$k]-1] - ${$PDB[3]}[$i-1]);
                 $r = sqrt($dist[0]**2 + $dist[1]**2 + $dist[2]**2);
		 if (grep {$_ eq ${$PDB[4]}[$hmlinks[$k]-1]} (keys %$default)){
		 	push(@H_pos,${$PDB[1]}[$hmlinks[$k]-1] - ($default->{${$PDB[4]}[$hmlinks[$k]-1]}/$r)*$dist[0]);
                 	push(@H_pos,${$PDB[2]}[$hmlinks[$k]-1] - ($default->{${$PDB[4]}[$hmlinks[$k]-1]}/$r)*$dist[1]);
                 	push(@H_pos,${$PDB[3]}[$hmlinks[$k]-1] - ($default->{${$PDB[4]}[$hmlinks[$k]-1]}/$r)*$dist[2]);
		 	@dist = ();
	 	 } else {
			print color("red"), "\nOne of the boundery atoms is not C,N or O. Currently, this situation is not implemented.\n", color("reset");
			exit;
		 }
		 if (${$PDB[0]}[$i-1] eq $solvent_name){
			push(@{$sol_med_pdb2real[0]},$i);
			push(@{$sol_med_pdb2real[1]},$count);
		 }
	} elsif (grep {$_ == $i} @medium_layer){
		 $count++;
		 printf MAPM "%10d  %10d\n", $count, $crd_2_pdb{$i};
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  M\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
                 printf XYZ1 "%s    %12.6f  %12.6f  %12.6f\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
		 if (${$PDB[0]}[$i-1] eq $solvent_name){
	                push(@{$sol_med_pdb2real[0]},$i);
                        push(@{$sol_med_pdb2real[1]},$count);
                 }
	} elsif ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] ne $solvent_name) && (grep {$_ == $i} @llinks)){
                 my $k = 0;
                 while ($k <= $#llinks){
                       if ($i == $llinks[$k]){
                               last;
                       }
                       $k++;
                 }
                 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  L H  %10d\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1], $hllinks[$k];
                 push(@dist,${$PDB[1]}[$hllinks[$k]-1] - ${$PDB[1]}[$i-1]);
                 push(@dist,${$PDB[2]}[$hllinks[$k]-1] - ${$PDB[2]}[$i-1]);
                 push(@dist,${$PDB[3]}[$hllinks[$k]-1] - ${$PDB[3]}[$i-1]);
                 $r = sqrt($dist[0]**2 + $dist[1]**2 + $dist[2]**2);
                 if (grep {$_ eq ${$PDB[4]}[$hllinks[$k]-1]} (keys %$default)){
                        push(@H_pos,${$PDB[1]}[$hllinks[$k]-1] - ($default->{${$PDB[4]}[$hllinks[$k]-1]}/$r)*$dist[0]);
                        push(@H_pos,${$PDB[2]}[$hllinks[$k]-1] - ($default->{${$PDB[4]}[$hllinks[$k]-1]}/$r)*$dist[1]);
                        push(@H_pos,${$PDB[3]}[$hllinks[$k]-1] - ($default->{${$PDB[4]}[$hllinks[$k]-1]}/$r)*$dist[2]);
                        @dist = ();
                 } else {
                        print color("red"), "\nOne of the boundery atoms is not C,N or O. Currently, this situation is not implemented.\n", color("reset");
                        exit;
                 } 
	} elsif ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] ne $solvent_name)){
                 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  L\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
	} elsif ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] eq $solvent_name)){
		 next;
        } else {
                 print color("red"), "\nAtom $i${$PDB[4]}[$i-1] was not found in any layer!\n", color("reset");
		 exit;
        }
}
print MAPH "\n";
print MAPM "\n";
close(MAPH);
close(MAPM);
system("cat crd_2_H.map crd_2_M.map > crd_2_HM.map; rm crd_2_H.map crd_2_M.map");
for my $i (1..$#{$PDB[0]}+1){
	if ((grep {$_ == $i} @low_layer) && (${$PDB[0]}[$i-1] eq $solvent_name)){
		 printf OUT "%s    0    %12.6f  %12.6f  %12.6f  L\n", ${$PDB[4]}[$i-1], ${$PDB[1]}[$i-1], ${$PDB[2]}[$i-1], ${$PDB[3]}[$i-1];
	}
}
close(OUT);
#open (INP, "<real_layers.xyz") or die;
#$count = 0;
#while (<INP>){
#	$_ =~ s/^\s*//;
#	$_ =~ s/\s*$//;
#	@_ = split(/\s+/,$_);
#	if ($_[0] ne ${$PDB[4]}[$count]){
#		print color("red"), "\nSomething went wrong during sorting! Please compare atom sequence in real_layer.xyz and the .pdb file!\n", color("reset");
#		exit;
#	}
#	$count++;
#}
print "done!\n\n";
close(INP);
for my $i (0..(scalar(@hmlinks)+scalar(@hllinks)-1)){
	printf XYZ "%s    %12.6f  %12.6f  %12.6f\n", 'H', $H_pos[3*$i], $H_pos[3*$i+1], $H_pos[3*$i+2];
        printf PDB1 "ATOM  %5d %3s  LIG     1    %8.3f%8.3f%8.3f\n", $high_count+$i+1, 'H', $H_pos[3*$i], $H_pos[3*$i+1], $H_pos[3*$i+2];
}
close(XYZ);
close(XYZ1);
`cp model-H.xyz HighLayer.xyz`;
close(PDB1);
}

#####################################

sub prepare_nonstd{
my (@parm,@nonstd_types,@nonstd_param,@nonstd_crg,@nonstd_crgmethod,@nonstd_additional);
if ($nonstd_resid eq '0'){
       print "Only standard and/or pre-defined residues used throughout.\n";
} else {
        @nonstd_resid = split(/\s+/,$nonstd_resid);
	my $coord = read_input("coordinate file");
        my $nonstd_types = read_input("non-standard residues atom types");
        if ($nonstd_types eq '0'){
            for my $i (0..$#nonstd_resid){
                $nonstd_types[$i] = 'gaff';
            }
        } else {
		@nonstd_types = split(/\s+/,$nonstd_types);
        }
        my $nonstd_param = read_input("non-standard residues parameters");
        if ($nonstd_param eq '0'){
            for my $i (0..$#nonstd_resid){
		$nonstd_param[$i] = 'gaff.dat';
            }
        } else {
		@nonstd_param = split(/\s+/,$nonstd_param);
        }
        my $nonstd_crg = read_input("non-standard residues charges");
        if ($nonstd_crg eq '0'){
            for my $i (0..$#nonstd_resid){
                $nonstd_crg[$i] = 0.0;
            }
        } else {
		@nonstd_crg = split(/\s+/,$nonstd_crg);
	}
        my $nonstd_crgmethod = read_input("non-standard residues charge methods");
        if ($nonstd_crgmethod eq '0'){
            for my $i (0..$#nonstd_resid){
                $nonstd_crgmethod[$i] = 'bcc';
            }
        } else {
		@nonstd_crgmethod = split(/\s+/,$nonstd_crgmethod);
	}
        my $nonstd_additional = read_input("non-standard residues additional antechamber keys"); 
        if ($nonstd_additional eq '0'){
		$nonstd_additional = '';
        }

        printf "%2d non-standard residues will be generated\n",scalar(@nonstd_resid);
        print "Name  Type  Charge  Model\n";
        for my $i (0..$#nonstd_resid){
                if ($nonstd_resid[$i] =~ /\w{4,}/){
                        print color("red"), "Non-standard residue names must have a 3-letter code!\n", color("reset");
                        exit;
                }
                print "  $nonstd_resid[$i]    $nonstd_types[$i]    $nonstd_crg[$i]    $nonstd_crgmethod[$i]\n";
        }
        for my $i (0..$#nonstd_resid){
		my ($col3,$col5,$col7,$col9,$col10,$col11);
		my $tmp_number = 0;
		my @resid = ();
                open PDBI, ">$nonstd_resid[$i].pdb" or die;
		open PDB, "<$coord" or die;
		while (<PDB>){
                        if (/^\s*ATOM\s+(\d+)/){
                                @_ = split(//,$_);
                                $col3 = join('',@_[12..15]);    #atom name
                                $col3 =~ s/^\s*//;
                                $col3 =~ s/\s*$//;
                                $col3 =~ s/\d+//g;
                                $col5 = join('',@_[17..19]);    #residue name
                                $col5 =~ s/^\s*//;
                                $col5 =~ s/\s*$//;
                                $col7 = join('',@_[22..25]);    #residue number
                                $col7 =~ s/^\s*//;
                                $col7 =~ s/\s*$//;
                                $col9 = join('',@_[30..37]);    #x-coord
                                $col9 =~ s/^\s*//;
                                $col9 =~ s/\s*$//;
                                $col10 = join('',@_[38..45]);   #y-coord
                                $col10 =~ s/^\s*//;
                                $col10 =~ s/\s*$//;
                                $col11 = join('',@_[46..53]);   #z-coord
                                $col11 =~ s/^\s*//;
                                $col11 =~ s/\s*$//;
				if ($col5 eq $nonstd_resid[$i] && $col7 == $tmp_number){
					print PDBI $_;
                                } elsif ($col5 eq $nonstd_resid[$i] && $tmp_number == 0){
                                        $tmp_number = $col7;
					print PDBI $_;
                                } elsif (($col5 ne $nonstd_resid[$i] && $tmp_number != 0) || ($col5 eq $nonstd_resid[$i] && $col7 != $tmp_number)){
                                        last;
                                }
                        }
                }
		close(PDB);
		close(PDBI);
        }
	for my $i (0..$#nonstd_resid){
		print "Generating non-standard residue $i ... \n";
                print "$AMBERHOME/bin/antechamber -fi pdb -i $nonstd_resid[$i].pdb -fo mol2 -o $nonstd_resid[$i]\_out.mol2 -at $nonstd_types[$i] -nc $nonstd_crg[$i] -c $nonstd_crgmethod[$i] $nonstd_additional\n";
                `$AMBERHOME/bin/antechamber -fi pdb -i $nonstd_resid[$i].pdb -fo mol2 -o $nonstd_resid[$i]\_out.mol2 -at $nonstd_types[$i] -nc $nonstd_crg[$i] -c $nonstd_crgmethod[$i] $nonstd_additional\n`;
                print "$AMBERHOME/bin/parmchk2 -i $nonstd_resid[$i]\_out.mol2 -f mol2 -o $nonstd_resid[$i].parm -p $AMBERHOME/dat/leap/parm/$nonstd_param[$i]\n";
		`$AMBERHOME/bin/parmchk2 -i $nonstd_resid[$i]\_out.mol2 -f mol2 -o $nonstd_resid[$i].parm -p $AMBERHOME/dat/leap/parm/$nonstd_param[$i]`;
                open (PARM, "<".$nonstd_resid[$i].".parm") or die;
                while (<PARM>){
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
                        if (scalar(@_) > 3 && $_[-3] eq 'ATTN,' && $_[-2] eq 'need' && $_[-1] eq 'revision'){
                                @_ = split(/-/,$_[0]);
                                if (!grep {$_ eq 'AL'} @_){
                                        print color("red"), "One or more missing parameters detected for residue $nonstd_resid[$i]. Complete $nonstd_resid[$i].parm and re-run leap_$nonstd_resid[$i].inp, leap_model-H.inp and leap_real.inp!\n", color("reset");
                                }
                        }
                }
                close(PARM);
		open (IN, ">leap_".$nonstd_resid[$i].".inp") or die;
		print IN "logFile leap.log\n";
		for my $j (@ff){
			print IN "source $j\n";
		}
		print IN "loadamberparams $nonstd_resid[$i].parm\n$nonstd_resid[$i] = loadmol2 $nonstd_resid[$i]\_out.mol2\n";
		print IN "check $nonstd_resid[$i]\nsaveoff $nonstd_resid[$i] $nonstd_resid[$i].lib\nquit\n";
		close(IN);
		`if [ -e leap.log ]; then rm leap.log; fi; tleap -f leap_$nonstd_resid[$i].inp; cp leap.log leap_$nonstd_resid[$i].log`;
		open LOG, "<leap.log" or die;
		while (<LOG>){
		        if (/^\s*saveOff: Argument \#1 is type String must be of type: \[unit parameter_set list\]/){
		                print color("red"), "Could not generate non-standard library for $nonstd_resid[$i]. Check leap_$nonstd_resid[$i]_error.log\n", color("reset");
				`mv leap.log leap_$nonstd_resid[$i]_error.log`;
				exit;
		        } elsif (/^\s*Building topology/){
		                print "done!\n";
		                last;
		        }
		}
		close(LOG);
	}
}

if ($predef_resid ne '0'){
	@predef_resid = split(/\s+/,$predef_resid);
	for my $i (0..($#predef_resid-2)/3){
		print "\nChecking predefined residue $predef_resid[3*$i] ...\n";
		open (IN, ">leap_".$predef_resid[3*$i].".inp") or die;
		print IN "logFile leap.log\n";
	        for my $j (@ff){
			print IN "source $j\n";
 	        }
	        print IN "loadamberparams $predef_resid[3*$i+2]\n";
		$_ = (split(/\./,$predef_resid[3*$i+1]))[1];
		if ($_ eq 'in' ){
			print IN "loadamberprep $predef_resid[3*$i+1]\n";
		} elsif ($_ eq 'mol2'){
			print IN "$predef_resid[3*$i] = loadmol2 $predef_resid[3*$i+1]\n";
		} elsif ($_ eq 'pdb'){
                        print IN "$predef_resid[3*$i] = loadpdb $predef_resid[3*$i+1]\n";
                } elsif ($_ eq 'lib'){
                        print IN "loadOff $predef_resid[3*$i+1]\n";
                }
		if ($_ ne 'lib'){
                        print IN "check $predef_resid[3*$i]\nsaveoff $predef_resid[3*$i] $predef_resid[3*$i].lib\nsavemol2 $predef_resid[3*$i] $predef_resid[3*$i]\_out.mol2 1\nquit\n";
                } else {
                        print IN "check $predef_resid[3*$i]\nsavemol2 $predef_resid[3*$i] $predef_resid[3*$i]\_out.mol2 1\nquit\n";
                }
	        close(IN);
	        `if [ -e leap.log ]; then rm leap.log; fi; tleap -f leap_$predef_resid[3*$i].inp; cp leap.log leap_$predef_resid[3*$i].log`;
		open LOG, "<leap.log" or die;
		while (<LOG>){
                	if (/^\s*saveOff: Argument \#1 is type String must be of type: \[unit parameter_set list\]/){
                	        print color("red"), "Could not load pre-defined non-standard library for $nonstd_resid[$i]. Possibly missing parameters. Check leap_$predef_resid[$i]_error.log\n", color("reset");
                	        `mv leap.log leap_$predef_resid[$i]_error.log`;
				exit;
                	} elsif (/^\s*check: Argument \#1 is type String must be of type: \[unit molecule residue atom\]/){
				print color("red"), "Could not load pre-defined non-standard library for $nonstd_resid[$i]. Possibly missing parameters. Check leap_$predef_resid[$i]_error.log\n" , color("reset");
				exit;
			} elsif (/^\s*Building topology/){
                	        print "done!\n";
                	        last;
                	}       
        	}
        	close(LOG);
	}
}
}

#####################################

sub load_std{
my $residue = $_[0];
my $resid_length = $_[1];
my $resid_number = $_[2];
my $residue_old = $residue;
if (grep {$_ eq $residue} (keys %aminoacids)){
	if ($resid_length == $aminoacids{$residue}){
	} elsif ($resid_length eq $N_aminoacids{$residue}){
		print "Residue $residue is an N-terminal residue. N$residue will be loaded instead\n";
		$residue_old = $residue;
		$residue = "N".$residue;
		$term_resid_counter++;
	} elsif ($resid_length eq $C_aminoacids{$residue}){
		print "Residue $residue is an C-terminal residue. C$residue will be loaded instead\n";
		$residue_old = $residue;
	        $residue = "C".$residue;
		$term_resid_counter++;
	}
}
open (OUT, ">leap_".$residue.".inp") or die;
print OUT "logFile leap.log\n";
for my $j (@ff){
        print OUT "source $j\n";
}
print OUT "saveoff $residue $residue.lib\nsavemol2 $residue $residue\_out.mol2 1\nquit\n";
close(OUT);
`if [ -e leap.log ]; then rm leap.log; fi;  tleap -f leap_$residue.inp; cp leap.log leap_$residue.log`;
open LOG, "<leap.log" or die;
while (<LOG>){
	if (/^\s*saveOff: Argument \#1 is type String must be of type: \[unit parameter_set list\]/){
		print color("yellow"), "Could not load library for $residue. Possibly an unspecified non-standard library!\n", color ("reset");
		print color("yellow"), "This is not a problem if the $residue.lib is present, otherwise please provide a .lib or include in the non-standard residues or pre-defined non-standard residues\n", color("reset");
	} elsif (/^\s*Building topology/){
		print "Loading library for residue $residue.\n";
		last;
	}
}
close(LOG);
if (grep {$_ eq $residue_old} (keys %aminoacids)){
	if (($resid_length eq $N_aminoacids{$residue_old}) || ($resid_length eq $C_aminoacids{$residue_old})){
		my $new_resid_name = $alphabet[-$term_resid_counter].$alphabet[-$term_resid_counter].$alphabet[-$term_resid_counter];
		system("sed 's/$residue/$new_resid_name/g' $residue.lib > tmp.lib; mv tmp.lib $residue.lib");
		system("mv $residue.lib $new_resid_name.lib");
		my $coord = read_input("coordinate file");
		my $old_pdb_inp = sprintf("%3s  %4d", $residue_old, $resid_number);
		my $new_pdb_inp = sprintf("%3s  %4d", $new_resid_name, $resid_number); 
		system("sed 's/$old_pdb_inp/$new_pdb_inp/g' $coord > tmp.pdb ; mv tmp.pdb $coord");
		$residue = $new_resid_name;
	}
}
return $residue;
}

#####################################

sub check_pdb{
my $coord = read_input("coordinate file");
open (INP, "<".$coord) or die;
open my $out, ">tmp" or die;
my (@tmp,@tmp_atom_names,@tmp_resid,$old_resid_number,$old_resid_name,$pdb_altered,$reorder,$count);
my ($col2,$col2_old,$col3,$col5,$col7,@tmp_equivalent_names,@lib_names,$resid_length);
$old_resid_number = 1;
$old_resid_name = '';
$count = 0;
$col2_old = 1;
$pdb_altered = 'no';
my $keep = 0;
while (my $tmp = <INP>){
        if ($tmp =~ /^\s*ATOM\s+(\d+)/){
                @_ = split(//,$tmp);
		$col2 = join('',@_[6..10]);     #atom number
                $col2 =~ s/^\s*//;
                $col2 =~ s/\s*$//;
                $col3 = join('',@_[12..15]);    #atom name
                $col3 =~ s/^\s*//;
                $col3 =~ s/\s*$//;
                $col5 = join('',@_[17..19]);    #residue name
                $col5 =~ s/^\s*//;
                $col5 =~ s/\s*$//;
		$col7 = join('',@_[22..25]);    #residue number
                $col7 =~ s/^\s*//;
                $col7 =~ s/\s*$//;
		if ($old_resid_number != $col7){
			$count = $col2 - scalar(@tmp_resid);
			$resid_length = $col2 - $col2_old;
			$col2_old = $col2;
			$old_resid_number = $col7;
			if ( (!(grep {$_ eq $old_resid_name} (@nonstd_resid,@predef_resid,@std_resid))) || ((grep {$_ eq $old_resid_name} (keys %aminoacids)) && $resid_length != $aminoacids{$old_resid_name}) ){
				#if (!-e "$BASEDIR/$old_resid_name.lib"){
					$old_resid_name = load_std($old_resid_name,$resid_length,$col7-1);
				#}
				push(@std_resid,$old_resid_name);
				#if one adds $old_resid_name to &aminoacids and $old_resid_name is WAT (which can happen if WAT is part of user-defined Medium or High layer, tleap will be executed for each water molecule in the .pdb
				#if ( !(grep {$_ eq $old_resid_name} (keys %aminoacids)) ){ 
				#	$aminoacids{$old_resid_name} = $resid_length;
				#}
			}
			@tmp_equivalent_names = generate_equivalent_names(@tmp_atom_names);
			($reorder,@lib_names) = check_order($old_resid_name,@tmp_equivalent_names);
			if ($reorder eq 'yes'){
				reorder_pdb($out,\@tmp_equivalent_names,\@lib_names,\@tmp_resid,$count);
				$pdb_altered = 'yes';
                                if ($keep == 1){
					print $out "TER\n";
                                        $keep = 0;
				}
			} else {
				for my $i (0..$#tmp_resid){
					print $out $tmp_resid[$i];
					$crd_2_pdb{$count + $i} = $count + $i;
				} 
                                if ($keep == 1){
                                	print $out "TER\n";
                                	$keep = 0;
				}
			}
			@tmp_resid = ();
			push(@tmp_resid,$tmp);
			@tmp_atom_names = ();
			push(@tmp_atom_names,$col3);
			$old_resid_name = $col5;
		} else {
			$old_resid_name = $col5;
			push(@tmp_resid,$tmp);
			push(@tmp_atom_names,$col3);
		}
	} elsif ($tmp =~ /^\s*TER/){
		#print $out $tmp;
                $keep = 1;
	}
}
$count = $col2 - scalar(@tmp_resid);
$resid_length = $col2 - $col2_old + 1;
$old_resid_number = $col7;
if ( (!(grep {$_ eq $old_resid_name} (@nonstd_resid,@predef_resid,@std_resid))) || ((grep {$_ eq $old_resid_name} (keys %aminoacids)) && $resid_length != $aminoacids{$old_resid_name}) ){
	#if (!-e "$BASEDIR/$old_resid_name.lib"){
        	$old_resid_name = load_std($old_resid_name,$resid_length,$col7);
	#}
        push(@std_resid,$old_resid_name);
}
@tmp_equivalent_names = generate_equivalent_names(@tmp_atom_names);
($reorder,@lib_names) = check_order($old_resid_name,@tmp_equivalent_names);
if ($reorder eq 'yes'){
        reorder_pdb($out,\@tmp_equivalent_names,\@lib_names,\@tmp_resid,$count);
        $pdb_altered = 'yes';
} else {
	for my $i (0..$#tmp_resid){
	        print $out $tmp_resid[$i];
                $crd_2_pdb{$count+$i+1} = $count+$i+1;
        }
}
print $out "TER\nEND\n";
close(INP);
close($out);
if ($pdb_altered eq 'yes'){
	system("mv $coord $coord\_orig; mv tmp $coord");
	print color("red"), "The atom order for one or more residues has to be changed in the .pdb file! Please inspect the new .pdb file and modify the definition of the High and Medium layers, as well as of the links in the cobram.parm accordingly", color("reset");
	exit;
} else {
	system("rm tmp");
}
%pdb_2_crd = reverse %crd_2_pdb;
}

#####################################

sub generate_equivalent_names{
my ($iname,@tmp);
my @equivalent_names = ();
my @atom_names = @_;
for my $i (0..$#atom_names){
	push(@{$equivalent_names[$i]},$atom_names[$i]);
	@tmp = split(//,$atom_names[$i]);
        if (($tmp[0] =~ /\d/ || $tmp[0] =~ /\'/) && scalar(@tmp) == 3){
                $iname = join('',$tmp[1],$tmp[2],$tmp[0]);
                push(@{$equivalent_names[$i]},$iname);
        } elsif (($tmp[0] =~ /\d/ || $tmp[0] =~ /\'/) && scalar(@tmp) == 4){
                $iname = join('',$tmp[1],$tmp[2],$tmp[3],$tmp[0]);
                push(@{$equivalent_names[$i]},$iname);
	}
        if (scalar(@tmp) == 3 && ($tmp[2] =~ /\d/ || $tmp[2] =~ /\'/) && $tmp[0] !~ /\d/){
                $iname = join('',$tmp[2],$tmp[0],$tmp[1]);
                push(@{$equivalent_names[$i]},$iname);
        } elsif (scalar(@tmp) == 4 && ($tmp[3] =~ /\d/ || $tmp[3] =~ /\'/) && $tmp[0] !~ /\d/){
                $iname = join('',$tmp[3],$tmp[0],$tmp[1],$tmp[2]);
                push(@{$equivalent_names[$i]},$iname);
	}
	my @tmp_equivalent_names = @{$equivalent_names[$i]};
	for my $k (@tmp_equivalent_names){
		for my $j (0..$#equivalence){
			if (grep {$_ eq $k} @{$equivalence[$j]}){
				for my $l (@{$equivalence[$j]}){
				        @tmp = split(//,$l);
                                        if (($tmp[0] =~ /\d/ || $tmp[0] =~ /\'/) && scalar(@tmp) == 3){
                                                $iname = join('',$tmp[1],$tmp[2],$tmp[0]);
                                                push(@{$equivalent_names[$i]},$iname);
                                        } elsif (($tmp[0] =~ /\d/ || $tmp[0] =~ /\'/) && scalar(@tmp) == 4){
                                                $iname = join('',$tmp[1],$tmp[2],$tmp[3],$tmp[0]);
                                                push(@{$equivalent_names[$i]},$iname);
                                        }
                                        if (scalar(@tmp) == 3 && ($tmp[2] =~ /\d/ || $tmp[2] =~ /\'/) && $tmp[0] !~ /\d/){
                                                $iname = join('',$tmp[2],$tmp[0],$tmp[1]);
                                                push(@{$equivalent_names[$i]},$iname);
                                        } elsif (scalar(@tmp) == 4 && ($tmp[3] =~ /\d/ || $tmp[3] =~ /\'/) && $tmp[0] !~ /\d/){
                                                $iname = join('',$tmp[3],$tmp[0],$tmp[1],$tmp[2]);
                                                push(@{$equivalent_names[$i]},$iname);
                                        }
				}
				push(@{$equivalent_names[$i]},@{$equivalence[$j]});
        	        }
		}
	}
	my %dummy;
	@{$equivalent_names[$i]} = grep { ! $dummy{ $_ }++ } @{$equivalent_names[$i]};
	if (scalar(@{$equivalent_names[$i]} > 1) && $print eq 'yes'){
		#print "Alternative names for $atom_names[$i]: @{$equivalent_names[$i]}\n";
	}
}
return @equivalent_names;
}

#####################################

sub check_order{
my $resid = shift(@_);
my @equivalent_names = @_;
my @lib_names = ();
my $reorder = 'no';
#print "Resid $resid\n";
open LIB, "<".$resid.".lib" or die;
while (<LIB>){
        if (/!entry\.$resid\.unit\.atoms table/){
                my $count = 0;
                while (<LIB>){
                        if (/!entry\.$resid\.unit\.atomspertinfo/){
                                last;
                        }
			$_ =~ s/^\s*//;
			$_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
			$_[0] =~ s/\"//g;
			push(@lib_names,$_[0]);
                }
                last;
        }
}
#print "Lib $resid: @lib_names\n";
for my $i (0..$#lib_names){
	if (!(grep {$_ eq $lib_names[$i]} @{$equivalent_names[$i]})){
		printf "Atom name %4s of residue %3s could not be found among atom names @{$equivalent_names[$i]}. Order in the .pdb and .lib differs. Will try to reorder the atoms in the .pdb to match the order of the .lib! The original pdb file will be stored as.pdb_orig.\n", $lib_names[$i], $resid;
		$reorder = 'yes';
		last;
	}
}
return($reorder,@lib_names);
}

#####################################

sub reorder_pdb{
my $out = $_[0];
my @equivalent_names = @{$_[1]};
my @lib_names = @{$_[2]};
my @resid = @{$_[3]};
my $count = $_[4]-1;
my $crd_count = $count;
for my $i (0..$#lib_names){
	for my $j (0..$#equivalent_names){
		if (grep {$_ eq $lib_names[$i]} @{$equivalent_names[$j]}){
			$count++;
			substr($resid[$j], 6, 5) = sprintf("%5d", $count); # replace the atom number
			print $out $resid[$j];
			$crd_2_pdb{$count} = $crd_count+$j+1; # map new order ($count, key) to original .crd order ($crd_count + $j + 1, value)
			last;
		}
	}
}
#print ".pdb order => .crd order\n";
#for my $key ( keys %crd_2_pdb ) {
#	my $value = $crd_2_pdb{$key};
#        print "$key => $value\n";
#}
}

#####################################

sub prepare_model{
my ($correct_type,$correct_charge);
print "Generating model-H.top ...\n";
#`babel -ixyz model-H.xyz -omol2 model-H.mol2`;
#`$AMBERHOME/bin/antechamber -fi pdb -i model-H1.pdb -fo mol2 -o model-H.mol2 -at amber -nc 0 -c bcc`;
`$AMBERHOME/bin/antechamber -fi pdb -i model-H1.pdb -fo mol2 -o model-H.mol2 -dr no`;
#`rm model-H1.pdb`;
open XYZ, "<model-H.xyz" or die;
$_ = <XYZ>;
$_ = <XYZ>;
open MOL1, "<model-H.mol2" or die;
open MOL2, ">tmp.mol2" or die;
while (<MOL1>){
    $_ =~ s/^\s*//;
    $_ =~ s/\s*$//;
    @_ = split(/\s+/,$_);
    if ($_[0] && $_[0] =~/\d+/ && $_[2] && $_[2] =~/\d+\.\d+/ && $_[3] && $_[3] =~/\d+\.\d+/ && $_[4] && $_[4] =~/\d+\.\d+/ && $_[5] && $_[5] =~/\w+/){
            my $tmp = <XYZ>;
            $tmp =~ s/^\s*//;
            $tmp =~ s/\s*$//;
            my @tmp = split(/\s+/,$tmp);
            $_ = join('  ',$_[0],$tmp[0],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7],$_[8]);
    }
    print MOL2 "$_\n";
}
close(XYZ);
close(MOL1);
close(MOL2);
`mv tmp.mol2 model-H.mol2 `;
open (INP, "<model-H.mol2") or die;
open (OUT, ">tmp") or die;
my $curr_resid_num = 0;
my $prev_resid_num = 0;
my $number = 0;
my $generic_name;
while(<INP>){
	$_ =~ s/^\s*//;
        @_ = split(/\s+/,$_);
	if ($_[0] && $_[0] =~/\d+/ && $_[0] <= scalar(@high_layer) && $_[5] && $_[5] =~ /\w+/){
		my $model_num = $_[0];
		#print "Number of atom $_[1] in the model-H: $model_num\n";
		my $real_num = $high_layer[$_[0]-1];
		#print "Number of atom $_[1] in real_layers: $real_num\n";
		my $curr_resid = ${$PDB[0]}[$real_num-1];
		$curr_resid_num = ${$PDB[5]}[$real_num-1];
		if ($curr_resid_num ne $prev_resid_num){
			$number++;
			$prev_resid_num = $curr_resid_num;
			$generic_name = $_[7];
		} 
		push(@model_resid,$curr_resid);
		my $res_num = 0;
		for my $i (0..$real_num-1){
			if (${$PDB[0]}[$i] eq $curr_resid){  
				$res_num++;
			}
		}
		#print "Number of atom $_[1] in residue $curr_resid: $res_num\n";
		open (IN, "<".$curr_resid."_out.mol2") or die;
		while (my $tmp = <IN>){
			$tmp =~ s/^\s*//;
		        my @tmp = split(/\s+/,$tmp);
			if ($tmp[0] && $tmp[0] !~ /^\d+$/ ){
				next;
			}
			if ($tmp[0] && $tmp[0] =~ /^\d+$/ && $tmp[1] =~ /^\d+$/ && $tmp[2] =~ /^\d+$/ && $tmp[3] =~ /^\d+$/ && $tmp[4] =~ /^\d+$/){
				#print "$tmp\n";
				$res_num = recursive_substraction($res_num,\@tmp);
			#print "Recursive Res_num $res_num\n";
			}
			if ($tmp[0] && $tmp[0] == $res_num && $tmp[5] && $tmp[5] =~ /^\w+/){
				$correct_type = $tmp[5];
				#charge taken from .lib file
				###$correct_charge = $tmp[8];
				#print "Correct_type $correct_type\n Correct_charge $correct_charge\n\n";
				last;
			}
		}
		close(IN);
		$_ = join('  ',$_[0],$_[1],$_[2],$_[3],$_[4],$correct_type,$number,$generic_name,'0.00',"\n");
	}
	if ($_[0] && $_[0] =~/\d+/ && $_[0] > scalar(@high_layer) && (split(//,$_[1]))[0] eq 'H'){
		$_[5] = 'AL';
		$_ = join('  ',$_[0],$_[1],$_[2],$_[3],$_[4],$_[5],$_[6],$_[7],'0.00',"\n");
	} 
	print OUT $_;
}
close(INP);
close(OUT);
`mv tmp model-H_out.mol2`;
my %dummy = ();
@_ = grep { ! $dummy{ $_ }++ } @model_resid;
@model_resid = @_;
open (IN, ">leap_model.inp") or die;
print IN "logFile leap.log\n";
for my $j (@ff){
	print IN "source $j\n";
}
for my $i (@model_resid){
	print IN "loadoff $i.lib\n";
	if (grep {$_ eq $i} @predef_resid){
		for my $j (0..($#predef_resid-2)/3){
			if ($i eq $predef_resid[3*$j]){
				print IN "loadamberparams $predef_resid[3*$j+2]\n";
			}
		}
	} elsif (grep {$_ eq $i} @std_resid){
	} else {
		$_ = `echo -n $i.parm`;
		#if ($_ eq $i.".parm"){
                if (! -e "$BASEDIR/$i.parm"){
		      	print color("red"), "There is no .parm file for residue $i\n", color("reset");
		        exit;	
		} elsif (scalar(split(/\s+/,$_)) > 1){
			print color("red"), "There is more than one .parm file for residue $i\n", color("reset");
			exit;
		} else {
			print IN "loadamberparams $_\n";
		}
	}
}
print IN "loadamberparams frcmod.ionsjc_tip3p\n";
print IN "model = loadmol2 model-H_out.mol2\n";
if (read_input("tleap instructions for the high layer")){
       $_ = read_input("tleap instructions for the high layer");
       my @instructions = split(/;/,$_);
       for my $i (@instructions){
               print IN $i."\n";
       }
}
print IN "check model\nsaveamberparm model model-H.top model-H.crd\nsavepdb model model-H.pdb\nquit\n";
close(IN);
`if [ -e leap.log ]; then rm leap.log; fi; tleap -f leap_model.inp; cp leap.log leap_model.log`;
my $link_atoms = read_input("links");
if ($link_atoms eq '0'){
        return;
} else {
	open (OUT, ">model-H.parm") or die;
	print OUT "remark goes here\nMASS\n";
	@_ = ();
	print OUT "AL 0.000         0.000               ATTN, need revision\n\n";
	open (IN, "<leap.log") or die;
	print OUT "BOND\n";
	while (<IN>){
		if (/^\s*Could not find bond parameter for:\s+(\w.*)/){
			my @tmp = split(/\s+/,$1);
			if (grep {$_ eq 'AL'} @tmp){
			} else {
				print color("red"), "No bond parameter was found for $1 during generation of model-H.top. Correct this and execute tleap -f leap_model.inp again.\n", color("reset");
			}
			if (grep {$_ eq $1} @_){
			} else {
				print OUT "$1    0.00   0.000       ATTN, need revision\n";
				push(@_,$1);
			}
		}
	}
	close(IN);
	open (IN, "<leap.log") or die;
	print OUT "\nANGLE\n";
	while (<IN>){
	        if (/^\s*Could not find angle parameter:\s+(\w.*)/){
			my @tmp = split(/\s+/,$1);
                        if (grep {$_ eq 'AL'} @tmp){
                        } else {
                                print color("red"), "No angle parameter was found for $1 during generation of model-H.top. Correct this and execute tleap -f leap_model.inp again.\n", color("reset");
                        }
			if (grep {$_ eq $1} @_){
			} else {
	                	print OUT "$1   0.00   0.000       ATTN, need revision\n";
				push(@_,$1);
			}
	        }
	}
	close(IN);
	open (IN, "<leap.log") or die;
	print OUT "\nDIHE\n\n";
	while (<IN>){
                if (/^\s*Could not find dihedral parameter:\s+(\w.*)/){
                        my @tmp = split(/\s+/,$1);
                        if (grep {$_ eq 'AL'} @tmp){
                        } else {
                                print color("red"), "No dihedral parameter was found for $1 during generation of model-H.top. Correct this and execute tleap -f leap_model.inp again.\n", color("reset");
                        }
		}
	}
	close(IN);
	print OUT "IMPROPER\n\n";
	print OUT "NONBON\n  AL          0.0000  0.0000             ATTN, need revision\n"; 
	close(OUT);
	open (IN, "<leap_model.inp") or die;
	open (LEAP, ">tmp") or die;
	while (<IN>){
		if (/check model/){
			print LEAP "loadamberparams model-H.parm\n";
			print LEAP "check model\n";
		} else {
			print LEAP "$_";
		}
	}
	close(LEAP);
	`if [ -e leap.log ]; then rm leap.log; fi; mv tmp leap_model.inp; tleap -f leap_model.inp; cp leap.log leap_model.log`;
	open (IN, "<leap.log") or die;
	@_=();
	while (<IN>){
                if (/^.*No torsion terms for\s+(\w.*)/){
                        my @tmp = split(/-/,$1);
                        if (grep {$_ eq 'AL'} @tmp){
                        } else {
                                print color("red"), "No torsion parameter was found for $1 during generation of model-H.top. Correct this and execute tleap -f leap_model.inp again.\n", color("reset");
                        }
                        if (grep {$_ eq $1} @_){
                        } else {
                                `sed "s/DIHE/DIHE\\n$1 0.0 0.0 0.0 ATTN, need revision/" model-H.parm > tmp1; mv tmp1 model-H.parm`;
                                push(@_,$1);
                        }
                }
        }
	close(IN);
	`tleap -f leap_model.inp; cp leap.log leap_model.log`;
}
open (IN, "<leap.log") or die;
while (<IN>){
	if (/^\s*Added missing heavy atom/){
        	print color("red"), "AMBER $_.\n", color("reset");
                exit;
	}
}
close(IN);
print "... done!\n";
}

######################################

sub recursive_substraction{

my $res_num = $_[0];	
my $tmp = $_[1];
my @tmp = @$tmp;
if ($res_num > $tmp[0]){
	$res_num = $res_num-$tmp[0];
	recursive_substraction($res_num,\@tmp);
} elsif ($res_num <= $tmp[0]){
	return $res_num;
}
}

######################################

sub check_residue{
my $iatom = $_[0];
foreach my $i (@{$RES[${$PDB[5]}[$iatom-1]-1]}){
	if (!grep {$_ eq $i} @high_layer){
		return 'HM';
	}
}
return 'H';
}

######################################

sub prepare_real{
my ($col1,$col2,$col3,$col4,$col5,$col6,$col7,$col8,$col9,$col10,$col11);
my $count = 0;
for my $i (@hmlinks,@hllinks){
        my $check=check_residue($i);
        if ($check eq 'H'){
        	next;
        } else {
		if ((grep {$_ eq ${$PDB[0]}[$i-1]} @UNIQUE_RES) && (grep {$_ eq ${$PDB[5]}[$i-1]} @UNIQUE_RES_NUM)){
			push (@REASSIGNED_RES,@{$RES[${$PDB[5]}[$i-1]-1]});
	        } else {
			push(@UNIQUE_RES,${$PDB[0]}[$i-1]);
			push(@UNIQUE_RES_NUM,${$PDB[5]}[$i-1]);
	                $new_resid->{${$PDB[5]}[$i-1]} = $alphabet[$count].$alphabet[$count].$alphabet[$count];
	                $count++;
			push (@REASSIGNED_RES,@{$RES[${$PDB[5]}[$i-1]-1]});  # $i is pdb atom number (starts at 1), $PDB[5]}[$i-1] is pdb resid number (starts at 1)
	        }
	}
}
#for my $i (keys %$new_resid){
#        print "Residue: $i, generic residue name: $new_resid->{$i}\n";
#}
#print "Unique residues: @UNIQUE_RES\n";
#print "Unique residue numbers: @UNIQUE_RES_NUM\n";
my %dummy = ();
my @REASSIGNED_RES = grep { ! $dummy{ $_ }++ } @REASSIGNED_RES;
#print "ATOMS IN RESIDUES TO BE REASSIGNED: @REASSIGNED_RES\n";
#print "UNIQUE RESIDUES FOR CHARGE REASSIGNMENT: @UNIQUE_RES\n";
foreach my $i (@UNIQUE_RES_NUM){
        reassign_charges($i);
}
my $coord = read_input("coordinate file");
print "\nBuilding real.pdb file ... ";
open (INP, "<".$coord) or die;
open (OUT, ">real.pdb") or die;
my @tmp;
my $kilo = 0;
my $TER = 0;
while (my $tmp = <INP>){
	if ($tmp =~ /^\s*ATOM\s+(\d+)/){
	       if ((grep {$_ eq $1} @REASSIGNED_RES)){
			@_ = split(//,$tmp);
		        $col1 = join('',@_[0..5]);
	                $col2 = join('',@_[6..10]);
	                $col2 =~ s/^\s*//;
        	        $col2 =~ s/\s*$//;
        	        $col3 = join('',@_[12..15]);
        	        $col4 = $_[16];
        	        $col5 = join('',@_[17..19]);
			$col5 =~ s/^\s*//;
                        $col5 =~ s/\s*$//;
        	        $col6 = $_[21];
        	        $col7 = join('',@_[22..25]);
                	$col8 = $_[26];
                	$col9 = join('',@_[30..37]);
                	$col10 = join('',@_[38..45]);
                	$col11 = join('',@_[46..53]);
			if ($col6 eq '@'){
	                        $kilo = 10000;
	                }
                        my $col3_tmp = $col3;
                        $col3_tmp =~ s/^\s*//;
                        $col3_tmp =~ s/\s*$//;
                        my @tmp = split(//,$col3_tmp);
        		if (($tmp[0] =~ /\d/ || $tmp[0] =~ /\'/) && scalar(@tmp) == 4){
                		$col3 = join('',$tmp[1],$tmp[2],$tmp[3],$tmp[0]);
        		}
			printf OUT "%6s%5d %4s%s%3s %s%4d%s   %8.3f%8.3f%8.3f\n",$col1,$col2,$col3,$col4,$new_resid->{$col7+$kilo},$col6,$col7,$col8,$col9,$col10,$col11;
                        $TER = 1;
	       } elsif ((grep {$_ == $1} @high_layer,@medium_layer) || ((grep {$_ == $1} @low_layer) && (${$PDB[0]}[$1-1] ne $solvent_name))) {
                        @_ = split(//,$tmp);
                        $col3 = join('',@_[12..15]);
                        $col3 =~ s/^\s*//;
                        $col3 =~ s/\s*$//;
                        my @tmp = split(//,$col3);
                        if (($tmp[0] =~ /\d/ || $tmp[0] =~ /\'/) && scalar(@tmp) == 4){
                                @_[12..15] = ($tmp[1],$tmp[2],$tmp[3],$tmp[0]);
                        }
                        $tmp = join('',@_);
			print OUT "$tmp";
                        $TER = 1;
	       } else {
                        $TER = 0;
			next;
	       }
	} elsif ($TER == 1 && $tmp =~ /^\s*TER/){
        	print OUT "$tmp";
        }
}
close(INP);
print "done!\n";
$TER = 0;
open (INP, "<".$coord) or die;
while (my $tmp = <INP>){
	if ($tmp =~ /^\s*ATOM\s+(\d+)/){
		if ((grep {$_ == $1} @low_layer) && (${$PDB[0]}[$1-1] eq $solvent_name)){
			print OUT "$tmp";
                        $TER = 1;
		} else {
			$TER = 0;
		}
	} elsif ($TER == 1 && $tmp =~ /^\s*TER/){
                print OUT "$tmp";
        }
}
close(INP);
close(OUT);
#pdb_to_mol2();
open (IN, ">leap_real.inp") or die;
print IN "logFile leap.log\n";
for my $j (@ff){
        print IN "source $j\n";
}
foreach my $i (@nonstd_resid){
	print IN "loadoff $i.lib\n";
	$_ = `echo -n $i.parm`;
	#if ($_ eq $i.".parm"){
        if (! -e "$BASEDIR/$i.parm"){
	        print color("red"), "There is no .parm file for residue $i\n", color("reset");
		exit;
	} elsif (scalar(split(/\s+/,$_)) > 1){
		print color("red"), "There is more than one .parm file for residue $i\n", color("reset");
		exit;
	} else {
		print IN "loadamberparams $_\n";
	}
}
foreach my $i (@std_resid){
	print IN "loadoff $i.lib\n";
}
foreach my $i (0..($#predef_resid-2)/3){        
	print IN "loadoff $predef_resid[3*$i].lib\n";
	print IN "loadamberparams $predef_resid[3*$i+2]\n";
}
print IN "loadamberparams frcmod.ionsjc_tip3p\n";
foreach my $i (@UNIQUE_RES_NUM){
        print IN "loadoff $new_resid->{$i}.lib\n";
}
#print IN "real = loadmol2 real_out.mol2\n"; #tleap does not take residue types from .mol2 
print IN "real = loadpdb real.pdb\n";  #old way: use pdb -> rearranges the atoms according to the libraries
#if (read_input("additional bond definitions")){
#	$_ = read_input("additional bond definitions");
#	my @addbonds = split(/\s+/,$_);
#	for my $i (@addbonds){
#		@_ = split(/-/,$i);
#		print IN "bond real.$_[0] real.$_[1]\n";
#	}
#}
if (read_input("tleap instructions for the whole system")){
       $_ = read_input("tleap instructions for the whole system");
       my @instructions = split(/;/,$_);
       for my $i (@instructions){
               print IN $i."\n";
       }
}
print IN "saveamberparm real real.top real.cord\n";
print IN "savepdb real real_libs.pdb\nquit\n";
close(IN);
`if [ -e leap.log ]; then rm leap.log; fi; tleap -f leap_real.inp; cp leap.log leap_real.log`;
open (IN, "<leap.log") or die;
while (<IN>){
        if (/^\s*Unknown residue:\s+(\w+)/){
        	print color("red"), "Unknown residue $1 detected while preparing the real.top. Correct this and execute tleap -f leap_real.inp or repeat the whole procedure.\n", color("reset");
        }
}
close(IN);
open (IN, "<leap.log") or die;
while (<IN>){
	if (/^\s*Could not find bond parameter for:\s+(\w.*)/){
	        my @tmp = split(/\s+/,$1);
	        if (grep {$_ eq 'AL'} @tmp){
	        } else {
	                print color("red"), "No bond parameter was found for $1 during generation of real.top. Correct this and execute tleap -f leap_real.inp or repeat the whole procedure.\n", color("reset");
	        }
	}
}
close(IN);
open (IN, "<leap.log") or die;
while (<IN>){
	if (/^\s*Could not find angle parameter:\s+(\w.*)/){
	        my @tmp = split(/\s+/,$1);
	        if (grep {$_ eq 'AL'} @tmp){
	        } else {
	                print color("red"), "No angle parameter was found for $1 during generation of real.top. Correct this and execute tleap -f leap_real.inp or repeat the whole procedure.\n", color("reset");
	        }
	}
}
close(IN);
open (IN, "<leap.log") or die;
while (<IN>){
	if (/^\s*Could not find dihedral parameter:\s+(\w.*)/){
	        my @tmp = split(/\s+/,$1);
	        if (grep {$_ eq 'AL'} @tmp){
	        } else {
	                print color("red"), "No dihedral parameter was found for $1 during generation of real.top. Correct this and execute tleap -f leap_real.inp or repeat the whole procedure.\n", color("reset");
	        }
	}
}
close(IN);
open (IN, "<leap.log") or die;
while (<IN>){
        if (/^\s*Added missing heavy atom/){
                print color("red"), "AMBER $_.\n", color("reset");
                exit;
	}
}
close(IN);
}

######################################

sub reassign_charges{

my %dummy=();
my $resid_num_i = $_[0];
my ($resid,@full_resid_high_layer);
for my $i (@high_layer){
        if (${$PDB[5]}[$i-1] == $resid_num_i){
                $resid = ${$PDB[0]}[$i-1];
                last;
        }
}
print "\nReassigning charges for residue: $resid -> $new_resid->{$resid_num_i} ...\n";
my @resid_high_layer = ();
my $resid_charge;
open (IN, ">leap_crg.inp") or die;
print IN "logFile leap.log\n";
for my $j (@ff){
        print IN "source $j\n";
}
print IN "loadoff $resid.lib\n";
print IN "charge $resid\nquit\n";
close(IN);
`if [ -f leap.log ]; then rm leap.log; fi; tleap -f leap_crg.inp`;
open (LEAP, "<leap.log") or die;
while (<LEAP>){
	if (/^\s*Total unperturbed charge:\s+(-{0,1}\d+\.\d+)/){
		$resid_charge = $1;
	}
}
#print "Re-assigning $resid\n";
foreach my $i (@high_layer){
        if (${$PDB[5]}[$i-1] eq $resid_num_i){
                my $count = $i;
                while (1){
                        if (${$PDB[5]}[$count-1] && ${$PDB[5]}[$count-1] eq $resid_num_i){
                                push(@full_resid_high_layer,$count);
                                $count++;
                        } else {
                                last;
                        }
                }
                $count = $i;
                while (1){
                        if (${$PDB[5]}[$count-1] && $count > 0 && ${$PDB[5]}[$count-1] eq $resid_num_i){
                                push(@full_resid_high_layer,$count);
                                $count--;
                        } else {
                                last;
                        }
                }
        }
}
@full_resid_high_layer = grep { ! $dummy{ $_ }++ } @full_resid_high_layer;
@full_resid_high_layer = sort {$a <=> $b} @full_resid_high_layer;
foreach my $i (@high_layer){
	#print "High layer atom: $i\n";
        if (${$PDB[5]}[$i-1] eq $resid_num_i){
		#print "Resid ID: ${$PDB[0]}[$i-1], Resid #: ${$PDB[5]}[$i-1]\n";
		my ($index) = grep {$full_resid_high_layer[$_] eq $i} 0..$#full_resid_high_layer;
		#print "Index: $index\n";
                push(@resid_high_layer,$index);
        }
}
#print "Full resid_high_layer: @full_resid_high_layer\n"; 
#print "Resid_high_layer: @resid_high_layer\n";
#print "Collecting residual charge.\n";
open (LIB, "<".$resid.".lib") or die;
my $model_charge = 0;
my $QM_model_charge = 0;
$_ = read_input("high layer charges");
$_ =~ s/^\s*//;
$_ =~ s/\s*$//;
@_ = split(/\s+/,$_);
my $check = 0;
for my $i (0..$#_/2){
	if ($_[2*$i+1] eq $resid_num_i){
		$QM_model_charge = $_[2*$i];
                $check = 1;
	}
}
if ($check == 0){
	print color("red"), "Residue $resid_num_i ($resid) is shared between two layers but the High Layer charge (should be integer) has not be specified by the user in cobram.parm under 'high layer charges'!\n", color("reset");
        exit;
}
print "QM charge: $QM_model_charge\n";
my $charge = 0;
my $abs_charge = 0;
while (<LIB>){
	if (/!entry\.$resid\.unit\.atoms table/){
		my $count = 0;
		while (my $tmp = <LIB>){
			if ($tmp =~ /!entry\.$resid\.unit\.atomspertinfo/){  
                                last;
                        }
                        @_ = split(/\s+/,$tmp);
			if (grep {$_ == $count} @resid_high_layer){
				$model_charge += $_[8];
			} elsif ($_[1] eq '"P"') {
			        $charge += $_[8];
			} else {
				$charge += $_[8];
				$abs_charge += abs($_[8]);
			}
			$count++;
		}
		last;
	}
}
$model_charge = $model_charge - $QM_model_charge;
printf "Total charge the residue should have: %9.6f\n", $resid_charge;
printf "Total charge the MM atoms have before redistribution of the residual charge: %9.6f\n", $charge;
printf "Total residual charge to redistribute among MM atoms: %9.6f\n", $model_charge;
$charge = 0;
close(LIB);
#$abs_charge -= abs($resid_charge);
open (INP, "<".$resid.".lib") or die;
open (LIB, ">".$new_resid->{$resid_num_i}.".lib") or die;
while (<INP>){
	if (/!entry\.$resid\.unit\.atoms table/){
		my $count = 0;
		print LIB "!entry.$new_resid->{$resid_num_i}.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg\n";
		while (my $tmp = <INP>){
			if ($tmp =~ /!entry\.$resid\.unit\.atomspertinfo/){
				print LIB "!entry.$new_resid->{$resid_num_i}.unit.atomspertinfo table  str pname  str ptype  int ptypex  int pelmnt  dbl pchg\n";
                                last;
                        }
			@_ = split(/\s+/,$tmp);
			if (grep {$_ == $count} @resid_high_layer){
                                print LIB "$_[0] $_[1] $_[2] $_[3] $_[4] $_[5] $_[6] $_[7] 0.000000\n";
                        } elsif ($_[1] eq '"P"') {
				printf LIB "$_[0] $_[1] $_[2] $_[3] $_[4] $_[5] $_[6] $_[7] $_[8]\n";
				$charge += $_[8];
			} else {
				printf LIB "$_[0] $_[1] $_[2] $_[3] $_[4] $_[5] $_[6] $_[7] %9.6f\n", $_[8]+abs($_[8])*$model_charge/$abs_charge;
				$charge += $_[8]+abs($_[8])*$model_charge/$abs_charge;
                        }
                        $count++;
		}
	} elsif (/(.*)$resid(.*)/){
		print LIB "$1$new_resid->{$resid_num_i}$2\n";
	} else {
		print LIB $_;
	}
}
close(INP);
close(LIB);
printf "Total charge the MM atoms have after redistribution of the residual charge: %9.6f\n", $charge;
print "... done!\n";
}

########################################

sub check_atom_types{
my (@model_types,$model_types,@real_types,@pos,@parm);
my ($line,$field,$sed,@tmp);
print "\nComparing atom types in real.top and model-H.top ...\n";
foreach my $i (@high_layer){
	$line = floor(($i-1)/20);
        $field = ($i-1)%20;
	open (INP, "<real.top") or die;
	while (<INP>){
	        if (/^\s*%FLAG AMBER_ATOM_TYPE/){
			$_ = <INP>;
			for my $j (0..$line){
				$_ = <INP>;
			}
			chomp;
			$_ =~ s/^\s*//;
			$_ =~ s/\s*$//;
			@_ = split(/\s+/,$_);
			if (scalar(split(//,$_[$field])) == 1){
				$_[$field] = $_[$field]." ";
			} elsif ((split(//,$_[$field]))[1] eq '*'){
				@tmp = split(//,$_[$field]);
	                        $_[$field] = $tmp[0]."\\".$tmp[1];
			}
			push(@real_types,$_[$field]);
			last;
		}
	}
	close(INP);
}
print "Atom types in real.top   : @real_types\n";
$model_types = types_from_model();
print "Atom types in model-H.top: @$model_types\n";
($sed,@pos) = compare_real_model(\@real_types,$model_types);
}	

#####################################

sub compare_real_model{
my $sed = 'no';
my @real_types = @{$_[0]};
my @model_types = @{$_[1]};
my @pos = ();
for my $i (0..$#model_types-scalar(@hmlinks)-scalar(@hllinks)){
	if ($real_types[$i] ne $model_types[$i]){
		printf color("yellow"), "real.top and model-H.top differ in atom type for atom $high_layer[$i](real)/%d(high): $real_types[$i] vs. $model_types[$i]\n", $i+1;
		printf color("reset");
		push(@pos,$i);
		$sed = 'yes';
	}
}
if (!@real_types){
	print color("red"), "Real atom types were not assign. Probably real.top was not created. Script terminates.\n";
	exit; 
}
if (!@model_types){
        print color("red"), "Model-H atom types were not assign. Probably model-h.top was not created. Script terminates.\n";
        exit;
}
if ($sed eq 'no'){
	print "Atom types of real.top and model-H.top identical.\n";
}
return ($sed,@pos);
}

#####################################

sub real_vs_libs{
my (@real_name,@libs_name,$tmp_real,$tmp_libs,@xyz,$index,$rearranged_links,@tmp_equivalent_names);
$print = 'no';
my $rearrange = 0; 
open REAL, "<real.pdb" or die;
open LIBS, "<real_libs.pdb" or die;
open XYZ, "<real_layers.xyz" or die;	
while(1){
	for my $resid (0..$#RES){
		for my $i (0..$#{$RES[$resid]}){
			$_ = <REAL>;
                        if ($_ =~ /^\s*ATOM/){
				$tmp_real = $_;
			} else {
                                $_ = <REAL>;
                                $tmp_real = $_;
                        }
			$_ = <LIBS>;
                        if ($_ =~ /^\s*ATOM/){
                                $tmp_libs = $_;
                        } else {
                                $_ = <LIBS>;
                                $tmp_libs = $_;
                        }
                        $_ = <XYZ>;
                        push(@{$xyz[$resid]},$_);			
			@_ = split(//,$tmp_real);
			$_ = join('',@_[12..15]);
			$_ =~ s/^\s*//;
	                $_ =~ s/\s*$//;
			@tmp_equivalent_names = generate_equivalent_names($_);
                        push(@{${$real_name[$resid]}[$i]},@{$tmp_equivalent_names[0]});
                        @_ = split(//,$tmp_libs);
			$_ = join('',@_[12..15]);
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        push(@{$libs_name[$resid]},$_);
		}
	}
	last;
}
close(REAL);
close(LIBS);
close(XYZ);
my $count = 0;
for my $resid (0..$#RES){
        my $res_count = 0;
        for my $i (0..$#{$RES[$resid]}){
                $res_count++;
                if (!(grep {$_ eq ${$libs_name[$resid]}[$i]} @{${$real_name[$resid]}[$i]})){
                        printf "real.pdb and real_libs.pdb differ for residue %d, atom %d: @{${$real_name[$resid]}[$i]} vs. ${$libs_name[$resid]}[$i]\n", $resid+1, $i+1;
			print color("red"), "This should not be the case, bacause the .pdb file was already reordered according to the .libs in the beginning. The script will terminate!\n";
			printf color("reset");
			exit;
                }
	}
}
}

#####################################

sub read_equivalent{
push(@{$equivalence[0]},0);
my %dummy;
my $found = 0;
for my $ff (@ff){
	if ( -e "$BASEDIR/$ff" ){
		open FF, "<$BASEDIR/$ff" or die "Check permissions for $ff";
	} elsif ( -e "$AMBERHOME/dat/leap/cmd/$ff" ){
		open FF, "<$AMBERHOME/dat/leap/cmd/$ff" or die "Check permissions for $AMBERHOME/dat/leap/cmd/$ff";
	} else {
		print color("red"), "Force field $ff could not be found in $BASEDIR or in $AMBERHOME/dat/leap/cmd/", color("reset");
		exit;
	}
	while (<FF>){
		if (/\s*addPdbAtomMap/){
			while(<FF>){
				$_ =~ s/^\s*//;
	                	$_ =~ s/\s*$//;
				@_ = split(/\s+/,$_);
				if ($_[0] eq '}'){
					last;
				} else {
					$_[1] =~ s/\"//g;
	                                $_[2] =~ s/\"//g;
					for my $i (0..$#equivalence){
						if ((grep {$_ eq $_[1]} @{$equivalence[$i]}) || (grep {$_ eq $_[2]} @{$equivalence[$i]})){
							push (@{$equivalence[$i]},$_[1],$_[2]);		# add equivalent notations to already existing notations
							$found = 1;
						} 
					}
					if ($found == 0){
						push (@{$equivalence[$#equivalence+1]},$_[1],$_[2]); # create new set of equivalent notations
					} else {
						$found = 0;
					}
				}
			}
		}
	}
}
for my $i (0..$#equivalence){
	@{$equivalence[$i]} = grep { ! $dummy{ $_ }++ } @{$equivalence[$i]};
}
shift(@equivalence);
print "List of equivalent notations:\n";
for my $i (0..$#equivalence){
	for my $j (0..$#{$equivalence[$i]}){
		print "${$equivalence[$i]}[$j] "; 
	}
	print "\n";
}
}

#####################################

sub types_from_model{
my (@model_types,@tmp,@_i);
open (INP, "<model-H.top") or die;
while (<INP>){
	if (/^\s*%FLAG AMBER_ATOM_TYPE/){
		$_ = <INP>;
		for my $i (0..floor((scalar(@high_layer)+scalar(@hmlinks)+scalar(@hllinks))/20.0001)){
			@_ = ();
			$_ = <INP>;
			chomp;
			@_i = split(//,$_);         
			for my $j (0..19){
				if ($_i[4*$j]){				
					$_[$j] = join('',@_i[4*$j..4*$j+3]);
					$_[$j] =~ s/^\s*//;
					$_[$j] =~ s/\s*$//;
				}
			}			
			for my $j (0..$#_){
				if (scalar(split(//,$_[$j])) == 1){
                                	$_[$j] = $_[$j]." ";
                        	} elsif ((split(//,$_[$j]))[1] eq '*'){
                                	@tmp = split(//,$_[$j]);
                                	$_[$j] = $tmp[0]."\\".$tmp[1];
                        	}
			}
			push(@model_types,@_);
		}
		last;
	}
}
close(INP);
return(\@model_types);
}

######################################

sub keep_files{
	print "\nSaving files to directory cobramm_files!\n";
	if ($switch == 1){
		my $coord = read_input("coordinate file");
        	system("if [ ! -d \$PWD/cobramm_files ]; then mkdir cobramm_files; fi; cp -f real_layers.xyz real.top model-H.top cobram.parm $coord crd_2_HM.map HighLayer.xyz HighMediumLayer.xyz cobramm_files/.");
		if (-e "$BASEDIR/$coord\_orig"){
			system("cp -f $coord\_orig cobramm_files/.");
		}
	} elsif ($switch == 3){
		my $coord = read_input("coordinate file");
		system("if [ ! -d \$PWD/cobramm_files ]; then mkdir cobramm_files; fi; cp -f real_layers.xyz real.top model-H.top cobram.parm $coord velocity.dat HighLayer.xyz HighMediumLayer.xyz cobramm_files/.");
		if (-e "$BASEDIR/$coord\_orig"){
                        system("cp -f $coord\_orig cobramm_files/.");
                }
		if (-e "$BASEDIR/rattle"){
                        system("cp -f rattle cobramm_files/.");
                }
	}	
}

#####################################

sub generate_log{
open (OUT, ">geometry.com") or die;
print OUT "%mem=1000MB\n %chk=geometry.chk\n #P freq(HPModes,ReadFC) geom=check guess=read nosymm\n\ntest\n\n1 1\n";
close(OUT);
system("$GAUSS_EXE < geometry.com > geometry.log");
}

#####################################

sub generate_molf{
my (@freqs,@redmass,@ints,@modes,@el,@pos,@x,@y,@z);
open (GAU, "<geometry.log") or die;
open (MOLF, ">molden.molf") or die;
print MOLF "ATOMIC GEOMETRYS(AU)\n";
my $nModes = 0;
while (<GAU>){
        if (/^\s+Z-Matrix orientation/ || /^\s+Input orientation/ || /^\s+Standard orientation/){
                $_ = <GAU>; $_ = <GAU>; $_ = <GAU>; $_ = <GAU>;
		while (1){
			$_ = <GAU>;
                        if (/------/){
				last;	
			}
			$_ =~ s/^\s*//;
        		$_ =~ s/\s*$//;
 			@_ = split(/\s+/,$_);
			$nAt = $_[0];
                        my ($el) = grep { $PTE{$_}[0] eq $_[1] } keys %PTE;
			push(@el,$el); push(@pos,$_[1]); push(@x,$_[3]); push(@y,$_[4]); push(@z,$_[5]); 
                }
		print MOLF " $nAt\n";
		for my $i (0..$nAt-1){
			printf MOLF "%5s %9d %15.8f %15.8f %15.8f\n", $el[$i], $pos[$i], $x[$i]/$au2ang, $y[$i]/$au2ang, $z[$i]/$au2ang;
		}
	}
        if (/^\s*Frequencies --/){
                $_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
                @_ = split(/\s+/,$_);
                my $length = $#_-1;
                for my $i (0..$length-1){
			push(@freqs,$_[$i+2]);
		}
	}
        if (/^\s*Red\. masses --/ || /^\s*Reduced masses --/){
		$_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
                @_ = split(/\s+/,$_);
		my $length = $#_-2;
                for my $i (0..$length-1){
                        push(@redmass,$_[$i+3]);
                }
	} 
	if (/^\s*IR Inten/){
                $_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
                @_ = split(/\s+/,$_);
                my $length = $#_-2;
                for my $i (0..$length-1){
                        push(@ints,$_[$i+3]);
                }
	}
	if (/^\s*Coord Atom Element/){
                my $length;
		for my $i (0..3*$nAt-1){
			$_ = <GAU>;
			$_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
			$length = $#_-2;
			for my $j (0..$length-1){
				if ($_[0] == 1){
					push(@{${$modes[$nModes+$j]}[0]},$_[$j+3]);
				} elsif ($_[0] == 2){
                                        push(@{${$modes[$nModes+$j]}[1]},$_[$j+3]);
                                } elsif ($_[0] == 3){
                                        push(@{${$modes[$nModes+$j]}[2]},$_[$j+3]);
                                }
			}
		}
		$nModes += $length;
		if ($nModes == 3*$nAt-6){
			last;
		}
	}
	if (/^\s*Atom  AN/ || /^\s*Atom AN/){
                my $length;
                for my $i (0..$nAt-1){
                        $_ = <GAU>;
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
                        $length = ($#_-1)/3;
                        for my $j (0..$length-1){
                        	push(@{${$modes[$nModes+$j]}[0]},$_[3*$j+2]);
                        	push(@{${$modes[$nModes+$j]}[1]},$_[3*$j+3]);
                        	push(@{${$modes[$nModes+$j]}[2]},$_[3*$j+4]);
                        }
                }
                $nModes += $length;
                if ($nModes == 3*$nAt-6){
                        last;
                }
        }
}
printf MOLF "%10d %9d\n", scalar(@freqs), scalar(@el);
print MOLF "FREQUENCY(cm**-1)\n";
for my $i (0..$#freqs){
	printf MOLF "%15.8f\n", $freqs[$i];
}
print MOLF "REDUCED MASS(AMU)\n";
for my $i (0..$#freqs){
        printf MOLF "%15.8f\n", $redmass[$i];
}
print MOLF "NORMAL MODE\n";
for my $i (1..$nModes){
        for my $j (0..$nAt-1){
		printf MOLF "%15.8f %15.8f %15.8f\n", ${${$modes[$i-1]}[0]}[$j]/sqrt($redmass[$i-1]), ${${$modes[$i-1]}[1]}[$j]/sqrt($redmass[$i-1]), ${${$modes[$i-1]}[2]}[$j]/sqrt($redmass[$i-1]);	
	}
}
$numModes=scalar(@freqs);
close(MOLF);
close(GAU);
}

#####################################

sub gen_init_cond{
my $iGeom = $_[0];
if ($sampling eq "Wigner"){
	#system("cp $BASEDIR/molden.molf $dirname");
	#generate_Wigner_input();
	#system("$COBRAMM_PATH/util/sampling.exe <wigner.inp >wigner.out");
	system("cp $BASEDIR/real_layers.xyz $dirname");
	if (-e "$BASEDIR/real.crd"){
		system("cp $BASEDIR/real.crd $dirname");
	}
	system("cp $BASEDIR/cobram.parm $dirname");
	system("cp $BASEDIR/geometry.fchk $dirname");
	$inact_modes = list_active_modes();
	$temp = read_input("temperature");
	$geomdir = sprintf("%s%s%03d", $BASEDIR, '/geom_', $iGeom);
	print "\n";
	system("echo \" WIGNER SAMPLING\"");
	if (scalar(@$inact_modes) == 0){
		system("echo \"$COBRAMM_PATH/util/cobramm-wignersampling.py prepare -n 1 -fd $geomdir -T $temp  > output.dat\"");
	} else {
		system("echo \"$COBRAMM_PATH/util/cobramm-wignersampling.py prepare -n 1 -fd $geomdir -T $temp -FM @$inact_modes > output.dat\"");
	}
	print "\n";
	system("$COBRAMM_PATH/util/cobramm-wignersampling.py prepare -n 1 -fd $geomdir -T $temp -FM @$inact_modes > output.dat");
	system("cp $geomdir/sampling/sample_0000/velocity.dat $geomdir/velocity.dat");
	system("cp $geomdir/sampling/sample_0000/geometry.xyz $geomdir/x_Q_ang.dat");
	system("rm -rf $geomdir/sampling $geomdir/geometry.fchk");
} elsif ($sampling eq "ZPE"){
	system("cp $BASEDIR/geometry.chk $dirname");
	my $rand = int(rand(1000000));
	open (OUT, ">gaussian.com") or die;
	print OUT "%mem=1000MB\n %chk=geometry.chk\n #P Nonstd
1/6=1,7=11,8=2500,9=1,10=5,18=40,29=2,38=1,42=1,44=$rand,70=3,80=1/1,18;\n99//99;\n\ntest\n\n1 1\n";
	close(OUT);
	print "\nGenerating the initial conditions ... ";
	system("$GAUSS_EXE < gaussian.com > gaussian.log");
	system("rm geometry.chk");
	print "done!\n";
}
}

#####################################

sub get_NM{
if (-e "$BASEDIR/geometry.chk.gz"){
	system("gunzip geometry.chk.gz");
	system("$GAUSS_DIR/formchk geometry.chk geometry.fchk");
} elsif (-e "$BASEDIR/geometry.chk"){
	system("$GAUSS_DIR/formchk geometry.chk geometry.fchk");
}
open (CHK, "<geometry.fchk") or die;
while (<CHK>){
	if (/^\s*Number of Normal Modes/){
                $_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
                @_ = split(/\s+/,$_);
                $_ = $_[5];
		last;
	}
}
close(CHK);
return $_;
}

#####################################

sub list_active_modes{
my (@range,@act_modes,@inact_modes,$temp);
$_ = read_input("active modes");
@_ = split(/\s+/,$_);
for my $i (0..$#_){
        @range = split(/-/,$_[$i]);
        if (!$range[1]){
                push(@act_modes,$range[0]);
        } else {
                for my $j ($range[0]..$range[1]){
                        push(@act_modes,$j);
                }
        }
        undef @range;
}
my %dummy = ();
@act_modes = grep { ! $dummy{ $_ }++ } @act_modes;
@act_modes = sort {$a <=> $b} @act_modes;
my $max_mode = $act_modes[-1];
print "Highest active mode: $max_mode\n";
@_ = ();
for my $i (1..$numModes){
        push(@_,$i);
}
my %diff;
@diff{ @act_modes }= ();
@inact_modes = grep !exists($diff{$_}), @_;
print "Inactive modes: ";
for my $i (0..$#inact_modes){
        print "$inact_modes[$i] ";
}
print "\n";

return(\@inact_modes);
}

#####################################
sub generate_Wigner_input{
my (@range,@act_modes,@inact_modes,$temp);
$_ = read_input("active modes");
@_ = split(/\s+/,$_);
for my $i (0..$#_){
        @range = split(/-/,$_[$i]);
        if (!$range[1]){
                push(@act_modes,$range[0]);
        } else {
                for my $j ($range[0]..$range[1]){
                        push(@act_modes,$j);
                }
        }
        undef @range;
}
my %dummy = ();
@act_modes = grep { ! $dummy{ $_ }++ } @act_modes;
@act_modes = sort {$a <=> $b} @act_modes;
my $max_mode = $act_modes[-1];
print "Highest active mode: $max_mode\n";
@_ = ();
for my $i (1..$numModes){
	push(@_,$i);
}
my %diff;
@diff{ @act_modes }= ();
@inact_modes = grep !exists($diff{$_}), @_;
print "Inactive modes: ";
for my $i (0..$#inact_modes){
	print "$inact_modes[$i] ";
}
print "\n";
$temp = read_input("temperature");
printf "Temperature in sampling procedure: %d\n", $temp;
open (MOLF, "<molden.molf") or die;
$_ = <MOLF>; $_ = <MOLF>;
$_ =~ s/^\s*//;
$_ =~ s/\s*$//;
$nAt = $_;
close(MOLF);
open (WIG, ">wigner.inp") or die;
print WIG "$nAt          read (*,*) n_atom
$numModes         read (*,*) n_mode
1           read (*,*) label_random
10000       read (*,*) nr
50          read (*,*) nbin
1           read (*,*) label_read_vib
3           read (*,*) label_es_output
molden.molf read (*,*) filename_es_output
1           read (*,*) label_displacement
1           read (*,*) n_geom\n";
if ($temp == 0){
	print WIG "1           read (*,*) label_method\n";
} else { 
	print WIG "2           read (*,*) label_method
$temp       read (*,*) t_k\n";
}
if (scalar(@inact_modes) == 0){
        print WIG "0           read (*,*) label_frozen
1           read (*,*) number_frozen 
1           read (*,*) list_frozen\n";
} else {
        printf WIG "1           read (*,*) label_frozen
%d           read (*,*) number_frozen\n", scalar(@inact_modes);
        for my $i (0..$#inact_modes){
                print WIG "$inact_modes[$i] ";
        }
        print WIG "         read (*,*) list_frozen\n";
}
close(WIG);
}

#####################################

sub replace_MM_HL{
my $iGeom = $_[0];
undef(%HL);
undef(%ML);
@LL=();
@HML=();
if ($heat == 0){
	$rst = read_input("MM geometry/velocity file");
} else {
	$rst = 'real.crd';
}
system("cp $BASEDIR/$rst $dirname");
open CRD, "<$rst" or die;
$_ = <CRD>;
$_ = <CRD>;
$_ =~ s/^\s*//;
$_ =~ s/\s*$//;
@_ = split(/\s+/,$_);
my $numAt = $_[0];
close(CRD);
#read_wigner();
if ($sampling eq "Wigner"){
	read_wigner();
} elsif ($sampling eq "ZPE"){
	read_gaussian();
}
if ($heat == 0){
	open (INP, "<".$BASEDIR."/crd_2_HM.map") or die;
	$_ = <INP>;
	while (<INP>){
		chomp;
		if ($_ eq ''){
			last;
		}
	        $_ =~ s/^\s*//;
	        $_ =~ s/\s*$//;
		@_ = split(/\s+/,$_);
		$HL{$_[0]} = $_[1];
	}
	$_ = <INP>;
        while (<INP>){
		chomp;
		if ($_ eq ''){
                        last;
                }
		$_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
		@_ = split(/\s+/,$_);
                $ML{$_[0]} = $_[1];
	}
	close(INP);
	@_ = keys %HL;
	push(@_,keys %ML);
	for my $i (1..$numAt){
		if ( !(grep {$_ == $i} keys @_) ){   # the keys of the HL, ML hashes are the positions of the H and M atoms in the gaussian.log (starting at 1), the values are the positions in the original .crd
			push (@LL,$i);
		}
	}
} else {
	my $count = 0;
	my $countHM = 0;
	open (INP, "<".$BASEDIR."/real_layers.xyz") or die;
	while (<INP>){
		$count++;
		$_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
		@_ = split(/\s+/,$_);
		if ($_[5] eq 'H'){
			$countHM++;
			$HL{$countHM} = $count;
		} elsif ($_[5] eq 'M'){
			$countHM++;
			$ML{$countHM} = $count;
		} elsif ($_[5] eq 'L') {
			push (@LL,$count);
		}
	}
	close(INP);
}
#print "HL: @HL\n";
#print "ML: @ML\n";
%HML = (%HL,%ML);
#%rHML = reverse %HML;
for my $key ( keys %HML ) {
	#print "$key -> $HML{$key}\n";
	push(@HML,$HML{$key});
}
@HML = sort {$a <=> $b} @HML; #this array will be used for the constraints on the H and M layer atoms during MD (order is as in the original .crd)
#print "HML: @HML\n";
close(INP);
$rst=replace($rst);
}

#####################################

sub replace{
my $dummy;
my @_i;
my $repl_vel = 0;
my $lrst = $_[0];
open (INP, "<".$lrst) or die;
my $count = 0;
open (OUT, ">tmp") or die;
$_ = <INP>;
print OUT "$_";
$_ = <INP>;
print OUT "$_";
$_ =~ s/^\s*//;
$_ =~ s/\s*$//;
($numAt,$dummy) = split(/\s+/,$_);
while (<INP>){
	@_ = ();
        $_ =~ s/\s*$//;
	@_i = split(//,$_);
	for my $i (0..5){
		if ($_i[12*$i]){
			$_[$i] = join('',@_i[12*$i..12*$i+11]);
			$_[$i] =~ s/^\s*//;
			$_[$i] =~ s/\s*$//;
		}
	}
	if ($count == $numAt && $heat == 0){
                $count = 0;
                $repl_vel = 1;
	}
	$count++;
	if (grep {$HML{$_} == $count} keys %HML){
		if ($repl_vel == 0){
			my ($key) = grep {$HML{$_} == $count} keys %HML;
			#print "key: $key\n";
			#my $value = $rHML{$key};
			#print "value: $value\n"; 
			$_ = sprintf("%12.7f%12.7f%12.7f", ${$xyz[$key-1]}[0],${$xyz[$key-1]}[1],${$xyz[$key-1]}[2]);
		} elsif ($repl_vel == 1) {
			$_ = sprintf("%12.7f%12.7f%12.7f", 0, 0, 0);
		}
	} else {
		$_ = sprintf("%12.7f%12.7f%12.7f", $_[0],$_[1],$_[2]);
	}
	if (!$_[3]){
		print OUT "$_\n";
		next;
	}
	$count++;
	if (grep {$HML{$_} == $count} keys %HML){
		if ($repl_vel == 0){
			my ($key) = grep { $HML{$_} == $count } keys %HML; 
			#print "key: $key\n";
			#my $value = $rHML{$key};
			#print "value: $value\n";
                	$_ = sprintf("$_%12.7f%12.7f%12.7f", ${$xyz[$key-1]}[0],${$xyz[$key-1]}[1],${$xyz[$key-1]}[2]);
		} elsif ($repl_vel == 1){
			$_ = sprintf("$_%12.7f%12.7f%12.7f", 0, 0, 0);
		}	
        } else {
                $_ = sprintf("$_%12.7f%12.7f%12.7f", $_[3],$_[4],$_[5]);
        }
	print OUT "$_\n";
}
close(INP);
close(OUT);
system("mv $lrst $lrst\_old; mv tmp $lrst");
system("head -n -1 $lrst > tmp; tail -n 1 $lrst\_old >> tmp; mv tmp $lrst");
return($lrst);
}

#####################################

sub read_wigner{
$mass = 0;
@at = ();
@xyz = ();
if ($mode eq 'gas-phase'){
	open (XYZ, ">real_layers.xyz") or die;
	open (XYZref, "<$BASEDIR/real_layers.xyz") or die;
}
print "XYZ:\n";
open (INP, "<x_Q_ang.dat") or die;
$_ = <INP>;
$_ =~ s/^\s*//;
$_ =~ s/\s*$//;
$nAt = $_;
$_ = <INP>;
for my $i (0..$nAt-1){
	$_ = <INP>;
	$_ =~ s/^\s*//;
	$_ =~ s/\s*$//;
	@_ = split(/\s+/,$_);
        ($at[$i],${$xyz[$i]}[0],${$xyz[$i]}[1],${$xyz[$i]}[2]) = ($_[0],$_[1],$_[2],$_[3]);
	printf "%2s %12.6f %12.6f %12.6f\n", $at[$i], ${$xyz[$i]}[0], ${$xyz[$i]}[1], ${$xyz[$i]}[2];
        if ($mode eq 'gas-phase'){
		$_ = <XYZref>;
		$_ =~ s/^\s*//;
                $_ =~ s/\s*$//;
		@_ = split(/\s+/,$_);
		if (scalar(@_) == 6){
			printf XYZ "%2s 0 %12.6f %12.6f %12.6f  %s\n", $at[$i], ${$xyz[$i]}[0],${$xyz[$i]}[1],${$xyz[$i]}[2], $_[5];
		} elsif (scalar(@_) == 8){
			printf XYZ "%2s 0 %12.6f %12.6f %12.6f  %s %s %10d \n", $at[$i], ${$xyz[$i]}[0],${$xyz[$i]}[1],${$xyz[$i]}[2], $_[5],  $_[6],  $_[7];
		}
	}
}
close(INP);
if ($mode eq 'gas-phase'){
	close(XYZ); close(XYZref);
        system("mkdir input; cp real_layers.xyz velocity.dat input/.");
}

print "Veloc:\n";
open (INP, "<velocity.dat") or die;
open (VEL, ">v_V_au.dat") or die;
printf VEL "$nAt\n\n";
for my $i (0..$nAt-1){
        $_ = <INP>;
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @_ = split(/\s+/,$_);
	(${$vel[$i]}[0],${$vel[$i]}[1],${$vel[$i]}[2]) = ($_[0],$_[1],$_[2]);
	printf "%2s %15.10f %15.10f %15.10f\n",$at[$i], ${$vel[$i]}[0], ${$vel[$i]}[1], ${$vel[$i]}[2];
	printf VEL "%2s %15.10f %15.10f %15.10f\n",$at[$i], ${$vel[$i]}[0], ${$vel[$i]}[1], ${$vel[$i]}[2];
}
close(INP);
close(VEL);
system("rm velocity.dat")
}

#####################################

sub read_gaussian{
$mass = 0;
@at = ();
@xyz0 = ();
open (INP, "<gaussian.log") or die;
while(<INP>){
	if (/^\s*Cartesian coordinates read from the checkpoint file/){
		$_ = <INP>;
		$_ = <INP>;
		while(<INP>){
			$_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
			if ($_[0] eq 'Recover'){
				$_ = <INP>;
	                	$_ =~ s/^\s*//;
				$_ =~ s/\s*$//;
	               	 	@_ = split(/\s+/,$_);
	               	 	$nAt = $_[1];
				last;
			} else {
				push(@at,$_[0]); push(@{$xyz0[0]},$_[1]); push(@{$xyz0[1]},$_[2]); push(@{$xyz0[2]},$_[3]);
				$mass += ${$PTE{$_[0]}}[1];
			}
		}
		print "Total mass: $mass\n";
		center_of_mass2();
        }
        if (/^\s*Symbolic Z-matrix/){
		$_ = <INP>;
                while(<INP>){
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
                        if (!$_[0]){
                                $_ = <INP>;
                                $_ =~ s/^\s*//;
                                $_ =~ s/\s*$//;
                                @_ = split(/\s+/,$_);
                                $nAt = $_[1];
                                last;
                        } else {
                                push(@at,$_[0]); push(@{$xyz0[0]},$_[2]); push(@{$xyz0[1]},$_[3]); push(@{$xyz0[2]},$_[4]);
                                $mass += ${$PTE{$_[0]}}[1];
                        }
                }
                print "Total mass: $mass\n";
                center_of_mass2();
	}
	if (/^\s*Cartesian coordinates: \(bohr\)/){
		print "XYZ:\n";
		for my $i (0..$nAt-1){
			$_ = <INP>;
			$_ =~ s/^\s*//;
        	        $_ =~ s/\s*$//;
	                @_ = split(/\s+/,$_);
			$_[3] =~ s/D/E/; $_[5] =~ s/D/E/; $_[7] =~ s/D/E/;
			(${$xyz[$i]}[0],${$xyz[$i]}[1],${$xyz[$i]}[2]) = ($_[3]*$au2ang+$cm[0],$_[5]*$au2ang+$cm[1],$_[7]*$au2ang+$cm[2]);
		 	printf "%2s %12.6f %12.6f %12.6f\n", $at[$i], ${$xyz[$i]}[0], ${$xyz[$i]}[1], ${$xyz[$i]}[2];
		}
	}
	if (/^\s*MW cartesian velocity: \(sqrt\(amu\)\*bohr\/sec\)/){
		print "Veloc:\n";
                for my $i (0..$nAt-1){
                        $_ = <INP>;
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
			$_[3] =~ s/D/E/; $_[5] =~ s/D/E/; $_[7] =~ s/D/E/;
                        (${$vel[$i]}[0],${$vel[$i]}[1],${$vel[$i]}[2]) = ($_[3]*$gau2au*sqrt(${$PTE{$at[$i]}}[1])/${$PTE{$at[$i]}}[1],$_[5]*$gau2au*sqrt(${$PTE{$at[$i]}}[1])/${$PTE{$at[$i]}}[1],$_[7]*$gau2au*sqrt(${$PTE{$at[$i]}}[1])/${$PTE{$at[$i]}}[1]);
			printf "%2s %12.6f %12.6f %12.6f\n",$at[$i], ${$vel[$i]}[0], ${$vel[$i]}[1], ${$vel[$i]}[2];
                }
        }
}
}

#####################################

sub center_of_mass2{
@cm = qw/0 0 0/;
for my $i (0..2){
	for my $j (0..$#{$xyz0[$i]}){
		$cm[$i] += ${$xyz0[$i]}[$j]*${$PTE{$at[$j]}}[1]/$mass;
		#print "@cm\n";
	}
}
print "Center of mass: @cm\n";
}

#####################################

sub analyze_sampling{
my ($totEn_modes,$avPotE_wig,$avKinE_wig,$avKinE_vel);

print "\n######## ANALYSIS OF THE WIGNER SAMPLING #########\n";

#obtain vibrational energy from normal 
my $freqdir = sprintf("%s%s", $BASEDIR, '/geom_001');
open (INP, "<".$freqdir."/output.dat") or die;
my $itotEn = 0;
while(<INP>){
	if (/frequencies used in the Wigner sampling/){
		while(1){
                	$_ = <INP>;
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
			if ($_[0] eq 'Total'){
				$totEn_modes = 1240*($itotEn/4)/10000000;
				printf "Energy in all active modes from frequency calculation: %7.2f eV\n", $totEn_modes;
				last;
			} else {
				for my $i (0..$#_){
					$itotEn += $_[$i];
				}
			}
		}
	}
}
close(INP);

#obtain average potential and kinetic energy from output.dat in geom_XXX folders
my $nGeoms = read_input("number of geometries");
$itotEn = 0;
my $jtotEn = 0;
for my $iGeom (1..$nGeoms){
	my $geomdir = sprintf("%s%s%03d", $BASEDIR, '/geom_', $iGeom);
	open (INP, "<".$geomdir."/output.dat") or die;
	while(<INP>){
	        if (/^\s*Total ePot/){
	        	$_ =~ s/^\s*//;
	                $_ =~ s/\s*$//;
	                @_ = split(/\s+/,$_);
			$itotEn += $_[3];
	        } 
		if (/^\s*Total eKin/){
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/\s+/,$_);
                        $jtotEn += $_[3];
                }
	}
	close(INP);
}
$avPotE_wig = $itotEn / $nGeoms;
$avKinE_wig = $jtotEn / $nGeoms;
printf "Average potential energy in all active modes from Wigner sampling: %7.2f eV\n", $avPotE_wig;
printf "Average kinetic energy in all active modes from Wigner sampling: %7.2f eV\n", $avKinE_wig;
#obtain average kinetic energy from v_V_au.dat in geom_XXX folders

$itotEn = 0;
for my $iGeom (1..$nGeoms){
        my $geomdir = sprintf("%s%s%03d", $BASEDIR, '/geom_', $iGeom);
        open (INP, "<".$geomdir."/v_V_au.dat") or die;
	$_ = <INP>;
	$_ = <INP>;
        while(<INP>){
        	$_ =~ s/^\s*//;
        	$_ =~ s/\s*$//;
        	@_ = split(/\s+/,$_);
        	$itotEn +=  0.5*($_[1]**2+$_[2]**2+$_[3]**2)*$amu2e*$PTE{$_[0]}[1]*$au2eV;
        }
        close(INP);
}
$avKinE_vel = $itotEn / $nGeoms;
printf "Average kinetic energy in all active modes from velocities: %7.2f eV\n\n", $avKinE_vel;

#compare all 4 values (should be identical;
my $error =  0;
if (abs($avPotE_wig - $avKinE_wig) > 0.1){
	print color("yellow"),  "Values for average kinetic and potential energy from Wigner sampling differ too much!\n", color("reset");
	$error = 1;	
}
if (abs($totEn_modes - $avPotE_wig) > 0.1 || abs($totEn_modes - $avKinE_wig) > 0.1){
	print color("yellow"), "Values for average kinetic / potential energies differ too much from expected value (1/4 quant per active mode)! This might not be a problem if sampling was done at finite temperature (in which case sampled values must be bigger)!\n", color("reset"); 
	$error = 1;
}
if (abs($avKinE_wig - $avKinE_vel) > 0.1){
	print color("yellow"), "Kinetic energy stored in velocities differs too much from average kinetic energy from Wigner sampling! These two values should be identical!\n", color("reset");
	$error = 1;
}
if ($error == 1){
	print color("yellow"), "Possible remedy: increase number of samples!\n", color("reset");
}

}

#####################################

sub run_MM_heat{
my @boxdim;
my $mmmd = read_input("MM input file");
my $top = read_input("MM topology file");
system("cp $BASEDIR/$mmmd $dirname/template.in");
system("cp $BASEDIR/$top $dirname");
open (INP, "<template.in") or die;
open my $out, ">tmp" or die;
my $ntr = 0;
my $ntt = 0;
while(<INP>){
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        if ($_ eq '&cntrl'){
                print $out "Heating\n&cntrl\n";
                while(<INP>){
                        $_ =~ s/^\s*//;
                        $_ =~ s/\s*$//;
                        @_ = split(/,/,$_);
                        for my $i (0..$#_){
                                $_[$i] =~ s/^\s*//;
                        }
                        if ($_[0] eq '/'){
				if ($ntr == 0 && $constrain == 1){ 
					print $out "ibelly = 1,\n"; 
					$ntr=1;
					bellymask($out);
					if ($ntt == 0 && $constrain == 1){
						print $out "ntt = 0,\n";  # if user didn't specify ntt in the .in file, then we use the default ntt=0
						$ntt = 1;
					}
				} elsif ($ntr == 0 && $constrain == 0){
					print $out "ntr = 1,\n"; 
					$ntr = 1;
				}
                                print $out "/\n";
                                last;
                        }
                        for my $i (@_){
                                my @tmp = split(/=/,$i);
                                $tmp[0] =~ s/\s*$//;
				if ($tmp[0] eq 'ntr'){
					if ($constrain == 1){
						print $out "ibelly = 1,\n"; 
						$ntr=1; 
						bellymask($out);
					} else {
						print $out "ntr = 1,\n"; $ntr = 1;
					}
				} elsif ($tmp[0] eq 'ntb'){
					print $out "ntb = 0,\nigb = 0,\n";
				} elsif ($tmp[0] eq 'ntp' || $tmp[0] eq 'pres0' || $tmp[0] eq 'taup'){
				} elsif ($tmp[0] eq 'nstlim'){
					print $out "nstlim = 5000,\n"
                                } elsif ($tmp[0] eq 'dt'){
					print $out "dt = 0.001,\n";
                                } elsif ($tmp[0] eq 'tempi'){
					print $out "tempi = 0.0,\n";
                                } elsif ($tmp[0] eq 'temp0'){
					print $out "temp0 = 50.0,\n";
                                } elsif ($tmp[0] eq 'ntx'){
					print $out "ntx = 1,\n";
                                } elsif ($tmp[0] eq 'irest'){
					print $out "irest = 0,\n";
                                } elsif ($tmp[0] eq 'ntpr'){
					print $out "ntpr = 100,\n";
                                } elsif ($tmp[0] eq 'ntt'){
					if ( $constrain == 1 && ($tmp[1] == 2 || $tmp[1] == 3) ){
                                                print $out "ntt = 2,\nig = -1,\n";
                                                $ntt = 1;
                                        } elsif ( $constrain == 1 && $tmp[1] != 2 && $tmp[1] != 3 ){
                                                print $out "$i,\n"; 
						$ntt =1;
                                        } elsif ( $constrain == 0 && ($tmp[1] == 2 || $tmp[1] == 3) ){
                                                print $out "$i,\nig = -1,\n";
                                                $ntt = 1;
                                        } elsif ( $constrain == 0 && $tmp[1] != 2 && $tmp[1] != 3 ){
                                                print $out "$i,\n";
                                                $ntt = 1;
                                        }
				} elsif ($tmp[0] eq 'ig' || $tmp[0] eq 'igb'){	
                                } elsif ($tmp[0] eq 'gamma_ln' && $constrain == 0){
					print $out "$i,\n";
				} else {
                                        print $out "$i,\n";
                                }
                        }
                }
        }
}
close(INP);
system("rm template.in");
if ($constrain == 0){
	print $out "Restrains\n100.0\nATOM $HML[0] ";
        my $index = -1;
        my $first = $HML[0];
        increment($index,$first,$out);
}
system("mv tmp heat0-50.in; cp heat0-50.in heat50-100.in; cp heat0-50.in heat100-150.in; cp heat0-50.in heat150-200.in; cp heat0-50.in heat200-250.in; cp heat0-50.in heat250-300.in");
for (my $i=0; $i<=200; $i+=50){
	my $j = $i+50;
	my $k = $i+100;
	system("sed 's/tempi = 0.0/tempi = $j.0/' heat$j-$k.in > tmp; mv tmp heat$j-$k.in");
	system("sed 's/temp0 = 50.0/temp0 = $k.0/' heat$j-$k.in > tmp; mv tmp heat$j-$k.in");
	system("sed 's/ntx = 1/ntx = 5/' heat$j-$k.in > tmp; mv tmp heat$j-$k.in");
	system("sed 's/irest = 0/irest = 1/' heat$j-$k.in > tmp; mv tmp heat$j-$k.in");
}
#open TOP, "<"."$BASEDIR/$top" or die;
#while (<TOP>){
#	if (/^\s*\%FLAG BOX_DIMENSIONS/){
#		$_ = <TOP>;
#		$_ = <TOP>;
#		$_ =~ s/^\s*//;
#	        $_ =~ s/\s*$//;
#		@_ = split(/\s+/,$_);
#		@boxdim = ($_[1],$_[2],$_[3]);
#		last;
#	}
#}
#close(TOP);
#open CRD, ">>".$rst or die;
#if ($numAt % 2 == 0){
#	printf CRD "%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", @boxdim, 90.0, 90.0, 90.0;
#} else { 
#	printf CRD "\n%12.7f%12.7f%12.7f%12.7f%12.7f%12.7f\n", @boxdim, 90.0, 90.0, 90.0;
#}
#close(CRD);
print "\nHeating the system to 300K in 6 stages of 5 ps ... (using $nproc cores) \n";
if ($constrain == 0){
	print "No periodic boundery conditions will be used! Harmonic restrains (100 kcal/mol) will be used to freeze the H and M layer atoms.\n";
        print color("yellow"), "Alternatively, H and M freezing through ibelly can be used. For this option set \$constrain = 1 in the script. Beware, that Langevin dynamics is incompatible with ibelly = 1 and ntt = 2 (Anderson coupling scheme) will be used instead.\n", color("reset");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i    heat0-50.in -p $top -c            $rst -r    heat0-50.rst -o    heat0-50.out -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i    heat0-50.in -p $top -c            $rst -r    heat0-50.rst -o    heat0-50.out -ref $rst");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i  heat50-100.in -p $top -c    heat0-50.rst -r  heat50-100.rst -o  heat50-100.out -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i  heat50-100.in -p $top -c    heat0-50.rst -r  heat50-100.rst -o  heat50-100.out -ref $rst");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat100-150.in -p $top -c  heat50-100.rst -r heat100-150.rst -o heat100-150.out -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat100-150.in -p $top -c  heat50-100.rst -r heat100-150.rst -o heat100-150.out -ref $rst");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat150-200.in -p $top -c heat100-150.rst -r heat150-200.rst -o heat150-200.out -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat150-200.in -p $top -c heat100-150.rst -r heat150-200.rst -o heat150-200.out -ref $rst");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat200-250.in -p $top -c heat150-200.rst -r heat200-250.rst -o heat150-200.out -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat200-250.in -p $top -c heat150-200.rst -r heat200-250.rst -o heat150-200.out -ref $rst");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat250-300.in -p $top -c heat200-250.rst -r heat250-300.rst -o heat250-300.out -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat250-300.in -p $top -c heat200-250.rst -r heat250-300.rst -o heat250-300.out -ref $rst");
} else {
	print "No periodic boundery conditions will be used! A bellymask will be used to freeze the H and M layer atoms. Beware, that Langevin dynamics is incompatible with ibelly = 1 and ntt = 2 (Anderson coupling scheme) will be used instead.\n";
	print color("yellow"), "Alternatively, restrains (100 kcal/mol, via ntr = 1) can be used. For this option set \$constrain = 0 in the script\n", color("reset");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i    heat0-50.in -p $top -c            $rst -r    heat0-50.rst -o    heat0-50.out");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i    heat0-50.in -p $top -c            $rst -r    heat0-50.rst -o    heat0-50.out");
        system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i  heat50-100.in -p $top -c    heat0-50.rst -r  heat50-100.rst -o  heat50-100.out");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i  heat50-100.in -p $top -c    heat0-50.rst -r  heat50-100.rst -o  heat50-100.out");
        system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat100-150.in -p $top -c  heat50-100.rst -r heat100-150.rst -o heat100-150.out");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat100-150.in -p $top -c  heat50-100.rst -r heat100-150.rst -o heat100-150.out");
        system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat150-200.in -p $top -c heat100-150.rst -r heat150-200.rst -o heat150-200.out");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat150-200.in -p $top -c heat100-150.rst -r heat150-200.rst -o heat150-200.out");
        system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat200-250.in -p $top -c heat150-200.rst -r heat200-250.rst -o heat150-200.out");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat200-250.in -p $top -c heat150-200.rst -r heat200-250.rst -o heat150-200.out");
        system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat250-300.in -p $top -c heat200-250.rst -r heat250-300.rst -o heat250-300.out");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i heat250-300.in -p $top -c heat200-250.rst -r heat250-300.rst -o heat250-300.out");
}
$rst = 'heat250-300.rst';
print "done!\n";
}

#####################################

sub bellymask{
my $pid = $_[0];
print $pid "bellymask = '\@$LL[0]";
my $prev = $LL[0];
my $pprev = -1;
my ($curr);
for my $i (1..$#LL){
	$curr = $LL[$i];
	if ($curr-$prev == 1){
		if ($i == $#LL){
                        print $pid "-$curr";
                } else {
                        $pprev = $prev;
                        $prev = $curr;
                }
	} elsif ($curr-$prev != 1 && $prev-$pprev == 1) {
		print $pid "-$prev,$curr";
		$pprev = $prev;
		$prev = $curr;
	} elsif ($curr-$prev != 1 && $prev-$pprev != 1){
		print $pid ",$curr";
                $pprev = $prev;
                $prev = $curr;
	}
}	
print $pid "',\n";
}

#####################################

sub run_MM_constr{
my $mmmd = read_input("MM input file");
my $top = read_input("MM topology file");
system("cp $BASEDIR/$mmmd $dirname");
system("cp $BASEDIR/$top $dirname");
open (INP, "<".$mmmd) or die;
open my $out, ">tmp" or die;
my $ntr = 0;
my $ig = 0;
my $ntt = 0;
while(<INP>){
	$_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
	if ($_ eq '&cntrl'){
		print $out "Heating\n&cntrl\n";
		while(<INP>){
			$_ =~ s/^\s*//;
        	        $_ =~ s/\s*$//;
			@_ = split(/,/,$_);
			for my $i (0..$#_){
				$_[$i] =~ s/^\s*//;
			}
			if ($_[0] eq '/'){
                                if ($ntr == 0 && $constrain == 1){
                                        print $out "ibelly = 1,\n";
                                        $ntr=1;
                                        bellymask($out);
					if ($ntt == 0 && $constrain == 1){
                                                print $out "ntt = 0,\n"; # default 
                                                $ntt = 1;
                                        }
                                } elsif ($ntr == 0 && $constrain == 0){
                                        print $out "ntr = 1,\n";
                                        $ntr = 1;
                                }
				print $out "/\n";
				last;
			}
			for my $i (@_){
				my @tmp = split(/=/,$i);
				$tmp[0] =~ s/\s*$//;
				if ($tmp[0] eq 'ntr'){
					if ($constrain == 1){
						print $out "ibelly = 1,\n"; 
						$ntr=1; 
						bellymask($out);
					} else {
						print $out "ntr = 1,\n"; 
						$ntr = 1;
					}
				} elsif ($tmp[0] eq 'ntb'){
					print $out "ntb = 0,\nigb = 0,\n";
				} elsif ($tmp[0] eq 'ntp' || $tmp[0] eq 'pres0' || $tmp[0] eq 'taup'){
				} elsif ($tmp[0] eq 'nstlim'){
					print $out "nstlim = 5000,\n";
                                } elsif ($tmp[0] eq 'dt'){
					print $out "dt = 0.002,\n";
                                } elsif ($tmp[0] eq 'ntwx'){
					print $out "ntwx = 100,\n";
                                } elsif ($tmp[0] eq 'ntpr'){
					print $out "ntpr = 100,\n";
                                } elsif ($tmp[0] eq 'ntt'){
					if ( $constrain == 1 && ($tmp[1] == 2 || $tmp[1] == 3) ){
                                                print $out "ntt = 2,\nig = -1,\n";
                                                $ntt = 1;
                                        } elsif ( $constrain == 1 && $tmp[1] != 2 && $tmp[1] != 3 ){
                                                print $out "$i,\n"; 
						$ntt =1;
                                        } elsif ( $constrain == 0 && ($tmp[1] == 2 || $tmp[1] == 3) ){
                                                print $out "$i,\nig = -1,\n";
                                                $ntt = 1;
                                        } elsif ( $constrain == 0 && $tmp[1] != 2 && $tmp[1] != 3 ){
                                                print $out "$i,\n";
                                                $ntt = 1;
                                        }
				} elsif ($tmp[0] eq 'ig' || $tmp[0] eq 'igb'){
                                } elsif ($tmp[0] eq 'gamma_ln' && $constrain == 0){
					print $out "$i,\n";
				} elsif ($tmp[0] eq 'ioutfm'){
                                        print $out "ioutfm = 0,\n";
			        } elsif ($tmp[0] eq 'ntwv' || $tmp[0] eq 'ntwf'){
				} elsif ($tmp[0] eq 'ntxo'){
                                        print $out "ntxo = 1,\n";
                                } else {
					print $out "$i,\n";
				}
			}
		}
	}
}
close(INP);
if ($constrain == 0){
	print $out "Restrains\n100.0\nATOM $HML[0] ";
	my $index = -1;
	my $first = $HML[0];
	increment($index,$first,$out);
}
system("mv tmp $mmmd");
print "\nEquilibrating system for 10ps (using $nproc cores) ... \n";
if ($constrain == 0){
	print "No periodic boundery conditions will be used! Periodic restrains (100 kcal/mol) will be used to freeze the H and M layer atoms.\n";
        print color("yellow"), "Alternatively, H and M freezing through ibelly can be used. For this option set \$constrain = 1 in the script. Beware, that Langevin dynamics is incompatible with ibelly = 1 and ntt = 2 (Anderson coupling scheme) will be used instead.\n", color("reset");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i $mmmd -p $top -c $rst -r $eqrst -o equilibration_10ps.out -x equilibration.mdcrd -ref $rst");
	system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i $mmmd -p $top -c $rst -r $eqrst -o equilibration_10ps.out -x equilibration.mdcrd -ref $rst");
} else {
	print "No periodic boundery conditions will be used! A bellymask will be used to freeze the H and M layer atoms. Beware, that Langevin dynamics is incompatible with ibelly = 1 and ntt = 2 (Anderson coupling scheme) will be used instead.\n";
        print color("yellow"), "Alternatively, restrains (100 kcal/mol, via ntr = 1) can be used. For this option set \$constrain = 0 in the script\n", color("reset");
	system("echo $MPIrun $AMBERHOME/bin/sander$MPI -O -i $mmmd -p $top -c $rst -r $eqrst -o equilibration_10ps.out -x equilibration.mdcrd");
        system("     $MPIrun $AMBERHOME/bin/sander$MPI -O -i $mmmd -p $top -c $rst -r $eqrst -o equilibration_10ps.out -x equilibration.mdcrd");
}
print "done!\n";
if ($constrain == 0){
	$eqrst=replace($eqrst);
}
print "\nGenerating a .pdb file ... ";
make_pdb();
system("$AMBERHOME/bin/cpptraj $BASEDIR/$top < ptraj_pdb.in > ptraj_pdb.out");
remake_real_layers();
#system("mv $eqpdb.1 $eqpdb");
print "done!\n";
}

#####################################

sub make_pdb{
open my $inp, ">ptraj_pdb.in" or die;
print $inp "trajin $eqrst\ntrajout $eqpdb pdb\n";
reimage_mask($inp);
print $inp "image";
close($inp);
}

#####################################

sub remake_real_layers{
@xyz = ();
open (XYZ, ">real_layers.xyz") or die;
open (XYZref, "<$BASEDIR/real_layers.xyz") or die;
open (PDB, "<".$eqpdb) or die;
while(<PDB>){
    $_ =~ s/^\s*//;
    $_ =~ s/\s*$//;
    @_ = split(/\s+/,$_);
    if ($_[0] ne 'TER' && $_[0] ne 'END' && $_[5] && $_[6] && $_[7]){
	($xyz[0],$xyz[1],$xyz[2]) = ($_[5],$_[6],$_[7]);
        $_ = <XYZref>;
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @_ = split(/\s+/,$_);
        if (scalar(@_) == 6){
            printf XYZ "%2s     0 %12.6f %12.6f %12.6f  %s\n", $_[0], $xyz[0], $xyz[1], $xyz[2], $_[5];
        } elsif (scalar(@_) == 8){
            printf XYZ "%2s     0 %12.6f %12.6f %12.6f  %s %s %10d \n", $_[0], $xyz[0], $xyz[1], $xyz[2], $_[5], $_[6], $_[7];
        }
    }
}
close(XYZ);
close(XYZref);
close(PDB);
}

#####################################

sub reimage_mask{
my $pid = $_[0];
print $pid "center \@$HML[0]";
my $prev = $HML[0];
my $pprev = -1;
my ($curr);
for my $i (1..$#HML){
        $curr = $HML[$i];
        if ($curr-$prev == 1){
                if ($i == $#HML){
                        print $pid "-$curr";
                } else {
                        $pprev = $prev;
                        $prev = $curr;
                }
        } elsif ($curr-$prev != 1 && $prev-$pprev == 1) {
                print $pid "-$prev,$curr";
                $pprev = $prev;
                $prev = $curr;
        } elsif ($curr-$prev != 1 && $prev-$pprev != 1){
                print $pid ",$curr";
                $pprev = $prev;
                $prev = $curr;
        }
}
print $pid "\n";
}

#####################################

sub increment{
my $index = shift;
my $first = shift;
my $out = shift;
$index++;
if ($HML[$index+1] && $HML[$index] == $HML[$index+1]-1){
	increment($index,$first,$out);
} elsif ($HML[$index+1]) {
	if ($HML[$index] == $first){
		print $out "\nEND\nRestrains\n100.0\nATOM $HML[$index+1] ";
	} else {
		print $out "$HML[$index]\nEND\nRestrains\n100.0\nATOM $HML[$index+1] ";
	}
	$first = $HML[$index+1];
	increment($index,$first,$out);
} else {
	if ($HML[$index] == $first){
                print $out "\nEND\nEND";
        } else {
                print $out "$HML[$index]\nEND\nEND";
        }
	return;
}
}

#####################################

sub create_velocity{
open (INP, "<".$eqrst) or die;
my $count = 0;
print "Generating velocity.dat ... ";
open (OUT, ">velocity.dat") or die;
$_ = <INP>;
$_ = <INP>;
for my $i (1..($numAt+1)/2){
	$_ = <INP>;
}
while (<INP>){
        $_ =~ s/^\s*//;
        $_ =~ s/\s*$//;
        @_ = split(/\s+/,$_);
	if ($count == $numAt){
		last;
	}
        $count++;
        if (grep {$_ eq $count} @HML){
                my ($index) = grep { $HML[$_] eq $count } 0..$#HML; #search for index of value in array
                printf OUT "%12.7f%12.7f%12.7f\n", ${$vel[$index]}[0],${$vel[$index]}[1],${$vel[$index]}[2];
        } else {
                printf OUT "%12.7f%12.7f%12.7f\n", $amb2au*$_[0],$amb2au*$_[1],$amb2au*$_[2];
        }
        if (!$_[3]){
                next;
        }
        $count++;
        if (grep {$_ eq $count} @HML){
                my ($index) = grep { $HML[$_] eq $count } 0..$#HML; #search for index of value in array
                printf OUT "%12.7f%12.7f%12.7f\n", ${$vel[$index]}[0],${$vel[$index]}[1],${$vel[$index]}[2];
        } else {
                printf OUT "%12.7f%12.7f%12.7f\n", $amb2au*$_[3],$amb2au*$_[4],$amb2au*$_[5];
        }
}
print "done!\n";
}

#####################################

sub edit_velocity{
print "Editting velocity.dat ... ";
open (INP, "<velocity.dat") or die;
open (OUT, ">tmp") or die;
my $count = 0;
while(<INP>){
	$count++;
	if ((grep {$_ eq $count} @high_layer) || (grep {$_ eq $count} @medium_layer)){
		print OUT $_;
	}
}
close(INP);
for my $i (1..scalar(@hmlinks)+scalar(@hllinks)){
	print OUT "   0.0000000   0.0000000   0.0000000\n";
}
close(OUT);
system("mv tmp velocity.dat");
print "done!\n";
}

#####################################

sub rattle_h2o_cl3{ 
my @mlist = ();
my @molecule_pdb;
my @molecule_real;
print "Generating rattle ... ";
open (OUT, ">rattle") or die;
print OUT "!RATTLE\n";
for my $i (@{$sol_med_pdb2real[0]}){
	if (${$PDB[0]}[$i-1] eq "WAT" || ${$PDB[0]}[$i-1] eq "CL3"){
		@molecule_pdb = complete_solvent($i);
		@molecule_real = ();
		for my $k (@molecule_pdb){
			my ($index) = grep { ${$sol_med_pdb2real[0]}[$_] eq $k } 0..$#{$sol_med_pdb2real[0]};
			push(@molecule_real,${$sol_med_pdb2real[1]}[$index]);
		}
		for my $j (@molecule_real){
			if (grep {$_ eq $j} @mlist){
				last; 
			} elsif (${$PDB[0]}[$i-1] eq "WAT"){
				print OUT "$molecule_real[0] $molecule_real[1]\n$molecule_real[1] $molecule_real[2]\n$molecule_real[0] $molecule_real[2]\n";
				push(@mlist,@molecule_real);
			} elsif (${$PDB[0]}[$i-1] eq "CL3"){
                                print OUT "$molecule_real[0] $molecule_real[1]\n$molecule_real[1] $molecule_real[2]\n$molecule_real[1] $molecule_real[3]\n$molecule_real[1] $molecule_real[4]\n$molecule_real[0] $molecule_real[2]\n$molecule_real[2] $molecule_real[3]\n$molecule_real[3] $molecule_real[4]\n$molecule_real[0] $molecule_real[3]\n$molecule_real[2] $molecule_real[4]\n";
                                push(@mlist,@molecule_real);
                        }
		}
	}
}
print OUT "?RATTLE\n";
print "done!\n";
}

########################################

sub prepare_switch3{
`mkdir input; cp velocity.dat ../cobram.parm $eqpdb input/.; cd input; sed "/coordinate file/ccoordinate file = $eqpdb" cobram.parm > tmp; mv tmp cobram.parm`;
}

