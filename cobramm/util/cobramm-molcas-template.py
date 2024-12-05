#!/usr/bin/env python3
# coding=utf-8

#    COBRAMM
#    Copyright (c) 2019 ALMA MATER STUDIORUM - Università di Bologna

#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.

#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#####################################################################################################
import os
import sys         # System-specific parameters and functions
# imports of local modules
try:
    sys.path.append(os.path.join(os.environ['COBRAM_PATH'], 'cobramm'))
except KeyError:
    raise RuntimeError('cannot import cobramm module.\n'
                       'Please check if $COBRAM_PATH is defined in your environment, if not run:\n'
                       'export COBRAM_PATH=/path/to/COBRAMM')
import logwrt
import questions
###################################################################

def main():
    logwrt.cobramstart()

    print("\n\n ===============================================================================================\n\n"
              "             GENERATOR OF MOLCAS TEMPLATE TO BE INCLUDED INTO THE COBRAM.COMMAND               \n\n"
        " ==============================================================================================="
        "\n\n")

    print("\n\nInsert the options to write a template for COBRAMM/MOLCAS calculations\n\n")


    ### write the commands for the reference calculation
    reference_string = ""
    reference = questions.ask("\nReference wavefunction [options: casscf, rasscf]: ", qtype="string", constrain=True, accepted=["casscf", "rasscf"])
    chosen_reference = (reference.upper().replace("-", ""))
    spin = questions.ask("\nMultiplicity:  ", qtype="int")
    actEl = questions.ask("\nNumber of active electrons: ", qtype="int")
    ciroot = questions.ask("\nNumber of CI roots: ", qtype="int")
    Inactive = questions.ask("\nNumber of inactive orbitals: ", qtype="int")
    if chosen_reference == "CASSCF":
        ras1 = 0
        ras2 = questions.ask("\nNumber of active orbitals: ", qtype="int")
        ras3 = 0
        hole = 0
        part = 0
    elif chosen_reference == "RASSCF":
        ras1 = questions.ask("\nNumber of RAS1 orbitals: ", qtype="int")
        ras2 = questions.ask("\nNumber of RAS2 orbitals: ", qtype="int")
        ras3 = questions.ask("\nNumber of RAS3 orbitals: ", qtype="int")
        hole = questions.ask("\nMaximum number of holes in RAS1: ", qtype="int")
        part = questions.ask("\nMaximum number of electrons in RAS3: ", qtype="int")
    reference_string += "&RASSCF\nSymm=1\nSpin={0}\nNactel={1} {7} {8}\nCIRoot={2} {2} {0}\n" \
                    "Inacti={3}\nRas1={4}\nRas2={5}\nRas3={6}\nLumorb\n\n".format(spin, actEl, ciroot, Inactive, ras1, ras2, ras3, hole, part)

    ### write the commands for the pt2 calculation

    pt2 = questions.ask("\nType of PT2 correction: [options: none, ss-pt2, xms-pt2, ms-pt2, rms-pt2]: ", qtype="string",
                          constrain=True, accepted=["none", "ss-pt2", "xms-pt2", "ms-pt2", "rms-pt2"])

    chosen_pt2 = (pt2.upper().replace("-", ""))

    pt2string = ""
    if chosen_pt2 == "XMSPT2" or chosen_pt2 == "MSPT2" or chosen_pt2 == "RMSPT2" or chosen_pt2 == "SSPT2":
        imag = questions.ask("\nImaginary level shift: ", qtype="float")
        ipea = questions.ask("\nIPEA level shift: ", qtype="float")
        nroots = questions.ask("\nNumber of PT2 roots: ", qtype="int")
        if chosen_pt2 == "XMSPT2":
            ms = "xmultistate={0} {1}".format(nroots, [state for state in range(1, nroots+1)]).replace("[", "").replace("]", "").replace(",", "")
        elif chosen_pt2 == "MSPT2":
            ms = "multistate={0} {1}".format(nroots, [int(state) for state in range(1, nroots+1)]).replace("[", "").replace("]", "").replace(",", "")
        elif chosen_pt2 == "RMSPT2":
            ms = "rmultistate={0} {1}".format(nroots, [int(state) for state in range(1, nroots+1)]).replace("[", "").replace("]", "").replace(",", "")
        elif chosen_pt2 == "SSPT2":
            ms = "multistate={0} {1}\nNOMULT".format(nroots, [int(state) for state in range(1, nroots+1)]).replace("[", "").replace("]", "").replace(",", "")

        pt2string += "&CASPT2\nimag={0}\nipea={1}\n{2}\nPROP\n".format(imag, ipea, ms)

    transientstring = ""

    transient = questions.ask("\nDo you want to setup a pump-probe calculation?: ", qtype="bool")
    if transient:
        nofstate = questions.ask("\nHow many state do you want to include?: ", qtype="int")
        transientstring += "\n############# SECOND CALCULATION. DO NOT MODIFY MANUALLY\n\n!molcas\n\n\n&RASSCF\nSymm=1\nSpin={0}\nNactel={1} {2} {3}" \
                          "\nCIRoot={4} {4} 1\nInacti={5}\nRas1={6}\nRas2={7}\nRas3={8}\nLumorb\n" \
                          ">>>COPY molcas.JobIph expand.JobIph\n\n&caspt2\nimag={9}\nipea={10}\nmaxiter=200\n" \
                          "multistate={4} {11}\nnomult\nPROP\n\n>>>COPY expand.JobIph JOB001\n&rassi\nNr of JobIphs" \
                          "\n1 \nALL\nONEL\nMEIN\nProperties=3; 'Mltpl 1' 1 'Mltpl 1' 2 'Mltpl 1' 3\n\n" \
                          ">>>COPY ../vanilla.JobMix JOB001\n>>>COPY expand.JobIph JOB002\n\n&rassi\nNr of JobIphs\n2\nALL\nSTOVERLAP\n?molcas\n\n"\
            .format(spin, actEl, hole, part, nofstate, Inactive, ras1, ras2, ras3, imag, ipea, [int(state) for state in range(1, nofstate+1)]).replace("[", "").replace("]", "").replace(",", "")

    ### write the command for ANO basis sets
    basis_set_string = ""
    basis_set = questions.ask("\nDo you want to use Atomic Natural Orbitals (ANO) basis set?: ", qtype="bool")
    if basis_set:
        first_contraction = questions.ask("\nPlease enter the contraction for the first set of atoms \n\n("
                                  "format: {atom number} {type of ANO b.s.} {contraction})\n\n eg. 3,5,9,12 ANO-L 2s1p  or 1-10 ANO-L 3s2p1d\n\n", qtype="string")
        first_contraction += "\n"
        other_contractions = ""
        more_contraction = True
        while more_contraction:
            contraction = questions.ask("\nPlease enter the contraction for the next set of atoms "
                                       "(press STOP for adding no more):\n", qtype="string")
            if contraction == "STOP":
                more_contraction = False
            else:
                other_contractions += "{}\n".format(contraction)
        basis_set_string = "!basisset\n{0}{1}?basisset".format(first_contraction, other_contractions)

    template_str = "!molcas\n\n" + reference_string + pt2string + " " + "\n?molcas\n\n"+ transientstring + basis_set_string

    ### write the file
    template = "command_molcas_template"
    if os.path.isfile(template) == True:
        overwrite = questions.ask("\nOverwrite 'command_molcas_template' file?: ", qtype="bool")
        if overwrite == False:
            template = questions.ask("\nname of the new file: ", qtype="string")
    with open(template, "w") as templ:
        templ.write(template_str)

    print("\n{} file written\n".format(template))

##############################################################################################################
if __name__ == '__main__':
    main()

