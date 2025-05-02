#! /usr/bin/env python3
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

# import statements of module from python standard library

import math  # mathematical functions
import copy
# imports of local modules

from kids_log import Timing, Log  # manages log file output + start/end procedures

# math libraries

import numpy as np  # numpy library for scientific computation


#####################################################################################################


class Charge:
    """ the class implements the object to work with QM/MM charges during a COBRAMM calculation.
        The object stores the atom charges (for both the QM and MM parts) and defines methods to work
        with them: access the data in various forms, update them during the calculation,
        process charges to redistribute the excess charge of atoms at the H/M and H/L boundary
        according to the atom link scheme """

    ###########################################################################################
    #                                   CONSTRUCTOR
    ###########################################################################################

    def __init__(self, geometry, CRG_real, CRG_model_H):    
        """ constructor of the charge class, takes as input ... and ... """
        # store a local copy of the atom lists
        self.__list_HIGH = geometry.list_HIGH
        self.__list_MEDIUM = geometry.list_MEDIUM
        self.__list_LOW = geometry.list_LOW
        self.__list_MM = geometry.list_MM
        self.__list_QM = geometry.list_QM
        self.__atomLabel = geometry.atomLabel

        self.CRG_real = CRG_real # temporary saved for comparision

        TotCRG_modelH_tmp = np.sum(CRG_model_H)
        if TotCRG_modelH_tmp != round(TotCRG_modelH_tmp):
            Log.writeLog("Total from model-H top is not an integer!\n")
            if (TotCRG_modelH_tmp-round(TotCRG_modelH_tmp)) < 0.005:
                Log.writeLog("Try to rescale charges from modelH from {} to {}!\n".format(
                    (TotCRG_modelH_tmp), round(TotCRG_modelH_tmp)
                    ))
                scale = round(TotCRG_modelH_tmp) / (TotCRG_modelH_tmp)
                CRG_model_H = [q*scale for q in CRG_model_H]
            else:
                Log.fatalError("cannot rescale the model-H charges\n")

        # first define CRG_real, after redistribution and scaling
        CRG_real_notscaled = self.embeddingCharges(geometry, CRG_real, CRG_model_H)

        self.CRG_real = np.array(self.scaledCharges(geometry, CRG_real_notscaled))

        # store initial modelH charges
        self.CRG_model_H = np.array(CRG_model_H)
        # exit(0)

        
    ###########################################################################################
    #                          GET METHODS FOR CHARGES OF REAL
    ###########################################################################################

    @property
    def CRG_emb(self):
        """
        property method for embedding charges
        :return: a numpy array with the charges of the embedding
        """
        return np.array([self.CRG_real[j - 1] for j in self.__list_MM])

    @property
    def CRG_model(self):
        """
        property method for charges of the QM region
        :return: a numpy array with the charges of the QM region
        """
        return np.array([self.CRG_real[j - 1] for j in self.__list_QM])

    @property
    def CRG_pod(self):
        """
        property method for charges of the MM region
        :return: a numpy array with the charges of the MM region
        """
        return np.array([self.CRG_real[j - 1] for j in self.__list_MM])

    @property
    def CRG_HIGH(self):
        """
        property method for charges of the H region
        :return: a numpy array with the charges of the H region
        """
        return np.array([self.CRG_real[j - 1] for j in self.__list_HIGH])

    @property
    def CRG_MEDIUM(self):
        """
        property method for charges of the M region
        :return: a numpy array with the charges of the M region
        """
        return np.array([self.CRG_real[j - 1] for j in self.__list_MEDIUM])

    @property
    def CRG_LOW(self):
        """
        property method for charges of the L region
        :return: a numpy array with the charges of the L region
        """
        return np.array([self.CRG_real[j - 1] for j in self.__list_LOW])

    @property
    def TotCRG_real(self):
        """
        property method for the total charge of real
        :return: a floating point number with the sum of CRG_real
        """
        return np.add.reduce(self.CRG_real)

    @property
    def TotCRG_emb(self):
        """
        property method for the total charge of the embedding
        :return: a floating point number with the sum of the embedding charges
        """
        return sum([self.CRG_real[j - 1] for j in self.__list_MM])

    @property
    def TotCRG_model(self):
        """
        property method for the total charge of the QM part
        :return: a floating point number with the sum of the QM charges
        """
        return sum([self.CRG_real[j - 1] for j in self.__list_QM])

    @property
    def TotCRG_pod(self):
        """
        property method for the total charge of the MM part
        :return: a floating point number with the sum of the MM charges
        """
        return sum([self.CRG_real[j - 1] for j in self.__list_MM])

    @property
    def TotCRG_HIGH(self):
        """
        property method for the total charge of the H part
        :return: a floating point number with the sum of the H charges
        """
        return sum([self.CRG_real[j - 1] for j in self.__list_HIGH])

    @property
    def TotCRG_MEDIUM(self):
        """
        property method for the total charge of the M part
        :return: a floating point number with the sum of the M charges
        """
        return sum([self.CRG_real[j - 1] for j in self.__list_MEDIUM])

    @property
    def TotCRG_LOW(self):
        """
        property method for the total charge of the L part
        :return: a floating point number with the sum of the L charges
        """
        return sum([self.CRG_real[j - 1] for j in self.__list_LOW])

    ###########################################################################################
    #                      GET METHODS FOR CHARGES OF MODEL-H
    ###########################################################################################

    @property
    def CRG_model_H_noHlink(self):
        """
        property method for modelH charges without H
        :return: a numpy array with the charges of the QM region with modelH charges
        """
        return self.CRG_model_H[:len(self.__list_QM)]

    @property
    def CRG_model_H_onlyHlink(self):
        """
        property method for modelH charges only H
        :return: a numpy array with the charges of the H atoms from modelH charges
        """
        return self.CRG_model_H[len(self.__list_QM):]

    @property
    def TotCRG_model_H(self):
        """
        property method for the total charge of modelH
        :return: a floating point number with the sum of modelH charges
        """
        return np.add.reduce(self.CRG_model_H)

    @property
    def TotCRG_model_H_noHlink(self):
        """
        property method for the total charge of the QM region for modelH charges
        :return: a floating point number with the sum of the QM charges from modelH
        """
        return np.add.reduce(self.CRG_model_H[:len(self.__list_QM)])

    @property
    def TotCRG_Hlink(self):
        """
        property method for the total charge of the H atoms charges for modelH
        :return: a floating point number with the sum of the H atom charges from modelH
        """
        return np.add.reduce(self.CRG_model_H[len(self.__list_QM):])

    ###########################################################################################

    # @staticmethod
    def embeddingCharges(self, geometry, CRG_real, CRG_model_H): # CRG_model_H is intrinsically from model-H.top
        """
        This function initialize the charges for the QM / MM calculation. Two different set of charges are given:
        CRG_real is the list of the charges for the real model part, which comes from the MM calculation setup
        CRG_model_H is the list of charges for the modelH part, which comes from the results of the QM calculation
        In this method, the two set of charges CRG_real and CRG_model_H are combined to define a unique
        set of charges, with the chosen model for redistribution at the boundary.

        In the following scheme, Q1 is an atom of the QM part, M1 is a MM atom, covalently bound to Q1.
        M1 is surrounded by M2 atoms, which are MM atoms that are chosen based on their distance from M1.

            QM           |             MM
                         |     M2
                         |    /
                Q1 --(H)-|-- M1--M2
                         |    \
                         |     M2

        This function takes the actual charge of the full system, and redefines them to avoid
        overpolarization at the boundary between QM and MM, when M1 is substituted with H.

        The behavior of the function is given by the variable whichemb, that takes values:
          0          : no redistribution of the charges
         10          : charges of M1 atoms are zeroed
         22,23,24,25 : charges of M1 are redistributed to M2 with different scaling factors
         """

        # Parameters defining the redistribution strategy
        whichemb = 25
        min_shell = 0.2
        max_shell = 1.6
        debug = True

        #print(CRG_real, CRG_model_H)
        NatomQM = len(geometry.list_QM)
        CRG_model_H_onlyHlink = np.array(CRG_model_H[NatomQM:])
        CRG_real_bak = copy.deepcopy(self.CRG_real)
        if np.sum(CRG_model_H_onlyHlink) != 0:
            self.rallyCharges(CRG_model_H) # only QM part charge is modified
            CRG_model_H[0:NatomQM] = [self.CRG_real[iatom - 1] for iatom in geometry.list_QM]
        # print(CRG_model_H, NatomQM, len(CRG_model_H))
        for i in range(NatomQM, len(CRG_model_H)):
            CRG_model_H[i] = 0.0 # virtual H always charged with 0
        self.CRG_real = copy.deepcopy(CRG_real_bak)
        # print(CRG_model_H)
        # exit()

        # define set of charges for the QM atoms only, coming from both sets of charges
        if len(geometry.list_QM) > 0:
            CRG_model = [CRG_real[iatom - 1] for iatom in geometry.list_QM]
            # print(CRG_model)
            # exit()
            CRG_model_H_noHlink = CRG_model_H[:len(geometry.list_QM)]
            Delta_CRG_model = np.add.reduce(np.array(CRG_model) - np.array(CRG_model_H_noHlink))
        else:
            CRG_model, CRG_model_H_noHlink = [], []
            Delta_CRG_model = 0.0

        # When requested, print debug info on charge distribution
        if debug:
            Log.writeLog('Using \'embeddingCharges\' to create the charge embedding used in QM calculation\n')
            if geometry.NsubH == 0:
                Log.writeLog('There are no atom-links, so the embedding charges will not be modified\n')
            else:
                Log.writeLog('There are {0} atom-links\n'.format(geometry.NsubH))
                Log.writeLog('The embedding scheme requested by is: {0}\n'.format(whichemb))
                if 0 <= whichemb < 10:
                    Log.writeLog('Embedding charges will not be modified\n')
                    Log.writeLog('Use this option with care: it will result in the loss of charge!\n')
                elif 10 <= whichemb < 20:
                    Log.writeLog('Charges of MM frontier atoms will be zeroed\n')
                    Log.writeLog('Use this option with care: it will result in the loss of charge!\n')
                elif 20 <= whichemb < 30:
                    Log.writeLog('Charges of MM frontier atoms will be distributed among the closest MM atoms\n')

        # Print debug info on model comparison
        if debug:
            Log.writeLog('\n' + '=' * 80 + '\n')
            Log.writeLog('Comparison between charges of "model" (as a part of "real")    \n')
            Log.writeLog('                          and "model" (as a part of "model-H") \n')
            Log.writeLog(' (1)    (2)        (3)               (4)              (3)-(4)  [(3)-(4)]*100   \n')
            Log.writeLog('Number Atom   model(from real)  model(from model-H)    Delta        %          \n')
            for i in range(len(geometry.list_QM)):
                idx = geometry.list_QM[i]
                diff = CRG_model[i] - CRG_model_H[i]
                perc = str(round((CRG_model[i] - CRG_model_H[i]) / ((CRG_model[i]) + 1.E-6) * 100.))
                Log.writeLog('{0:6d} {1:2s}     {2:10f}        {3:10f}         {4:10f}    {5:5s}%\n'.format(
                    idx, geometry.atomLabel[idx-1], CRG_model[i], CRG_model_H[i], diff, perc))
            Log.writeLog('                                                  _____________\n')
            Log.writeLog('                                         Delta_tot = ' + str(Delta_CRG_model) + "\n")
            Log.writeLog('=' * 80 + '\n\n')
            if whichemb == 25:
                Log.writeLog('Delta_tot will not be distributed\n\n')
            else:
                Log.writeLog('Delta_tot will be equally distributed among the ' + str(geometry.NsubH)
                                + ' atom-links\n\n')

        # no redistribution of the charges
        if 0 <= whichemb < 10 or geometry.NsubH == 0:
            CRG_real_temporary = CRG_real  # the new array if charges is identical (linked) to the old

        # zero charges for MM atoms at the boundary (only for debug purposes)
        elif 10 <= whichemb < 20:
            CRG_real_temporary = CRG_real[:]  # the new array is an hard copy of the old one

            # loop over the MM atoms involved in the atoms links, and set their charge to zero
            for iAtom in geometry.atomLink_BA[0]:
                CRG_real_temporary[iAtom - 1] = 0.0

        # redistribute charges
        elif 20 <= whichemb < 30:

            # loop over the boundary atoms, and for each boundary atom define the list of the close MM atoms
            vect_nM2 = []  # vect_nM2 contains the number of M2 atoms for each M1 atom
            matr_M2 = []  # list of M2 index lists, one list for each M1 atom
            crg_matr_M2 = []  # list of M2 charge lists, one list for each M1 atom
            abs_crg_matr_M2 = []  # list of M2 abs(charge) lists, one list for each M1 atom
            AT_matr_M2 = []  # list of M2 label lists, one list for each M1 atom
            for iAtomM1 in geometry.atomLink_BA[0]:
                # define the list of atoms that belong to the shell between min and max value
                atomM2List = [iAtomMM for iAtomMM in geometry.list_MM
                              if min_shell <= geometry.distance(iAtomM1 - 1, iAtomMM - 1) <= max_shell]
                # store the number of M2 atoms
                vect_nM2.append(len(atomM2List))
                # store the list of M2 atoms
                matr_M2.append(atomM2List)
                # store the list of charges of the M2 atoms
                crg_matr_M2.append([CRG_real[iAtomM2-1] for iAtomM2 in atomM2List])
                # store the list of abs values of the charge of the M2 atoms
                abs_crg_matr_M2.append([abs(CRG_real[iAtomM2-1]) for iAtomM2 in atomM2List])
                # store the list of atomic labels of the M2 atoms
                AT_matr_M2.append([geometry.atomLabel[iAtomM2-1] for iAtomM2 in atomM2List])

            # now compute a vector that stores the sum of M2 charges (and their abs values) for each M1 atom
            vect_crg_M2 = [np.add.reduce(atomM2Charges) for atomM2Charges in crg_matr_M2]
            abs_vect_crg_M2 = [np.add.reduce(atomM2AbsCharges) for atomM2AbsCharges in abs_crg_matr_M2]
            # and check if there is a M2 group that has zero net charge, which is incompatible with whichemb 23
            for iAtomM1, atomM2Charge in zip(geometry.atomLink_BA[0], vect_crg_M2):
                if atomM2Charge <= 1.E-6 and whichemb == 23:
                    Log.fatalError('The sum of the charges of the MM atoms surrounding the boundary atom nr.' +
                                      '{0} is {1}\nCharge is too small: cannot apply embedding scheme nr. 23!'.format(
                                          iAtomM1, atomM2Charge))

            # compute the charge to redistribute
            # for whichemb 25, redistribute only the charges of the M1 atoms
            if whichemb == 25:
                CRG_to_distribute = [CRG_real[iAtom - 1] for iAtom in geometry.atomLink_BA[0]]
            # for other whichemb, redistribute the charge of the M1 atoms + the difference between real
            # and modelH in the high level region
            else:
                CRG_to_distribute = (np.array([CRG_real[iAtom - 1] for iAtom in geometry.atomLink_BA[0]]) +
                                     (Delta_CRG_model / geometry.NsubH))

            # compute the scaling factor to use in the charge redistribution
            scaling_factor_matr_M2 = []
            if whichemb == 22:  # divide the excess charge in equal parts #USE THIS MAYBE
                for nM2Atoms in vect_nM2:
                    scaling_factor_matr_M2.append([1.0 / nM2Atoms] * nM2Atoms)
            elif whichemb == 23:  # divide the excess charge proportionally to the actual atom charge
                for M2ChargeList, M2TotalCharge in zip(crg_matr_M2, vect_crg_M2):
                    scaling_factor_matr_M2.append(np.array(M2ChargeList) / M2TotalCharge)
            elif whichemb == 24 or whichemb == 25:  # divide the excess charge proportionally to the abs(atom charge)
                for M2ChargeList, M2TotalCharge in zip(abs_crg_matr_M2, abs_vect_crg_M2):
                    scaling_factor_matr_M2.append(np.array(M2ChargeList) / M2TotalCharge)
            else:
                Log.fatalError('Invalid option for charge redistribution: ' +
                                  'scheme nr. {} does not exist\n'.format(whichemb))

            # merges the charges of all atoms in a list with the correct order
            # QM atom charges are taken from CRG_model_H_noHlink
            # MM normal atom charges are taken from CRG_real
            # charge of MM atom which are connected to QM is zeroed
            CRG_real_temporary = []
            nModelH = 0
            for iatom in range(1, geometry.atomNum + 1):
                if iatom in geometry.atomLink_BA[0]:
                    CRG_real_temporary.append(0.0)
                elif iatom in geometry.list_MM:
                    CRG_real_temporary.append(CRG_real[iatom - 1])
                elif iatom in geometry.list_QM:
                    CRG_real_temporary.append(CRG_model_H_noHlink[nModelH])
                    nModelH += 1
                else:
                    Log.fatalError('Incorrect atom index in charges.embeddingCharges!')

            # now loop over the atom links, and add the charge to redistribute to the neighbor atoms
            for i in range(geometry.NsubH):
                for k, scaling in zip(matr_M2[i], scaling_factor_matr_M2[i]):
                    CRG_real_temporary[k - 1] += CRG_to_distribute[i] * scaling

            # Print debug info on charge redistribution
            if debug and 20 <= whichemb < 30:
                # initial information and legend
                if whichemb == 25:
                    redistrModelCharge = 0.0
                else:
                    redistrModelCharge = Delta_CRG_model / geometry.NsubH
                Log.writeLog(' ka = Delta_tot / {0} = {1}\n'.format(geometry.NsubH, redistrModelCharge))
                Log.writeLog(' kb = (charge of each M1) + ka  \n')
                Log.writeLog(' kc = kb * sf \n')
                Log.writeLog(' new_crg = old_crg + kc \n')
                Log.writeLog(' diff% = (new_crg - old_crg)*100 \n\n')
                # table with charge redistribution data
                Log.writeLog('=' * 80 + '\n')
                Log.writeLog(
                    '       M1       M2     old_crg        kb       sf      kc        new_crg   diff% \n')
                Log.writeLog('-' * 80 + '\n')
                for i in range(geometry.NsubH):
                    atomM1 = geometry.atomLink_BA[0][i]
                    kb = CRG_real[atomM1-1] + redistrModelCharge
                    tmp1 = ' %2c%6u            %10f  %f \n' % (geometry.atomLabel[atomM1-1], atomM1,
                                                                  CRG_real[atomM1-1], kb)
                    Log.writeLog(tmp1)
                    for j in range(vect_nM2[i]):
                        if whichemb == 22:
                            string_sf = '1/' + str(vect_nM2[i])
                        else:
                            string_sf = str(round(scaling_factor_matr_M2[i][j], 4))
                        kc = kb * scaling_factor_matr_M2[i][j]
                        atomid = matr_M2[i][j]
                        perc = (CRG_real_temporary[atomid - 1] - crg_matr_M2[i][j])*100.
                        tmp2 = '          %2c%6u   %10f          %8s  %10f  %10f %5s%c\n' \
                               % (AT_matr_M2[i][j], matr_M2[i][j], crg_matr_M2[i][j], string_sf, kc,
                                  CRG_real_temporary[atomid - 1], str(round(perc, 1)), '%')
                        Log.writeLog(tmp2)
                    Log.writeLog('=' * 80 + '\n\n')

        else:  # in other cases, do not touch CRG_real_temporary
            CRG_real_temporary = CRG_real  # the new array if charges is identical (linked) to the old

        return CRG_real_temporary

    ###########################################################################################

    @staticmethod
    def scaledCharges(geometry, CRG_real):
        """ Scale the embedding charges according to the chosen scheme. """

        # Parameters defining the scaling function
        useScaling = 3  # 1 = use scaling function for emb charges; 3 = do not use scaling
        scalingFunction = 1  # type of scaling function: 0 = all charges set to zero; 1 = 1/r; 2,3 = sigmoidal

        if useScaling == 3:  # do not scale

            # return embedding charges without any modification
            CRG_real_scaled = CRG_real

        elif useScaling == 1:  # do scale

            # initialize the return list with the same charges of CRG_real
            CRG_real_scaled = CRG_real[:]

            # loop over all the MM atom, and scale their charges according to the minimum distance from the QM atoms
            for i in geometry.list_MM:

                # compute minimum distance of the MM atom from all the QM atoms
                r = min([geometry.distance(j, i) for j in geometry.list_QM])

                # now compute scaling function depending on the chosen scalingFunction
                sf = 1.0
                if scalingFunction == '0':  # charges are zeroed
                    sf = 0.0
                elif scalingFunction == '1':  # charges are scaled with 1 / distance
                    sf = 1 / r
                elif scalingFunction == '2':  # charges are scaled with the first sigmoidal function
                    C, B, T, M = 1.0, -1.0, 5.0, 8.0
                    sf = C / ((1 + T * math.exp(-B * (r - M))) ** (1 / T))
                elif scalingFunction == '3':  # charges are scaled with the second sigmoidal function
                    C, B, T, M = 1.0, -3.0, 10.0, 4.0
                    sf = C / ((1 + T * math.exp(-B * (r - M))) ** (1 / T))
                else:
                    Log.fatalError('Incorrect scaling funct ({0}) in charges.scaledCharges!'.format(scalingFunction))

                # store scaled charge to the return list
                CRG_real_scaled[i - 1] *= sf

        else:  # do not scale

            # return embedding charges without any modification
            CRG_real_scaled = CRG_real

        # return list with results
        return CRG_real_scaled

    ###########################################################################################

    def rallyCharges(self, CRG_model_H):


        whichemb = 23
        debug = True

        NatomQM = len(self.__list_QM)
        NsubH = len(CRG_model_H) - NatomQM

        CRG_model_H_noHlink = np.array(CRG_model_H[:NatomQM])
        CRG_model_H_onlyHlink = np.array(CRG_model_H[NatomQM:])

        if debug:
            Log.writeLog('\nUpdate of charges for the H layer with the QM calculated ones.\n')
            Log.writeLog('The requested scheme is nr. {0}: \n'.format(whichemb))
            if NsubH == 0:
                Log.writeLog('There are not atom-links, so the QM charges will not be modified.\n')
            else:
                Log.writeLog('There are {0} atom-links.\n'.format(NsubH))
                if 0 <= whichemb < 10:
                    Log.writeLog('No QM charge modification will be performed.\n')
                    Log.writeLog('Use this option with care: it will results in the loss of charge!\n')
                elif 10 <= whichemb < 20:
                    Log.writeLog('Charges from previous step will be kept.\n')
                    Log.writeLog('Use this option with care: it is intended only for debugging purposes!\n')
                elif 20 <= whichemb < 30:
                    Log.writeLog('Charges of H atoms from modelH will be distributed among the other QM atoms.\n')
                Log.writeLog('\n')

        dict_CRG_real_temporary = {}
        for j in self.__list_MM:  # MM region the charges are the same of real ones
            dict_CRG_real_temporary[j - 1] = self.CRG_real[j - 1]

        if 0 <= whichemb < 10 or NsubH == 0:  # QM charges are not modified, useful for no-atom-links calculations

            for i in range(NatomQM):  # QM region the charges are the same of QM calculation
                dict_CRG_real_temporary[self.__list_QM[i] - 1] = CRG_model_H[i]

        elif 10 <= whichemb < 20:  # the new array will be a hard copy of the old one (only for debug purposes)

            for i in self.__list_QM:  # QM region the charges are the same of initial charge list
                dict_CRG_real_temporary[i - 1] = self.CRG_real[i - 1]

        elif 20 <= whichemb < 30:

            # compute the charge to redistribute
            CRG_to_distribute = np.add.reduce(CRG_model_H_onlyHlink)

            # compute the scaling factors
            if whichemb == 22:  # partition equally
                scaling_factor = [1.0 / NatomQM] * NatomQM
            elif whichemb == 23:  # proportionally to the absolute value of the charge
                abs_CRG_model_H_noHlink = abs(CRG_model_H_noHlink)
                Tot_abs_CRG_model_H_noHlink = np.add.reduce(abs_CRG_model_H_noHlink)
                scaling_factor = [abs_CRG_model_H_noHlink[i] / Tot_abs_CRG_model_H_noHlink for i in range(NatomQM)]
            else:  # do not partition charge (excess charge will be lost!)
                scaling_factor = [0.0] * NatomQM

            # assign to QM region the modified charges
            Tot_new_CRG_model = 0.0
            for i in range(NatomQM):
                dict_CRG_real_temporary[self.__list_QM[i] - 1] = CRG_model_H[i] + scaling_factor[i] * CRG_to_distribute
                Tot_new_CRG_model += CRG_model_H[i] + scaling_factor[i] * CRG_to_distribute
            # compute difference between old and new QM total charge (it should be zero!)
            check1 = Tot_new_CRG_model - self.TotCRG_model

            if debug:
                Log.writeLog(" Old 'model' charges sum  = {0:.6f}\n".format(self.TotCRG_model))
                Log.writeLog(" QM  'model' charges sum  = {0:.6f}\n".format(sum(CRG_model_H_noHlink)))
                Log.writeLog(" CRG to distribute        = {0:.6f}\n".format(CRG_to_distribute))
                Log.writeLog(" New 'model' charges sum  = {0:.6f}\n\n".format(Tot_new_CRG_model))
                Log.writeLog('Comparison between old charges of "model"  \n')
                Log.writeLog('               and new charges of "model"  \n')
                Log.writeLog('=' * 80 + '\n')
                Log.writeLog(' (1)    (2)       (3)             (4)            (5)                (5)-(3)     \n')
                Log.writeLog('Number Atom  model(old)      model(QM)      model(new)          Delta       %   \n')
                for j in range(len(self.__list_QM)):
                    i = self.__list_QM[j]  # get atom nr in the full structure
                    deltaCharge = dict_CRG_real_temporary[i-1] - self.CRG_real[i-1]
                    Log.writeLog(
                        "{0:6d} {1:2s}    {2:10f}     {3:10f}      {4:10f}      {5:10f}    {6:5.1f}%\n".format(
                            i, self.__atomLabel[i-1], self.CRG_real[i-1], CRG_model_H[j], dict_CRG_real_temporary[i-1],
                            deltaCharge, deltaCharge*100))
                Log.writeLog('=' * 80 + '\n\n')
                Log.writeLog('  (1) Old total charge of real = {0:12f}\n'.format(self.TotCRG_model))
                Log.writeLog('  (2) New total charge of real = {0:12f}\n'.format(Tot_new_CRG_model))
                Log.writeLog('                       (1)-(2) = {0:12f} (should be zero)\n'.format(check1))

            if (abs(check1 - 0.0)) > 1.E-6:
                Log.writewarning("difference between old and new QM total charge is not zero!")
            Log.writeLog('\n')

        # convert the dictonary to a list
        crg_real_temporary = [dict_CRG_real_temporary[i] for i in range(len(self.CRG_real))]

        # now store the new charge vector as attribute of the charge instance
        self.CRG_real = np.array(crg_real_temporary)
        self.CRG_model_H = np.array(CRG_model_H)

    ###########################################################################################

    def checkConsistency(self):

        rnd = 10  # round value
        tlr1 = 1E-6  # tolerance of checks on charge electro-constance

        Log.writeLog('\n')
        check1 = self.TotCRG_real - (self.TotCRG_model + self.TotCRG_pod)
        check3 = (self.TotCRG_pod - self.TotCRG_emb)
        check4 = self.TotCRG_real - (self.TotCRG_HIGH + self.TotCRG_MEDIUM + self.TotCRG_LOW)

        Log.writeLog('  (1) TotCRG_real            = %12f  \n' % (round(self.TotCRG_real, rnd)))
        Log.writeLog('\n')
        Log.writeLog('  (2) TotCRG_model           = %12f  \n' % (round(self.TotCRG_model, rnd)))
        Log.writeLog('  (3) TotCRG_pod             = %12f  \n' % (round(self.TotCRG_pod, rnd)))
        if (abs(check1 - 0.0)) <= tlr1:
            Log.writeLog('(ck1)       {(1)-[(2)+(3)]}  = %12f ok: should be zero \n' % (round(check1, rnd)))
            Log.writeLog('                             it differs from zero less than %e \n' % tlr1)
        else:
            Log.writeLog(
                '(ck1)       {(1)-[(2)+(3)]}  = %12f warning: should be zero rather than %e\n' % (check1, check1))
            Log.writeLog('                             it differs from zero more than %e \n' % tlr1)
        Log.writeLog('\n')

        Log.writeLog('  (7) TotCRG_emb             = %12f  \n' % (round(self.TotCRG_emb, rnd)))
        if (abs(check3 - 0.0)) <= tlr1:
            Log.writeLog('(ck3)       {[(3)-(7)]}      = %12f ok: should be zero \n' % (round(check3, rnd)))
            Log.writeLog('                             it differs from zero less than %e \n' % tlr1)
        else:
            Log.writeLog('(ck3)       {[(3)-(7)]}      = %12f warning: should be zero rather than %e\n' % (check3, check3))
            Log.writeLog('                             it differs from zero more than %e \n' % tlr1)
        Log.writeLog('\n')
        Log.writeLog('  (8) TotCRG_HIGH            = %12f  \n' % (round(self.TotCRG_HIGH, rnd)))
        Log.writeLog('  (9) TotCRG_MEDIUM          = %12f  \n' % (round(self.TotCRG_MEDIUM, rnd)))
        Log.writeLog(' (10) TotCRG_LOW             = %12f  \n' % (round(self.TotCRG_LOW, rnd)))
        if (abs(check4 - 0.0)) <= tlr1:
            Log.writeLog('(ck4)  {(1)-[(8)+(9)+(10)]}  = %12f ok: should be zero \n' % (round(check4, rnd)))
            Log.writeLog('                             it differs from zero less than %e \n' % tlr1)
        else:
            Log.writeLog('(ck4)  {(1)-[(8)+(9)+(10)]}  = %12f warning: should be zero rather than %e\n' % (check4, check4))
            Log.writeLog('                             it differs from zero more than %e \n' % tlr1)
        Log.writeLog('\n')
