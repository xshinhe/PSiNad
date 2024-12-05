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

def ask(question, qtype:str, constrain=False, accepted=[]):

    
    ## STRING
    if qtype == "string":

        answered = False
        while answered == False:
            answer = input(question)
            if constrain:
                if answer in accepted:
                    answered = True
                else:
                    print("ERROR: Please enter one of the following: {}\n".format([allowed for allowed in accepted]))
                    continue
            else:
                answered = True


    ## INTEGER

    elif qtype == "int":

        answered = False
        while answered == False:
            try:
                answer = int(input(question))
                answered = True
            except ValueError:
                print("ERROR: Please enter an integer\n")
                continue

    ## FLOAT

    elif qtype == "float":

        answered = False
        while answered == False:
            try:
                answer = float(input(question))
                answered = True
            except ValueError:
                print("ERROR: Please enter a float\n")
                continue

    ## BOOLEAN

    elif qtype == "bool":

        answered = False
        yes = ["yes", "YES", "Yes"]
        no = ["no", "NO", "No"]
        while answered == False:
            given_answer = input(question)
            if given_answer in yes:
                answer = True
                answered = True
            elif given_answer in no:
                answer = False
                answered = True
            else:
                print("\nERROR: Please enter yes or no\n")
                continue

    return answer

###################################################################

