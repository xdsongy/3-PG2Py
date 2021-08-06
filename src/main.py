#-------------------------------------------------------------------------------
# Name:        main
# Purpose:     main function entry for 3-PG2Py
#
#
#
# Author:      Xiaodong Song
#
# Created:     08/08/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import numpy as np
import scipy.stats as ss
import sys, argparse, os
import sqlite3
from argparse import ArgumentParser, ArgumentTypeError, Action
from copy import *
import pandas as pd

from sobol_lib import *

# Global variables
import globalVariables

# Initialize soil type data into a dictionary, data retrieved from sheet "3PG_Parameters"
import getSoilTypeData 

# Initialize 3-PG model parameters retrieved from sheet "3PG_Parameters"
import getModelInitialParameterValues

# Retrieve site series name into a lis
import getSiteSeries

# Retrieve single site parameters 
import getSingleSiteParameters

# Contains the main entrance of the 3PGPy --- run3PG_2() function
import main_module, main_module_EnKF, main_module_spatial

# dual state ensemble Kalman Filter
import EnKF, Spatial

import time
from time import ctime, sleep
import threading


def main():
    print ("\n------------------------------------------------")
    print ("3-PG2Py terminated.")
        
def varianceBasedSA(sitesCollection, endage):
    """
    Sensitivity Analysis using Variance-based method, variable definition
    """

    # 35 parameters, exclude gammaF0, MaxIntcptn, LAImaxIntcptn, edgeEffect and Y
    # sa_parameters = ['alphaCx',
    #                 'MaxCond',
    #                 'tWaterMax',
    #                 'k',
    #                 'nWS',
    #                 'FR',
    #                 'pRn',
    #                 'pFS20',
    #                 'fN0',
    #                 'pRx',
    #                 'CoeffCond',
    #                 'SLA1',
    #                 'Topt',
    #                 'LAIgcx',
    #                 'Tmin',
    #                 'Tmax',
    #                 'gammaF1',
    #                 'fullCanAge',
    #                 'aWS',
    #                 'rAge']

    sa_parameters = ['alphaCx', 'MaxCond']

    if len(sa_parameters) < 2:
        print('Parameter number should not be less than 2 in SA.\n')
        exit()

    # Define the total number of the Monte Carlo simulation
    # N=2**m, m=5...14, refer Y. Tang et al.: lumped model sensitivity analysis, P800
    m = range(1,15) #5...14
    N = 2**6   # 2**14 should be enough

    # Define the variance extent
    var_extent = 0.3

    # Calculating the number of 'k' (number of selected parameters)
    k = len(sa_parameters)

    # Define empty matrices A, B to hold parameters' values
    global A, B, C, D       # revised 20210705
    A, B = np.array([]), np.array([])
    C, D = [], [] # C, D each has total k members

    # 3PG2 output values
    yA, yB, yCi, yDi= np.zeros(N), np.zeros(N), np.zeros(N), np.zeros(N)
    yC, yD = [], [] # collection of yCi and yDi

    # first-order (S) and total-effect (ST) sensitivity indices
    S, ST = [], []
    S_sobol, ST_sobol = [], [] # sobol's variance based method, refer J. Nossent, Environmental Modelling & Software xxx (2011)

    # Count number of simulation times of 3PG2
    iCount = 0

    #===================================================================================
    # Vary +/-n% either side of the selected parameters
    #===================================================================================
    def RandomizeParaValues():
        global A, B, C

        # Generate Sobol Quasi-random sequences with k dimensions
        dim_num = k # dimension

        # skips are very important in Sobol algorithm. For certain dim_num and N, the Sobol
        # will generate totally the same sequence quasi-random number. So to keep A and B be
        # different, skip2 must skip an additional number N on the base of skip1
        skip = 100

        #=======================================================================
        # Create two quasi-random vectors
        #=======================================================================
        # Note: the purpose is to generate a (N,2k) array, but sobol is (2k, N) array, so transpose sobol
        sobol= i4_sobol_generate(2*dim_num, N, skip)
        sobol = sobol.transpose()

        for i in range(0, N):
            if i==0: # 1st row
                r_A = sobol[i][0:dim_num]
                r_B = sobol[i][dim_num:2*dim_num]
            else: # the left rows
                r_A = np.vstack((r_A, sobol[i][0:dim_num]))
                r_B = np.vstack((r_B, sobol[i][dim_num:2*dim_num]))

        # Iterate through to generate one random number from a uniform distribution,
        # the range of the uniform distribution is [reference_value -/+30%*reference_value]
        for index, para in enumerate(sa_parameters):
            if para in InitialParameters_3PG.paraDict.keys() or para == 'FR': # to ensure the para name is correct
                if para != 'FR':
                    left_bound  = InitialParameters_3PG.paraDict[para] - var_extent*InitialParameters_3PG.paraDict[para]
                    right_bound = InitialParameters_3PG.paraDict[para] + var_extent*InitialParameters_3PG.paraDict[para]

                # Note: FR is not in the InitialParameters_3PG.paraDict[??]
                if para == 'FR':
                    for key in sitesCollection.keys():
                        plot = sitesCollection[key]
                        FR= plot.FR
                        left_bound  = FR - var_extent*FR
                        right_bound = FR + var_extent*FR
                        if right_bound > 1.0:
                            right_bound = 1.0

                # Generate the random para value from the uniform distribution constrained by the left_ and right_ bounds
                # a and b are arrays with only one row, but note that r_A[:, index] and r_B[:, index] are column arrays (just 1 column)
                a = left_bound + (right_bound - left_bound)*r_A[:, index]
                b = left_bound + (right_bound - left_bound)*r_B[:, index]

                # Attention: To put A, B and C into variance-based method, they have to be transposed,
                # refer A. Saltelli, 2008, "Global Sensitivity Analysis: The Primer, p164-165"
                if index==0:
                    A = np.array(a)
                    B = np.array(b)
                else:
                    A = np.vstack((A, a))
                    B = np.vstack((B, b))
            else:
                print ("para name does not exist.")

        # tanspose A and B (this is the A and B we actually need)
        A = A.transpose()
        B = B.transpose()

        # Create matrix C
        GenerateMatrixC()

    #===========================================================================
    # Generate matrix C based on A and B
    #===========================================================================
    def GenerateMatrixC():
        global A, B, C

        # Generate C, major part is A, each column is replaced by corresponding B column
        for i in range(0, k):
            Ci_Replacement = np.copy(A) # make a copy of A
            Ci_Replacement[:, i] = B[:, i]    # replace row i using row i in B
            C.append(Ci_Replacement)   # push into C list

        # Generate D, major part is B, each column is replaced by corresponding A column
        for i in range(0, k):
            Ci_Replacement = np.copy(B) # make a copy of B
            Ci_Replacement[:, i] = A[:, i]    # replace row i using row i in A
            D.append(Ci_Replacement)   # push into D list

    #===========================================================================
    # Change selected parameter (in sa_parameters[]) values indexed by 0 to N-1,
    # the values are stored in matrices A, B and C
    # @para index_N: from 0 to N-1, Monte Carlo simulation times
    # @para index_matrix: 'A', 'B' and 'C', indicating different matrices
    # @para index_Ci: list index in C[]
    #===========================================================================
    def ChangeParaValues(index_N, index_matrix, index_Ci):
        global A, B, C

        if index_matrix == 'A':
            for i in range(0, k):
                para = sa_parameters[i]
                InitialParameters_3PG.paraDict[para] = A[index_N][i]
                # Change FR value
                if para == 'FR':
                    for key in sitesCollection.keys():
                        plot = sitesCollection[key]
                        plot.FR = A[index_N][i]

            # Check if Topt > Tmax (not making sense), exchange these two variables
            if InitialParameters_3PG.paraDict['Topt'] > InitialParameters_3PG.paraDict['Tmax']:
                tmax = InitialParameters_3PG.paraDict['Tmax']
                InitialParameters_3PG.paraDict['Tmax'] = InitialParameters_3PG.paraDict['Topt']
                InitialParameters_3PG.paraDict['Topt'] = tmax

            for key in sitesCollection.keys():
                plot = sitesCollection[key]
                plot.assignTrueParameters()

        if index_matrix == 'B':
            for i in range(0, k):
                para = sa_parameters[i]
                InitialParameters_3PG.paraDict[para] = B[index_N][i]
                # Change FR value
                if para == 'FR':
                    for key in sitesCollection.keys():
                        plot = sitesCollection[key]
                        plot.FR = B[index_N][i]

            # Check if Topt > Tmax (not making sense), exchange these two variables
            if InitialParameters_3PG.paraDict['Topt'] > InitialParameters_3PG.paraDict['Tmax']:
                tmax = InitialParameters_3PG.paraDict['Tmax']
                InitialParameters_3PG.paraDict['Tmax'] = InitialParameters_3PG.paraDict['Topt']
                InitialParameters_3PG.paraDict['Topt'] = tmax

            for key in sitesCollection.keys():
                plot = sitesCollection[key]
                plot.assignTrueParameters()

        if index_matrix == 'C' or index_matrix == 'D':
            if index_matrix == 'C':
                Ci = C[index_Ci]
            elif index_matrix == 'D':
                Ci = D[index_Ci]

            for i in range(0, k):
                para = sa_parameters[i]
                InitialParameters_3PG.paraDict[para] = Ci[index_N][i]
                # Change FR value
                if para == 'FR':
                    for key in sitesCollection.keys():
                        plot = sitesCollection[key]
                        plot.FR = Ci[index_N][i]

            # Check if Topt > Tmax (not making sense), exchange these two variables
            if InitialParameters_3PG.paraDict['Topt'] > InitialParameters_3PG.paraDict['Tmax']:
                tmax = InitialParameters_3PG.paraDict['Tmax']
                InitialParameters_3PG.paraDict['Tmax'] = InitialParameters_3PG.paraDict['Topt']
                InitialParameters_3PG.paraDict['Topt'] = tmax

            for key in sitesCollection.keys():
                plot = sitesCollection[key]
                plot.assignTrueParameters()


    #===========================================================================
    # Single-thread simulation
    #===========================================================================
    t0 = time.time()

    # Generate the selected parameters' values using Monte Carlo method
    RandomizeParaValues()

    # Calculate yA and yB, the output variables we selected the follows:
    # avDBH, BasArea, LAI, StandVol, WR, WF, WS, ET, Transp, fASW
    for index_N in range(0, N):
        for key in sitesCollection.keys():
            ChangeParaValues(index_N, 'A', 0) # here index_Ci=0 is just a dummy variable
            plot = sitesCollection[key]
            plot.endAge = endage
            plot = main_module.run3PG_2_class(sitesCollection[key], False) # Run 3PG2
            yA[index_N] = plot.site.WS
            ChangeParaValues(index_N, 'B', 0) # here index_Ci=0 is just a dummy variable
            plot = main_module.run3PG_2_class(sitesCollection[key], False)
            yB[index_N] = plot.site.WS

            iCount += 2
            print (iCount)

    # Calculate yC, the output variable we select 'BasArea'
    for index_Ci in range(0, k):
        for index_N in range(0, N):
            ChangeParaValues(index_N, 'C', index_Ci)
            plot = main_module.run3PG_2_class(sitesCollection[key], False)
            yCi[index_N] = plot.site.WS

            iCount += 1
            print (iCount)

        yC.append(np.copy(yCi))

    #=======================================================================================
    # Calculate S (first-order sensitivity indices) and ST (total-effect sensitivity indices)
    #=======================================================================================
    # Saltelli method, refer Computer Physics Communications 181 (2010) 259-270
    S_total, ST_total = [], [] # Containers of S and ST, member number is len(m)
    for index in range(0, len(m)):
        S, ST = [], [] # clear empty
        sample_size = 2**m[len(m)-index-1] # m from 14 to 5

        for i in range(0, k):
            # Saltelli method
            f0_2 = (sum(yA[0:sample_size])/sample_size)**2
            vY = np.dot(yA[0:sample_size],yA[0:sample_size])/sample_size - f0_2
            Si  = np.dot(yB[0:sample_size], yC[i][0:sample_size]-yA[0:sample_size])/(vY*sample_size)   # refer eq. (2) and Table 2-(b)
            STi = sum(np.dot(yA[0:sample_size]-yC[i][0:sample_size], yA[0:sample_size]-yC[i][0:sample_size]))/(vY*2.0*sample_size)  # refer eq (4) and Table 2-(f)

            S.append(Si)
            ST.append(STi)

        S_total.append(copy(S))
        ST_total.append(copy(ST))

    # SA output file
    for key in sitesCollection.keys():
        sa_output_filename = key + ".txt"

        if sys.platform == 'win32':
            sa_output_path = "..\\Output\\" + sa_output_filename
        else:
            sa_output_path = r"../Output/" + sa_output_filename

        f = open(sa_output_path, "w")

        for index in range(0, len(m)):
            f.write("sample_exponent=%s\n" % m[len(m)-index-1])
            Si = S_total[index]
            STi = ST_total[index]
            for i in range(0, k):
                f.write("%-15s %15.6f %15.6f\n" % (sa_parameters[i], Si[i], STi[i]))

        f.close()

    #===========================================================================
    # Write file of yA, yB and yC (N row * k column)
    #===========================================================================
    for key in sitesCollection.keys():
        sa_output_filename = key + ".txt"

        if sys.platform == 'win32':
            ya = "..\\Output\\ya_" + sa_output_filename
            yb = "..\\Output\\yb_" + sa_output_filename
            yc = "..\\Output\\yc_" + sa_output_filename
        else:
            ya = r"../Output/ya_" + sa_output_filename
            yb = r"../Output/yb_" + sa_output_filename
            yc = r"../Output/yc_" + sa_output_filename

        f_ya = open(ya, "w")
        f_yb = open(yb, "w")
        f_yc = open(yc, "w")

        # Write yA, yB and yC
        for i in range(0,len(yA)):
            f_ya.write("%-15.6f\n" % yA[i])
            f_yb.write("%-15.6f\n" % yB[i])

            # write yC
            line = ''
            for j in range(0 ,k):
                line += "%-15.6f" % yC[j][i]
            f_yc.write(line + "\n")

        f_ya.close()
        f_yb.close()
        f_yc.close()

    #===========================================================================
    # Computing Si and STi. Write two files, one is Si, the other is STi,
    # the structure is(N row * k column), each column is one parameter
    #===========================================================================
    for key in sitesCollection.keys():
        sa_output_filename = key + ".txt"

        if sys.platform == 'win32':
            si_output_path = "..\\Output\\si_" + sa_output_filename
            sti_output_path = "..\\Output\\sti_" + sa_output_filename
        else:
            si_output_path = r"../Output/si_" + sa_output_filename
            sti_output_path = r"../Output/sti_" + sa_output_filename

        f_si  = open(si_output_path, "w")
        f_sti = open(sti_output_path, "w")

        # write file header using parameter names
        for para in sa_parameters:
            f_si.write("%-15s" % para)
            f_sti.write("%-15s" % para)

        f_si.write("\n")
        f_sti.write("\n")

        for sample_size in range(2, N):
            line_si, line_sti = '', '' # empty lines
            for i in range(0, k):
                # Saltelli method
                f0_2 = (sum(yA[0:sample_size])/sample_size)**2
                vY = np.dot(yA[0:sample_size],yA[0:sample_size])/sample_size - f0_2
                Si  = np.dot(yB[0:sample_size], yC[i][0:sample_size]-yA[0:sample_size])/(vY*sample_size)   # refer eq. (2) and Table 2-(b)
                STi = sum(np.dot(yA[0:sample_size]-yC[i][0:sample_size], yA[0:sample_size]-yC[i][0:sample_size]))/(vY*2.0*sample_size)  # refer eq (4) and Table 2-(f)
                line_si  += "%-15.6f" % Si
                line_sti += "%-15.6f" % STi

            # write S and ST into f_si and f_sti, respectively
            f_si.write(line_si + "\n")
            f_sti.write(line_sti + "\n")

        f_si.close()
        f_sti.close()


    print ("\nSensitivity analysis execution time is: %4.3f s" % (time.time()-t0))

    main()


def FAST(endage):
    # ===========================================================================
    # FAST method
    # ===========================================================================

    # input parameters
    # input_parameters = ['SLA1',
    #                     'nWS',
    #                     'alphaCx',
    #                     'MaxCond',
    #                     'tWaterMax',
    #                     'pFS20',
    #                     'gammaF1',
    #                     'k',
    #                     'pRn',
    #                     'pRx',
    #                     'CoeffCond',
    #                     'fN0',
    #                     'Topt',
    #                     'aWS',
    #                     'Tmax',
    #                     'Tmin',
    #                     'LAIgcx',
    #                     'fullCanAge',
    #                     'cVPD0',
    #                     'rAge']

    input_parameters = ['SLA1',
                        'nWS']

    # using the following format is you're using command line mode to run FAST
    output_filename = "fast_%s" % endage + ".txt"
    # output_filename = "fast_example.txt"

    path = ""
    if sys.platform == 'win32':
        path = "..\\Output\\" + output_filename
    else:
        path = r"../Output/" + output_filename

    output = open(path, "w")

    for para in input_parameters:
        output.write("%-15s" % para)
    output.write("\n")

    print("Note: 3-PG2Py is a Python version of 3-PGxl(VBA based).\n")
    print("=================Initialization=================")

    count = 0
    for times in range(0, 30):
        # Generate global variables
        GlobalVar = globalVariables.globalVariablesClass()

        # set the name of the Excel sheet, in which the plot name for FAST SA is given
        GlobalVar.setSiteSeries("Site")

        # Seed excel file & sheet names into Class readSoilTypeInfoClass(), then fill the soil type info
        SoilType = getSoilTypeData.readSoilTypeInfoClass(GlobalVar.excel_file,
                                                         GlobalVar.parameters_3PG)

        # Seed 3PG initial parameters using Class readInitialParametersClass(), then fill the parameter dictionary
        InitialParameters_3PG = getModelInitialParameterValues.readInitialParametersClass(GlobalVar.excel_file,
                                                                                          GlobalVar.parameters_3PG)

        # Get site series name list
        SiteSeriesList = getSiteSeries.getSiteSeriesClass(GlobalVar.excel_file,
                                                          GlobalVar.siteSeries)

        # Get each single site's parameters
        SingleSiteParameters = getSingleSiteParameters.assembleSingleSite(GlobalVar.excel_file,
                                                                          SiteSeriesList.siteSeries,
                                                                          SoilType.soilDict,
                                                                          InitialParameters_3PG.paraDict,
                                                                          GlobalVar.sitesParameterCollection)

        print("\nTotal number of sites: %d\n" % len(SiteSeriesList.siteSeries))

        sitesCollection = GlobalVar.sitesParameterCollection

        # Change age parameter
        for key in sitesCollection.keys():
            plot = sitesCollection[key]
            plot.endAge = endage
            # print(plot.endAge)

        # Check whether the output file "Py3PG2Output.sqlite" exists. If exist, clear it.
        try:
            # Create output file. If it exists already, clear it's content
            file = open(GlobalVar.outputFilename, 'w')
            file.close()
        except IOError:
            print(e)

        print("===================Simulation===================")

        ##################################################

        # global variables
        global X

        # Total model run times, must be odd number
        N = 201    # 2001
        count += 1

        # frequency and harmonics
        wi = 1
        M = 6

        # input parameter number
        kk = len(input_parameters)

        # create total N points over (-pi, pi)
        s = []
        for j in range(1, N + 1):
            s.append(-np.pi + np.pi / N + (2 * np.pi / N) * (j - 1))  # refer Xu 2008 RESS eq.(5)

        # permute s
        S = []
        for i in range(0, kk):
            S.append(list(np.random.permutation(s)))

        # for the k input parameters, permute the N points (s), respectively
        X = deepcopy(S)  # container of all the Xi(permutation of s)
        for i in range(0, kk):
            for j in range(0, N):
                # using formula xi=1/2 + (1/pi)*arcsin(sin(wi*s)), refer Saltelli 1999 eq. (19)
                X[i][j] = 0.5 + (1.0 / np.pi) * np.arcsin(np.sin(wi * X[i][j]))

        # change Xi from 0-1 to real input parameter values, the range is varied 30% at either side.
        var_extent = 0.3

        # ===================================================================================
        # Vary +/-n% either side of the input parameters
        # ===================================================================================
        def RandomizeParaValues():
            global X

            # Iterate through to generate one random number from a uniform distribution,
            # the range of the uniform distribution is [reference_value -/+30%*reference_value]
            for index, para in enumerate(input_parameters):
                if para in InitialParameters_3PG.paraDict.keys() or para == 'FR':  # to ensure the para name is correct
                    if para != 'FR':
                        left_bound = InitialParameters_3PG.paraDict[para] - var_extent * InitialParameters_3PG.paraDict[
                            para]
                        right_bound = InitialParameters_3PG.paraDict[para] + var_extent * \
                                      InitialParameters_3PG.paraDict[para]

                    # Note: FR is not in the InitialParameters_3PG.paraDict[??]
                    if para == 'FR':
                        for key in sitesCollection.keys():
                            plot = sitesCollection[key]
                            FR = plot.FR
                            left_bound = FR - var_extent * FR
                            right_bound = FR + var_extent * FR
                            if right_bound > 1.0:
                                right_bound = 1.0

                    # Generate the random para value from the uniform distribution constrained by the left_ and right_ bounds
                    # a and b are arrays with only one row, but note that r_A[:, index] and r_B[:, index] are column arrays (just 1 column)
                    for i in range(0, N):
                        X[index][i] = left_bound + (right_bound - left_bound) * X[index][i]

                else:
                    print("para name does not exist.")

        # ===========================================================================
        # Change selected parameter (in input_parameters[]) values indexed by 0 to k-1
        # @para index_N: from 0 to N-1
        # ===========================================================================
        def ChangeParaValues(index_N, flag=0):
            global X

            if flag == 1:
                for key in sitesCollection.keys():
                    plot = sitesCollection[key]
            #                    plot.FR = 0.2

            for i in range(0, kk):
                para = input_parameters[i]
                InitialParameters_3PG.paraDict[para] = X[i][index_N]
                # Change FR value
                if para == 'FR':
                    for key in sitesCollection.keys():
                        plot = sitesCollection[key]
                        plot.FR = X[i][index_N]

                        # Check if Topt > Tmax (not making sense), exchange these two variables
            if InitialParameters_3PG.paraDict['Topt'] > InitialParameters_3PG.paraDict['Tmax']:
                tmax = InitialParameters_3PG.paraDict['Tmax']
                InitialParameters_3PG.paraDict['Tmax'] = InitialParameters_3PG.paraDict['Topt']
                InitialParameters_3PG.paraDict['Topt'] = tmax

            for key in sitesCollection.keys():
                plot = sitesCollection[key]
                plot.assignTrueParameters()

                # Generate the input parameters' values

        RandomizeParaValues()

        # model output
        Y = []
        Y_FR = []
        # begin the simulations
        iCount = 0  # Count number of simulation times of 3PG2
        for index_N in range(0, N):
            for key in sitesCollection.keys():
                ChangeParaValues(index_N)
                plot = main_module.run3PG_2_class(sitesCollection[key], False)  # Run 3PG2
                # print("endAge", plot.site.endAge)
                Y.append(plot.site.LAI)
                iCount += 1
                print(iCount, count)

        # =================================================================================
        # Reord the output Y so that Xi[0,...,N-1], i=0,...,k-1 are in increasing order
        # =================================================================================
        Y_reorder = []
        S_reorder = []
        for i in range(0, kk):
            rank = ss.rankdata(
                X[i]) - 1  # ranked in increasing order, rank(1,...,N), so must minus 1
            Y_new = [0] * N  # reordered output corresponding to X[i]
            s_new = [0] * N
            for j in range(0, N):
                Y_new[int(rank[j])] = Y[j]  # minus 1 to guarantee the index between 0~N-1
                s_new[int(rank[j])] = S[i][j]
            Y_reorder.append(Y_new)
            S_reorder.append(s_new)

        # ===========================================================================
        # Calculate Fourier spectrum
        # ===========================================================================
        Vi_collection = []
        for para_index in range(0, kk):
            Vi = 0.0  # partial variance, refer Xu 2008 RESS, eq.(9)
            for w in range(1, M + 1):  # or w
                A, B = 0.0, 0.0
                for sample_index in range(0, N):
                    A += Y_reorder[para_index][sample_index] * np.cos(S_reorder[para_index][sample_index] * w)
                    B += Y_reorder[para_index][sample_index] * np.sin(S_reorder[para_index][sample_index] * w)
                Ak = A * 2 / N
                Bk = B * 2 / N
                Vi += 0.5 * (Ak ** 2 + Bk ** 2)
            Vi_collection.append(Vi)

        print(Vi_collection)

        # ===========================================================================
        # Calculate total variance
        # ===========================================================================
        V_collection = []  # total variance
        for para_index in range(0, kk):
            V = 0.0
            for w in range(1, int((N-1)/2) + 1):
                A, B = 0.0, 0.0
                for sample_index in range(0, N):
                    A += Y[sample_index] * np.cos(s[sample_index] * w)
                    B += Y[sample_index] * np.sin(s[sample_index] * w)
                Ak = A * 2 / N
                Bk = B * 2 / N
                V += 0.5 * (Ak ** 2 + Bk ** 2)
            V_collection.append(V)

        Vi = []
        for i in range(0, kk):
            Vi.append(Vi_collection[i] / V_collection[i])

        print(Vi)

        for ele in Vi:
            output.write("%-15.6f" % ele)
        output.write("\n")

    time.ctime()
    output.close()


class spatialSimulationDemo():

    def __init__(self, sitesCollection):
        # pixel size, for example 10
        self.pixelNumber = 10

        plot = sitesCollection[list(sitesCollection.keys())[0]]
        plot.endAge = args.endage
        ensembleContainer = Spatial.spatialSimulation(plot, self.pixelNumber)

        for ensembleID in range(self.pixelNumber):
            ensembleContainer.subscribe(main_module_spatial.run3PG2(plot, ensembleID))

        # change initial parameter values
        parameters_list = []    # add any parameter needs to be changed among different pixels
        ensembleContainer.changeInitialParameterValues(parameters_list)

        standAge = float(plot.InitialAge[0]) + float(plot.InitialAge[1]) / 10

        timer = 0
        while True:
            print('Simulation time = ', timer)

            # Main body of model simulation
            ensembleContainer.notifySubscribers()

            timer += 1
            standAge += 1 / 12.0
            if standAge - plot.endAge > 0 and abs(standAge - plot.endAge) > 0.01:
                break

class dualStateEnKF():

    def __init__(self, sitesCollection, parameterSpaceSamplingIndex):

        # ensemble size
        self.ensembleSize = 30

        # To calculate the Kalman gains in A(3) and B(3), the last term is about the variance of observation/prediction.
        # We will use std_scaleFactor multiple the average of predictions at time t to represent the std of prediction.
        self.std_scaleFactor = 0.1

        plot = sitesCollection[list(sitesCollection.keys())[0]]
        plot.endAge = args.endage
        ensembleContainer_A = EnKF.dualStateEnKF(plot, self.ensembleSize)

        # open sqlite connection
        DB = sqlite3.connect(GlobalVar.outputFilename, check_same_thread=False)

        for ensembleID in range(self.ensembleSize):
            ensembleContainer_A.subscribe(main_module_EnKF.run3PG_2_class(plot, ensembleID))

        # Change values of the parameters and forcing factors and collect state variables' values at each step.
        # Assimilate observed variables' values at specific stand ages, note that the observedStateVariables[]
        # is only a subset of stateVariables[]
        forcingFactors_enkf = ['Tx', 'Tn', 'VPD', 'Rain', 'SolarRad']
        parameters_enkf = ['MaxCond', 'nWS', 'alphaCx']
        parameter_bounds = {'MaxCond':[0.01, 0.05], 'nWS':[1, 4], 'alphaCx':[0.01, 0.2]}
        stateVariables_enkf = ['WF', 'WS', 'WR', 'StandVol', 'forestSW']  # 'StemNo',
        observedStateVariables_enkf = ['avDBH']      # e.g., 'avDBH', 'LAI', 'StandVol'

        # Ref: getSingleSiteParameters.py Line 240-243
        standAge = float(plot.InitialAge[0]) + float(plot.InitialAge[1]) / 10

        # Collect pseudo observed values corresponding to the output variables listed in observedStateVariables_enkf[]
        # using default site parameter values. Note that if true observation data is available, i.e., there should be many
        # missing values compared with this pseudo scenario, the deliverPrerunning() function should be redesigned.
        self.pseudoObservedVariableValues = ensembleContainer_A.deliverPrerunning(DB, standAge, plot.endAge, observedStateVariables_enkf)

        # change initial parameter values
        ensembleContainer_A.changeInitialParameterValues(parameters_enkf, parameter_bounds)

        # to store all the values of parameters_enkf, stateVariables_enkf and observedStateVariables_enkf by line
        cols = parameters_enkf + stateVariables_enkf + observedStateVariables_enkf
        cols.append('Observed')
        df = pd.DataFrame(columns = cols)

        timer = 0
        while True:
            print('Simulation time = ', timer)

            # create perturbed observed variable's values
            if timer in self.pseudoObservedVariableValues['time']:
                observedValues = self.createObservedVariableValues(self.pseudoObservedVariableValues,
                                                                   observedStateVariables_enkf,
                                                                   self.ensembleSize,
                                                                   timer)

            # Get the month index in order to get the disturbed forcing data for each ensemble member
            monthIndex = divmod(plot.InitialMonth + timer, 12)[1]  # zero based month index, i.e., from 0 to 11

            # Deepcopy ensembleContainer_A. Keeping ut (forcing data) unchanged as in Step A, the deepcopied
            # ensembleContainer_A will use the updated parameter values (and ut) to get the updated state and observed
            # variable values In Step B(1-2). Note that in A(1) and B(1), the two f(.) (or state of the plot instance)
            # are completely the same (state variable values are the same), except values of the ensemble of thetas.
            ensembleContainer_B = deepcopy(ensembleContainer_A)

            # Step C: Kernel smoothing of parameter samples at each step
            thetas_forecast = ensembleContainer_A.deliverChangeParameterValues(timer, parameters_enkf, parameter_bounds)

            # 1. Main body of model simulation, including climate forcing data perturbation.
            # 2. Obtain model state trajectories (state forecast) by propagating the replicates of forcing data and
            # parameters through the forward model
            forcingDataPerturbations = ensembleContainer_A.fixForcingData(monthIndex, forcingFactors_enkf)
            ensembleContainer_A.deliverModelRunning(monthIndex, forcingFactors_enkf, forcingDataPerturbations)

            # Step A(1-2): collect forecasted state variables (stateVariables_enkf[]) and
            # observations (observedStateVariables_enkf[]).
            states_forecast_A, predictions_forecast_A = ensembleContainer_A.deliverCollectState_ObservationValues(stateVariables_enkf,
                                                                                                             observedStateVariables_enkf)
            # Step A(3): Update the ensemble of parameters according to the standard Kalman equation.
            self.updatedParameters = self.updateParameters(timer, predictions_forecast_A, thetas_forecast, observedValues)

            # Step B: re-run the simulation at time t using ensembleContainer_B and self.updatedParameters. Note that
            # the forcing data is unchanged.
            # Step B(1-2)
            ensembleContainer_B.changeParameterValues(self.updatedParameters)
            ensembleContainer_B.deliverModelRunning(monthIndex, forcingFactors_enkf, forcingDataPerturbations)
            states_forecast_B, predictions_forecast_B = ensembleContainer_B.deliverCollectState_ObservationValues(stateVariables_enkf,
                                                                                                                  observedStateVariables_enkf)

            # Step B(3): Update the ensemble of model states according to the standard Kalman equation.
            # In this step, the state variable(s) are updated. Then at the next time t+1, the state(s) and thetas in f(.)
            # shall be updated for the next simulation (the thetas should be first smoothed).
            self.updatedStates = self.updateStateVariables(timer, predictions_forecast_B, states_forecast_B, observedValues)
            ensembleContainer_B.changeStateValues(self.updatedStates)

            if timer in self.pseudoObservedVariableValues['time']:
                p = self.pseudoObservedVariableValues['time'].index(timer)
                print('observed = ', self.pseudoObservedVariableValues['avDBH'][p])
                print('predicted = ', np.mean(predictions_forecast_B['avDBH']))

            # Replace the states listed in stateVariables_enkf{} in ensembleContainer_B by self.updatedStates,
            # then replace ensembleContainer_A by ensembleContainer_B
            # ensembleContainer_B.changeStateValues(self.updatedStates)
            ensembleContainer_A = deepcopy(ensembleContainer_B)

            # write self.updatedParameters, self.updatedStates and predictions_forecast_B to pandas df
            self.writeDataFrame(df, timer, self.updatedParameters, self.updatedStates, predictions_forecast_B, self.pseudoObservedVariableValues)

            timer += 1
            standAge += 1 / 12.0
            if standAge - plot.endAge > 0 and abs(standAge - plot.endAge) > 0.01:
                # check if the sub-directory "EnKF" under "Output" exists; if not, then create it.
                if not os.path.isdir('..\\Output\\EnKF'):
                    os.mkdir('..\\Output\\EnKF')
                # write pandas df to txt file
                fname = '..\\Output\\EnKF\\' + str(parameterSpaceSamplingIndex) + '.txt'
                df.to_csv(fname)

                break

        # close sqlite connection
        DB.commit()
        DB.cursor().close()
        DB.close()


    def createObservedVariableValues(self, pseudoObservedVariableValues, observedStateVariables_enkf, ensembleSize, timer):
        observedValues = []
        for i in range(ensembleSize):
            index = pseudoObservedVariableValues['time'].index(timer)
            data = pseudoObservedVariableValues[observedStateVariables_enkf[0]][index]
            val = data + np.random.normal(0, 0.3 * data)
            while val / data < 0.5 or val / data > 1.5 and val > 0:
                val = data + np.random.normal(0, 0.3 * data)

            observedValues.append(val)

        return observedValues


    def writeDataFrame(self, df, timer, updatedParameters, updatedStates, predictions_forecast_B, pseudoObservedVariableValues):
        temp = []

        for key in updatedParameters.keys():
            temp.append(np.mean(updatedParameters[key]))

        for key in updatedStates.keys():
            temp.append(np.mean(updatedStates[key]))

        for key in predictions_forecast_B.keys():
            temp.append(np.mean(predictions_forecast_B[key]))
            if timer in self.pseudoObservedVariableValues['time']:
                index = self.pseudoObservedVariableValues['time'].index(timer)
                temp.append(self.pseudoObservedVariableValues[key][index])
            else:
                temp.append(0)

        df.loc[timer] = temp


    def cov(self, series1, series2):
        """
        Calculate the covariance
        """
        # mean values of series1 and series2
        mean_series1, mean_series2 = np.mean(series1), np.mean(series2)

        sum = 0.0
        for i in range(len(series1)):
            sum += (series1[i] - mean_series1)*(series2[i] - mean_series2)

        return sum/(len(series1)-1)


    def updateStateVariables(self, timer, predictions_forecast, states_forecast, observedValues):
        """
        Step B(3): Update the ensemble of state variable(s) according to the standard Kalman equation.
        :return: Updated ensemble of states.
        """
        #------------------ Ref B(3), first we calculate the Kalman gain ------------------#
        # Step 1: calculate the cross covariance of state ensemble and prediction ensemble
        # forecast error covariance
        cov_predictions = {}
        for key in predictions_forecast.keys():
            # cov_predictions[key] = self.cov(predictions_forecast[key], predictions_forecast[key])
            cov_predictions[key] = np.var(predictions_forecast[key])

        # state
        cov_state_prediction = {}
        for key_state in states_forecast.keys():
            # In current stage, the predictions_forecast{} only contains one member. For compatibility in the future,
            # we still use the following loop in case of there are more than one prediction variables.
            for key_prediction in predictions_forecast.keys():
                cov_state_prediction[key_state] = self.cov(predictions_forecast[key_prediction], states_forecast[key_state])

        # prediction variance
        # In current stage, the varPrediction{} only contains one member.
        varPrediction = {}
        for key in predictions_forecast.keys():
            # varPrediction[key] = self.std_scaleFactor*np.mean(predictions_forecast[key])
            varPrediction[key] = np.var(predictions_forecast[key])

        # Step 2: calculate K_states in B(3)
        K_states = {}
        for key_state in states_forecast.keys():
            for key_prediction in cov_predictions.keys():
                # K_states[key_state] = cov_state_prediction[key_state]*(1./(cov_predictions[key_prediction] + varPrediction[key_prediction]))
                K_states[key_state] = cov_state_prediction[key_state] * (
                            1. / (cov_predictions[key_prediction] + np.var(observedValues)))

        # Step 3: calculate the updated states. If at time t there is one observation, then update the parameters/thetas;
        # else keeping the updated parameters/thetas same as the thetas_forecast.
        if timer in self.pseudoObservedVariableValues['time']:  # it means that there is one observation at current time
            for key_state in states_forecast.keys():
                for key_prediction in predictions_forecast.keys():
                    for index in range(self.ensembleSize):
                        # observedValue = self.pseudoObservedVariableValues[key_prediction][timer] + np.random.normal(0,varPrediction[key_prediction]**0.5)
                        states_forecast[key_state][index] += K_states[key_state]*(observedValues[index] - predictions_forecast[key_prediction][index])

        return states_forecast


    def updateParameters(self, timer, predictions_forecast, thetas_forecast, observedValues):
        """
        Step A(3): Update the ensemble of parameters according to the standard Kalman equation.
        :return: Updated ensemble of parameters.
        """
        #------------------ Ref A(3), first we calculate the Kalman gain ------------------#
        # Step 1: calculate the cross covariance of parameter ensemble and prediction ensemble
        # forecast error covariance
        cov_predictions = {}
        for key in predictions_forecast.keys():
            # cov_predictions[key] = self.cov(predictions_forecast[key], predictions_forecast[key])
            cov_predictions[key] = np.var(predictions_forecast[key])

        # parameter
        cov_theta_prediction = {}
        for key_theta in thetas_forecast.keys():
            # In current stage, the predictions_forecast{} only contains one member. For compatibility in the future,
            # we still use the following loop in case of there are more than one prediction variables.
            for key_prediction in predictions_forecast.keys():
                cov_theta_prediction[key_theta] = self.cov(predictions_forecast[key_prediction], thetas_forecast[key_theta])

        # Step 2: calculate K_thetas in A(3)
        K_thetas = {}
        for key_theta in thetas_forecast.keys():
            for key_prediction in cov_predictions.keys():
                # K_thetas[key_theta] = cov_theta_prediction[key_theta]*(1./(cov_predictions[key_prediction] + varPrediction[key_prediction]))
                K_thetas[key_theta] = cov_theta_prediction[key_theta] * (
                            1. / (cov_predictions[key_prediction] + np.var(observedValues)))

        # Step 3: calculate the updated thetas. If at time t there is one observation, then update the parameters/thetas;
        # else keeping the updated parameters/thetas same as the thetas_forecast.
        if timer in self.pseudoObservedVariableValues['time']:  # it means that there is one observation at current time
            for key_theta in thetas_forecast.keys():
                for key_prediction in predictions_forecast.keys():
                    for index in range(self.ensembleSize):
                        thetas_forecast[key_theta][index] += K_thetas[key_theta]*(observedValues[index] - predictions_forecast[key_prediction][index])

        return thetas_forecast


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="===================3-PG2Py Simulation===================")
    parser.add_argument('-m', '--mode', dest='mode', type=str, required=True, choices=['r', 'spatial', 'VB', 'FAST', 'EnKF'], \
                                                                    help='execution mode: r (3-PG2Py simulation), \
                                                                                        spatial (spatial simulation), \
                                                                                        VB (Variance-based SA), \
                                                                                        FAST (FAST SA), \
                                                                                        EnKF (ensemble Kalman Filter)')
    parser.add_argument('-a', '--endage', dest='endage', type=int, required=True, help='endAge of the simulation and \
                                                                                       should be less than 200 years')
    print(parser.description)
    args = parser.parse_args()

    if args.endage > 200:
        raise argparse.ArgumentTypeError("%s is an invalid positive int value for the endage" % args.endage)

    # Generate global variables
    GlobalVar = globalVariables.globalVariablesClass()

    if args.mode == 'r':    # run 3-PG2Py directly
        GlobalVar.setSiteSeries("SiteSeries")
    else:     # For VB, FAST, and EnKF, we only use the site listed in the "Site" Sheet in 3PGxl_Parameters.xls
        GlobalVar.setSiteSeries("Site")
    
    # Seed excel file & sheet names into Class readSoilTypeInfoClass(), then fill the soil type info 
    SoilType = getSoilTypeData.readSoilTypeInfoClass(GlobalVar.excel_file, 
                                                     GlobalVar.parameters_3PG)
   
    # Seed 3PG initial parameters using Class readInitialParametersClass(), then fill the parameter dictionary 
    InitialParameters_3PG = getModelInitialParameterValues.readInitialParametersClass(GlobalVar.excel_file, 
                                                                                      GlobalVar.parameters_3PG)
    
    # Get site series name list
    SiteSeriesList = getSiteSeries.getSiteSeriesClass(GlobalVar.excel_file, 
                                                      GlobalVar.siteSeries)
    
    print ("Note: 3-PG2Py is a Python version of 3-PGxl(VBA based).\n")
    print ("=================Initialization=================")
    
    # Get each single site's parameters 
    SingleSiteParameters = getSingleSiteParameters.assembleSingleSite(GlobalVar.excel_file, 
                                                                      SiteSeriesList.siteSeries, 
                                                                      SoilType.soilDict,
                                                                      InitialParameters_3PG.paraDict,
                                                                      GlobalVar.sitesParameterCollection)

    sitesCollection = GlobalVar.sitesParameterCollection

    # Check whether the output file "Py3PG2Output.sqlite" exists. If exist, clear it.
    try:
        # Create output file. If it exists already, clear its content
        file = open(GlobalVar.outputFilename, 'w')
        file.close()
    except IOError:
        print('File open error.')

    if args.mode == 'r':    # run 3-PG2Py directly
        t0 = time.time()
        for key in sitesCollection.keys():
            plot = sitesCollection[key]
            plot.endAge = args.endage
            main_module.run3PG_2_class(plot)
            print('{} ... is finished.\n'.format(key))
        print("\nTotal execution time is: %4.3f s" % (time.time() - t0))
    elif args.mode == 'spatial':  # spatial simulation
        print('spatial simulation mode')
        spatialSimulationDemo(sitesCollection)
    elif args.mode == 'VB':    # variance-based SA
        GlobalVar.setSiteSeries("Site")
        varianceBasedSA(sitesCollection, args.endage)
    elif args.mode == 'FAST':    # FAST SA
        FAST(args.endage)
    elif args.mode == 'EnKF':    # dual EnKF. For example: python main.py -m EnKF -a 50
        # Total number of ensemble running, and each ensemble running is initiated with a different set of
        # parameter values in order to quantify the uncertainties of the parameters and state variables.
        # Ref: Line 1 on P144 in Moradkhani 2005.
        parameterSpaceSamplingNum = 500   # 500?
        for i in range(parameterSpaceSamplingNum):
            print('Parameter space sampling number: ', i+1, '/', str(parameterSpaceSamplingNum))
            EnKF_Series = dualStateEnKF(sitesCollection, i)
