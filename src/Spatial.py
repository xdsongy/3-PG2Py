# -------------------------------------------------------------------------------
# Name:        EnKF
# Purpose:     Spaital simulation interface
# Author:      Xiaodong Song
#
# Created:     1/08/2021
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2021
# -------------------------------------------------------------------------------
# !/usr/bin/env python

import numpy as np
import sys

import getSingleSiteParameters

# Contains the main entrance of the 3PGPy --- run3PG_2() function
import main_module_EnKF


class spatialSimulation():

    def __init__(self, plot, pixelNumber):

        # pixel number
        self.n = pixelNumber

        # end age of the simulation
        self.T = plot.endAge

        # self.site = plot

        # container for the ensemble
        self.plotInstanceObservers = []


    # multiple deliver* type functions need to be defined in order to deliver different kinds of callings during the
    # simulation to each ensemble member. Accordingly, the Observer class (main_module_spatial) needs to 'subscribe' the
    # following given deliver* functions.
    def notifySubscribers(self):
        for index, siteInstance in enumerate(self.plotInstanceObservers):
            print('pixel ID = ', index)
            siteInstance.run3PG2_Spatial()


    def changeInitialParameterValues(self, parameters_list):
        print('changeInitialParameterValues()')
        for index, siteInstance in enumerate(self.plotInstanceObservers):
            for key in parameters_list:
                str = 'siteInstance.site.' + key + '=' + 'one specific value'
                exec(str)


    def subscribe(self, observer):
        self.plotInstanceObservers.append(observer)
