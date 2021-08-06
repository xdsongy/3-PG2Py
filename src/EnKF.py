# -------------------------------------------------------------------------------
# Name:        EnKF
# Purpose:     dual state ensemble Kalman Filter
# Ref:         Moradkhani H, et al. Dual state-parameter estimation of hydrological models using ensemble Kalman filter.
#               Advances in Water Resources, 2005, 28(2):135-147.
# Author:      Xiaodong Song
#
# Created:     12/07/2021
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2021
# -------------------------------------------------------------------------------
# !/usr/bin/env python

import numpy as np
import sys

import getSingleSiteParameters

# Contains the main entrance of the 3PGPy --- run3PG_2() function
import main_module_EnKF


class dualStateEnKF():

    def __init__(self, plot, ensembleSize):

        # ensemble size of the dual state Kalman filter
        self.n = ensembleSize

        # end age of the assimilation
        self.T = plot.endAge

        self.site = plot

        # container for the ensemble
        self.siteInstanceObservers = []

        # get the variances of each selected climate factor
        self.forcingDataVariances = {}
        self.getForcingDataVariances()

        # self.runEnKF()

    # multiple deliver* type functions need to be defined in order to deliver different kinds of callings during the
    # simulation to each ensemble member. Accordingly, the Observer class (main_module_EnKF) needs to 'subscribe' the
    # following given deliver* functions.
    def deliverModelRunning(self, monthIndex, forcingFactors_enkf, forcingDataPerturbations):
        for siteInstance in self.siteInstanceObservers:
            siteInstance.run3PG2_EnsembleVersion(monthIndex, forcingFactors_enkf, forcingDataPerturbations)


    def fixForcingData(self, monthIndex, forcingFactors_enkf):
        """
        Fix the forcing data perturbations at each time, because in A(1) and B(1) two identical f(.) will be run so
        it's necessary to make sure the forcing data in A(1) and B(1) are the same.
        """
        forcingDataVariances = {}
        for variable in forcingFactors_enkf:
            # forcingDataVariances[variable] = np.random.normal(0, self.forcingDataVariances[variable][monthIndex])
            forcingDataVariances[variable] = self.forcingDataVariances[variable][monthIndex]

        return forcingDataVariances



    def deliverCollectState_ObservationValues(self, stateVariables_enkf, observedStateVariables_enkf):
        """
        Step A(1-2): collect forecasted state variables (stateVariables_enkf[]) and
        observations (observedStateVariables_enkf[]).
        """

        # collect state variable values
        stateVariableValues = {}

        for stateVariable in stateVariables_enkf:
            stateVariableValues[stateVariable] = []

        for siteInstance in self.siteInstanceObservers:
            for stateVariable in stateVariables_enkf:
                stateVariableValues[stateVariable].append(eval('siteInstance.site.' + stateVariable))

        # collect observation variable values
        observationVariableValues = {}

        for observationVariable in observedStateVariables_enkf:
            observationVariableValues[observationVariable] = []

        for siteInstance in self.siteInstanceObservers:
            for observationVariable in observedStateVariables_enkf:
                observationVariableValues[observationVariable].append(eval('siteInstance.site.' + observationVariable))

        return stateVariableValues, observationVariableValues


    def collectParameterValues(self, time, parameters_enkf):
        """
        Collect parameters as listed in parameters_enkf for all the plot instances in self.siteInstanceObservers[].
        """
        parameterValues = {}

        for parameter in parameters_enkf:
            parameterValues[parameter] = []

        for siteInstance in self.siteInstanceObservers:
            for parameter in parameters_enkf:
                parameterValues[parameter].append(eval('siteInstance.site.' + parameter))

        return parameterValues


    def kernelSmoothingParameters(self, parameterValues, parameter_bounds):
        """
        Kernel smoothing parameter values, follow eq.(30)
        """

        print('---------------------------------------------------------------------------------------')
        # Ref Eq.(8) in: Chen M,  Liu S,  Tieszen L. State-Parameter Estimation of Ecosystem Models Using a
        # Smoothed Ensemble Kalman Filter[J]. preprint, 2006.
        a = 0.49
        h2 = 1 - a**2  # h square
        for key in parameterValues.keys():
            theta_mean = np.mean(parameterValues[key])
            std = theta_mean*0.15
            tau = np.random.normal(0, std)

            for i in range(len(parameterValues[key])):
                while a*parameterValues[key][i] + (1-a)*theta_mean + h2*tau < parameter_bounds[key][0] or \
                    a * parameterValues[key][i] + (1 - a) * theta_mean + h2 * tau > parameter_bounds[key][1]:
                    tau = np.random.normal(0, std)
                parameterValues[key][i] = a*parameterValues[key][i] + (1-a)*theta_mean + h2*tau

        return parameterValues


    def changeInitialParameterValues(self, parameters_enkf, parameter_bounds):

        for key in parameters_enkf:
            parVal = np.random.uniform(parameter_bounds[key][0], parameter_bounds[key][1])
            for siteInstance in self.siteInstanceObservers:
                str = 'siteInstance.site.' + key + '=' + 'parVal'
                exec(str)


    def deliverChangeParameterValues(self, time, parameters_enkf, parameter_bounds):

        parameterValues = self.collectParameterValues(time, parameters_enkf)

        # Kernel smoothing parameter values
        changedParameterValues = self.kernelSmoothingParameters(parameterValues, parameter_bounds)

        # change parameter values for all the plot instances using changedParameterValues{}
        for index, siteInstance in enumerate(self.siteInstanceObservers):
            for key in changedParameterValues.keys():
                str = '{0}{1} = {2}[\'{3}\'][{4}]'.format('siteInstance.site.', key, 'changedParameterValues', key, index)
                exec(str)

        return changedParameterValues


    def changeParameterValues(self, updatedParameters):
        for index, siteInstance in enumerate(self.siteInstanceObservers):
            for key in updatedParameters.keys():
                str = '{0}{1} = {2}[\'{3}\'][{4}]'.format('siteInstance.site.', key, 'updatedParameters', key, index)
                exec(str)


    def changeStateValues(self, updatedStates):
        for index, siteInstance in enumerate(self.siteInstanceObservers):
            for key in updatedStates.keys():
                str = '{0}{1} = {2}[\'{3}\'][{4}]'.format('siteInstance.site.', key, 'updatedStates', key, index)
                exec(str)



    def deliverPrerunning(self, DB, standAge, endAge, observedStateVariables_enkf):
        """
        Collect model outputs (as listed in observedStateVariables_enkf[]) using default parameter values of the site.
        We use these outputs as the pseudo observed variable values.
        :return: A dict{} object with lists of outputs.
        """
        # Select plot instance with ensembleID of 0 to collect the pseudo observed data
        ensembleID = 0

        plot = self.siteInstanceObservers[ensembleID].site
        # print('MaxCond', plot.MaxCond)

        observedVariables = {}
        observedVariables['time'] = []
        for variable in observedStateVariables_enkf:
            observedVariables[variable] = []

        timer = 0
        while True:
            self.siteInstanceObservers[ensembleID].run3PG_2()
            # print(timer, self.siteInstanceObservers[0].site.GPP)

            # 12 means one observation per year
            if divmod(timer, 6)[1] == 0:
                for variable in observedStateVariables_enkf:
                    observedVariables[variable].append(eval('plot.' + variable))

                # add time record
                observedVariables['time'].append(timer)

            timer += 1
            standAge += 1 / 12.0
            if standAge - endAge > 0 and abs(standAge - endAge) > 0.01:
                break

        # Re-initialize the plot instance with ensembleID of 0.
        self.siteInstanceObservers[0] = main_module_EnKF.run3PG_2_class(plot, ensembleID, DB)

        return observedVariables


    def deliverForcing_Parameter_State_Observed_Values(self, forcings, parameters, stateVariables, observedStateVariables):
        pass


    def deliverPostExecutionProcesses(self):
        pass

    def subscribe(self, observer):
        self.siteInstanceObservers.append(observer)


    def runEnKF(self):
        for ensembleID in range(self.n):
            main_module_EnKF.run3PG_2_class(self.site, ensembleID)


    def getForcingDataVariances(self):
        """
        Get the variance of each selected forcing variable, the variances will be added to the respective forcing data
        in run3PG_2() during the running of the model, i.e., the disturbances of the forcing data are dynamically added
        to the original data per month.
        """

        self.forcingDataVariances['Tx'] = self.forcingDataVariance(self.site.mTx)
        self.forcingDataVariances['Tn'] = self.forcingDataVariance(self.site.mTn)
        self.forcingDataVariances['Tav'] = self.forcingDataVariance(self.site.mTav)
        self.forcingDataVariances['VPD'] = self.forcingDataVariance(self.site.mVPD)
        self.forcingDataVariances['Rain'] = self.forcingDataVariance(self.site.mRain)
        self.forcingDataVariances['SolarRad'] = self.forcingDataVariance(self.site.mRad)


    def forcingDataVariance(self, climateFactor):
        """
        Calculate the variance of each forcing data (weather data).
        :return: a list contains the variance of each month, with length of 12 from Jan to Dec.
        """

        if divmod(len(climateFactor), 12)[1] > 0.00001:     # not divisible by 12
            raise Exception('The total length of the climate data is not divisible by 12.')

        # number of years of the climate data
        nYear = divmod(len(climateFactor), 12)[0]

        # standard deviation
        stds = np.zeros(12)
        for i in range(12):
            temp = np.zeros(nYear)
            for j in range(nYear):
                temp[j] = climateFactor[i + j*12]
            stds[i] = np.std(temp)

        return stds