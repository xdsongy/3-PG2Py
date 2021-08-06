#-------------------------------------------------------------------------------
# Name:        main_module_EnKF
# Purpose:     Main entrance of the Py3PG2 for the dual state ensemble Kalman filter
#
#
# Author:      Xiaodong Song
#
# Created:     14/07/2021
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2021
#-------------------------------------------------------------------------------
#!/usr/bin/env python

# Retrieve single site parameters
import getSingleSiteParameters

# Global variables
import globalVariables

# Water balance module
import water_module

import numpy as np
import math
import re
import sqlite3
import time
import copy

#
# Refer run3PG_2() in VBA. 
# This is the 3-PG main routine, called whenever the model is to be run.
# All detailed 3-PG calculations are in do3PGmonthlyStep(). Here the model
# is initialised, the metdata and silvicultural events are assigned, and
# cumulative results are handled.
#
class run3PG2():
    
    def __init__(self, singleSiteInstance, ensembleID, bWriteOutput=True):
       
        # Generate global variables
        self.GlobalVar = globalVariables.globalVariablesClass()
       
        # Initiate self.site variable with the getSingleSiteParametersClass class
        # self.site = getSingleSiteParameters.getSingleSiteParametersClass()
        self.site = copy.deepcopy(singleSiteInstance)

        # In VB and FAST, the output table (sqlite) should be not written, i.e., self.bWriteOutput = False
        # self.bWriteOutput = bWriteOutput

        # ensemble ID
        self.ensembleID = ensembleID

        # self.run3PG_2(ensembleID)

        # Initialize the stand age
        self.initialiseStandAge()

        self.initialiseParameters()

        self.assignVaryData()

        self.initialiseStand()

        # The entrance to function initialiseSoilWaterBalance()
        self.waterModule = water_module.water_module_class(self.site)

        self.assignSilviculturalEvents(self.site.standAge)

        # set all output variables to be zero
        self.initialiseOutputData()

        # Note: the simulation starts with the month after planting as the first month of growth
        self.site.monthOfYear = self.site.InitialMonth + 1
        self.site.metMonth = self.site.monthOfYear
        self.site.metyear = self.site.InitialYear
        self.site.monthCounter = 1

        self.getAgeDependentfactors()

        self.getStandCharacteristics()

        # Headings and variables flagged as having initial values are output here
        # self.write3PGresults(self.GlobalVar.opStart, 0)

        # open the connection to sqlite DB, only need to be open in the 1st ensemble member
        # self.connection = sqlite3.connect(self.GlobalVar.outputFilename, check_same_thread=False)
        # self.connection = DB
        # self.cursor = self.connection.cursor()

        # Create a table in sqlite
        # self.writeOutput(0)

        # month counter, write into sqlite
        self.monthIncrease = 0


    def run3PG_2(self):

        # Do loop over monthly calculations
        # while self.site.standAge <= self.site.endAge - 0.01:
        # increase this counter
        self.monthIncrease += 1

        # Climate variables are zero based arrays, so self.site.metMonth shall minus 1
        self.assignMetdata(self.site.metMonth-1, self.site.metyear)
        # print(self.site.metMonth)

        # Here is the 3-PG run for the current month
        self.do3PGmonthlyStep(self.site.monthOfYear, self.site.standAge)

        # Now update stand age and deal with age-related factors, silvicultural events
        # and time-variant factors
        self.site.standAge += 1.0/12.0

        self.site.uLAI = self.fnUnderstoryLAI(self.site.standAge)
        self.site.pLAI = self.fnPastureLAI(self.site.standAge)

        self.getAgeDependentfactors()
        self.getStandCharacteristics()
        self.assignSilviculturalEvents(self.site.standAge)
        self.assignVaryData()

        # Update and handle output
        self.updateOutputData()

        if self.isOutputMonth(self.site.standAge, self.site.monthCounter):
            self.write3PGresults(self.GlobalVar.opEndMonth, self.site.monthCounter)

            if self.site.standAge < self.site.endAge - 0.01:
                self.initialiseOutputData()

        if self.site.monthOfYear < 12:
            self.site.monthOfYear += 1
        else:
            self.site.monthOfYear = 1

        if self.site.monthCounter < 12:
            self.site.monthCounter += 1
        else:
            self.site.monthCounter = 1

            # self.writeOutput(1)

    def run3PG2_Spatial(self):
        # increase this counter
        self.monthIncrease += 1

        # Climate variables are zero based arrays, so self.site.metMonth shall minus 1
        self.assignMetdata(self.site.metMonth - 1, self.site.metyear)

        # Here is the 3-PG run for the current month
        self.do3PGmonthlyStep(self.site.monthOfYear, self.site.standAge)

        # Now update stand age and deal with age-related factors, silvicultural events
        # and time-variant factors
        self.site.standAge += 1.0 / 12.0

        self.site.uLAI = self.fnUnderstoryLAI(self.site.standAge)
        self.site.pLAI = self.fnPastureLAI(self.site.standAge)

        self.getAgeDependentfactors()
        self.getStandCharacteristics()
        self.assignSilviculturalEvents(self.site.standAge)
        self.assignVaryData()

        # Update and handle output
        self.updateOutputData()

        if self.isOutputMonth(self.site.standAge, self.site.monthCounter):
            self.write3PGresults(self.GlobalVar.opEndMonth, self.site.monthCounter)

            if self.site.standAge < self.site.endAge - 0.01:
                self.initialiseOutputData()

        if self.site.monthOfYear < 12:
            self.site.monthOfYear += 1
        else:
            self.site.monthOfYear = 1

        if self.site.monthCounter < 12:
            self.site.monthCounter += 1
        else:
            self.site.monthCounter = 1

        # self.writeOutput(1)


    #
    # This procedure gets the starting month and initial stand age
    # It is assumed that InitialYear and InitialMonth are the calendar year
    # and month of the year at which the stand is initialised, and that
    # PlantedYear and PlantedMonth are the calendar year and month of the
    # year at which the stand was planted. If InitialYear is less than
    # PlantedYear, it is assumed to be referenced back to the PlantedYear
    # Refer function initialiseStandAge() in VBA
    def initialiseStandAge(self):

        if self.site.InitialYear < self.site.PlantedYear:
            self.site.InitialYear += self.site.PlantedYear
        
        self.site.standAge = (self.site.InitialYear + self.site.InitialMonth/12.0) - (self.site.PlantedYear + self.site.PlantedMonth/12.0)
        
        self.site.startAge = int(self.site.standAge)
        
        if self.site.standAge < 0:
            print ("Error: standAge of the site %s is error" % self.site.siteName)
            
        if self.site.standAge > self.site.endAge:
            print ("Error: standAge is greater than endAge in site %s" % self.site.siteName)
            

    #
    # Perform various parameter initialisations required before the model can be run
    # Refer function initialiseParameters() in VBA
    #
    def initialiseParameters(self):
        
        # assign the texture-dependent soil parameters for this soil class
        self.site.calculateSoilParameters()
            
        # check fN(FR) for no effect: fNn = 0 ==> fN(FR)=1 for all FR   
        if self.site.fNn == 0:
            self.site.fN0 = 1
        
        # derive the pFS-allometric parameters
        self.site.pfsPower = np.log(self.site.pFS20 / self.site.pFS2) / np.log(20 / 2)
        self.site.pfsConst = (self.site.pFS2 / 2) ** self.site.pfsPower
        
        # other checks
        self.site.poolFractn = max(0, min(1, self.site.poolFractn))  
        
            
    #
    # Assign values to Vary block data
    # Refer function assignVaryData() in VBA 
    #
    def assignVaryData(self):
        
        # I put this function here just for consistent with the 
        # coding in VBA version 3_PG, it does nothing in current stage
        
        if self.site.nVary == 0:
            pass        
          

    #
    # Perform various stand-related intialisations assign initial state of stand
    #
    def initialiseStand(self):
        
        # assign initial state of stand
        self.site.WLi = 0
        self.site.WS = self.site.WSi
        self.site.WF = self.site.WFi
        self.site.WR = self.site.WRi
        self.site.WL = self.site.WLi
        self.site.TotalW = self.site.WSi + self.site.WFi + self.site.WRi
        self.site.StemNo = self.site.StemNoi    
        self.site.zEffect = 1.0 + (self.site.edgeTrees / self.site.StemNo) * (self.site.edgeEffect / 100.0)
            
        self.site.uLAI = self.fnUnderstoryLAI(self.site.standAge)   
        self.site.pLAI = self.fnPastureLAI(self.site.standAge)        
        self.site.noRootZone = (self.site.spRootVol == 0) or (self.site.soilDepth == 0) or (self.site.soilDepth == self.GlobalVar.MissingValue)
 
        # sundry stuff
        self.site.thinEventNo = 1
        self.site.dfolEventNo = 1
        self.site.MAIx = 0
        self.site.LAIx = 0
        self.site.avStemGR = 0
        
        # set output frequency
        if len(self.site.outputDates[0]) == 0:
            self.site.noOutputDates = 0
        else:
            self.site.noOutputDates = len(self.site.outputDates)
        
        if self.site.noOutputDates > 0:
            self.site.opFrequency = self.GlobalVar.opfMonthly
        
        self.site.outputDateNo = 1
        
    #
    # Get understory LAI - a stand age determined variable
    # Refer function fnUnderstoryLAI() in VBA 
    #
    def fnUnderstoryLAI(self, age):
        
        return self.fnGenericLAI(age, 
                          self.site.uAge1, self.site.uLAI1, 
                          self.site.uAge2, self.site.uLAI2, 
                          self.site.uAge3, self.site.uLAI3, 
                          self.site.uIndex)
        
        
    #
    # Get pasture LAI - a stand age determined variable
    # Refer function fnPastureLAI() in VBA
    #
    def fnPastureLAI(self, age):
        
        return self.fnGenericLAI(age, 
                          self.site.pAge1, self.site.pLAI1, 
                          self.site.pAge2, self.site.pLAI2, 
                          self.site.pAge3, self.site.pLAI3, 
                          self.site.pIndex)
        
  
    #
    # The following deal with canopy and understory cover, and light interception
    # Refer function fnGenericLAI() in VBA 
    #
    def fnGenericLAI(self, 
                     age, 
                     Age0, L0,
                     Age1, L1,
                     Age2, L2, 
                     n):
       
        # Get understory cover
        y = 0
        
        if n <= 0:
            if age <= Age0:
                y = L0
            elif age <= Age1:
                y = L1
            else:
                y = L2 
        else:
            if age <= Age0:
                y = L0
            elif age <= Age1:
                a = (2 ** (n - 1)) * (L1 - L0) / (Age1 - Age0) ** n
                if age < 0.5 * (Age1 + Age0):
                    y = L0 + a * (age - Age0) ** n
                else:
                    y = L1 - a * (Age1 - age) ** n
            elif age < Age2:
                a = (2 ** (n - 1)) * (L2 - L1) / (Age2 - Age1) ** n
                if age < 0.5 * (Age2 + Age1):
                    y = L1 + a * (age - Age1) ** n
                else:
                    y = L2 - a * (Age2 - age) ** n
            else:
                y = L2
        
        return max(y, 0)
        

    #
    # Do age-dependent silvicultural events
    # Refer function assignSilviculturalEvents() in VBA
    #
    def assignSilviculturalEvents(self, age):
        
        # Assign values to silvicultural events
        if self.site.nFR > 0:
            self.site.FR = self.fnLookup(age, self.site.nFR, self.site.FRages, self.site.FRVals)
        
        if self.site.nMinASW > 0:
            self.site.minASW = self.fnLookup(age, self.site.nMinASW, self.site.minASWages, self.site.minASWVals)
 
        if self.site.nIrrig > 0:
            self.site.applIrrig = self.fnLookup(age, self.site.nIrrig, self.site.irrigAges, self.site.irrigVals)
 
 
    #
    # Perform a look-up in a table to determine a value of a function.
    # The boolean "interpolateLookups" is a global parameter that determines
    # if the look-ups are based on a histogram function or linear interpolation
    # Refer function fnLookup() in VBA 
    #
    def fnLookup(self, x, nTable, xTable, yTable):
        
        if self.GlobalVar.interpolateLookups:
            return self.fnInterpolateValue(x, nTable, xTable, yTable)
        else:
            return self.fnBinValue(x, nTable, xTable, yTable)
 
 
    #
    # Perform a linear interpolation between values in a lookup table to determine
    # the value of "table" at "x".
    # Refer function fnInterpolateValue() in VBA
    # @param x type(int) Age
    # @param nTable type(int) Events number
    # @param xTable type(list) Age list
    # @param yTable type(list) Corresponding value list vs. Age
    #
    def fnInterpolateValue(self, x, nTable, xTable, yTable):
        
        if x <= xTable[0]:
            return yTable[0]
        elif x >= xTable[nTable-1]:
            return yTable[nTable-1]
        else:
            i = 0
            while x >= xTable[i]:
                i += 1
            
            slope = (yTable[i] - yTable[i - 1]) / (xTable[i] - xTable[i - 1])
            f = yTable[i - 1] + slope * (x - xTable[i - 1])
            return f
 
 
    #
    # Perform a look-up in a table of bin-values. The values in the bins apply
    # for values of "x" up to and including the "x"-values in "table"
    # Refer function fnBinValue() in VBA
    # @param x type(int) Age
    # @param nTable type(int) Events number
    # @param xTable type(list) Age list
    # @param yTable type(list) Corresponding value list vs. Age
    #
    def fnBinValue(self, x, nTable, xTable, yTable):
        
        i = 0
        while x > xTable[i] or i < nTable-1:
            i += 1
        
        if i <= nTable-1:
            return yTable[i]
        else:
            return yTable[nTable]
        
 
    #
    # Initialise and update the cumulative output variables
    # Refer function initialiseOutputData() in VBA
    #
    def initialiseOutputData(self):
        
        self.site.cSelfThin = 0
        self.site.cMortality = 0
        self.site.cRadInt = 0
        self.site.cGPP = 0
        self.site.cNPP = 0
        self.site.cCVI = 0
        self.site.cLitter = 0
        self.site.cStemDM = 0
        self.site.cfc_Transp = 0
        self.site.cpc_Transp = 0
        self.site.cfu_Transp = 0
        self.site.cfs_Evap = 0
        self.site.cps_Evap = 0
        self.site.cET = 0
        self.site.cTransp = 0
        self.site.cEsoil = 0
        self.site.cRainInt = 0
        self.site.cEsoil = 0
        self.site.cf_RainInt = 0
        self.site.cp_RainInt = 0
        self.site.cRunoff = 0
        self.site.cf_Runoff = 0
        self.site.cp_Runoff = 0
        self.site.cSupIrrig = 0            
            
            
    #
    # The following deal with age-dependent and stand factors
    # Refer function getAgeDependentfactors() in VBA 
    #
    def getAgeDependentfactors(self):
       
        # Determine the age-dependent factors
        self.site.SLA = self.fnExp(self.site.standAge, self.site.SLA0, self.site.SLA1, self.site.tSLA, 2)
        self.site.fracBB = self.fnExp(self.site.standAge, self.site.fracBB0, self.site.fracBB1, self.site.tBB, 1)
        self.site.rho = self.fnExp(self.site.standAge, self.site.rho0, self.site.rho1, self.site.tRho, 1)
        self.site.CanopyCover = self.fnCanopyCover(self.site.standAge)
          
 
    #
    # This is a standardised function used for various effects in 3-PG
    # Refer function fnExp() in VBA 
    #
    def fnExp(self, x, g0, gx, tg, ng):
        
        if tg != 0:
            y = gx + (g0 - gx) * np.exp(-np.log(2) * (x / tg) ** ng)
        else:
            y = gx
            
        return y
 

    #
    # Get canopy cover - a stand age determined variable
    # Refer function fnCanopyCover() in VBA 
    #
    def fnCanopyCover(self, age):
        
        if (age < self.site.fullCanAge) and (self.site.fullCanAge > 0):
            y = (age + 0.01) / self.site.fullCanAge
        else:
            y = 1
            
        return y
            
    
    #
    # Determine the stand characteristics
    # Refer function getStandCharacteristics() in VBA 
    #
    def getStandCharacteristics(self):
        
        self.site.LAI = self.site.WF * self.site.SLA * 0.1
        
        if self.site.CanopyCover * self.site.treeCover == 0:
            self.site.fLAI = 0
        else:
            self.site.fLAI = self.site.WF * self.site.SLA * 0.1 / (self.site.CanopyCover * self.site.treeCover)
            
        self.site.avWS = self.site.WS * 1000 / self.site.StemNo
        self.site.avDBH = (self.site.avWS / self.site.aWs) ** (1 / self.site.nWs)
        self.site.BasArea = (((self.site.avDBH / 200) ** 2) * np.pi) * self.site.StemNo
        self.site.Height = self.site.aH * (self.site.avDBH ** self.site.nHB) * (self.site.StemNo ** self.site.nHN)
        
        if self.site.aV > 0:
            self.site.StandVol = self.site.aV * self.site.avDBH ** self.site.nVB * self.site.StemNo ** self.site.nVN
        else:
            self.site.StandVol = self.site.WS * (1 - self.site.fracBB) / self.site.rho
        
        self.site.CVI = self.site.StandVol - self.site.oldVol
        self.site.oldVol = self.site.StandVol
        
        if self.site.standAge > 0:
            self.site.MAI = self.site.StandVol / self.site.standAge
        else:
            self.site.MAI = 0.0   
    
    
    #
    # This is the primary output routine called by run3PG. It uses the value of
    # "runType" to determine which output procedure to call.
    #
    # The value of "action" depends on the stage in a run at which the output
    # procedure was called, and hence determines what is to be done at that stage.
    #
    def write3PGresults(self, action, month):
        pass
        
 
    #
    # Assign metdata for current month m in year y
    # Refer function assignMetdata() in VBA 
    # @param m type(int) Month
    # @param y type(int) Year
    #    
    def assignMetdata(self, m, y):
        
        # Here the coding is a little different with the original VBA code,
        # we omit the first part of the "If isDailyMetdata Then", because
        # isDailyMetdata is always False in current development stage 
        if not self.site.isDailyMetdata:
            self.site.SolarRad = self.site.mRad[m]
            self.site.Tx = self.site.mTx[m]
            self.site.Tn = self.site.mTn[m]
            self.site.Tav = self.site.mTav[m]
            self.site.VPD = self.site.mVPD[m]            
            self.site.FrostDays = math.ceil(self.site.mFrostDays[m])           
            self.site.RainDays = math.ceil(self.site.mRainDays[m])
            self.site.DayLength = self.site.mDayLength[m]
            self.site.Rain = self.site.mRain[m]

            # print('m=',m)
            # print(self.site.Tx)
            # print(self.site.Tn)
            # raise Exception('print (m)')

            if self.site.Rain == 0:
                self.site.RainDays = 1
                
            if self.site.metMonth < 12 * self.site.mYears - 1:
                self.site.metMonth += 1
            else:
                self.site.metMonth = 0
                
            #self.site.metMonth = self.fnNextValue(m, 12 * self.site.mYears)
            

    #
    # Increment n by 1, wrapping back to 1 after max
    #
    def fnNextValue(self, n, max):
        
        if n < max-1:
            n += 1
        else:
            n = 0
            
        return n
        
       
    #
    # This is the monthly, i.e. inner-most, loop in the 3PG model
    # Do the monthly calculations - i.e. all growth & mortality processes
    # Refer function do3PGmonthlyStep() in VBA
    # @param month type(int) month of a year, between 1 to 12
    # @param age type(float) standAge
    #    
    def do3PGmonthlyStep(self, month, age):
        
        # Get growth modifiers and light interception data
        self.getGrowthModifiers()

        # Do the soil water balance
        #... this will be either monthly or daily depending on selected mode
        # Note: we only implement the monthly mode in current stage
        self.doSoilWaterBalance(month, age)
  
        # Determine biomass increments and losses and then update biomass pools
        self.getNPP(month)
        self.getAllocationRates()
        self.getLossRates()
        self.doBiomassAllocation()
        
        # Determine age and stress-related stem mortality, and do self-thinning
        self.doMortality(age)
        self.doSelfThinning()
        
        if self.site.StemNo <= 0:
            print ("Error: the StemNo is zero for plot: %s" % self.site.siteName)
            return
        
        # Update tree and stand data at the end of this time period,
        # taking thinning or defoliation into account
        # Note: minus 1 means the index of the event list begin from zero
        if self.site.thinEventNo <= self.site.nThin:
            self.doThinning(age, self.site.thinEventNo-1) 
        
        if self.site.dfolEventNo <= self.site.nDfol:
            self.doDefoliation(age, self.site.dfolEventNo-1)
        
        self.site.TotalW = self.site.WF + self.site.WR + self.site.WS
        
        self.site.WUE = self.fnEfficiency(self.site.NPP, self.site.ET)
        
        self.finaliseSoilwaterBalance()
        
   
    #
    # Calculate all growth modifiers
    # Refer function getGrowthModifiers() in VBA 
    #
    def getGrowthModifiers(self):
        
        if self.GlobalVar.use3PG1modifiers:
            self.site.gmTemp = self.growthModifierTemp(self.site.Tav)
        else:
            self.site.gmTemp = self.growthModifierTemp(self.site.Tx)
 
        self.site.gmVPD = self.growthModifierVPD(self.site.VPD)
        self.site.gmSW = self.growthModifierSW(self.site.rzRelASW)
        self.site.gmSalt = self.growthModifierSalt(self.site.EC)
        self.site.gmNut = self.growthModifierNut(self.site.FR)
        self.site.gmFrost = self.growthModifierFrost(self.site.FrostDays)
        self.site.gmAge = self.growthModifierAge(self.site.standAge)
        self.getCO2modifiers()
        self.site.physMod = min(self.site.gmVPD, self.site.gmSW) * self.site.gmAge * self.site.gmSalt
        
 
    #
    # Get the temperature-dependent growth modifier
    #
    def growthModifierTemp(self, t):
        
        p = ((self.site.Tmax - self.site.Topt) / (self.site.Topt - self.site.Tmin))
        
        if t <= self.site.Tmin or t >= self.site.Tmax:
            gm = 0.0
        else:
            gm = ((t - self.site.Tmin) / (self.site.Topt - self.site.Tmin)) * ((self.site.Tmax - t) / (self.site.Tmax - self.site.Topt)) ** p
 
        return gm
            
            
    #
    # Get the VPD-dependent growth modifier
    #
    def growthModifierVPD(self, VPD):
        
        return np.exp(-self.site.CoeffCond * VPD)
            
            
    #
    # Get the soilwater-dependent growth modifier
    #
    def growthModifierSW(self, rASW):
     
        c = self.site.cTheta ** self.site.nTheta
        w = (1 - rASW) ** self.site.nTheta  
        
        return  (1 - w) / (1 + (1 - 2 * c) * w / c)
            
            
    #
    # Get the EC-dependent growth modifier
    #
    def growthModifierSalt(self, EC):
            
        if EC <= self.site.EC0:
            gm = 1.0
        elif EC < self.site.EC1:
            gm = 1 - ((EC - self.site.EC0) / (self.site.EC1 - self.site.EC0)) ** self.site.ECn
        else:
            gm = 0.0
            
        return gm
  
  
    #
    # Get the nutrition-dependent growth modifier
    #
    def growthModifierNut(self, FR):
        
        if self.site.fNn == 0:
            return 1
        else:
            return 1 - (1 - self.site.fN0) * (1 - FR) ** self.site.fNn
  
  
    #
    # Get the frost-dependent growth modifier
    #
    def growthModifierFrost(self, t):
        
        return 1 - self.site.kF * (t / 30.0)
  
  
    #
    # Get the age-dependent growth modifier
    #
    def growthModifierAge(self, age):
        
        gm = 1
        
        if self.site.nAge > 0:
            relAge = age / self.site.MaxAge
            gm = (1.0 / (1 + (relAge / self.site.rAge) ** self.site.nAge))
            
        return gm


    def getCO2modifiers(self):
        
        self.site.gmCO2NPP = 1.0
        self.site.gmCO2gC  = 1.0


    #
    # 
    #
    def doSoilWaterBalance(self, month, age):
        
        if self.site.isDailyMetdata:
            self.waterModule.doDailySoilWaterBalance(month, age)
        else:
            self.waterModule.doMonthlySoilWaterBalance(month, age)


    #
    # Determine NPP ...
    #
    def getNPP(self, month):
        
        # Get radiation use efficiency and intercepted radiation
        if self.GlobalVar.use3PG1modifiers:
            partMod = self.site.physMod
        else:
            partMod = self.site.gmSW * self.site.gmAge * self.site.gmSalt

        self.site.alphaC = self.site.alphaCx * self.site.gmNut * self.site.gmTemp * self.site.gmFrost * self.site.gmCO2NPP * partMod
        self.site.epsilon = self.site.gDM_mol * self.site.molPAR_MJ * self.site.alphaC
        self.getCanopyRadiationInterception (self.site.SolarRad)
        self.site.RadInt = self.GlobalVar.daysInMonth[month-1] * self.site.CanopyCover * self.site.treeCover * self.site.fRadInt

        # fracET corrects for low SW, 100 ==> tDm/ha
        self.site.GPP = self.site.fracET * self.site.epsilon * self.site.RadInt * self.site.zEffect / 100
        self.site.NPP = self.site.GPP * self.site.y
        

    #
    # Get radiation intercepted by forest canopy
    #
    def getCanopyRadiationInterception(self, Rad):

        self.site.fFracInt = self.fnFractionIntercepted(self.site.k, self.site.fLAI)
        self.site.fRadInt = self.site.fFracInt * Rad
        

    #
    # Get fraction of incident light intercepted by canopy of effective LAI = L
    #
    def fnFractionIntercepted(self, k, L):
        
        return 1 - math.exp(-k * L)


    #
    # Get biomas allocation coefficients
    #
    def getAllocationRates(self):
        
        self.site.pFS = self.site.pfsConst * self.site.avDBH ** self.site.pfsPower
        self.site.pR = self.site.pRx * self.site.pRn / (self.site.pRn + (self.site.pRx - self.site.pRn) * self.fnRootAlloc())
        self.site.pS = (1 - self.site.pR) / (1 + self.site.pFS)
        self.site.pF = 1 - self.site.pR - self.site.pS
        

    #
    # Get the root allocation modifier
    #
    def fnRootAlloc(self):
        
        m = self.site.m0 + (1 - self.site.m0) * self.site.FR
        
        if self.GlobalVar.use3PG1modifiers:
            return m * self.site.physMod
        else:
            return m * self.site.gmSW

        
    #
    # Get loss rates
    #
    def getLossRates(self):
        
        self.site.gammaF = self.fnGammaF(self.site.standAge)
        self.site.gammaR = self.fnGammaR()


    #
    # Get literfall rate. Currently this is simply age dependent
    #
    def fnGammaF(self, age):

        if self.site.tgammaF * self.site.gammaF1 == 0:
            gF = self.site.gammaF1
        else:
            kF = 12 * math.log(1 + self.site.gammaF1 / self.site.gammaF0) / self.site.tgammaF
            gF = self.site.gammaF1 * self.site.gammaF0 / (self.site.gammaF0 + (self.site.gammaF1 - self.site.gammaF0) * math.exp(-kF * age))
        
        return gF


    #
    # Get root-turnover rate - currently a constant parameter
    #
    def fnGammaR(self):
        
        return self.site.gammaR


    #
    # Do the biomass allocation and take losses into account
    #
    def doBiomassAllocation(self):

        self.site.incrWF = self.site.NPP * self.site.pF
        self.site.incrWR = self.site.NPP * self.site.pR
        self.site.incrWS = self.site.NPP * self.site.pS
        self.site.lossWF = self.site.gammaF * self.site.WF        
        self.site.lossWR = self.site.gammaR * self.site.WR
        self.site.WF += self.site.incrWF - self.site.lossWF
        self.site.WR += self.site.incrWR - self.site.lossWR
        self.site.WS += self.site.incrWS
        self.site.WL += self.site.lossWF


    #
    # Determine stems removed by age- & stress-related mortality
    # gammaN is % per annum
    # fnStressMortality is fraction dying per month
    #
    def doMortality(self, age):
        
        self.site.gammaN = self.fnExp(age, self.site.gammaN0, self.site.gammaN1, self.site.tgammaN, self.site.ngammaN)
        gN = self.site.gammaN / 12.0 / 100.0

        self.doUpdateStems(gN, self.site.mortality)


    #
    # Update biomass pools to take mortality into account, where
    #  gN = fraction of stems that have died
    #  delStems = number that die
    #
    def doUpdateStems(self, gN, delStems):
        
        delStems = gN * self.site.StemNo
        self.site.WF -= self.site.mF * gN * self.site.WF
        self.site.WR -= self.site.mr * gN * self.site.WR
        self.site.WS -= self.site.mS * gN * self.site.WS
        self.site.StemNo -= delStems


    #
    # Calculate self-thinning mortality
    #
    def doSelfThinning(self):

        self.site.selfThin = 0
        if self.site.StemNo > 0:
            self.site.wSmax = self.fnMaxWs(self.site.StemNo)  
            self.site.avWS = self.site.WS * 1000.0 / self.site.StemNo
            
            if self.site.wSmax < self.site.avWS:  # apply self-thinning
                gN = self.fnSelfthinningMortality(self.site.StemNo, self.site.WS)   
                self.doUpdateStems(gN, self.site.selfThin)     
                self.site.wSmax = self.fnMaxWs(self.site.StemNo)   
                self.site.avWS = self.site.WS * 1000.0 / self.site.StemNo
            

    #
    # Get the maximum stem mass allowed by selfthinning when stockign is N
    #
    def fnMaxWs(self, n):
        
        return self.site.wSx1000 * (1000.0 / n) ** self.site.thinPower


    #
    # Determine fraction of stems removed by self-thinning by applying the
    # Newton-Rhapson method to solve for the new stem numbers
    #
    def fnSelfthinningMortality(self, oldN, oldW):

        # accuracy of convergence
        accuracy = 1 / 1000.0         
        
        # maximum number of iterations
        maxItern = 5                

        n = oldN / 1000.0
        x1 = 1000 * self.site.mS * oldW / oldN
        i = 0
        
        # apply Newton-Rhapson method 
        dN = 0.0
        while math.fabs(dN) >= accuracy or i <= maxItern:
            i += 1
            x2 = self.site.wSx1000 * n ** (1 - self.site.thinPower)
            y = x2 - x1 * n - (1 - self.site.mS) * oldW
            dy = (1 - self.site.thinPower) * x2 / n - x1
            dN = -y / dy
            n = n + dN 
            if n < 0.0:
                print ("Error: the stem number has become negative.")
            
        return (oldN - 1000.0 * n) / oldN
        

    #
    # If it is time to do a thinning, carry out the thinning (if there are
    # stems to remove) and update the thinnning event number "eventNo"
    # @param age type(float) stand age
    # @param eventNo type(int) the index of the events, zero based
    #    
    def doThinning(self, age, eventNo):        
        if math.fabs(age - self.site.thinAges[eventNo]) < 0.01:
            if self.site.StemNo > self.site.thinVals[eventNo]:
                delN = (self.site.StemNo - self.site.thinVals[eventNo]) / self.site.StemNo
                self.site.StemNo = self.site.StemNo * (1 - delN)
                self.site.WF = self.site.WF * (1 - delN * self.site.thinWF[eventNo])
                self.site.WR = self.site.WR * (1 - delN * self.site.thinWR[eventNo])
                self.site.WS = self.site.WS * (1 - delN * self.site.thinWS[eventNo])
                
            self.site.thinEventNo += 1
                
        
    #
    # If it is the time of a defoliation, carry out the defoliation
    # and update the defoliation event number "eventNo"
    #
    def doDefoliation(self, age, eventNo):
        
        if age >= self.site.defolAges[eventNo]:
            self.site.WF = self.site.WF * self.site.defolVals[eventNo]
            self.site.dfolEventNo += 1


    #
    # Get an efficiency
    #
    def fnEfficiency(self, qOut, qIn):

        if qIn == 0:
            return 0
        else:
            return 100.0 * qOut / qIn


    #
    # Little things that need to be done prior to the next time step
    # 
    def finaliseSoilwaterBalance(self):

        self.doGrowRootZone()

        self.site.rzRelASW = self.waterModule.fnRelASW(self.site.rzSW, self.site.rzVol)
        self.site.nrRelASW = self.waterModule.fnRelASW(self.site.nrSW, self.site.soilVol - self.site.rzVol)


    #
    # As roots grow the volume of soil they access changes. This routine determines the
    # new volume of the root zone due to root growth, with an upper limit being the full
    # soil profile, and then takes into account the soil water that is thus moved from the
    # root-free zone to the root zone, or conversely if the root zone volume decreases.
    #
    def doGrowRootZone(self):

        if self.site.noRootZone:
            return
        
        oldRootDepth = self.site.rootDepth
        self.site.rootDepth = self.waterModule.fnRootDepth(self.site.WR)    
        self.site.rzVol = self.site.treeCover * self.site.rootDepth
        delDepth = self.site.rootDepth - oldRootDepth
        
        if self.site.soilDepth - oldRootDepth == 0:
            delSW = 0
        elif delDepth > 0:
            delSW = delDepth * self.site.nrSW / (self.site.soilDepth - oldRootDepth)
        else:
            delSW = delDepth * self.site.rzSW / oldRootDepth
        
        self.site.rzSW += delSW
        self.site.nrSW -= delSW   
        self.site.rootFrac = self.site.rootDepth / self.site.soilDepth


    #
    #
    #
    def updateOutputData(self):
        
        self.site.cSelfThin  += self.site.selfThin
        self.site.cMortality += self.site.mortality
        self.site.cRadInt    += self.site.RadInt
        self.site.cGPP       += self.site.GPP
        self.site.cNPP       += self.site.NPP
        self.site.cCVI       += self.site.CVI
        self.site.cLitter    += self.site.lossWF
        self.site.cStemDM    += self.site.incrWS
        self.site.cfc_Transp += self.site.fTransp
        self.site.cpc_Transp += self.site.pTransp
        self.site.cfu_Transp += self.site.uTransp
        self.site.cET        += self.site.ET
        self.site.cTransp    += self.site.Transp
        self.site.cfs_Evap   += self.site.fEsoil
        self.site.cps_Evap   += self.site.pEsoil
        self.site.cEsoil     += self.site.Esoil
        self.site.cf_RainInt += self.site.fRainInt
        self.site.cp_RainInt += self.site.pRainInt
        self.site.cRainInt   += self.site.RainInt
        self.site.cf_Runoff  += self.site.fRunoff
        self.site.cp_Runoff  += self.site.pRunoff
        self.site.cRunoff    += self.site.RunOff
        self.site.cSupIrrig  += self.site.SupIrrig
        
        # Calculate current light & water uses efficiencies
        self.site.cEpsilonWS = self.fnEfficiency(self.site.cStemDM, self.site.cRadInt)
        self.site.cEpsilon = self.fnEfficiency(self.site.cGPP, self.site.cRadInt)
        self.site.cWUE = self.fnEfficiency(self.site.cNPP, self.site.cET)
        
        # Update peak LAI & MAI and age at peaks
        # Refer function updateMaximum() in VBA
        if self.site.fLAI > self.site.LAIx:
            self.site.LAIx = self.site.fLAI
            self.site.ageLAIx = self.site.standAge        
        
        if self.site.MAI > self.site.MAIx:
            self.site.MAIx = self.site.MAI
            self.site.ageMAIx = self.site.standAge


    #
    # Is output required this month? If so, the results of this month will be outputed
    #
    def isOutputMonth(self, age, month):
            
        yes = False
        
        if self.site.noOutputDates > 0:
            yes = self.isOutputDate(age)
        else:
            yes = ((self.site.opFrequency == self.GlobalVar.opfMonthly) or 
                   ((self.site.opFrequency == self.GlobalVar.opfAnnual) and (month == 12)) or 
                   ((self.site.opFrequency == self.GlobalVar.opfRotation) and (age >= self.site.endAge - 0.001)))

        return yes


    #
    # Check to see if the current stand age corresponds to an output date
    #
    def isOutputDate(self, age):
        
        if self.site.outputDateNo <= self.site.noOutputDates:
            # similar with '3y11m'
            date = self.site.outputDates[self.site.outputDateNo-1]              
            parsedAge = re.split('[y,m]', date)
            y = int(parsedAge[0])   # year
            m = int(parsedAge[1])   # month
            d = y + m / 12.0
            yes = math.fabs(age - d) < 0.01
            
            if yes:
                self.site.outputDateNo += 1
        else:
            yes = False
            
        return yes
        

    #
    # For each site, create one sqlite table using siteName, and insert records into it
    #
    def writeOutput(self, dbWriteFlag):
        pass
        # Generate table name using self.site.siteName.
        # Note: some siteName contain space, like "Barrett 11", it's not allowed to contain 
        # space in sqlite table name, so here we replace space(" ") using "_"
        # tableName = self.site.siteName.replace(' ', '_') + '_' + str(self.ensembleID)
        # tableName = tableName.replace('+', 'P')   # 'P' means '+', e.g. Cann River PF + Mac
        # tableName = tableName.replace('-', 'M')   # 'M' means '-', e.g. Cann River PF - Mac
        #
        # # create table
        # if dbWriteFlag == 0:
        #     # delete table if exists
        #     sqlClause = 'DROP TABLE IF EXISTS ' + tableName
        #     self.cursor.execute(sqlClause)
        #     self.connection.commit()
        #
        #     sqlClause = "CREATE TABLE " + tableName + " (" + \
        #     "SimTime integer," + \
        #     "standAge float," + \
        #     "Stems integer," + \
        #     "MeanDBH float," + \
        #     "BasalArea float," + \
        #     "StandVolume float," + \
        #     "MAI float," + \
        #     "PlotlevelLAI float," + \
        #     "FoliageDM float," + \
        #     "StemDM float," + \
        #     "RootDM float," + \
        #     "BarkBranchFraction float" + \
        #      ");"
        #
        #     self.cursor.execute(sqlClause)
        #
        # # insert records
        # if dbWriteFlag == 1:
        #     sqlClause = "INSERT INTO " + tableName + " values(" + \
        #     str(self.monthIncrease) + "," + \
        #     str(self.site.standAge) + "," + \
        #     str(self.site.StemNo) + "," + \
        #     str(self.site.avDBH) + "," + \
        #     str(self.site.BasArea) + "," + \
        #     str(self.site.StandVol) + "," + \
        #     str(self.site.MAI) + "," + \
        #     str(self.site.LAI) + "," + \
        #     str(self.site.WF) + "," + \
        #     str(self.site.WS) + "," + \
        #     str(self.site.WR) + "," + \
        #     str(self.site.fracBB) + \
        #     ");"
        #
        #     self.cursor.execute(sqlClause)