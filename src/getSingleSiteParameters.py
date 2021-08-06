#-------------------------------------------------------------------------------
# Name:        getSingleSiteParameters
# Purpose:     Retrieve the parameters for a single site (for each plot in single excel sheet), 
#              the site series can be get from module "getSiteSeries"
#
#
# Author:      Xiaodong Song
#
# Created:     09/08/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python

#
# Module for excel file reading 
#
import xlrd

import numpy as np

#
# Global variables
#
from globalVariables import globalVariablesClass

class getSingleSiteParametersClass():
    
    def __init__(self, soilDB={}, paraDict={}, excel_file='', site=''):
        
        if excel_file != '' and site != '':
            
            # Refer function readSingleSiteData() in VBA
            self.Initialise3PGxlInput()
            self.initialiseSiteData()
            self.initialiseStandData()
            self.assignDefaultParameters()
            self.assignHiddenParameters()
            
            # Assign parameter dictionary to local variable self.paraDict
            self.paraDict = paraDict
            
            # Using true parameters' value to overlap values given by function assignDefaultParameters()
            self.assignTrueParameters() 
                        
            # Set the initial excel file path and sheet name
            self.excel_file = excel_file
            self.siteName = site
            
            # Assign soil database to single site local variable
            self.soilDB = soilDB
            
            # Call readSiteParaemters() to initialize single site parameters
            self.readParaemters()
         
    
    #
    # Fill the dictionary with data read from excel "Standard 3PGxl soil types" block
    # Refer funtion readParaemters() in VBA
    def readParaemters(self):
        
        # Open the workbook
        self.wb = xlrd.open_workbook(self.excel_file, 'rb')
        
        # Open the sheet which contains the certain site parameters
        self.sh = self.wb.sheet_by_name(self.siteName)  
        
        self.readOutputAges()
        
        self.readSiteFactors()
        
        self.readStandData()
        
        self.getChangedParameterValues()
        
        self.readVaryBlocks()
        
        self.assignInitialBiomassPools() 
        
        self.readUnderstoryData()
        
        self.readPastureData()
        
        self.readDailyWaterBalance()
        
        self.readClimateData()
               
        self.readSilvicultualEvents()
        
        self.readObservedData()
        
        # print ("Site: %-15s ......initialized." % self.siteName)
        
        
    #
    # Get site-specific changes to current 3-PG parameter values, refer function getChangedParameterValues() in VBA 
    # Note: I don't think site specific parameters are useful, so I will just omit this fuction and left it as a 
    # reminder only
    def getChangedParameterValues(self):
        pass
        
        
    #
    # This code reads the data for Vary sub-blocks. The data is stored in one dimensional
    # arrays and is parsed from these when it is assigned to the paramters or factors
    # in the tPG_Model module.
    # Refer function readVaryBlocks() in VBA. Because the "Vary block" is missing in our data(excel file),
    # so I will just omit this fuction and left it as a reminder only
    #
    def readVaryBlocks(self):
        pass
        
    
    #
    # Determine and assing the initial biomass pools, refer function assignInitialBiomassPools() in VBA 
    #
    def assignInitialBiomassPools(self):
        
        # derived from seedling mass and stand density
        if self.SeedlingMass > 0:
            self.WFi = (self.fracWF / (self.fracWF + self.fracWR + self.fracWS)) * self.SeedlingMass * self.StemNoi / 10 ** 6
            self.WRi = (self.fracWR / (self.fracWF + self.fracWR + self.fracWS)) * self.SeedlingMass * self.StemNoi / 10 ** 6
            self.WSi = (self.fracWS / (self.fracWF + self.fracWR + self.fracWS)) * self.SeedlingMass * self.StemNoi / 10 ** 6
        # derived from total stand biomass
        elif self.StandMass > 0:
            self.WFi = (self.fracWF / (self.fracWF + self.fracWR + self.fracWS)) * self.StandMass
            self.WRi = (self.fracWR / (self.fracWF + self.fracWR + self.fracWS)) * self.StandMass
            self.WSi = (self.fracWS / (self.fracWF + self.fracWR + self.fracWS)) * self.StandMass
        else:
            # nothing to do - the initial pool biomass data are read directly
            pass 
        
        # print "SeedlingMass = %s, StandMass = %s" % (self.SeedlingMass, self.StandMass)
     
           
    #
    # Get "Output ages" corresponding with the observed data
    #
    def readOutputAges(self):
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        self.locateCell("Output ages", cellLocation)
        self.outputDates = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value.split(',')
        
        
    #
    # Get "Site data" block parameters (refer sub routine 'readSiteFactors()' in VBA)
    # 
    def readSiteFactors(self):        
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        self.locateCell("Latitude =", cellLocation)
        self.siteLat = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Fertility rating =", cellLocation)
        self.FR = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Soil class =", cellLocation)
        self.soilClass = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Maximum ASW =", cellLocation)
        self.maxASW = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        if self.maxASW == '':
            self.maxASW = globalVariablesClass.MissingValue
        
        self.locateCell("Salinity =", cellLocation)
        self.EC = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Ambient CO2 =", cellLocation)
        self.CO2 = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Rain intensity =", cellLocation)
        self.rainIntensity = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.isSoilProfile = (self.maxASW == globalVariablesClass.MissingValue)
        
        if self.isSoilProfile:
            self.locateCell("Soil depth =", cellLocation)
            self.soilDepth = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value 
            
            if not self.locateCell("Watertable depth", cellLocation):
                self.wtblDepth = globalVariablesClass.MissingValue
            else:
                self.wtblDepth = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            
            self.locateCell("%Stones =", cellLocation)
            self.pStones = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value 
        else:
            self.locateCell("Minimum ASW =", cellLocation)
            self.minASW = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            
        # read soil class data
        if self.soilClass != '?':
            self.readSoilDatabase(self.soilClass)
        
        self.readSoilData()
            
    
    #
    # Read the soil texture data base and retrieve the various texture-dependent
    # soil parameters for the texture class 'tClass'
    # Refer function readSoilDatabase() in VBA 
    # @param tClass type(str) Soil class id
    #
    def readSoilDatabase(self, tClass):
        
        for key in self.soilDB.keys():
            if key == tClass:                
                self.SWsat  = self.soilDB[key]['SWsat']              
                self.SWfcap = self.soilDB[key]['SWfcap']
                self.SWwilt = self.soilDB[key]['SWwilt']
                self.cTheta = self.soilDB[key]['cTheta']
                self.nTheta = self.soilDB[key]['nTheta']
                self.ESoil1 = self.soilDB[key]['Esoil1']
                self.ESoil2 = self.soilDB[key]['Esoil2']
                self.kDrain = self.soilDB[key]['kDrain']
                self.kSCond = self.soilDB[key]['kSCond']               
 
            
    #
    # Get "Stand data" block parameters
    #
    def readStandData(self):
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        self.locateCell("Species =", cellLocation)
        self.species = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Date planted =", cellLocation)        
        self.PlantedDate  = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value.split('/')
        self.PlantedYear  = int(self.PlantedDate[0])
        self.PlantedMonth = int(self.PlantedDate[1])
        
        # Attention: 0.5(0 year, May), 0.12(0 year, December)
        self.locateCell("Initial age =", cellLocation)        
        self.InitialAge   = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value.split('.')
        self.InitialYear  = self.PlantedYear + int(self.InitialAge[0])
        self.InitialMonth = self.PlantedMonth + int(self.InitialAge[1])
        if self.InitialMonth >= 12:
            self.InitialYear += 1
            self.InitialMonth -= 12
        
        self.locateCell("End age =", cellLocation)
        self.endAge = float(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
         
        # Judge whether this planting is seedling since the simulation
        if self.PlantedYear == self.InitialYear and self.PlantedMonth == self.InitialMonth:
            self.isSeedlings = True 
        else:     
            self.isSeedlings = False
         
        # Get initial biomass data, refer getInitialBiomassData() function in VBA   
        if self.isSeedlings:
            if self.locateCell("Seedling mass", cellLocation):
                self.SeedlingMass = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            else:
                self.SeedlingMass = 0
                if globalVariablesClass.ERROR_MESSAGE_SHOW:
                    print ("Input error, SeedlingMass = %s is not a valid seedling mass" % self.SeedlingMass)
      
        else:
            if self.locateCell("Stand mass", cellLocation):
                self.StandMass = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            else:
                self.StandMass = 0
                self.locateCell("Initial WF =", cellLocation)
                self.WFi = float(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
                self.locateCell("Initial WS =", cellLocation)
                self.WSi = float(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
                self.locateCell("Initial WR =", cellLocation)
                self.WRi = float(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
                
        # Reads fraction of total mass in a biomas pool, refer getBiomassFractions() functoin in VBA
        if self.locateCell("WF fraction", cellLocation):
            x = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        else:
            x = 0
            
        if x > 0:
            self.fracWF = x
            
            if self.locateCell("WR fraction", cellLocation):
                self.fracWR = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            else:
                self.fracWR = 0
            
            if self.locateCell("WS fraction", cellLocation):
                self.fracWS = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            else:
                self.fracWS = 0
        
        self.locateCell("Initial stocking =", cellLocation)
        self.StemNoi = float(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
        
        self.locateCell("Edge trees =", cellLocation)
        self.edgeTrees = int(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
        
        self.locateCell("Initial ASW =", cellLocation)
        self.SWi = float(self.sh.cell(cellLocation['row'], cellLocation['col']+1).value)
        
        self.locateCell("%Tree cover =", cellLocation)
        self.treeCover = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value / 100.0
 
              
    #
    # Get "Soil data" block parameters, refer function "readSoilData()" in VBA
    #
    def readSoilData(self):
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        if self.isSoilProfile:
            if self.locateCell("SWsat", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
                self.SWsat  = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            if self.locateCell("SWfcap", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
                self.SWfcap = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            if self.locateCell("SWwilt", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
                self.SWwilt = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            if self.locateCell("kDrain", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
                self.kDrain = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
            if self.locateCell("kSCond", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
                self.kSCond = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        if self.locateCell("cTheta", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
            self.cTheta  = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        if self.locateCell("nTheta", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
            self.nTheta  = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        if self.locateCell("Esoil1", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
            self.ESoil1  = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        if self.locateCell("Esoil2", cellLocation) and self.sh.cell(cellLocation['row'], cellLocation['col']+1).value != '':
            self.ESoil2  = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value    


    #
    # Get "Understory data" block parameters, refer function "readUnderstoryData()" in VBA
    #        
    def readUnderstoryData(self):
       
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        # Get the location of the cell with value of "Understory data"
        self.locateCell("xUnderstory data", cellLocation)
        
        # Max gC
        x = self.sh.cell(cellLocation['row']+1, cellLocation['col']+1).value        
        if x != '':
            self.uMaxCond = x            
        
        # Index
        x = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value        
        if x != '':
            self.uIndex = x
            
        # uAge1 & uLAI1
        x = self.sh.cell(cellLocation['row']+4, cellLocation['col']).value        
        if x != '':
            self.uAge1 = x
            
        x = self.sh.cell(cellLocation['row']+4, cellLocation['col']+1).value        
        if x != '':
            self.uLAI1 = x
            
        # uAge2 & uLAI2
        x = self.sh.cell(cellLocation['row']+5, cellLocation['col']).value        
        if x != '':
            self.uAge2 = x
            
        x = self.sh.cell(cellLocation['row']+5, cellLocation['col']+1).value        
        if x != '':
            self.uLAI2 = x
            
        # uAge3 & uLAI3
        x = self.sh.cell(cellLocation['row']+6, cellLocation['col']).value        
        if x != '':
            self.uAge3 = x
            
        x = self.sh.cell(cellLocation['row']+6, cellLocation['col']+1).value        
        if x != '':
            self.uLAI3 = x
        
        
    #
    # Get "Pasture data" block parameters, refer function "readPastureData()" in VBA
    #        
    def readPastureData(self):
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        # Get the location of the cell with value of "Understory data"
        self.locateCell("xPasture data", cellLocation)        
        
        # Max gC
        x = self.sh.cell(cellLocation['row']+1, cellLocation['col']+1).value        
        if x != '':
            self.pMaxCond = x
        
        # Index
        x = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value        
        if x != '':
            self.pIndex = x
            
        # uAge1 & uLAI1
        x = self.sh.cell(cellLocation['row']+4, cellLocation['col']).value        
        if x != '':
            self.pAge1 = x
            
        x = self.sh.cell(cellLocation['row']+4, cellLocation['col']+1).value        
        if x != '':
            self.pLAI1 = x
            
        # uAge2 & uLAI2
        x = self.sh.cell(cellLocation['row']+5, cellLocation['col']).value        
        if x != '':
            self.pAge2 = x
            
        x = self.sh.cell(cellLocation['row']+5, cellLocation['col']+1).value        
        if x != '':
            self.pLAI2 = x
            
        # uAge3 & uLAI3
        x = self.sh.cell(cellLocation['row']+6, cellLocation['col']).value        
        if x != '':
            self.pAge3 = x
            
        x = self.sh.cell(cellLocation['row']+6, cellLocation['col']+1).value        
        if x != '':
            self.pLAI3 = x


    #
    # Get "Daily water balance?" block parameters, this module seems still unused in 3-PG2  
    #        
    def readDailyWaterBalance(self):
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        self.locateCell("Daily water balance", cellLocation)
        x = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        if x == 'n':        
            self.dailyWaterBalance = False
        else:
            self.dailyWaterBalance = True
        
        
    #
    # Get "Climate data" block parameters, refer funtion "getSinglesiteClimateData()" in VBA 
    #        
    def readClimateData(self):
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
         
        self.locateCell("Climate data is", cellLocation)
        x = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        if x != 'Daily' and x == 'Monthly':
            self.isDailyMetdata = False
        else:
            print ("Plot %s climate data attribute error" % self.siteName)
            
        self.locateCell("Met station :", cellLocation)
        self.MetStation = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        self.locateCell("Worksheet :", cellLocation)
        self.MetdataSheet = self.sh.cell(cellLocation['row'], cellLocation['col']+1).value
        
        # Get the meta data for the metStation "self.MetStation" into a dictionary structure
        self.getSingleSiteMetData()
        
    
    #
    # Define metData collection for this site, following with the call to getSingleSiteMetaDataClass
    # This dict will have the following structure:
    #    [[12*7 sequential meta data(numpy array)],
    #     [12*7 sequential meta data(numpy array)],
    #     ...
    #     [12*7 sequential meta data(numpy array)]]
    #
    def getSingleSiteMetData(self):
        
        # Collection of the meta data for this site
        self.metDataCollection = []
        
        # List of the years of the met data for this station
        self.metYears = []        
        
        # Open the sheet(Silo metdata) which contains the certain site parameters        
        self.sh_metaData = self.wb.sheet_by_name(self.MetdataSheet)  
        
        self.getClimateVarNames()        
        
        # Location dictionary contains the [row, column] position of a certain cell
        rowRange = {'start': 0, 'end': 0}
        
        # Get the range of the meta records for this station
        self.locateMetaStationCell(self.MetStation, rowRange)

        # Search the 1st column of the "self.MetdataSheet" excel sheet from top to bottom,
        # match the cell value equal "self.MetStation", then store the following 12*n 
        # meta variables into one list (metaVarList)      
        for row in range(rowRange['start'], rowRange['end'] + 1):
            
            year = int(self.sh_metaData.cell(row, 1).value)
            
            if year >= self.InitialYear and year <= self.InitialYear + self.endAge:
            
                self.metYears.append(int(self.sh_metaData.cell(row, 1).value))
            
                # Get 12*n meta variables, n is the climate variables' number
                metaVarList = []
            
                # len(self.climateVarNames)-1)*12 is the total number of climate data in one line
                for col in range(0, (len(self.climateVarNames)-1)*12):
                    metaVarList.append(self.sh_metaData.cell(row, col+2).value)
                      
                self.metDataCollection.append(metaVarList) 
        
        # Total years climate records of this met station
        self.mYears = len(self.metYears)
         
        self.setMetDataArrays()
        
        self.computeMetdataVbls()
        
        
    #
    # Compute metdata variables not supplied as inputs
    # Refer function computeMetdataVbls() in VBA 
    #
    def computeMetdataVbls(self):
        
        self.mDayLength = np.zeros(self.mYears*12)
        
        for month in range(0, self.mYears*12):
            m = month % 12 
            self.mDayLength[month] = 86400 * self.getDayLength(globalVariablesClass.midMonthDay[m])
            
            if self.computeTav:
                self.mTav[month] = (self.mTn[month] + self.mTx[month]) / 2
            else:
                self.mTx[month] = self.mTav[month]
                self.mTn[month] = self.mTav[month]
                                
            if self.computeVPD:
                self.mVPD[month] = self.getVPD(self.mTx[month], self.mTn[month])
 
               
    #
    # Gets daily "mean" VPD in mBar - based on daily max and min temperatures
    # Refer function getVPD() in VBA
    #
    def getVPD(self, Tx, Tn):
        
        VPDx = 6.1078 * np.exp(17.269 * Tx / (237.3 + Tx))
        VPDn = 6.1078 * np.exp(17.269 * Tn / (237.3 + Tn))
    
        return (VPDx - VPDn) / 2
    
    
    #
    # gets fraction of day when sun is "up"
    # Refer function getDayLength() in VBA
    #
    def getDayLength(self, dayOfYear):
        
        SLAt = np.sin(np.pi * self.siteLat / 180)
        cLat = np.cos(np.pi * self.siteLat / 180)
    
        sinDec = 0.4 * np.sin(0.0172 * (dayOfYear - 80))        
        cosH0 = -sinDec * SLAt / (cLat * np.sqrt(1 - (sinDec) ** 2))
                
        if cosH0 > 1:
            return 0
        elif cosH0 < -1:
            return 1
        else:
            return np.arccos(cosH0) / np.pi

        
    #
    # Using the members in self.climateVarNames list(except 'YEAR'), construct corresponding 
    # arrays (1-dimension) and fill met data retrieved from self.metDataCollection
    #
    def setMetDataArrays(self):
            
        # Initialize cliamte variables    
        self.mTx        = np.zeros(len(self.metDataCollection)*12)
        self.mTn        = np.zeros(len(self.metDataCollection)*12)
        self.mRain      = np.zeros(len(self.metDataCollection)*12)
        self.mEpan      = np.zeros(len(self.metDataCollection)*12)
        self.mRad       = np.zeros(len(self.metDataCollection)*12)
        self.mRainDays  = np.zeros(len(self.metDataCollection)*12)
        self.mFrostDays = np.zeros(len(self.metDataCollection)*12)
        self.mTav       = np.zeros(len(self.metDataCollection)*12)
        self.mVPD       = np.zeros(len(self.metDataCollection)*12)
        
        for climateVar in self.climateVarNames:
            if "TMAX" == climateVar:
                self.fillArray(self.mTx, self.climateVarNames.index(climateVar))                
            if "TMIN" == climateVar:
                self.fillArray(self.mTn, self.climateVarNames.index(climateVar))                
            if "RAIN" == climateVar:
                self.fillArray(self.mRain, self.climateVarNames.index(climateVar))
            if "EVAP" == climateVar:
                self.fillArray(self.mEpan, self.climateVarNames.index(climateVar))
            if "RADTN" == climateVar:
                self.fillArray(self.mRad, self.climateVarNames.index(climateVar))
            if "RAIN DAYS" == climateVar:
                self.fillArray(self.mRainDays, self.climateVarNames.index(climateVar))
            if "FROST DAYS" == climateVar:
                self.fillArray(self.mFrostDays, self.climateVarNames.index(climateVar))
            if "TAV" == climateVar:
                self.fillArray(self.mTav, self.climateVarNames.index(climateVar))
            if "VPD" == climateVar:
                self.fillArray(self.mVPD, self.climateVarNames.index(climateVar))
      
        
    #
    # Fill climate array using values retrieved from self.metDataCollection
    # @param array type(np.array) The climate variable to be filled
    # @param index type(int) Indicate the position of the climate variable in 
    #    self.climateVarNames list (Year Tmax Tmin Rain Evap Radtn Rain days  Frost days)
    #
    def fillArray(self ,array, index):
        
        # The position interval of the climate variable        
        index_range = [(index-1)*12, index*12]
        
        # Array index (1-dimension)
        j = 0
        
        for numYear in range(0, len(self.metYears)):
            for i in range(index_range[0], index_range[1]):                
                array[j] = self.metDataCollection[numYear][i]
                j += 1
        

    #
    # Scan the row after "Climate data" cell in "Silo metdata" sheet to 
    # read the metdata input variables' names:
    # e.g. "Year Tmax Tmin Rain Evap Radtn Rain days  Frost days"
    #
    def getClimateVarNames(self):
        
        # Define the string to be found
        strFinder = "Climate data:"
        
        # Climate variable names in "Silo metdata"
        self.climateVarNames = []
        
        self.hasYear    = False
        self.computeVPD = True
        self.computeTav = True
        
        row_position, col_position = 0, 0
        
        # Locate the strFinder
        for row in range(0, self.sh_metaData.nrows):
            for col in range(0, self.sh_metaData.ncols):                
                if strFinder in str(self.sh_metaData.cell(row, col).value):
                    row_position, col_position = row, col
                    
        # The cells after  "Climate data:" and "Year" are the variables we want to keep, so plus 2 
        for offset in range(1, self.sh_metaData.ncols):
            varName = str(self.sh_metaData.cell(row_position, col_position+offset).value).upper()
            if varName != '':
                self.climateVarNames.append(varName)
        
        if 'YEAR' in self.climateVarNames:
            self.hasYear = True
        if 'VPD' in self.climateVarNames:
            self.computeVPD = False
        if 'TAV' in self.climateVarNames:
            self.computeTav = False
        
            
    #
    # Get "Silvicultural events" block parameters, refer funtion "readSilviculturalEvents()" in VBA  
    #        
    def readSilvicultualEvents(self):
        
        # Initialize number of envents
        self.nFR     = 0
        self.nMinASW = 0
        self.nIrrig  = 0
        self.nThin   = 0
        self.nDfol   = 0
                
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        # Check whether the "Silvicultural events" tag exists, if not, exist this function
        if not self.locateCell("Silvicultural events", cellLocation):
            return 
        
        # Read fertility block        
        self.FRages = []
        self.FRVals = []
        
        if self.locateCell("Fertility", cellLocation) or (self.locateCell("Fertility/", cellLocation) and globalVariablesClass.FertilityState):         
            x  = self.sh.cell(cellLocation['row']+2, cellLocation['col']).value # below 'Age'
            x1 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value # below 'FR'
            while x != '':
                self.nFR += 1
                
                # Add 'Age' values into self.FRages 
                self.FRages.append(float(x))            
                x = self.sh.cell(cellLocation['row']+2+self.nFR, cellLocation['col']).value
                
                # Add 'FR' values into self.FRVals
                self.FRVals.append(float(x1))
                x1 = self.sh.cell(cellLocation['row']+2+self.nFR, cellLocation['col']+1).value
                
        # Read irrigation block        
        self.irrigAges = []
        self.irrigVals = []
        
        if self.locateCell("Irrigation", cellLocation) or (self.locateCell("Irrigation/", cellLocation) and globalVariablesClass.IrrigationState):        
            x  = self.sh.cell(cellLocation['row']+2, cellLocation['col']).value # below 'Age'
            x1 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value # below 'FR'
            while x != '':
                self.nIrrig += 1
                
                # Add 'Age' values into self.irrigAges 
                self.irrigAges.append(float(x))            
                x = self.sh.cell(cellLocation['row']+2+self.nIrrig, cellLocation['col']).value
                
                # Add 'FR' values into self.irrigVals
                self.irrigVals.append(float(x1))
                x1 = self.sh.cell(cellLocation['row']+2+self.nIrrig, cellLocation['col']+1).value
                
        # Read thinning block
        self.locateCell("Thinning", cellLocation)
        self.thinAges = []
        self.thinVals = []
        self.thinWF   = []
        self.thinWR   = []
        self.thinWS   = []
        
        x  = self.sh.cell(cellLocation['row']+2, cellLocation['col']).value # below 'Age'
        x1 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value # below 'Stocking'
        x2 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+2).value # below 'F'
        x3 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+3).value # below 'R'
        x4 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+4).value # below 'S'

        while x != '':
            self.nThin += 1
            
            # Add 'Age' values into self.thinAges 
            self.thinAges.append(float(x))            
            x = self.sh.cell(cellLocation['row']+2+self.nThin, cellLocation['col']).value
            
            # Add 'Stocking' values into self.thinVals
            self.thinVals.append(float(x1))
            x1 = self.sh.cell(cellLocation['row']+2+self.nThin, cellLocation['col']+1).value
            
            # Add 'F' values into self.thinWF
            self.thinWF.append(float(x2))
            x2 = self.sh.cell(cellLocation['row']+2+self.nThin, cellLocation['col']+2).value
            
            # Add 'R' values into self.thinWR
            self.thinWR.append(float(x3))
            x3 = self.sh.cell(cellLocation['row']+2+self.nThin, cellLocation['col']+3).value
            
            # Add 'S' values into self.thinWS
            self.thinWS.append(float(x4))
            x4 = self.sh.cell(cellLocation['row']+2+self.nThin, cellLocation['col']+4).value   
            
        # Percentage correction
        self.correctFractions(self.thinWF)
        self.correctFractions(self.thinWR)    
        self.correctFractions(self.thinWS)     
        #print self.thinAges, self.thinVals, self.thinWF, self.thinWR, self.thinWS 
        
        # Read MinASW block
        if self.locateCell("minASW", cellLocation):
            self.minASWages = []
            self.minASWVals = []
            
            x  = self.sh.cell(cellLocation['row']+2, cellLocation['col']).value # below 'Age'
            x1 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value # below 'minASW value'
        
            while x != '':
                self.nMinASW += 1
            
                # Add 'Age' values into self.minASWages 
                self.minASWages.append(float(x))            
                x = self.sh.cell(cellLocation['row']+2+self.nMinASW, cellLocation['col']).value
                
                # Add 'minASW' values into self.minASWVals
                self.minASWVals.append(float(x1))
                x1 = self.sh.cell(cellLocation['row']+2+self.nMinASW, cellLocation['col']+1).value
        
        # Read defoliation block
        if self.locateCell("Defoliation", cellLocation):
            self.defolAges = []
            self.defolVals = []
        
            x  = self.sh.cell(cellLocation['row']+2, cellLocation['col']).value # below 'Age'
            x1 = self.sh.cell(cellLocation['row']+2, cellLocation['col']+1).value # below 'Defoliation value'
        
            while x != '':
                self.nDfol += 1
            
                # Add 'Age' values into self.defolAges 
                self.defolAges.append(float(x))            
                x = self.sh.cell(cellLocation['row']+2+self.nDfol, cellLocation['col']).value
                
                # Add 'Defoliation' values into self.defolVals
                self.defolVals.append(float(x1))
                x1 = self.sh.cell(cellLocation['row']+2+self.nDfol, cellLocation['col']+1).value
        
            self.correctFractions(self.defolVals)
     
            
    #
    # Get "Observed data" block parameters 
    # Notice: Observed data are species specific, so this function may be changed due to different species and data 
    # Planting: Clad/Mac  
    # Observed parameters(the others can be inferred from these parameters): 
    #     in Excel                           in Coding
    #     @Stand age                         standAge
    #     @Initial stocking(stems/ha)        initialStocking
    #     @Stocking(stems/ha)                stocking
    #     @Height(m)                         height
    #     @DBH(cm)                           DBH
    #     @BA(m2/ha)                         BA
    #          
    def readObservedData(self):
        
        # Number of observed data record
        self.observedDataNum = 0
        
        # Define the collection of all the observed data using a dictionary structure
        self.observedData = {}
        
        standAge, initialStocking, stocking, height, DBH, BA = [], [], [], [], [], []
        
        # Location dictionary contains the [row, column] position of a certain cell
        cellLocation = {'row': 0, 'col': 0}
        
        self.locateCell("Observed data", cellLocation)
        # add 3 to arrive the actual position below "Stand age"
        x = self.sh.cell(cellLocation['row']+3, cellLocation['col']).value 
        
        while x != '':
            standAge.append(x)
            self.observedDataNum += 1
            x = self.sh.cell(cellLocation['row']+3+self.observedDataNum, cellLocation['col']).value
        
        for columnIndex in range(0, self.sh.ncols):
            # Get the observed paramters' name
            paraName = self.sh.cell(cellLocation['row']+2, cellLocation['col'] + columnIndex).value
            
            if "Initial stocking" in paraName:
                for rowIndex in range(0, self.observedDataNum):
                    x = self.sh.cell(cellLocation['row']+rowIndex+3, cellLocation['col']+columnIndex).value
                    initialStocking.append(x)
        
            if "Stocking" in paraName:
                for rowIndex in range(0, self.observedDataNum):
                    x = self.sh.cell(cellLocation['row']+rowIndex+3, cellLocation['col']+columnIndex).value
                    stocking.append(x)
        
            if "Height" in paraName:
                for rowIndex in range(0, self.observedDataNum):
                    x = self.sh.cell(cellLocation['row']+rowIndex+3, cellLocation['col']+columnIndex).value
                    height.append(x)
        
            if "DBH" in paraName:
                for rowIndex in range(0, self.observedDataNum):
                    x = self.sh.cell(cellLocation['row']+rowIndex+3, cellLocation['col']+columnIndex).value
                    DBH.append(x)
        
            if "BA" in paraName:
                for rowIndex in range(0, self.observedDataNum):
                    x = self.sh.cell(cellLocation['row']+rowIndex+3, cellLocation['col']+columnIndex).value
                    BA.append(x)
        
        # Add the 6 parameter list into dictionary self.observedData
        self.observedData['standAge'] = standAge
        self.observedData['initialStocking'] = initialStocking
        self.observedData['stocking'] = stocking
        self.observedData['height'] = height
        self.observedData['DBH'] = DBH
        self.observedData['BA'] = BA
        
        
    # Make sure some variable are in the correct range. If not, change them into percentage numbers
    # Refer funciton "correctFractions()" in VBA
    # @param varList type(list) A list of variables
    def correctFractions(self, varList):
        
        isPercentage = True
        
        for i in range(0, len(varList)):
            if varList[i] > 0 and varList[i] < 2:
                isPercentage = False
                
        if isPercentage:
            for i in range(0, len(varList)):
                varList[i] = varList[i] / 100.0
        
      
    #
    # Iterate across the range of the excel sheet ("self.sh") to locate the cell 
    # which contains the "searchStr", and return the location of the cell through "cellLocation"
    # @param searchStr type(str) the cell value to be located
    # @param cellLocation type(dictionary) return the location of the cell with the value of 'searchStr'
    #    
    def locateCell(self, searchStr, cellLocation):
        
        # Reset cellLocation to (0,0)
        cellLocation['row'] = 0
        cellLocation['col'] = 0
        
        # Iterate across the range of the excel sheet ("self.sh") to locate the cell 
        # which contains the "searchStr", and return the location of the cell through "cellLocation"
        for row in range(0, self.sh.nrows):
            for col in range(0, self.sh.ncols):      
                cellVal = str(self.sh.cell(row, col).value).upper()
                cellVal = cellVal.strip()   
                if searchStr.upper() == cellVal:
                    cellLocation['row'] = row
                    cellLocation['col'] = col                                                  
        
        if cellLocation['row'] == 0 and cellLocation['col'] == 0:
            return False
        else:
            return True
       
                
    #
    # Iterate across the 1st column of the metaData sheet, and return the range of the 
    # meta records for a certain station
    # @param stationName type(str) Station name to be located
    # @param rowRange type(dict) Range of the record lines of the meta station 
    #    
    def locateMetaStationCell(self, stationName, rowRange):
        
        # Reset rowRange to (0,0)
        rowRange['start'], rowRange['end'] = 0, 0
                 
        # Define flags 
        startFind = True
        
        # Iterate across the range(1st column) of the excel sheet ("Silo metdata") to locate 
        # the cells which contain the stationName, and return the location of the range through "rowRange"
        for row in range(0, self.sh_metaData.nrows):                               
            if stationName.upper() in str(self.sh_metaData.cell(row, 0).value).upper() and startFind:
                rowRange['start'] = row
                rowRange['end'] = rowRange['start'] - 1
                startFind = False
            
            if stationName.upper() in str(self.sh_metaData.cell(row, 0).value).upper():
                rowRange['end'] += 1 
     
        if rowRange['start'] <= rowRange['end'] == 0:
            print ("Meta station: %s data error" % stationName)
            return False
        else:
            return True
    

    #
    # Calculate maximum available soil water content: depth in m, SWfcap in m3/m3
    # Refer function calculateSoilParameters() in VBA
    # 
    def calculateSoilParameters(self):
        
        # Assign the texture-dependent soil parameters for this soil class
        if self.isSoilProfile:
            self.maxASW = 1000 * (1 - self.pStones / 100) * self.soilDepth * (self.SWfcap - self.SWwilt)
            self.maxSW  = 1000 * (1 - self.pStones / 100) * self.soilDepth * self.SWfcap

        
    #
    # Give all 3PGxl inputs and control data built-in default values, refer Initialise3PGxlInput() in VBA
    #
    def Initialise3PGxlInput(self):
        
        self.siteName = globalVariablesClass.Unknown
        self.MetStation = globalVariablesClass.Unknown
        self.MetdataSheet = globalVariablesClass.Unknown
        self.speciesName = globalVariablesClass.Unknown
        self.mYears = 0
        self.nFR = 0
        self.nMinASW = 0
        self.nIrrig = 0
        self.nThin = 0
        self.nDfol = 0
        self.nVary = 0
        self.PlantedDate = ""
        self.PlantedYear = 0
        self.PlantedMonth = 0
        self.InitialDate = ""
        self.InitialYear = 0
        self.InitialMonth = 0
        self.InitialAge = ""
        self.endAge = 0
        self.opFrequency = 0
        self.clearOutput = False
        self.isDailyMetdata = False
 
            
    #
    # Refer the function initialiseSiteData() in VBA
    #
    def initialiseSiteData(self):
      
        # Initialise the site data
        self.siteLat = globalVariablesClass.MissingValue
        self.FR = globalVariablesClass.MissingValue
        self.EC = 0
        self.soilDepth = globalVariablesClass.MissingValue
        self.wtblDepth = globalVariablesClass.MissingValue
        self.pStones = 0
        self.CO2 = 370
        self.rainIntensity = 99999
        self.treeCover = 1
        self.maxASW = globalVariablesClass.MissingValue
        self.minASW = 0
          
        # Initiliase texture-related soil factors
        self.SWsat  = globalVariablesClass.MissingValue
        self.SWfcap = globalVariablesClass.MissingValue
        self.SWwilt = globalVariablesClass.MissingValue
        self.kSCond = globalVariablesClass.MissingValue
        self.kDrain = globalVariablesClass.MissingValue
        self.cTheta = globalVariablesClass.MissingValue
        self.nTheta = globalVariablesClass.MissingValue
        self.ESoil1 = globalVariablesClass.MissingValue
        self.ESoil2 = globalVariablesClass.MissingValue
        self.soilClass = "?"


    #
    # Initialise the stand to nonsense: these values are checked before 3-PG is run
    # to ensure all initialisations have been completed. Refer function initialiseStandData() in VBA 
    #
    def initialiseStandData(self):
        
        self.startAge = globalVariablesClass.MissingValue
        self.SeedlingMass = 0
        self.StandMass = 0
        self.fracWF = 0.5
        self.fracWR = 0.25
        self.fracWS = 0.25
        self.WFi = globalVariablesClass.MissingValue
        self.WRi = globalVariablesClass.MissingValue
        self.WSi = globalVariablesClass.MissingValue
        self.SWi = 1E+99
        self.StemNoi = globalVariablesClass.MissingValue
        self.edgeTrees = 0
        self.uLAI = 0
        self.pLAI = 0


    #
    # Default parameter values - for E. globulus, taken from Sands & Landsberg (2002)
    # Refer function assignDefaultParameters() in VBA
    #
    def assignDefaultParameters(self):

        self.MaxAge = 50          #Determines rate of "physiological decline" of forest
        self.nAge = 4             #Parameters in age-modifier
        self.rAge = 0.95
        self.fullCanAge = 0       #Age at full canopy cover
        self.k = 0.5              #Radiation extinction coefficient
        self.cVPD0 = 5            #LAI at which VPD in canopy reduced to 50%
        self.gammaR = 0.015       #Root turnover rate per month
        self.tWaterMax = 0.25     #Max thickness of water retained on leaf
        self.MaxIntcptn = 0.15    #Max proportion of rainfall intercepted by canopy
        self.LAImaxIntcptn = 0    #LAI required for maximum rainfall interception
        self.MaxCond = 0.02       #Maximum canopy conductance (gc, m/s)
        self.LAIgcx = 3.33        #LAI required for maximum canopy conductance
        self.gAc = 0.2            #Canopy aerodynamic conductance, assumed constant
        self.gAs = 0.2            #Soil aerodynamic conductance, assumed constant
        self.CoeffCond = 0.05     #Determines response of canopy conductance to VPD
        self.y = 0.47             #Assimilate use efficiency
        self.Tmax = 40            #"Critical" biological temperatures: max, min
        self.Tmin = 8.5           #  and optimum. Reset if necessary/appropriate
        self.Topt = 16
        self.kF = 1               #Number of days production lost per frost days
        self.pFS2 = 1             #Foliage:stem partitioning ratios for D = 2cm
        self.pFS20 = 0.15         #  and D = 20cm
        self.aWs = 0.095          #Stem allometric parameters
        self.nWs = 2.4
        self.pRx = 0.8            #maximum root biomass partitioning
        self.pRn = 0.25           #minimum root biomass partitioning
        self.spRootVol = 10       #specific root volume (m3 soil/kgDM)
        self.m0 = 0               #Value of m when FR = 0
        self.fN0 = 1              #Value of fN when FR = 0
        self.fNn = 0              #Power of (1-FR) in fN
        self.EC0 = 5              #Lower threshold for salinity modifier
        self.EC1 = 15             #Upper threshold for salinity modifier
        self.ECn = 2              #Power for salinity modifier
        self.alphaCx = 0.06       #Canopy quantum efficiency
        self.edgeEffect = 20      #Growth enhancement of edge trees
        self.poolFractn = 0       #Determines fraction of excess water that remains on site
       
        #Stem mortality
        self.gammaN1 = 0          #Coefficients in stem mortality rate
        self.gammaN0 = 0
        self.tgammaN = 2
        self.ngammaN = 1
        self.wSx1000 = 300        #Max tree stem mass (kg) likely in mature stands of 1000 trees/ha
        self.thinPower = 3 / 2    #Power in self-thinning law
        self.mF = 0               #Leaf mortality fraction
        self.mr = 0.2             #Root mortality fraction
        self.mS = 0.2             #Stem mortality fraction
        self.mortality = 0        # ?????????
        
        #Litterfall rate
        self.gammaF1 = 0.027      #Coefficients in monthly litterfall rate
        self.gammaF0 = 0.001
        self.tgammaF = 2
        
        #Specific leaf area
        self.SLA0 = 11            #specific leaf area at age 0 (m^2/kg)
        self.SLA1 = 4             #specific leaf area for mature trees (m^2/kg)
        self.tSLA = 2.5           #stand age (years) for SLA = (SLA0+SLA1)/2
        
        #Branch & bark fraction
        self.fracBB0 = 0.75       #branch & bark fraction at age 0 (m^2/kg)
        self.fracBB1 = 0.15       #branch & bark fraction for mature trees (m^2/kg)
        self.tBB = 2              #stand age (years) for fracBB = (fracBB0+fracBB1)/2
        
        #Basic density
        self.rho0 = 0.45          #basic density for young trees (t/m3)
        self.rho1 = 0.45          #basic density for old trees (t/m3)
        self.tRho = 4             #age at which rho = average of old and young values
        
        #Mean stem height allometric relationship
        self.aH = 0
        self.nHB = 0
        self.nHN = 0
        
        #Stand volume allometric relationship
        self.aV = 0
        self.nVB = 0
        self.nVN = 0
        
        #Spare parameters
        self.param1 = 0
        self.param2 = 0
        self.param3 = 0
        
        #Conversion factors
        self.Qa = -90             #intercept of net v. solar radiation relationship (W/m2)
        self.Qb = 0.8             #slope of net v. solar radiation relationship
        self.gDM_mol = 24         #conversion of mol to gDM
        self.molPAR_MJ = 2.3      #conversion of MJ to PAR
        
        #Understory LAI parameters
        self.uMaxCond = 0         #maximum understory canopy conductance
        self.uIndex = 2           #index for understory model
        self.uAge1 = 0            #age and LAI at first reference point
        self.uLAI1 = 0
        self.uAge2 = 0            #age and LAI at second reference point
        self.uLAI2 = 0
        self.uAge3 = 0            #age and LAI at third reference point
        self.uLAI3 = 0
        
        #Pasture LAI parameters
        self.pMaxCond = 0         #maximum understory canopy conductance
        self.pIndex = 2           #index for understory model
        self.pAge1 = 0            #age and LAI at first reference point
        self.pLAI1 = 0
        self.pAge2 = 0            #age and LAI at second reference point
        self.pLAI2 = 0
        self.pAge3 = 0            #age and LAI at third reference point
        self.pLAI3 = 0
        
        # in function getStandCharacteristics()
        self.oldVol = 0.0
      
        
    #
    # Assign hidden parameters & controls, refer function assignHiddenParameters() in VBA 
    #
    def assignHiddenParameters(self):
        self.poolFractn = 0
        self.param1 = 0
        self.param2 = 0
        self.param3 = 0
       

    #
    # Using true parameters' value to overlap values given by function assignDefaultParameters() 
    #
    def assignTrueParameters(self):
        
        for key in self.paraDict.keys():
            # Allometric relationships & partitioning
            if key == "pFS2":
                self.pFS2 = self.paraDict[key] 
                continue
            if key == "pFS20":
                self.pFS20 = self.paraDict[key] 
                continue            
            if key == "aWS":
                self.aWs = self.paraDict[key] 
                continue
            if key == "nWS":
                self.nWs = self.paraDict[key] 
                continue                        
            if key == "pRx":
                self.pRx = self.paraDict[key] 
                continue 
            if key == "pRn":
                self.pRn = self.paraDict[key] 
                continue
            if key == "spRootVol":
                self.spRootVol = self.paraDict[key] 
                continue
            
            # Litterfall & root turnover
            if key == "gammaF1":
                self.gammaF1 = self.paraDict[key] 
                continue
            if key == "gammaF0":
                self.gammaF0 = self.paraDict[key] 
                continue            
            if key == "tgammaF":
                self.tgammaF = self.paraDict[key] 
                continue
            if key == "gammaR":
                self.gammaR = self.paraDict[key]
                
                continue                        
            
            # Temperature modifier (gmTemp)
            if key == "Tmin":
                self.Tmin = self.paraDict[key] 
                continue            
            if key == "Topt":
                self.Topt = self.paraDict[key] 
                continue
            if key == "Tmax":
                self.Tmax = self.paraDict[key] 
                continue  
                        
            # Frost modifier (gmFrost)
            if key == "kF":
                self.kF = self.paraDict[key] 
                continue  

            # Fertitlity effects (gmNutr)
            if key == "m0":
                self.m0 = self.paraDict[key] 
                continue            
            if key == "fN0":
                self.fN0 = self.paraDict[key] 
                continue
            if key == "fNn":
                self.fNn = self.paraDict[key] 
                continue  

            # Salinity effects (gmSalt)
            if key == "EC0":
                self.EC0 = self.paraDict[key] 
                continue            
            if key == "EC1":
                self.EC1 = self.paraDict[key] 
                continue
            if key == "ECn":
                self.ECn = self.paraDict[key] 
                continue  

            # Age modifier (gmAge)
            if key == "MaxAge":
                self.MaxAge = self.paraDict[key] 
                continue            
            if key == "nAge":
                self.nAge = self.paraDict[key] 
                continue
            if key == "rAge":
                self.rAge = self.paraDict[key] 
                continue  

            # Stem mortality & self-thinning
            if key == "gammaN1":
                self.gammaN1 = self.paraDict[key] 
                continue            
            if key == "gammaN0":
                self.gammaN0 = self.paraDict[key] 
                continue
            if key == "tgammaN":
                self.tgammaN = self.paraDict[key] 
                continue  
            if key == "ngammaN":
                self.ngammaN = self.paraDict[key] 
                continue            
            if key == "wSx1000":
                self.wSx1000 = self.paraDict[key] 
                continue
            if key == "thinPower":
                self.thinPower = self.paraDict[key] 
                continue  
            if key == "mF":
                self.mF = self.paraDict[key] 
                continue            
            if key == "mR":
                self.mR = self.paraDict[key] 
                continue
            if key == "mS":
                self.mS = self.paraDict[key] 
                continue  

            # Specific leaf area
            if key == "SLA0":
                self.SLA0 = self.paraDict[key] 
                continue            
            if key == "SLA1":
                self.SLA1 = self.paraDict[key] 
                continue
            if key == "tSLA":
                self.tSLA = self.paraDict[key] 
                continue  

            # Light interception & VPD attenuation
            if key == "k":
                self.k = self.paraDict[key] 
                continue            
            if key == "fullCanAge":
                self.fullCanAge = self.paraDict[key] 
                continue
            if key == "cVPD0":
                self.cVPD0 = self.paraDict[key] 
                continue  

            # Rainfall interception
            if key == "tWaterMax":
                self.tWaterMax = self.paraDict[key] 
                continue            
            if key == "MaxIntcptn":
                self.MaxIntcptn = self.paraDict[key] 
                continue
            if key == "LAImaxIntcptn":
                self.LAImaxIntcptn = self.paraDict[key] 
                continue  

            # Production and respiration
            if key == "alphaCx":
                self.alphaCx = self.paraDict[key] 
                continue            
            if key == "edgeEffect":
                self.edgeEffect = self.paraDict[key] 
                continue
            if key == "Y":
                self.Y = self.paraDict[key] 
                continue  

            # Conductance
            if key == "MaxCond":
                self.MaxCond = self.paraDict[key] 
                continue            
            if key == "LAIgcx":
                self.LAIgcx = self.paraDict[key] 
                continue
            if key == "CoeffCond":
                self.CoeffCond = self.paraDict[key] 
                continue  
            if key == "gAc":
                self.gAc = self.paraDict[key] 
                continue            
            if key == "gAs":
                self.gAs = self.paraDict[key] 
                continue
            
            # Branch and bark fraction (fracBB)
            if key == "fracBB0":
                self.fracBB0 = self.paraDict[key] 
                continue            
            if key == "fracBB1":
                self.fracBB1 = self.paraDict[key] 
                continue
            if key == "tBB":
                self.tBB = self.paraDict[key] 
                continue  

            # Basic Density
            if key == "rho0":
                self.rho0 = self.paraDict[key] 
                continue            
            if key == "rho1":
                self.rho1 = self.paraDict[key] 
                continue
            if key == "tRho":
                self.tRho = self.paraDict[key] 
                continue  

            # Stem height
            if key == "aH":
                self.aH = self.paraDict[key] 
                continue            
            if key == "nHB":
                self.nHB = self.paraDict[key] 
                continue
            if key == "nHN":
                self.nHN = self.paraDict[key] 
                continue  

            # Stem volume
            if key == "aV":
                self.aV = self.paraDict[key] 
                continue            
            if key == "nVB":
                self.nVB = self.paraDict[key] 
                continue
            if key == "nVN":
                self.nVN = self.paraDict[key] 
                continue  

            # Conversion factors
            if key == "Qa":
                self.Qa = self.paraDict[key] 
                continue            
            if key == "Qb":
                self.Qb = self.paraDict[key] 
                continue
            if key == "gDM_mol":
                self.gDM_mol = self.paraDict[key] 
                continue  
            if key == "molPAR_MJ":
                self.molPAR_MJ = self.paraDict[key] 
                continue   
 
        
# 
# Assemble single site parameters using site series list got from module "getSiteSeries" and, 
# single site parameters got from class "getSingleSiteParametersClass". After execution, one 
# dictionary structure as following will be formed in the memory:
#    [site_1 name: getSingleSiteParametersClass instance 1]
#    [site_2 name: getSingleSiteParametersClass instance 2]
#      ......
#    [site_n name: getSingleSiteParametersClass instance n]
#  
# @param excel_file type(str) Excel file name and path defined in module "globalVariables"
# @param siteSeriesList type(list) One list of site series, refer module "getSiteSeries"
# @param soilDB type(dictionary) Soil database from sheet "3PG_Parameters" - Standard 3PGxl soil types
# @param paraDict type(dictionary) Parameters read from excel sheet "3PG_Parameters" for certain species
# @param sitesParameterCollection type(dictionary) Reference of the global variable "sitesParameterCollection" 
#        defined in module "globalVariables"
#
def assembleSingleSite(excel_file, siteSeriesList, soilDB, paraDict, sitesParameterCollection):
    
    # Iterate the site series
    for site in siteSeriesList:        
        # Using class getSingleSiteParametersClass() instance fill the site specific parameters
        sitesParameterCollection[site] = getSingleSiteParametersClass(soilDB, paraDict, excel_file, site)
        
    




