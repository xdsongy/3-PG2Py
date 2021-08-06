#-------------------------------------------------------------------------------
# Name:        globalVariables
# Purpose:     collection of all the global variables  
#
# Author:      Xiaodong Song
#
# Created:     08/08/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import numpy as np
import sys

class globalVariablesClass():
    
    # output sqlite file name, here we use sqlite DB to contain the output ot Py3PG2
    if sys.platform == 'win32':
        outputFilename = r"..\Output\Py3PG2Output.sqlite"   # Win32
    else:    
        outputFilename = r"../Output/Py3PG2Output.sqlite"   # Linux, sys.platform = 'linux2'
    
    # in Excel sheet "3PG_Parameters", indicate the start and end row number of the parameters' list
    paraStartRow, paraEndRow = 14, 99
    
    # in Excel sheet "3PG_Parameters", indicate the column numbers of the parameters' name and value (eg. Clad/Mac)
    paraNameCol, paraValueCol = 2, 4
    
    # Define the row and column ranges of soil data block
    soilRowStart, soilRowEnd = 11, 23
    soilColStart, soilColEnd = 17, 29 
    
    # A switch to control whether to show the error or debugging message 
    ERROR_MESSAGE_SHOW = True
    
    # Given missing value
    MissingValue = -999
    
    Unknown = "Unknown"
    
    # Constants for output frequency
    opfNone = 0
    opfRotation = 1
    opfAnnual = 2
    opfMonthly = 3
    
    # Constants for stage in run at which write3PGResults called
    opStart = 0
    opEndMonth = 1            
    opEndRun = 2
    
    # Growth modifer
    use3PG1modifiers = False
    
    #===========================================================================
    # Silvicultural events
    #===========================================================================
    FertilityState  = False # Indicate when we meet "Fertility/", whether or not to read 
                            # the "Age" and "FR" variables below. If there is no back-slash 
                            # hehind "Fertility", always read the following "Age" and "FR" values   
    
    IrrigationState = False # same as above 
  
    # Number of day in each moonth
    daysInMonth = np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31])
    
    # Middle day number of each month
    midMonthDay = np.array([16, 44, 75, 105, 136, 166, 197, 228, 258, 289, 319, 350])

    def setSiteSeries(self, sheetName):

        self.siteSeries = sheetName

    
    def __init__(self):
        
        #=======================================================================
        # Parameter file info (Excel file & sheets) 
        #=======================================================================
        
        # Given the Excel file name and path
        if sys.platform == 'win32':
            self.excel_file = r"..\Data\3PGxl_Parameters.xls"
        else:
            self.excel_file = r"../Data/3PGxl_Parameters.xls"   # Linux, sys.platform = 'linux2'
        
        # Set the 3-PG model initial parameters' sheet name
        self.parameters_3PG = u"3PG_Parameters"  
        
        # Set the "SiteSeries" sheet name
        # self.siteSeries = u"Site"
        # self.siteSeries = u"SiteSeries"
        
        #=======================================================================
        # Standard 3PGxl options (reference Excel file)
        #=======================================================================
        
        # Output frequency (None, Rotation, Annual or Monthly)
        self.outputFrequency = 'Annual'
                
        # Interpolate lookups? (True/False)
        self.interpolateLookups = True
        
        # Output data (3-PG names from 3PGxl_Outputs)
        self.outputData = ['StemNo', 'WF', 'WR', 'WS', 'standVol', 'LAI', 'MAI', 'avDBH']
        
        #=======================================================================
        # Assembly of single sites' parameters as a dictionary 
        # (refer module 'getSingleSiteParameters')
        #=======================================================================
        
        self.sitesParameterCollection = {}