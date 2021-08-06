#-------------------------------------------------------------------------------
# Name:        getSoilTypeData
# Purpose:     Extract soil type data from the Excel sheets: "3PG_Parameters" into memory,  
#              make "Soil class" column as dictionary keys, the other 12 soil attributes as 
#              the corresponding values. Just like the following structure:
#
#                {'C':  {'Sand': 30, 'Clay': 50, ...},
#                 'CL': {'Sand': 33, 'Clay': 34, ...},
#                 ... }
#
#
#
# Author:      Xiaodong Song
#
# Created:     08/08/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python

#
# Module for excel file reading 
#
import xlrd

# Global variables
import globalVariables

class readSoilTypeInfoClass():
    
    def __init__(self, excel_file_path, excel_sheet):
       
        # Generate global variables
        self.GlobalVar = globalVariables.globalVariablesClass()
       
        # Set the initial excel file path and sheet name
        self.excel_file_path = excel_file_path
        self.excel_sheet =excel_sheet
        
        # Call this function to read and fill soilDict{} 
        self.fillDict()
        
        
    # Fill the dictionary with data read from excel "Standard 3PGxl soil types" block
    def fillDict(self):
        
        # Open the workbook
        wb = xlrd.open_workbook(self.excel_file_path, 'rb')
        
        # Open the sheet which contains the 3PG model initial parameters
        sh = wb.sheet_by_name(self.excel_sheet)
        
        # Define the soil type data dictionary, all soil type info will be stored into this dict
        self.soilDict = {}
        
        # Set the range of soil data block 
        row_start, row_end = self.GlobalVar.soilRowStart, self.GlobalVar.soilRowEnd  # e.g. 11, 23
        col_start, col_end = self.GlobalVar.soilColStart, self.GlobalVar.soilColEnd  # e.g. 17, 29
        
        # define key words for the soil dictionary (excel block header)
        keys = []
        
        for i in range(col_start+1, col_end+1):
            # Append the attributes' name and romove the "% " before the first 3 attributes
            keys.append(sh.cell(row_start, i).value.replace('% ', ''))
        
        for i in range(row_start+1, row_end+1):
            # Number values for each soil type (one line)
            values = []
            
            for j in range(col_start+1, col_end+1):
                values.append(sh.cell(i, j).value)
            
            # Assembly each line into one dictionary
            sub_dict = {}
            
            for index, items in enumerate(keys):
                sub_dict[items] = values[index]
                
            # Push into the self.soilDict
            self.soilDict[sh.cell(i, col_start).value] = sub_dict

