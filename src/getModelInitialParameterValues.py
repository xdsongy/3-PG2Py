#-------------------------------------------------------------------------------
# Name:        getModelInitialParameterValues
# Purpose:     Extract model initial parameter values from "Standard 3PG parameters" section
#              in excel sheet "3PG_Parameters"  
#
#
# Author:      Xiaodong Song
#
# Created:     08/08/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python


# Module for excel file reading 
import xlrd

# Global variables
import globalVariables

class readInitialParametersClass():
    
    def __init__(self, excel_file_path, excel_sheet):
       
        # Generate global variables
        self.GlobalVar = globalVariables.globalVariablesClass()
       
        # Set the initial excel file path and sheet name
        self.excel_file_path = excel_file_path
        self.excel_sheet =excel_sheet
        
        # Call this function to read and fill initialParametersDict{} 
        self.fillDict()
        
        
    # Fill the dictionary with data read from excel "Standard 3PGxl soil types" block
    def fillDict(self):
        
        # Open the workbook
        wb = xlrd.open_workbook(self.excel_file_path, 'rb')
        
        # Open the sheet which contains the 3PG model initial parameters
        sh = wb.sheet_by_name(self.excel_sheet)
        
        # 3-PG initial parameters dictionary
        self.paraDict = {}
                
        #=======================================================================
        # Set the location / limitation (rows & columns) of the parameters' name and value
        #=======================================================================        
        # paraStartRow: beginning row number of the parameter list
        # paraEndRow: ending row number of the parameter list
        row_start, row_end = self.GlobalVar.paraStartRow, self.GlobalVar.paraEndRow 
                                                                                    
        # paraNameCol: indicate the column number of parameter name
        # paraValueCol: indicate the column of planting 'Clad/Mac'
        col_paraName, col_paraValue = self.GlobalVar.paraNameCol, self.GlobalVar.paraValueCol  
                                                                                        
        
        for row in range(row_start, row_end+1):
            # Read parameters' name and value
            para_name  = sh.cell(row, col_paraName).value
            para_value = sh.cell(row, col_paraValue).value
            
            if len(para_name) > 0:
                self.paraDict[para_name] = para_value
                
