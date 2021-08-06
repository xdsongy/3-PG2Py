#-------------------------------------------------------------------------------
# Name:        getSiteSeries
# Purpose:     Get site series name list from excel sheet "SiteSeries", and store the 
#              sites' names into a list 
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

class getSiteSeriesClass():
    
    def __init__(self, excel_file_path, excel_sheet):
       
        # Set the initial excel file path and sheet name
        self.excel_file_path = excel_file_path
        self.excel_sheet =excel_sheet
        
        # Call this function to read and fill initialParametersDict{} 
        self.fillList()
        
        
    # Fill the dictionary with data read from excel "Standard 3PGxl soil types" block
    def fillList(self):
        
        # Open the workbook
        wb = xlrd.open_workbook(self.excel_file_path, 'rb')
        
        # Open the sheet which contains the 3PG model initial parameters
        sh = wb.sheet_by_name(self.excel_sheet)
        
        # Site series list, contain planting plots' names (same with excel sheets' names)
        self.siteSeries = []
        
        row_num = 0
        
        while len(sh.cell(row_num, 0).value) > 0:                     
            self.siteSeries.append(sh.cell(row_num, 0).value)
            
            # has arrived the limit of the row number of the sheet
            if sh.nrows == row_num + 1:
                break
            else:
                row_num += 1
                
            