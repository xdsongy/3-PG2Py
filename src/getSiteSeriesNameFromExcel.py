#-------------------------------------------------------------------------------
# Name:        getSiteSeriesNameFromExcel.py
# Purpose:     Get site series name list from excel file "3PGxl_Parameters.xls"
#              
#
#
# Author:      Xiaodong Song
#
# Created:     26/09/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python

import xlrd

excel_file = r"..\Data\3PGxl_Parameters.xls"
output_file = r"..\Output\excel_sheets_namelist.txt"

wb = xlrd.open_workbook(excel_file)

# Open the output txt file to store the excel sheet names
fp = open(output_file, 'w')

sheets = wb._sheet_list

for sheet in sheets:
    fp.write(sheet.name + "\n")


fp.close()

