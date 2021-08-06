#===============================================================================
# 
# Name: SA_Cluster_Scheduler
#
# Purpose: Using BootStrap method to        
#
# Created:     08/10/2011
# 
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
# 
#===============================================================================

import shutil, os, subprocess, time, itertools

# Node number of CL*-WRON, here "*" is represented by the nodeNumList  
nodeNumList = ['01','02','03','04','05','06','07','08','09','10','15','16','17','18','19','20','21','22']

# Make sure no network connections
os.system('NET USE * /delete /yes')

def copyData(nodeNum):
    
    # Connect to network drive
    networkPath = '\\\\wron\\cluster\\cl' + nodeNum
    os.system('NET USE ' + networkPath + ' /User:son027 XD@csiro')
    
    # Clean up
    filePath = '/carbon_hotspots/GIS_Data/carbon_viability/viable_summaries_pnas.txt'
    if os.path.exists(networkPath + filePath): os.remove(networkPath + filePath)
    filePath = '/carbon_hotspots/GIS_Data/carbon_viability/error.txt'
    if os.path.exists(networkPath + filePath): os.remove(networkPath + filePath)
    filePath = '/carbon_hotspots/GIS_Data/carbon_viability/output.txt'
    if os.path.exists(networkPath + filePath): os.remove(networkPath + filePath)
    
    # Copy AML file to get latest version
    filePath = '/carbon_hotspots/GIS_Data/carbon_viability/calc_viab_pnas_parallel.aml'
    if os.path.exists(networkPath + filePath): os.remove(networkPath + filePath)
    shutil.copy2('h:' + filePath, networkPath + filePath)
    
    # Copy data
    # if os.path.exists(networkPath + '/carbon_hotspots'): os.remove(networkPath + '/carbon_hotspots')
    # shutil.copytree('h:/carbon_hotspots', networkPath + '/carbon_hotspots')
    
    # Close connection    
    os.system('NET USE ' + networkPath + ' /delete')

# Copy and/or clean up data, note POWERAPP data copied manually to c:\hpc folders
# for nodeNum in nodeNumList:
    # copyData(nodeNum)







if __name__ == '__main__':
    pass