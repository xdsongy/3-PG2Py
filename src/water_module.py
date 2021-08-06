#-    ------------------------------------------------------------------------------
# Name:        water_module
# Purpose:     code for the 3-PG2 Water Balance model  
#
#
#
# Author:      Xiaodong Song
#
# Created:     17/08/2011
# Copyright:   (c) Xiaodong Song <xdsongy@gmail.com> 2011
#-------------------------------------------------------------------------------
#!/usr/bin/env python


#
# Global variables
#
from globalVariables import globalVariablesClass

#
# Retrieve single site parameters 
#
import getSingleSiteParameters

import math

class water_module_class():
    
    def __init__(self, singleSiteInstance):
        
        self.site = getSingleSiteParameters.getSingleSiteParametersClass()
        self.site = singleSiteInstance
        
        self.initialiseSoilWaterBalance()
    
    
    #
    # Initial SW must be between min and max SW
    # Refer function initialiseSoilWaterBalance() in VBA 
    #
    def initialiseSoilWaterBalance(self):
        
        if self.site.minASW > self.site.maxASW:
            self.site.minASW = self.site.maxASW
            
        if self.site.SWi <= self.site.minASW:
            self.site.SWi = self.site.minASW    
        
        if self.site.isSoilProfile:
            if self.site.SWi >= self.site.maxSW:
                self.site.SWi = self.site.maxSW
        else:
            if self.site.SWi >= self.site.maxASW:
                self.site.SWi = self.site.maxASW
        
        self.site.applIrrig = 0
        self.site.fPooledSW = 0
        self.site.pPooledSW = 0
        self.site.fEsum0 = 0
        self.site.fEsum1 = 0
        self.site.pEsum = 0
        self.site.soilVol = self.site.treeCover * self.site.soilDepth
        self.site.rootDepth = self.fnRootDepth(self.site.WR)        
        self.site.rzVol = self.site.treeCover * self.site.rootDepth
        self.site.forestSW = self.site.treeCover * self.site.SWi
        self.site.pastureSW = (1.0 - self.site.treeCover) * self.site.SWi
        self.site.totalSW = self.site.forestSW + self.site.pastureSW
    
        if self.site.isSoilProfile:
            self.site.fASW = self.fnWaterContent(self.site.forestSW, self.site.SWwilt)            
        else:
            self.site.fASW = self.site.forestSW
        
        self.site.rootFrac = self.site.rootDepth / self.site.soilDepth
        self.site.rzSW = self.site.rootFrac * self.site.forestSW
        self.site.nrSW = (1.0 - self.site.rootFrac) * self.site.forestSW
        self.site.rzRelASW = self.fnRelASW(self.site.rzSW, self.site.rzVol)
        self.site.nrRelASW = self.fnRelASW(self.site.nrSW, self.site.soilVol - self.site.rzVol)
        
        
    def fnRootDepth(self, WR):
        
        if self.site.spRootVol == 0:
            z = self.site.soilDepth
        else:
            z = (0.1 * WR * self.site.zEffect * self.site.spRootVol) / (self.site.treeCover * (1 - self.site.pStones / 100.0))
            
        return min(self.site.soilDepth, z)


    #
    # Routines for soil water balance model, returns water content of soil above critSW
    # Refer function fnWaterContent() in VBA
    #
    def fnWaterContent(self, SW, critSW):
        
        if self.site.isSoilProfile:
            vol = max(0, SW - 1000 * (1 - self.site.pStones / 100) * self.site.soilVol * critSW)
        else:
            vol = SW
            
        return vol


    #
    # Get relative available soil water content for a volume vol of soil
    # Refer function fnRelASW() in VBA
    #
    def fnRelASW(self, SW, vol):
        
        if vol == 0 or SW == 0:
            rASW = 0.0
        elif self.site.isSoilProfile:            
            rASW = (0.001 * SW / ((1 - self.site.pStones / 100) * vol) - self.site.SWwilt) / (self.site.SWfcap - self.site.SWwilt)
        else:
            if self.site.treeCover > 0:
                rASW = SW / (self.site.maxASW * self.site.treeCover)
            else:
                rASW = 0.0
                
        return max(0, min(1, rASW))


    #
    # Do the soil water balance for both forest and pasture
    # Daily soil water balance model
    # Note: we didn't implement daily water balance model, maybe later...
    #
    def doDailySoilWaterBalance(self, month, age):
        pass


    #
    # Do the soil water balance for both forest and pasture
    # @param month type(int) month of a year, between 1 to 12
    # @param age type(float) standAge
    #
    def doMonthlySoilWaterBalance(self, month, age):
        
        self.site.mDays = globalVariablesClass.daysInMonth[month-1]

        self.initialiseMonthlySums()
        
        # Get components of SW balance
        self.getInterceptedRadiation(self.site.SolarRad)
        
        self.getWithinCanopyVPDs(self.site.VPD)
        
        self.getConductances(self.site.VPD, self.site.rzRelASW)
        
        self.getRainfallInterception(self.site.Rain, self.site.RainDays)
        
        self.getTranspiration(self.site.mDays)
        
        self.getSoilEvaporation(self.site.Rain, int(self.site.RainDays), self.site.mDays / self.site.RainDays)
        
        self.getMaxSWcontents()

        # Do the forest & pasture SW balance
        self.doForestWaterBalance()
        self.doPastureWaterBalance()
                
        # Tidy up loose ends
        self.doWatertableAccess()
        self.getPlotSWdetails()


    def initialiseMonthlySums(self):
        
        self.site.RainInt = 0
        self.site.RunOff = 0
        self.site.Drainage = 0
        self.site.Transp = 0
        self.site.Esoil = 0
        self.site.rzRecharge = 0


    #
    # Get radiation intercepted each component
    #
    def getInterceptedRadiation(self, Rad):
        
        # ... forest canopy
        self.getCanopyRadiationInterception(Rad)

        # ... fractional light interceptions
        self.site.uFracInt0 = self.fnFractionIntercepted(self.site.k, self.site.uLAI)
        self.site.uFracInt1 = self.site.uFracInt0
        self.site.pFracInt = self.fnFractionIntercepted(self.site.k, self.site.pLAI)

        # ... intercepted radiation
        self.site.fuRadInt0 = self.site.uFracInt0 * Rad
        self.site.fuRadInt1 = self.site.uFracInt1 * (1 - self.site.fFracInt) * Rad
        self.site.fsRadInt0 = (1 - self.site.uFracInt0) * Rad
        self.site.fsRadInt1 = (1 - self.site.uFracInt1) * (1 - self.site.fFracInt) * Rad
        self.site.pcRadInt = self.site.pFracInt * Rad
        self.site.psRadInt = (1 - self.site.pFracInt) * Rad


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
    # Get VPD at various places in the canopy
    #
    def getWithinCanopyVPDs(self, VPD):
        
        self.site.fcVPD  = VPD
        self.site.fuVPD0 = VPD
        self.site.fuVPD1 = VPD * math.exp(-math.log(2) * self.site.fLAI / self.site.cVPD0)
        self.site.fsVPD0 = VPD * math.exp(-math.log(2) * self.site.uLAI / self.site.cVPD0)
        self.site.fsVPD1 = VPD * math.exp(-math.log(2) * (self.site.fLAI + self.site.uLAI) / self.site.cVPD0)
        self.site.pcVPD  = VPD
        self.site.psVPD  = VPD * math.exp(-math.log(2) * self.site.pLAI / self.site.cVPD0)

    
    #
    # Get the conductances for the current day
    #
    def getConductances(self, VPD, relASW):
        
        # Get the SW modifier as it changes daily with relASW
        self.site.gmVPD = self.growthModifierVPD(VPD)
        self.site.gmSW = self.growthModifierSW(relASW)

        if globalVariablesClass.use3PG1modifiers:
            self.site.physMod = min(self.site.gmVPD, self.site.gmSW) * self.site.gmAge * self.site.gmSalt * self.site.gmCO2gC 
        else:
            self.site.physMod = self.site.gmTemp * self.site.gmVPD * self.site.gmSW * self.site.gmAge * self.site.gmSalt * self.site.gmCO2gC
            
        # get conductances
        self.site.fCond = self.fnCanCond(self.site.MaxCond * self.site.physMod, self.site.fLAI)
        self.site.uCond = self.fnCanCond(self.site.uMaxCond * self.site.physMod, self.site.uLAI)
        self.site.pCond = self.fnCanCond(self.site.pMaxCond, self.site.pLAI)
        

    #
    # Get the canopy conductance for a canopy of LAI = L
    #
    def fnCanCond(self, gCx, L):
        
        return gCx * min(1, L / self.site.LAIgcx)


    #
    # Get the VPD-dependent growth modifier
    #
    def growthModifierVPD(self, VPD):
        
        return math.exp(-self.site.CoeffCond * VPD)


    #
    # Get the soilwater-dependent growth modifier
    #
    def growthModifierSW(self, rASW):
     
        c = self.site.cTheta ** self.site.nTheta
        w = (1 - rASW) ** self.site.nTheta  
        
        return  (1 - w) / (1 + (1 - 2 * c) * w / c)


    #
    # Get rainfall interception by forest and pasture using the leaf water-retention
    # model. Total rainfall is tRain in rEvents rain events
    #
    def getRainfallInterception(self, tRain, rEvents):

        # forest rainfall interception
        e0 = self.fnPenman(self.site.fRadInt, self.site.fcVPD, self.site.DayLength, self.site.gAc, 100000)
        y = self.fnDailyRainInterception(self.site.fLAI + self.site.uLAI, tRain / rEvents, e0)
        
        self.site.fRainInt = rEvents * self.site.treeCover * y
        
        # pasture rainfall interception
        y = 0
        e0 = self.fnPenman(self.site.pcRadInt, self.site.pcVPD, self.site.DayLength, self.site.gAc, 100000)
        y = self.fnDailyRainInterception(self.site.pLAI, tRain / rEvents, e0)
        self.site.pRainInt = rEvents * (1 - self.site.treeCover) * y


    #
    # Apply the Penman-Monteith formula (Landsberg & Gower, 1997) to compute
    # transpiration by the stand in kg/m2/day (i.e. mm/day)
    #
    def fnPenman(self, q, VPD, h, gBL, gC):
        
        e20 = 2.2          # rate of change of saturated VP with T at 20C
        rhoAir = 1.2       # density of air, kg/m3
        lambda_1 = 2460000   # latent heat of vapourisation of H2O (J/kg)
        VPDconv = 0.000622 # VPD ---> saturation deficit = 18/29/1000
        
        netRad = self.site.Qa + self.site.Qb * (q * 10 ** 6 / h)
        Etransp = gC * (e20 * netRad + rhoAir * lambda_1 * (VPDconv * VPD) * gBL) / (gC * (1 + e20) + gBL)
        fnPenman = (Etransp / lambda_1) * h            # now kg/m2/day
        
        return fnPenman 


    #
    # Modelling rainfall interception by canopy
    #
    def fnDailyRainInterception(self, L, Rain, e0):
        
        # Compute daily rainfall interception by a canopy of leaf area index L
        if self.site.rainIntensity == 0:
            e = 0
        else:
            e = e0 / (24 * self.site.rainIntensity)
            
        if Rain <= L * self.site.tWaterMax / (1 - e):
            y = Rain
        else:
            y = L * self.site.tWaterMax + e * Rain
        
        return y


    #
    # Get transpiration accumlated over a period of nDays days at the current
    # ambient climatic conditions
    #
    def getTranspiration(self, nDays):
        
        # forest canopy
        fPenman = self.fnPenman(self.site.fRadInt, self.site.fcVPD, self.site.DayLength, self.site.gAc, self.site.fCond)
        self.site.fTransp = nDays * fPenman * self.site.CanopyCover * self.site.treeCover * self.site.zEffect        
        
        # understory - separate into full sun and shaded components
        uPenman0 = self.fnPenman(self.site.fuRadInt0, self.site.fuVPD0, self.site.DayLength, self.site.gAc, self.site.uCond)
        uPenman1 = self.fnPenman(self.site.fuRadInt1, self.site.fuVPD1, self.site.DayLength, self.site.gAc, self.site.uCond)
        uTransp0 = nDays * uPenman0 * (1 - self.site.CanopyCover) * self.site.treeCover
        uTransp1 = nDays * uPenman1 * self.site.CanopyCover * self.site.treeCover
        self.site.uTransp = uTransp0 + uTransp1
        
        # pasture transpiration
        pPenman = self.fnPenman(self.site.pcRadInt, self.site.pcVPD, self.site.DayLength, self.site.gAc, self.site.pCond)
        self.site.pTransp = nDays * pPenman * (1 - self.site.treeCover)


    #
    # Get the soil evaporation using the Ritchie two-phase model.
    # Total rain fall is tRain in rEvents rainfall events, nDays apart
    #
    def getSoilEvaporation(self, tRain, rEvents, nDays):

        gSoil = 100000
        
        # forest region (separately without and with tree canopy cover)
        self.site.fsPenman0 = self.fnPenman(self.site.fsRadInt0, self.site.fsVPD0, self.site.DayLength, self.site.gAs, gSoil)
        self.site.fsPenman1 = self.fnPenman(self.site.fsRadInt1, self.site.fsVPD1, self.site.DayLength, self.site.gAs, gSoil)
        
        val = self.fnSoilEvap(self.site.fsPenman0, self.site.fEsum0, tRain / rEvents, rEvents, nDays)
        self.site.fEsoil0 = val[0] * self.site.treeCover * (1 - self.site.CanopyCover)
        self.site.fEsum0 = val[1]
        
        val = self.fnSoilEvap(self.site.fsPenman1, self.site.fEsum1, tRain / rEvents, rEvents, nDays)
        self.site.fEsoil1 = val[0] * self.site.treeCover * self.site.CanopyCover
        self.site.fEsum1 = val[1]
        
        self.site.fEsoil = self.site.fEsoil0 + self.site.fEsoil1
        self.site.fEsum = self.site.fEsum0 + self.site.fEsum1
        
        # pasture region
        self.site.psPenman = self.fnPenman(self.site.psRadInt, self.site.psVPD, self.site.DayLength, self.site.gAs, gSoil)
        
        val = self.fnSoilEvap(self.site.psPenman, self.site.pEsum, tRain / rEvents, rEvents, nDays)
        self.site.pEsoil = val[0] * (1 - self.site.treeCover)
        self.site.pEsum = val[1]
        
        
    #
    #
    #
    def fnSoilEvap(self, epot, Esum, Rain, rEvents, nDays):

        y = 0
        if epot > 0:
            for n in range(0, rEvents):
                Esum = max(0, Esum - Rain)
                t0 = self.fnEvapTime(epot, Esum)
                e  = self.fnEvapSum(epot, t0 + nDays) - Esum
                Esum = Esum + e
                y = y + e
                
        # Return a turple, 1st is the sum, 2nd is the changed reference variable (named fEsum0, fEsum1 and pEsum)
        return (y, Esum) 


    #
    # Get time from last wetting event required to evaporate Esum
    #
    def fnEvapTime(self, epot, Esum):
        
        if Esum <= self.site.ESoil1:
            y = Esum / epot
        else:
            tS1 = self.site.ESoil1 / epot
            if self.site.ESoil2 > 0:
                y = tS1 + (self.site.ESoil2 / epot / 2) * ((1 + (Esum - self.site.ESoil1) / self.site.ESoil2) ** 2 - 1)
            else:
                y = tS1
                
        return y

    
    #
    # Get accumatled soil evaporation since last wetting event
    # using the Ritchie two-phase model
    #
    def fnEvapSum(self, epot, t):
        
        tS1 = self.site.ESoil1 / epot
        if t <= tS1:
            y = t * epot
        else:
            if self.site.ESoil2 > 0: 
                y = self.site.ESoil1 + self.site.ESoil2 * (math.sqrt(1 + 2 * (epot / self.site.ESoil2) * (t - tS1)) - 1)
            else:
                y = self.site.ESoil1
                
        return y


    #
    # 
    #
    def getMaxSWcontents(self):
        
        # Get volumetric soil water contents of forest zones
        if self.site.isSoilProfile:
            self.site.rzSWsat = 1000 * (1 - self.site.pStones / 100) * self.site.rzVol * self.site.SWsat
            self.site.nrSWsat = 1000 * (1 - self.site.pStones / 100) * (self.site.soilVol - self.site.rzVol) * self.site.SWsat
            self.site.rzSWfcap = 1000 * (1 - self.site.pStones / 100) * self.site.rzVol * self.site.SWfcap
            self.site.nrSWfcap = 1000 * (1 - self.site.pStones / 100) * (self.site.soilVol - self.site.rzVol) * self.site.SWfcap
        else:
            self.site.rzSWsat = self.site.treeCover * self.site.maxASW
            self.site.nrSWsat = 0
            self.site.rzSWfcap = self.site.rzSWsat
            self.site.nrSWfcap = 0
        
        # Get volumetric soil water content of pasture zone
        if self.site.isSoilProfile:
            self.site.pSWsat = 1000 * (1 - self.site.treeCover) * (1 - self.site.pStones / 100) * self.site.soilDepth * self.site.SWsat
            self.site.pSWfcap = 1000 * (1 - self.site.treeCover) * (1 - self.site.pStones / 100) * self.site.soilDepth * self.site.SWfcap
        else:
            self.site.pSWsat = (1 - self.site.treeCover) * self.site.maxASW
            self.site.pSWfcap = self.site.pSWsat


    #
    # Do the soilwater balance in the forest
    #
    def doForestWaterBalance(self):
        
        # distribute effective rain between the root and root-free zones
        effctvRain = self.site.treeCover * self.site.Rain - self.site.fRainInt
        self.site.rzSW = self.site.rzSW + self.site.rootFrac * effctvRain
        self.site.nrSW = self.site.nrSW + (1 - self.site.rootFrac) * effctvRain
        
        # distribute soil evaporation between the two zones (check SW doesn't go -ve)
        rzEsoil = min(self.site.rzSW, self.site.rootFrac * self.site.fEsoil)
        nrEsoil = min(self.site.nrSW, (1 - self.site.rootFrac) * self.site.fEsoil)
        self.site.rzSW = self.site.rzSW - rzEsoil
        self.site.nrSW = self.site.nrSW - nrEsoil
        self.site.fEsoil = rzEsoil + nrEsoil
        
        # get runoff
        rzSWxs = max(self.site.rzSW - self.site.rzSWsat, 0)
        nrSWxs = max(self.site.nrSW - self.site.nrSWsat, 0)
        self.site.rzSW = self.site.rzSW - rzSWxs
        self.site.nrSW = self.site.nrSW - nrSWxs
        self.site.fRunoff = rzSWxs + nrSWxs
        
        # distribute understory transpiration between soil zones
        rzTransp = self.site.fTransp + self.site.rootFrac * self.site.uTransp
        nrTransp = (1 - self.site.rootFrac) * self.site.uTransp
        self.site.rzSW = self.site.rzSW - rzTransp
        self.site.nrSW = self.site.nrSW - nrTransp
        
        # allow the root zone to wet up
        self.site.rzRechg = self.fnRootZoneFlow(self.site.mDays)            
        self.site.rzSW = self.site.rzSW + self.site.rzRechg
        self.site.nrSW = self.site.nrSW - self.site.rzRechg
        
        # determine drainage
        rzSWxs = max(self.site.rzSW - self.site.rzSWfcap, 0)
        nrSWxs = max(self.site.nrSW - self.site.nrSWfcap, 0)
        kD = 1 - math.exp(-self.site.mDays * self.site.kDrain)
        self.site.rzSW = self.site.rzSW - kD * rzSWxs
        self.site.nrSW = self.site.nrSW - kD * nrSWxs
        self.site.fDrainage = kD * (rzSWxs + nrSWxs)

        # update relative soil water
        self.site.rzRelASW = self.fnRelASW(self.site.rzSW, self.site.rzVol)  
        self.site.nrRelASW = self.fnRelASW(self.site.nrSW, self.site.soilVol - self.site.rzVol)

        # correct forest transpiration for insufficent soil water to support it
        # the same correction factor (fracET) is later applied to NPP
        self.site.fracET = 1
        
        self.site.fTransp = self.site.fracET * self.site.fTransp
        self.site.uTransp = self.site.fracET * self.site.uTransp


    #
    # Modelling water flow between root and non-root zones
    #   
    def fnRootZoneFlow(self, nDays):

        # Analytical solution for amount of water noving over a period nDays days
        # from non-root zone to root zone
        y = 0
        
        if self.site.isSoilProfile and (self.site.soilVol != self.site.rzVol) and (self.site.kSCond != 0):
            nrVol = self.site.soilVol - self.site.rzVol
            tS0 = self.site.rzVol * nrVol / (self.site.kSCond * self.site.soilVol)
            SWa = (self.site.rzSW * nrVol - self.site.nrSW * self.site.rzVol) / self.site.soilVol
            SWb = self.site.rzVol * (self.site.rzSW + self.site.nrSW) / self.site.soilVol
            y = SWa * math.exp(-nDays / tS0) + SWb - self.site.rzSW

        return y


    #
    # Do the soilwater balance for the pasture
    #
    def doPastureWaterBalance(self):

        if self.site.treeCover == 1:
            self.site.pRunoff = 0.0
            self.site.pDrainage = 0.0            
            return
        
        # ... add rain
        effctvRain = (1 - self.site.treeCover) * self.site.Rain - self.site.pRainInt
        self.site.pastureSW += effctvRain

        # ... check for runoff
        self.site.pRunoff = max(self.site.pastureSW - self.site.pSWsat, 0)
        self.site.pastureSW -= self.site.pRunoff
        
        # ... take out soil evaporation
        self.site.pEsoil = min(self.site.pastureSW, self.site.pEsoil)
        self.site.pastureSW -= self.site.pEsoil        
        
        # ... update accumulated evaporation
        self.site.pEsum = max(0, self.site.pEsum + self.site.pEsoil - effctvRain)

        # ... take out pasture transpiration
        self.site.pTransp = min(self.site.pastureSW, self.site.pTransp)
        self.site.pastureSW -= self.site.pTransp        
        
        # ... determine drainage
        self.site.pDrainage = self.site.kDrain * max(self.site.rzSW - self.site.rzSWfcap, 0)
        self.site.pastureSW -= self.site.pDrainage


    #
    # Deal with water table access
    # ... currently this is only relevant to the case where maxASW is specified
    #
    def doWatertableAccess(self):
        
        self.site.SupIrrig = 0
        
        if self.site.isSoilProfile:
            # this has not yet been implemented - would be via capilary rise
            return 
        else:
            rzMinASW = self.site.treeCover * self.site.minASW
            pMinASW = (1 - self.site.treeCover) * self.site.minASW
            fSupIrrig = 0
            if self.site.rzSW < rzMinASW:
                fSupIrrig = rzMinASW - self.site.rzSW
                self.site.rzSW = rzMinASW
                self.site.nrSW = 0
            
            pSupIrrig = 0
            if self.site.pastureSW < pMinASW:
                pSupIrrig = pMinASW - self.site.pastureSW
                self.site.pastureSW = pMinASW
            
            self.site.SupIrrig = fSupIrrig + pSupIrrig

    
    #
    # Update total water pools and flows these are current pools
    #
    def getPlotSWdetails(self):
        
        self.site.forestSW = self.site.rzSW + self.site.nrSW
        self.site.totalSW = self.site.forestSW + self.site.pastureSW
        self.site.fASW = self.fnWaterContent(self.site.forestSW, self.site.SWwilt)
        self.site.pASW = self.fnWaterContent(self.site.pastureSW, self.site.SWwilt)
  
        # these fluxes are updated as running totals to give monthly values
        self.site.RainInt  = self.site.RainInt + self.site.fRainInt + self.site.pRainInt
        self.site.RunOff   = self.site.RunOff + self.site.fRunoff + self.site.pRunoff
        self.site.Drainage = self.site.Drainage + self.site.fDrainage + self.site.pDrainage
        self.site.Transp   = self.site.Transp + self.site.fTransp + self.site.uTransp + self.site.pTransp
        self.site.Esoil    = self.site.Esoil + self.site.fEsoil + self.site.pEsoil
        self.site.rzRecharge = self.site.rzRecharge + self.site.rzRechg

        # all water losses
        self.site.ET = self.site.RainInt + self.site.Esoil + self.site.Transp










