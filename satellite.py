'''

Satellite Object Nameing / Iter convention:

Satellite.SubsetParams['Param'] - returns value list

Satellite.SubsetParams['Param'][0] - returns raw value
Satellite.SubsetParams['Param'][1] - returns unit of value

Subset List:

init: all parameters initially given by satellte developers; from mission requirements
orbital: all parameters relating to the orbital mechanics
optical: all parameters relating to the optical payload (payload subsets may exist w/ multiple payloads eg. optical, IR)
adcs: all parameters relating to the adcs



'''

import numpy as np
import scipy.constants
from IPython import embed

class Satellite(object):

    def __init__(self, initData, ):
        setattr(self,'initParams',initData)
        setattr(self,'subsetList', initData.keys())
        for subset in initData:
            setattr(self, subset, initData[subset])
        self.earthR = [6378.14,'km'] #[m] - Earth Radius using Spherical Model
        self.sunIncidentAngleDeg = [23.5, 'deg']
        self.earthGravConst = 3.986e14 #Earth-Gravity Constant in m3/s2
        self.lightSpeed = 299792458.0 # speed of light in [m/s]
        self.earthMagneticMoment = 7.96e15 #tesla/m^3

        '''
        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])
        '''
    def set(self, *initial_data, **kwargs):
        for dictionary in initial_data:
            for key in dictionary:
                setattr(self, key, dictionary[key])
        for key in kwargs:
            setattr(self, key, kwargs[key])
        return self

    def get(self,param):
        #TODO: Fix for case of multiple param names in different subsets
        for subset in self.subsetList:
            try:
                value = getattr(self,subset)[param][0] 
                return value
            except:
                pass
        return None
        #quickly get param raw value from name
    
    def getU(self,param):
        #TODO: Fix for case of multiple param names in different subsets
        for subset in self.subsetList:
            try:
                value = getattr(self,subset)[param][1] 
                return value
            except:
                pass
        return None
        #quickly get param raw value from name

    def merge_dicts(self, x, y):
        z = x.copy()   # start with x's keys and values
        z.update(y)    # modifies z with y's keys and values & returns None
        return z

    def setSubset(self,subsetName,subsetData):
        try:
            setattr(self, subsetName, self.merge_dicts(getattr(self,subsetName),subsetData))
        except AttributeError:
            setattr(self, subsetName, subsetData)
        return self

    def calculateOpticalParameters(self):
        self.calculateSensorViewingParams()
        self.calculatePixelDataParams()
        self.calculateSensorIntegrationParams()
        self.calculateSensorOptics()

    def calculateOrbitalParameters(self,subsetName='orbital',results = {}):
        #from IPython import embed; embed()
        if self.getU('altitude') != 'km':
            raise Exception('altitude not in km: %s' % self.getU('altitude'))
        results['orbitalRadius'] = [(self.earthR[0]+self.get('altitude'))*1e3,'m'] #earth orbital radius
        results['period'] = [(1.658669e-4)*((self.earthR[0]+self.get('altitude')))**(1.5), 'min'] #P [mins]
        results['angularVelocity'] = [np.deg2rad(6.0/results['period'][0]), 'rad/s'] #omega [deg/s] ; <= 0.071 deg/s (max angular vel for circular orbit)
        results['groundTrackVelocity'] = [2*np.pi*self.earthR[0]/(results['period'][0]*60.0), '%s/s'%self.earthR[1]] #V_g [m/s] ; <= 7905.0 m/s for circular orbit
        results['nodeShift'] = [(results['period'][0]/1436.0)*360.0, 'deg'] #dL [deg] - spacing between sucessive node crossings on the equator
        results['earthAngularRadius'] = [np.arcsin(self.earthR[0]/(self.earthR[0]+self.get('altitude'))),'rad'] #p [rad] - angular radius of spherical earth wrt spacecraft pov
        results['eclipseMax'] = [(results['earthAngularRadius'][0]/np.pi)*results['period'][0],'min']
        results['daylightMax'] = [results['period'][0]-results['eclipseMax'][0],'min']
        self.setSubset(subsetName, results)
        print 'Subset: {}: {}\n'.format(subsetName,results)
        return self

    def calculateSensorViewingParams(self, subsetName='optical',results = {}):
        #TODO: This func assumes spherical model of earth, eventually will use earthOblatenessModel()
        results['maxHorizonDistance'] = [self.earthR[0]*np.tan((np.pi/2)-self.get('earthAngularRadius')), self.earthR[1]] #D_max [m] - distance to horizon
        results['elevationAngleMin'] = [np.arccos(np.sin(np.deg2rad(self.get('nadirAngleMaxDeg'))) / np.sin(self.get('earthAngularRadius'))),'rad'] #epsilon_min [rad] - at target between spacecraft and local horizontal
        results['incidenceAngleMax'] = [(np.pi/2) - results['elevationAngleMin'][0],'rad']
        results['earthCentralAngleMax'] = [(np.pi/2) - np.deg2rad(self.get('nadirAngleMaxDeg')) - results['elevationAngleMin'][0],'rad'] #lambda [rad] - at center of earth from subsatellite point to nadirMax
        results['distanceToOffNadirMax'] = [self.earthR[0]*np.sin(results['earthCentralAngleMax'][0])/np.sin(np.deg2rad(self.get('nadirAngleMaxDeg'))),self.earthR[1]] #Dn_max [m] - i.e. Slant range; distance from satellite to max off-nadir target
        results['swathWidthAngle'] = [2*results['earthCentralAngleMax'][0],'rad'] #[rad] - determines coverage
        self.setSubset(subsetName, results)
        print 'Subset: {}: {}\n'.format(subsetName,results)
        return self
    
    def calculatePixelDataParams(self, subsetName='optical',results = {}):
        results["IFOV"] = [self.get('alongTrackGSD_ECAMax') / (self.get('distanceToOffNadirMax')*1e3),'rad'] #IFOV [rad] - Instantaneous Field of View; one pixel width
        results["crossTrackPixelRes_ECAMax"] = [self.get('alongTrackGSD_ECAMax') / np.cos(self.get('incidenceAngleMax')),self.getU('alongTrackGSD_ECAMax')] #X_max [m] - max cross-track pixel resolution @ ECAMax
        results["crossTrackGPR_Nadir"] = [results['IFOV'][0] * self.get('altitude'),self.getU('altitude')] #X [m] - cross-track Ground Pixel Resolution @ Nadir
        results["alongTrackGPR_Nadir"] = [results['IFOV'][0] * self.get('altitude'),self.getU('altitude')] #Y [m] - along-track Ground Pixel Resolution @ Nadir
        results["crossTrackPixelNum"] = [2*np.deg2rad(self.get('nadirAngleMaxDeg'))/results['IFOV'][0],'num'] #Z_c [num] - num of cross-track pixels
        results["alongTrackSwathNumRate"] = [self.get('groundTrackVelocity') / results['alongTrackGPR_Nadir'][0], 'num/s'] #Z_a [num/s] - num of swaths recorded per sec
        results["pixelRecordRate"] = [results['crossTrackPixelNum'][0] * results['alongTrackSwathNumRate'][0],'num/s'] #Z [num/s] - num of pixels recorded per sec
        results["dataRate"] = [results['pixelRecordRate'][0] * self.get('pixelBitEncodeNum'),'%sps'%self.getU('pixelBitEncodeNum')] #DR [Mbps] - data rate
        self.setSubset(subsetName, results)
        print 'Subset: {}: {}\n'.format(subsetName,results)
        return self

    def calculateSensorIntegrationParams(self, subsetName='optical',results = {}):
        #Verify detector time constant, T_det, is < pixelIntegrationPeriod in table 9-12 in SMAD
        results["pixelIntegrationPeriod"] = [self.get('alongTrackGPR_Nadir') * self.get('pixelWhiskbroomInstNum') / (self.get('groundTrackVelocity') * self.get('crossTrackPixelNum')), 'sec'] #T_i [sec]
        results["pixelReadOutFreq"] = [1.0/results["pixelIntegrationPeriod"][0],'Hz'] #F_p [Hz]
        self.setSubset(subsetName, results)
        return self

    def calculateSensorOptics(self, subsetName='optical',results = {}):
        results["focalLength"] = [self.get('altitude')*self.get('squareDetectorWidth')/self.get('crossTrackGPR_Nadir'),'m'] #f [m]
        results["defractionLimitedApertureDiameter"] = [2.44*self.get('operatingWavelength')*results['focalLength'][0]*self.get('imagingQualityFactor')/self.get('squareDetectorWidth'), 'm']
        results["opticsFNum"] = [results['focalLength'][0]/results["defractionLimitedApertureDiameter"][0],'num'] #F# [num] - typical range: 4-6 
        results["opticalFOV"] = [self.get('IFOV') * self.get('pixelWhiskbroomInstNum'),'rad'] #FOV [rad] - FOV for the array of pixels
        results["cutoffFreq"] = [results['defractionLimitedApertureDiameter'][0] / (self.get('operatingWavelength')*self.get('altitude')),'num/m'] #F_c [line pairs / m] - referred to nadir
        results["crossTrackNyquistFreq"] = [1./(2*self.get('crossTrackGPR_Nadir')),'num/m'] #F_nc [line pairs / m]
        results["alongTrackNyquistFreq"] = [1./(2*self.get('alongTrackGPR_Nadir')),'num/m'] #F_na [line pairs / m]
        results["nyquistFreqRelative"] = [[results["crossTrackNyquistFreq"][0]/results["cutoffFreq"][0],results["alongTrackNyquistFreq"][0]/results["cutoffFreq"][0]],"[%,%]"] #[F_qc,F_qa] [%,%] - % of the cutoff freq used in this case
        #results["opticsPSF"] = 
        self.setSubset(subsetName, results)
        return self

    def calculateCustomOpticalParams(self, subsetName='csOptical',results={}):
        '''
        csSensorSize
        csAspectRatio
        csLensFormat
        csFocalLength
        '''
        results['csImageRadius'] = [(self.get('csSensorSize')/2.0),'m'] #imaging radius in m(of sensor)
        results['csGroundRadius'] = [(results['csImageRadius'][0]*(self.get('altitude')))/self.get('csFocalLength'),'m'] #max ground radius in m (cicular)
        results['csDectectorArea'] = [np.pi*(results['csImageRadius'][0]**2),'m^2'] #area of sensor
        results['csFovDiameter'] = [2*np.arctan(results['csImageRadius'][0]/self.get('csFocalLength')),'rad'] #angular diameter of FOV in radians
        #areaG = pi*(groundRadius**2) #ground object FOV in m squared (circular)
        results['csDiagonalLength'] = [(self.get('csSensorSize')/self.get('csLensFormat')*(2*(results['csGroundRadius'][0]))),'m'] #ground radius adjusted for sensor size and lens radius
        results['csArea2SensorHypotenuseRatio'] = [results['csDiagonalLength'][0]/(self.get('csSensorSize')),'num'] #finding ratio of area to sensor hypotenous
        results['csSensorDiagonal'] = [np.linalg.norm(self.get('csAspectRatio')),'num'] #finding sensor diagonal
        results['csSensorRatio'] = [self.get('csSensorSize')/results['csSensorDiagonal'][0],'num'] #the ratio of the given sensor size to the ratio based on its aspects
        results['csImageArea'] = [(results['csArea2SensorHypotenuseRatio'][0]*(np.array(self.get('csAspectRatio'))*results['csSensorRatio'][0])).tolist(),'[m,m]']
        self.setSubset(subsetName, results)
        return self

    def calculateGroundStationParams(self, subsetName='groundStation', results={}):
        #SMAD - pg. 121
        #TODO verify results
        latPole = np.deg2rad(90.0 - self.get('inclination'))
        longPole = np.deg2rad(self.get('lNode') - 90.0)
        deltaLong = np.abs(longPole - np.deg2rad(self.get('gsLong')))
        results['gsEarthCentralAngleMin'] = [np.arcsin(np.sin(latPole)*np.sin(np.deg2rad(self.get('gsLat'))) + np.cos(latPole)*np.cos(np.deg2rad(self.get('gsLat')))*np.cos(deltaLong)),'rad']
        results['gsNadirAngleMin'] = [np.arctan((np.sin(self.get('earthAngularRadius'))*np.sin(results['gsEarthCentralAngleMin'][0]))/(1 - np.sin(self.get('earthAngularRadius'))*np.cos(results['gsEarthCentralAngleMin'][0]))),'rad']
        results['gsElevationAngleMax'] = [(np.pi/2) - results['gsEarthCentralAngleMin'][0] - results['gsNadirAngleMin'][0],'rad']
        results['gsDistanceMin'] = [self.earthR[0]*(np.sin(results['gsEarthCentralAngleMin'][0])/np.sin(results['gsNadirAngleMin'][0])),'m']
        results['gsAngularRateMax'] = [np.deg2rad((2*np.pi*(self.earthR[0]+self.get('altitude'))) / (self.get('period')*results['gsDistanceMin'][0]))/60.0,'rad/s']
        results['gsAzimuthRange'] = [2*np.arccos(np.tan(results['gsEarthCentralAngleMin'][0]) / np.tan(self.get('earthCentralAngleMax'))),'rad']
        results['gsViewTime'] = [(self.get('period') / 180.0)*np.arccos(np.cos(self.get('earthCentralAngleMax'))/results['gsEarthCentralAngleMin'][0]),'min']
        results['gsViewTimeMax'] = [self.get('period')*(self.get('earthCentralAngleMax')/np.pi),'min']
        self.setSubset(subsetName, results)
        return self

    def calculateCommsParams(self, subsetName='comms', results={}):
        F = (1/self.get('earthCentralAngleMax'))*np.arccos(np.cos(self.get('earthCentralAngleMax'))/np.cos(self.get('gsEarthCentralAngleMin')))
        F *= 0.8 #80% for Circular LEO, and 86% or more of all passes will have F > 0.5
        results["dataQuantity"] = self.get('dataRate')*(F*self.get('gsViewTimeMax') - self.get('gsInitCommTime'))
        #energyPerBit2noiseDensityRatio of 5-10 adaquate for recieving binary data with low probability of error with some forward error correction
        results["energyPerBit2noiseDensityRatio"] = self.get('transmitterPower')*self.get('transmitter2antennaGainLineLoss')*self.get('transmitAntennaGain')*\
            self.get('spaceLoss')*self.get('transmissionPathLoss')*self.get('receiveAntennaGain')/(scipy.constants.Boltzmann*self.get('systemNoiseTemperature')*self.get('dataRate'))

    def calculateEpsParams(self, subsetName='eps', results={}):
        results['solarArrayPowerOutputDaylight'] = [(self.get('avgTotalPower')*((self.get('eclipseMax')/self.get('effSA2Batt2Load'))+(self.get('daylightMax')/self.get('effSA2Load'))))/self.get('daylightMax'),'W'] # P_sa [W] - amount of power that must be produced by arrays
        results['powerOutputSunNormalEst'] = [self.get('solarCellEff')*self.get('solarIlluminationIntensity'),'W/m^2'] # P_o - power output with sun normal to the surface of the sun
        results['powerProductionBOL'] = [results['powerOutputSunNormalEst'][0]*self.get('solarArrayInherentDegradation')*np.cos(np.deg2rad(self.sunIncidentAngleDeg[0])),'W/m^2'] # P_BOL - power production at beginning of life
        results['lifeDegradation'] = [(1.0 - self.get('solarCellPerformanceDegradation'))**self.get('lifetimeNominal'),'num'] # L_d - life degradation
        results['powerProductionEOL'] = [results['powerProductionBOL'][0]*results['lifeDegradation'][0],'W/m^2'] # P_EOL - power production at the end of life
        results['solarArrayAreaEOL'] = [results['solarArrayPowerOutputDaylight'][0]/results['powerProductionEOL'][0],'m^2'] # A_sa - area of solar arrays required to produce the necessary power (P_sa)
        results['solarArrayMassEst'] = [results['solarArrayPowerOutputDaylight'][0]*(1.0/(self.get('solarArraySpecificPerformance'))),'kg']
        self.setSubset(subsetName, results)
        print 'Subset: {}: {}\n'.format(subsetName,results)
        return self

    def calculateADCSParams(self, subsetName='adcs', results={}):
        results['disturbanceTorque_gravityGradient'] = [((3*self.earthGravConst)/(2*(self.get('orbitalRadius')**3))) * np.abs(np.max(self.get('PMOI')) - np.min(self.get('PMOI'))) * (np.sin(np.deg2rad(2*self.get('nadirAngleMaxDeg')))),'Nm']
        solarRad_F = (self.get('solarIlluminationIntensity')/self.lightSpeed) * self.get('solarSurfaceArea') * (1 + self.get('reflectanceFactor')) * (np.cos(np.deg2rad(self.get('adcsIncidenceAngle'))))
        results['disturbanceTorque_solarRadiation'] = [solarRad_F*(self.get('centreSolarPress')-self.get('COM')),'Nm']
        aero_F = 0.5*self.get('atmosphericDensity')*self.get('aeroDragCoeff')*self.get('aeroSurfaceArea')*(self.get('groundTrackVelocity')*1e3)**2 #km to m
        results['disturbanceTorque_aerodynamic'] = [aero_F*(self.get('centreAeroDrag')-self.get('COM')),'Nm']
        earthMagneticField = (2*self.earthMagneticMoment/(self.get('orbitalRadius'))**3)
        results['disturbanceTorque_magneticField'] = [self.get('residualDipole')*earthMagneticField,'Nm']
        maxDisturbanceTorque = np.max([results['disturbanceTorque_gravityGradient'][0],results['disturbanceTorque_solarRadiation'][0],results['disturbanceTorque_aerodynamic'][0],results['disturbanceTorque_magneticField'][0]])
        results['reactionWheelDisturbanceRejection'] = [self.get('disturbanceMarginFactor')*maxDisturbanceTorque,'Nm']
        results['slewTorque'] = [4*np.deg2rad(self.get('nadirAngleMaxDeg'))*np.max(self.get('PMOI'))/(self.get('slewTime_nadirAngleMaxDeg'))**2,'Nm']
        results['momentumStorage'] = [maxDisturbanceTorque*self.get('period')*60/4.*0.707,'Nms']

        results['reactionWheelTorque'] = [results['momentumStorage'][0]*self.get('componentSelectionMargin'),'Nms']
        results['magnetTorquerTorque'] = [(maxDisturbanceTorque/earthMagneticField)*self.get('componentSelectionMargin'),'A/m^2']
        
        embed()
        self.setSubset(subsetName, results)
        print 'Subset: {}: {}\n'.format(subsetName,results)
        return self

    def calculateEpsBat(self, subsetName='eps', results={}):
        results['BattCapacity'] = [self.get('avgTotalPower')*self.get('eclipseMax')/(60*self.get('depthDischarge')*self.get('BattQuant')*self.get('effBatt2Load')),'Whr']
        self.setSubset(subsetName, results)
        print 'Subset: {}: {}\n'.format(subsetName,results)
        return self

    def earthOblatenessModeling(satVectorEFF):
        f = 1/298.257 # f - earth flattening factor ; f = (R_e - R_p) / R_e : R_e (equatorial radius) ~ 6378.140, R_p (polar radius)
        a = np.sqrt(np.dot(satVectorEFF**2,[1,1,(1 - earthFlattening)**2])) #a - equatorial radius, c - polar radius
        #lat (lambda), long (phi) = geocentric lat, long of the observers (sats) position
        #azi = azimuth angle of the horizon vector
        R = a*(1-f) / np.sqrt(1 - (2 - f)*f*np.cos(lat)**2) #Earth Radius to local surface
        d = np.linalg.norm(satVectorEFF) #distance from earth center to satellite in EFF
        term1 = (((d**2 - R**2)/a**2)*(1 + ((2-f)*f*R**2*np.cos(lat)**2*np.sin(azi)**2 / (1-f)**2*a**2)))**0.5
        earthAngularRadius = 1.0 / np.arctan(term1 + ((2-f)*f*R**2*np.sin(azi) / 2*(1-f)**2)*a**2)
        '''
        H = np.sqrt(np.sum(satVectorEFF**2)-a**2) #init estimate (assumes spherical earth)
        #iteratively compute max distance to horizon
        for n in range(10):
            H = np.sqrt(np.sum([]))
        '''




