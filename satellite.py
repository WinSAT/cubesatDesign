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

class Satellite(object):

    def __init__(self, initData, subsetList=['optical','orbital']):
        setattr(self,'initParams',initData)
        setattr(self,'subsetList', subsetList)
        for subset in initData:
            setattr(self, subset, initData[subset])
        self.earthR = 6378.14e3 #[m] - Earth Radius using Spherical Model
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
        for subset in self.subsetList:
            try:
                value = getattr(self,subset)[param][0] 
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
        results['period'] = [(1.658669e-4)*(self.earthR+self.get('altitude'))**(1.5), 'mins'] #P [mins]
        results['angularVelocity'] = [np.deg2rad(6.0/results['period'][0]), 'rad/s'] #omega [deg/s] ; <= 0.071 deg/s (max angular vel for circular orbit)
        results['groundTrackVelocity'] = [2*np.pi*self.earthR/results['period'][0], 'm/s'] #V_g [m/s] ; <= 7905.0 m/s for circular orbit
        results['nodeShift'] = [(results['period'][0]/1436.0)*360.0, 'deg'] #dL [deg] - spacing between sucessive node crossings on the equator
        self.setSubset(subsetName, results)
        return self

    def calculateSensorViewingParams(self, subsetName='optical',results = {}):
        #TODO: This func assumes spherical model of earth, eventually will use earthOblatenessModel()
        results['earthAngularRadius'] = [np.arcsin(self.earthR/(self.earthR+self.get('altitude'))),'rad'] #p [rad] - angular radius of spherical earth wrt spacecraft pov
        results['maxHorizonDistance'] = [self.earthR*np.tan(np.deg2rad(90.0)-results['earthAngularRadius'][0]),'m'] #D_max [m] - distance to horizon
        results['elevationAngleMin'] = [np.arccos(np.sin(np.deg2rad(self.get('nadirAngleMaxDeg'))) / np.sin(results['earthAngularRadius'][0])),'rad'] #epsilon_min [rad] - at target between spacecraft and local horizontal
        results['incidenceAngleMax'] = [(np.pi/2) - results['elevationAngleMin'][0],'rad']
        results['earthCentralAngle'] = [(np.pi/2) - np.deg2rad(self.get('nadirAngleMaxDeg')) - results['elevationAngleMin'][0],'rad'] #lambda [rad] - at center of earth from subsatellite point to nadirMax
        results['distanceToMaxOffNadir'] = [self.earthR*np.sin(results['earthCentralAngle'][0])/np.sin(np.deg2rad(self.get('nadirAngleMaxDeg'))),'m'] #Dn_max [m] - i.e. Slant range; distance from satellite to max off-nadir target
        results['swathWidthAngle'] = [2*results['earthCentralAngle'][0],'rad'] #[rad] - determines coverage
        self.setSubset(subsetName, results)
        return self
    
    def calculatePixelDataParams(self, subsetName='optical',results = {}):
        results["IFOV"] = [self.get('alongTrackGSD_ECAMax') / self.get('distanceToMaxOffNadir'),'rad'] #IFOV [rad] - Instantaneous Field of View; one pixel width
        results["crossTrackPixelRes_ECAMax"] = [self.get('alongTrackGSD_ECAMax') / np.cos(self.get('incidenceAngleMax')),'m'] #X_max [m] - max cross-track pixel resolution @ ECAMax
        results["crossTrackGPR_Nadir"] = [results['IFOV'][0] * self.get('altitude'),'m'] #X [m] - cross-track Ground Pixel Resolution @ Nadir
        results["alongTrackGPR_Nadir"] = [results['IFOV'][0] * self.get('altitude'),'m'] #Y [m] - along-track Ground Pixel Resolution @ Nadir
        results["crossTrackPixelNum"] = [2*np.deg2rad(self.get('nadirAngleMaxDeg'))/results['IFOV'][0],'num'] #Z_c [num] - num of cross-track pixels
        results["alongTrackSwathNumRate"] = [self.get('groundTrackVelocity') / results['alongTrackGPR_Nadir'][0], 'num/s'] #Z_a [num/s] - num of swaths recorded per sec
        results["pixelRecordRate"] = [results['crossTrackPixelNum'][0] * results['alongTrackSwathNumRate'][0],'num/s'] #Z [num/s] - num of pixels recorded per sec
        results["dataRate"] = [results['pixelRecordRate'][0] * self.get('pixelBitEncodeNum'),'Mbps'] #DR [Mbps] - data rate
        self.setSubset(subsetName, results)
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




