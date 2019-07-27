import numpy as np
import satellite
from IPython import embed

FireSatParams = {
	'orbital': {
		'lifetimeNominal': [5, 'yr'],
		'altitude': [700.0,'km'], #[m]
		'inclination': [50.0,'deg'], #WinSAT
		'lNode': [0.0, 'deg'] #TODO ?? - Instantaneous Ascending Node, RAAN? 
	},
	'optical': {
		'nadirAngleMaxDeg': [57.9,'deg'], #eta [deg] - max target range from sat to off-nadir target
		'alongTrackGSD_ECAMax': [68.0,'m'], #Y_max [m] - max along-track Ground Sampling Distance; design param
		'pixelBitEncodeNum': [8,'b'], #B [bits] - num of bits used to encode each pixel
		'pixelWhiskbroomInstNum': [256,'num'], #N_m [num] - must be large enough to allow sufficient integration time; design param
		'squareDetectorWidth': [30e-6,'m'], #d [m]; design param
		'imagingQualityFactor': [1.1,'num'], #Q [num] - 0.5<Q<2 (1.1 for good img quality); design param
		'operatingWavelength': [4.2e-6,'m'], #lambda [m] - based on subject trades; design param
	},
	'csOptical': {
		'csSensorSize': [4.0e-3,'m'], #sensor size in mm
		'csAspectRatio': [[4.0,3.0],'[num,num]'], #aspect ratio of sensor mm
		'csLensFormat': [6.0e-3,'m'], #lens format in mm
		'csFocalLength': [16.0e-3,'m']
	},
	'groundStation': {
		'gsLat': [42.304524,'deg'], #WinSAT - UWindsor CEI 42.304524, -83.062185
		'gsLong': [-83.062185,'deg'],
		'gsInitCommTime': [2.0, 'min'],
		'gsMissedPassesMargin': [2.5, 'num'] #~2-3
	},
	'eps': {
		'avgTotalPower': [110.0, 'W'],
		'effSA2Batt2Load': [0.6, 'num'], #X_e - efficiency of Solar Arrays to Batt to Load - 0.6 for PPT, 0.65 for Direct Energy
		'effSA2Load': [0.8, 'num'], #X_e - efficiency of Solar Arrays to Load directly - 0.8 for PPT, 0.85 for Direct Energy
		'solarCellEff': [0.148, 'num'], #0.148 for Si, 0.185 for GaAs, 0.22 for Multijunction 
		'solarArrayInherentDegradation': [0.77, 'num'], #nominal 0.77, range 0.49-0.88
		'solarCellPerformanceDegradation': [0.0375, 'num/yr'], #3.75% for Si, 2.75% for GaAs, 0.5% for Multijunction
		'solarArraySpecificPerformance': [25.0,'W/kg'],
		'BattQuant': [3,'num'], # N - number of batteries
		'effBatt2Load':[0.9,'num'], # n - Transmission efficiency between the battery and the load
		'depthDischarge':[.2,'num'], # DOD - Depth of Dischanrge
	}
}

WinSATParams = {
	'orbital': { #also general parameters
		'altitude': [500.0,'km'], #[m]
		'inclination': [51.6,'deg'], #WinSAT
		'lNode': [0.0, 'deg'], #TODO ?? - Instantaneous Ascending Node, RAAN? 
		'solarIlluminationIntensity': [1367.0, 'W/m^2'],
		'PMOI': [[5451073.76e-9, 5430996.67e-9, 1358707.37e-9],'[kgm^2,kgm^2,kgm^2]'], #kgm^2 - x,y,z - nadir, orthonormal, velcity vector
	},
	'optical': {
		'nadirAngleMaxDeg': [25.0,'deg'], #eta [deg] - max target range from sat to off-nadir target
		'alongTrackGSD_ECAMax': [50.0,'m'], #Y_max [m] - max along-track Ground Sampling Distance; design param
		'pixelBitEncodeNum': [8,'b'], #B [bits] - num of bits used to encode each pixel
		'pixelWhiskbroomInstNum': [256,'num'], #N_m [num] - must be large enough to allow sufficient integration time; design param
		'squareDetectorWidth': [30e-6,'m'], #d [m]; design param
		'imagingQualityFactor': [1.1,'num'], #Q [num] - 0.5<Q<2 (1.1 for good img quality); design param
		'operatingWavelength': [4.2e-6,'m'], #lambda [m] - based on subject trades; design param
	},
	'csOptical': {
		'csSensorSize': [4.0e-3,'m'], #sensor size in mm
		'csAspectRatio': [[4.0,3.0],'[num,num]'], #aspect ratio of sensor mm
		'csLensFormat': [6.0e-3,'m'], #lens format in mm
		'csFocalLength': [16.0e-3,'m']
	},
	'groundStation': {
		'gsLat': [42.304524,'deg'], #WinSAT - UWindsor CEI 42.304524, -83.062185
		'gsLong': [-83.062185,'deg'],
		'gsInitCommTime': [2.0, 'min'],
		'gsMissedPassesMargin': [2.5, 'num'] #~2-3
	},
	'adcs': {
		'reflectanceFactor': [0.6, 'num'],
		'adcsIncidenceAngle': [0, 'deg'],
		'centreSolarPress': [0.05, 'm'],
		'solarSurfaceArea': [0.3405*0.1,'m^2'],
		'COM': [0, 'm'], #centre of mass
		
		'centreAeroDrag': [0.3405/2, 'm'],
		'atmosphericDensity' : [1e-13,'kg/m^3'],
		'aeroDragCoeff': [2.2,'num'],
		'aeroSurfaceArea': [0.01, 'm^2'],

		'residualDipole' : [1.0, 'A/m^2'],
		'slewTime_nadirAngleMaxDeg': [60.0, 's'],

		'disturbanceMarginFactor': [5, 'num'],
		'componentSelectionMargin': [9, 'num']

	}

}

sat = satellite.Satellite(WinSATParams)
#sat = satellite.Satellite(WinSATParams)

sat.calculateOrbitalParameters()
sat.calculateADCSParams()
sat.calculateSensorViewingParams()
sat.calculatePixelDataParams()
sat.calculateEpsParams()
sat.calculateEpsBat()