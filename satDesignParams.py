import numpy as np
import satellite

FireSatParams = {
	'orbital': {
		'altitude': [700.0e3,'m'], #[m]
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
	'groundStation': {
		'gsLat': [42.304524,'deg'], #WinSAT - UWindsor CEI 42.304524, -83.062185
		'gsLong': [-83.062185,'deg'],
	}
}

WinSATParams = {
	'altitude': 500.0e3, #[m]
	'nadirAngleMaxDeg': 25.0, #eta [deg] - max target range from sat to off-nadir target
	'alongTrackGSD_ECAMax': 50.0, #Y_max [m] - max along-track Ground Sampling Distance @ ECAMax; design param
	'pixelBitEncodeNum': 8, #B [num] - num of bits used to encode each pixel
	'pixelWhiskbroomInstNum': 256, #N_m [num] - must be large enough to allow sufficient integration time; design param
	'squareDetectorWidth': 30e-6, #d [m]; design param
	'imagingQualityFactor': 1.1, #Q [num] - 0.5<Q<2 (1.1 for good img quality); design param
	'operatingWavelength': 4.2e-6, #lambda [m] - based on subject trades; design param
}

FireSat = satellite.Satellite(FireSatParams)
FireSat.calculateOrbitalParameters()
FireSat.calculateOpticalParameters()
FireSat.calculateGroundStationParams()

from IPython import embed; embed()