function compute_fbb, nuvmag, anglearcmin
		ft0 = (nuvmag-19.68)*(-0.9d-2) > (-8.3d-3) < 0
		ft1 =  (nuvmag-20.35)*(-1) > (nuvmag-20.35)*(0.2) - 1.4
		return, ft0*(anglearcmin^ft1)		
end