#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later


function /wave quick_avg2(int wavenum, string dataset, int view) // /WAVE lets your return a wave

// averaging over total sweeps without centering or removing any data
	
	variable i
	string w2d
	string w2x
	string avg_name
	int nc 
	int nr 
	
	w2d="dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	avg_name = "dat"+num2str(wavenum)+"quickavg" //new array
	
	wave wavenm = $w2d
	
	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)
	
	
	
	make /n=(nc, 1) /o $avg_name
	wave avgwave = $avg_name
	
	avgwave = wavenm[0][p]                       //first row
	
	for (i=1; i < nr ; i+=1)
		avgwave += wavenm[i][p]                 //sums all the rows
	endfor
	
	avgwave = avgwave / nr                     //divide by total rows
	
	
	
	matrixtranspose wavenm						//needs to be transposed so the correct scaling is used
	copyscales wavenm avgwave                   //setting x-scaling
	
	if (view == 1)
	
		display avgwave
		Label bottom "voltage"
    	Label left "current"
    	ModifyGraph fSize=24
    	ModifyGraph gFont="Gill Sans Light"
    	ModifyGraph width={Aspect,1.62},height=300
    	ModifyGraph mode=2,lsize=2,rgb=(21845,21845,21845)
		Legend/C/N=text0/J/A=RT "\\Z14\\Z16\\s(avgcurrent) quick average of dat" + num2str(wavenum)
	
	endif
	
	return avgwave
	
end
