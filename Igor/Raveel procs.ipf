#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function avg_raveel(wavenum)

// averaging over total sweeps without centering.

	variable wavenum
	variable i
	string w2d
	int nc 
	int nr 
	
	w2d="dat"+num2str(wavenum)+"cscurrent_2d"
	wave wavenm = $w2d
	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)
	
	make /n=(nc, 1 ) /o avgcscurrent
	avgcscurrent = wavenm[0][p] //first row
	
	for (i=1; i < nr ; i+=1)
		avgcscurrent += wavenm[i][p] //sums all the rows
	endfor
	
	avgcscurrent = avgcscurrent / nr
	
	display avgcscurrent
	
	//return avgcscurrent
	
end



