#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function /WAVE avg_raveel(int wavenum, string dataset) // /WAVE lets your return a wave

// averaging over total sweeps without centering or removing any data
//
//inputs = wavenum: datnum
//         dataset: "cscurrent_2d","dotcurrent_2d"
	
	variable i
	string w2d
	int nc 
	int nr 
	
	w2d="dat"+num2str(wavenum)+dataset
	wave wavenm = $w2d
	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)
	
	make /n=(nc, 1 ) /o avgcurrent
	avgcurrent = wavenm[0][p] //first row
	
	for (i=1; i < nr ; i+=1)
		avgcurrent += wavenm[i][p] //sums all the rows
	endfor
	
	avgcurrent = avgcurrent / nr
	display avgcurrent
	
	Label bottom "voltage"
    Label left dataset
	
	return avgcurrent
	
end


	
function plot2d_raveel(num,dataset)
//plots the repeats against the sweeps for the specified dataset
//
//inputs =  num: datnum
//          dataset: "cscurrent_2d","dotcurrent_2d"

	variable num 
	string dataset
	string wvname
	
	wvname="dat"+num2str(num)+dataset
	
	wave wavenm = $wvname
			
	if (dimsize(wavenm,0)<151)
		matrixtranspose wavenm
	endif
	
	display; //start with empty graph
	appendimage wavenm //append image of data
	ModifyImage $wvname ctab= {*,*,Turbo,0} //setting color (idk why it prefers the pointer)
	ColorScale /A=RC /E width=20 //puts it on the right centre, /E places it outside

    Label bottom "voltage"
    Label left dataset

    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300

end




// plotted thetas from fitting each sweep = good and bad for each line
	
	//need to figure out fitting in Igor
	//gradient descent helped figure out a specific parameter, need to look at dat analysis to recall
	// amplitude was just the max and min value subtracted, theta was set to 5
	
	

// plotted the bad thetas only (above some standard deviation)
	//nice to see the errors
	

// plotted the avg of the good thetas with the avg theta value
	//I belive this averaging takes the good thetas, the data is centred then averaged.
	// I averaged over all data, 


//Savitsky golay smoothing exists in Igor
	//https://www.wavemetrics.com/products/igorpro/dataanalysis/signalprocessing/smoothing