#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function /wave avg_raveel(int wavenum, string dataset) // /WAVE lets your return a wave

// averaging over total sweeps without centering or removing any data
//
//inputs = wavenum: datnum
//         dataset: "cscurrent_2d","dotcurrent_2d"
	
	variable i
	string w2d
	string w2x
	int nc 
	int nr 
	
	w2d="dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	
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
	
	avgcurrent = avgcurrent / nr // divide by total rows
	
	display avgcurrent[][0] vs $w2x
	Label bottom "voltage"
    Label left dataset
    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
	
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


//Function to get initial paramaters

function /wave get_initial_params(sweep, x_array)

// for a given sweep returns a guess of initial parameters for the fit function: Charge transiton
	wave sweep //have to declare function parameters first
	wave x_array
//	
//	
//	
	variable amp = wavemax(sweep) - wavemin(sweep) //might be worthwile looking at a maximum/ minimum near the middle
	variable const = mean(sweep)
	variable theta = 5
	
	duplicate /o sweep sweepsmooth
	Smooth/S=4 101, sweepsmooth ;DelayUpdate
	
	differentiate sweepsmooth
    extract/INDX sweepsmooth, extractedwave, sweepsmooth == wavemin(sweepsmooth)
    variable mid = x_array[extractedwave[0]]

	
	//smoothing and differentiating //could be useful for amp guess and linear term guess
	//Duplicate/O avgcurrent,avgcurrent_smth;DelayUpdate
	//Smooth/S=4 121, avgcurrent_smth;DelayUpdate
	//differentiate avgcurrent_smth

	
	variable lin = 0.001  // differentiated value of flat area?
	
	Make /D/N=5/O W_coef
	W_coef[0] = {amp,const,theta,mid,lin}
	
	killwaves extractedwave, sweepsmooth
	return W_coef

end

	
	
	
//function to perform a fit

function /wave fit_transition(current_array,x_array)
	
	wave current_array
	wave x_array

	wave W_coef
	
	//string wnamex
	
	//wnamex = "dat"+num2str(datnum)+"x_array"
	
	
	//Make/D/N=5/O W_coef
	//W_coef = get_initial_params(sweep)
	
	get_initial_params(current_array, x_array)
	
	//smoothing and differentiating //could be useful for amp guess and linear term guess
	Duplicate/O current_array, smooth_current;DelayUpdate
	Smooth/S=4 201, smooth_current;DelayUpdate
	
	
	FuncFit Chargetransition W_coef smooth_current[][0] /X= x_array /D 
	
	//return W_coef
end
	


//function to perform all the fits

function /wave get_fit_params(int wavenum, string dataset)
//
	variable i
	string w2d
	string w2x
	int nc 
	int nr
	wave fit_params
	wave temp_wave
	wave W_coef
	wave W_sigma
	 
//	
	w2d="dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
//	
	wave wavenm = $w2d
//	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
//	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)
//	
//	
	make /N= (2 * nr, 5) /O fit_params
	make /N= (nc) /O temp_wave
	  
	for (i=0; i < nr ; i+=1)
		//create wave for for each i'th row
		//run fit_transition(current_array, $w2x)
		//store in corresponding 2d array
		
	  temp_wave = wavenm[i][p]
      fit_transition(temp_wave, $w2x)
      fit_params[1 * i] = W_coef[q]
      fit_params[(nr/2) + i] = W_sigma[q] 
//		
	endfor
	
	return fit_params
	
end
	
	

//function plot_thetas(int wavenum, string dataset)
//	
//	wave fit_params
//	int nr
//	get_fit_params(wavenum, dataset)
//	nr = dimsize(fit_params,0)
//	
//	display fit_params[0,nr/2 - 1][2]
//	
//	ModifyGraph fSize=24
//    ModifyGraph gFont="Gill Sans Light"
//    ModifyGraph width={Aspect,1.62},height=300
//    ModifyGraph mode=3,rgb=(0,0,0,32768)
//    
//    Label bottom "repeat"
//    Label left "theta values"
//    
//    
//    
//	
//end


function plot_thetas(int wavenum, string dataset)
	
	wave fit_params
	variable thetamean
	variable thetastd
	variable i
	int nr
	
	get_fit_params(wavenum, dataset)
	nr = dimsize(fit_params,0)
	nr = nr/2-1               //only the first half are associated with theta values,
	                          //the second have are the uncertainties 
	 
	duplicate /O/R =[0,nr][2] fit_params thetas
	
	thetamean = mean(thetas)
	thetastd = sqrt(variance(thetas))
	
	make /o/n =(nr) meanwave
	make /o/n =(nr) stdwave
	make /o/n =(nr) stdwave2
	make /o/n = 0 goodthetas
	make /o/n = 0 goodthetasx
	make /o/n = 0 badthetas
	make /o/n = 0 badthetasx
	
	
	meanwave = thetamean
	stdwave = thetamean - 2 * thetastd
	stdwave2 = thetamean + 2 * thetastd
	
	
	//display thetas, meanwave, stdwave, stdwave2
	
	
	for (i=0; i < nr ; i+=1)
		
		if (abs(thetas[i] - thetamean) < (2 * thetastd))
			
			insertPoints /v = (thetas[i]) nr, 1, goodthetas
			insertpoints /v = (i) nr, 1, goodthetasx
			
		else
		
			insertPoints /v = (thetas[i]) nr, 1, badthetas
			insertpoints /v = (i) nr, 1, badthetasx
			
		endif
		
	endfor
			
		
	
	
	display meanwave, stdwave, stdwave2
	appendtograph goodthetas vs goodthetasx
	appendtograph badthetas vs badthetasx	
	
	
	ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    ModifyGraph lstyle(meanwave)=3,rgb(meanwave)=(17476,17476,17476)
    ModifyGraph lstyle(stdwave)=3,rgb(stdwave)=(30583,30583,30583)
    ModifyGraph lstyle(stdwave2)=3,rgb(stdwave2)=(30583,30583,30583)
    ModifyGraph mode(goodthetas)=3,lsize(goodthetas)=2
    ModifyGraph mode(badthetas)=3,rgb(badthetas)=(29524,1,58982)

    
    
    
    Label bottom "repeat"
    Label left "theta values"
    
    
    
	
end



function plot_badthetas(int wavenum, string dataset)

	int i 
	int nr
	wave badthetasx
	string w2d
	string w2x
	
	w2d= "dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	
	wave wavenm = $w2d
	duplicate /o wavenm, wavenmcopy
	
	plot_thetas(wavenum, dataset)
	nr = dimsize(badthetasx,0)
	
	
	if (dimsize(wavenmcopy,1)<151)
		matrixtranspose wavenmcopy
	endif
	
	display
	
	for(i=0; i < nr; i +=1)
		appendtograph wavenmcopy[badthetasx[i]][] vs $w2x
	endfor

	QuickColorSpectrum()
	
	ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
	Label bottom "voltage"
    Label left dataset
	
end 
	
	
function full_procedure(int wavenum, string dataset)

	avg_raveel(wavenum, dataset)
	plot_badthetas(wavenum, dataset)
	plot2d_raveel(wavenum,dataset)

end



//from: https://www.wavemetrics.com/forum/igor-pro-wish-list/automatically-color-traces-multi-trace-graph

Function QuickColorSpectrum()                            // colors traces with 12 different colors
    String Traces    = TraceNameList("",";",1)               // get all the traces from the graph
    Variable Items   = ItemsInList(Traces)                   // count the traces
    Make/FREE/N=(11,3) colors = {{65280,0,0}, {65280,43520,0}, {0,65280,0}, {0,52224,0}, {0,65280,65280}, {0,43520,65280}, {0,15872,65280}, {65280,16384,55552}, {36864,14592,58880}, {0,0,0},{26112,26112,26112}}
    Variable i
    for (i = 0; i <DimSize(colors,1); i += 1)
        ModifyGraph rgb($StringFromList(i,Traces))=(colors[0][i],colors[1][i],colors[2][i])      // set new color offset
    endfor
End





// plotted thetas from fitting each sweep = good and bad for each line

//avg_raveel(3920,"cscurrent_2d")
//ModifyGraph rgb=(0,0,0)
//Make/D/N=5/O W_coef
//W_coef[0] = {0.02,1.17,-0.1,0,1}
//FuncFit Chargetransition W_coef avgcurrent[][0] /X=dat3920x_array /D 


	
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
	
	
	
	//duplicate /O/R = fit_params[0,nr/2 - 1][2] fit_params, wave2