#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later


function /wave quick_avg2(int wavenum, string dataset, int view) // /WAVE lets your return a wave

//  averaging over total sweeps without centering or removing any data
// wave returned is named "dat" + "wavenum" + "quickavg" eg. dat3320quickavg


//    Inputs: wavenum is the dat number
//            dataset implies either "cscurrent_2d" or "dotcurrent_2d"
//            view if set to 1, will display the quick_avg2 plot

	
	variable i
	string w2d
	string w2x
	string avg_name
	int nc 
	int nr 
	
	w2d="dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	avg_name = "dat"+num2str(wavenum)+"quickavg" //new array
	
	duplicate /o $w2d wavenm
	
	
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


function plot2d_heatmap2(num,dataset)

//plots the repeats against the sweeps for dataset cscurrent_2d

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
    Label left "repeats"

    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 raw 2D plot of dat" + num2str(num)

end


function /wave get_initial_params2(sweep)

// for a given sweep returns a guess of initial parameters for the fit function: Charge transiton

	wave sweep 
	duplicate /o sweep x_array
	x_array = x
		
	variable amp = wavemax(sweep) - wavemin(sweep) //might be worthwile looking for a maximum/minimum with differentiation
	//variable amp = 0.001
	variable const = mean(sweep)
	variable theta = 5
	
	duplicate /o sweep sweepsmooth
	Smooth/S=4 201, sweepsmooth ;DelayUpdate
	
	differentiate sweepsmooth
    extract/INDX sweepsmooth, extractedwave, sweepsmooth == wavemin(sweepsmooth)
    variable mid = x_array[extractedwave[0]]

	//extract/INDX sweepsmooth, extractedwave, sweepsmooth == 0 //new
	//variable amp = sweep[extractedwave[0]] - sweep[extractedwave[1]] // new
	
	
	variable lin = 0.001  // differentiated value of flat area?
	
	Make /D/N=5/O W_coef
	W_coef[0] = {amp,const,theta,mid,lin}
	
	killwaves extractedwave, sweepsmooth
	return W_coef

end



function /wave fit_transition2(current_array, condition)
// fits the current_array, If condition is 0 it will get initial params, If 1:
// define a variable named W_coef_guess = {} with the correct number of arguments

	
	wave current_array
	int condition

	wave W_coef_guess
	wave W_coef
	wave smooth_current
	
	//duplicate /o current_array x_array
	//x_array = x
	
	
	if(condition == 0)
		get_initial_params2(current_array)
	else
		duplicate W_coef_guess, W_coef
	endif
	
	//smoothing and differentiating //could be useful for amp guess and linear term guess
	
	Duplicate /o current_array smooth_current;DelayUpdate
	Smooth/S=4 51, smooth_current;DelayUpdate                           //double check all the smoothing
	FuncFit Chargetransition W_coef smooth_current[][0] /D    //removed the x_array
	
end




function /wave get_fit_params2(int wavenum, string dataset, int condition)

	variable i
	string w2d
	string fit_params_name
	int nc 
	int nr
	wave fit_params
	wave temp_wave
	wave W_coef
	wave W_sigma
	 

	w2d="dat"+num2str(wavenum)+dataset //current 2d array
	fit_params_name = "dat"+num2str(wavenum)+"fit_params" //new array
	wave wavenm = $w2d
	
	
	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)

	make /N= (2* nr , 5) /o $fit_params_name
	wave fit_params = $fit_params_name
	
	
	make /N= (nc) /o temp_wave
	  
	for (i=0; i < nr ; i+=1)
	
		//create wave for for each i'th row
		//run fit_transition(current_array, $w2x)
		//store in corresponding 2d array
		
	  temp_wave = wavenm[i][p]
	  
	  	matrixtranspose wavenm						//needs to be transposed so the correct scaling is used
	  	copyscales wavenm temp_wave
	  	matrixtranspose wavenm
	  
      fit_transition2(temp_wave, condition)
      fit_params[1 * i] = W_coef[q]
      fit_params[(nr/2) + i] = W_sigma[q] 
		
	endfor
	
	return fit_params
	
end







function chargetransition_procedure2(int wavenum, int condition)

  // if condition is 0, It gets a guess of initial params from each sweep
  // If condition is 1, It uses the fit paramaters based on the quickavg wave
	
	
	string dataset = "cscurrent_2d"
	string quickavg = "dat" + num2str(wavenum) + "quickavg"
	string xvals = "dat" + num2str(wavenum) + "x_array"
	wave W_coef
	
	
	quick_avg2(wavenum, dataset, 1) // quick average plot
	fit_transition2($quickavg, 0)
	
	duplicate /o W_coef W_coef_guess //guess based on the quick average
	
	
//	plot_badthetas(wavenum, dataset, condition) // thetas vs repeat plot and bad theta sweep plot
//	centering(wavenum, dataset, condition) // centred plot and average plot
//	plot2d_heatmap(wavenum,dataset) // raw 2d plot

end