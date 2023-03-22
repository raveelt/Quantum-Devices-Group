#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

function /wave avg_wave(wave waved, wave wavex)

// averaging 2D wave down to 1D
	
	variable i
	int nc 
	int nr 
	
	duplicate /o waved wavenm
	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)
	
	make /n=(nc, 1 ) /o avg_current
	avg_current = wavenm[0][p] //first row
	
	for (i=1; i < nr ; i+=1)
		avg_current += wavenm[i][p] //sums all the rows
	endfor
	
	avg_current = avg_current / nr // divide by total rows
	
	display avg_current[][0] vs wavex
	Label bottom "voltage"
    Label left "current"
    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    
   
	
	return avg_current
	
end




function /wave quick_avg(int wavenum, string dataset) // /WAVE lets your return a wave

// averaging over total sweeps without centering or removing any data
	
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
    Label left "current"
    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    ModifyGraph mode=2,lsize=2,rgb=(21845,21845,21845)
	Legend/C/N=text0/J/A=RT "\\Z14\\Z16\\s(avgcurrent) quick average of dat" + num2str(wavenum)
	
	return avgcurrent
	
end


	
function plot2d_heatmap(num,dataset)

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



function /wave get_initial_params(sweep, x_array)

// for a given sweep returns a guess of initial parameters for the fit function: Charge transiton

	wave sweep 
	wave x_array
	
		
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
	
	//killwaves extractedwave, sweepsmooth
	return W_coef

end

	
	
	
//function to perform a fit

function /wave fit_transition(current_array,x_array, condition)
	
	wave current_array
	wave x_array
	int condition
	wave avgcurrent
	
	wave W_coef
	
	if(condition == 0)
		get_initial_params(current_array, x_array)
	else
		//get_initial_params(avgcurrent, x_array)
		fit_transition(avgcurrent, x_array, 0) // probably would be faster to calculate it once then duplicate each time.
	endif
	
	//smoothing and differentiating //could be useful for amp guess and linear term guess
	
	Duplicate/O current_array, smooth_current;DelayUpdate
	Smooth/S=4 51, smooth_current;DelayUpdate                           //double check all the smoothing
	FuncFit Chargetransition W_coef smooth_current[][0] /X= x_array /D 
	
end
	


//function to perform all the fits

function /wave get_fit_params(int wavenum, string dataset, int condition)

	variable i
	string w2d
	string w2x
	int nc 
	int nr
	wave fit_params
	wave temp_wave
	wave W_coef
	wave W_sigma
	 

	w2d="dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	
	wave wavenm = $w2d
	
	if (dimsize(wavenm,1)<151)
		matrixtranspose wavenm
	endif
	
	nr = dimsize(wavenm,0) //number of rows (total sweeps)
	nc = dimsize(wavenm,1) //number of columns (data points)
	

	make /N= (2 * nr, 5) /O fit_params
	make /N= (nc) /O temp_wave
	  
	for (i=0; i < nr ; i+=1)
		//create wave for for each i'th row
		//run fit_transition(current_array, $w2x)
		//store in corresponding 2d array
		
	  temp_wave = wavenm[i][p]
      fit_transition(temp_wave, $w2x, condition)
      fit_params[1 * i] = W_coef[q]
      fit_params[(nr/2) + i] = W_sigma[q] 
		
	endfor
	
	return fit_params
	
end
	
	

function plot_thetas(int wavenum, string dataset, int condition)
	
	wave fit_params
	variable thetamean
	variable thetastd
	variable i
	int nr
	
	get_fit_params(wavenum, dataset, condition)
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
			
			insertPoints /v = (thetas[i]) nr, 1, goodthetas // value of theta
			insertpoints /v = (i) nr, 1, goodthetasx        // the repeat
			
		else
		
			insertPoints /v = (thetas[i]) nr, 1, badthetas // value of theta
			insertpoints /v = (i) nr, 1, badthetasx        // repeat
			
		endif
		
	endfor
			
		
	
	
	display meanwave, stdwave, stdwave2
	appendtograph goodthetas vs goodthetasx
	appendtograph badthetas vs badthetasx	
	
	
	ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    ModifyGraph lstyle(meanwave)=3,rgb(meanwave)=(17476,17476,17476)
    ModifyGraph lstyle(stdwave)=3,rgb(stdwave)=(52428,1,1)
    ModifyGraph lstyle(stdwave2)=3,rgb(stdwave2)=(52428,1,1)
    ModifyGraph mode(goodthetas)=3,lsize(goodthetas)=2, rgb(goodthetas)=(2,39321,1)
    ModifyGraph mode(badthetas)=3
    Legend/C/N=text0/J/A=RT "\\s(meanwave) mean\r\\s(stdwave) 2*std\r\\s(goodthetas) good\r\\s(badthetas) outliers"
	TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 thetas of dat" + num2str(wavenum)
    
    
    
    Label bottom "repeat"
    Label left "theta values"
    
    
    
	
end



function plot_badthetas(int wavenum, string dataset, int condition)

	int i 
	int nr
	wave badthetasx
	string w2d
	string w2x
	
	w2d= "dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	
	wave wavenm = $w2d
	duplicate /o wavenm, wavenmcopy
	
	plot_thetas(wavenum, dataset, condition)
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
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 bad thetas of dat" + num2str(wavenum)
	
end 


function centering(int wavenum, string dataset, int condition)

	wave fit_params
	wave avg_current
	wave w_coef
	wave w_sigma
	wave goodthetasx
	string w2d
	string w2x
	int i
	int nr
	int nc
	int nnr
	
	w2d= "dat"+num2str(wavenum)+dataset //current 2d array
	w2x = "dat"+num2str(wavenum)+"x_array" //voltage array
	
	wave wavex = $w2x
	wave waved = $w2d
	
	get_fit_params(wavenum, dataset, condition)
	
	duplicate /o $w2d wavecopy
	duplicate /O $w2d centered_2dx
	//duplicate /o $w2d new2dwave
	
	if (dimsize(centered_2dx,1)<151)
		matrixtranspose centered_2dx
		matrixtranspose wavecopy
	endif
	
	nr = dimsize(centered_2dx,0)
	nc = dimsize(centered_2dx,1)
	nnr = dimsize(goodthetasx,0)
	
	centered_2dx = 0
	
	duplicate /o/r = [0,nr][3] fit_params mids
	duplicate /o/r =[nc/10, nc - nc/10] $w2x new_x
	
	
	make /o/n = (nnr, (dimsize(new_x,0))) new2dwave
	
	
	for(i = 0; i < nnr; i += 1)
	
		duplicate /o wavex wavex2
		matrixtranspose wavex2
		
		wavex2 -= mids[goodthetasx[i]]
		wavex2 += mean(mids)
		centered_2dx[1 * goodthetasx[i]] = wavex2[q] //this forloop collects the centred x data
		
		duplicate /o/r=[goodthetasx[i]][0,nc] wavecopy sweep
		
		Interpolate2/T=2/E=2/Y=new_y/X=new_x/I=3 wavex2, sweep
		
		matrixtranspose new_y
		
		new2dwave[1*i] = new_y[q]
		 	
	endfor
	
	display; //start with empty graph
	
	matrixtranspose new2dwave
	
	appendimage new2dwave //append image of data
	ModifyImage new2dwave ctab= {*,*,Turbo,0} //setting color (idk why it prefers the pointer)
	ColorScale /A=RC /E width=20 //puts it on the right centre, /E places it outside

    Label left "repeats"
    //Label left dataset

    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 Centred good thetas of dat" + num2str(wavenum)
    
     
    avg_wave(new2dwave, new_x) //centred and averaged 2D data, returns a wave called avg_current
    ModifyGraph mode=2,lsize=3,rgb=(48059,48059,48059)
    
    
    fit_transition(avg_current, new_x, condition) // get fit transition
    make /o/n=(dimsize(new_x,0)) fit
    fit = w_coef[0]*tanh((new_x - w_coef[3])/(2*w_coef[2])) + w_coef[4]*new_x + w_coef[1] //theres already a function, I dont need to make it like this.
    appendToGraph fit vs new_x
    Legend/C/N=text0/J/E=2 "\\s(avg_current) average\r\\s(fit) fit"
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 dat" + num2str(wavenum) + " average"
    TextBox/C/N=text2/A=MC/E=2 "\\Z14\\Z16 theta = " + num2str(w_coef[2]) + "+/-" +  num2str(W_sigma[2])

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



function chargetransition_procedure(int wavenum, int condition)

	string dataset
	dataset = "cscurrent_2d"
	
	quick_avg(wavenum, dataset) // quick average plot
	plot_badthetas(wavenum, dataset, condition) // thetas vs repeat plot and bad theta sweep plot
	centering(wavenum, dataset, condition) // centred plot and average plot
	plot2d_heatmap(wavenum,dataset) // raw 2d plot

end


/////////////////// NOTES /////////////////////////////////////////////////////////////////////////////


// bounds for fit parameters
// setting amplitude by differentiation seems to not work
// lin term could maybe be set by picking portions of the start of the wave or the differentiation wave


//from meeting:
			//DONE// use fit parameters from quick avg for all sweeps // DONE, But i dont love the layout
			//DONE// fixed how the fit parameters are layed out W_coef then W_sigma
			

//		doesnt work for conductance peaks that are shifting (maybe have an on/off option?)

// 		    //DONE// all plots in one window

// 		    //DONE//fix x - scaling - suprisingly time consuming.

// 		deal with some interlacing (dividing current set into subsets of sweeps)
						//(would likely just need another function)

// 		     //DONE// keep all avg_current waves (name them accordingly)
							//dat3320quickavg
							//dat3320centavg
							//dat3320fit_params

// 		centering option (yes/no)


// other:
//       centered data is being halfed after x-scaling fix - need to check whats going on
//       multi-graph layout axis need to be adjusted


// waves created example: dat3320quickavg, dat3320fit_params, dat3320centavg




// run charge transition by doing the following:

//chargetransition_procedure2(3914, 1);
//MultiGraphLayout(WinList("*", ";", "WIN:1"), 3, 20, "AllGraphLayout");

//reduce matrix size? igor procedure file?
	//


///////////////////////////////////////////////////////////////////////////////////////////////////////////////