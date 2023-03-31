#pragma TextEncoding = "UTF-8"
#pragma rtGlobals=3				// Use modern global access method and strict wave access
#pragma DefaultTab={3,20,4}		// Set default tab width in Igor Pro 9 and later

#include <Reduce Matrix Size>


function /wave quick_avg2(int wavenum, string dataset, int view) // /WAVE lets your return a wave

//  averaging over total sweeps without centering or removing any data
// wave returned is named "dat" + "wavenum" + "quickavg" eg. dat3320quickavg


//    Inputs: wavenum is the dat number
//            dataset implies either "cscurrent_2d" or "dotcurrent_2d"
//            view if set to 1, will display the quick_avg2 plot

	
	variable i
	string w2d
	string avg_name
	int nc 
	int nr 
	
	w2d="dat"+num2str(wavenum)+dataset //current 2d array
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
		duplicate /o W_coef_guess, W_coef
	endif
	
	//smoothing and differentiating //could be useful for amp guess and linear term guess
	
	Duplicate /o current_array smooth_current;DelayUpdate
	Smooth/S=4 51, smooth_current;DelayUpdate                           //double check all the smoothing
	FuncFit Chargetransition W_coef smooth_current[][0] /D   //removed the x_array
	
end




function /wave get_fit_params2(int wavenum, string dataset, int condition)
	// returns wave with the name wave "dat"+ wavenum +"fit_params" eg. dat3320fit_params
	
	//If condition is 0 it will get initial params, If 1:
    // define a variable named W_coef_guess = {} with the correct number of arguments
	

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

	make /N= (nr , 10) /o $fit_params_name
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
      fit_params[1 * i][,4] = W_coef[q]
      fit_params[1 * i][5,] = W_sigma[q-5]         //I genuinely cant believe this worked 
		    											// i dont think the q-5 does anything, should double check
	endfor
	
	return fit_params
	
end


function plot_thetas2(int wavenum, string dataset, int condition)

	//If condition is 0 it will get initial params, If 1:
    // define a variable named W_coef_guess = {} with the correct number of arguments
	
	string fit_params_name = "dat"+num2str(wavenum)+"fit_params"
	variable thetamean
	variable thetastd
	variable i
	int nr
	
	
	get_fit_params2(wavenum, dataset, condition)
	wave fit_params = $fit_params_name
	nr = dimsize(fit_params,0)
	
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



function plot_badthetas2(int wavenum, string dataset, int condition)

	int i 
	int nr
	wave badthetasx
	string w2d
	
	w2d= "dat"+num2str(wavenum)+dataset //current 2d array
	
	wave wavenm = $w2d
	duplicate /o wavenm, wavenmcopy
	
	plot_thetas2(wavenum, dataset, condition)
	nr = dimsize(badthetasx,0)
	
	
	if (dimsize(wavenmcopy,1)<151)
		matrixtranspose wavenmcopy
	endif
	
	display
	
	for(i=0; i < nr; i +=1)
		appendtograph wavenmcopy[badthetasx[i]][]
	endfor

	QuickColorSpectrum2()
	
	ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
	Label bottom "voltage"
    Label left dataset
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 bad thetas of dat" + num2str(wavenum)
	
end 



function centering2(int wavenum, string dataset, int condition)

	string fit_params_name = "dat"+num2str(wavenum)+"fit_params"
	string centered_avg_name = "dat"+num2str(wavenum)+"centavg"
	wave avg_current
	wave fit_smooth_current
	wave w_coef
	wave w_sigma
	wave goodthetasx
	string w2d
	int i
	int nr
	int nc
	int nnr
	
	
	wave fit_params = $fit_params_name
	w2d= "dat"+num2str(wavenum)+dataset //current 2d array
	wave waved = $w2d
	
	get_fit_params2(wavenum, dataset, condition)
	
	duplicate /o $w2d wavecopy
	duplicate /O $w2d centered_2dx
	
	if (dimsize(centered_2dx,1)<151)
		matrixtranspose centered_2dx
		matrixtranspose wavecopy
	endif
	
	duplicate /o /r = [0][] wavecopy wavex
	matrixtranspose wavex
	wavex = x
	
	
	
	nr = dimsize(centered_2dx,0)
	nc = dimsize(centered_2dx,1)
	nnr = dimsize(goodthetasx,0)
	
	centered_2dx = 0
	
	duplicate /o/r = [0,nr][3] fit_params mids
	duplicate /o/r =[nc/10, nc - nc/10] wavex new_x
	
	
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
	copyscales wavex new2dwave
	
	appendimage new2dwave //append image of data
	ModifyImage new2dwave ctab= {*,*,Turbo,0} //setting color (idk why it prefers the pointer)
	ColorScale /A=RC /E width=20 //puts it on the right centre, /E places it outside

    Label left "repeats"
    //Label left dataset

    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 Centred good thetas of dat" + num2str(wavenum)
    
     
    avg_wave2(new2dwave) //centred and averaged 2D data, returns a wave called avg_current
    
    duplicate /o avg_current $centered_avg_name 
    ModifyGraph mode=2,lsize=3,rgb=(48059,48059,48059)
    
    
    fit_transition2(avg_current, 1) // get fit transition
    make /o/n=(dimsize(new_x,0)) fit
    //fit = w_coef[0]*tanh((new_x - w_coef[3])/(2*w_coef[2])) + w_coef[4]*new_x + w_coef[1] //theres already a function, I dont need to make it like this.
    //appendToGraph fit vs new_x
    appendToGraph fit_smooth_current
    Legend/C/N=text0/J/E=2 "\\s(avg_current) average\r\\s(fit) fit"
    TextBox/C/N=text1/A=MT/E=2 "\\Z14\\Z16 dat" + num2str(wavenum) + " average"
    TextBox/C/N=text2/A=MC/E=2 "\\Z14\\Z16 theta = " + num2str(w_coef[2]) + "+/-" +  num2str(W_sigma[2])

end


function /wave avg_wave2(wave waved)

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
	
	matrixtranspose wavenm
	copyscales wavenm avg_current
	
	display avg_current[][0]
	Label bottom "voltage"
    Label left "current"
    ModifyGraph fSize=24
    ModifyGraph gFont="Gill Sans Light"
    ModifyGraph width={Aspect,1.62},height=300
    
   
	
	return avg_current
	
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
	
	
	
	plot_badthetas2(wavenum, dataset, condition) // thetas vs repeat plot and bad theta sweep plot
	centering2(wavenum, dataset, condition) // centred plot and average plot
	plot2d_heatmap2(wavenum,dataset) // raw 2d plot
	
	

end







//function chargetransition_procedure2m(int wavenum, int condition)
//
//	chargetransition_procedure2(wavenum, condition)
//	MultiGraphLayout(WinList("*", ";", "WIN:1"), 3, 20, "AllGraphLayout")
//
//end


//from: https://www.wavemetrics.com/forum/igor-pro-wish-list/automatically-color-traces-multi-trace-graph

Function QuickColorSpectrum2()                            // colors traces with 12 different colors
    String Traces    = TraceNameList("",";",1)               // get all the traces from the graph
    Variable Items   = ItemsInList(Traces)                   // count the traces
    Make/FREE/N=(11,3) colors = {{65280,0,0}, {65280,43520,0}, {0,65280,0}, {0,52224,0}, {0,65280,65280}, {0,43520,65280}, {0,15872,65280}, {65280,16384,55552}, {36864,14592,58880}, {0,0,0},{26112,26112,26112}}
    Variable i
    for (i = 0; i <DimSize(colors,1); i += 1)
        ModifyGraph rgb($StringFromList(i,Traces))=(colors[0][i],colors[1][i],colors[2][i])      // set new color offset
    endfor
End







//from:
// https://www.wavemetrics.com/code-snippet/stacked-plots-multiple-plots-layout

function MultiGraphLayout(GraphList, nCols, spacing, layoutName)
    string GraphList        // semicolon separated list of graphs to be appended to layout
    variable nCols      // number of graph columns
    string layoutName   // name of the layout
    variable spacing        // spacing between graphs in points!
   
    // how many graphs are there and how many rows are required
    variable nGraphs = ItemsInList(GraphList)
    variable nRows = ceil(nGraphs / nCols)
    variable LayoutWidth, LayoutHeight 
    variable gWidth, gHeight
    variable maxWidth = 0, maxHeight = 0
    variable left, top
    variable i, j, n = 0

    string ThisGraph

    // detect total layout size from individual graph sizes; get maximum graph size as column/row size
    for(i=0; i<nGraphs; i+=1)
       
        ThisGraph = StringFromList(i, GraphList)
        GetWindow $ThisGraph gsize
        gWidth = (V_right - V_left)
        gHeight = (V_bottom - V_top)
       
        // update maximum
        maxWidth = gWidth > maxWidth ? gWidth : maxWidth
        maxHeight = gHeight > maxHeight ? gHeight : maxHeight  
    endfor
   
    // calculate layout size
    LayoutWidth = maxWidth * nCols + ((nCols + 1) * spacing)
    LayoutHeight = maxHeight * nRows + ((nRows +1) * spacing)
   
    // make layout; kill if it exists
    DoWindow $layoutName
    if(V_flag)
        KillWindow $layoutName
    endif
   
    NewLayout/N=$layoutName/K=1/W=(517,55,1451,800)
    LayoutPageAction size=(LayoutWidth, LayoutHeight), margins=(0,0,0,0)
    ModifyLayout mag=0.75
   
    //append graphs
    top = spacing
    for(i=0; i<nRows; i+=1)
   
        // reset vertical position for each column
        left = spacing
       
        for (j=0; j<    nCols; j+=1)
       
            ThisGraph = StringFromList(n, GraphList)
            if(strlen(ThisGraph) == 0)
                return 0
            endif
           
            GetWindow $ThisGraph gsize
            gWidth = (V_right - V_left)
            gHeight = (V_bottom - V_top)
           
            AppendLayoutObject/F=0 /D=1 /R=(left, top, (left + gWidth), (top + gHeight)) graph $ThisGraph
       
            // shift next starting positions to the right
            left += maxWidth + spacing
           
            // increase plot counter
            n += 1             
        endfor  
       
        // shift next starting positions dwon
        top += maxHeight + spacing
    endfor
   
    return 1
end



// https://www.wavemetrics.com/code-snippet/stacked-plots-multiple-plots-graph


function MakeStackedGraph(yWaveList, xWaveList, nCols, spacing, GraphName, mirror)
    string yWaveList    // semicolon separated list containing wave names to be plotted
    string xWaveList    // semicolon separated list containg corresonding x data; if list items are empty, y-waves will be plotted against x-scaling
    variable nCols      // number of columns within the stacked graph
    variable spacing    // spacing between plots in terms of fraction of total plot area
    string GraphName    // name of the stacked graph
    variable mirror     // mirror axis on = 1, or off = 0
   
    variable nGraphs
    variable nRows
    variable nGaps
    variable yLength
    variable xLength
    variable y0, x0, y1, x1
    variable i, j, n = 0
   
    string yWave, xWave
    string yAxisName, xAxisName, yMirrorName, xMirrorName
       
    // how many graphs are there and how many rows are required
    nGraphs = ItemsInList(yWaveList)
    nRows = ceil(nGraphs / nCols)
   
    // calculate length of axis from given spacing
    nGaps = nCols - 1
    xLength = (1 - spacing * nGaps) / nCols    
   
    nGaps = nRows - 1
    yLength = (1 - spacing * nGaps) / nRows
   
    // Display empty window; kill if it exists
    DoWindow/F $GraphName
    if(V_flag)
        KillWindow $GraphName
    endif
    Display/K=1 /N= $Graphname 
   
    // append traces
    for(i=0; i<nRows; i+=1)
   
        // reset vertical axis position
        y0 = 0
        x1 = 0

        for (j=0; j<    nCols; j+=1)
           
            // get wave names from lists
            yWave = StringFromList(n, yWaveList)
            yAxisName = "yAxis" + num2str(n)               
            xWave = StringFromList(n, xWaveList)
            xAxisName = "xAxis" + num2str(n)
           
            if(strlen(xWave) == 0)
                // if x-string is empty; plot against wave scaling
                AppendToGraph/L=$yAxisName /B=$xAxisName $yWave
            else
                AppendToGraph/L=$yAxisName /B=$xAxisName $yWave vs $xWave
            endif
           
            // set lower left position of y and x axis
            ModifyGraph freePos($yAxisName)={y0,kwFraction}
            ModifyGraph freePos($xAxisName)={x0,kwFraction}
           
            // set length of the axis
            ModifyGraph axisEnab($yAxisName)={y1,(y1+yLength)}
            ModifyGraph axisEnab($xAxisName)={x1,(x1+xLength)}
           
            // do some formatting
            ModifyGraph tick($yAxisName) = 2, tick($xAxisName) = 2
            ModifyGraph btlen($yAxisName) = 4, btlen($xAxisName) = 4
       
            // append mirror axis
            if(mirror)
                yMirrorName = "yMirror" + num2str(n)
                NewFreeAxis/L $yMirrorName
                ModifyFreeAxis $yMirrorName master = $yAxisName
                ModifyGraph freePos($yMirrorName) ={(y0+xLength), kwFraction}
                ModifyGraph axisEnab($yMirrorName) ={y1, (y1+ylength)}
                ModifyGraph noLabel($yMirrorName)=2
                ModifyGraph tick($yMirrorName) = 0, btlen($yMirrorName) = 4
               
                xMirrorName = "xMirror" + num2str(n)
                NewFreeAxis/B $xMirrorName
                ModifyFreeAxis $xMirrorName master = $xAxisName
                ModifyGraph freePos($xMirrorName) ={(x0+yLength), kwFraction}
                ModifyGraph axisEnab($xMirrorName) ={x1, (x1+xlength)}
                ModifyGraph noLabel($xMirrorName)=2
                ModifyGraph tick($xMirrorName) = 0, btlen($xMirrorName) = 4
            endif
           
            // shift next starting positions to the right
            y0 += xLength + spacing
            x1 += xlength + spacing
           
            // increase plot counter
            n += 1             
        endfor  
       
        // shift next starting positions up
        x0 += yLength + spacing
        y1 += yLength + spacing
    endfor
   
    return 1
end







/////// Dealing Interlacing ////////



// improvements on this function
//			let it take an argument of names for all the waves created


//     		an option or a new function all together that seperates the waves by grouping
//									i.e grouping x number of rows in a m by n matrix creating 
//                                      a total of m/x waves, also takes an argument for naming?
//          it could group based on amount of splits e.g split 2D wave into 4
//          it could group based on number of rows indicated. 


function periodic_seperation(string interlaced, int period)
	//interlaced = wavename
	// output = "interlaced" + "period" + "0"
	
	int i
	int j

	duplicate /o $interlaced interlaced_c
	
	if (dimsize(interlaced_c,1)<151)    //this might not be great depending on the matrix size
		matrixtranspose interlaced_c    //ask tim/johann for an example dat to test it out?
	endif
	
	
	
	int nr = dimsize(interlaced_c,1) //number of rows (total sweeps)
	int nc = dimsize(interlaced_c,0) //number of columns (data points)
	
	
	// makes waves based on the period //
	// make list of all new string names? // create this as an option to include its own names?
	// loop through a name, point to that wave, assign the corresponding row
	
	
	make /t /o /n = (period) wave_names   //wavename list (text wave)
	
	for (i=0; i < period ; i+=1)
	
		string wave_name = interlaced + "period" + num2str(i)
		wave_names[i] = wave_name
		make /n=(nr/period, nc) /o $wave_name
		
	endfor
	
	
	for (i=0; i < nr/period ; i+=1)	
		for (j=0; j < period ; j+=1)
		
			wave tempwave = $(wave_names[j]) 
			tempwave[1*i] = interlaced_c[q][(period * i + j)]
		
		endfor
	endfor

end





////$somewave + "_somestring" // this is wrong
////$(somewave + "_somestring") // this is correct


