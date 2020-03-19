PRO Error_Estimate_Plot

   ; Example data.
   data = cgDemoData(1)
   time = cgScaleVector(Findgen(N_Elements(data)), 0, 6)
   high_error = (data + cgScaleVector(RandomU(seed, N_Elements(data)), 3, 7)) < 35
   low_error = (data - cgScaleVector(RandomU(seed, N_Elements(data)), 2, 6)) > (-5)
   
   ; Set up variables for the plot. Normally, these values would be 
   ; passed into the program as positional and keyword parameters.
   xtitle = 'Time'
   ytitle = 'Signal Strength'
   title = 'Error Estimate Plot'
   position = [0.125, 0.125, 0.9, 0.925]
   
   ; Set up a "window" for the plot. The PostScript output will have
   ; the same aspect ratio as the graphics window on the display.
   cgDisplay, 600, 500, Title='Error Estimate Plot'
      
   ; Draw the line plot with no data
   cgPlot, time, data, Title=title, XTitle=xtitle, YTitle=ytitle, $
      Position=position, /NoData, YRange=[-5, 35], YStyle=1
      
   ; Fill in the error estimates.
   cgColorFill, [time, Reverse(time), time[0]], $
       [high_error, Reverse(low_error), high_error[0]], $
       Color='sky blue'

   ; Draw the line plot with no data
   cgPlotS, time, data, Color='red', PSym=-16, SymColor='olive', $
      SymSize=1.0, Thick=2
   stop
END ;*****************************************************************

; ; This main program shows how to call the program and produce
; ; various types of output.

;   ; Display the line plot in a graphics window.
;   cgDisplay, 600, 500
;   Error_Estimate_Plot

  
;   ; Display the line plot in a resizeable graphics window.
;   cgWindow, 'Error_Estimate_Plot'
  
;   ; Create a PostScript file.
;   cgPS_Open, Filename='error_estimate_plot.ps'
;   Error_Estimate_Plot
;   cgPS_Close, /PNG

; END