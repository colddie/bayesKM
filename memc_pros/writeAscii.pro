Pro writeAscii, parms, filename

    size  = n_elements(parms)
    parms = reform(parms,size)

    for i = 0, size-1 do begin
        openw,lun,filename,/get_lun,/append
        printf,lun,parms[i]
        close,lun 
        free_lun,lun
    endfor

End