pskeywords = psconfig_gdw( cancelled=cancelled)
; if not cancelled then begin
if pskeywords.toprinter then pskeywords.filename = getenv("TMPDIR") + '/idltemp' + strtrim(string(long(systime(1))),2) + ".ps"
thisDevice = !D.Name
set_plot, 'PS'
device, _EXTRA=pskeywords
plot, findgen(11)
device, /close
Set_Plot, thisDevice
if pskeywords.toprinter then begin
    command = pskeywords.command + ' ' + pskeywords.filename
    spawn, command, exit=status
    file_delete, pskeywords.filename
endif

end
