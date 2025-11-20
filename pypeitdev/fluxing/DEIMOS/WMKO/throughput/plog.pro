;;------------------------------------------------------------------------
pro plog, logunit, msg
;;------------------------------------------------------------------------
;; print a message to terminal and, with timestamp, to logfile
;;------------------------------------------------------------------------
message, msg, /info, level=-1
now = gdwtimestamp()
buf = '[' + now + '] '+msg
printf, logunit, buf
flush, logunit
end
