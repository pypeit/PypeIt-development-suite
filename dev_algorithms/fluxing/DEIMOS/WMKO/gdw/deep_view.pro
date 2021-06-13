; $Id: xplot3d.pro,v 1.16 2000/09/06 22:43:43 paulcs Exp $
;
;  Copyright (c) 1997-2000, Research Systems, Inc. All rights reserved.
;       Unauthorized reproduction prohibited.
;
;+
; NAME:
;       xplot3d
;
; PURPOSE:
;       Provide an interactive version of IDL's plot_3dbox
;       command.
;
; CATEGORY:
;       Widgets, Object Graphics.  Plotting.
;
; CALLING SEQUENCE:
;       xplot3d, x, y, z
;
; REQUIRED INPUTS:
;       None.  Fake data will be used if no data is supplied in call.
;
; OPTIONAL INPUTS
;
;       x: A vector of X data values.
;
;       y: A vector of Y data values.
;
;       z: A vector of Z data values.
;
; OPTIONAL KEYWORD PARAMETERS:
;
;       GROUP: Widget group_leader.  When the group leader
;       is destroyed, this program will be destroyed.
;       MODAL: If set, block other IDL widgets from receiving events.
;       BLOCK: If set, block IDL command line.
;       TITLE: String.  A title for the plot.
;       [XYZ]TITLE: String.  Title for x, y or z axis.
;       NAME: A string.  This is a name for the data curve being plotted.
;           The name is displayed on xplot3d's toolbar if/when
;           the curve is selected with the mouse.  (To select
;           the curve with the mouse, xplot3d must
;           be in select mode.  You can put xplot3d in select mode
;           by clicking on the rightmost button in xplot3d's
;           toolbar.)
;       [XYZ]RANGE: Two element array specifying axis range.  [min, max]
;       OVERPLOT: Draw the curve in the most recently created view.
;           Title keywords, range keywords and the modal keyword are
;           ignored if we are overplotting.
;
;       LINESTYLE: Same as IDLgrPolyline::init keyword.
;       THICK: Same as IDLgrPolyline::init keyword.
;       COLOR: [r, g, b] triplet.  Color of curve.  Same as IDLgrPolyline.
;       SYMBOL: Same as IDLgrPolyline::init keyword.
;       TEST: If set, do not require input arguments.  Plot an example
;           sinusoidal curve instead.
;
; EXAMPLE 1:
;       IDL> x = indgen(20)
;       IDL> y = sin(x/3.)
;       IDL> z = x
;       IDL> xplot3d, x, y, z, thick=2, name='interesting curve'
;
; EXAMPLE 2:
;       IDL> x = indgen(20)
;       IDL> y = sin(x/3.)
;       IDL> z = x
;       IDL> oOrb = obj_new('orb', color=[0, 0, 255])
;       IDL> oOrb->Scale, .75, .1, .5
;       IDL> oSymbol = obj_new('IDLgrSymbol', oOrb)
;       IDL> xplot3d, x, y, z, thick=2, symbol=oSymbol
;
; KNOWN PROBLEMS:
;       Symbols (specified via the SYMBOL keyword) are not drawn in
;           a suitably muted color when projected (via the
;           "2D Projection..." menu functionality).
;       Heap variables may not be cleaned up properly if an error
;           occurs in procedure xplot3d.
;
;-
pro xplot3d__set_prop_on_children, oContainer, _extra=extra
compile_opt hidden
;
;Purpose: set properties on all immediate children contained
;in a given container.  The properties to be set are given by
;keyword.
;
oChildren = oContainer->Get(/all, count=count)
for i=0L,count-1 do begin
    oChildren[i]->SetProperty, _extra=extra
    endfor

end
;--------------------------------------------------------------------
pro xplot3d__update_button_bmps, state
compile_opt hidden

fname_prefixes = ['pan',       'rotate', 'zoom',  'select']
modes =          ['translate', 'rotate', 'scale', 'select']
wButtons = [state.wPan, state.wRotate, state.wZoom, state.wSelect]

indx = where(strlowcase(state.mode) eq modes)
indx = indx[0]

for i=0,n_elements(wButtons)-1 do begin
    filename = filepath( $
        fname_prefixes[i] $
            + (i eq indx ? '_active.bmp':'.bmp'), $
        subdirectory=state.bitmap_subdir $
        )
    widget_control, wButtons[i], /bitmap, set_value=filename
    endfor

end
;--------------------------------------------------------------------
pro xplot3d__free, var
compile_opt idl2, hidden

on_error, 2

case size(var, /tname) of
    'OBJREF': obj_destroy, var
    'POINTER': ptr_free, var
    'LONG': widget_control, var, /destroy
    else:
    endcase

end
;--------------------------------------------------------------------
pro xplot3d__draw, oWindow, oPicture, vector=vector
compile_opt hidden
;
;Purpose: On some platforms, when IDLgrWindow::Draw is invoked, math
;errors (e.g. "% Program caused arithmetic error: Floating illegal
;operand") are printed.  xplot3d__draw suppresses the
;printing of these errors.
;
;Flush and print any accumulated math errors
;
void = check_math(/print)
;
;Silently accumulate any subsequent math errors.
;
orig_except = !except
!except = 0
;
;Draw.
;
if n_elements(vector) gt 0 then begin
    oWindow->Draw, oPicture, vector=vector
    endif $
else begin
    oWindow->Draw, oPicture
    endelse
;
;Silently flush any accumulated math errors.
;
void = check_math()
;
;Restore original math error behavior.
;
!except = orig_except
end
;--------------------------------------------------------------------
pro xplot3d__set_axes_prop, state, _extra=extra

for i=0,n_elements(state.oAxes)-1 do $
    state.oAxes[i]->SetProperty, _extra=extra

end
;--------------------------------------------------------------------
pro xplot3d__set_backg_color, state, color
compile_opt hidden

on_error, 2

if size(color, /tname) ne 'STRING' then $
    message, 'second argument must be a string.'

case strlowcase(color) of
    'black': clr = [0, 0, 0]
    'gray': clr = [127, 127, 127]
    'white': clr = [255, 255, 255]
    endcase

xplot3d__set_prop_on_children, state.oScene, color=clr

widget_control, state.wBlackButton, sensitive=strlowcase(color) ne 'black'
widget_control, state.wWhiteButton, sensitive=strlowcase(color) ne 'white'
widget_control, state.wGrayButton,  sensitive=strlowcase(color) ne 'gray'

end
;--------------------------------------------------------------------
pro deep_view_event, event

widget_control, event.top, get_uvalue=pState
;
;If mouse buttons are down, only process draw widget events.
;
if (*pState).btndown eq 1 then begin
    if event.id eq (*pState).wDraw then begin
        if event.type eq 0 then begin ; Button down event?...
            return ; ...ignore it.  A mouse button is already down.
            end
        end $
    else begin
        return
        end
    end
;
;Handle resize events.
;
if tag_names(event, /structure_name) eq 'WIDGET_BASE' then begin
    widget_control, /hourglass
    toolbar_geom = widget_info((*pState).wToolbarBase, /geometry)
    fudge = 6
    xsize = event.x - fudge
    ysize = event.y - toolbar_geom.scr_ysize - fudge
    widget_control, (*pState).wDraw, xsize=xsize, ysize=ysize
    end
;
case event.id of
    (*pState).wReset: begin
        (*pState).oView2->ResetTransforms
        (*pState).oView2->SetViewVolume, $
            (*pState).oWindow, $
            /quiet, $
            /isotropic, $
            xrange=[0, 1], $
            yrange=[0, 1], $
            zrange=[0, 1]

        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end
    (*pState).wRotate: begin
        widget_control, (*pState).wLabel, set_value=' '

        (*pState).mode = 'rotate'
        xplot3d__update_button_bmps, *pState
        (*pState).oView2->SetProperty, mode=(*pState).mode
        end
    (*pState).wZoom: begin
        widget_control, (*pState).wLabel, set_value=' '

        (*pState).mode = 'scale'
        xplot3d__update_button_bmps, *pState
        (*pState).oView2->SetProperty, mode=(*pState).mode
        end
    (*pState).wPan: begin
        widget_control, (*pState).wLabel, set_value=' '

        (*pState).mode = 'translate'
        xplot3d__update_button_bmps, *pState
        (*pState).oView2->SetProperty, mode=(*pState).mode
        end
    (*pState).wSelect: begin

        (*pState).mode = 'select'
        xplot3d__update_button_bmps, *pState
        (*pState).oView2->SetProperty, mode=(*pState).mode
        end
    (*pState).wQuitButton : begin
        widget_control, event.top, /destroy
        return
        end
    (*pState).wVRMLButton : begin
        file = dialog_pickfile( $
            /write, $
            file='untitled.wrl', $
            group=event.top, $
            filter='*.wrl' $
            )
        if (file ne '') then begin
            widget_control, /hourglass
            (*pState).oWindow->GetProperty, $
                dimensions=dimensions, $
                resolution=resolution,$
                units=units, $
                color_model=color_model, $
                n_colors=n_colors
            oVRML = obj_new('IDLgrVRML', $
                dimensions=dimensions, $
                resolution=resolution, $
                units=units, $
                color_model=color_model, $
                n_colors=n_colors $
                )
            oVRML->SetProperty, filename=file
            xplot3d__draw, oVRML, (*pState).oView2
            obj_destroy,oVRML
            end
        end
    (*pState).wExportImageButton: begin
        geo = widget_info((*pState).wDraw, /geometry)
        oBuff = obj_new('IDLgrBuffer', $
            dimensions=[geo.scr_xsize, geo.scr_ysize] $
            )
        oBuff->Draw, (*pState).oScene
        oBuff->GetProperty, image_data=image_data
        void = dialog_write_image(image_data, dialog_parent=event.top)
        obj_destroy, oBuff
        end

    (*pState).wClipboardButton: begin
        (*pState).oWindow->GetProperty, $
            dimension=wdims, $
            resolution=res, $
            color_model=cm, $
            units=units, $
            n_colors=icolors
        oClipboard = obj_new('IDLgrClipboard', $
            dimensions=wdims, $
            resolution=res, $
            color_model=cm, $
            units=units, $
            n_colors=icolors $
            )
        xplot3d__draw, oClipboard, (*pState).oScene, /vector
        obj_destroy, oClipboard
        end
    (*pState).wPrintButton: begin
        oPrinter = obj_new('IDLgrPrinter')
        if dialog_printersetup( $
            oPrinter, $
            dialog_parent=event.top $
            ) $
        then begin
            if dialog_printjob( $
                oPrinter, $
                dialog_parent=event.top $
                ) $
            then begin
                widget_control, /hourglass
;
;               Convert views to inches.
;
                oViews = (*pState).oScene->Get(/all, count=count)

                default_arr = bytarr(count > 1)
                unit_arr = bytarr(count > 1)
                dim_arr = fltarr(2, count > 1)
                loc_arr = fltarr(2, count > 1)

                for i=0,count-1 do begin
                    oViews[i]->GetProperty, $
                        units=units, $
                        dimensions=dim, $
                        location=loc
                    unit_arr[i] = units
                    dim_arr[0, i] = dim
                    loc_arr[0, i] = loc

                    dim = oViews[i]->ViewportDimensions( $
                        (*pState).oWindow, $
                        location=loc, $
                        /inches, $
                        defaulted=defaulted $
                        )
                    default_arr[i] = defaulted

                    oViews[i]->SetProperty,$
                        location=loc, $
                        dimensions=dim, $
                        units=1
                    end
;
;               Print.
;
                xplot3d__draw, oPrinter, (*pState).oScene, /vector
                oPrinter->NewDocument
;
;               Restore views to their original units.
;
                for i=0,count-1 do begin
                    oViews[i]->SetProperty,$
                        location=default_arr[i] ? [0, 0] : loc_arr[*, i], $
                        dim=default_arr[i] ? [0, 0] : dim_arr[*, i], $
                        units=unit_arr[i]
                    endfor
                endif
            endif
        obj_destroy,oPrinter
        end
    (*pState).wBlackButton: begin
        xplot3d__set_backg_color, *pState, 'black'
        xplot3d__set_axes_prop, *pState, color=[255, 255, 255]
        if n_elements(*(*pState).pTitle) eq 1 then begin
            *(*pState).pTitle->SetProperty, color=[255, 255, 255]
            end
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end
    (*pState).wWhiteButton: begin
        xplot3d__set_backg_color, *pState, 'white'
        xplot3d__set_axes_prop, *pState, color=[0, 0, 0]
        if n_elements(*(*pState).pTitle) eq 1 then begin
            *(*pState).pTitle->SetProperty, color=[0, 0, 0]
            end
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end
    (*pState).wGrayButton: begin
        xplot3d__set_backg_color, *pState, 'gray'
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wAllShadow: begin
        (*pState).oShadowXYModel->SetProperty, hide=0
        (*pState).oShadowYZModel->SetProperty, hide=0
        (*pState).oShadowXZModel->SetProperty, hide=0
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wXYShadow: begin
        (*pState).oShadowXYModel->SetProperty, hide=0
        (*pState).oShadowYZModel->SetProperty, hide=1
        (*pState).oShadowXZModel->SetProperty, hide=1
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wYZShadow: begin
        (*pState).oShadowXYModel->SetProperty, hide=1
        (*pState).oShadowYZModel->SetProperty, hide=0
        (*pState).oShadowXZModel->SetProperty, hide=1
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wXZShadow: begin
        (*pState).oShadowXYModel->SetProperty, hide=1
        (*pState).oShadowYZModel->SetProperty, hide=1
        (*pState).oShadowXZModel->SetProperty, hide=0
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wNoShadow: begin
        (*pState).oShadowXYModel->SetProperty, hide=1
        (*pState).oShadowYZModel->SetProperty, hide=1
        (*pState).oShadowXZModel->SetProperty, hide=1
        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wBoxAxes: begin
        xplot3d__set_axes_prop, *pState, hide=0

        (*pState).oXaxis->SetProperty, ticklen=1.0
        (*pState).oYaxis->SetProperty, ticklen=1.0
        (*pState).oZaxis->SetProperty, ticklen=1.0

        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wSimpleAxes: begin
        xplot3d__set_axes_prop, *pState, hide=1

        (*pState).oXaxis->SetProperty, hide=0, ticklen=0
        (*pState).oYaxis->SetProperty, hide=0, ticklen=0
        (*pState).oZaxis->SetProperty, hide=0, ticklen=0

        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wNoAxes: begin
        xplot3d__set_axes_prop, *pState, hide=1

        xplot3d__draw, (*pState).oWindow, (*pState).oScene
        end

    (*pState).wDraw: begin
        case event.type of
            4 : begin ; Expose.
                xplot3d__draw, $
                    (*pState).oWindow, $
                    (*pState).oScene
                end
            0 : begin ; Button press
                case event.press of
                    4 : begin ; Right mouse-button.
                        end
                    2 : begin ; Middle mouse-button.
                        end
                    1 : begin ; Left mouse button.
                        oSelectedView = (*pState).oWindow->Select( $
                            (*pState).oScene, $
                            [event.x, event.y] $
                            )
                        oSelectedView = oSelectedView[0]
                        if obj_valid(oSelectedView) then begin
                            if oSelectedView eq (*pState).oView2 then begin
                                void = (*pState).oView2->Update(event)

                                if (*pState).mode eq 'select' then begin
                                    (*pState).oView2->GetProperty, $
                                        selected=oSelected

                                    if obj_valid(oSelected[0]) then begin
                                        oSelected[0]->GetProperty, name=name
                                        if name eq '' then $
                                            name = obj_class(oSelected[0])
                                        endif $
                                    else begin
                                        name = ''
                                        endelse

                                    widget_control, $
                                        (*pState).wLabel, $
                                        set_value=name
                                    endif $
                                else begin
                                    (*pState).we_are_manipulating = 1b
                                    endelse

                                end
                            end
                        end
                    else:
                    endcase
                end
            2 : begin ; Button motion.
                if (*pState).we_are_manipulating then begin
                    if (*pState).oView2->Update(event) then begin
                        xplot3d__draw, (*pState).oWindow, (*pState).oScene
                        endif
                    endif
                end
            1 : begin ; Button release.
                if (*pState).we_are_manipulating then begin
                    if (*pState).oView2->Update(event) then begin
                        xplot3d__draw, (*pState).oWindow, (*pState).oScene
                        endif
                    (*pState).we_are_manipulating = 0b
                    endif
                end
            else:
            endcase
        end
    else:
    endcase

end
;--------------------------------------------------------------------
Pro xplot3d_cleanup, wID
compile_opt idl2, hidden

widget_control, wID, get_uvalue=pState

if (*pState).group_leader_is_fabricated then begin
    if widget_info(*(*pState).pGroupLeader, /valid_id) then begin
        widget_control, *(*pState).pGroupLeader, /destroy
        end
    end
;
;Clean up heap variables.
;
for i=0,n_tags(*pState)-1 do begin
    case size((*pState).(i), /tname) of
        'POINTER': $
            ptr_free, (*pState).(i)
        'OBJREF': $
            obj_destroy, (*pState).(i)
        else:
        endcase
    end
ptr_free, pState
end
;--------------------------------------------------------------------
pro deep_view, x, y, z, $
    _extra=extra, $
    xtitle=xtitle, $
    ytitle=ytitle, $
    ztitle=ztitle, $
    title=title, $
    xrange=xrng, $
    yrange=yrng, $
    zrange=zrng, $
    overplot=overplot, $
    double_view=double_view, $
    group=group_leader, $   ; IN: (opt)
    block=block, $          ; IN: (opt) block command line.
    test=test, $
    modal=modal

on_error, 2 ; Return to caller on error.

common xplot3d, oView2, xs, ys, zs, xaxrange, yaxrange, zaxrange, $
    oShadowXYModel, oShadowXZModel, oShadowYZModel, oWindow, oScene, $
    oPolylineModel

widget_control, /hourglass

if keyword_set(test) then begin
    title = "Example Data"
    y = findgen(101)
    y = sin(y/5) / exp(y/50)
    x = findgen(n_elements(y))
    z = x
    name = 'Example Data'
    end

if n_elements(x) eq 0 or n_elements(y) eq 0 or n_elements(z) eq 0 then $
    message, 'Requires 3 non-keyword arguments (x, y and z).'

if n_elements(x) ne n_elements(y) $
or n_elements(x) ne n_elements(z) then $
    message, $
        'Invalid Input: X, Y, and Z should have the same number ' + $
        'of elements.'

oPolyline = obj_new('IDLgrPolyline',  $
    x, $
    y, $
    z, $
    _extra=extra, $
    name=name, $
    color=[255, 0, 0] $
    )

if (obj_valid(oView2) eq 0) or (keyword_set(overplot) eq 0) then begin
    doing_initial_plot = 1b

    oHelvetica14pt = obj_new('IDLgrFont', 'Helvetica', size=14)
    oHelvetica10pt = obj_new('IDLgrFont', 'Helvetica', size=10)

    if n_elements(xtitle) eq 0 then xtitle='X Axis'
    if n_elements(ytitle) eq 0 then ytitle='Y Axis'
    if n_elements(ztitle) eq 0 then ztitle=' '

    oScene = obj_new('IDLgrScene')

    if n_elements(title) gt 0 then begin
        oView1 = obj_new('IDLexInscribingView', $
            location=[0, .9], $
            dimensions=[1, .1], $
            units=3 $
            )
        oTitle = obj_new('IDLgrText', $
            title , $
            alignment=0.5, $
            font=oHelvetica14pt, $
            recompute=2 $
            )

        oModel = obj_new('IDLgrModel')
        oModel->Add, oTitle

        oView1->Add, oModel

        oView2 = obj_new('IDLexObjview', $
            location=[0, 0], $
            dimensions=[1, .9], $
            double=double_view, $
            units=3 $
            )
        oScene->Add, [oView1, oView2]
        end $
    else begin
        oView2 = obj_new('IDLexObjview', double=double_view)
        oScene->Add, oView2
        end

    oXtitle = obj_new('IDLgrText', xtitle, recompute=0)
    oYtitle = obj_new('IDLgrText', ytitle, recompute=0)
    oZtitle = obj_new('IDLgrText', ztitle, recompute=0)

    oXaxis = obj_new('IDLgrAxis', $
        0, $
        ticklen=0.1, $
        minor=0, $
        title=oXtitle, $
        /exact $
        )
    oXaxis->GetProperty, ticktext=oXaxisText
    oXaxisText->SetProperty, font=oHelvetica10pt

    oXaxis2 = obj_new('IDLgrAxis', 0, minor=0, notext=1, /exact, ticklen=0)
    oXaxis3 = obj_new('IDLgrAxis', 0, minor=0, notext=1, /exact, ticklen=0)

    oYaxis = obj_new('IDLgrAxis', $
        1, $
        ticklen=0.1, $
        minor=0, $
        title=oYtitle, $
        /exact $
        )
    oYaxis->GetProperty, ticktext=oYaxisText
    oYaxisText->SetProperty, font=oHelvetica10pt

    oYaxis2 = obj_new('IDLgrAxis', 1, minor=0, notext=1, /exact, ticklen=0)
    oYaxis3 = obj_new('IDLgrAxis', 1, minor=0, notext=1, /exact, ticklen=0)

    oZaxis = obj_new('IDLgrAxis', $
        2, $
        ticklen=0.1, $
        minor=0, $
        title=oZtitle, $
        /exact $
        )
    oZaxis->GetProperty, ticktext=oZaxisText
    oZaxisText->SetProperty, font=oHelvetica10pt

    oZaxis2 = obj_new('IDLgrAxis', 2, minor=0, notext=1, /exact, ticklen=0)
    oZaxis3 = obj_new('IDLgrAxis', 2, minor=0, notext=1, /exact, ticklen=0)

    dummy = -1
    if n_elements(xrng) eq 0 then $
        xrange = float([min(x), max(x)]) $
    else $
        xrange = float(xrng)

    if n_elements(yrng) eq 0 then $
        yrange = float([min(y), max(y)]) $
    else $
        yrange = float(yrng)

    if n_elements(zrng) eq 0 then $
        zrange = float([min(z), max(z)]) $
    else $
        zrange = float(zrng)

    if xrange[0] ge xrange[1] then begin
        xrange[1] = xrange[0] + 1
        end
    xs = norm_coord(xrange)
    oXaxis->SetProperty, range=xrange,  location=[dummy, 0, 0], $
        xcoord_conv=xs
    oXaxis2->SetProperty, range=xrange, location=[dummy, 1, 0], $
        xcoord_conv=xs
    oXaxis3->SetProperty, range=xrange, location=[dummy, 1, 1], $
        xcoord_conv=xs

    if yrange[0] ge yrange[1] then begin
        yrange[1] = yrange[0] + 1
        end
    ys=norm_coord(yrange)
    oYaxis->SetProperty, range=yrange,  location=[0, dummy, 0], $
        ycoord_conv=ys
    oYaxis2->SetProperty, range=yrange, location=[1, dummy, 0], $
        ycoord_conv=ys
    oYaxis3->SetProperty, range=yrange, location=[1, dummy, 1], $
        ycoord_conv=ys

    if zrange[0] ge zrange[1] then begin
        zrange[1] = zrange[0] + 1
        end
    zs=norm_coord(zrange)
    oZaxis->SetProperty, range=zrange,  location=[0, 1, dummy], $
        zcoord_conv=zs
    oZaxis2->SetProperty, range=zrange, location=[1, 0, dummy], $
        zcoord_conv=zs
    oZaxis3->SetProperty, range=zrange, location=[1, 1, dummy], $
        zcoord_conv=zs

    oXaxis->GetProperty, xrange=xaxrange, tickvalues=xtickvalues
    oYaxis->GetProperty, yrange=yaxrange, tickvalues=ytickvalues
    oZaxis->GetProperty, zrange=zaxrange, tickvalues=ztickvalues
;
;   Create gridlines.
;
    zNumlines = n_elements(ztickvalues)
    oGridz = objarr(zNumlines)
    for i=0,zNumlines-1 do begin
        oGridz[i]=obj_new('IDLgrPolyline', $
            [xaxrange[1], xaxrange[1]], $
            [yaxrange[1], yaxrange[0]],$
            [ztickvalues[i], ztickvalues[i]], $
            xcoord_conv=xs, $
            ycoord_conv=ys, $
            zcoord_conv=zs, $
            linestyle=([6, 1])[ $
                ztickvalues[i] ne zaxrange[0] and $
                ztickvalues[i] ne zaxrange[1] $
                ] $
            )
        endfor

    xNumlines = n_elements(xtickvalues)
    oGridx = objarr(xNumlines)
    for i=0,xNumlines-1 do begin
        oGridx[i]=obj_new('IDLgrPolyline', $
            [xtickvalues[i], xtickvalues[i]], $
            [yaxrange[1], yaxrange[1]],$
            [zaxrange[0], zaxrange[1]], $
            xcoord_conv=xs, $
            ycoord_conv=ys, $
            zcoord_conv=zs, $
            linestyle=([6, 1])[ $
                xtickvalues[i] ne xaxrange[0] and $
                xtickvalues[i] ne xaxrange[1] $
                ] $
            )
        endfor

    yNumlines = n_elements(ytickvalues)
    oGridy = objarr(yNumlines)
    for i=0,yNumlines-1 do begin
        oGridy[i]=obj_new('IDLgrPolyline', $
            [xaxrange[1], xaxrange[1]], $
            [ytickvalues[i], ytickvalues[i]],$
            [zaxrange[1], zaxrange[0]], $
            xcoord_conv=xs, $
            ycoord_conv=ys, $
            zcoord_conv=zs, $
            linestyle=([6, 1])[ $
                ytickvalues[i] ne yaxrange[0] and $
                ytickvalues[i] ne yaxrange[1] $
                ] $
            )
        endfor
;
    oAxes = [ $
        oXaxis, $
        oXaxis2, $
        oXaxis3, $

        oYaxis, $
        oYaxis2, $
        oYaxis3, $

        oZaxis, $
        oZaxis2, $
        oZaxis3, $

        oGridx, $
        oGridy, $
        oGridz $
        ]

    oPolylineModel = obj_new('IDLgrModel')
    oShadowXYModel = obj_new('IDLgrModel', /hide)
    oShadowXZModel = obj_new('IDLgrModel', /hide)
    oShadowYZModel = obj_new('IDLgrModel', /hide)
    end

shadow_color = [150, 150, 150]
shadow=intarr(n_elements(z))
oShadowXY = obj_new('IDLgrPolyline',  $
    x, $
    y, $
    shadow+zaxrange[0], $
    _extra=extra $
    )
oShadowXY->SetProperty, color=shadow_color

shadow=intarr(n_elements(y))
oShadowXZ = obj_new('IDLgrPolyline',  $
    x, $
    shadow+yaxrange[1], $
    z, $
    _extra=extra $
    )
oShadowXZ->SetProperty, color=shadow_color

shadow=intarr(n_elements(x))
oShadowYZ = obj_new('IDLgrPolyline',  $
    shadow+xaxrange[1], $
    y, $
    z, $
    _extra=extra $
    )
oShadowYZ->SetProperty, color=shadow_color

oPolylineModel->Add, oPolyline
oShadowXYModel->Add, oShadowXY
oShadowXZModel->Add, oShadowXZ
oShadowYZModel->Add, oShadowYZ

if keyword_set(doing_initial_plot) then begin
    oModel = obj_new('IDLgrModel')
    oModel->Translate, -.5, -.5, -.5
    oModel->Rotate, [0,0,1],  30
    oModel->Rotate, [1,0,0], -60
    oModel->Translate, .5, .5, .5
    oModel->Add, [ $
        oPolylineModel, $
        oShadowXYModel, $
        oShadowYZModel, $
        oShadowXZModel, $
        oAxes $
        ]

    oView2->SetProperty, subject=oModel
    end

oPolyline->SetProperty,xcoord_conv=xs, ycoord_conv=ys, zcoord_conv=zs
oShadowXY->SetProperty,xcoord_conv=xs, ycoord_conv=ys, zcoord_conv=zs
oShadowYZ->SetProperty,xcoord_conv=xs, ycoord_conv=ys, zcoord_conv=zs
oShadowXZ->SetProperty,xcoord_conv=xs, ycoord_conv=ys, zcoord_conv=zs

if keyword_set(doing_initial_plot) then begin
    oXaxis -> SetProperty, ticklen=1, gridstyle=1
    oYaxis -> SetProperty, ticklen=1, gridstyle=1
    oZaxis -> SetProperty, ticklen=1, gridstyle=1

    if n_elements(group_leader) ne 0 then begin
        if not widget_info(group_leader, /valid_id) then begin
            message, 'Specified Group Leader is not valid.'
            end
        end $
    else begin
        if keyword_set(modal) then begin
;
;           Modal widgets require a group leader.  A group leader was not
;           specified, so fabricate an invisible one.

            group_leader = widget_base(map=0)
            group_leader_is_fabricated = 1b
            end
        end

    if keyword_set(modal) then begin
        tlb = widget_base( $
            title='XPlot3D', $
            column=1, $
            tlb_size_events=1, $
            /modal, $
            group_leader=group_leader $
            )
        mbar = widget_base(tlb, /row, /frame)
        end $
    else begin
        tlb = widget_base( $
            title='XPlot3D', $
            column=1, $
            tlb_size_events=1, $
            group_leader=group_leader, $
            mbar=mbar $
            )
        end
;
;   Create the menu bar.
;
    wFileButton = widget_button(mbar, value='File', /menu)
        wExportImageButton = widget_button( $
            wFileButton, $
            value='Export Image...' $
            )
        wVRMLButton = widget_button( $
            wFileButton, $
            value="Export VRML..." $
            )
        wPrintButton = widget_button(wFileButton, value='Print...')
        wQuitButton = widget_button(wFileButton, /separator, value='Quit')

    wEditButton = widget_button(mbar, /menu, value='Edit')
        wClipboardButton = widget_button( $
            wEditButton, $
            value='Copy Scene to Clipboard' $
            )

    wViewButton = widget_button(mbar, value='View', /menu)
        wBackgroundButton = widget_button( $
            wViewButton, $
            value='Background', $
            /menu $
            )
            wBlackButton = widget_button(wBackgroundButton, value='Black')
            wWhiteButton = widget_button( $
                wBackgroundButton, $
                value='White', $
                sensitive=0 $
                )
            wGrayButton = widget_button(wBackgroundButton, value='Gray')

    if lmgr(/demo) then begin
        widget_control, wPrintButton, sensitive=0
        widget_control, wExportImageButton, sensitive=0
        widget_control, wClipboardButton, sensitive=0
        widget_control, wVRMLButton, sensitive=0
        end
;
;   Create toolbar.
;
    bitmap_subdir = ['resource', 'bitmaps']
    wToolbarBase = widget_base(tlb, /row, /frame)
        wReset = widget_button( $
            wToolbarBase, $
            value=filepath('reset.bmp', subdirectory=bitmap_subdir), $
            /bitmap $
            )
        void = widget_base(wToolbarBase, /row, scr_xsize=8)
        wRotate = widget_button( $
            wToolbarBase, $
            value=filepath('rotate_active.bmp', subdirectory=bitmap_subdir), $
            /bitmap $
            )
        wZoom = widget_button( $
            wToolbarBase, $
            value=filepath('zoom.bmp', subdirectory=bitmap_subdir), $
            /bitmap $
            )
        wPan = widget_button( $
            wToolbarBase, $
            value=filepath('pan.bmp', subdirectory=bitmap_subdir), $
            /bitmap $
            )
        wSelect = widget_button( $
            wToolbarBase, $
            value=filepath('select.bmp', subdirectory=bitmap_subdir), $
            /bitmap $
            )
        wLabel = widget_label(wToolbarBase, value=' ', /dynamic_resize)

    wDraw = widget_draw(tlb, $
        xsize=400, $
        ysize=400, $
        graphics_level=2, $
        retain=0, $
        /expose_events, $
        /button_events, $
        /motion_events $
        )
;
;   Shadow controls.
;
    wProjection = widget_button(wViewButton, value='2D Projection', /menu)
    wAllShadow = widget_button(wProjection, Value='All On')
    wXYShadow = widget_button(wProjection, Value='XY Plane')
    wYZShadow = widget_button(wProjection, Value='YZ Plane')
    wXZShadow = widget_button(wProjection, Value='XZ Plane')
    wNoShadow = widget_button(wProjection, Value='All Off')
;
;   Axes controls.
;
    wAxes = widget_button(wViewButton, value='Axes', /menu)
    wSimpleAxes = widget_button(wAxes, value='Simple Axes')
    wBoxAxes = widget_button(wAxes, value='Box Axes')
    wNoAxes = widget_button(wAxes, value='No Axes')

    widget_control, tlb, /realize
;
;   Draw our scene in the window.
;
    widget_control, wDraw, Get_Value=oWindow

    if n_elements(oView1) eq 1 then begin
        oView1->SetViewVolume, oWindow, /quiet, /isotropic
        end
    oView2->SetViewVolume, oWindow, /isotropic, /quiet, $
        xrange=[0,1], $
        yrange=[0,1], $
        zrange=[0,1]

    oWindow->Draw, oScene

    oXaxis->GetProperty, ticktext=oTicktext
    oTickText->SetProperty, recompute_dimensions=0
    oYaxis->GetProperty, ticktext=oTicktext
    oTickText->SetProperty, recompute_dimensions=0
    oZaxis->GetProperty, ticktext=oTicktext
    oTickText->SetProperty, recompute_dimensions=0
;
    widget_control, tlb, set_uvalue=ptr_new({ $
        btndown: 0b, $
        pGroupLeader: ptr_new(group_leader), $
        group_leader_is_fabricated: $
            keyword_set(group_leader_is_fabricated), $
        mode: 'rotate', $ ; ; 'translate', 'scale', 'select'
        we_are_manipulating: 0b, $
        bitmap_subdir: bitmap_subdir, $
        wDraw: wDraw, $
        wQuitButton: wQuitButton, $
        wVRMLButton: wVRMLButton, $
        wExportImageButton: wExportImageButton, $
        wClipboardButton: wClipboardButton, $
        wPrintButton: wPrintButton, $
        wReset: wReset, $
        wRotate: wRotate, $
        wZoom: wZoom, $
        wPan: wPan, $
        wSelect: wSelect, $
        wLabel: wLabel, $
        wBlackButton: wBlackButton, $
        wWhiteButton: wWhiteButton, $
        wGrayButton: wGrayButton, $
        wAllShadow: wAllShadow, $
        wXYShadow: wXYShadow, $
        wYZShadow: wYZShadow, $
        wXZShadow: wXZShadow, $
        wNoShadow: wNoShadow, $
        wSimpleAxes: wSimpleAxes, $
        wBoxAxes: wBoxAxes, $
        wNoAxes: wNoAxes, $
        wToolbarBase: wToolbarBase, $
        oWindow: oWindow, $
        oModel: oModel, $
        oPolylineModel: oPolylineModel, $
        oShadowXYModel: oShadowXYModel, $  ; 2D projection.
        oShadowXZModel: oShadowXZModel, $  ; 2D projection.
        oShadowYZModel: oShadowYZModel, $  ; 2D projection.
        oXaxis: oXaxis, $
        oYaxis: oYaxis, $
        oZaxis: oZaxis, $
        oAxes: oAxes, $
        oHelvetica10pt: oHelvetica10pt, $
        oHelvetica14pt: oHelvetica14pt, $
        oXtitle: oXtitle, $
        oYtitle: oYtitle, $
        oZtitle: oZtitle, $
        pTitle: ptr_new(oTitle), $
        oScene: oScene, $
        oView2: oView2 $
        })

    xmanager, $
        'xplot3d', $
        tlb, $
        cleanup='xplot3d_cleanup', $
        no_block=keyword_set(block) eq 0
    end $
else begin
    xplot3d__draw, oWindow, oScene
    end

if keyword_set(group_leader_is_fabricated) then begin
;
;   Leave GROUP_LEADER parameter like we found it: undefined.
;
    ptr_free, ptr_new(group_leader, /no_copy)
    end

end

