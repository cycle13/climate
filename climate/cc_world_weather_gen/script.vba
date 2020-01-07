Sub ProgressBar()

Dim FName, FNameDisplay, FParam, Datareturn, ChangePoint, AvData As Variant
Dim ensemblestart, ensembleend, startpoint, datapoint, CurrentLat, CurrentLong, stepcounter, loadend As Integer
Dim scenariotime, ensemble As Variant
Dim PctDone As Single
Dim CarryOver As Variant

On Error GoTo Terminate

' Sets the scenario timeframe
If Worksheets("Convert File").Opt2020.Value = True Then scenariotime = 2020
If Worksheets("Convert File").Opt2050.Value = True Then scenariotime = 2050
If Worksheets("Convert File").Opt2080.Value = True Then scenariotime = 2080
Worksheets("HelpCalc").Range("C13") = scenariotime
Worksheets("Convert File").Range("B10").Value = "Selected scenario: A2 scenario ensemble for the " & scenariotime & "'s"
Worksheets("Convert File").Range("I37").Value = "    A2 scenario for the " & scenariotime & "'s"

' Clears old morphed data
Application.Cursor = xlWait
Worksheets("HadCM3").Range("A1:P9000").Clear
Worksheets("Morphed Weather").Range("A1:BA9000").Clear
Worksheets("CC EPW").Range("A1:A9000").Clear
Worksheets("CC TMY2").Range("A1:A9000").Clear
Worksheets("Convert File").Range("E13:Q21") = 0
Worksheets("Scenario Store").Range("D3:P11") = 0
Worksheets("Convert File").Range("B50").Value = "    No morphed weather file"
Worksheets("Convert File").Range("B51").Value = " "

' Determines the HadCM3 parameter to be loaded
For c = 1 To 9
stepcounter = 0
If c = 1 Then
    d = 4
    FParam = "TEMP"
    ensemblestart = 1
    ensembleend = 3
End If
If c = 2 Then                           '2050's A2b is missing
    d = d + 3
    FParam = "TMAX"
    If scenariotime = 2050 Then
        ensemblestart = 1
        ensembleend = 2
    Else
        ensemblestart = 1
        ensembleend = 3
    End If
End If
If c = 3 Then
    d = d + 3
    FParam = "TMIN"
    If scenariotime = 2050 Then
        ensemblestart = 1
        ensembleend = 2
    Else
        ensemblestart = 1
        ensembleend = 3
    End If
End If
If c = 4 Then
    d = d + 3
    FParam = "DSWF"
    ensemblestart = 1
    ensembleend = 3
End If
If c = 5 Then
    d = d + 3
    FParam = "TCLW"
    ensemblestart = 2
    ensembleend = 3
End If
If c = 6 Then
    d = d + 3
    FParam = "PREC"
    ensemblestart = 1
    ensembleend = 6
End If
If c = 7 Then
    d = d + 3
    FParam = "RHUM"
    ensemblestart = 1
    ensembleend = 3
End If
If c = 8 Then
    d = d + 3
    FParam = "MSLP"
    ensemblestart = 1
    ensembleend = 1
End If
If c = 9 Then
    d = d + 3
    FParam = "WIND"
    ensemblestart = 1
    ensembleend = 6
End If
For x = ensemblestart To ensembleend
    If x = 1 Then
        Sheet_Select = "Data 2a"
        ensemble = "a"
    End If
    If x = 2 Then                           '2050's A2b is missing
        If scenariotime = 2050 And FParam = "TMAX" Or scenariotime = 2050 And FParam = "TMIN" Then
            Sheet_Select = "Data 2c"
            ensemble = "c"
        Else
            Sheet_Select = "Data 2b"
            ensemble = "b"
        End If
    End If
    If x = 3 Then
        Sheet_Select = "Data 2c"
        ensemble = "c"
    End If
    If x = 4 Then                           'For parameters that require mean values (WIND, PREC)
        Sheet_Select = "Data 2a"
        ensemble = "a"
    End If
    If x = 5 Then                           'For parameters that require mean values (WIND, PREC)
        Sheet_Select = "Data 2b"
        ensemble = "b"
    End If
    If x = 6 Then                           'For parameters that require mean values (WIND, PREC)
        Sheet_Select = "Data 2c"
        ensemble = "c"
    End If

If x < 4 Then                               'Loads the change data with respect to the 1980's (dif data not mea data)
'Reads the HadCM3 file header info for conversion
    ThisWorkbook.Worksheets(Sheet_Select).Range("A1:GZ2000").Clear
    FName = Worksheets("Convert File").Data_Path & "\HADCM3_A2" & ensemble & "_" & FParam & "_" & scenariotime & ".dif"
    FNameDisplay = "HADCM3_A2" & ensemble & "_" & FParam & "_" & scenariotime & ".dif"
    Workbooks.OpenText Filename:=FName, DataType:=xlDelimited, Tab:=False, Space:=False, ConsecutiveDelimiter:=False
    startpoint = 0
    datapoint = 0
    For i = 1 To 12
        For k = 1 To 6
            Datareturn = ActiveWorkbook.ActiveSheet.Cells(startpoint + k, 1)
            ThisWorkbook.Worksheets(Sheet_Select).Cells(datapoint + k, 1) = Datareturn
        Next k
        If c = 9 Then                           ' for WIND data
            startpoint = startpoint + 698
            datapoint = datapoint + 78
        Else                                    ' for all other data
            startpoint = startpoint + 707
            datapoint = datapoint + 79
        End If
    Next i
    ActiveWorkbook.Close SaveChanges:=False
    Application.ScreenUpdating = True

'Reads the HadCM3 data for conversion
    Workbooks.OpenText Filename:=FName, DataType:=xlDelimited, Tab:=True, Space:=True, ConsecutiveDelimiter:=True
    ActiveWorkbook.ActiveSheet.Columns(1).Delete
    k = 7
    M = 7
    n = 1
    o = 1
    For g = 1 To 12
        ' Displays the progress bar
        PctDone = g / 12
        With Status
            .CurrentFile.Text = FNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Gets the montly data in the correct format
        If c = 9 Then                       ' for WIND data
            loadend = 6912
        Else                                ' for all other data
            loadend = 7008
        End If
        For i = 1 To loadend
            Datareturn = ActiveWorkbook.ActiveSheet.Cells(k, n)
            ThisWorkbook.Worksheets(Sheet_Select).Cells(M, o) = Datareturn
            n = n + 1
            o = o + 1
            If n / 11 = 1 Then
                k = k + 1
                n = 1
            End If
            If i / 96 = Int(i / 96) Then
                M = M + 1
                o = 1
            End If
        Next i
        k = k + 7
        M = M + 6
        n = 1
        o = 1
    Next g
    ActiveWorkbook.Close SaveChanges:=False

' Returns the grid point data to the HadCM3 sheet
    ' Inserts the header information
    With Worksheets("HadCM3")
        .Range("A1") = Worksheets("HelpCalc").Range("C9")
        .Range("A1").Font.Bold = True
        .Range("A1").Font.Size = 12
        .Range("A2") = Worksheets("Convert File").Range("B10")
        ' Inserts individual file data for the four datapoints and formats the table
        For y = 1 To 4
            If c = 9 Then                                               ' for WIND data
                CurrentLat = Worksheets("HelpCalc").Cells(y + 1, 16)
                CurrentLong = Worksheets("HelpCalc").Cells(y + 1, 17)
            Else                                                        ' for all other data
                CurrentLat = Worksheets("HelpCalc").Cells(y + 1, 7)
                CurrentLong = Worksheets("HelpCalc").Cells(y + 1, 8)
            End If
            .Cells(d, 1).Font.Bold = True
            .Cells(d, 1) = "File:"
            .Cells(d, 3).Font.Bold = True
            .Cells(d, 3) = FNameDisplay
            .Cells(d, 8) = "Data format: " & Worksheets(Sheet_Select).Cells(5, 1)
            .Cells(d + 2, 1) = "Grid Lat/Long"
            .Cells(d + 2, 3) = "Jan"
            .Cells(d + 2, 4) = "Feb"
            .Cells(d + 2, 5) = "Mar"
            .Cells(d + 2, 6) = "Apr"
            .Cells(d + 2, 7) = "May"
            .Cells(d + 2, 8) = "Jun"
            .Cells(d + 2, 9) = "Jul"
            .Cells(d + 2, 10) = "Aug"
            .Cells(d + 2, 11) = "Sep"
            .Cells(d + 2, 12) = "Oct"
            .Cells(d + 2, 13) = "Nov"
            .Cells(d + 2, 14) = "Dec"
            .Cells(d + 2 + y, 1) = CurrentLat & "/" & CurrentLong
            .Cells(d + 7, 1) = "Average"
            CurrentLat = CurrentLat + 6
            For t = 1 To 12
                .Cells(d + 2 + y, t + 2) = Worksheets(Sheet_Select).Cells(CurrentLat, CurrentLong)
                If c = 9 Then
                    CurrentLat = CurrentLat + 78
                Else
                    CurrentLat = CurrentLat + 79
                End If
            Next t
            ' Transforms MSLP from Pa to hPa
            If c = 8 Then
                For P = 1 To 12
                    .Cells(d + 2 + y, P + 2) = .Cells(d + 2 + y, P + 2) / 100
                Next P
            End If
        Next y
        ' Calculates the average values for the four datapoints and formats the table
        For Z = 1 To 12
            ChangePoint = Worksheets("HelpCalc").Cells(1, Z + 2).AddressLocal(RowAbsolute:=False, ColumnAbsolute:=False)
            ChangePoint = Left(ChangePoint, 1)
            AvData = ChangePoint & d + 3 & ":" & ChangePoint & d + 6
            .Cells(d + 7, Z + 2) = Application.WorksheetFunction.Average(Worksheets("HadCM3").Range(AvData))
            .Cells(d + 2, Z + 2).HorizontalAlignment = xlRight
        Next Z
        AvData = "A" & d + 2 & ":N" & d + 2
        .Range(AvData).Borders(xlEdgeBottom).LineStyle = xlContinuous
        .Range(AvData).Font.Bold = True
        AvData = "A" & d + 7 & ":N" & d + 7
        .Range(AvData).Borders(xlEdgeTop).LineStyle = xlContinuous
        .Range(AvData).Font.Bold = True
        d = d + 10
    End With
    stepcounter = stepcounter + 1

Else    'Loads baseline 1980's data for WIND and PREC (mea data not dif data)

'Reads the HadCM3 file header info for conversion
    FName = Worksheets("Convert File").Data_Path & "\HADCM3_A2" & ensemble & "_" & FParam & "_" & "1980.mea"
    FNameDisplay = "HADCM3_A2" & ensemble & "_" & FParam & "_" & "1980.mea"
    Workbooks.OpenText Filename:=FName, DataType:=xlDelimited, Tab:=False, Space:=False, ConsecutiveDelimiter:=False
    startpoint = 0
    If c = 9 Then
        datapoint = 939                         ' for WIND data
    Else
        datapoint = 951                         ' for PREC data
    End If
    For i = 1 To 12
        For k = 1 To 5
            Datareturn = ActiveWorkbook.ActiveSheet.Cells(startpoint + k, 1)
            ThisWorkbook.Worksheets(Sheet_Select).Cells(datapoint + k, 1) = Datareturn
        Next k
        If c = 9 Then                           ' for WIND data
            startpoint = startpoint + 697
            datapoint = datapoint + 77
        Else                                    ' for PREC data
            startpoint = startpoint + 706
            datapoint = datapoint + 78
        End If
    Next i
    ActiveWorkbook.Close SaveChanges:=False
    Application.ScreenUpdating = True

'Reads the HadCM3 data for conversion
    Workbooks.OpenText Filename:=FName, DataType:=xlDelimited, Tab:=True, Space:=True, ConsecutiveDelimiter:=True
    ActiveWorkbook.ActiveSheet.Columns(1).Delete
    k = 6
    If c = 9 Then
        M = 945                               ' for WIND data
    Else
        M = 957                               ' for PREC data
    End If
    n = 1
    o = 1
    For g = 1 To 12
        ' Displays the progress bar
        PctDone = g / 12
        With Status
            .CurrentFile.Text = FNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Gets the montly data in the correct format
        If c = 9 Then                       ' for WIND data
            loadend = 6912
        Else                                ' for PREC data
            loadend = 7008
        End If
        For i = 1 To loadend
            Datareturn = ActiveWorkbook.ActiveSheet.Cells(k, n)
            ThisWorkbook.Worksheets(Sheet_Select).Cells(M, o) = Datareturn
            n = n + 1
            o = o + 1
            If n / 11 = 1 Then
                k = k + 1
                n = 1
            End If
            If i / 96 = Int(i / 96) Then
                M = M + 1
                o = 1
            End If
        Next i
        k = k + 6
        M = M + 5
        n = 1
        o = 1
    Next g
    ActiveWorkbook.Close SaveChanges:=False

' Returns the mean value adapted grid point data to the HadCM3 sheet (WIND, PREC)
        ' Updates the individual file data for the four datapoints
        If x = 4 Then d = d - 30
        Worksheets("HadCM3").Cells(d, 3) = Worksheets("HadCM3").Cells(d, 3) & " & 1980.mea"
        If c = 6 Then Worksheets("HadCM3").Cells(d, 8) = "Data format: Total Precipitation Change (%)"
        If c = 9 Then Worksheets("HadCM3").Cells(d, 8) = "Data format: Wind Speed Change (%)"
        For y = 1 To 4
            If c = 9 Then                                               ' for WIND data
                CurrentLat = Worksheets("HelpCalc").Cells(y + 1, 16)
                CurrentLong = Worksheets("HelpCalc").Cells(y + 1, 17)
            Else                                                        ' for all other data - PREC
                CurrentLat = Worksheets("HelpCalc").Cells(y + 1, 7)
                CurrentLong = Worksheets("HelpCalc").Cells(y + 1, 8)

            End If
            If c = 6 Then CurrentLat = CurrentLat + 956
            If c = 9 Then CurrentLat = CurrentLat + 944
            For t = 1 To 12
                If Worksheets(Sheet_Select).Cells(CurrentLat, CurrentLong) = 0 Then
                    Worksheets("HadCM3").Cells(d + 2 + y, t + 2) = 0
                Else
                    Worksheets("HadCM3").Cells(d + 2 + y, t + 2) = Worksheets("HadCM3").Cells(d + 2 + y, t + 2) / Worksheets(Sheet_Select).Cells(CurrentLat, CurrentLong) * 100
                    If c = 9 Then
                        CurrentLat = CurrentLat + 77
                    Else
                        CurrentLat = CurrentLat + 78
                    End If
                End If
            Next t
        Next y
        ' Calculates the average values for the four datapoints and returns them to the table
        For Z = 1 To 12
            ChangePoint = Worksheets("HelpCalc").Cells(1, Z + 2).AddressLocal(RowAbsolute:=False, ColumnAbsolute:=False)
            ChangePoint = Left(ChangePoint, 1)
            AvData = ChangePoint & d + 3 & ":" & ChangePoint & d + 6
            Worksheets("HadCM3").Cells(d + 7, Z + 2) = Application.WorksheetFunction.Average(Worksheets("HadCM3").Range(AvData))
        Next Z
        d = d + 10
End If

Next x


' Calculates the average values for the three model runs and formats the table
With Worksheets("HadCM3")
    AvData = "A" & d & ":N" & d
    .Range(AvData).Font.Bold = True
    .Range(AvData).Borders(xlEdgeBottom).LineStyle = xlDouble
    .Cells(d, 1) = FParam & " ensemble:"
    ' Returns the data to the main worksheet
    For b = 1 To 12
        If stepcounter = 1 Then .Cells(d, b + 2) = .Cells(d - 3, b + 2)
        If stepcounter = 2 Then .Cells(d, b + 2) = (.Cells(d - 3, b + 2) + .Cells(d - 13, b + 2)) / 2
        If stepcounter = 3 Then .Cells(d, b + 2) = (.Cells(d - 3, b + 2) + .Cells(d - 13, b + 2) + .Cells(d - 23, b + 2)) / 3
        If c = 5 Then
            Worksheets("Convert File").Cells(12 + c, b + 4) = .Cells(d, b + 2) * 100
        Else
            Worksheets("Convert File").Cells(12 + c, b + 4) = .Cells(d, b + 2)
        End If
    Next b
    AvData = "E" & 12 + c & ":P" & 12 + c
    Worksheets("Convert File").Cells(12 + c, 17) = Application.WorksheetFunction.Average(Worksheets("Convert File").Range(AvData))
End With

Next c

' Copies the scenario information to a safe sheet to avoid mistakes
CarryOver = Worksheets("Convert File").Range("E13:Q21")
Worksheets("Scenario Store").Range("D3:P11") = CarryOver

' Unloads the progress bar
Application.ScreenUpdating = True
Application.Cursor = xlDefault
Unload Status

End


' On error loading the HadCM3 data the program is terminated and any data on the main data sheet deleted
Terminate:
Application.Cursor = xlDefault
Application.ScreenUpdating = True
Unload Status
Worksheets("Convert File").Range("B10").Value = "No scenario selected"
Worksheets("Convert File").Activate
MsgBox "An error has occured while loading the file HADCM3_A2" & ensemble & "_" & FParam & "_" & scenariotime & ".dif! Please check whether the file has been downloaded and is located in the correct directory! Please also note that some files require renaming after downloading as detailed in the manual." & Chr(10) & Chr(10) & "If the problem persists please try to close and reopen the climate change weather file generator.", vbExclamation
Exit Sub

End Sub

Sub ProgressBarMorphing()

Dim scenariotime, Celldaypoint, DayYear As Integer
Dim Scalingfactor As Double
Dim Datareturn, Datacheck As Variant
Dim FormText, StringCopy, ParaNameDisplay As String
Dim Ordernumber, StringPos, CharacterNo, StringReturn As Integer
Dim PctDone As Single

On Error GoTo Terminate

' Resets the data on the main sheet in case it was changed
 Datacheck = Worksheets("Selection Data").Range("R33:S35")
 Worksheets("Convert File").Range("F35:G37") = Datacheck

' Clears previous data in the worksheet and jumps to the morphed weather sheet, sets the cursor hourglass
    Worksheets("Morphed Weather").Range("A1:BA9000").Clear
    Worksheets("CC EPW").Range("A1:A9000").Clear
    Worksheets("CC TMY2").Range("A1:A9000").Clear
    Worksheets("Convert File").Range("B50").Value = "    No morphed weather file"
    Worksheets("Convert File").Range("B51").Value = " "
    Application.Goto reference:=Worksheets("Morphed Weather").Range("A1")
    Application.Cursor = xlWait

' Copies and returns the header information
    With Worksheets("Morphed Weather")
        .Range("A1") = Worksheets("EPW").Range("A1")
        .Range("A2") = "DESIGN CONDITIONS,0"
        .Range("A3") = "TYPICAL/EXTREME PERIODS,0"
        .Range("A5") = Worksheets("EPW").Range("A5")
        .Range("A6").Value = "COMMENTS 1, This climate change adapted weather file, which bases on HadCM3 GCM A2 ensemble data, has been generated using the CCWorldWeatherGen tool V.1.9. The original weather file used for generating this climate change adapted weather data may be copyrighted material. Therefore, generated weather files can only be used by persons or entities who possess the corresponding licensed weather file. DISCLAIMER OF WARRANTIES: The entire risk as to the quality and performance of the calculated climate change weather data in this file is with you. In no event will the authors of the weather file generation tool be liable to you for any damages, including without limitation any lost profits, lost savings, or other incidental or consequential damages arising out of the use or inability to use this data."
        .Range("A8").Value = "DATA PERIODS,1,1,Data,Sunday,1/1,12/31"
        .Range("A10:AF10").Value = Worksheets("EPW").Range("A10:AF10").Value
        .Range("F11:F8770").HorizontalAlignment = xlRight
        .Range("A1").Font.Bold = True
        .Range("A1").Font.Size = 12
        .Range("A10:AG10").Font.Bold = True
        .Range("A10:AG10").HorizontalAlignment = xlRight
        .Range("A10:AG10").Borders(xlEdgeBottom).LineStyle = xlContinuous
        .Range("AG10") = "DayYear"

    ' Corrects incorrect time zones
        If .Range("A1").Value = "LOCATION,Tenerife,-,ESP,SWEC,600200,28.47,-16.25,1.0,46.0" Then
            .Range("A1").Value = "LOCATION,Tenerife,-,ESP,SWEC,600200,28.47,-16.25,0.0,46.0"
            Worksheets("Selection Data").Range("R36") = 0
        End If
        If .Range("A1").Value = "LOCATION,Las Palmas,-,ESP,SWEC,600300,27.93,-15.38,1.0,25.0" Then
            .Range("A1").Value = "LOCATION,Las Palmas,-,ESP,SWEC,600300,27.93,-15.38,0.0,25.0"
            Worksheets("Selection Data").Range("R36") = 0
        End If

    ' Generates the second comment line from the original line
        FormText = ThisWorkbook.Worksheets("EPW").Range("A6")
        StringPos = InStr(1, FormText, ",")
        CharacterNo = ThisWorkbook.Worksheets("EPW").Range("A6").Characters.Count
        StringReturn = CharacterNo - StringPos
        StringCopy = Right(FormText, StringReturn)
        .Range("A7") = "COMMENTS 2, The original data use for generating this file is: " & StringCopy

    End With

' Inserts the basic weather file data
    If Worksheets("Convert File").Opt2020.Value = True Then scenariotime = 2020
    If Worksheets("Convert File").Opt2050.Value = True Then scenariotime = 2050
    If Worksheets("Convert File").Opt2080.Value = True Then scenariotime = 2080
    With Worksheets("Morphed Weather")
        .Range("A11:A8770").Value = scenariotime
        .Range("B11:E8770").Value = Worksheets("EPW").Range("B11:E8770").Value
        .Range("F11:F8770").Value = "*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?*?"
        ' Wind direction: Remains unchanged
        .Range("U11:U8770").Value = Worksheets("EPW").Range("U11:U8770").Value
    End With

' Generates the day of the year
    ParaNameDisplay = "Day of the Year"
    Celldaypoint = 11
    DayYear = 1
    For i = 1 To 365
        ' Displays the progress bar
        PctDone = i / 365
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        'Performs the day of the year calculations
        For k = Celldaypoint To Celldaypoint + 23
            ThisWorkbook.Worksheets("Morphed Weather").Cells(k, 33) = DayYear
        Next k
        Celldaypoint = Celldaypoint + 24
        DayYear = DayYear + 1
    Next i

' Morphing procedures

' Dry bulb temperature: Derived of UKCIP TEMP, TMAX and TMIN (DryT / TEMP)
    Dim DryT_TEMP As Double
    Dim Monthjump As Integer
    Dim Cellshift As Integer
    Dim SmoothTMAX, SmoothTMIN, SmoothAV, SmoothMAX, SmoothMIN As Double
    Dim TMAXStep, TMINStep, AVStep, MAXStep, MINStep As Double
    Dim Temprange_Copy As Variant

    ' Generates the mean / maximum / minimum temps in the Selection Data sheet
    Temprange_Copy = Worksheets("EPW").Range("G11:G8770")
    Worksheets("Selection Data").Range("D12:D8771").Value = Temprange_Copy

    ' Starts the morphing procedure for dry bulb temperature
    ParaNameDisplay = "Dry Bulb Temperature (°C)"
    P = 0
    For j = 1 To 12
        If j = 1 Then
            Cellshift = 10
            Monthjump = 768                                     '744 hours + 24 hours for smoothing
        End If
        If j = 2 Then
            Cellshift = Cellshift + Monthjump
            Monthjump = 672
        End If
        If j = 3 Or j = 5 Or j = 7 Or j = 8 Or j = 10 Then
            Cellshift = Cellshift + Monthjump
            Monthjump = 744
        End If
        If j = 4 Or j = 6 Or j = 9 Or j = 11 Or j = 12 Then     '720 hours for December - 24 hours for smoothing
            Cellshift = Cellshift + Monthjump
            Monthjump = 720
        End If
        For i = Cellshift + 1 To Cellshift + Monthjump
            ' Displays the progress bar
            PctDone = P + i / 8760
            With Status_M
                .CurrentFile.Text = ParaNameDisplay
                .FrameProgress.Caption = Format(PctDone, "0%")
                .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
            End With
            DoEvents
            ' Performs the morphing
            Scalingfactor = (Worksheets("Convert File").Cells(14, j + 4) - Worksheets("Convert File").Cells(15, j + 4)) / (Worksheets("Selection Data").Cells(7, j + 1) - Worksheets("Selection Data").Cells(9, j + 1))
            DryT_TEMP = Worksheets("EPW").Cells(i, 7) + Worksheets("Convert File").Cells(13, j + 4) + (Scalingfactor * (Worksheets("EPW").Cells(i, 7) - Worksheets("Selection Data").Cells(5, j + 1)))
            If DryT_TEMP > 50 Then DryT_TEMP = 50
            If DryT_TEMP < -50 Then DryT_TEMP = -50
            Worksheets("Morphed Weather").Cells(i, 7) = DryT_TEMP
            If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                Smooth = (Worksheets("Convert File").Cells(13, j + 5) - Worksheets("Convert File").Cells(13, j + 4)) / 48
                SmoothTMAX = (Worksheets("Convert File").Cells(14, j + 5) - Worksheets("Convert File").Cells(14, j + 4)) / 48
                SmoothTMIN = (Worksheets("Convert File").Cells(15, j + 5) - Worksheets("Convert File").Cells(15, j + 4)) / 48
                SmoothAV = (Worksheets("Selection Data").Cells(5, j + 2) - Worksheets("Selection Data").Cells(5, j + 1)) / 48
                SmoothMAX = (Worksheets("Selection Data").Cells(7, j + 2) - Worksheets("Selection Data").Cells(7, j + 1)) / 48
                SmoothMIN = (Worksheets("Selection Data").Cells(9, j + 2) - Worksheets("Selection Data").Cells(9, j + 1)) / 48
                For k = 1 To 48
                    Smoothstep = Smooth * k
                    TMAXStep = SmoothTMAX * k
                    TMINStep = SmoothTMIN * k
                    AVStep = SmoothAV * k
                    MAXStep = SmoothMAX * k
                    MINStep = SmoothMIN * k
                    Scalingfactor = ((Worksheets("Convert File").Cells(14, j + 4) + TMAXStep) - (Worksheets("Convert File").Cells(15, j + 4) + TMINStep)) / ((Worksheets("Selection Data").Cells(7, j + 1) + MAXStep) - (Worksheets("Selection Data").Cells(9, j + 1) + MINStep))
                    DryT_TEMP = Worksheets("EPW").Cells(i - 24 + k, 7) + (Worksheets("Convert File").Cells(13, j + 4) + Smoothstep) + (Scalingfactor * (Worksheets("EPW").Cells(i - 24 + k, 7) - (Worksheets("Selection Data").Cells(5, j + 1) + AVStep)))
                    If DryT_TEMP > 50 Then DryT_TEMP = 50
                    If DryT_TEMP < -50 Then DryT_TEMP = -50
                    Worksheets("Morphed Weather").Cells(i - 24 + k, 7) = DryT_TEMP
                Next k
                i = i + 24
            End If
        Next i
        P = 0
    Next j


' Total Sky Cover & Opaque Sky Cover: Simple scaling (TotSky & OpSky / TCLW)
    Dim TotSky_TCLW, OpSky_TCLW As Double
    ParaNameDisplay = "Total & Opaque Sky Cover (1/10)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Performs the morphing
        For j = 1 To 12
            If Worksheets("Morphed Weather").Cells(i, 2).Value = j Then
                Scalingfactor = Worksheets("Convert File").Cells(17, j + 4) / 10
                TotSky_TCLW = Worksheets("EPW").Cells(i, 23) + Scalingfactor
                If Worksheets("EPW").Cells(i, 23) = 0 Then
                    OpSky_TCLW = TotSky_TCLW / 2
                Else
                    OpSky_TCLW = (TotSky_TCLW * Worksheets("EPW").Cells(i, 24)) / Worksheets("EPW").Cells(i, 23)
                End If
                If OpSky_TCLW > TotSky_TCLW Then OpSky_TCLW = TotSky_TCLW
                If TotSky_TCLW > 10 Then TotSky_TCLW = 10
                If TotSky_TCLW < 0 Then TotSky_TCLW = 0
                If OpSky_TCLW > 10 Then OpSky_TCLW = 10
                If OpSky_TCLW < 0 Then OpSky_TCLW = 0
                Worksheets("Morphed Weather").Cells(i, 23) = TotSky_TCLW
                Worksheets("Morphed Weather").Cells(i, 24) = OpSky_TCLW
                If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                    Smooth = (Worksheets("Convert File").Cells(17, j + 5) - Worksheets("Convert File").Cells(17, j + 4)) / 48
                    For k = 1 To 48
                        Smoothstep = Smooth * k
                        Scalingfactor = (Worksheets("Convert File").Cells(17, j + 4) + Smoothstep) / 10
                        TotSky_TCLW = Worksheets("EPW").Cells(i - 24 + k, 23) + Scalingfactor
                        If Worksheets("EPW").Cells(i - 24 + k, 23) = 0 Then
                            OpSky_TCLW = TotSky_TCLW / 2
                        Else
                            OpSky_TCLW = (TotSky_TCLW * Worksheets("EPW").Cells(i - 24 + k, 24)) / Worksheets("EPW").Cells(i - 24 + k, 23)
                        End If
                        If TotSky_TCLW > 10 Then TotSky_TCLW = 10
                        If TotSky_TCLW < 0 Then TotSky_TCLW = 0
                        If OpSky_TCLW > 10 Then OpSky_TCLW = 10
                        If OpSky_TCLW < 0 Then OpSky_TCLW = 0
                        Worksheets("Morphed Weather").Cells(i - 24 + k, 23) = TotSky_TCLW
                        Worksheets("Morphed Weather").Cells(i - 24 + k, 24) = OpSky_TCLW
                    Next k
                i = i + 24
                End If
            End If
        Next j
    Next i

' Relative humidity: Simple shifting (RelHum - RHUM)
    Dim RelHum_RHUM As Double
    ParaNameDisplay = "Relative Humidity (%)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Performs the morphing
        For j = 1 To 12
            If Worksheets("Morphed Weather").Cells(i, 2).Value = j Then
                RelHum_RHUM = Worksheets("EPW").Cells(i, 9) + Worksheets("Convert File").Cells(19, j + 4)
                If RelHum_RHUM > 100 Then RelHum_RHUM = 100
                If RelHum_RHUM < 1 Then RelHum_RHUM = 1
                Worksheets("Morphed Weather").Cells(i, 9) = RelHum_RHUM
                If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                    Smooth = (Worksheets("Convert File").Cells(19, j + 5) - Worksheets("Convert File").Cells(19, j + 4)) / 48
                    For k = 1 To 48
                        Smoothstep = Smooth * k
                        RelHum_RHUM = Worksheets("EPW").Cells(i - 24 + k, 9) + Worksheets("Convert File").Cells(19, j + 4) + Smoothstep
                        If RelHum_RHUM > 100 Then RelHum_RHUM = 100
                        If RelHum_RHUM < 1 Then RelHum_RHUM = 1
                        Worksheets("Morphed Weather").Cells(i - 24 + k, 9) = RelHum_RHUM
                    Next k
                    i = i + 24
                End If
            End If
        Next j
    Next i

' Dew point temperature (according to ASHRAE handbook) &
' Horizontal infrared radiation from the sky for EPW file generation following Crawford and Duchon
    Dim DryT_Round As Integer
    Dim Pws_Value, Pws_Step, Pws_StepPoint As Double
    Dim WaterSaturationPressure, WaterVapPressure As Double
    Dim DewPoiT, alpha As Double
    Dim CloudCover, DryTKelvin, AtmoEmis, IRRad As Double
    ParaNameDisplay = "Dew Point Temp (°C) & Infrared Radiation (Wh/m²)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Performs the morphing
        ' Determines the partial pressure Pw of the water vapour (based on RH)
        ' a) Determine the temperature for looking up the correct saturation pressure, interpolate, return value
        DryT_Round = Round(ThisWorkbook.Worksheets("Morphed Weather").Cells(i, 7))
        If ThisWorkbook.Worksheets("Morphed Weather").Cells(i, 7) - DryT_Round < 0 Then
            DryT_Round = DryT_Round - 1
        End If
        Pws_Value = ThisWorkbook.Worksheets("Psychrometrics").Cells(DryT_Round + 64, 2)
        Pws_Step = (ThisWorkbook.Worksheets("Psychrometrics").Cells(DryT_Round + 65, 2) - ThisWorkbook.Worksheets("Psychrometrics").Cells(DryT_Round + 64, 2)) / 10
        Pws_StepPoint = (Round(ThisWorkbook.Worksheets("Morphed Weather").Cells(i, 7), 1) - DryT_Round) * 10
        ' b) Calculate saturation pressure
        WaterSaturationPressure = Pws_Value + (Pws_Step * Pws_StepPoint)
        ' c) Calculate vapour pressure
        WaterVapPressure = ThisWorkbook.Worksheets("Morphed Weather").Cells(i, 9) / 100 * WaterSaturationPressure
        ' d) Calculates the dew point temperature
        alpha = Log(WaterVapPressure)
        If WaterVapPressure >= 0.61115 Then
            DewPoiT = 6.54 + (14.526 * alpha) + (0.7389 * alpha ^ 2) + (0.09486 * alpha ^ 3) + (0.4569 * WaterVapPressure ^ 0.1984)
        Else
            DewPoiT = 6.09 + (12.608 * alpha) + (0.4959 * alpha ^ 2)
        End If
        If DewPoiT > 30 Then DewPoiT = 30
        If DewPoiT < -60 Then DewPoiT = -60
        ThisWorkbook.Worksheets("Morphed Weather").Cells(i, 8) = DewPoiT
        ' e) Calculates the IR radiation
        CloudCover = Worksheets("Morphed Weather").Cells(i, 23) / 10
        DryTKelvin = Worksheets("Morphed Weather").Cells(i, 7) + 273.15
        AtmoEmis = CloudCover + ((1 - CloudCover) * 1.24 * ((WaterVapPressure * 10) / DryTKelvin) ^ (1 / 7))
        IRRad = AtmoEmis * 0.0000000567 * DryTKelvin ^ 4
        Worksheets("Morphed Weather").Cells(i, 13) = IRRad
    Next i

' Atmospheric pressure: Simple monthly shifting (Press / MSLP)
    Dim Press_MSLP As Double
    ParaNameDisplay = "Atmospheric Pressure (Pa)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Performs the morphing
        For j = 1 To 12
            If Worksheets("Morphed Weather").Cells(i, 2).Value = j Then
                Press_MSLP = Worksheets("EPW").Cells(i, 10) + (Worksheets("Convert File").Cells(20, j + 4) * 100)
                If Press_MSLP > 110000 Then Press_MSLP = 110000
                If Press_MSLP < 70000 Then Press_MSLP = 70000
                Worksheets("Morphed Weather").Cells(i, 10) = Press_MSLP
                If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                    Smooth = (Worksheets("Convert File").Cells(20, j + 5) - Worksheets("Convert File").Cells(20, j + 4)) * 100 / 48
                    For k = 1 To 48
                        Smoothstep = Smooth * k
                        Press_MSLP = Worksheets("EPW").Cells(i - 24 + k, 10) + (Worksheets("Convert File").Cells(20, j + 4) * 100) + Smoothstep
                        If Press_MSLP > 110000 Then Press_MSLP = 110000
                        If Press_MSLP < 70000 Then Press_MSLP = 70000
                        Worksheets("Morphed Weather").Cells(i - 24 + k, 10) = Press_MSLP
                    Next k
                i = i + 24
                End If
            End If
        Next j
    Next i

' Calculates the extraterrestrial radiation parameters
    For i = 11 To 8770
    ParaNameDisplay = "Extraterrestrial Radiation Parameters (Wh/m²)"
    ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
    ' Performs the morphing
    ' Calculates the solar baseline parameters
        Dim LongitudeCalc, LocalTime, TimeZone As Double
        Dim VarE, VarB As Double
        Dim HourAngle, Declination, sin_SolarAlt, SolarAlt, SolarAltR, PII As Double
        Dim DayofYear As Double
        Dim ExtHorRad_spot, ExtHorRad_sum As Double
        Dim nosunminutes As Integer
        PII = Application.WorksheetFunction.Pi
        LongitudeCalc = Worksheets("Selection Data").Range("R34")
        TimeZone = Worksheets("Selection Data").Range("R36")
        DayofYear = Worksheets("Morphed Weather").Cells(i, 33)
        VarB = DayofYear * 360 / 365.25
        VarE = -0.128 * Sin((VarB - 2.8) * PII / 180) - 0.165 * Sin((2 * VarB + 19.7) * PII / 180)
        ExtHorRad_sum = 0
    ' Calculates and checks whether it is night or day for Extraterrestrial Direct Normal Radiation
        Dim CorrectionFac As Double
        CorrectionFac = 1 + 0.03344 * Cos((VarB - 2.8) * PII / 180)
        ExtDirRad = 1367 * CorrectionFac
        nosunminutes = 0
    ' Calculates Extraterrestrial Horizontal Radiation
        For k = 1 To 60
            LocalTime = (Worksheets("Morphed Weather").Cells(i, 4) - 1 + k / 60) + (LongitudeCalc - (TimeZone * 15)) / 15 + VarE
            HourAngle = 15 * (LocalTime - 12)
            Declination = Application.WorksheetFunction.Asin(0.3978 * Sin((VarB * PII / 180) - 1.4 + 0.0355 * Sin((VarB * PII / 180) - 0.0489)))
            sin_SolarAlt = Sin(Declination) * Sin(Worksheets("Selection Data").Range("R33") * PII / 180) + _
                Cos(Declination) * Cos(Worksheets("Selection Data").Range("R33") * PII / 180) * _
                Cos(HourAngle * PII / 180)
            SolarAlt = Application.WorksheetFunction.Asin(sin_SolarAlt) * 180 / PII
            ExtHorRad_spot = Sin(SolarAlt * PII / 180) * ExtDirRad
            If ExtHorRad_spot < 0 Then
                ExtHorRad_spot = 0
                nosunminutes = nosunminutes + 1
            End If
            ExtHorRad_sum = ExtHorRad_sum + ExtHorRad_spot
        Next k
        ExtHorRad = ExtHorRad_sum / 60
        ExtDirRad = ExtDirRad * ((60 - nosunminutes) / 60)
    ' Writes back Extraterrestrial Horizontal Radiation & Extraterrestrial Direct Normal Radiation
        If ExtHorRad > 1415 Then ExtHorRad = 1415
        Worksheets("Morphed Weather").Cells(i, 11) = ExtHorRad
        If ExtDirRad > 1415 Then ExtDirRad = 1415
        Worksheets("Morphed Weather").Cells(i, 12) = ExtDirRad
    Next i

' Checks the solar radiation data integety
    Dim Recheckvalue, kstep, solarhour, solarhour_calc, solarhourcheck, solarhourfinal, solarhourcheck_calc, hourshiftup, hourshiftdown, overridecheck, TimeShiftStep, RadCalcFail As Integer
    Dim TimeShift, TimeShiftAdj, AdjustRad, DaySunrise, DaySunset As Double
    Worksheets("Selection Data").Range("Z11:AD8770").Clear
    Worksheets("EPW").Range("N11:N8770").Copy Destination:=Worksheets("Selection Data").Range("AC11:AC8770")
    Recheckvalue = 0
    TimeShiftStep = 0
    solarhourfinal = 0
    With Worksheets("Radiation Log")
            .Range("A1:E9000").Clear
            .Range("A1").Font.Bold = True
            .Range("A1").Font.Size = 12
            .Range("A1") = "Solar Radiation Inconsistencies Log File, CCWorldweatherGen, V.1.9 for Site: " & Worksheets("HelpCalc").Range("C11") & ", " & Worksheets("HelpCalc").Range("C12")
            .Range("A2") = "This file lists all hours where the global horizontal raditation of the original EPW file is equal to or exceeds the calculated extraterrestrial horizontal radiation. Where this frequently occurs at a significant level,"
            .Range("A3") = "i.e. the extraterrestrial radiation is 0 and the global horiziontal radiation exceeds 10 Wh/m2 this may indicate that the EPW file follows the solar time convention instead of the commonly used local standard time convention."
            .Range("A5:E5").Font.Bold = True
            .Range("A5:E5").HorizontalAlignment = xlRight
            .Range("A5:E5").Borders(xlEdgeBottom).LineStyle = xlContinuous
            .Range("A5") = "Month"
            .Range("B5") = "Day"
            .Range("C5") = "Hour"
            .Range("D5") = "ExHoRad(calc)"
            .Range("E5") = "GlHoRad(EPW)"
    End With
Recheck:
    solarhour = 0
    solarhourcheck = 0
    hourshiftdown = 0
    hourshiftup = 0
    overridecheck = 0
    RadCalcFail = 0
    If Recheckvalue = 0 Then
        ParaNameDisplay = "Solar Radiation Parameter Consistency Check"
    Else
        ParaNameDisplay = "Re-examining Solar Radiation Consistency"
        Worksheets("Radiation Log").Cells(solarhourfinal + 2, 1) = "The following list provides the remaining solar radiation inconsistencies after re-adjustment."
    End If
    For i = 11 To 8770
    ' Displays the progress bar
        PctDone = i / 17520
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
    ' Checks whether Extraterrestrial Horizontal Radiation is larger than EPW Global Horizontal Radiation
        If Worksheets("Morphed Weather").Cells(i, 11) < Worksheets("EPW").Cells(i, 14) And Worksheets("EPW").Cells(i - 1, 14) = 0 Then solarhour = solarhour + 1
        If Worksheets("Morphed Weather").Cells(i, 11) < Worksheets("EPW").Cells(i, 14) And Worksheets("EPW").Cells(i + 1, 14) = 0 Then solarhour = solarhour + 1
        If Round(Worksheets("Morphed Weather").Cells(i, 11), 0) <= Worksheets("EPW").Cells(i, 14) And Worksheets("EPW").Cells(i, 14) > 0 And Recheckvalue = 0 Then
            solarhourcheck = solarhourcheck + 1
            Worksheets("Radiation Log").Cells(solarhourcheck + 5, 1) = Worksheets("EPW").Cells(i, 2)
            Worksheets("Radiation Log").Cells(solarhourcheck + 5, 2) = Worksheets("EPW").Cells(i, 3)
            Worksheets("Radiation Log").Cells(solarhourcheck + 5, 3) = Worksheets("EPW").Cells(i, 4)
            Worksheets("Radiation Log").Cells(solarhourcheck + 5, 4) = Round(Worksheets("Morphed Weather").Cells(i, 11), 0)
            Worksheets("Radiation Log").Cells(solarhourcheck + 5, 5) = Worksheets("EPW").Cells(i, 14)
            If Worksheets("EPW").Cells(i, 14) - Round(Worksheets("Morphed Weather").Cells(i, 11), 0) > 15 Then overridecheck = overridecheck + 1
        End If
        If Round(Worksheets("Morphed Weather").Cells(i, 11), 0) <= Round(Worksheets("EPW").Cells(i, 14), 0) And Round(Worksheets("EPW").Cells(i, 14), 0) > 0 And Recheckvalue = 1 Then
            solarhourcheck = solarhourcheck + 1
            Worksheets("Radiation Log").Cells(solarhourfinal + solarhourcheck + 2, 1) = Worksheets("EPW").Cells(i, 2)
            Worksheets("Radiation Log").Cells(solarhourfinal + solarhourcheck + 2, 2) = Worksheets("EPW").Cells(i, 3)
            Worksheets("Radiation Log").Cells(solarhourfinal + solarhourcheck + 2, 3) = Worksheets("EPW").Cells(i, 4)
            Worksheets("Radiation Log").Cells(solarhourfinal + solarhourcheck + 2, 4) = Round(Worksheets("Morphed Weather").Cells(i, 11), 0)
            Worksheets("Radiation Log").Cells(solarhourfinal + solarhourcheck + 2, 5) = Round(Worksheets("EPW").Cells(i, 14), 0)
            If Round(Worksheets("EPW").Cells(i, 14), 0) - Round(Worksheets("Morphed Weather").Cells(i, 11), 0) > 15 Then overridecheck = overridecheck + 1
        End If
    Next i
    If Recheckvalue = 0 Then solarhourfinal = solarhourcheck + 5
    If Recheckvalue = 1 Then solarhourfinal = solarhourfinal + solarhourcheck + 2
    For i = 11 To 8770
    ' Displays the progress bar
        PctDone = (i + 8760) / 17520
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        kstep = 1
    ' Checks whether the Global horizontal radiation data is shifted down by one row, i.e. appearing earlier than expected
        If Worksheets("Morphed Weather").Cells(i, 11) < Worksheets("EPW").Cells(i, 14) Then
            Do While Worksheets("Morphed Weather").Cells(i + kstep, 11) > 0
                If Worksheets("EPW").Cells(i + kstep, 14) = 0 And Worksheets("Morphed Weather").Cells(i + kstep, 11) > 10 Then hourshiftdown = hourshiftdown + 1
                kstep = kstep + 1
            Loop
            i = i + kstep
        End If
    ' Checks whether the Global horizontal radiation data is shifted up by one row, i.e. appearing later than expected
        If Worksheets("Morphed Weather").Cells(i, 11) > 5 And Worksheets("EPW").Cells(i, 14) > 5 Then
            Do While Worksheets("EPW").Cells(i + kstep, 14) > 0
                If Worksheets("Morphed Weather").Cells(i + kstep, 11) = 0 And Worksheets("EPW").Cells(i + kstep, 14) > 10 Then hourshiftup = hourshiftup + 1
                kstep = kstep + 1
            Loop
            i = i + kstep
        End If
    Next i

    ' Returns the radiation check result if any inconsistencies are found (Code for internal checking only, therefore commented out)
    ' If solarhour < 10 Then solarhour_calc = 0
    ' If solarhourcheck < 10 Then solarhourcheck_calc = 0
    ' If solarhour + solarhourcheck + hourshiftdown + hourshiftup > 0 Then
    '    MsgBox "There are solar radiation data inconsistencies in your file. These are as follows:" & Chr(10) & Chr(10) & "Hours global horizontal radiation appears earlier than expected: " & hourshiftdown & Chr(10) & "Hours global horizontal radiation appears later than expected: " & hourshiftup & Chr(10) & "Hours where the first/last hour of the day is above extraterrestrial radiation: " & solarhour & Chr(10) & "Total number of hours above extraterrestrial radiation: " & solarhourcheck & Chr(10) & "Total number of hours with difference > 15 Wh/m²: " & overridecheck
    ' End If

    ' Fix for data early by one hour, i.e. shift data down
    If hourshiftdown > 5 And hourshiftup = 0 And Recheckvalue = 0 Then
        Attention_shift_down.Show
    ' Exits morphing routine for individual assessment by user if cancel was selected
        If Worksheets("Selection Data").Range("Q43") = 0 Then
            Application.ScreenUpdating = True
            Application.Cursor = xlDefault
            Unload Status_M
            GoTo Logfile
        End If
    ' Tries to amend data inconsistencies
        If Worksheets("Selection Data").Range("Q43") = 1 Then
            MsgBox "Warning: This action will change the relation of the solar radiation data to the remaining data of your EPW file. The EPW solar radition data will be recalculated to match the local standard time convention. However, the change from solar time to local standard time is an assumption. The reasons for the inconsistencies between the calculated extraterrestrial horizontal radiation and your particular EPW file are not known. It is your own responsibility to check whether the file remains meteorologically consistent." & Chr(10) & Chr(10) & "Solar radiation data consistency will be re-examined after this action. It is possible that the inconsistencies cannot be resolved. Please beware that if major inconsistencies remain morphing will not be possible.", vbExclamation
            TimeShift = (LongitudeCalc / 15) - TimeZone
            If TimeShift > 0 Then TimeShift = TimeShift * -1
            TimeShiftStep = Int(TimeShift) * -1 - 1
            TimeShiftAdj = TimeShift + TimeShiftStep
        ' Starts the amendment
        For i = 11 To 8770
            ParaNameDisplay = "Global Horizontal Radiation Adjustment (Wh/m²)"
        ' Displays the progress bar
            PctDone = i / 8760
            With Status_M
                .CurrentFile.Text = ParaNameDisplay
                .FrameProgress.Caption = Format(PctDone, "0%")
                .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
            End With
            DoEvents
        ' Performs the adjustment
        ' Writes back the LST corresponding to the SOT of the original EPW global horizontal radiation
            Worksheets("Selection Data").Cells(i, 26) = Worksheets("Selection Data").Cells(i, 25) - TimeShift
            If Worksheets("Selection Data").Cells(i, 26) > 24 Then Worksheets("Selection Data").Cells(i, 26) = Worksheets("Selection Data").Cells(i, 26) - 24
            If Worksheets("Selection Data").Cells(i, 26) < 0 Then Worksheets("Selection Data").Cells(i, 26) = Worksheets("Selection Data").Cells(i, 26) + 24
        ' Calculates the solar baseline parameters
            DayofYear = Worksheets("Morphed Weather").Cells(i, 33)
            VarB = DayofYear * 360 / 365.25
            VarE = -0.128 * Sin((VarB - 2.8) * PII / 180) - 0.165 * Sin((2 * VarB + 19.7) * PII / 180)
            nosunminutes = 0
        ' Calculates sunrise and sunset in LST
            For k = 1 To 60
                LocalTime = (Worksheets("Morphed Weather").Cells(i, 4) - 1 + k / 60) + (LongitudeCalc - (TimeZone * 15)) / 15 + VarE
                HourAngle = 15 * (LocalTime - 12)
                Declination = Application.WorksheetFunction.Asin(0.3978 * Sin((VarB * PII / 180) - 1.4 + 0.0355 * Sin((VarB * PII / 180) - 0.0489)))
                sin_SolarAlt = Sin(Declination) * Sin(Worksheets("Selection Data").Range("R33") * PII / 180) + _
                    Cos(Declination) * Cos(Worksheets("Selection Data").Range("R33") * PII / 180) * _
                    Cos(HourAngle * PII / 180)
                SolarAlt = Application.WorksheetFunction.Asin(sin_SolarAlt) * 180 / PII
                If SolarAlt < 0 Then
                    nosunminutes = nosunminutes + 1
                End If
            Next k
            Worksheets("Selection Data").Cells(i, 27) = nosunminutes / 60
            If Worksheets("Selection Data").Cells(i, 27) = 1 Then Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i, 25)
            If Worksheets("Selection Data").Cells(i, 27) = 0 Then Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i, 25)
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i, 27) < 1 And Worksheets("Selection Data").Cells(i - 1, 27) = 1 Then
                 Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i - 1, 25) + Worksheets("Selection Data").Cells(i, 27)
            End If
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i, 27) < 1 And Worksheets("Selection Data").Cells(i - 1, 27) = 0 Then
                 Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i, 25) - Worksheets("Selection Data").Cells(i, 27)
            End If
        ' Adjusts the radation data
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i + 1, 27) = 0 Then
                If Worksheets("Selection Data").Cells(i, 27) = 1 Then
                    DaySunrise = 1
                Else
                    DaySunrise = 1 - Worksheets("Selection Data").Cells(i, 27)
                End If
            End If
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i - 1, 27) = 0 Then
                If Worksheets("Selection Data").Cells(i, 27) = 1 Then
                    DaySunset = 0
                Else
                    DaySunset = 1 - Worksheets("Selection Data").Cells(i, 27)
                End If
            End If
            If i < 14 Then
                AdjustRad = 0
            Else
                AdjustRad = (Worksheets("EPW").Cells(i - TimeShiftStep, 14) * (TimeShiftAdj + 1)) + (Worksheets("EPW").Cells(i - 1 - TimeShiftStep, 14) * (TimeShiftAdj * -1))
            End If
            If AdjustRad > 0 And Worksheets("Selection Data").Cells(i - 1, 30) = 0 Then
                AdjustRad = AdjustRad * DaySunrise
            End If
            Worksheets("Selection Data").Cells(i, 30) = AdjustRad
            If AdjustRad = 0 And i > 13 And Worksheets("Selection Data").Cells(i - 1, 30) > 0 Then
                AdjustRad = Worksheets("Selection Data").Cells(i - 1, 30) * DaySunset
                Worksheets("Selection Data").Cells(i - 1, 30) = AdjustRad
            End If
        Next i
            Worksheets("Selection Data").Range("AD11:AD8770").Copy Destination:=Worksheets("EPW").Range("N11:N8770")
        End If
        Recheckvalue = 1
        GoTo Recheck
    End If

    ' Fix for data late by one hour, i.e. shift data up
    If hourshiftdown = 0 And hourshiftup > 5 And Recheckvalue = 0 Then
        Attention_shift_up.Show
    ' Exits morphing routine for individual assessment by user if cancel was selected
        If Worksheets("Selection Data").Range("Q43") = 0 Then
            Application.ScreenUpdating = True
            Application.Cursor = xlDefault
            Unload Status_M
            GoTo Logfile
        End If
    ' Tries to amend data inconsistencies
        If Worksheets("Selection Data").Range("Q43") = 1 Then
            MsgBox "Warning: This action will change the relation of the solar radiation data to the remaining data of your EPW file. The EPW solar radition data will be recalculated to match the local standard time convention. However, the change from solar time to local standard time is an assumption. The reasons for the inconsistencies between the calculated extraterrestrial horizontal radiation and your particular EPW file are not known. It is your own responsibility to check whether the file remains meteorologically consistent." & Chr(10) & Chr(10) & "Solar radiation data consistency will be re-examined after this action. It is possible that the inconsistencies cannot be resolved. Please beware that if major inconsistencies remain morphing will not be possible.", vbExclamation
            TimeShift = (LongitudeCalc / 15) - TimeZone
            If TimeShift < 0 Then TimeShift = TimeShift * -1
            TimeShiftStep = Int(TimeShift) * -1 - 1
            TimeShiftAdj = TimeShift + TimeShiftStep
        ' Starts the amendment
        For i = 11 To 8770
            ParaNameDisplay = "Global Horizontal Radiation Adjustment (Wh/m²)"
        ' Displays the progress bar
            PctDone = i / 8760
            With Status_M
                .CurrentFile.Text = ParaNameDisplay
                .FrameProgress.Caption = Format(PctDone, "0%")
                .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
            End With
            DoEvents
        ' Performs the adjustment
        ' Writes back the LST corresponding to the SOT of the original EPW global horizontal radiation
            Worksheets("Selection Data").Cells(i, 26) = Worksheets("Selection Data").Cells(i, 25) - TimeShift
            If Worksheets("Selection Data").Cells(i, 26) > 24 Then Worksheets("Selection Data").Cells(i, 26) = Worksheets("Selection Data").Cells(i, 26) - 24
            If Worksheets("Selection Data").Cells(i, 26) < 0 Then Worksheets("Selection Data").Cells(i, 26) = Worksheets("Selection Data").Cells(i, 26) + 24
        ' Calculates the solar baseline parameters
            DayofYear = Worksheets("Morphed Weather").Cells(i, 33)
            VarB = DayofYear * 360 / 365.25
            VarE = -0.128 * Sin((VarB - 2.8) * PII / 180) - 0.165 * Sin((2 * VarB + 19.7) * PII / 180)
            nosunminutes = 0
        ' Calculates sunrise and sunset in LST
            For k = 1 To 60
                LocalTime = (Worksheets("Morphed Weather").Cells(i, 4) - 1 + k / 60) + (LongitudeCalc - (TimeZone * 15)) / 15 + VarE
                HourAngle = 15 * (LocalTime - 12)
                Declination = Application.WorksheetFunction.Asin(0.3978 * Sin((VarB * PII / 180) - 1.4 + 0.0355 * Sin((VarB * PII / 180) - 0.0489)))
                sin_SolarAlt = Sin(Declination) * Sin(Worksheets("Selection Data").Range("R33") * PII / 180) + _
                    Cos(Declination) * Cos(Worksheets("Selection Data").Range("R33") * PII / 180) * _
                    Cos(HourAngle * PII / 180)
                SolarAlt = Application.WorksheetFunction.Asin(sin_SolarAlt) * 180 / PII
                If SolarAlt < 0 Then
                    nosunminutes = nosunminutes + 1
                End If
            Next k
            Worksheets("Selection Data").Cells(i, 27) = nosunminutes / 60
            If Worksheets("Selection Data").Cells(i, 27) = 1 Then Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i, 25)
            If Worksheets("Selection Data").Cells(i, 27) = 0 Then Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i, 25)
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i, 27) < 1 And Worksheets("Selection Data").Cells(i - 1, 27) = 1 Then
                 Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i - 1, 25) + Worksheets("Selection Data").Cells(i, 27)
            End If
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i, 27) < 1 And Worksheets("Selection Data").Cells(i - 1, 27) = 0 Then
                 Worksheets("Selection Data").Cells(i, 28) = Worksheets("Selection Data").Cells(i, 25) - Worksheets("Selection Data").Cells(i, 27)
            End If
        ' Adjusts the radation data
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i + 1, 27) = 0 Then
                If Worksheets("Selection Data").Cells(i, 27) = 1 Then
                    DaySunrise = 1
                Else
                    DaySunrise = 1 - Worksheets("Selection Data").Cells(i, 27)
                End If
            End If
            If Worksheets("Selection Data").Cells(i, 27) > 0 And Worksheets("Selection Data").Cells(i - 1, 27) = 0 Then
                If Worksheets("Selection Data").Cells(i, 27) = 1 Then
                    DaySunset = 0
                Else
                    DaySunset = 1 - Worksheets("Selection Data").Cells(i, 27)
                End If
            End If
            If i < 14 Then
                AdjustRad = 0
            Else
                AdjustRad = (Worksheets("EPW").Cells(i - TimeShiftStep, 14) * (TimeShiftAdj + 1)) + (Worksheets("EPW").Cells(i - 1 - TimeShiftStep, 14) * (TimeShiftAdj * -1))
            End If
            If AdjustRad > 0 And Worksheets("Selection Data").Cells(i - 1, 30) = 0 Then
                AdjustRad = AdjustRad * DaySunrise
            End If
            Worksheets("Selection Data").Cells(i, 30) = AdjustRad
            If AdjustRad = 0 And i > 13 And Worksheets("Selection Data").Cells(i - 1, 30) > 0 Then
                AdjustRad = Worksheets("Selection Data").Cells(i - 1, 30) * DaySunset
                Worksheets("Selection Data").Cells(i - 1, 30) = AdjustRad
            End If
        Next i
            Worksheets("Selection Data").Range("AD11:AD8770").Copy Destination:=Worksheets("EPW").Range("N11:N8770")
        End If
        Recheckvalue = 1
        GoTo Recheck
    End If

    ' Terminates the procedures if the re-examination fails or if there is a shift up as well as down
    If hourshiftdown > 14 And hourshiftup = 0 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown > 3 And hourshiftup = 1 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown > 3 And hourshiftup = 2 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown > 3 And hourshiftup = 3 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown = 0 And hourshiftup > 14 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown = 1 And hourshiftup > 3 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown = 2 And hourshiftup > 3 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown = 3 And hourshiftup > 3 And Recheckvalue = 1 Then RadCalcFail = 1
    If hourshiftdown > 0 And hourshiftup > 0 And Recheckvalue = 0 Then RadCalcFail = 1
    If hourshiftdown > 3 And hourshiftup > 3 And Recheckvalue = 1 Then RadCalcFail = 1
    If overridecheck > 19 Then RadCalcFail = 1
    If RadCalcFail = 1 Then
        Attention_terminate.Show
        Application.ScreenUpdating = True
        Application.Cursor = xlDefault
        Unload Status_M
        Worksheets("Convert File").Activate
        Worksheets("Convert File").Range("B50").Value = "    No morphed weather file"
        Worksheets("Convert File").Range("B51").Value = " "
        Worksheets("Morphed Weather").Range("A1:BA9000").Clear
        Worksheets("Selection Data").Range("AC11:AC8770").Copy Destination:=Worksheets("EPW").Range("N11:N8770")
        Worksheets("Selection Data").Range("Q43") = 0
        GoTo Logfile
    End If

    ' Explanantion on corrections to global horizontal radiation on the first / last hour of the day
    If solarhour > 9 And Recheckvalue = 0 Then
        Attention_radiation_firstlast.Show
    ' Exits morphing routine for individual assessment by user if cancel was selected
        If Worksheets("Selection Data").Range("Q43") = 0 Then
            Application.ScreenUpdating = True
            Application.Cursor = xlDefault
            Unload Status_M
            GoTo Logfile
        End If
        If Worksheets("Selection Data").Range("Q43") = 1 Then
            Do
                b = 1
                A = MsgBox("Warning: This action will amend global horizontal radiation data during morphing where it exceeds the extraterrestrial horizontal radiation. This was found to be the case on a number of occasions for the first/last hour of daylight in your original EPW data. It is your own responsibility to check whether these changes are acceptable and whether the morphed weather file remains meteorologically consistent." & Chr(10) & Chr(10) & "Do you want to generate a radiation log file now highlighting the hours with inconsistencies (recommended)?", vbExclamation + vbYesNo)
                If A = 7 Then
                    GoTo TotHorRad_Continue
                Else
                    GoTo Logfile
                End If
            Loop While b > 1
        End If
    End If

    'Lets the user save the radiation log file in case an amendment from SOT to LST has been successfully conducted
    If Recheckvalue = 1 Then
        Do
            b = 1
            A = MsgBox("Shifting the solar radiation parameters of the EPW file helped to resolve the inconsistencies highlighted before. However, there are still a number of occasions where the re-calculated global horizontal radiation exceeds the extraterrestrial horizontal radiation. This will be amended during the morphing process. It is your own responsibility to check whether these changes are acceptable and whether the morphed weather file remains meteorologically consistent." & Chr(10) & Chr(10) & "Do you want to generate a radiation log file now highlighting the initial and the current hours with inconsistencies (recommended)?", vbExclamation + vbYesNo)
            If A = 7 Then
                GoTo TotHorRad_Continue
            Else
                GoTo Logfile
            End If
        Loop While b > 1
    End If

TotHorRad_Continue:
' Total horizontal solar irradiance: Scaling using the baseline climate as scaling factor (GlRad / DSWF)
    Dim GlRad_DSWF, MonthDays As Double
    ParaNameDisplay = "Global Horizontal Radiation (Wh/m²)"
    For i = 11 To 8770
    ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
    'Performs the morphing
        For j = 1 To 12
            If Worksheets("Morphed Weather").Cells(i, 2).Value = j Then
            MonthDays = 24 * Worksheets("Selection Data").Cells(j + 16, 19)
                Scalingfactor = 1 + (Worksheets("Convert File").Cells(16, j + 4) / (Worksheets("Selection Data").Cells(j + 16, 18) * 1000 / MonthDays))
                GlRad_DSWF = Scalingfactor * Worksheets("EPW").Cells(i, 14)
                If GlRad_DSWF > 1200 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i, 11) * 0.87
                If GlRad_DSWF < 0 Then GlRad_DSWF = 0
                If Worksheets("Morphed Weather").Cells(i, 11) < 0.5 Then
                    GlRad_DSWF = 0
                    Worksheets("Morphed Weather").Cells(i, 11) = 0
                    Worksheets("Morphed Weather").Cells(i, 12) = 0
                End If
            ' Correction to account for incorrect global horizontal radiation data in the original EPW file
                If Worksheets("Morphed Weather").Cells(i, 11) < GlRad_DSWF And Worksheets("Morphed Weather").Cells(i, 4) > 10 And Worksheets("Morphed Weather").Cells(i, 4) < 16 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i, 11) * 0.85
                If Worksheets("Morphed Weather").Cells(i, 11) < GlRad_DSWF And Worksheets("Morphed Weather").Cells(i, 4) < 11 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i, 11) * 0.6
                If Worksheets("Morphed Weather").Cells(i, 11) < GlRad_DSWF And Worksheets("Morphed Weather").Cells(i, 4) > 15 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i, 11) * 0.6
                Worksheets("Morphed Weather").Cells(i, 14) = GlRad_DSWF
                If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                    Smooth = (Worksheets("Convert File").Cells(16, j + 5) - Worksheets("Convert File").Cells(16, j + 4)) / 48
                    For k = 1 To 48
                        Smoothstep = Smooth * k
                        Scalingfactor = 1 + ((Worksheets("Convert File").Cells(16, j + 4) + Smoothstep) / (Worksheets("Selection Data").Cells(j + 16, 18) * 1000 / MonthDays))
                        GlRad_DSWF = Scalingfactor * Worksheets("EPW").Cells(i - 24 + k, 14)
                        If GlRad_DSWF > 1200 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i - 24 + k, 11) * 0.87
                        If GlRad_DSWF < 0 Then GlRad_DSWF = 0
                        If Worksheets("Morphed Weather").Cells(i - 24 + k, 11) = 0 Then
                            GlRad_DSWF = 0
                        End If
                        If Worksheets("Morphed Weather").Cells(i - 24 + k, 11) < GlRad_DSWF And Worksheets("Morphed Weather").Cells(i - 24 + k, 4) > 10 And Worksheets("Morphed Weather").Cells(i - 24 + k, 4) < 16 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i - 24 + k, 11) * 0.85
                        If Worksheets("Morphed Weather").Cells(i - 24 + k, 11) < GlRad_DSWF And Worksheets("Morphed Weather").Cells(i - 24 + k, 4) < 11 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i - 24 + k, 11) * 0.6
                        If Worksheets("Morphed Weather").Cells(i - 24 + k, 11) < GlRad_DSWF And Worksheets("Morphed Weather").Cells(i - 24 + k, 4) > 15 Then GlRad_DSWF = Worksheets("Morphed Weather").Cells(i - 24 + k, 11) * 0.6
                       Worksheets("Morphed Weather").Cells(i - 24 + k, 14) = GlRad_DSWF
                    Next k
                i = i + 24
                End If
            End If
        Next j
    Next i

' Diffuse horizontal solar irradiance: BRL model by Boland, Ridley and Laurent (DiRad)
    Dim DiRad, ClearIndexHour, ClearIndexDay, ClearPersistence, DiRad_Check_A, DiRad_Check_B, sin_SolarAlt_step, sin_SolarAlt_sum As Double
    Dim DayId, minutes As Integer
    ParaNameDisplay = "Diffuse Horizontal Radiation (Wh/m²)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        'Performs the morphing
            If Worksheets("Morphed Weather").Cells(i, 11) > 0 Then
                ClearIndexHour = Worksheets("Morphed Weather").Cells(i, 14) / Worksheets("Morphed Weather").Cells(i, 11)
                DayId = Worksheets("Morphed Weather").Cells(i, 33)
                ClearIndexDay = Worksheets("Selection Data").Cells(DayId + 11, 23)
                PII = Application.WorksheetFunction.Pi
                LongitudeCalc = Worksheets("Selection Data").Range("R34")
                TimeZone = Worksheets("Selection Data").Range("R36")
                DayofYear = Worksheets("Morphed Weather").Cells(i, 33)
                VarB = DayofYear * 360 / 365.25
                VarE = -0.128 * Sin((VarB - 2.8) * PII / 180) - 0.165 * Sin((2 * VarB + 19.7) * PII / 180)
                Declination = Application.WorksheetFunction.Asin(0.3978 * Sin((VarB * PII / 180) - 1.4 + 0.0355 * Sin((VarB * PII / 180) - 0.0489)))
                'new
                sin_SolarAlt_sum = 0
                minutes = 0
                For k = 1 To 60
                    LocalTime = (Worksheets("Morphed Weather").Cells(i, 4) - 1 + k / 60) + (LongitudeCalc - (TimeZone * 15)) / 15 + VarE
                    HourAngle = 15 * (LocalTime - 12)
                    sin_SolarAlt_step = Sin(Declination) * Sin(Worksheets("Selection Data").Range("R33") * PII / 180) + _
                        Cos(Declination) * Cos(Worksheets("Selection Data").Range("R33") * PII / 180) * _
                        Cos(HourAngle * PII / 180)
                    If sin_SolarAlt_step > 0 Then
                        sin_SolarAlt_sum = sin_SolarAlt_sum + sin_SolarAlt_step
                        minutes = minutes + 1
                    End If
                Next k
                If minutes = 0 Then
                    SolarAlt = 0
                Else
                    sin_SolarAlt = sin_SolarAlt_sum / minutes
                    SolarAlt = Application.WorksheetFunction.Asin(sin_SolarAlt) * 180 / PII
                End If
                If SolarAlt < 0 Then SolarAlt = 0
                ' new end
                ' calculates the clearness persistence, i.e. inertia
                    If Worksheets("Morphed Weather").Cells(i - 1, 11) = 0 Then ClearPersistence = Worksheets("Morphed Weather").Cells(i + 1, 14) / Worksheets("Morphed Weather").Cells(i + 1, 11)
                    If Worksheets("Morphed Weather").Cells(i + 1, 11) = 0 Then ClearPersistence = Worksheets("Morphed Weather").Cells(i - 1, 14) / Worksheets("Morphed Weather").Cells(i - 1, 11)
                    If Worksheets("Morphed Weather").Cells(i - 1, 11) > 0 And Worksheets("Morphed Weather").Cells(i + 1, 11) > 0 Then
                        ClearPersistence = (Worksheets("Morphed Weather").Cells(i - 1, 14) / Worksheets("Morphed Weather").Cells(i - 1, 11) + _
                            Worksheets("Morphed Weather").Cells(i + 1, 14) / Worksheets("Morphed Weather").Cells(i + 1, 11)) / 2
                    End If
                DiRad = Worksheets("Morphed Weather").Cells(i, 14) * (1 / (1 + Exp(-5.38 + 6.63 * ClearIndexHour + 0.006 * LocalTime - 0.007 * SolarAlt + 1.75 * ClearIndexDay + 1.31 * ClearPersistence)))
            Else
                DiRad = 0
            End If
            If Worksheets("Selection Data").Cells(i - TimeShiftStep, 29) > 0 And Worksheets("Selection Data").Cells(i - TimeShiftStep, 29) = Worksheets("EPW").Cells(i - TimeShiftStep, 16) Then
                DiRad_Check_A = DiRad + Worksheets("EPW").Cells(i - TimeShiftStep - 1, 16)
                If DiRad_Check_A = DiRad Then DiRad = Worksheets("Morphed Weather").Cells(i, 14)
                DiRad_Check_B = DiRad + Worksheets("EPW").Cells(i - TimeShiftStep + 1, 16)
                If DiRad_Check_B = DiRad Then DiRad = Worksheets("Morphed Weather").Cells(i, 14)
            End If
            If Worksheets("Morphed Weather").Cells(i, 11) = 0 Then
                DiRad = 0
            End If
            If DiRad > 700 Then DiRad = 700
            If DiRad < 0 Then DiRad = 0
            Worksheets("Morphed Weather").Cells(i, 16) = DiRad
    Next i

' Calculates the lighting parameters
    For i = 11 To 8770
    ParaNameDisplay = "Direct Normal Radiation & Lighting Parameters"
    ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
    ' Performs the morphing
    ' Calculates Direct Normal Radiation Radiation from solar elevation and irradiance data
        Dim DirRad, NorRad As Double
        PII = Application.WorksheetFunction.Pi
        LongitudeCalc = Worksheets("Selection Data").Range("R34")
        TimeZone = Worksheets("Selection Data").Range("R36")
        DayofYear = Worksheets("Morphed Weather").Cells(i, 33)
        VarB = DayofYear * 360 / 365.25
        VarE = -0.128 * Sin((VarB - 2.8) * PII / 180) - 0.165 * Sin((2 * VarB + 19.7) * PII / 180)
        sin_SolarAlt_sum = 0
        minutes = 0
        For k = 1 To 60
            LocalTime = (Worksheets("Morphed Weather").Cells(i, 4) - 1 + k / 60) + (LongitudeCalc - (TimeZone * 15)) / 15 + VarE
            HourAngle = 15 * (LocalTime - 12)
            Declination = Application.WorksheetFunction.Asin(0.3978 * Sin((VarB * PII / 180) - 1.4 + 0.0355 * Sin((VarB * PII / 180) - 0.0489)))
            sin_SolarAlt_step = Sin(Declination) * Sin(Worksheets("Selection Data").Range("R33") * PII / 180) + _
                Cos(Declination) * Cos(Worksheets("Selection Data").Range("R33") * PII / 180) * _
                Cos(HourAngle * PII / 180)
            If sin_SolarAlt_step > 0 Then
                sin_SolarAlt_sum = sin_SolarAlt_sum + sin_SolarAlt_step
                minutes = minutes + 1
            End If
        Next k
        If minutes = 0 Then
            NorRad = 0
            SolarAlt = 0
            SolarAltR = 0
        Else
            sin_SolarAlt = sin_SolarAlt_sum / minutes
            SolarAlt = Application.WorksheetFunction.Asin(sin_SolarAlt) * 180 / PII
            DirRad = Worksheets("Morphed Weather").Cells(i, 14) - Worksheets("Morphed Weather").Cells(i, 16)
            NorRad = (minutes / 60) * (DirRad / sin_SolarAlt)
            If NorRad > 1100 Then NorRad = 1100
            If NorRad < 0 Then NorRad = 0
            If Worksheets("Morphed Weather").Cells(i, 14) < 0.5 Then NorRad = 0
        End If
        Worksheets("Morphed Weather").Cells(i, 15) = NorRad
    ' Calcualtes the illuminance parameters and zenith luminance
        Dim GlHorIll, DirNorIll, DifHorIll, ZenLum As Double
        Dim Fudgerange As Integer
        Dim HorDifRad, ZenAng, SkyClear, SkyBright, PreciWaterCont As Double
        Dim AirMass, CoeffiA, CoeffiB, CoeffiC, CoeffiD As Double
        If Worksheets("Morphed Weather").Cells(i, 14) < 0.5 Then
            GlHorIll = 0
            DirNorIll = 0
            DifHorIll = 0
            ZenLum = 0
        Else
            HorDifRad = Worksheets("Morphed Weather").Cells(i, 16)
            ZenAng = (90 - (minutes / 60) * SolarAlt) * PII / 180
            SkyClear = Round(((HorDifRad + NorRad) / HorDifRad + 1.041 * (ZenAng) ^ 3) / (1 + 1.041 * (ZenAng) ^ 3), 3)
            SolarAltR = Round((minutes / 60) * SolarAlt, 1)
            If SolarAltR < 0 Then SolarAltR = 0
            AirMass = Worksheets("Illuminance").Cells(SolarAltR * 10 + 2, 9)
            SkyBright = HorDifRad * AirMass / Worksheets("Morphed Weather").Cells(i, 12)
            PreciWaterCont = Exp(0.07 * Worksheets("Morphed Weather").Cells(i, 8) - 0.075)
            If SkyClear < 1.065 Then
                Fudgerange = 1
            End If
            If SkyClear > 1.064 And SkyClear < 1.23 Then
                Fudgerange = 2
            End If
            If SkyClear > 1.229 And SkyClear < 1.5 Then
                Fudgerange = 3
            End If
            If SkyClear > 1.499 And SkyClear < 1.95 Then
                Fudgerange = 4
            End If
            If SkyClear > 1.949 And SkyClear < 2.8 Then
                Fudgerange = 5
            End If
            If SkyClear > 2.799 And SkyClear < 4.5 Then
                Fudgerange = 6
            End If
            If SkyClear > 4.499 And SkyClear < 6.2 Then
                Fudgerange = 7
            End If
            If SkyClear > 6.199 Then
                Fudgerange = 8
            End If
        ' Global Horizontal Illuminance
            CoeffiA = Worksheets("Illuminance").Cells(Fudgerange + 2, 2)
            CoeffiB = Worksheets("Illuminance").Cells(Fudgerange + 2, 3)
            CoeffiC = Worksheets("Illuminance").Cells(Fudgerange + 2, 4)
            CoeffiD = Worksheets("Illuminance").Cells(Fudgerange + 2, 5)
            GlHorIll = Worksheets("Morphed Weather").Cells(i, 14) * (CoeffiA + CoeffiB * PreciWaterCont + CoeffiC * Cos(ZenAng) + CoeffiD * Application.WorksheetFunction.Ln(SkyBright))
            If GlHorIll > 130000 Then GlHorIll = 130000
            If GlHorIll < 0 Then GlHorIll = 0
        ' Direct Normal Illuminance
            CoeffiA = Worksheets("Illuminance").Cells(Fudgerange + 22, 2)
            CoeffiB = Worksheets("Illuminance").Cells(Fudgerange + 22, 3)
            CoeffiC = Worksheets("Illuminance").Cells(Fudgerange + 22, 4)
            CoeffiD = Worksheets("Illuminance").Cells(Fudgerange + 22, 5)
            DirNorIll = NorRad * (CoeffiA + CoeffiB * PreciWaterCont + CoeffiC * Exp(5.73 * ZenAng - 5) + CoeffiD * SkyBright)
            If DirNorIll > 110000 Then DirNorIll = 110000
            If DirNorIll < 0 Then DirNorIll = 0
        ' Diffuse Horizontal Illuminance
            CoeffiA = Worksheets("Illuminance").Cells(Fudgerange + 12, 2)
            CoeffiB = Worksheets("Illuminance").Cells(Fudgerange + 12, 3)
            CoeffiC = Worksheets("Illuminance").Cells(Fudgerange + 12, 4)
            CoeffiD = Worksheets("Illuminance").Cells(Fudgerange + 12, 5)
            DifHorIll = HorDifRad * (CoeffiA + CoeffiB * PreciWaterCont + CoeffiC * Cos(ZenAng) + CoeffiD * Application.WorksheetFunction.Ln(SkyBright))
            If DifHorIll > 80000 Then DifHorIll = 80000
            If DifHorIll < 0 Then DifHorIll = 0
        ' Zenith Luminance
            CoeffiA = Worksheets("Illuminance").Cells(Fudgerange + 32, 2)
            CoeffiB = Worksheets("Illuminance").Cells(Fudgerange + 32, 3)
            CoeffiC = Worksheets("Illuminance").Cells(Fudgerange + 32, 4)
            CoeffiD = Worksheets("Illuminance").Cells(Fudgerange + 32, 5)
            ZenLum = HorDifRad * (CoeffiA + CoeffiB * Cos(ZenAng) + CoeffiC * Exp(-3 * ZenAng) + CoeffiD * SkyBright)
            If ZenLum > 70000 Then ZenLum = 70000
            If ZenLum < 0 Then ZenLum = 0
        End If
        If DifHorIll > GlHorIll Then DifHorIll = GlHorIll
        If DirNorIll = 0 Then DifHorIll = GlHorIll
        If Round(DifHorIll, 0) = Round(GlHorIll, 0) Then DirNorIll = 0
        Worksheets("Morphed Weather").Cells(i, 17) = GlHorIll
        Worksheets("Morphed Weather").Cells(i, 18) = DirNorIll
        Worksheets("Morphed Weather").Cells(i, 19) = DifHorIll
        Worksheets("Morphed Weather").Cells(i, 20) = ZenLum
    Next i
    Worksheets("Selection Data").Range("AC11:AC8770").Copy Destination:=Worksheets("EPW").Range("N11:N8770")

' Wind speed: Simple monthly stretching (WindSp / WIND)
    Dim WS_WIND As Double
    ParaNameDisplay = "Wind Speed (m/s)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Performs the morphing
        For j = 1 To 12
            If Worksheets("Morphed Weather").Cells(i, 2).Value = j Then
                WS_WIND = (1 + Worksheets("Convert File").Cells(21, j + 4) / 100) * Worksheets("EPW").Cells(i, 22)
                If WS_WIND > 40 Then WS_WIND = 40
                Worksheets("Morphed Weather").Cells(i, 22) = WS_WIND
                If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                    Smooth = (Worksheets("Convert File").Cells(21, j + 5) - Worksheets("Convert File").Cells(21, j + 4)) / 48
                    For k = 1 To 48
                        Smoothstep = Smooth * k
                        WS_WIND = (1 + (Worksheets("Convert File").Cells(21, j + 4) + Smoothstep) / 100) * Worksheets("EPW").Cells(i - 24 + k, 22)
                        If WS_WIND > 40 Then WS_WIND = 40
                        Worksheets("Morphed Weather").Cells(i - 24 + k, 22) = WS_WIND
                    Next k
                i = i + 24
                End If
            End If
        Next j
    Next i

' Precipitation: Simple monthly stretching (Prec / PREC)
    Dim Prec_PREC As Double
    If Application.WorksheetFunction.Sum(Worksheets("EPW").Range("AC11:AC8770")) = 0 Or Application.WorksheetFunction.Sum(Worksheets("EPW").Range("AC11:AC8770")) = 867240 Or Application.WorksheetFunction.Sum(Worksheets("EPW").Range("AC11:AC8770")) = 8751240 Then
        Worksheets("Morphed Weather").Range("AC11:AC8770") = 0
    Else
    ParaNameDisplay = "Precipitation (mm)"
    For i = 11 To 8770
        ' Displays the progress bar
        PctDone = i / 8760
        With Status_M
            .CurrentFile.Text = ParaNameDisplay
            .FrameProgress.Caption = Format(PctDone, "0%")
            .LabelProgress.Width = PctDone * (.FrameProgress.Width - 16)
        End With
        DoEvents
        ' Performs the morphing
        For j = 1 To 12
            If Worksheets("Morphed Weather").Cells(i, 2).Value = j Then
                Prec_PREC = (1 + Worksheets("Convert File").Cells(18, j + 4) / 100) * Worksheets("EPW").Cells(i, 29)
                If Prec_PREC > 100 Then Prec_PREC = 100
                Worksheets("Morphed Weather").Cells(i, 29) = Prec_PREC
                If Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) < 0 And Worksheets("Morphed Weather").Cells(i + 1, 3) - Worksheets("Morphed Weather").Cells(i, 3) > -31 Then
                    Smooth = (Worksheets("Convert File").Cells(18, j + 5) - Worksheets("Convert File").Cells(18, j + 4)) / 48
                    For k = 1 To 48
                        Smoothstep = Smooth * k
                        Prec_PREC = (1 + (Worksheets("Convert File").Cells(18, j + 4) + Smoothstep) / 100) * Worksheets("EPW").Cells(i - 24 + k, 29)
                        If Prec_PREC > 100 Then Prec_PREC = 100
                        Worksheets("Morphed Weather").Cells(i - 24 + k, 29) = Prec_PREC
                    Next k
                i = i + 24
                End If
            End If
        Next j
    Next i
    End If

' Missing and impossible fields
    With Worksheets("Morphed Weather")
        .Range("Y11:Y8770") = 9999                  ' Visibility
        .Range("Z11:Z8770") = 99999                 ' Ceiling height
        .Range("AA11:AA8770") = 9                   ' PWO
        .Range("AB11:AB8770") = 999999999           ' PWC
        .Range("AD11:AD8770") = 999                 ' Aerosol optical depth
        .Range("AE11:AE8770") = 999                 ' Snow depth
        .Range("AF11:AF8770") = 99                  ' Days since last snowfall
    End With

' Calculates the ground temperature at 3 depths for the file header
    Dim Tamb_mean, Tamb_amp, Tamb_min_day, PhaseLag, GroundTemp, BetaP, Tmaxi, Tmini, var_y, var_z As Double
    Dim GroundTempJan, GroundTempFeb, GroundTempMar, GroundTempApr, GroundTempMai, GroundTempJun, GroundTempJul, GroundTempAug, GroundTempSep, GroundTempOct, GroundTempNov, GroundTempDec As String
    Dim GroundT_int As Integer
    Dim GrDepth, GroundT_raw As Double
    Dim GroundT_Text, GroundT_dec, Temprange As String

    ' Sets the main variables required for the calculation
    GroundT_Text = "GROUND TEMPERATURES,3"
    PII = Application.WorksheetFunction.Pi
    Temprange_Copy = Worksheets("Morphed Weather").Range("G11:G8770")
    Worksheets("Selection Data").Range("I12:I8771") = Temprange_Copy
    If Worksheets("Selection Data").Range("Q12") = 1 Then Tamb_min_day = 15
    If Worksheets("Selection Data").Range("Q12") = 2 Then Tamb_min_day = 46
    If Worksheets("Selection Data").Range("Q12") = 3 Then Tamb_min_day = 74
    If Worksheets("Selection Data").Range("Q12") = 4 Then Tamb_min_day = 95
    If Worksheets("Selection Data").Range("Q12") = 5 Then Tamb_min_day = 135
    If Worksheets("Selection Data").Range("Q12") = 6 Then Tamb_min_day = 166
    If Worksheets("Selection Data").Range("Q12") = 7 Then Tamb_min_day = 196
    If Worksheets("Selection Data").Range("Q12") = 8 Then Tamb_min_day = 227
    If Worksheets("Selection Data").Range("Q12") = 9 Then Tamb_min_day = 258
    If Worksheets("Selection Data").Range("Q12") = 10 Then Tamb_min_day = 288
    If Worksheets("Selection Data").Range("Q12") = 11 Then Tamb_min_day = 319
    If Worksheets("Selection Data").Range("Q12") = 12 Then Tamb_min_day = 349
    PhaseLag = Tamb_min_day * 0.017214 + 0.341787
    Tmaxi = Application.WorksheetFunction.Max(Worksheets("Selection Data").Range("T17:T28"))
    Tmini = Application.WorksheetFunction.Min(Worksheets("Selection Data").Range("T17:T28"))
    Tamb_amp = (Tmaxi - Tmini) / 2
    Tamb_mean = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range("T17:T28"))
    ' Sets the different depths and calculates the temperatures
    For k = 1 To 3
        If k = 1 Then GrDepth = 0.5
        If k = 2 Then GrDepth = 2
        If k = 3 Then GrDepth = 4
        BetaP = GrDepth * (PII / (0.055741824 * 365)) ^ 0.5
        var_y = ((Exp(-BetaP)) ^ 2 - 2 * Exp(-BetaP) * Cos(BetaP) + 1) / (2 * BetaP ^ 2)
        var_z = (1 - Exp(-BetaP) * (Cos(BetaP) + Sin(BetaP))) / (1 - Exp(-BetaP) * (Cos(BetaP) - Sin(BetaP)))
        For i = 1 To 365
        'Calculates the ground temperature
            GroundTemp = Tamb_mean - Tamb_amp * Cos(2 * (PII / 365) * i - PhaseLag - Atn(var_z)) * var_y ^ 0.5
            Worksheets("Selection Data").Cells(i + 11, 15) = GroundTemp
        Next i
    ' Calculates the monthly averages and returns them to the headertext
    ' January
    Temprange = "O12:O42"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempJan = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempJan = GroundT_int & "." & GroundT_dec
    End If
    ' February
    Temprange = "O43:O70"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempFeb = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempFeb = GroundT_int & "." & GroundT_dec
    End If
    ' March
    Temprange = "O71:O101"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempMar = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempMar = GroundT_int & "." & GroundT_dec
    End If
    ' April
    Temprange = "O102:O131"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempApr = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempApr = GroundT_int & "." & GroundT_dec
    End If
    ' May
    Temprange = "O132:O162"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempMai = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempMai = GroundT_int & "." & GroundT_dec
    End If
    ' June
    Temprange = "O163:O192"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempJun = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempJun = GroundT_int & "." & GroundT_dec
    End If
    ' July
    Temprange = "O193:O223"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempJul = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempJul = GroundT_int & "." & GroundT_dec
    End If
    ' August
    Temprange = "O224:O254"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempAug = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempAug = GroundT_int & "." & GroundT_dec
    End If
    ' September
    Temprange = "O255:O284"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempSep = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempSep = GroundT_int & "." & GroundT_dec
    End If
    ' October
    Temprange = "O285:O315"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempOct = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempOct = GroundT_int & "." & GroundT_dec
    End If
    ' November
    Temprange = "O316:O345"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempNov = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempNov = GroundT_int & "." & GroundT_dec
    End If
    ' December
    Temprange = "O346:O376"
    GroundT_raw = Application.WorksheetFunction.Average(Worksheets("Selection Data").Range(Temprange))
    GroundT_int = Fix(Round(GroundT_raw, 1))
    GroundT_dec = Right(Round(GroundT_raw * 100), 2)
    If GroundT_raw < 0 And GroundT_int = 0 Then
        GroundTempDec = "-" & GroundT_int & "." & GroundT_dec
    Else
        GroundTempDec = GroundT_int & "." & GroundT_dec
    End If
    GroundT_Text = GroundT_Text & "," & GrDepth & ",,,," & GroundTempJan & "," & GroundTempFeb & "," & GroundTempMar & "," & GroundTempApr & "," & GroundTempMai & "," & GroundTempJun & "," & GroundTempJul & "," & GroundTempAug & "," & GroundTempSep & "," & GroundTempOct & "," & GroundTempNov & "," & GroundTempDec
    Next k
    Worksheets("Morphed Weather").Range("A4").Value = GroundT_Text

' Unloads the progress bar and returns to the main worksheet
Application.ScreenUpdating = True
Application.Cursor = xlDefault
Unload Status_M
Worksheets("Convert File").Range("B50").Value = "    Morphed EPW file for: " & Worksheets("HelpCalc").Range("C10")
Worksheets("Convert File").Range("B51").Value = "    HadCM3 A2 emissions senario ensemble for the " & _
    Worksheets("Morphed Weather").Range("A11") & "'s"
Worksheets("Convert File").Activate

End

Logfile:

' Lets the user decide whether he wants to save the radiation log file and if yes opens a new workbook

    City = Worksheets("HelpCalc").Range("C11")
    Country = Worksheets("HelpCalc").Range("C12")
    SimYear = Worksheets("HelpCalc").Range("C13")
    Do
        A = MsgBox("Do you want to save the radiation log file for " & City & "?", vbQuestion + vbYesNo)
        b = 1
        If A = 7 Then
            b = MsgBox("Your radiation log file will not be saved!", vbExclamation + vbOKCancel)
            If b = 1 Then
                If Worksheets("Selection Data").Range("Q43") = 1 Then
                    GoTo TotHorRad_Continue
                Else
                    Worksheets("Convert File").Activate
                End If
                Exit Sub
            End If
        Else
            Datareturn = Worksheets("Radiation Log").Range("A1:E8768")
            Set AddBook = Workbooks.Add
            Set AddSheet = Worksheets.Add
            With AddSheet
                .Name = City
                .Range("A1:E8768") = Datareturn
            End With
        End If
    Loop While b > 1

    ' Lets the user save the file to a specified location
    Do
        FNameString = Country & "_" & City & "_radiation_log"
        FName = Application.GetSaveAsFilename _
            (FNameString, filefilter:="LOG Files (*.log), *.log")
        If FName <> False Then
            A = MsgBox("Save as " & FName & " ?", vbQuestion + vbYesNo)
            If A = 7 Then
                MsgBox "The file was not saved!", vbExclamation + vbOKOnly
                ActiveWorkbook.Close SaveChanges:=False
                If Worksheets("Selection Data").Range("Q43") = 1 Then
                    GoTo TotHorRad_Continue
                Else
                    Worksheets("Convert File").Activate
                End If
                Exit Sub
            Else
                FileNo = FreeFile
                Open FName For Output As #FileNo
                With ActiveWorkbook.ActiveSheet
                    For x = 1 To 4
                        Print #FileNo, .Cells(x, 1).Value
                    Next x
                    For x = 5 To solarhourfinal
                        Print #FileNo, .Cells(x, 1).Value & "," & .Cells(x, 2).Value & "," & .Cells(x, 3).Value & "," & .Cells(x, 4).Value & "," & .Cells(x, 5).Value
                    Next x
                End With
                Close #FileNo
                MsgBox "The file was successfully saved!", vbExclamation + vbOKOnly
                ActiveWorkbook.Close SaveChanges:=False
                If Worksheets("Selection Data").Range("Q43") = 1 Then
                    GoTo TotHorRad_Continue
                Else
                    Worksheets("Convert File").Activate
                End If
                Exit Sub
            End If
        Else
            b = MsgBox("Your radiation log file will not be saved!", vbExclamation + vbOKCancel)
            If b = 1 Then
                MsgBox "The file was not saved!", vbExclamation + vbOKOnly
                ActiveWorkbook.Close SaveChanges:=False
                If Worksheets("Selection Data").Range("Q43") = 1 Then
                    GoTo TotHorRad_Continue
                Else
                    Worksheets("Convert File").Activate
                End If
                Exit Sub
            End If
        End If
    Loop While b > 1

End

' On error morphing the file the program is terminated and the data on the morphed weather sheet deleted
Terminate:
Application.Cursor = xlDefault
Application.Goto reference:=Worksheets("Morphed Weather").Range("A1")
Worksheets("Convert File").Activate
MsgBox "An error has occured while morphing your data! Please try again!" & Chr(10) & Chr(10) & "If the problem persists please try to close and reopen the climate change weather file generator.", vbExclamation
Worksheets("Convert File").Range("B50").Value = "    No morphed weather file"
Worksheets("Convert File").Range("B51").Value = " "
Worksheets("Morphed Weather").Range("A1:BA9000").Clear
Exit Sub

End Sub
