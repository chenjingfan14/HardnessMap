; AutoIT script to control Emcotest/Struers Durascan 'ecosworkflow' software.
; Confirmed to work through a complete analysis once initial mesh has been completed
; Should work right from the start of a meshed component
; To run from command line: AutoIt3.exe autoWorkflow.au3
;(c) M Roy 2016
#include <File.au3>
#include <ButtonConstants.au3>
#include <GUIConstantsEx.au3>
#include <EditConstants.au3>
#include <WindowsConstants.au3>
#include <WinAPI.au3>
#include <Misc.au3>
#include <GUIListBox.au3>
#include <StaticConstants.au3>
#Region ### START Koda GUI section ### Form=D:\Dropbox\AutoIt\autoWorkflowGUI.kxf
$autoWorkflowGUI = GUICreate("autoWorkflow", 390, 487, 198, 114)
$RunButton = GUICtrlCreateButton("Run", 17, 19, 113, 25)
$ExitButton = GUICtrlCreateButton("Stop", 257, 19, 113, 25)
GUICtrlSetFont(-1, 8, 800, 0, "MS Sans Serif")

$WorkingDirectory = GUICtrlCreateInput("D:\Documents\Hardness\Fix", 18, 76, 353, 21)
$ProgramPath=GUICtrlCreateInput("D:\Documents\Hardness\Fix\Programs", 18, 128, 353, 21)
$OutlinePath = GUICtrlCreateInput("D:\Documents\Hardness\Fix\Results", 16, 176, 353, 21)
$RelocatePath = GUICtrlCreateInput("D:\Documents\Hardness\Fix\Results\refindents.txt", 16, 232, 353, 21)

$Prefix= GUICtrlCreateInput("something", 16, 272, 145, 21)
$Iter = GUICtrlCreateInput("7", 184, 272, 81, 21)
$nRp = GUICtrlCreateInput("120", 288, 272, 81, 21)

$XMLDefault = GUICtrlCreateInput("D:\Specimen", 16, 456, 353, 21)

$Path1 = GUICtrlCreateLabel("Program path", 16, 104, 67, 17)
$Label2 = GUICtrlCreateLabel("User-defined results && outlines path", 16, 152, 175, 17)
GUICtrlCreateLabel("Prefix", 16, 256, 30, 17)
$Iterations = GUICtrlCreateLabel("Iterations", 184, 256, 47, 17)
$Label1 = GUICtrlCreateLabel("F9-Kill F10-Pause", 155, 24, 86, 17)
$Label3 = GUICtrlCreateLabel("XML path (Change under settings > General Settings > File Locations)", 16, 432, 332, 17)
$Label4 = GUICtrlCreateLabel("Working directory", 18, 52, 87, 17)
$Label5 = GUICtrlCreateLabel("Full path to re-orientation points *.txt file", 16, 208, 188, 17)
$Label6 = GUICtrlCreateLabel("Refining points", 288, 256, 74, 17)

GUICtrlSetColor(-1, 0x000000)
GUISetState(@SW_SHOW)
#EndRegion ### END Koda GUI section ###

HotKeySet("{F9}", "Kill")
HotKeySet("{F10}", "Pause")
Global $paused = False
Global $runCount = 0


While 1
	$nMsg = GUIGetMsg()
	Switch $nMsg
		Case $GUI_EVENT_CLOSE
			Exit
		Case $RunButton
			Global $refLoc=GUICtrlRead($RelocatePath)
			Global $maindir=GUICtrlRead($WorkingDirectory)
			Global $NumRefiningPnts=GUICtrlRead($nRp)
			$pp_gui = GUICtrlRead($ProgramPath)
			$op_gui = GUICtrlRead($OutlinePath)
			$pref_gui=GUICtrlRead($Prefix)
			Global $it = GUICtrlRead($Iter)
			Global $XML = GUICtrlRead($XMLdefault)
			
			; Check if the corresponding results file is there
			Global $runCount = 0
			$filenames = _FileListToArray( $op_gui, "*.spe")
			$check = 0
			For $i = 1 to $filenames[0]
			    if StringLeft($filenames[$i], StringLen($pref_gui & " - ")) = $pref_gui & " - " Then
						$check = 1
						$fnameFound=$filenames[$i]
				ConsoleWrite("Found results for run: "& $runCount & "   " & $filenames[$i] & @CRLF)
				EndIf
			Next
			
			
			if $check = 1 Then

				
			Else
				ConsoleWrite("Did not find results for run: "& $runCount & @CRLF)
				RunEcos($pp_gui,$op_gui,$pref_gui)
				
			EndIf
			
			; Now check for other results files

			For $i = 1 to $filenames[0]
				for $j = 1 to $it
				if StringLeft($filenames[$i], StringLen($pref_gui & "_" & String($j))) = $pref_gui & "_" & String($j) Then
					$runCount+=1
					;Need to fix this!!!
					; ConsoleWrite("Compared " & StringLeft($filenames[$i], StringLen($pref_gui & "_" & String($j))) & "   " & $pref_gui & "_" & String($j) & @CRLF)
					ConsoleWrite("Found results for run: "& $runCount & "   " & $filenames[$i] & @CRLF)
				
				EndIf
				Next
			Next
			
			

			
		
			if $it-$runCount > 0 Then

			$runstohappen = $it - $runCount
		
			
			ConsoleWrite("There will be this many further runs: " & $runstohappen & @CRLF)
			
			For $i = 1 to $runstohappen

				
				$spath =_PathMake("", $maindir, "Mcommands", ".m") 
				$spath = StringTrimLeft ( $spath, 1 )
				ConsoleWrite("Running Matlab refinement routine . . ." & @CRLF)

				$hFileOpen = FileOpen($spath, 2)
				FileWrite($hFileOpen, "NumRefiningPnts=" & $NumRefiningPnts & ";" & @CRLF)
				FileWrite($hFileOpen, "Prefix='" & $pref_gui & "';" & @CRLF)
				FileWrite($hFileOpen, "speOut='" & $pp_gui & "\" & $pref_gui & "';" & @CRLF)
				if $runCount = 0 Then
					FileWrite($hFileOpen, "LastRunWorkspace='" & $maindir & "\" & $pref_gui &  "_setup.mat';" & @CRLF)
				Else
				 	FileWrite($hFileOpen, "LastRunWorkspace='" & $maindir & "\" &  $pref_gui & "_"&$runCount & "_setup.mat';" & @CRLF)
				EndIf
				$filenames = _FileListToArray( $op_gui, "*.spe")
				$check = 0
				For $j = 1 to $filenames[0]
					
					if $runCount =0 Then
						if StringLeft($filenames[$j], StringLen($pref_gui & " - ")) = ($pref_gui & " - ") Then
							FileWrite($hFileOpen, "LastRunResults='" & $op_gui & "\" & $filenames[$j] & "';" & @CRLF)
						EndIf
					Else
						if StringLeft($filenames[$j], StringLen($pref_gui & "_" & String($runCount))) = ($pref_gui & "_" & String($runCount)) Then
							FileWrite($hFileOpen, "LastRunResults='" & $op_gui & "\" & $filenames[$j] & "';" & @CRLF)
						EndIf
					EndIf
				Next

				FileWrite($hFileOpen, "RefLoc='" & $refLoc & "';" & @CRLF)
				FileClose($hFileOpen)
					
				;then run the following command
				$cmd= "cd /d " & $maindir & " & matlab -wait -minimize -nosplash -nodesktop -log file matcmdLog.log -r ""run('MapRefineAuto.m');exit;"" "
					
				RunWait(@ComSpec & " /c " & $cmd )
					
				RunEcos($pp_gui,$op_gui,$pref_gui & "_" & $runCount +1)
						
				$runCount = $runCount+1				

			
			Next

			EndIf ; it is greater than 0
			
			

			
			ConsoleWrite("Iterations satisfied; program complete. " & @CRLF & "***********************************" & @CRLF)

			
		Case $ExitButton
			Exit

	EndSwitch
WEnd



Func RunEcos($pp,$op,$pref)
	Global $hWnd = WinGetHandle(" ecos Workflow", "")
	$ecosRunning = WinExists($hWnd, "")
	if $ecosRunning = 1 Then

		$spath =_PathMake("", $pp, $pref, ".spe")
		$spath = StringTrimLeft ( $spath, 1 )
		$check = _FilePathIsValid($spath)
		; ConsoleWrite("pathvalid:" & $spath & $check & @CRLF)
		if $check = 1 Then
			WinSetState($hWnd, "", @SW_RESTORE )
			$success = ControlClick ($hWnd, "", "[NAME:btnLoadPattern]", "Left")
			ConsoleWrite("Started ecos Workflow & clicked first button (1 successful) " & $success & @CRLF)
			WinWait("Open","",1500)

			$oWnd = WinGetHandle("Open", "")
			ControlSend($oWnd, "", "[CLASS:Edit; INSTANCE:1]", $spath & "{ENTER}")
			Sleep(1000)
			Send("{ENTER}")

			Global $CurrentRunResultsFileName = ControlGetText ( $hWnd, "", "[NAME:tbSpecimenName]" )
			$success = ControlClick ($hWnd, "", "[NAME:pbNextTPSpecimen]", "Left")
			ConsoleWrite("SpecimenLoad: " & $success & @CRLF)
			; ConsoleWrite("CurrentResultsFname: " & $CurrentRunResultsFileName & @CRLF)
			Sleep(1000)
			$success = ControlClick ($hWnd, "", "[NAME:pbNextTPMethod]", "Left")
			ConsoleWrite("method: " & $success & @CRLF)
			Sleep(1000)
			$success = ControlClick ($hWnd, "", "[NAME:pbStartMeasurementTPPosition]", "Left")
			ConsoleWrite("RUNNING . . . " & @CRLF)
			
			;; Do a winwait for the messagebox identifying that the test has completed ;;
			WinWait("Info","","")
			$success = ControlClick ("Info", "","Button1", "Left")
			ConsoleWrite("info: " & $success & @CRLF)
			;; Check what the iteration count is ;;
			;; Push to history ;;
			Sleep(1000)
			WinSetState($hWnd, "", @SW_RESTORE )
			; WindowsForms10.Window.8.app.0.f96fc5_r14_ad178
			$success = ControlClick ($hWnd, "", "[NAME:pbMoveSpecToHist]", "Left")
			ConsoleWrite("HistoryMove: " & $success & @CRLF)
			
			;; Copy the results .spe file to the directory specified
			Local $rpath =_PathMake("", $XML, $CurrentRunResultsFileName, ".spe")
			Local $dpath =_PathMake("", $op, $CurrentRunResultsFileName, ".spe")
			$success = FileCopy (StringTrimLeft ( $rpath, 1 ), StringTrimLeft ( $dpath, 1 ))
			ConsoleWrite("FINISHED. Copying: "  & StringTrimLeft ( $rpath, 1 ) & @CRLF & "to" & @CRLF & StringTrimLeft ( $dpath, 1 ) & @CRLF)

			
			;; Get the interface ready for the next run by tabbing back to 'Specimen'
			Do
				ControlCommand ($hWnd, "", "WindowsForms10.SysTabControl32.app.0.f96fc5_r14_ad11", "TabLeft", "")
				$TabNo = ControlCommand ($hWnd, "", "WindowsForms10.SysTabControl32.app.0.f96fc5_r14_ad11", "CurrentTab", "")
			Until $TabNo = 1
			
		EndIf
	Else
		Tooltip('autoWorkflow cannot find an instance of ecos Workflow . . .',0,0)
	EndIf ;==>ecosRunning

EndFunc ;==> RunEcos

Func RunEcosTest($pp,$op,$pref)
	; Test Function to replace RunEcos. Pauses/resumes with escape once results are copied over
	ConsoleWrite("RunEcosTest on " & $pp & "  " & $op & "  " & $pref & @CRLF)
	Local $hDLL = DllOpen("user32.dll")
	Do
		Sleep(1)
	Until _IsPressed("1B", $hDLL)
EndFunc


Func Kill()
	; with fire
	Exit 0
EndFunc

Func Pause()
	; pause and tell the user the script is paused
	$paused = NOT $paused
	While $paused
		sleep(100)
		Tooltip('autoWorkflow paused . . .',0,0)
	Wend
	Tooltip("")
EndFunc

Func _FilePathIsValid($sPath)
    Local $sInvalid_Chars = '*?|:<>"/'
    Local $sPattern = '(?i)^[a-z]:([^\Q' & $sInvalid_Chars & '\E]+)?(\\[^\Q' & $sInvalid_Chars & '\E\\]+)?(\.[^\Q' & $sInvalid_Chars & '\E\\\.]+|\\)?$'
    Local $iPathIsValid = StringRegExp($sPath, $sPattern)
    
    Return $iPathIsValid
EndFunc

Func dirUp($s_path, $i_times = 1)
    If $i_times = 0 Then Return $s_path
    $s_path &= '\..'
    For $i = 1 To $i_times
        $s_path &= '\..'
    Next
    Return $s_path
EndFunc  ;==>dirUp