#include <File.au3>
#include <FileConstants.au3>
#include <Array.au3>

$array = _FileListToArray(@ScriptDir, "*.csv")
$outFile = FileOpen("results_combined.csv", $FO_OVERWRITE)
FileWriteLine($outFile, "Strat,F(x),Average,Std. Deviation,Range,Median,Min")

For $i = 1 to $array[0]
  processFile($array[$i], $outFile)
Next

FileClose($outFile)

MsgBox(0, "Done", "Done")
Exit

Func processFile($fileName, $fileHnd)
   $fileNameSplit = StringSplit($fileName, '-')

   if $fileNameSplit[1] <> "results" Then
	  Return
   EndIf

   $strat = $fileNameSplit[2]
   $f = $fileNameSplit[3]

   FileWrite($fileHnd, $strat & "," & $f & ",")

   $curFile = FileOpen($fileName, $FO_READ)
   $lastLine = FileReadLine($curFile, -1)
   $lastLineSplit = StringSplit($lastLine, ",")

   $avg = StringFormat("%f", calcAverage($lastLineSplit))
   $stdDev = StringFormat("%f", calcStdDev($lastLineSplit))
   $range = StringFormat("%f", calcRange($lastLineSplit))
   $median = StringFormat("%f", calcMedian($lastLineSplit))
   $min = StringFormat("%f", calcMin($lastLineSplit))

   FileWrite($fileHnd, $avg & ",")
   FileWrite($fileHnd, $stdDev & ",")
   FileWrite($fileHnd, $range & ",")
   FileWrite($fileHnd, $median & ",")
   FileWriteLine($fileHnd, $min)

EndFunc

Func calcAverage($arr)
   $sum = 0.0

   For $i = 1 to $arr[0]
	  $sum += Number($arr[$i])
   Next

   Return  $sum / $arr[0]
EndFunc

Func calcStdDev($arr)
   $avg = calcAverage($arr)
   $tmp = $arr

   For $i = 1 to $tmp[0]
	  $tmp[$i] =  Number($tmp[$i]) - $avg
	  $tmp[$i] = $tmp[$i] * $tmp[$i]
   Next

   $avg2 = calcAverage($tmp)

   Return Sqrt($avg2)
EndFunc

Func calcRange($arr)
   $min = Number($arr[1])
   $max = Number($arr[1])

   For $i = 2 to $arr[0]
	  $cur = Number($arr[$i])

	  If $cur < $min Then
		 $min = $cur
	  EndIf

	  If $cur > $max Then
		 $max = $cur
	  EndIf
   Next

   Return $max - $min
EndFunc

Func calcMedian($arr)
   $tmp = $arr

   _ArrayDelete($tmp, 0)
   _ArraySort($tmp)

   If Mod(UBound($tmp), 2) <> 0 Then
	  $index = Int(UBound($tmp) / 2)
	  Return $tmp[$index]
   Else
	  $lower = Int(UBound($tmp) / 2)
	  $upper = $lower + 1
	  Return ($tmp[$lower] + $tmp[$upper]) / 2
   EndIf
EndFunc

Func calcMin($arr)
   $min = Number($arr[1])

   For $i = 2 to $arr[0]
	  $cur = Number($arr[$i])

	  If $cur < $min Then
		 $min = $cur
	  EndIf
   Next

   Return $min
EndFunc
