#include <File.au3>
#include <FileConstants.au3>

$array = _FileListToArray(@ScriptDir, "*.csv")
$outFile = FileOpen("results_combined.csv", $FO_OVERWRITE)
FileWriteLine($outFile, "f/strat,Average,Std. Deviation,Range,Median,Min")

For $i = 1 to $array[0]
  processFile($array[$i], $outFile)
Next

FileClose($outFile)

Func processFile($fileName, $fileHnd)
   $fileNameSplit = StringSplit($fileName, '-')

   if $fileNameSplit[1] <> "results" Then
	  Return
   EndIf

   $alg = $fileNameSplit[2]

   if $alg <> "DE" Then
	  Return
   EndIf

   $strat = $fileNameSplit[3]
   $f = $fileNameSplit[4]

   FileWrite($fileHnd, $strat & " " & $f & ",")

   $curFile = FileOpen($fileName, $FO_READ)
   $lastLine = FileReadLine($curFile, -1)
   $lastLineSplit = StringSplit($lastLine, ",")

   $avg = calcAverage($lastLineSplit)

   FileWriteLine($fileHnd, $avg)

EndFunc

Func calcAverage($arr)
   $sum = 0

   For $i = 1 to $arr[0]
	  $sum += Number($arr[$i])
   Next

   Return StringFormat("%f", $sum / $arr[0])
EndFunc