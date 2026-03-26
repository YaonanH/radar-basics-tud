echo ""

$port= new-Object System.IO.Ports.SerialPort COM9,115200,None,8,one
$port.open()

$port.writeline("flushCfg")
$port.readline()
$port.readline() > $null

$port.writeline("dfeDataOutputMode 1")
$port.readline()
$port.readline() > $null

$port.writeline("channelCfg 15 5 0")
$port.readline()
$port.readline() > $null

$port.writeline("adcCfg 2 1")
$port.readline()
$port.readline() > $null

$port.writeline("adcbufCfg 0 0 0 1")
$port.readline()
$port.readline() > $null

$port.writeline("profileCfg 0 77 429 7 57.13 0 0 70 1 240 4884 0 0 30")
$port.readline()
$port.readline() > $null

$port.writeline("chirpCfg 0 0 0 0 0 0 0 1")
$port.readline()
$port.readline() > $null

$port.writeline("chirpCfg 1 1 0 0 0 0 0 4")
$port.readline()
$port.readline() > $null

$port.writeline("frameCfg 0 1 16 1 100 1 0")
$port.readline()
$port.readline() > $null

$port.writeline("lowPower 0 0")
$port.readline()
$port.readline() > $null

$port.writeline("sensorStart")
$port.readline()
$port.readline() > $null

echo "... running ..."
Start-Sleep -s 1

$port.writeline("sensorStop")
$port.readline() > $null
echo "CaptureDemo:/>sensorStop"
$port.readline() > $null

echo " (data now available) "

$port.close()
exit
