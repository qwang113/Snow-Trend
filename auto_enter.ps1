$INTERVAL = 10

Add-Type -AssemblyName System.Windows.Forms

while ($true) {

    $proc = Get-Process | Where-Object {
        $_.ProcessName -eq "cmd" -and $_.MainWindowTitle -like "*weekly*"
    }

    if ($proc) {
        $wshell = New-Object -ComObject WScript.Shell
        $wshell.AppActivate($proc[0].Id)

        Start-Sleep -Milliseconds 200
        [System.Windows.Forms.SendKeys]::SendWait("{ENTER}")

        Write-Output "Enter sent"
    } else {
        Write-Output "No matching window"
    }

    Start-Sleep -Seconds $INTERVAL
}