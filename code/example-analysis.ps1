# Define the responses in the order the Python script expects them.
$responses = @(
    "example-analysis",
    "y",
    "DBr",
    "y",
    "384",
    "y",
    "79.904",
    "2.01410",
    "y",
    "5218fa",
    "y",
    "y",
    "DCl",
    "y",
    "480",
    "y",
    "35.453",
    "2.01410",
    "y",
    "1E90FF",
    "y",
    "y",
    "HBr",
    "y",
    "382",
    "y",
    "79.904",
    "1.00784",
    "y",
    "b700b1",
    "y",
    "y",
    "HCl",
    "y",
    "478",
    "y",
    "35.453",
    "1.00784",
    "y",
    "2e6577",
    "y",
    "n",
    "1685.03",
    "1831.078",
    "y",
    "0.30897",
    "8",
    "y",
    "1847.9",
    "1951.383",
    "y",
    "0.29508",
    "5",
    "y",
    "1960.049",
    "2079.589",
    "n",
    "1959.539",
    "2080.099",
    "y",
    "0.33888",
    "10",
    "y",
    "2100.234",
    "2229.46",
    "y",
    "0.2669",
    "6",
    "y",
    "2371.175",
    "2541.946",
    "y",
    "0.28632",
    "6",
    "y",
    "2572.787",
    "2698.444",
    "y",
    "0.29537",
    "6",
    "y",
    "2726.481",
    "2864.883",
    "y",
    "0.28859",
    "20",
    "y",
    "2904.135",
    "3058.339",
    "y",
    "0.26541",
    "16",
    "y",

)

# Join responses with newline characters.
$inputData = $responses -join "`n"

# Start the Python process and send the responses as standard input.
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = "python"
$psi.Arguments = "analyze.py"
$psi.RedirectStandardInput = $true
$psi.RedirectStandardOutput = $true
$psi.UseShellExecute = $false
$psi.CreateNoWindow = $true

$proc = New-Object System.Diagnostics.Process
$proc.StartInfo = $psi
$proc.Start() | Out-Null

# Write the prepared input data to the Python process.
$proc.StandardInput.WriteLine($inputData)
$proc.StandardInput.Close()

# Read the output from the Python process.
$output = $proc.StandardOutput.ReadToEnd()
$proc.WaitForExit()

# Display the Python output.
Write-Output $output
