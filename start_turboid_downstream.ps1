param($1)
if ($1 -eq "process_data") {
    Write-Host ""
    Write-Host "Welcome to the Vinogradova-Lab TurboID analysis downsteam Pipeline!"
    Write-Host "You are about to process your files!"
    Write-Host "You will be asked for folder paths, please provide the full folder path!"
    Write-Host ""
    $INPUTFOLDER = Read-Host "Enter your input folder path"
    Write-Host ""
    $OUTPUTFOLDER = Read-Host "Enter your output folder path"
    Write-Host ""
    $FASTATABLE = Read-Host "Enter your fasta table file path (file name: main_fasta_table_without_signal_p.csv):"
    Write-Host ""
    docker run --platform linux/amd64 -m 6g -t -i -v ${INPUTFOLDER}:/work_dir/input_folder -v ${OUTPUTFOLDER}:/work_dir/output_folder -v ${FASTATABLE}:/work_dir/fasta_table_tp_fp_list turboid_downstream $1
}
elseif ($1 -eq "create_lists_for_yuvals_group") {
    Write-Host ""
    Write-Host "Welcome to the Vinogradova-Lab script which creates lists for Yuval's group!"
    Write-Host "You are about to create your reactivity changes file!"
    Write-Host "You will be asked for folder paths, please provide the full folder path!"
    Write-Host ""
    $OUTPUTFOLDER = Read-Host "Enter your output folder path"
    Write-Host ""
    $FASTATABLE = Read-Host "Enter your fasta table file path (file name: main_fasta_table_without_signal_p.csv):"
    Write-Host ""
    docker run --platform linux/amd64 -m 6g -t -i -v ${OUTPUTFOLDER}:/work_dir/output_folder -v ${FASTATABLE}:/work_dir/fasta_table_tp_fp_list turboid_downstream $1
}
else {
    Write-Host "Please provide the correct program name"
}