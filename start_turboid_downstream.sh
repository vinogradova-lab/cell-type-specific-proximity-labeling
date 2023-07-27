#!/bin/zsh
if [ $1 = "process_data" ]
then
  echo ""
  echo "Welcome to the Vinogradova-Lab TurboID analysis downsteam Pipeline!"
  echo "You are about to process your files!"
  echo "You will be asked for folder paths, please provide the full folder path!"
  echo ""
  echo "Enter your input folder path:"
  read INPUTFOLDER
  echo ""
  echo "Enter your output folder path:"
  read OUTPUTFOLDER
  docker run -t -i -v $INPUTFOLDER:/work_dir/input_folder -v $OUTPUTFOLDER:/work_dir/output_folder turboid_downstream $1
elif [ $1 = "create_fasta_file" ]
then
  echo ""
  echo "Welcome to the Vinogradova-Lab script which creates the main fasta table for the turboid project!"
  echo "You are about to create the main fasta file table for the turboid project!"
  echo "You will be asked for folder paths, please provide the full folder path!"
  echo ""
  echo "Enter your input folder path which contains both reactivity and whole proteome files:"
  read INPUTFOLDER
  echo ""
  echo "Enter your output folder path:"
  read OUTPUTFOLDER
  docker run -t -i -v $INPUTFOLDER:/work_dir/input_folder -v $OUTPUTFOLDER:/work_dir/output_folder turboid_downstream $1
else
  echo "Please provide the correct program name"
fi
