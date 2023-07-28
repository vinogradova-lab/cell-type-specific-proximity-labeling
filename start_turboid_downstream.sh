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
  echo ""
  echo "Enter the folder path to the fasta table file (please provide folder where the file is located)"
  read FASTATABLE
  echo ""
  docker run -t -i -v $INPUTFOLDER:/work_dir/input_folder -v $OUTPUTFOLDER:/work_dir/output_folder -v $FASTATABLE:/work_dir/fasta_table_tp_fp_list turboid_downstream $1
elif [ $1 = "create_lists_for_yuvals_group" ]
then
  echo ""
  echo "Welcome to the Vinogradova-Lab script which creates lists for Yuval's group!"
  echo "You will be asked for folder paths, please provide the full folder path!"
  echo ""
  echo "Enter your output folder path:"
  read OUTPUTFOLDER
  echo ""
  echo "Enter the folder path to the fasta table file (please provide folder where the file is located)"
  read FASTATABLE
  echo ""
  docker run -t -i -v $OUTPUTFOLDER:/work_dir/output_folder -v $FASTATABLE:/work_dir/fasta_table_tp_fp_list turboid_downstream $1
else
  echo "Please provide the correct program name"
fi