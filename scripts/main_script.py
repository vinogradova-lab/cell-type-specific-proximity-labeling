import fire 
import pathlib
import questionary
import os 
from readin_funcs import get_input_output_folder
from turboid_analysis_pipeline import turbo_id_analysis
from create_lists_for_yuval import get_lists_for_yuvals_group

def process_data(): 
    
    workdir_path = pathlib.Path.cwd()
    input_folder_path = get_input_output_folder(workdir_path, "input_folder")
    output_folder_path = get_input_output_folder(workdir_path, "output_folder")
    fasta_table_path = get_input_output_folder(workdir_path, "fasta_table_tp_fp_list")

    input_folder_path = pathlib.Path(input_folder_path)
    input_folder_path = input_folder_path / "01_processed_files"

    output_folder_path = pathlib.Path(output_folder_path)
    output_dir_name = (questionary.text("How do you want to label this run? - make sure to mention if it is 1 pep or 2 pep per protein").ask())
    output_folder_path = output_folder_path / output_dir_name
    os.mkdir(output_folder_path)

    fasta_table_path = pathlib.Path(fasta_table_path)

    turbo_id_analysis(input_folder_path, output_folder_path, fasta_table_path)

    goodbye_message = (
            questionary.select(
                "The program completed succesfully! Thank you for your patience.",
                choices=['Goodbye!']
            ).ask()
        )

    return

def create_lists_for_yuvals_group(): 
    workdir_path = pathlib.Path.cwd()
    output_folder_path = get_input_output_folder(workdir_path, "output_folder")
    fasta_table_path = get_input_output_folder(workdir_path, "fasta_table_tp_fp_list")

    output_folder_path = pathlib.Path(output_folder_path)
    fasta_table_path = pathlib.Path(fasta_table_path)

    get_lists_for_yuvals_group(output_folder_path, fasta_table_path)

    goodbye_message = (
            questionary.select(
                "The program completed succesfully! Thank you for your patience.",
                choices=['Goodbye!']
            ).ask()
        )

if __name__ == "__main__":
    fire.Fire()
