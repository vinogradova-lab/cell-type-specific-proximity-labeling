from csv import reader
import json
import questionary 
import pathlib 

def get_input_output_folder(workdir_path, folder_name):

    """ This function identifies if the program was started with provided input and output folders
        i.e this is helpful if program is started without docker - then it will ask for folder paths
        otherwise if it is started in docker we usually provide input and output folder path through the -v parameter.
        In another scenario if you have input and output folder already in the workdir_path prgram will also not ask for folders
        but will ask for confirmation if you are sure that you want to continue and use these.

        Input: workdir_path (as pathlib object/ or not it will get converted to one anyway) and folder name
        Output: folder_path
    """

    workdir_path_obj = pathlib.Path(workdir_path)
    folder_path = workdir_path_obj / folder_name
    query = folder_path.exists() #it exists if started with docker, it does not exists if started without

    if query == True:
        answer = (
            questionary.select(
                "An {} was provided or already exists in your working directory, if either is true confirm with yes to continue:".format(folder_name),
                choices=['yes','no'],
            ).ask()
        )
        if answer == "yes":
            folder_path = folder_path
        else:
            print("Remove the folder called {} from your working directory or please change the name and try again!".format(folder_name))
            exit()
    else:
        folder_path = (questionary.path("Select folder (hit / and direct to your folder, full path required)").ask())
    return folder_path 


def get_channel_name_dict(input_folder_path, list_of_file_names):
    meta_data_list = []
    with open(input_folder_path / "metadata_col.csv", 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            meta_data_list.append(row)

    file_channel_dict = {}
    # file_turbo_id_dict = {}
    for file_name in list_of_file_names:
        channel_dict = {}
        # turboid_dict = {}
        for meta_data_item in meta_data_list[1:]:
            if file_name in meta_data_item:
                channel_dict[meta_data_item[1]] = meta_data_item[2]
                # turboid_dict[meta_data_item[2]] = float(meta_data_item[4])
        file_channel_dict[file_name] = channel_dict
    # file_turbo_id_dict[file_name] = turboid_dict

    # save as dict as json file
    with open("../results/file_channel_details.json", "w") as outfile:
        json.dump(file_channel_dict, outfile)
    
    return file_channel_dict


def get_cond_info(input_folder_path, list_of_file_names):
    condition_list = []
    with open(input_folder_path / "conditions_metadata.csv", 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            condition_list.append(row)

    file_condition_dict = {}
    for file_name in list_of_file_names:
        condition_dict = {}
        for file_row in condition_list[1:]:
            if file_name in file_row:
                condition_dict["conditions"] = file_row[1].replace(" ", "").split(",")
                condition_dict["control_labelling"] = file_row[2]
                condition_dict["treatment_labelling"] = file_row[3]
                condition_dict["file_type"] = file_row[4]
        file_condition_dict[file_name] = condition_dict
    
    # save as dict as json file
    with open("../results/file_condition_details.json", "w") as outfile:
        json.dump(file_condition_dict, outfile)
    
    return file_condition_dict 