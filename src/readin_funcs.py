from csv import reader
import json

def get_channel_name_dict(input_folder_path, list_of_file_names):
    meta_data_list = []
    with open(input_folder_path / "metadata_col.csv", 'r') as read_obj:
        csv_reader = reader(read_obj)
        for row in csv_reader:
            meta_data_list.append(row)

    file_channel_dict = {}
    #file_turbo_id_dict = {}
    for file_name in list_of_file_names:
        channel_dict = {}
        #turboid_dict = {}
        for meta_data_item in meta_data_list[1:]:
            if file_name in meta_data_item:
                channel_dict[meta_data_item[1]] = meta_data_item[2]
                #turboid_dict[meta_data_item[2]] = float(meta_data_item[4])
        file_channel_dict[file_name] = channel_dict
    #file_turbo_id_dict[file_name] = turboid_dict

    #save as dict as json file
    with open("file_channel_details.json", "w") as outfile:
        json.dump(file_channel_dict, outfile)
    
    return(file_channel_dict)

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
        file_condition_dict[file_name] = condition_dict
    
    #save as dict as json file
    with open("file_condition_details.json", "w") as outfile:
        json.dump(file_condition_dict, outfile)
    
    return(file_condition_dict)