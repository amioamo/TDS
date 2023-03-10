import pandas as pd
import os
import pickle
def data_loader(dataset, data_path):
    if dataset == 'human1000':
        gen = pd.read_csv(data_path, index_col='IID')

    elif dataset == 'opensnp':
        gen = pd.read_csv(data_path)
    else:
        print('not default dataset')

    return gen


def partition(dataset, data, scenario, seed):
    case_ = data['y'] == 1
    control_ = data['y'] == 0

    if dataset == 'human1000':
        if scenario == '3' or scenario == '5':
            ## when both parties have imbalanced datasets (case:control= 300:900, vs 900:300)
            data_A_case = data[case_].sample(300, random_state=seed)
            data_A_control = data[control_].sample(900, random_state=seed)

            data_B_case = data[case_].loc[~ data[case_].index.isin(data_A_case.index)]
            data_B_control = data[control_].loc[~ data[control_].index.isin(data_A_control.index)]
        elif scenario == '4':
            data_A_case = data[case_].sample(450, random_state=seed)
            data_A_control = data[control_].sample(750, random_state=seed)

            data_B_case = data[case_].loc[~ data[case_].index.isin(data_A_case.index)]
            data_B_control = data[control_].loc[~ data[control_].index.isin(data_A_control.index)]
        else:
            ## when both are balanced or when party A (case:control= 600:600) and party B (case:control=300:900)
            data_A_case = data[case_].sample(600, random_state=seed)
            data_A_control = data[control_].sample(600, random_state=seed)

            data_B_case = data[case_].loc[~ data[case_].index.isin(data_A_case.index)]
            data_B_control = data[control_].loc[~ data[control_].index.isin(data_A_control.index)]

        data_A = pd.concat([data_A_case, data_A_control]).sample(len(data_A_case) + len(data_A_control))
        data_B = pd.concat([data_B_case, data_B_control]).sample(len(data_B_case) + len(data_B_control))
    elif dataset == 'opensnp':
        data_A_case = data[case_].sample(200, random_state=seed)
        data_A_control = data[control_].sample(200, random_state=seed)

        data_B_case = data[case_].loc[~ data[case_].index.isin(data_A_case.index)]
        data_B_control = data[control_].loc[~ data[control_].index.isin(data_A_control.index)]

        data_A = pd.concat([data_A_case, data_A_control]).sample(len(data_A_case) + len(data_A_control))
        data_B = pd.concat([data_B_case, data_B_control]).sample(len(data_B_case) + len(data_B_control))

    else:
        print('not default dataset')

    return data_A, data_B


def save_results(res, data_path, file_name):
    if not os.path.isdir(data_path):
        os.makedirs(data_path)

    if isinstance(res, pd.DataFrame):
        res.to_csv(data_path + file_name + '.csv', index=False)
    else:
        with open(data_path + file_name + '.pkl', 'wb') as f:
            pickle.dump(res, f)
