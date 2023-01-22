import argparse

from model import *
from utils import (
    data_loader, save_results
)
import time
import numpy as np

def main(args):
    data_A = data_loader(args['dataset'], args['d_A'])
    data_B = data_loader(args['dataset'], args['d_B'])
    gt_pval = data_loader(args['dataset'], args['GT'])

    seeds = [0, 100, 210, 970, 3506, 23525, 32451, 2466, 1242]

    iters = args['t']
    batch_size = args['b']

    thresholds = {1: args['thre'], 0: 1}

    for i_iter in range(len(iters)):
        ### Starting Phase 1 #########

        step1_path = args['out_dir'] + str(i_iter) + '/step1/'
        step1 = process_1(data_A, data_B)
        p_A, p_B = step1.fit()
        p_A_label, p_B_label, insig_list, column_list = step1.predict(p_A, p_B, thresholds)

        save_results(p_A, step1_path, 'p_A')
        save_results(p_B, step1_path, 'p_B')

        save_results(p_A_label, step1_path, 'p_A_label')
        save_results(column_list, step1_path, 'sig_list')
        save_results(insig_list, step1_path, 'insig_list')

        ###### Starting Phase 2 #######
        step2_path = args['out_dir'] + str(i_iter) + '/step2/'

        nss = [batch_size for _ in range(iters)]

        step2 = process_2(data_A, data_B, insig_list, nss, seeds)
        p_value_uni, selected_idx = step2.fit()
        pval, iter_label, vote_label = step2.predict(p_value_uni, column_list, thresholds)
        all_pred = step2.evaluate_overall(gt_pval, vote_label)

        save_results(p_value_uni, step2_path, 'p_value_uni')
        save_results(selected_idx, step2_path, 'selected_idx')

        save_results(pval, step2_path, 'pval')
        save_results(iter_label, step2_path, 'iter_label')
        save_results(vote_label, step2_path, 'vote_label')
        save_results(all_pred, step2_path, 'all_pred')

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--dataset')
    parser.add_argument('--d_A')
    parser.add_argument('--d_B')
    parser.add_argument('--out_dir',  default='../results')
    parser.add_argument('--GT')
    parser.add_argument('--thre', default = 0.3)
    parser.add_argument('--t', default = 5)
    parser.add_argument('--b', default = 300)

    args_dict = vars(parser.parse_args())
    #args = argparse.Namespace(**args_dict)
    main(args_dict)
