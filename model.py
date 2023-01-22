import pandas as pd
import numpy as np
import statsmodels.api as sm
import warnings
from sklearn.metrics import confusion_matrix
warnings.filterwarnings('ignore')



def get_pvalues(data):
    p_value=[]
    for i in range(data.shape[1]-2):
        df = pd.DataFrame({'x': data.iloc[:,i], 'y': data['y']})
        df['intecept']=1
        try:
            logit_model=sm.Logit(df['y'],df[['x','intecept']],missing = 'drop')
            result=logit_model.fit(disp=0)
            p_value.append(result.pvalues[0])
        except:
            p_value.append(2)
    return np.array(p_value)


def obtain_label(data, thresholds):
    """
    Given the p-values for each iteration and the thresholds
    Return the class label for each iteration
    """
    data_new = pd.DataFrame(columns=data.columns)
    data_new['SNP'] = data['SNP']
    for col in list(data_new.columns)[1:]:
        data_new[col] = pd.cut(data[col], bins=[-0.000001] + list(thresholds.values())[:-1] + [float('inf')],
                               labels=list(thresholds.keys()))
    return data_new



class process_1():
    def __init__(self, data_A, data_B):
        self.data_A = data_A
        self.data_B = data_B

    def fit(self):
        p_A = pd.DataFrame({'SNP': data_A.columns[:-1]})
        p_B = pd.DataFrame({'SNP': data_B.columns[:-1]})

        p_A['p_value'] = get_pvalues(data_A)
        p_B['p_value'] = get_pvalues(data_B)

        return p_A, p_B

    def predict(self, p_A, p_B, thresholds):
        p_A_label = obtain_label(p_A, thresholds)
        p_B_label = obtain_label(p_B, thresholds)

        insig_A = set(p_A_label[p_A_label['p_value'] == 0]['SNP'])
        insig_B = set(p_B_label[p_B_label['p_value'] == 0]['SNP'])

        insig_list = list(insig_A.intersection(insig_B))
        sig_list = list(p_A_label[~p_A_label.set_index('SNP').index.isin(insig_list)]['SNP'])

        p_A_label['predict_label'] = [1 for _ in range(len(p_A_label))]
        p_A_label['predict_label'].mask(p_A_label.set_index('SNP').index.isin(insig_list), 0, inplace=True)

        p_B_label['predict_label'] = [1 for _ in range(len(p_B_label))]
        p_B_label['predict_label'].mask(p_B_label.set_index('SNP').index.isin(insig_list), 0, inplace=True)

        return p_A_label, p_B_label, insig_list, sig_list

    def evaluate(self, gt_pval, p_A_label):
        true_label = gt_pval['Label']
        predict_label = p_A_label['predict_label']

        tn, fp, fn, tp = confusion_matrix(true_label, predict_label).ravel()

        return fn, tn


class process_2():

    def __init__(self, data_A, data_B, insig_list, nss, seeds):

        self.data_A = data_A
        self.data_B = data_B
        self.insig_list = insig_list
        self.nss = nss
        self.seeds = seeds

    def iterative_uni(self, data1, data2):
        """
        uniformly sampling records and get inference
        nss: iteratiton rounds
        seed: random state
        """
        p_value = np.zeros((data1.shape[1] - 1, len(self.nss)))

        selected_uni = []

        for idn, n in enumerate(self.nss):
            df_uni1 = data1.sample(n, random_state=self.seeds[idn])
            df_uni2 = data2.sample(n, random_state=self.seeds[idn])

            df_uni = pd.concat([df_uni1, df_uni2], axis=0)
            df_uni = df_uni.sample(len(df_uni))

            p_value[:, idn] = get_pvalues(df_uni)

            selected_uni.append(list(df_uni1.index))

        return p_value, selected_uni

    def ensembleVote(self, data):
        """
        Given the class label for each iteration
        Return the class label via majority voting
        """
        data_new = pd.DataFrame(columns=data.columns)
        data_new['SNP'] = data['SNP']
        # data_new['p_value'] = data['p_value']
        # data_new['group'] = data['group']

        for idn in range(len(self.nss)):
            if idn == 0:
                data_new[str(idn + 1)] = data[str(idn + 1)]
            else:
                iters = [str(i + 1) for i in range(idn + 1)]
                u = data[iters].mode(axis=1).iloc[:, 0]
                data_new[str(idn + 1)] = u

        return data_new

    def fit(self):

        """
        uniformly sampling records and get inference
        nss: list of subset size
        seed: random state
        """

        # filter out the insig snps from step 1
        data_A_sub = self.data_A.loc[:, ~self.data_A.columns.isin(self.insig_list)]
        data_B_sub = self.data_B.loc[:, ~self.data_B.columns.isin(self.insig_list)]

        return self.iterative_uni(data_A_sub, data_B_sub)

    def predict(self, p_value_uni, column_list, thresholds):
        # get label for each iteration

        p_values = pd.DataFrame({'SNP': column_list})

        for idn, n in enumerate(nss):
            p_uni_sample = dict(zip(column_list, p_value_uni[:, idn]))

            p_values[str(idn + 1)] = p_values['SNP'].map(p_uni_sample)

        iter_label = obtain_label(p_values, thresholds)

        # get the label at each iteration via voting
        vote_label = self.ensembleVote(iter_label)
        return p_values, iter_label, vote_label

    def evaluate_process2(self, gt_pval, vote_label):
        true_label = gt_pval[~gt_pval.set_index('SNP').index.isin(self.insig_list)]['Label']
        tp_ = []
        fn_ = []

        for idn in range(len(nss)):
            predict_label = vote_label[str(idn + 1)]
            tn, fp, fn, tp = confusion_matrix(true_label, predict_label).ravel()
            tp_.append(tp)
            fn_.append(fn)

        return tp_, fn_

    def evaluate_overall(self, gt_pval, vote_label):
        true_label = gt_pval['Label']

        predict_label = pd.DataFrame({'SNP': self.data_A.columns[:-1]})

        for idn in range(len(nss)):
            label_iter = np.zeros((len(self.data_A.columns[:-1])))
            idx_list = [list(self.data_A.columns[:-1]).index(snp) for snp in self.insig_list]

            mask = np.ones(len(label_iter), dtype=bool)
            mask[idx_list] = False

            label_iter[mask] = vote_label[str(idn + 1)]

            predict_label[str(idn + 1)] = label_iter

        return predict_label