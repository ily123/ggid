import numpy as np
import pandas as pd
import sklearn.metrics

import diffusion

# from sklearn.metrics import roc_auc_score


class LOOValitation:
    """
    Leave-one-out cross validation of input labels.
    """

    def __init__(self, network, labels):
        """
        Inits validation experiment.
        """
        self.network = network
        self.labels = labels
        if len(self.labels) < 2:
            raise ("need at least 2 labels for cross-validation")

    def run_validation(self):
        """
        Runs Leave-One-Out cross-validation.
        """

        label_score = {}
        score_vectors = []
        merged = None
        for left_out_label in self.labels:
            input_labels = [i for i in self.labels if i != left_out_label]
            dif = diffusion.Diffusion(self.network, input_labels)
            result = dif.diffuse()
            score_vectors.append(result.final_labels)
            label_score[left_out_label] = result.get_result_for_protein(left_out_label)
            print(left_out_label, result.get_result_for_protein(left_out_label))
        avg_result = self.average_results(score_vectors)
        avg_result = self.format_result_as_pd_df(avg_result, self.network.proteins)
        avg_result = self.replace_labels(avg_result, label_score)

        self.avg_result = avg_result
        return avg_result

    @staticmethod
    def average_results(score_vectors):
        """
        Averages results of several diffusion experiments.
        """

        agg = score_vectors[0]
        for result in score_vectors[1:]:
            agg = np.add(agg, result)
        avg = agg / len(score_vectors)
        return avg

    @staticmethod
    def format_result_as_pd_df(avg_result, proteins):
        """
        Converts result vector into pandas data frame.
        """

        return pd.DataFrame({"protein_name": proteins, "avg_label": avg_result})

    @staticmethod
    def replace_labels(avg_result, label_score):
        """
        Insert left-out label scores into final label result.
        """
        print(label_score)
        avg_result["label"] = 0
        for protein, score in label_score.items():
            avg_result.loc[avg_result.protein_name == protein, "avg_label"] = score
            avg_result.loc[avg_result.protein_name == protein, "label"] = 1
        return avg_result

    def get_roc(self):
        # self.avg_result
        labels = self.avg_result.label
        pred_scores = self.avg_result.avg_label
        fpr, tpr, _ = sklearn.metrics.roc_curve(
            y_true=labels, y_score=pred_scores, pos_label=1
        )
        auc = sklearn.metrics.auc(fpr, tpr)
        auc2 = sklearn.metrics.roc_auc_score(y_true=labels, y_score=pred_scores)
        print(auc)
        print(auc2)
        return [tpr, fpr, auc]

    # @staticmethod
    # def get_auc():
