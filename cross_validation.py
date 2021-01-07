"""Runs LOO cross-validation diffusion experiment.

Usage:

    loo_experiment = LOOValidation(network, input_nodes)
    result = loo_experiment.run_validation()
    print(result.head())
    tpr, fpr, auc = loo_experiment.get_roc()

"""

import numpy as np
import pandas as pd
import scipy.stats
import sklearn.metrics

import diffusion


class LOOValitation:
    """Diffusion leave-one-out cross validation experiment."""

    def __init__(self, network, input_nodes):
        """Inits validation experiment with network and input nodes."""
        self.network = network
        self.input_nodes = input_nodes
        if len(self.input_nodes) < 2:
            raise ValueError("need at least 2 input nodes for cross-validation")
        self.result = None

    def run_validation(self):
        """Runs Leave-One-Out cross-validation with the input set."""
        left_out_score = {}
        post_diffusion_scores = []
        for left_out in self.input_nodes:
            input_sans_one = [i for i in self.input_nodes if i != left_out]
            dif = diffusion.Diffusion(self.network, input_sans_one)
            dif_result = dif.diffuse()
            post_diffusion_scores.append(dif_result.final_state)
            left_out_score[left_out] = dif_result.get_result_for_protein(left_out)
        avg_scores = self.average_results(post_diffusion_scores)
        result = self.format_result_as_pd_df(avg_scores, self.network.proteins)
        result = self.insert_leftout_scores(result, left_out_score)
        # preset post-diffusion scores as zscores and ranks
        result["zscore"] = scipy.stats.zscore(result.final_state)
        result["rank"] = result.zscore.rank(ascending=False)
        result.sort_values(by="zscore", ascending=False, inplace=True)
        self.result = result
        return result

    @staticmethod
    def average_results(score_vectors):
        """Averages results of several diffusion experiments."""
        agg = score_vectors[0]
        for result in score_vectors[1:]:
            agg = np.add(agg, result)
        avg = agg / len(score_vectors)
        return avg

    @staticmethod
    def format_result_as_pd_df(avg_scores, proteins):
        """Converts result vector into pandas data frame with protein names column."""
        result = pd.DataFrame(
            {
                "protein": proteins,
                "initial_state": [0] * len(proteins),
                "final_state": avg_scores,
            }
        )
        return result

    @staticmethod
    def insert_leftout_scores(result, left_out_score):
        """Insert left-out scores into final result table."""
        for protein, score in left_out_score.items():
            result.loc[result.protein == protein, "final_state"] = score
            result.loc[result.protein == protein, "initial_state"] = 1
        return result

    def get_roc(self):
        """Calculate ROC."""
        class_labels = self.result.initial_state
        pred_scores = self.result.final_state
        fpr, tpr, _ = sklearn.metrics.roc_curve(
            y_true=class_labels, y_score=pred_scores, pos_label=1
        )
        auc = sklearn.metrics.auc(fpr, tpr)
        return [tpr, fpr, auc]
