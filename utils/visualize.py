import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as tick
from matplotlib.patches import Rectangle
import seaborn as sns
import warnings

from config import exp_cfg, formal_method_names
from typing import Union, List
import pickle
import os


class PerformanceVisualizer:
    def __init__(self, task: str):
        assert task in ('assign', 'cluster'), ValueError(f"{task} is an invalid argument.")
        self.record_path = exp_cfg.record_path
        self.task = task
        self.record_names = [file for file in os.listdir(self.record_path)
                             if file.endswith('pkl') and file.split('-')[1] == task]
        self.task_records = self.load_record(self.record_names)  # {dataset-n_gene: record}
        self.save_dir = exp_cfg.figure_path

        self.datasets = np.unique([file[:-4].split('-')[0] for file in self.record_names])
        self.n_genes = sorted(np.unique([file[:-4].split('-')[2] for file in self.record_names]), key=lambda x: int(x))

        if task == 'assign':
            self.task_methods = self.task_records[self.datasets[0] + '-' + self.n_genes[0]].singleR_F1.columns.tolist()
        else:
            self.task_methods = self.task_records[self.datasets[0] + '-' + self.n_genes[0]].seurat_ARI.columns.tolist()

    def load_record(self, file_name: Union[str, List[str]]):
        if isinstance(file_name, str):
            with open(self.record_path + file_name, 'rb') as f:
                return pickle.load(f)
        else:
            if len(file_name) > 0:
                return {name[:-4].split('-')[0] + '-' + name[:-4].split('-')[2]: self.load_record(name) for name in
                        file_name}
            else:
                warnings.warn("file_name doesn't have an element. Return null dict.")
                return dict()

    def plot_lines(self, metric: str = 'markers_found'):
        fig, ax = plt.subplots(figsize=(9 if self.task == 'assign' else 7, 5))
        drop_methods = ['scGeneFit', 'cv2'] if self.task == 'assign' else ['cv2']
        spec_gene_table = pd.DataFrame(np.zeros(shape=(len(self.n_genes), len(self.task_methods))),
                                       index=self.n_genes, columns=self.task_methods).drop(columns=drop_methods)
        for dataset in self.datasets:
            for n_gene in self.n_genes:
                record = self.task_records[dataset + '-' + n_gene]
                if metric == 'markers_found':
                    spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc[metric, :]
                else:
                    if self.task == 'assign':
                        spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc['singlecellnet_F1', :]
                        spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc['singleR_F1', :]
                        spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc['itclust_F1', :]
                    else:
                        spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc['seurat_ARI', :]
                        spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc['sc3_ARI', :]
                        spec_gene_table.loc[n_gene, :] += record.summary().drop(columns=drop_methods).loc['cidr_ARI', :]

        if metric == 'markers_found':
            spec_gene_table = spec_gene_table.div(len(self.datasets) * pd.Series(self.n_genes, index=spec_gene_table.index).astype(np.float_), axis=0)
            ylabel, title = 'The percentage of marker genes', None
        else:
            spec_gene_table = spec_gene_table.div(len(self.datasets) * 3)
            ylabel = 'The average F1-score' if self.task == 'assign' else 'The average ARI'
        spec_gene_table.rename(columns=formal_method_names, inplace=True)
        spec_gene_table.plot(marker='.', linestyle='--', linewidth=1.5, cmap='tab20', ax=ax,
                             xlabel='Number of selected genes', ylabel=ylabel)
        ax.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=6 if self.task == 'assign' else 4, frameon=False)
        if metric == 'markers_found':
            ax.yaxis.set_major_formatter(tick.PercentFormatter(xmax=1, decimals=0))
        fig.set_tight_layout(True)
        file_name = '-'.join([self.task, 'line']) + '.jpg'
        plt.savefig(self.save_dir + file_name, dpi=120, bbox_inches='tight')

    def plot_metric_heatmap(self, metric='MRR'):
        if self.task == 'assign':
            fig, axes = plt.subplots(2, 7, sharey='all', figsize=(10, 6))
        else:
            fig, axes = plt.subplots(1, 7, sharey='all', figsize=(10, 4))
            axes = axes[np.newaxis, :]
        fig.subplots_adjust(wspace=0.05)
        for i, method in enumerate(self.task_methods):
            spec_method_table = pd.DataFrame(np.zeros((len(self.datasets), len(self.n_genes))),
                                             index=self.datasets, columns=self.n_genes)

            for dataset in self.datasets:
                for n_gene in self.n_genes:
                    record = self.task_records[dataset + '-' + n_gene]
                    ratio_factor = record.n_marker_contain if metric == 'markers_found' else 1.
                    if hasattr(record, metric):
                        spec_method_table.loc[dataset, n_gene] = record.summary().loc[metric, method] / ratio_factor
                    else:
                        raise ValueError(f"record of {dataset + '-' + n_gene} doesn't have metric {metric}.")
            # start to plot
            if i != len(self.task_methods) - 1:
                sns.heatmap(spec_method_table, annot=True, fmt='.3f' if metric != 'markers_found' else '.1%',
                            ax=axes[i // 7, i % 7], cmap="YlGnBu", cbar=False, vmin=0,
                            vmax=1 if metric != 'MRR' else 0.1)
                axes[i // 7, i % 7].text(1, 4.5, formal_method_names[method], ha='center', va='center')
            else:
                sns.heatmap(spec_method_table.iloc[:, 0].to_frame(), annot=True,
                            fmt='.3f' if metric != 'markers_found' else '.1%',
                            ax=axes[i // 7, i % 7], cmap="YlGnBu", cbar=True, vmin=0,
                            vmax=1 if metric != 'MRR' else 0.1,
                            cbar_kws={
                                'format': tick.FuncFormatter(
                                    lambda x, y: '{:.0%}'.format(x)) if metric == 'markers_found' else None,
                                'shrink': 0.8,
                                'fraction': 0.3,
                                'anchor': (0.0, 0.0),
                                'panchor': (1.0, 1.0),
                                'use_gridspec': False
                            })
                axes[i // 7, i % 7].text(0.5, 4.5, formal_method_names[method], ha='center',
                                         va='center')
                # axes[i // 7, i % 7 + 1].set_visible(False)
            if i % 7 != 0:
                axes[i // 7, i % 7].tick_params(left=False)  # hide tick marks
            axes[i // 7, i % 7].set_ylim([0, 5])
            axes[i // 7, i % 7].fill_between(np.arange(0.0, 2.0, 0.01), 4, 5, facecolor='grey', alpha=0.5, zorder=0)
        file_name = '-'.join([self.task, 'heatmap', metric]) + '.jpg'
        plt.savefig(self.save_dir + file_name, dpi=150, bbox_inches='tight')
        plt.show()

    def plot_task_bar(self, shared_metric=False):
        """
        plot downstream metric in a single image.

        :return:
        """
        if self.task == 'assign':
            metrics = ['singlecellnet_F1', 'scmap_cell_F1', 'singleR_F1']
            drop_methods = ['cv2', 'scGeneFit']
            ylim, ylim_last = (0, 1.05), (0, 1.05)
        else:
            metrics = ['seurat_ARI', 'sc3_ARI']
            drop_methods = ['cv2']
            ylim, ylim_last = (0, 1), (0, 1)
        if shared_metric:
            metrics = ['markers_found', 'MRR']
            ylim, ylim_last = (0, 0.6), (0, 0.1)
        fig, axes = plt.subplots(len(metrics), 1, figsize=(8, 8))
        for i, metric in enumerate(metrics):
            spec_metric_table = pd.DataFrame(np.zeros((len(self.n_genes), len(self.task_methods))),
                                             index=self.n_genes, columns=self.task_methods)  # row method, # col n_gene
            for dataset in self.datasets:
                for n_gene in self.n_genes:
                    record = self.task_records[dataset + '-' + n_gene]
                    ratio_factor = record.n_marker_contain if metric == 'markers_found' else 1.0
                    spec_metric_table.loc[n_gene, :] += record.summary().loc[metric, :] / ratio_factor
            spec_metric_table /= len(self.datasets)

            # start to plot
            # spec_metric_table.drop(columns=drop_methods, inplace=True)
            spec_metric_table.rename(columns=formal_method_names, inplace=True)
            if i != len(metrics) - 1:
                spec_metric_table.plot(kind='bar', ax=axes[i], rot=0, title=metrics[i], cmap='tab20', width=0.8,
                                       ylim=ylim, legend=False)
            else:
                spec_metric_table.plot(kind='bar', ax=axes[i], rot=0, title=metrics[i], cmap='tab20', width=0.8,
                                       ylim=ylim_last)
                axes[i].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=4)
            # annotate number
            self.mark_number(axes[i], spec_metric_table)
        fig.set_tight_layout(True)
        partial_fig_name = '-bar-downstream_metrics.jpg' if not shared_metric else '-bar-markers_and_MRR.jpg'
        fig.savefig(self.save_dir + self.task + partial_fig_name, dpi=150, bbox_inches='tight')
        fig.show()

    def plot_rare_type_bar(self):
        rare_dataset = ['PBMC10000', 'xin']  # only use PBMC10000
        rare_f1 = ['singlecellnet_F1_rare', 'scmap_cell_F1_rare', 'singleR_F1_rare']
        fig, axes = plt.subplots(len(rare_f1), 1, figsize=(8, 8))

        for i, metric in enumerate(rare_f1):
            spec_metric_table = pd.DataFrame(np.zeros((len(self.n_genes), len(self.task_methods))),
                                             index=self.n_genes, columns=self.task_methods)
            for n_gene in self.n_genes:
                record = self.task_records[rare_dataset[0] + '-' + n_gene]
                spec_metric_table.loc[n_gene, :] += record.summary().loc[metric, :]

            # start to plot
            spec_metric_table.rename(columns=formal_method_names, inplace=True)
            if i != len(rare_f1) - 1:
                spec_metric_table.plot(kind='bar', ax=axes[i], rot=0, title=rare_f1[i], cmap='tab20', width=0.8,
                                       ylim=(0, 1.05),
                                       legend=False)
            else:
                spec_metric_table.plot(kind='bar', ax=axes[i], rot=0, title=rare_f1[i], cmap='tab20', width=0.8,
                                       ylim=(0, 1.05))
                axes[i].legend(loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=4)
            self.mark_number(axes[i], spec_metric_table)
        plt.tight_layout()
        plt.savefig(self.save_dir + 'assign-rare_type_f1.jpg', dpi=150, bbox_inches='tight')
        plt.show()

    def mark_number(self, ax, table: pd.DataFrame):
        rects = [ch for ch in ax.get_children() if isinstance(ch, Rectangle)]
        y_offset = 1e-2 if table.stack().max() > 0.1 else 1e-3
        for i, method in enumerate(table.columns):
            for j, n_gene in enumerate(self.n_genes):
                ax.text(rects[2 * i + j].get_x() + rects[2 * i + j].get_width() * 0.5,
                        rects[2 * i + j].get_height() + y_offset,
                        s='{:.2f}'.format(table.loc[n_gene, method]),
                        fontsize='xx-small' if self.task == 'assign' else 'x-small',
                        ha='center', va='baseline')

    def plot_box(self, metric: str):
        plt.style.use("ggplot")
        valid_method = [method for method in self.task_methods if method not in ('cv2', 'scGeneFit')]
        min_value = 1
        if self.task == 'assign':
            fig, axes = plt.subplots(2, 6, sharey='all', figsize=(10, 6))
            fig.subplots_adjust(wspace=0.05)
            for i, method in enumerate(valid_method):
                spec_metric_table = pd.DataFrame(np.zeros(shape=(5, len(self.n_genes))), columns=self.n_genes)
                for j, n_gene in enumerate(self.n_genes):
                    record = self.task_records['segerstolpe-' + n_gene]
                    spec_metric_table.loc[:, n_gene] = getattr(record, metric)[method]
                sns.boxplot(data=spec_metric_table, ax=axes[i // 6, i % 6])
                # axes[i // 6, i % 6].set_yscale('log')
                current_min = np.round(spec_metric_table.stack().min(), 1)
                min_value = current_min if current_min < min_value else min_value

            min_value -= 0.05
            interval = 0.05 if 1 - min_value > 0.1 else 0.01
            fig.text(0.04, 0.5, metric.replace('_', '-'), va='center', rotation='vertical')
            for k, method in enumerate(valid_method):
                if k % 6 != 0:
                    axes[k // 6, k % 6].tick_params(left=False)
                axes[k // 6, k % 6].text(0.5, 1.15 - 0.15 * min_value, formal_method_names[method],
                                         ha='center', va='center')
                axes[k // 6, k % 6].fill_between(np.arange(-0.5, 1.5, 0.01),  # x
                                                 1.05 - 0.05 * min_value, 1.3 - 0.3 * min_value,  # y
                                                 facecolor='grey', alpha=0., zorder=0)
                axes[k // 6, k % 6].set_yticks(np.arange(min_value, 1, interval))
                axes[k // 6, k % 6].set_ylim([min_value, 1.3 - 0.3 * min_value])
        else:
            fig, axes = plt.subplots(1, 6, sharey='all', figsize=(10, 4))
            fig.subplots_adjust(wspace=0.05)
            for i, method in enumerate(valid_method):
                spec_method_list = []
                for j, n_gene in enumerate(self.n_genes):

                    left_table = pd.DataFrame(np.zeros((10, 3)), columns=['value', 'n_gene', 'split'])
                    left_table['n_gene'], left_table['split'] = n_gene, 'S1'
                    right_table = pd.DataFrame(np.zeros((10, 3)), columns=['value', 'n_gene', 'split'])
                    right_table['n_gene'], right_table['split'] = n_gene, 'S2'

                    record = self.task_records['segerstolpe-' + n_gene]
                    for k in range(10):
                        single_run = eval(getattr(record, metric)[method].loc[k])
                        left_table.loc[k, 'value'], right_table.loc[k, 'value'] = single_run[0], single_run[1]
                    spec_method_list.append(pd.concat([left_table, right_table]).reset_index(drop=True))
                single_image_data = pd.concat(spec_method_list).reset_index(drop=True)
                current_min = np.round(single_image_data['value'].min(), 1)
                min_value = current_min if current_min < min_value else min_value
                sns.boxplot(x='n_gene', y='value', hue='split', data=single_image_data, ax=axes[i])
                if i == 0:
                    axes[i].set_ylabel(metric.replace('_', '-'))
                    axes[i].legend_.remove()
                elif i == len(valid_method) - 1:
                    axes[i].set_ylabel("")
                    axes[i].tick_params(left=False)
                    axes[i].legend(loc='center left', bbox_to_anchor=(1.1, 0.5),
                                   ncol=1, facecolor='white', frameon=False)
                else:
                    axes[i].set_ylabel("")
                    axes[i].legend_.remove()
                    axes[i].tick_params(left=False)
            min_value -= 0.1
            for k, method in enumerate(valid_method):
                axes[k].set_xlabel('')
                axes[k].text(0.5, 1.15 - 0.15 * min_value, formal_method_names[method], ha='center', va='center')
                axes[k].fill_between(np.arange(-0.7, 1.7, 0.01),  # x
                                     1.05 - 0.05 * min_value, 1.25 - 0.25 * min_value,  # y
                                     facecolor='grey', alpha=0.5, zorder=0)
                axes[k].set_yticks(np.arange(min_value, 1.05, 0.1))
                axes[k].set_ylim([min_value, 1.25 - 0.25 * min_value])
        file_name = '-'.join([self.task, 'box', metric]) + '.jpg'
        plt.savefig(self.save_dir + file_name, dpi=150, bbox_inches='tight')
        plt.show()

    def plot_markers_found_rate(self, min_quantile: float = 0.1):
        table = pd.DataFrame(np.zeros((len(self.n_genes), len(self.task_methods))), index=self.n_genes)
        for n_gene in self.n_genes:
            for dataset in self.datasets:
                record = self.task_records[dataset + '-' + n_gene]
                table.loc[n_gene, :] = record.markers_found.mean(axis=0) / record.n_marker_contain
        table /= len(self.datasets)
        table.plot()
        plt.show()

    def plot_all(self):
        shared_metrics = ['markers_found', 'MRR']
        if self.task == 'assign':
            task_spec_metrics = ['singlecellnet_F1', 'scmap_cell_F1', 'singleR_F1']
            self.plot_rare_type_bar()
        else:
            task_spec_metrics = ['seurat_ARI', 'sc3_ARI']
        for metric in shared_metrics + task_spec_metrics:
            self.plot_metric_heatmap(metric)
        for metric in task_spec_metrics:
            self.plot_box(metric)
        self.plot_task_bar(shared_metric=True)
        self.plot_task_bar()
