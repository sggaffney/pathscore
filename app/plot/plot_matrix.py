# import os
# from datetime import datetime

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.collections as collections
from matplotlib.transforms import Bbox

from utils import get_target_bbox


class GeneMatrix:
    """Holds patient-gene matrix info for a single pathway.
    Write matrix to file with call to export_matrix."""

    def __init__(self, genepatients_dict, annot_dict=None, 
                 exclusive_genes=None, hypermutated=None,
                 n_all_patients=None):
        """Add gene-patientList pairs into dictionary, and optional annot_dict.

        Args:
            genepatients_dict (dict): {gene: patient_list}
            annot_dict (dict): {(hugo, patient): annot (str)}
            exclusive_genes (list)
        """
        self.genepatients = genepatients_dict
        self.annot_dict = annot_dict if annot_dict is not None else {}
        self.exclusive_genes = exclusive_genes
        self.hypermutated = hypermutated
        self.n_all_patients = n_all_patients
        self.matrix = self._build_matrix()
        
    def _build_matrix(self):
        """ build gene-patient dataframe, sorted for display.
                  a      e      b      c      d      f
        KRAS      k   True      m   True   True  False
        TP53   True   True   True  False  False   True
        EGFR   True      e  False  False  False  False
        BRAF  False  False  False  False  False       
        """

        use_genes = set(self.genepatients)
        patients = set()
        for gene in self.genepatients:
            patients = patients.union(set(self.genepatients[gene]))

        annot_dict_use = dict([i for i in self.annot_dict.items() if i[0][0] in use_genes])
        patients = sorted(list(patients))
        # sort genes
        genes = sorted(list(use_genes))
        genes.sort(key=lambda g: g in self.exclusive_genes, reverse=True)
        genes.sort(key=lambda g: len(self.genepatients[g]), reverse=True)

        # This is original 'matrix' input to matlab
        m = pd.DataFrame(index=genes, columns=patients, data=False)
        for gene, patient_list in self.genepatients.items():
            m.loc[gene, patient_list] = True
        for gene, patient in annot_dict_use:
            m.loc[gene, patient] = annot_dict_use[(gene, patient)]

        # m_status = ~(m == False)
        m_status = m.applymap(lambda x: x is not False)

        # RE-SORT BY WEIGHT (powers of 2)
        self.n_genes = len(use_genes)
        twos = np.power(2, range(self.n_genes-1, -1, -1), dtype=np.float64)  # gene weights

        weights = m_status.apply(lambda s: sum(s*twos), axis=0)
        patient_order = weights.reset_index().\
            sort_values([0, 'index'], axis=0, ascending=[False, True])['index'].values

        matrix = m.reindex(columns=patient_order)
        return matrix

    # def export_matrix(self, outfile):
    #     """Writes tab-separated matrix file for patient/gene
    #     pairs in pathway."""
    #     pass

    def plot_matrix(self, fontsize=8, max_label=150, box_px=20, show_limits=False):
        mp = MatrixPlotter(self.matrix, n_all_patients=self.n_all_patients, 
                           exclusive_genes=self.exclusive_genes, 
                           hypermutated=self.hypermutated)
        hfig, ax = mp.draw_matrix(fontsize=fontsize, max_label=max_label,
                                  show_limits=show_limits, box_px=box_px)
        return hfig, ax

    def save_matrix(self, out_path=None, **plot_kwargs):
        # if not out_path:
        #     d = datetime.utcnow()
        #     maf_dir = os.path.dirname(self.maf_path) \
        #         if self.maf_path else os.getcwd()
        #     fig_name = 'matrix_' + d.strftime('%Y-%m-%d_%H%M%S%f') + '.pdf'
        #     out_path = os.path.join(maf_dir, fig_name)
        hfig, ax = self.plot_matrix(**plot_kwargs)
        hfig.savefig(out_path)  # from app.plot.plot_dendrogram import plot_dendrogram
        return out_path


class MatrixPlotter(object):

    def __init__(self, matrix_df, exclusive_genes=None, n_all_patients=None,
                 hypermutated=None):
        """Generate matrix plot using matrix dataframe.

        Args:
            matrix_df (pandas DataFrame): columns are hugo symbols, rows are
                patient_ids. Contents are boolean for mutation, or string for
                annotation.
        """
        self.matrix_df = matrix_df

        # self.m_status = ~(matrix_df == False)
        self.m_status = matrix_df.applymap(lambda x: x is not False)
        # assert isinstance(self.m_status, pd.DataFrame)
        self.use_genes = list(matrix_df.index)
        self.use_patients = list(matrix_df.columns)
        self.exclusive_genes = exclusive_genes
        self.hypermutated = hypermutated

        self.upi = float(96)  # units per inch
        
        self.n_genes = len(self.use_genes)
        self.n_patients = len(self.use_patients)

        # GENE LABELS
        if n_all_patients:
            p_counts = self.m_status.sum(axis=1)
            gene_labels = []
            format_str = "{gene} {pc:3.0f}%"
            for gene in self.use_genes:
                gene_labels.append(format_str.format(
                    gene=gene, pc=float(p_counts[gene]) / n_all_patients * 100))
            self.gene_labels = gene_labels
            n_mutated = sum(matrix_df.any(axis=0))
            pc_patients = n_mutated / float(n_all_patients) * 100
            patient_str = 'patient' if self.n_patients == 1 else 'patients'
            self.title = "{n} {patients} ({pc:.1f}%)".format(n=n_mutated,
                                                             patients=patient_str,
                                                             pc=pc_patients)
        else:
            self.gene_labels = self.use_genes
            self.title = "{n} patients".format(n=self.n_patients)
        patient_labels = []
        for i in self.use_patients:
            label = '{}*'.format(i) if i in self.hypermutated else i
            patient_labels.append(label)
        self.patient_labels = patient_labels

    def draw_matrix(self, fontsize=10, max_label=100, box_px=30, show_limits=False):
        """Generate matrix figure."""
        # figsize = [(self.ax_x['length'] + self.ax_x['padding']) / self.upi,
        #            (self.ax_y['length'] + self.ax_y['padding']) / self.upi]

        pbox_lgene = float(box_px)  # gene box width (data units)
        pbox_lpatient = float(box_px)  # patient box width (data units)

        self.ax_x = dict(length=pbox_lpatient * self.n_patients,
                         # padding=float(90),
                         num=self.n_patients,
                         labels=self.patient_labels,
                         boxlen=pbox_lpatient)
        self.ax_y = dict(length=pbox_lgene * self.n_genes,
                         # padding=float(70),
                         num=self.n_genes,
                         labels=self.gene_labels,
                         boxlen=pbox_lgene)

        figsize = [self.ax_x['length'] / self.upi,
                   self.ax_y['length'] / self.upi]
        subplot_kw = {'aspect': 'equal', 'position': [0, 0, 1, 1],
                      'xlim': [0, self.ax_x['length']],
                      'ylim': [0, self.ax_y['length']]}
        hfig, ax = plt.subplots(1, 1, figsize=figsize, subplot_kw=subplot_kw)
        ax.invert_yaxis()
        obj_list = []  # for calculating figure extents
        # ADD GRAY BACKGROUND
        rect = patches.Rectangle((0, 0), self.ax_x['length'],
                                 self.ax_y['length'], fc=[0.85, 0.85, 0.85],
                                 ec='none')
        patch = ax.add_patch(rect)
        obj_list.append(patch)
        # ADD MUTATION BOXES
        m_reindex = self.m_status.set_index(np.arange(0, self.n_genes))
        m_reindex.columns = np.arange(0, self.n_patients)
        s = m_reindex.stack()
        mut_inds = s[s].index.values
        for g_ind, p_ind in mut_inds:
            gene = self.use_genes[g_ind]
            # patient = self.use_patients[p_ind]
            if gene in self.exclusive_genes:
                self._add_box(p_ind, g_ind, ax, color='red')
            else:
                self._add_box(p_ind, g_ind, ax, color='blue')
            annot = self.matrix_df.iloc[g_ind, p_ind]
            if annot is not np.bool_(True) and annot is not True:
                self._add_annot(p_ind, g_ind, ax, annot, fontsize=fontsize,
                                fontweight='bold')
        # WHITE LINE GRID
        self._add_grid(ax)

        ax.set_frame_on(False)
        ax.axes.get_xaxis().set_visible(False)
        ax.axes.get_yaxis().set_visible(False)
        ax.tick_params(axis=u'both', which=u'both', length=0)

        labels = self._add_ax_labels(fontsize=fontsize)
        obj_list.extend(labels)

        # ADJUST FIGURE
        bb = get_target_bbox(obj_list, ax, hfig)
        # print(bb)

        # RIGHT WHITESPACE PADDING FOR KNOWN WIDTH
        if max_label is not None:
            title_overflow = max(bb.xmax - self.ax_x['length'], 0)
            r_padding = max([max_label + bb.xmin - title_overflow, 0])
        else:
            r_padding = 0
        # VERTICAL PADDING TO AVOID CROPPING TEXT
        v_padding = 5
        h_padding = 5
        
        # print("R padding: {}".format(r_padding))
        xmin = bb.xmin - h_padding
        xmax = bb.xmax + h_padding + r_padding
        ymin = bb.ymin - v_padding
        ymax = bb.ymax + v_padding

        bb_new = Bbox(np.array([[xmin, ymax], [xmax, ymin]]))
        # print(bb_new)
        if show_limits:
            rect = patches.Rectangle(bb_new.p0, bb_new.width, bb_new.height, alpha=0.2)
            ax.add_artist(rect)

        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymax, ymin])

        fwidth_in = (xmax - xmin) / self.upi
        fheight_in = (ymax - ymin) / self.upi
        hfig.set_figwidth(fwidth_in)
        hfig.set_figheight(fheight_in)
        
        return hfig, ax

    def _add_grid(self, ax):
        """Add white grid lines.

        One line 'after' each box excluding last one.
        ax_x/y has 'labels', 'padding', 'length', 'num', 'boxlen'

        Args:
            ax: Axis handle
        """
        n_x = int(self.ax_x['length'] / self.ax_x['boxlen'])
        n_y = int(self.ax_y['length'] / self.ax_y['boxlen'])
        # draw X grid
        if n_x > 1:
            xy_start = np.vstack((self.ax_x['boxlen'] *
                                  np.array(range(1, n_x, 1)),
                                  np.zeros(n_x - 1)))
            xy_end = np.vstack((self.ax_x['boxlen'] *
                                np.array(range(1, n_x, 1)),
                                self.ax_y['length'] * np.ones(n_x - 1)))
            segments = zip(zip(*xy_start), zip(*xy_end))
            lc_x = collections.LineCollection(segments, colors=(1, 1, 1, 1),
                                              linewidths=2)
            ax.add_collection(lc_x)
        # draw Y grid
        if n_y > 1:
            xy_start = np.vstack((np.zeros(n_y - 1),
                                  self.ax_y['boxlen'] *
                                  np.array(range(1, n_y, 1))))
            xy_end = np.vstack((self.ax_x['length'] * np.ones(n_y - 1),
                                self.ax_y['boxlen'] *
                                np.array(range(1, n_y, 1))))
            segments = zip(zip(*xy_start), zip(*xy_end))
            lc_y = collections.LineCollection(segments, colors=(1, 1, 1, 1),
                                              linewidths=2)
            ax.add_collection(lc_y)

        ax.add_patch(
            patches.Rectangle((0, 0), self.ax_x['length'], self.ax_y['length'],
                              fc='none', ec='white')
        )

    def _add_box(self, x_i, y_i, ax, color='red'):
        """Add rectangle for mutation using x-index and y-index."""
        l = x_i * self.ax_x['boxlen']
        b = y_i * self.ax_y['boxlen']
        w = self.ax_x['boxlen']
        h = self.ax_y['boxlen']
        rect = patches.Rectangle((l, b), w, h, fc=color, ec='none')
        ax.add_patch(rect)

    def _add_annot(self, x_i, y_i, ax, annot, fontsize=8, fontweight='bold', color='w'):
        """Add annotation string to mutation box using x,y index."""
        x = x_i * self.ax_x['boxlen'] + self.ax_x['boxlen']/2
        y = y_i * self.ax_y['boxlen'] + self.ax_y['boxlen']/2
        t = ax.text(x, y, annot, ha='center', va='center',
                    fontsize=fontsize, color=color, fontweight=fontweight)
        return t

    def _add_ax_labels(self, fontsize=10):
        """Add manual x and y labels"""
        labels = []
        for ind, ylabel in enumerate(self.ax_y['labels']):
            x = 0
            y = ind * self.ax_y['boxlen'] + self.ax_y['boxlen']/2
            t = plt.text(x, y, ylabel + ' ', ha='right', va='center', 
                         fontsize=fontsize)
            labels.append(t)
        for ind, xlabel in enumerate(self.ax_x['labels']):
            x = ind * self.ax_x['boxlen'] + self.ax_x['boxlen']/2
            y = self.ax_y['length']
            t = plt.text(x, y, xlabel + ' ', ha='right', va='top', 
                         fontsize=fontsize, rotation=45, 
                         rotation_mode="anchor")
            labels.append(t)
        # Add title
        x = 0 # self.ax_x['length']
        y = 0
        t = plt.text(x, y, self.title, ha='left', va='bottom', fontsize=fontsize)
        labels.append(t)
        return labels
