import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import collections as col
import itertools
import math

# helper function
def square_distance(a, b):
    s = 0
    for x, y in itertools.izip(a, b):
        d = x - y
        s += d * d
    return s

Node = col.namedtuple('Node', 'location axis left_child right_child')

class Analytics(object):
    def __init__(self, enviroment, k=3, distance_cutoff = 4):
        self.dist_cut = distance_cutoff
        self._enviroment = enviroment
        self.poly_poses_coords = self._enviroment.poly_poses_coords
        self.dim = k
        self._coord_list = self._enviroment.residue_list
        self.root = self._build_tree(list(self._coord_list))
        self.dataframe = self._build_dataframe()
        self.df_dist = self.dataframe.loc['Distance']
        self.multi_contact_list = self._find_multi_contacts(self.df_dist,
                self.dataframe)
        self._multi_concat = self._concatenate_multi_contact_list(
                self.multi_contact_list)

        # k equals the number of dimensions
    def _build_tree(self, coord_list, depth=0):
        if coord_list == []:
            return None
        else:
            # cycle through axis
            current_axis = depth % self.dim

            # sort the coordinate list according to the current axis
            coord_list.sort(key=lambda x: x[current_axis])
            median = len(coord_list) / 2

            # build the tree recursivly
            return Node(location = coord_list[median], axis = current_axis,
                        left_child=self._build_tree(coord_list[:median],
                            depth + 1),
                        right_child=self._build_tree(coord_list[median +1:],
                            depth + 1))

    def _nearest_neighbor(self, query):
        best_hit = [None, float('inf')]

        def recursive_search(position):
            if position is None:
                return 
            location, axis, left_child, right_child = position

            position_sd = square_distance(location, query)
            if position_sd < best_hit[1]:
                best_hit[:] = location, position_sd

            diff = query[axis] - location[axis]
            close, away = (left_child, right_child) if diff <= 0 else (
                    right_child, left_child)

            recursive_search(close)
            if diff ** 2 < best_hit[1]:
                recursive_search(away)

        recursive_search(self.root)
        return best_hit[0], math.sqrt(best_hit[1])

    def _calc_nearest_neighbor_list(self):
        nn_list = []
        for poly in range(len(self.poly_poses_coords)):
            for mono in self.poly_poses_coords[poly]:
                nn_list.append(self._nearest_neighbor(mono))
        return nn_list

    def _build_dataframe(self):
        pp_coords = self.poly_poses_coords
        nn_list = self._calc_nearest_neighbor_list()
        composite_nn_list = [nn_list[x:x+len(pp_coords[0])] for x in range(0,
            len(nn_list),len(pp_coords[0]))]
        df_list = []
        for i in range(len(pp_coords)):
            df = pd.DataFrame(columns = enumerate([residue.get_resname()
                for residue in self._enviroment._poly_poses[i].get_residues()], 1),
                              index = ['Monomer coord', 'Protein coord',
                                  'Residue Id', 'Distance'])
            df.columns.name = 'Pose %i' % (i+1)
            df.loc['Monomer coord'] = pp_coords[i]
            df.loc['Protein coord'] = [item[0][:3] for item in composite_nn_list[i]]
            df.loc['Residue Id'] = [item[0][3] for item in composite_nn_list[i]]
            df.loc['Distance'] = [item[1] for item in composite_nn_list[i]]
            df_list.append(df)
        return pd.concat(df_list)

    # Plot to show the probability of each monomer to bind to the
    # protein per snapshot / the probability of a monomer to be part
    # of a "multi-contact" snapshot.
    def probability_distribution_polymer(self, binding="multi"):
        """binding = single: Probability of each monomer to bind to the protein per snapshot,
        binding= multi: Probability of each monomer to be part of a snapshot with multiple binding events"""
        if binding == "single":
            self._plot_probability_distribution_polymer(self.df_dist)
        elif binding == 'multi':
            self._plot_probability_distribution_polymer(
                    self._multi_concat.loc['Distance'])
        else:
           raise ValueError('Only "multi" and "single" binding probability is supported')

    def _plot_probability_distribution_polymer(self, df):
        series = np.sum(df < self.dist_cut) / len(self.df_dist)
        series.plot.bar()
        plt.xlabel('Monomer')
        plt.ylabel('Probability')
        plt.show()

    def _find_multi_contacts(self, df_dist, dataframe):
        multi_contact_list = []
        for i in range(len(df_dist)):
            if sum(np.sum(df_dist.iloc[[i]] < self.dist_cut)) >= 2:
                multi_contact_list.append(dataframe.iloc[(i+1)*4-4:(i+1)*4])
        return multi_contact_list

    # This is neccesary so that if there are no multi-contacts/just
    # one for a run there wont be a pandas error that there is nothing
    # to concatenate.
    def _concatenate_multi_contact_list(self, multi_contact_list):
        if not multi_contact_list:
            multi_concat = pd.DataFrame(columns = enumerate(
                [residue.get_resname() for residue
                    in self._enviroment._poly_poses[0].get_residues()], 1),
                              index = 2*['Monomer coord', 'Protein coord',
                                  'Residue Id', 'Distance'])
                multi_concat.values.fill(0)
        elif len(multi_contact_list) == 1:
            multi_concat = multi_contact_list
        else:
            multi_concat = pd.concat(multi_contact_list)
        return multi_concat


    # Plot to show the probability of each residue to bind with a monomer per
    # snapshot / the probability of a residue to be part of a "multi-contact"
    # snapshot.
    def probability_distribution_protein(self, binding="multi"):
        """binding = single: Probability of a residue to bind with the polymer
        per snapshot, binding= multi: Probability of each residue to be
        part of a snapshot with multiple binding events"""
        if binding == "single":
            self._plot_probability_distribution_protein(self.dataframe)
        elif binding == "multi":
            self._plot_probability_distribution_protein(self._multi_concat)
        else:
           raise ValueError('Only "multi" and "single" binding probability is
                   supported')

    def _contact_frequencies_protein(self, df):
        res_list = []
        for i in range(len(df.loc['Distance'])):
            for mono in list(df):
                if df[mono].loc['Distance'].iloc[i] < self.dist_cut:
                    res_list.append(df[mono].loc['Residue Id'].iloc[i])
        list_count = col.Counter(res_list)
        for key in list_count.keys():
            list_count[key] = float(list_count[key]) / len(self.df_dist)
        return list_count

    def _plot_probability_distribution_protein(self, df):
        list_count = self._contact_frequencies_protein(df)
        for key in list_count.keys():
            list_count[key, self._enviroment.res_id_dict[key]] = list_count.pop(key)
        odict = col.OrderedDict(sorted(list_count.items()))
        plt.bar(range(len(odict)), list(odict.values()), align='center')
        plt.xticks(range(len(odict)), list(odict.keys()))
        plt.xlabel('Residue')
        plt.ylabel('Probability')
        fig = plt.gcf()
        fig.set_size_inches(8, 5)
        plt.show()

    # Create a dataframe of all residues and monomers taking part in multi
    # contact events and their relative frequencies of being part of such
    # an event.
    def multi_contact_frequencies(self, multi_concat):
        multi_concat = self._multi_concat
        df = pd.DataFrame(0, columns = [mono for mono in multi_concat],
                  index=set(itertools.chain.from_iterable(
                      multi_concat.loc['Residue Id'].values)))
        for mono in list(multi_concat):
            for i in range(len(multi_concat.loc['Distance'])):
                if multi_concat[mono].loc['Distance'].iloc[i] < self.dist_cut:
                    df.loc[multi_concat[mono].loc['Residue Id'][i]][mono] += 1
        for residue in df.index:
            if all(df.loc[residue] == 0):
                df = df.drop(residue)
        for mono in df.columns:
            if all(df[mono] == 0):
                df = df.drop(mono, axis=1)
        return(df.divide(len(self.df_dist)))
