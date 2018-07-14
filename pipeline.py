import PolyLibScan as pls
import MDlatticeAnalysisTool as mdl
from epitopsy import EnergyGrid as EG
from epitopsy import DXFile as dxf
import pathlib2 as pl
import collections as col
import matplotlib.pyplot as plt
import sklearn as sk
import shutil
import yaml

# Helpers
def project(project_folder):
        project = pls.Analysis.Project(project_folder)
        return project

def calc_run_average(counter_list):
    cc_protein = {}
    for i in range(len(counter_list)):
        for key in counter_list[i]:
            if key in cc_protein:
                cc_protein[key] = cc_protein[key] + counter_list[i][key]
            else:
                cc_protein[key] = counter_list[i][key]
    average_contact_counts_prot = {key: value / len(counter_list) for key,
            value in cc_protein.iteritems()}
    return average_contact_counts_prot


# Use PolyLibScan to create pdb/pqr files of all polymer poses/snapshots
# and the target protein.
class PolyLibScan_call(object):
    def __init__(self, project, project_folder, output_folder_full_path,
            job = 0, run = 0):
        self.run = run
        self.project = project
        self.sample = project.jobs[job][run]
        self.output_folder = output_folder_full_path
        self.protein_name = self.sample.meta()['protein']
        self.polymer_name = self.sample.meta()['polymer_name']
        self.sampling_rate = int(self.sample.meta()['sampling_rate'])
        self.pose_file = self.sample.pymol.pose.file_name(molecule = 'polymer',
                state = 'full')

    def create_pose_pdb(self, molecule='polymer', state='full'):
        gen = self.sample.pymol.pose.pdb_data()
        output_path = self._save_pdb(self.pose_file, gen)

    def _save_pdb(self, file_name, pdb_data):
        output_path = self.output_folder+file_name
        with open(output_path, 'w') as f:
            f.write(self.sample.pymol.pose.template.render(pymol=pdb_data))

    def create_pqr(self):
        pqr_obj = pls.helpers.to_pqr.to_pqr(self.sample)
        pqr_obj.save_pqr(self.output_folder+self.protein_name + '.pqr',
                molecule='protein', state='end')
        pqr_obj.save_pqr(self.output_folder+self.polymer_name + '.pqr',
                molecule='polymer', state='start')

    def copy_protein_pdb(self):
        shutil.copy2(self.sample.job.project.protein_path, self.output_folder)

# Use MDL to search for multi-contacts in the PLS results.
class MDL_call(PolyLibScan_call):
    def __init__(self, project, project_folder, output_folder_full_path,
            job = 0, run = 0):
        super(MDL_call, self).__init__(project, project_folder,
                output_folder_full_path, job, run)
        self.analysis = self.mdl_run()

    def mdl_run(self):
        env = mdl.enviroment.Enviroment(self.output_folder + self.protein_name
                + '.pqr', self.output_folder + self.pose_file)
        print('Building kd-Tree and searching for nearest neighbors for run %i, this may take a while'
                % self.run)
        ana = mdl.analysis.Analytics(env)
        return ana

# Use Epitopsy to create an occupancy isosurface.
class Epitopsy_call(MDL_call):
    def __init__(self, project, project_folder, output_folder_full_path,
            job = 0, run = 0, meshsize=[4,4,4]):
        super(Epitopsy_call, self).__init__(project, project_folder,
            output_folder_full_path, job, run)
        self.meshsize = meshsize

    def epitopsy_run(self):
        EG.electrostatics(self.output_folder + self.protein_name + '.pdb',
                self.output_folder + self.polymer_name+'.pqr',
                cubic_box=False, verbose=True, mesh_size = self.meshsize,
                box_center = tuple(self.analysis._enviroment.geometric_center),
                centering = 'center')
        EG.scan(self.output_folder + self.protein_name + '.pdb',
                self.output_folder + self.polymer_name + '.pqr',
                APBS_dx_path = self.output_folder)
        if not pl.Path.exists(pl.PosixPath(self.output_folder
            + self.protein_name + '_epi.dx')):
            shutil.move('./' + self.protein_name + '_epi.dx',
                    self.output_folder + self.protein_name + '_epi.dx')
        if not pl.Path.exists(pl.PosixPath(self.output_folder
            + self.protein_name + '_mic.dx.gz')):
            shutil.move('./' + self.protein_name + '_mic.dx.gz',
                    self.output_folder + self.protein_name + '_mic.dx.gz')
        # Epitopsy energy to occupancy.
        EG.energy2occupancy(self.output_folder + self.protein_name
                + '_epi.dx',self.output_folder + self.protein_name + '_mic.dx.gz',
            self.output_folder + self.protein_name + '_occ.dx', normalize=True)
        EG.isosurface_ramp(self.output_folder + self.protein_name + '_occ.dx')

    # This method transforms the cartesian coordinates of protein
    # atoms to the corresponding epitopsy box coordinates.
    def mdl_coords_to_box(self, coord_list):
        box = dxf.read_dxfile(self.output_folder + self.protein_name
                + '_epi.dx', 'esp')
        if coord_list:
            multi_c_list = box.transform_real_to_box_space(coord_list)
        else:
            multi_c_list = []
        return multi_c_list

# Visualize the PLS/MDL multi-contacts and the Epitopsy occupancy
# isosurfaces, pymol has to be started with the -R option.
class Pymol_call(MDL_call):
    def __init__(self, project, project_folder, output_folder_full_path,
            job = 0, run = 0):
        super(Pymol_call, self).__init__(project, project_folder,
                output_folder_full_path, job, run)
        self._connect = mdl.pymol.connect()
        self.pymol = self.pymol_load()

    def pymol_load(self):
        protein_pdb_path = self.output_folder+self.protein_name+ '.pdb' 
        pym = mdl.pymol.Pymol(self.analysis, protein_pdb_path)
        return pym


class Pipeline_run(object):
    def __init__(self, project_folder, output_folder_full_path, job = 0,
            pymol = 'True', epitopsy = 'True', meshsize = [4,4,4]):
        self._epitopsy = epitopsy
        self.meshsize = meshsize
        self.project_folder = project_folder
        self.output_folder = output_folder_full_path
        self.job = job
        self.project = project(self.project_folder)
        self._pls_call = PolyLibScan_call(self.project, self.project_folder,
                self.output_folder, job=self.job)
        self.coord_list = self.run()
        self._epi_call = self.epitopsy(self.meshsize)
        self._probability_plot_polymer = self.probability_plot_polymer()

    def run(self):
        self._pls_call.create_pqr()
        self._pls_call.copy_protein_pdb()
        counter_list = []
        coord_list = []
        coord_dict = col.OrderedDict()
        mono_list = []
        for i in range(len(self.project.jobs[self.job])):
            try:
                PolyLibScan_call(self.project, self.project_folder,
                        self.output_folder,job=self.job, run = i).create_pose_pdb()
                mdl_c = MDL_call(self.project, self.project_folder,
                        self.output_folder,job=self.job, run = i)
                for n in range(len(mdl_c.analysis.multi_contact_list)):
                    for mono in list(mdl_c.analysis.multi_contact_list[n]):
                        if mdl_c.analysis.multi_contact_list[n][mono].loc['Distance'] < 4:
                            if tuple(mdl_c.analysis.multi_contact_list[n][
                                mono].loc['Monomer coord']) not in coord_dict.keys():
                                coord_dict[tuple(
                                    mdl_c.analysis.multi_contact_list[n][
                                        mono].loc['Monomer coord'])] = 1
                            else:
                                coord_dict[tuple(
                                    mdl_c.analysis.multi_contact_list[n][
                                        mono].loc['Monomer coord'])] += 1
                            coord_list.append(mdl_c.analysis.multi_contact_list[n][
                                mono].loc['Monomer coord'])
                            mono_list.append(mono)
                counter_list.append(mdl_c.analysis._contact_frequencies_protein(
                    mdl_c.analysis._multi_concat))
            except (ValueError, AttributeError):
                # Both errors happen when the multi_concat list is empty,
                # if there was no multi-contact event.
                pass
        with open(self.output_folder + self._pls_call.protein_name
                + '_coords_dict_job' + str(self.job) + '.yml', 'w') as h_:
            yaml.dump(coord_dict, stream=h_)
        with open(self.output_folder + self._pls_call.protein_name
                + '_coords_job' + str(self.job) + '.yml', 'w') as g_:
            yaml.dump(coord_list, stream=g_)
        with open(self.output_folder + self._pls_call.protein_name
                + '_monos_job' + str(self.job) + '.yml', 'w') as f_:
            yaml.dump(mono_list, stream=f_)
        avg_cc_prot = calc_run_average(counter_list)
        if avg_cc_prot.has_key(0):
            del avg_cc_prot[0]
        for residue in avg_cc_prot.keys():
            if int(avg_cc_prot[residue]*100) == 0:
                del avg_cc_prot[residue]
        with open(self.output_folder + self._pls_call.protein_name + '_job'
                + str(self.job) + '.yml', 'w') as f_:
            yaml.dump(avg_cc_prot, stream=f_)
        f_.close()
        return coord_list

    def epitopsy(self, meshsize):
        if self._epitopsy == 'True':
            epi_c = Epitopsy_call(self.project, self.project_folder,
                    self.output_folder, job=self.job, meshsize = self.meshsize)
            epi_c.epitopsy_run()
            multi_c_list = epi_c.mdl_coords_to_box(self.coord_list)
            with open(self.output_folder + self._pls_call.protein_name
                    + '_box_coords_job' + str(self.job) + '.yml', 'w') as f_:
                yaml.dump(multi_c_list, stream=f_)

    def _mean_absolute_error(self, pls_c, mdl_c):
        if self._epitopsy == 'True':
            try:
                box = dxf.read_dxfile(self.output_folder + pls_c.protein_name
                        + '_occ.dx', 'esp')
                with open(self.output_folder + pls_c.protein_name
                        + '_box_coords_job' + str(self.job) + '.yml') as f_:
                    coords = yaml.load(f_)
                coords_tuple = []
                for i in range(len(coords)):
                    coords_tuple.append(tuple(coords.tolist()[i]))
                counter = col.Counter(coords_tuple)
                freq_occ_dict = col.OrderedDict()
                for key in counter:
                    freq_occ_dict[key] = counter[key], box.box[key]
                frequency = []
                occupancy = []
                for key in freq_occ_dict:
                    frequency.append(freq_occ_dict[key][0])
                    occupancy.append(freq_occ_dict[key][1])
                probability = [float(freq)/((len(mdl_c.analysis.df_dist))*int(
                    self.project.jobs[0].meta['sampling_rate']))
                    for freq in frequency]
                return sk.metrics.mean_absolute_error(probability, occupancy)
            except ValueError:
                return None

    def probability_plot_polymer(self):
        mdl_c = MDL_call(self.project, self.project_folder, self.output_folder,
                job=self.job)
        pls_c = PolyLibScan_call(self.project, self.project_folder,
                self.output_folder, job=self.job)
        with open(self.output_folder + pls_c.protein_name + '_monos_job'
                + str(self.job) + '.yml', 'r') as f_:
            monos = yaml.load(f_)
        c_monos = col.Counter(monos)
        prob_monos = col.OrderedDict()
        for key in c_monos:
            prob_monos[key] = float(c_monos[key])/(len(mdl_c.analysis.df_dist)
                    * pls_c.sampling_rate)
        prob_monos = col.OrderedDict(sorted(prob_monos.iteritems()))
        plt.bar(range(len(prob_monos)), list(prob_monos.values()),
                align='center')
        plt.xticks(range(len(prob_monos)), list(prob_monos.keys()))
        fig = plt.gcf()
        fig.set_size_inches(21, 5)
        plt.xlabel('Monomer (Position, Name)')
        plt.ylabel('Probability')
        plt.savefig(self.output_folder + self._pls_call.polymer_name
                + 'polymer_probability.png', dpi = 300)
        plt.clf()
        plt.close()

    def pymol(self, epitopsy = 'True'):
        pymol_c = Pymol_call(self.project, self.project_folder,
                self.output_folder, job=self.job)
        pymol_c.pymol.setup()
        if epitopsy == 'True':
            pymol_c.pymol.py_con.show_as('surface')
            pymol_c.pymol.py_con.bg_color('white')
            pymol_c.pymol.py_con.set('antialias', 2)
            pymol_c.pymol.py_con.set('surface_quality', 1)
            pymol_c.pymol.py_con.load(self.output_folder
                    + self._pls_call.protein_name + '_epi.dx')
            pymol_c.pymol.py_con.load(self.output_folder
                    + self._pls_call.protein_name + '_esp.dx')
            pymol_c.pymol.py_con.load(self.output_folder
                    + self._pls_call.protein_name + '_occ.dx')
            pymol_c.pymol.py_con.load(self.output_folder
                    + self._pls_call.protein_name + '_occ.pml')
            pymol_c.pymol.py_con.color('forest', '*occ*')
            pymol_c.pymol.py_con.isosurface('epi_neg_1kbT',
                    self._pls_call.protein_name + '_epi', -1)
            pymol_c.pymol.py_con.isosurface('epi_pos_1kbT',
                    self._pls_call.protein_name + '_epi', +1)
            pymol_c.pymol.py_con.isosurface('epi_neg_2kbT',
                    self._pls_call.protein_name + '_epi', -2)
            pymol_c.pymol.py_con.isosurface('epi_pos_2kbT',
                    self._pls_call.protein_name + '_epi', +2)
            pymol_c.pymol.py_con.set('transparency', 0.5, 'epi_neg_1kbT')
            pymol_c.pymol.py_con.set('transparency', 0.5, 'epi_pos_1kbT')
            pymol_c.pymol.py_con.color('firebrick', 'epi_neg*')
            pymol_c.pymol.py_con.color('skyblue', 'epi_pos*')
            pymol_c.pymol.py_con.disable('epi_pos*')
            pymol_c.pymol.py_con.disable('epi_neg*')
            pymol_c.pymol.py_con.isosurface('esp_neg_1kbT',
                    self._pls_call.protein_name + '_esp', -1)
            pymol_c.pymol.py_con.isosurface('esp_pos_1kbT',
                    self._pls_call.protein_name + '_esp', +1)
            pymol_c.pymol.py_con.color('salmon', 'esp_neg*')
            pymol_c.pymol.py_con.color('lightblue', 'esp_pos*')
            pymol_c.pymol.py_con.disable('esp_pos*')
            pymol_c.pymol.py_con.disable('esp_neg*')
        with open(self.output_folder + self._pls_call.protein_name
                + '_job' + str(self.job) + '.yml', 'r') as f_:
            avg_cc_prot = yaml.load(f_)
        pymol_c.pymol.show_bindings(counter_list = avg_cc_prot)
