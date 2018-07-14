import yaml
import collections as col
import PolyLibScan as pls
import pathlib2 as pl
import LammpsSubmit as LS
import epitopsy as epi
import MDlatticeAnalysisTool as mdl
import sklearn as sk
import operator
import random
import os
import time

class Genetic_algorithm(object):
    gene_list = ['BA', 'BP', 'Bor', 'CBS', 'Dan', 'Glu', 'NTA', 'Dod', 'C5-BP',
            'C6-BP', 'C6-CyH', 'C6-Sp', 'CBS-Boc', 'ED', 'Iso', 'NTA-Glu']
    generation = 0

    @classmethod
    def update_gen(cls):
        cls.generation += 1

    def __init__(self, local_project_folder, cluster_project_folder,
            best_sample = 5, lucky_few = 2, nr_children = 4,
            nr_generations = 6, mutation_rate = 0.02, meshsize = [4,4,4]):
        self.best_sample = best_sample
        self.lucky_few = lucky_few
        self.nr_children = nr_children
        self.nr_generations = nr_generations
        self.mutation_rate = mutation_rate
        self.meshsize = meshsize
        self.project_folder = local_project_folder
        self.cluster_folder = cluster_project_folder
        self._run = self.run(self.nr_generations, self.project_folder,
                self.cluster_folder, self.best_sample, self.lucky_few,
                self.nr_children)

    # This is the analysis part of the fitness function,
    # the fitness is calculated by calculating the mean absolute
    # error between the probability of a multi-contact event in a
    # given gridpoint and the predicted occupancy of the polymer on
    # this gridpoint and dividing the result by the total amount of
    # multi-contact events. This quantifies the difference between
    # the epitopsy prediction and the observed MD data for a
    # gridpoint while also accounting for how well in general the
    # polymer binds to the protein in MD simulations.
    # In this case, the lower the fitness the better
    def fitness(self, project, output_folder, job, protein_name, pose_file):
        try:
            box = epi.DXFile.read_dxfile(output_folder
                    + protein_name + '_occ.dx', 'esp')
            with open(output_folder + protein_name 
                    + '_box_coords_job' + str(job) + '.yml') as f_:
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
                        e = mdl.enviroment.Enviroment(output_folder
                                + protein_name + '.pdb', output_folder
                                + pose_file)
                        a = mdl.analysis.Analytics(e)
                        probability = [float(freq)/(len(a.df_dist)
                            * int(project.jobs[0].meta['sampling_rate']))
                            for freq in frequency]
            return sk.metrics.mean_absolute_error(probability,
                                occupancy) / sum(probability)
        except ValueError:
            return None

    # Compute performance of the individuals.
    def compPerformance(self, population, analysis_folder, generation_pls):
        pop_perf = {}
        for j in enumerate(population):
            ouput_folder = analysis_folder + j[1].name + '/'
            job =  j[0]
            # This could be a problem with the 1A2C / 1lya mixup
            protein_name = generation_pls.jobs[job].meta['protein']
            pose_file = generation_pls.jobs[job][0].pymol.pose.file_name(
                    molecule = 'polymer', state = 'full')
            pop_perf[j[1].name] = self.fitness(generation_pls, ouput_folder,
                    job, protein_name, pose_file)
        # Remove the jobs where no multi-contact event took place.
        for k,v in pop_perf.items():
            if v == None:
                del pop_perf[k]
        return sorted(pop_perf.items(), key = operator.itemgetter(1))

    # Select the indiviuals used to generate the next generation.
    def selePop(self, pop_sorted, best_sample, lucky_few):
        breeder = []
        for i in range(best_sample):
            breeder.append(pop_sorted[i][0])
        for i in range(lucky_few):
            breeder.append(random.choice(pop_sorted)[0])
            random.shuffle(breeder)
        return breeder

    # Load the polymer sequences (the 'genome') of each individual.
    def load_configs(self, configs):
        sequences = {}
        for i in range(len(configs)):
            with open(str(configs[i])) as f_:
                foo = yaml.load(f_)
                sequences[configs[i].parent.name] = foo[
                        'sim_parameter']['poly_sequence']
        return sequences

    # Recombination function.
    def create_child(self, parent1, parent2):
        child = []
        for i in range(len(parent1)):
            if (int(100 * random.random()) < 50):
                child.append(parent1[i])
            else:
                child.append(parent2[i])
        return child

    # Create the next generation.
    def next_pop(self, sequences, breeder, nr_children):
        next_pop = []
        for i in range(len(breeder)/2):
            for j in range(nr_children):
                l = self.create_child(sequences[breeder[i]],
                        sequences[breeder[len(breeder) -1 -i]])
                next_pop.append(l)
        return next_pop

    # Mutate a monomer ('gene') according to the mutation probability,
    # the possible genes are given by the gene list.
    def mutate(self, next_population):
        for i in range(len(next_population)):
            for j in range(len(next_population[i])):
                if random.random() <= self.mutation_rate:
                    next_population[i][j] = Genetic_algorithm.gene_list[int(
                        random.random() * len(Genetic_algorithm.gene_list))]
        return next_population

    # Save the current best individual, it's generation,
    # fitness and sequence.
    def currBest(self, pop_perf, current_best, history):
        for i in range(len(pop_perf)):
            if pop_perf[i][1] < current_best[2]:
                current_best = ('generation' + str(Genetic_algorithm.generation),
                        pop_perf[i][0], pop_perf[i][1],
                        history['generation' + str(
                            Genetic_algorithm.generation)][pop_perf[i][0]]['sequence'])
        return current_best

    # Check the status of jobs running on the cluster.
    def jobchecker(self, project):
        job_state = []
        for job in project:
            job.state.update()
            job_state.append(job.state.is_finished)
        return job_state

    def run(self, nr_generations, project_folder, cluster_folder, best_sample,
            lucky_few, nr_children):
        history = col.OrderedDict()
        current_best = ('generationX', 'test', float('inf'), 'No sequence')
        while Genetic_algorithm.generation < nr_generations:
            print Genetic_algorithm.generation
            if Genetic_algorithm.generation != 0:
                while True:
                    if all(self.jobchecker(project)):
                        break
                    else:
                        print "not all MD simulations are finished, going to sleep for 10 min"
                        time.sleep(600)
                print "MD simulation finished, downloading and analysing Data"
                os.mkdir(project_folder + 'generation'
                        + str(Genetic_algorithm.generation))
            gen_folder = project_folder + 'generation' + str(
                    Genetic_algorithm.generation)
            md_folder = project_folder + 'generation' + str(
                    Genetic_algorithm.generation) + '/MD/'
            ana_folder = project_folder + 'generation' + str(
                    Genetic_algorithm.generation) + '/analysis/'
            os.mkdir(ana_folder)
            if Genetic_algorithm.generation != 0:
                project.retrieve(gen_folder)
            # Population = list(pl.Path(md_folder + 'jobs/').glob('*'))
            configs = list(pl.Path(md_folder + 'jobs/').glob(
                '*/config_with_setup.yml'))

            # Getting our pls data and generating our population.
            generation_pls = pls.Analysis.Project(md_folder)
            population = [pl.Path(md_folder + 'jobs/'
                + job.db_path.parent.name) for job in generation_pls.jobs]

            # Getting our MDL & epitopsy data
            print "Retrieved MD data, starting Epitopsy analysis"
            for individual in enumerate(population):
                os.mkdir(ana_folder + '%s' % individual[1].name)
                mdl.pipeline.Pipeline_run(md_folder, ana_folder + '%s/'
                        % individual[1].name, job=individual[0],
                        meshsize=self.meshsize)

            # Calculating the Fitness and creating the next Generation.
            print "Calculating Fitness and creating the next Generation"
            pop_perf = self.compPerformance(population, ana_folder,
                    generation_pls)
            breeder = self.selePop(pop_perf, best_sample, lucky_few)
            sequences = self.load_configs(configs)
            next_population = self.mutate(self.next_pop(sequences, breeder,
                nr_children))
            history['generation' + str(Genetic_algorithm.generation)] = {}
            for i in range(len(pop_perf)):
                history['generation' + str(
                    Genetic_algorithm.generation)][pop_perf[i][0]] = {}
                history['generation' + str(Genetic_algorithm.generation)][
                        pop_perf[i][0]]['fitness'] = pop_perf[i][1]
                history['generation' + str(Genetic_algorithm.generation)][
                        pop_perf[i][0]]['sequence'] = sequences[pop_perf[i][0]]
            current_best = self.currBest(pop_perf, current_best, history)
            print "Preparing MD simulation runs"

            # Path to slurm config.
            slurm_config = pl.Path('slurm.cfg')
            if not slurm_config.exists():
                print slurm_config, 'does not exist'
            config_folder = pl.Path(md_folder + 'static/')

            # Create dict with all static files needed to run the simulation
            statics = {'config': config_folder.joinpath('parameters.yml'),
                    'active_site': config_folder.joinpath('active_sites.h5'),
                    'protein': config_folder.joinpath(
                        generation_pls.jobs[0].meta['protein']),
                    'script': config_folder.joinpath('sampling2.py'),
                    #'ic50': config_folder.joinpath('ic50.h5'),
                    'surface_db': config_folder.joinpath(
                        'hydrophobic_parameters.h5'),
                    'protonation_db': config_folder.joinpath('protonation.h5')}

            # Check if everything is in order.
            for file_ in statics.values():
                if not file_.exists():
                    print file_, 'does not exist'
            statics = dict(zip(statics.keys(), map(str, statics.values())))
            config_defaults = pl.Path(md_folder + '/jobs/abcd/config.yml')
            if not config_defaults.exists():
                raise FileNotFoundError

            # Create project.
            project = LS.Project(cluster_folder + 'generation'
                    + str(Genetic_algorithm.generation) + '/MD/',
                    defaults=config_defaults.as_posix(),statics=statics,
                    slurm_parameters=slurm_config.as_posix())

            # Creates folders and copies static data to cluster.
            project.create_environment()

            # Initialize parameters that are used to vary the job settings
            # like polymer type, timesteps, etc.
            parameters = {}
            # Simulation parameters.
            sim_parameter = {}
            # Lammps parameters.
            lmp_parameter = {}
            # Bundling them up.
            parameters['lammps_parameter'] = lmp_parameter
            parameters['sim_parameter'] = sim_parameter

            for genome in next_population:
                parameters['sim_parameter']['poly_sequence'] = genome
                project.create_job(parameters)

            print "Starting MD simulations"
            Genetic_algorithm.update_gen()
            print Genetic_algorithm.generation
            project.run_jobs()

        f_ = open(project_folder + 'history.yml', 'w')
        yaml.dump(history, stream = f_)
        f_.close()

        g_ = open(project_folder + 'best_hit.yml', 'w')
        yaml.dump(current_best, stream = g_)
        g_.close()
