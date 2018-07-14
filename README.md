The MDLatticeAnalysisTool (mdl) is a toolkit designed to identify instances of multi-contact-events (mce) in PolyLibScan (pls) datasets. A mce is a snapshot of a pls run where more than one bead is in close vicinity of the target protein (the standard threshold used is 4 A). The pipeline allows comparisons between the mdl results and occupancy predictions made by Epitopsy for the corresponding Polymer / Protein interactions. It also contains a genetic algorithm making use of the mdl pipeline to optimize polymer composition in respect to the amount of mces and the overlap between predicted occupancies and observed probabilities of mces in the pls simulations.

# Dependencies
- PolyLibScan
- Epitopsy
- pathlib2
- sklearn
- yaml
- pandas
- matplotlib
- Bio.PDB

# Genetic algorithm specific:
- LammpsSubmit

# MDLatticeAnalysisTool (mdl)
The mdl can be used to identify and quantify mces of pls runs. A example use case is shown below:

	import MDlatticeAnalysisTool as mdl

Create a Environment to parse pdb data and create a box around the protein

	env = mdl.environment.Environment('path/to/pls/project', 'path/to/polymer_poses.pdb')

Calculate the distance to the nearest protein for each monomer, this is done for each
snapshot of a single repetition (run) of a pls job. 

	ana = mdl.analysis.Analytics(env, distance_cutoff = 4)

Dataframe containing information for each snapshot on the monomer coordinate, coordinate of the
closest protein atom, distance (in A) and to which residue the protein atom belongs.

	ana.dataframe

Plot to show the probability of each monomer to bind to the protein per snapshot 
(binding ="single") / the probability of a monomer to be part of a "multi-contact" snapshot
(binding="multi")
		
	ana.probability_distribution_polymer(binding="multi")

Plot to show the probability of each protein residue to bind to a monomer per snapshot 
(binding="single") / the probability of a residue to be part of a "multi-contact" snapshot (binding="multi")

	ana.probability_distribution_protein(binding="multi")

Visualize the mces in Pymol (pymol has to be started with the -R flag)

	pym = mdl.Pymol.pymol(ana)

load the Protein into Pymol

	pym.setup()

Color protein residues according to their probability to bind to a monomer (binding='single') or to be part
in a mce (binding='multi') (between 1 and 5%: yellow; 5 to 10%:
orange; 10 to 20%: brown; over 20%: red)

	pym.show_bindings(binding = "multi")

# Pipeline

The pipeline incorporates pls analysis, mdl analysis and (optionally, but recommended) epitopsy analysis. It will automatically create all the files it needs from the pls project and will analyze all samples/repetitions/runs of a pls job. 
Input: 
	   
	    Pls project folder 
	    (empty) output folder
	    job id
Output: 
		 
		Pose pdb files of all runs
		.pdb and .pqr file of the protein
		.pqr file of the polymer
		Epitopsy esp, epi and occ files
		'protein'job'nr'.yml, containing the average mces for each residue
		'protein'_coords_job'nr'.yml, containing all coordinates where a mce took place
		'protein'_coords_dict_job'nr'.yml containing all coordinates where a mce took place and how often such an event took place
		'protein'_box_coords_job'nr'.yml containing all coordinates where a mce took place, translated to epitopsy gid coordinates
		a plot in .png format showing the probabilities of all monomers to be part of a mce

To start the pipeline the following commands are used:

	from MDlatticeAnalysisTool import pipeline
	pipeline.Pipeline_run('path/to/pls/project/', '/path/to/output/folder/', job=0, meshsize = [0.8, 0.8, 0.8])

job: Job nr of the pls job that you want to analyze

meshsize: Size of the Epitopsy mesh, smaller = higher resolution but longer calculation time


The output files can be further analyzed with the scripts rendering.pml and mean_absolute_error.py:

rendering.pml is a pymol script to visualize the epitopsy and mdl results in pymol, you have to adjust the path to the path of your 'protein'_job'nr'.yml in it and can then start it by typing in the terminal:

	
	pymol rendering.pml

mean_absolute_error.py is a python script to calculate the mean absolute error between the occupancies predicted by epitopsy and the mces of the pls simulation, it will create a yaml file called mean_absolute_error.yml, containg both the raw mae and the mae weighted by the sum of mce probabilities.

To use it, adjust:

	line 4 to the path of the folder containing the MDlatticeAnalysisTool
	line 10 to your pipeline output folder
	line 11 to your pls project folder
	and line 12 to the pls job nr

you can run it in the terminal with the command:

	python mean_absolute_error.py

# Genetic algorithm:
The genetic algorithm aims to improve the overlap between mce probabilities and epitopsy occupancies and the frequencies of mce by mimicking natural selection: Each monomer is seen as a 'gene' and each polymer as a 'individual', a group of polymers is seen as a 'population'. For each 'generation', the fitness for each indiviual is calculated by the mean absolute error between the the occupancies predicted by epitopsy and the mce of the pls simulation, weighted by the sum of the mce probabilities. The n individual with the highest fitness will be recombined for creating the next generation, but also some 'lucky few' randomly chosen from the whole population will be used for recombination, so that local minima are of less impact. Also, the genes can randomly mutate (being switched with another monomer) according to a predefined mutation rate. The algorithm switches between running pls simulations on the cluster and analyzing the data locally. A example use case is shown below:

	import MDlatticeAnalysisTool as mdl

	mdl.genetic_algorithm.Genetic_algorithm('/local/project/folder/', '/cluster/project/folder/', best_sample = 4, lucky_few = 2, nr_children = 4, nr_generations = 11, mutation_rate = 0.02, meshsize = [4,4,4])
local project folder: The folder containing your initial parent generation, the hierarchy has to be project_folder/generation0/MD/jobs/''parents'', the configuration file is taken from a parent called 'abcd', which needs to be in the jobs folder. the MD folder also has to contain a static folder, which can be copied from the pls project


cluster project folder: the project folder on the cluster, has to be empty


best_sample: The best n individuals used to create the next generation


lucky_few: n randomly chosen individuals from the whole population that will be used for recombination


nr_children: Number of children created from each parent pair


nr_generations: Number of generations the algorithm will run and analyse


mutation rate: probability for a gene to mutate, per gene per generation

meshsize: The size of the Epitopsy mesh, has to be adjusted depending on the available ressources

# output:
Each individual will have its own analysis folder, containing the same output as described in the Pipeline section. Furthermore, the local project folder will contain a 'best_hit.yml' file, containing the generation, jobname, fitness and sequence of the individual with the highest fitness. It will also contain a 'history.yml' file, containing the same information for all individuals, which allows for example plotting the mean difference in fitness between generations.