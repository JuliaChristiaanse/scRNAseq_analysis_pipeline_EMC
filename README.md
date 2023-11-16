Please make sure you have python 3.10.9 installed on your local machine before you follow these steps:

Can be ran with the following code:
- Create a project directory
- Clone this folder to your directory
- (important) move into the folder you just cloned
- run py -m pipenv install	(or py -m pipenv sync if you ran pipenv shell first)
- now 2 things will happen:
	1. a pipenv virtual environment will be created with an abstract name and a path
		to that pipenv. remember the name
	2. All dependencies will be installed from the pipfile.lock. This will take longer
		depending on the dependencies. IMPORTANT: this might take a while!
- select the pipenv as python interpreter and kernal for jupyter
- Install R (v 4.3.0)
- Install these libraries:
	BiocManager				1.30.22
	Seurat					4.3.0.1
	SingleR					2.2.0
	Scuttle					1.10.2
	Scran					1.28.2
	Reticulate				1.31
	Seuratdisk				0.0.0.9020
	Seuratdata				0.2.2
	Patchwork				1.1.3
	Sceasy					0.0.7
	SingleCellExperiment	1.22.0	
	Seuratwrappers			0.3.1
	Hdf5r					1.3.8
	Monocle3				1.3.4
	Tidyverse				2.0.0
	EnhancedVolcano			1.18.0
	Ggplot2					3.4.4
	Dplyr					1.1.2
- You now have all required installments
- Sometimes you need to make a settings.json file:
	{
    "python.analysis.extraPaths": ["./Change name to your working directory"]
	}
- Run code accordingly

Marker genes format:
- an example on how to format your markergenes can be found in markergenes.txt.
- The code accounts for a global marker genes "group" and a bunch of marker genes all seperated by a comma.

Happy analyzing :D