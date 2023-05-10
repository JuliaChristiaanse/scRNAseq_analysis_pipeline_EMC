Please make sure you have python 3.10.9 installed on your local machine before you follow these steps:

Can be ran with the following code:
- (optional) create a new folder
- run git clone "link of this project"
- (important) move into the folder you just cloned
- run py -m pipenv install	(or py -m pipenv sync if you ran pipenv shell first)
- now 2 things will happen:
	1. a pipenv virtual environment will be created with an abstract name and a path
		to that pipenv. remember the name
	2. All dependencies will be installed from the pipfile.lock. This will take longer
		depending on the dependencies. mine is += 4 minutes.
- select the pipenv as python interpreter and kernal for jupyter
- All packages should be installed
- Run code accordingly