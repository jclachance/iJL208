from setuptools import setup, find_packages

setup(name='iJL208',
      version='0.1.0',
      download_url = 'https://github.com/peterldowns/mypackage/archive/0.1.0.tar.gz',
      description='Bundle of code and data for the metabolic model of Mesoplasma florum',
      url='https://github.com/jclachance/iJL208/',
      author='Jean-Christophe Lachance',
      author_email='jean-christophe.lachance@usherbrooke.ca',
      license='MIT',
      packages=find_packages(),
      package_data={'iJL208':['data/*.*']},
      include_package_data=True,
      install_requires=['cobra>=0.11.0','pandas','numpy>=1.13',
			'BioPython','seaborn',
			'matplotlib'])
