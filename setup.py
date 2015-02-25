import os
import sys
from setuptools import setup

version_py = os.path.join(os.path.dirname(__file__), 'endseq', 'version.py')
version = open(version_py).read().strip().split('=')[-1].replace('"','')
long_description = """
trans-seq - analysis pipeline for trans-splicing of lola or mod(mdg4)
"""

setup(
        name="trans-seq",
        version=version,
        install_requires=['matplotlib >1.3.0',
                          'pandas >0.14.0',
                          'numpy > 1.8.0',
                          'HTSeq',
                          'seaborn',
                          'brewer2mpl'],
        requires = ['python (>=2.7, <3.0)'],
        packages=['trans-seq',
                  'trans-seq.scripts'],
        author="Mohan Bolisetty",
        description='A toolset for working with endseq data',
        long_description=long_description,
        url="",
        package_dir = {'trans-seq': "trans-seq"},
        package_data = {'trans-seq': []},
        zip_safe = False,
        include_package_data=True,
##        scripts = ['endseq/scripts/endseq-script'],
        entry_points = {
            'console_scripts' : [
                 'trans-seq = trans-seq.trans-seq_main:main', 
            ],
        },  
        author_email="mohanbolisetty@gmail.com",
        classifiers=[
            'Development Status :: 4 - Beta',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License (GPL)',
            'Topic :: Scientific/Engineering :: Bio-Informatics']
    )
