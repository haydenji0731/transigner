from setuptools import setup, Extension
import subprocess
import os

with open("README.md", "r") as f:
    long_description = f.read()

def get_em_dir():
    conda_prefix = os.getenv('CONDA_PREFIX')
    return os.path.join(conda_prefix, 'bin')

setup(
	name="transigner",
	version="1.1.3",
	author="Hyun Joo Ji",
	author_email="hji20@jh.edu",
	description="assign long RNA-seq reads to transcripts",
    long_description=long_description,
    long_description_content_type="text/markdown",
	url="https://github.com/haydenji0731/transigner",
	install_requires=[
        'pysam',
        'numpy',
        'pandas',
    ],
	python_requires='>=3.8',
	packages=['transigner'],
    data_files=[('bin', ['em'])],
	entry_points={'console_scripts': ['transigner = transigner.run_transigner:main'],},
)

# conda_prefix = os.getenv('CONDA_PREFIX')
# dest_dir = os.path.join(conda_prefix, 'bin')
# subprocess.call(f'cp em {dest_dir}', shell=True)
