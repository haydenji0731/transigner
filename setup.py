import setuptools


setuptools.setup(
	name="TranSigner",
	version="1.1.0",
	author="Hyun Joo Ji",
	author_email="hji20@jh.edu",
	description="Long read-to-transcript assignment creator",
	url="https://github.com/haydenji0731/transigner",
	install_requires=["pysam>=0.22.1", "tqdm>=4.66.5", "pyfastx>=2.1.0"],
	python_requires='>=3.8',
	packages=['transigner'],
	entry_points={'console_scripts': ['transigner = liftoff.run_liftoff:main'], },
)