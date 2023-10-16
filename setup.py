from setuptools import setup

setup(
	name="snoopsv",
	version="0.0.1",
	description="A software to detect structural variation signatures"
				" from long reads (PacBio/ONT)",
	author="Milad Mortazavi",
	author_email="smmortazavi@ucsd.edu",
	url="https://github.com/mldmort/snoopSV",
	packages=["snoopsv"],
	install_requires=[
		'pysam', 'pandas', 'numpy', 'pytest'
	],
	entry_points={
		"console_scripts": [
			"snoopsv = snoopsv.__main__:main"
		],
	},
	classifiers=[
		"License :: OSI Approved :: MIT License",
		"Programming Language :: Python :: 3",
	],
)
