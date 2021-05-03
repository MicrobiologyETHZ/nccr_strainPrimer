from setuptools import setup

requires = [
    'click',
    'biopython>=1.7',
    'pandas'

    ]

setup(
    name="strainPrimer",
    version="1.0",
    author="Anna Sintsova",
    author_email="ansintsova@ethz.ch",
    description=("Bioinformatic toolkit for barcode mapping and counting"),
    license="LICENSE",
    keywords="primers",
    install_requires=requires,
    packages=['strainPrimer'],
    entry_points={
        'console_scripts': ['strainPrimer=strainPrimer.main:main'],
    }
)
