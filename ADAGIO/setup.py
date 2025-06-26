#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['networkx', 'sklearn', 'pandas', 'dill', 'numpy']

test_requirements = ['pytest>=3', ]

setup(
    author="Mert Erden",
    author_email='mert.erden@tufts.edu',
    python_requires='>=3.8',
    classifiers=[
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
    description="ADAGIO Ranking of Genes",
    entry_points={
        'console_scripts': [
            't_map=t_map.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='t_map',
    name='t_map',
    packages=find_packages(include=['t_map', 't_map.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
