#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

requirements = ['pyopenms', 'pandas', 'matplotlib',
                'typing', 'requests', 'sqlalchemy',
                'pathlib', 'sqlalchemy_utils', 'click']

test_requirements = ['pytest>=3', ]

setup(
    author="Group 5",
    author_email='dingmingliu@outlook.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
    ],
    description="Group 5 Spectra Package",
    entry_points={
        'console_scripts': [
            'project_spectra=project_spectra.cli:main',
        ],
    },
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='project_spectra',
    name='project_spectra',
    packages=find_packages(include=['project_spectra', 'project_spectra.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.1.0',
    zip_safe=False,
)
