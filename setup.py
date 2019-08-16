#!/usr/bin/env python

import setuptools
from distutils.core import setup

if __name__ == '__main__':
    setup(
        name="assess_workflows",
        version="0.1",
        descriptions="Workflow execution for assess",
        author="Eileen Kuehn",
        author_email="eileen.kuehn@kit.edu",
        url="https://github.com/eileen-kuehn/assess_workflows",
        packages=setuptools.find_packages(),
        dependency_links=[],
        install_requires=[],
    )
