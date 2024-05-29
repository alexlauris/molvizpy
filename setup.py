from setuptools import setup, find_packages
import codecs
import os

VERSION = '0.0.1'
DESCRIPTION = 'Molvizpy'
LONG_DESCRIPTION = ""

# Setting up
setup(
    name="molvizpy",
    version=VERSION,
    author="molvizpy",
    author_email="<alexandre.lauris@epfl.ch>",
    description=DESCRIPTION,
    long_description=LONG_DESCRIPTION,
    packages=find_packages(),
    install_requires=["streamlit"],
    keywords=['python', 'chemistry', 'visualization'],
    entry_points={
        "console_scripts": [
            "molvizpy_app = molvizpy.app:main",
        ],
    },
)