import io
import os
import re

from setuptools import find_packages
from setuptools import setup


def read(filename):
    filename = os.path.join(os.path.dirname(__file__), filename)
    text_type = type(u"")
    with io.open(filename, mode="r", encoding='utf-8') as fd:
        return re.sub(text_type(r':[a-z]+:`~?(.*?)`'), text_type(r'``\1``'), fd.read())


setup(
    name="wcb_outflow",
    version="0.1.0",
    url="",
    license='MIT',

    description="Circulation in warm-conveyor belts",
    long_description=read("README.rst"),

    packages=find_packages(include=["wcb_outflow"], exclude=('tests',)),

    install_requires=[
        "numpy",
        "matplotlib",
        "scipy",
        "scitools-iris>=3.0",
        "cartopy",
        "pylagranto",
    ],

    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python',
        'Programming Language :: Python :: 3',
    ],
)
