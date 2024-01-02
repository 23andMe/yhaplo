"""Specify package code and data for yhaplo."""

from setuptools import find_packages, setup

setup(
    package_data={
        "yhaplo.data.tree": ["*"],
        "yhaplo.data.variants": ["*"],
    },
    packages=find_packages(exclude=["tests*"]),
    url="https://github.com/23andMe/yhaplo",
)
