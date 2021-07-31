from setuptools import setup, find_packages

setup(
    name="spatial-clustering",
    version="0.1",
    author="Alison Peard",
    author_email="alison.peard@gmail.com",
    description="Spatial clustering package",
    url="https://github.com/alisonpeard/spatial-clustering",

    packages=find_packages(),
    install_requires=[
        "networkx",
        "numpy",
        "scipy",
        "sklearn",
        "matplotlib",
        "pyperclip"
    ],

    scripts=['bin/community']
)
