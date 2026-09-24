from setuptools import setup, find_packages

setup(
    name="pepdesign",
    version="0.1.0",
    description="MEGA27-09 in-silico peptide design suite (CNN+GNN, real data, benchmarked)",
    packages=find_packages(include=["pepdesign*"]),
    python_requires=">=3.9",
)
