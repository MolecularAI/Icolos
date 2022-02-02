from setuptools import setup, find_namespace_packages

setup(
    name="icolos",
    maintainer="Christian Margreitter, Harry Moore",
    version="1.5.0",
    url="https://github.com/MolecularAI/Icolos",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    description="Icolos Workflow Manager",
    python_requires=">=3.8",
    entry_points={"console_scripts": ["icolos = icolos.scripts.executor:main"]},
)
