from setuptools import setup, find_namespace_packages
from pathlib import Path

long_description = (Path(__file__).parent / "README.md").read_text()

setup(
    name="icolos",
    maintainer="Christian Margreitter, Harry Moore",
    version="1.9.0",
    url="https://github.com/MolecularAI/Icolos",
    packages=find_namespace_packages(where="src"),
    package_dir={"": "src"},
    description="Icolos Workflow Manager",
    long_description=long_description,
    long_description_content_type="text/markdown",
    python_requires=">=3.8",
    entry_points={
        "console_scripts": [
            "icolos = icolos.scripts.executor:main",
            "validator = icolos.scripts.validator:main",
            "sdf2smi = icolos.scripts.sdf2smi:main",
        ]
    },
)
