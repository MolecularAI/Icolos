from setuptools import setup

setup(
    name="icolos",
    maintainer="Christian Margreitter, Harry Moore",
    version="1.4.0",
    url="https://github.com/MolecularAI/Icolos",
    packages=["icolos"],
    package_dir={"": "src"},
    # include_package_data=True,
    # package_dir={"config": "icolos/config"},
    # package_data={"icolos": ["src/icolos/config/logging/*.json"]},
    description="Icolos Workflow Manager",
    python_requires=">=3.8",
    entry_points={"console_scripts": ["icolos = icolos.scripts.executor:main"]},
)
