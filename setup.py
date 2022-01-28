from setuptools import setup, find_packages

setup(
    name="icolos",
    maintainer="Christian Margreitter, Harry Moore",
    version="1.4.0",
    packages=find_packages("."),
    include_package_data=True,
    package_dir={"config": "icolos/config"},
    package_data={"icolos": ["config/logging/*.json"]},
    description="Icolos Workflow Manager",
    entry_points="""
		[console_scripts]
		icolos=icolos.scripts.cli:entry_point
  	""",
    python_requires=">=3.8",
)
