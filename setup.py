from setuptools import setup,find_packages
setup(
    name="Py-Paired-Tag",
    version="1.0",
    description="Python version of Paired-Tag data processing pipeline",
    author="Sihao Huang",
    packages=find_packages(),
    include_package_data=True,
    scripts=["bin/ppt"],
    license="GPL 2.0"
)