import setuptools

with open('README_PIP.md', 'r') as fh:
    long_description = fh.read()

setuptools.setup(
    name='disteval',
    version='0.3',
    py_modules=['disteval'],
    author='ba-lab',
    author_email='adhikarib@umsl.edu',
    description='DISTEVAL: For inter-residue protein distance evaluation',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/ba-lab/disteval',
    packages=setuptools.find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License'
    ],
    install_requires=[
        "numpy", "scikit-learn"
    ],
    entry_points = {
        'console_scripts': [
            'disteval = disteval:main'
        ]
    }
)