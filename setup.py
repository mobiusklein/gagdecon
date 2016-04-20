from setuptools import setup, find_packages



required = []
with open('requirements.txt') as f:
    required = f.read()

setup(
    name="gagdecon",
    version="0.1.3",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "gagdecon-ms1 = gagdecon.app:main"
        ]
    },
    zip_safe=False,
    install_requires=required,
    maintainer='Joshua Klein',
    maintainer_email="jaklein@bu.edu",
    description="Identify Glycosaminoglycans from raw mass spectral profiles"
    )
