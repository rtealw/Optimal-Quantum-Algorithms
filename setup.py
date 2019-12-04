import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="quantum-query-optimizer",
    version="0.1.2",
    author="R. Teal Witter & Michael Czekanski",
    author_email="rtealwitter@gmail.com, michaeltczekanski@gmail.com",
    description="A toolkit to find the optimal quantum query complexity and query optimal quantum algorithm of arbitrary Boolean functions.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rtealw/QuantumQueryOptimizer",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ], 
    python_requires='>=3.6',
)
