import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
  long_description = fh.read()

setuptools.setup(
  name="focont", # Replace with your own username
  version="0.1",
  author="Okan Demir",
  author_email="demir@ee.bilkent.edu.tr",
  description="Fixed order controller design tool",
  long_description=long_description,
  long_description_content_type="text/markdown",
  url="https://github.com/dmrokan/focont",
  packages=setuptools.find_packages(),
  classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License",
    "Operating System :: OS Independent",
  ],
  python_requires='>=3.8',
)

