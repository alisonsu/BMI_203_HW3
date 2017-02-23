# BMI_203_HW3

[![Build
 Status](https://travis-ci.org/alisonsu/BMI_203_HW3.svg?branch=master)](https://travis-ci.org/alisonsu/BMI_203_HW3)

Skeleton for HW3


## structure

```
.
├── README.md
├── sequences
│   ...
├── hw2skeleton
│   ├── __init__.py
│   ├── __main__.py
│   ├── align.py
│   ├── io.py
│   └── utils.py
└── test
├── test_align.py
└── test_io.py
```

## usage

To use the package, first run

```
conda install --yes --file requirements.txt
```

to install all the dependencies in `requirements.txt`. In this homework,
the code for each part of the assignment can be accessed at the topmost 
branch of this repository.


## testing

Testing is as simple as running

```
python -m pytest
```

from the root directory of this project.


## contributors

Original design by Scott Pegg. Refactored and updated by Tamas Nagy.

