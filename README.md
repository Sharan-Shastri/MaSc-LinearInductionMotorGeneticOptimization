# MaSc-LinearInductionMotorGeneticOptimization
Title: **_2D Hybrid Magnetic Field Model Performance Optimization for Linear Induction Motors_**

This repository includes most of the relevant Python code for my Ma.Sc. motor simulation and optimization software.
Due to intellectual property, not all the code can be released.
An attached [PDF](https://github.com/MichaelThamm/MaSc-LinearInductionMotorGeneticOptimization/blob/main/ProjectExplanation_GitHub.pdf) file is provided in this repo that highlights the top-down foundation of the coding project.
Information on the code containing intellectual property can only be found in this .PPTX file

Python modules required to run the code:

* math [documentation](https://docs.python.org/3/library/math.html)
* cmath [documentation](https://docs.python.org/3/library/cmath.html)
* timeit [documentation](https://docs.python.org/3/library/timeit.html)
* contextlib [documentation](https://docs.python.org/3/library/contextlib.html)
* collections [documentation](https://docs.python.org/3/library/collections.html)
* tkinter [documentation](https://docs.python.org/3/library/tk.html)
* json [documentation](https://docs.python.org/3/library/json.html)
* numpy [documentation](https://numpy.org/doc/)
* matplotlib [documentation](https://matplotlib.org/)
* inquirer [documentation](https://python-inquirer.readthedocs.io/en/latest/)

Note: [GitHub__LIM_ShowFromJSON.py](https://github.com/MichaelThamm/MaSc-LinearInductionMotorGeneticOptimization/blob/main/GitHub__LIM_ShowFromJSON.py) is the main file that calls the other .py files in this order:

* [GitHub__LIM_SlotPoleCalculation.py](https://github.com/MichaelThamm/MaSc-LinearInductionMotorGeneticOptimization/blob/main/GitHub__LIM_SlotPoleCalculation.py)
* [GitHub__LIM_Grid.py](https://github.com/MichaelThamm/MaSc-LinearInductionMotorGeneticOptimization/blob/main/GitHub__LIM_Grid.py)
* [GitHub__LIM_Show.py](https://github.com/MichaelThamm/MaSc-LinearInductionMotorGeneticOptimization/blob/main/GitHub__LIM_Show.py)

Make sure to run GitHub__LIM_ShowFromJSON.py from a terminal or the user input functionality will not work
