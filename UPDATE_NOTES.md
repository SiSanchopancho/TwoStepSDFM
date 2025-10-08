# Upate notes for version 0.1.0

* Added a function to compute nowcasts for one or more quarterly target variables via prediction averaging (see README for more details).
* Added classes to the outputs of each function and respective generic plot and print methods.
* Added parameter checks for each function.
* ``twoStepSDFM()`` can now take zoo/ts objects as input for the data matrix.
* ``simFM()`` can now produce zoo/ts objects for the data, factors, and measurement errors.
* Several minor bug fixes.
* Structural project clean up.
