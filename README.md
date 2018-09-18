# Annotate_XL
OOP Python project to provide indepepndent annotations of a cross-link:spectrum match


Annotate_XL creates annotated spectra for a cross-link spectrum match. 
For a given cross-link identification and peak list containing m/z and intensity values Annotate_XL creates all possible fragment ions and matches them to observed ions with a mass that is within 10ppm of the calculated theoretical mass. 

Annotate_XL requires a cross-link ID of the form: "αSequence-βSequence-an-bn" where n represents the position of the cross-linker.
e.g. DTHKSEIAHR-FKDLGEEHFK-a4-b2

It also requires a deconvoluted CSV file containing the peak list in the following format:
m/z, intensity,
84.0838012695312,8.83298371701586,
86.0891036987305,6.32303576257556,
103.05290222168,2.62033348992199,
110.071502685547,26.9497413348613,
120.080596923828,22.4569751705114,
…

Please see test_peak_list.csv for an example.

Upon first execution:
- Open Annotate_xl.py and change the OBSERVED_BASE_DIR to the location of your Annotate_XL.py download and save.

To execute the code: 
- Place your peak_list.csv into the same folder as Annotate_XL.py. 
- Navigate to the Annotate_xl directory in your terminal. 
- Type python Annotate_XL.py cross-link-id peak_list.csv into your terminal.

In the same folder as your peak list Annotate_XL will generate a CSV file containing all annotations named “cross-link-id_annotatexl.csv and a print quality (300 dpi) PNG named “cross-link-id.png”

Annotate_XL comes with a test file to get you started. To run the example files:
- Ensure you have opened Annotate_xl.py and changed the OBSERVED_BASE_DIR to the location of your Annotate_XL.py download and save.
- Navigate to the Annotate_xl directory in your terminal. 
- Type python Annotate_XL.py DTHKSEIAHR-FKDLGEEHFK-a4-b2 test_peak_list.
- Two files are generated in you Annotate_XL directory: DTHKSEIAHR-FKDLGEEHFK-a4-b2_annotatexl.csv and DTHKSEIAHR-FKDLGEEHFK-a4-b2.png 

Annotate_XL has currently been tested on Linux/Unix operating systems. Stay tuned for the development of a future web portal…


===========================================================================
Author: Juliette M.B James March 2018
This code is offered under a GNU GPLv3 License. For more details on 
licensing terms please see: https://choosealicense.com/licenses/gpl-3.0/
