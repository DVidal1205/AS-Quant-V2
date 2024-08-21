AS-Quant generates the quantification level of different types of alternative splicing events between two biological conditions and assesses the significant of the events.

# CHANGES IN VERSION 2.0

## Generation of AS Events

To ensure the most accurate and dated results, the 7 alternative splicing events are now procedurally generated from a gtf file. Methods for generation are in methods/preprocess.py

## Paralellization

1. AS-Quant now supports parallel processing to speed up the quantification process. Users can specify the number of **cores** to use for parallel processing with the **cores** flag in config.ini. If no cores are specified, then the default is to run the code sequentially. Users can also use the MAX value to use all available cores on the machine.

## Change in Samtools

Debian distributions have issues with easily configuring execution permissions for the samtools binary. To address this issue, I have installed the samtools binary to my system in the /usr/local/bin directory. To accomodate this change, I had to modify samtools_directory in the preprocess.py file to /usr/local/bin/samtools. This change will allow the program to run on Debian distributions without any issues. (Note: the version of samtools downloaded on my system is identical to the one in the repository. I am looking into a more permanent solution for this issue to add in the documentation.)

# Installation

AS-Quant tool can be downloaded/cloned directly from github. Users need to have python3 installed on their machine. It can work both on Windows and Linux platform.

Users will need to run the as-quant.py file in order to demonstrate the AS events.
$python3 as_quant.py

To generate plot for a specific exon, users need to run make_plots.py and provide input in a specific pattern: chrom:geneName:start-end
$python3 make_plots.py

# User manual

User manual for AS-Quant is available on github https://github.com/compbiolabucf/AS-Quant/blob/main/AS-Quant_User_Manual.pdf. Note that a -c flag has been added to specify the number of cores to use for parallel processing. The default number of cores is the number of cores available on the machine.

# Sample data

AS-Quant example data for mouse(mouse mm10) is included in the folder https://github.com/compbiolabucf/AS-Quant/tree/main/sample_input_mouse. To run from scratch, users can use the sample data for quantization of AS events.

# Citation

Please use the following information to cite.

MDPI and ACS Style
Fahmi, N.A.; Nassereddeen, H.; Chang, J.; Park, M.; Yeh, H.; Sun, J.; Fan, D.; Yong, J.; Zhang, W. AS-Quant: Detection and Visualization of Alternative Splicing Events with RNA-seq Data. Int. J. Mol. Sci. 2021, 22, 4468. https://doi.org/10.3390/ijms22094468

AMA Style
Fahmi NA, Nassereddeen H, Chang J, Park M, Yeh H, Sun J, Fan D, Yong J, Zhang W. AS-Quant: Detection and Visualization of Alternative Splicing Events with RNA-seq Data. International Journal of Molecular Sciences. 2021; 22(9):4468. https://doi.org/10.3390/ijms22094468

Chicago/Turabian Style
Fahmi, Naima Ahmed, Heba Nassereddeen, Jaewoong Chang, Meeyeon Park, Hsinsung Yeh, Jiao Sun, Deliang Fan, Jeongsik Yong, and Wei Zhang. 2021. "AS-Quant: Detection and Visualization of Alternative Splicing Events with RNA-seq Data" International Journal of Molecular Sciences 22, no. 9: 4468. https://doi.org/10.3390/ijms22094468

# Contact the Author

Naima Ahmed Fahmi: fnaima@knights.ucf.edu
Wei Zhang: wzhang.cs@ucf.edu
