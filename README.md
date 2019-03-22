# shale_microcrack_growth
Matlab scripts that simulate microcrack growth in the shale during kerogen thermal maturation. The formulation details were discussed in [Yang and Mavko (2018)](https://doi.org/10.1306/05111817062) and [Yang (2019)](https://searchworks.stanford.edu/view/13122311) Ch.2. Two main files run the geological and laboratory settings separately.

When running `main_lab.m` and `main_geo.m`, one can set the parameter `var` to 'n2V' (lab), 'D' (geo), 'R', or 'KIc' to test the effect of the volume expansion parameter, pore size or fracture toughness, respectively. One can reproduce figures shown in the above references for each scenario. 
