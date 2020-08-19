# phage-antibiotics

This repositery contains all code written for "The Impact of Pf Bacteriophages on the Fitness of Pseudomonas aeruginosa: A Mathematical Modeling Approach" 

### There are three main files, which run the analysis and generate all figures in the paper and Supplementary Material:

  'main_graph.m' - generates figures 2-5 in the main text
  
  'sensitivity.m' - generates figures S1 and S2 in the Supplementary Material
  
  'sensitivity_seq.m' - generates figure S3 in the Supplementary Material
  
### In order to run the analysis, the additional files are required in the same folder:

  'growthRate.m' 
  
  'compareRegimens.m'
  
  'EmaxcompareNoResist.m'
  
  'competeRegimens.m'
  
  'EmaxcompeteNoResist.m' 
  
 ### For proper figure formatting, the following scripts should be downloaded from File Exchange:
 
  'tight_subplot.m' -  Pekka Kumpulainen (2020). tight_subplot(Nh, Nw, gap, marg_h, marg_w) (https://www.mathworks.com/matlabcentral/fileexchange/27991-tight_subplot-nh-nw-gap-marg_h-marg_w), MATLAB Central File Exchange. Retrieved August 19, 2020. 
  
  'format_ticks.m' - Alexander Hayes (2020). Format Tick Labels (https://www.mathworks.com/matlabcentral/fileexchange/15986-format-tick-labels), MATLAB Central File Exchange. Retrieved August 19, 2020. 
  
  'boxplotGroup.m' -  Adam Danz (2020). boxplotGroup (https://www.mathworks.com/matlabcentral/fileexchange/74437-boxplotgroup), MATLAB Central File Exchange. Retrieved August 19, 2020. 
