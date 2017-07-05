# Approximate Georeferencing

Social scientists have rapidly seized on the research potential of GIS data. However, the analyses that spatial data afford almost always require some amount of record linkage to join quantities of interest to the relevant spatial covariates. In locations with a nascent spatial data industry, researchers’ typical problems managing messy data are compounded with the difficulties of working with imprecise or unstandardized GIS files.

We outline methods of record linkage specialized to the context of spatial data. We use these techniques to assign geographic coordinates to all Kenyan polling stations from 1997 to 2013. Through a generalizable example of this automated geocoding process, we provide what we hope to be a helpful introduction to rigorous, systematic record linkage in the world of messy spatial data.

To read the main document, please download or clone this repository to your local machine and open `code/doc.html` in a web browser. 

[//]: # "You may also view the document [here](http://htmlpreview.github.com/?https://github.com/coletl/georef/blob/master/code/doc.html)."

To reproduce our results, you should
1) download/clone this repository, 
2) download the [replication data](https://drive.google.com/drive/folders/0B8K1PQKTPN42bS1SUXpvSUNvNUU?usp=sharing), and 
3) place the data files in a directory named `data` in the main project directory.

Opening `georef.Rproj` in RStudio will automatically set your working directory, so you can run the code without editing file paths.

Please direct correspondence to Cole Tanigawa-Lau at <coletl@nyu.edu>.
