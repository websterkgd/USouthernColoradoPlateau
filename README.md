# USouthernColoradoPlateau
===================

An analysis of Uranium in unregulated water sources from the southern Colorado Plateau. This repo contains unique code and files that are not hosted elsewhere. 

### Repo Contents

* **USCP_Data:** - Data files used in this research
	* AquifersColoradoPlateau.JPEG - A figure from Robson and Banta 1995
	* CoCoAq_0-1.KML - A KML file of the portion of the Coconino-De Chelly aquifer with 0 - 1k mg/L TDS
  	* CoCoAq_1-3.KML - A KML file of the portion of the Coconino-De Chelly aquifer with 1 - 3k mg/L TDS 
 	* CoCoAq_3-25.KML - A KML file of the portion of the Coconino-De Chelly aquifer with 3 - 25k mg/L TDS
 	* CoCoAq_g25.KML - A KML file of the portion of the Coconino-De Chelly aquifer with > 25k mg/L TDS
 	* CoCoAqSurface.KML - A KML file of the southern surface portion of the Coconino-De Chelly aquifer  
  	* CoconinoDeChellyTDS_Cropped.JPEG - A cropped figure of the TDS in the Coconino-DeChelly aquifer from Robson and Banta 1995
  	* CoconinoDeChellyTDS.JPEG - A figure of the TDS in the Coconino-DeChelly aquifer from Robson and Banta 1995
	* F1_PanelAwoMarkings.JPG - Satellite imagery for panel 1A
	* F1_PanelBwoMarkings.JPG - Satellite imagery for figures 1B and 3
	* Mspat.KML - KML file of mines in the study region
	* Mspatelev.txt - a txt file of mines in the region with elevation data
	* Ucaus.csv - a csv file of data extracted from various sources
	* Uspat.KML - a KML file of water source locations in this study
	* UWellList.csv - a csv file of water sources in the study area
  
* **USCP_Rcode:** - Analytical and figure generating code
	* CreateCoconino_20211017.r - Creates a kml file of the southern surface portion of the Coconino-De Chelly aquifer
	* CreateCoconinoTDS_0-1k_20211022.r - Creates a kml file of the portion of the Coconino-De Chelly aquifer with 0 - 1k mg/L TDS
	* CreateCoconinoTDS_1-3k_20211022.r - Creates a kml file of the portion of the Coconino-De Chelly aquifer with 1 - 3k mg/L TDS
	* CreateCoconinoTDS_g3k_20211022.r - Creates a kml files of the portions of the Coconino-De Chelly aquifer with > 3k mg/L TDS
	* CreateMinesKML_20220707.r - Creates a kml file of mines within the study region
	* CreateUKML_20220707.r - Creates a kml file of water sources in the study region
	* Figure1_20220706.r - Creates Figure 1
	* Figure2_20220706.r - Creates Figure 2
	* Figure3_20220707.r - Creates Figure 3
	* Ucausal_20220712.r - Analysis of points within and outside of the Chinle Fm. Also generates Figure 6
	* UMultiRegression.r - Multiregression analysis of U in waters of the Colorado Plateau. Also generates Figures 4 and 5
	* UvsCoconino_20220718.r - Analysis of U concentration vs varying regions of the Coconino-De Chelly aquifer
	* UvsMines_20220707.r - Analysis of U concentration vs the presence of mines upstream
	* WaterUvsChinle_20220707.r - Analysis of U concentration vs the presence of the Chinle Formation
	* WaterUvsMorrison_20220707.r - Analysis of U concentration vs the presence of the Morrison Formation

## Contributors

[Dr. Kevin Webster](https://websterkgd.com/): Assistant Professor of STEM, [Diné College](https://www.dinecollege.edu/academics/meet-our-faculty-stem/); Associate Research Scientist, [Planetary Science Institute](https://www.psi.edu/about/staffpage/webster).
