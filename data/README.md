# Ebola data
Epidemiological data from the Ebola epidemic in West Africa between 2013 and 2016.


### Data collection

These data are the result of the effort of many, many people.

- Raw data on country-level mortality rates were obtained from https://github.com/cmrivers/ebola/blob/master/country_timeseries.csv ;
- Raw data on health centres for Liberia and Guinea were obtained from http://ebola-data.un.org/Liberia/Liberia_WADC00086_OHDR_Health_Facilities_UNDP_LIBGov%202007/ and http://ebola-data.un.org/Guinea/GN_WADC00420_HXL_Health_Facilities/ ;
- Raw data on sub-national cases and mortality was obtained from [HDX](https://data.humdata.org/dataset/rowca-ebola-cases);
- Data from tables were extracted from PDFs using [Tabula](http://tabula.technology/);
- Raw data on Ebola treatment units (ETU) were obtained from https://data.humdata.org/dataset/ebola-treatment-centers .

### Data processing

All processing is done in [R](https://www.r-project.org/). The dependencies are [**tidyverse**](https://www.tidyverse.org/), [**maptools**](https://rdrr.io/rforge/maptools/), [**rgdal**](https://cran.r-project.org/web/packages/rgdal/index.html) and [**spdep**](https://github.com/r-spatial/spdep).
