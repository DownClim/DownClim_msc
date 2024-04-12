# Outputs
Apr 12, 2024

All outputs and intermediary files of the analyses:

- `bencmark.tsv`: aggregated benchmark data of the wokflow
- `tmf/*.tif`: tropical moist forests distribution per country (Ivory
  Coast, French Guiana, new Caldonia), product (CMIP6 or CORDEX) and
  experiment (SSP1-2.6, SSP5-8.5, RCP-2.6, RCP-8.5)

``` r
fs::dir_tree()
```

    .
    ├── README.md
    ├── README.qmd
    ├── README.rmarkdown
    ├── bencmark.tsv
    └── tmf
        ├── Côte-d'Ivoire_CMIP6_ssp126.tif
        ├── Côte-d'Ivoire_CMIP6_ssp585.tif
        ├── Côte-d'Ivoire_CORDEX_rcp26.tif
        ├── Côte-d'Ivoire_CORDEX_rcp85.tif
        ├── French-Guiana_CMIP6_ssp126.tif
        ├── French-Guiana_CMIP6_ssp585.tif
        ├── French-Guiana_CORDEX_rcp26.tif
        ├── French-Guiana_CORDEX_rcp85.tif
        ├── New-Caledonia_CMIP6_ssp126.tif
        ├── New-Caledonia_CMIP6_ssp585.tif
        ├── New-Caledonia_CORDEX_rcp26.tif
        └── New-Caledonia_CORDEX_rcp85.tif
