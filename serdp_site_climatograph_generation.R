
# Load required libraries
library(terra)
library(tidyverse)
library(here)

### TO DO 
# Get prism data
# move script to GEE package folder, replace filenames for pheno (twice) and guilds

# Load PRISM climate data
#   Source: 
tmean <- c(terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_01_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_02_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_03_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_04_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_05_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_06_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_07_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_08_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_09_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_10_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_11_bil.bil")),
           terra::rast(here::here("..","PRISM","30_year_normals","PRISM_tmean_30yr_normal_800mM4_all_bil","PRISM_tmean_30yr_normal_800mM4_12_bil.bil")))

precip <- c(terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_01_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_02_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_03_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_04_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_05_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_06_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_07_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_08_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_09_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_10_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_11_bil.bil")),
            terra::rast(here::here("..","PRISM","30_year_normals","PRISM_ppt_30yr_normal_800mM4_all_bil","PRISM_ppt_30yr_normal_800mM4_12_bil.bil")))

climate <- c(tmean, precip)
rm(tmean, precip)

# Get subset of climate data which includes each study site
subsetClimateData <- function(site_raster_filepath, site_name)
{
  # Load existing site-level raster 
  site_raster <- terra::rast(site_raster_filepath)
  # Aggregate study site to make projection more efficient 
  site_agg <- terra::aggregate(site_raster, 50)
  # Reproject to same CRS as climate data
  site_reproj <- terra::project(site_agg, terra::crs(climate))
  # Generate random points within boundary of site data file 
  samp_points <- terra::spatSample(site_reproj, 1000, xy=TRUE) %>%
    drop_na() %>%
    dplyr::select(c("x","y"))
  # Extract data from climate raster at each point 
  climate_points <- as.data.frame(terra::extract(climate, samp_points))
  names(climate_points) <- c("index", paste("tmean_",1:12,sep=""), paste("precip_",1:12,sep=""))
  climate_points$site <- site_name
  return(climate_points)
}

climate_mcbcp <- subsetClimateData("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/phenology/pendleton_2019_phenoseries.tif", "Pendleton")
climate_vsfb <- subsetClimateData("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/phenology/vandenberg_2019_phenoseries.tif", "Vandenberg")
climate_fh <- subsetClimateData("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/phenology/san_pedro_2021_phenoseries.tif", "Huachuca")
all_climate_data <- rbind(climate_mcbcp, climate_vsfb, climate_fh)
summarizeVariable <- function(variable_name)
{
  column_names <- names(all_climate_data)
  summarized <- all_climate_data %>% 
    dplyr::select(c(column_names[which(grepl(variable_name, column_names))], "site")) %>%
    pivot_longer(1:12, 
                 names_to="varname_str", 
                 values_to="value") %>% 
    mutate(month = as.numeric(substr(varname_str, nchar(variable_name)+2, 1000))) %>%
    group_by(site, month) %>% 
    summarize(mean_val = mean(value, na.rm=TRUE),
              stdev = sd(value, na.rm=TRUE))
  return(summarized)
}
temp_summary <- summarizeVariable("tmean")
names(temp_summary) <- c("site", "month", "temp_mean", "temp_stdev")
precip_summary <- summarizeVariable("precip")
names(precip_summary) <- c("site", "month", "precip_mean", "precip_stdev")
all_climate_summary <- merge(temp_summary, precip_summary, by=c("site", "month"))

color_scale = c("Pendleton" = "green3",
                  "Vandenberg" = "blue",
                  "Huachuca" = "orange")
climate_plot <- ggplot(all_climate_summary) + 
  geom_col(aes(x=month, y=precip_mean/3, fill=site, group=site), position="dodge") + 
  geom_line(aes(x=month, y=temp_mean, col=site), size=1) + 
  scale_color_manual(values = color_scale) + 
  xlab("Month") + 
  ylab("Average Daily Temperature (°C)") + 
  scale_x_continuous(expand=c(0,0),
                     breaks=(1:12)) + 
  scale_y_continuous(expand=c(0,0),
                     breaks=(1:5)*6,
                     sec.axis=sec_axis(~ .*3, name="Average Monthly Precipitation (mm)")) + 
  theme(legend.position = "none")
climate_plot
ggsave(here::here("climatograph.png"), climate_plot, width=4.5, height=3)
# add in a second axis here for precip / tmax

# Plot for phenology
#   Veg type data (deciduous riparian woodland class is 1 in all rasters)
riparian_pendleton <- terra::rast("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/class_outputs/pendleton_all_data.tif") == 1
riparian_vandenberg <- terra::rast("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/class_outputs/vandenberg_2019_test_2019_train.tif") == 1
riparian_huachuca <- terra::rast("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/class_outputs/san_pedro_all_data.tif") == 1
# mask Huachuca to just the San Pedro 

#   Phenology data
pheno_pendleton <- terra::rast("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/phenology/pendleton_2019_phenoseries.tif")
pheno_vandenberg <- terra::rast("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/phenology/vandenberg_2019_phenoseries.tif")
pheno_huachuca <- terra::rast("D:/SERDP/GEE_Classifier/phenology_classifier_drylands/phenology/san_pedro_2021_phenoseries.tif")
#   Get phenology data for riparian veg
pheno_pendleton_df <- as.data.frame(values(c(riparian_pendleton, pheno_pendleton))) %>%
  drop_na() %>%
  filter(class == 1) %>% 
  mutate(site = "Pendleton")
pheno_vandenberg_df <- as.data.frame(values(c(riparian_vandenberg, pheno_vandenberg))) %>%
  drop_na() %>%
  filter(class == 1) %>% 
  mutate(site = "Vandenberg")
pheno_huachuca_df <- as.data.frame(values(c(riparian_huachuca, pheno_huachuca))) %>%
  drop_na() %>%
  filter(class == 1) %>% 
  mutate(site = "Huachuca")
pheno_df <- rbind(pheno_pendleton_df, pheno_vandenberg_df, pheno_huachuca_df) %>%
  pivot_longer(2:13, names_to = "month_str", values_to = "NDVI") %>%
  mutate(month = as.numeric(substr(month_str, 1, nchar(month_str)-30))+1)
pheno_summary <- pheno_df %>%
  group_by(site, month) %>%
  summarize(mean_ndvi = mean(NDVI),
            std_ndvi = sd(NDVI))
#   Create Plot
ggplot(pheno_summary) + 
  geom_line(aes(x=month, y=mean_ndvi, group=site, col=site)) + 
  scale_x_continuous(expand=c(0,0), breaks=1:12) + 
  scale_y_continuous(expand=c(0,0), breaks=(1:5)/5) + 
  ggtitle("Phenology by Study Site") + 
  scale_color_manual(values = color_scale) + 
  xlab("Month") + 
  ylab("Average NDVI")


