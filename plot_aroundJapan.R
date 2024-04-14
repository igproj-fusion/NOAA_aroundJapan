


#####################################################
#####################################################
# NOAA NCEI Global Sea Surface Temperature 


setwd("C:/Users/myaji/Documents/NOAA")
pacman::p_load(
  tidyverse,
  here,
  ncdf4,
  future,
  furrr,
  raster,
  terra,
  ggokabeito,
  rnaturalearth)


## Functions

#' Get NCDF Files from NOAA
#'
#' @param url The endpoint URL of the AVHRR data, <https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/>
#' @param subdir The subdirectory of monthly data to get. A character string of digits, of the form "YYYYMM". No default.
#'
#' @return  A directory of NCDF files.
#' @export
#'
#'
get_nc_files <- function(url = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/",
                         subdir) {
  local <- here::here("raw/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/")
  
  localdir <- paste0(local, subdir)
  
  if(!fs::dir_exists(localdir)) fs::dir_create(localdir)
  
  files <- rvest::read_html(paste0(url,subdir)) |>
    rvest::html_elements("a") |>
    rvest::html_text2()
  files <- subset(files, str_detect(files, "nc"))
  
  full_urls <- paste0(url, subdir, "/", files)
  full_outpaths <- paste0(localdir, "/", files)
  
  walk2(full_urls, full_outpaths, \(x,y) httr::GET(x, httr::write_disk(y, overwrite = TRUE)))
  
}

#
#' Remove -prelim nc files if final nc version exists
#'
#' @param subdir A character vector. The name of the subdirectory of .nc files you want to clean up.
#'
#' @return Returns "No duplicates" if there are none; otherwise silently deletes dupes.
#' @export
#'
#' @examples
clean_prelims <- function(subdir) {
  local <- here("raw/www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1/access/avhrr/")
  path <- paste0(local, subdir)
  
  all_nc_files <- basename(fs::dir_ls(path, regexp = "[.]html", invert = TRUE))
  prelim_nc_files <- basename(fs::dir_ls(path, glob = "*_preliminary.nc"))
  final_nc_files <- all_nc_files[str_detect(all_nc_files, "_preliminary",
                                            negate = TRUE)]
  
  prelim_nc_dates <- str_extract(prelim_nc_files,
                                 paste0(subdir, "\\d{2}"))
  final_nc_dates <- str_extract(final_nc_files,
                                paste0(subdir, "\\d{2}"))
  
  dupes <- intersect(prelim_nc_dates, final_nc_dates)
  
  if(rlang::is_empty(dupes)) return(message("No duplicates"))
  
  
  ## Deletion
  dupes <- paste0(dupes, collapse = "|")
  deletion_candidates <- str_detect(prelim_nc_files, dupes)
  delete_these <- prelim_nc_files[deletion_candidates]
  if(!rlang::is_empty(delete_these)) fs::file_delete(paste0(path, "/", delete_these))
  
}


## For the filename processing
## This one gives you an unknown number of chunks each with approx n elements
chunk <- function(x, n) split(x, ceiling(seq_along(x)/n))

## This one gives you n chunks each with an approx equal but unknown number of elements
#chunk2 <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))

## The Terra way. Should be considerably faster. (And it is)
## Even faster with the chunked names, to maximize layered raster processing,
## Chunks of 25 elements or so seem to work quickly enough.
#layerinfo <- tibble(
#  num = c(1:4),
#  raw_name = c("anom_zlev=0", "err_zlev=0",
#               "ice_zlev=0", "sst_zlev=0"),
#  name = c("anom", "err",
#           "ice", "sst"))

process_raster <- function(fnames, crop_area = c(-80, 0, 0, 60), layerinfo = layerinfo) {
  
  tdf <- terra::rast(fnames) |>
    terra::rotate() |>   # Convert 0 to 360 lon to -180 to +180 lon
    terra::crop(crop_area) # Manually crop to a defined box.  Default is roughly N. Atlantic lat/lon box
  
  wts <- terra::cellSize(tdf, unit = "km") # For scaling
  
  # global() calculates a quantity for the whole grid on a particular SpatRaster
  # so we get one weighted mean per file that comes in
  out <- data.frame(date = terra::time(tdf),
                    means = terra::global(tdf, "mean", weights = wts, na.rm=TRUE))
  out$var <- rownames(out)
  out$var <- gsub("_.*", "", out$var)
  out <- reshape(out, idvar = "date",
                 timevar = "var",
                 direction = "wide")
  
  colnames(out) <- gsub("weighted_mean\\.", "", colnames(out))
  out
}


## Seasons for plotting
season <-  function(in_date){
  br = yday(as.Date(c("2019-03-01",
                      "2019-06-01",
                      "2019-09-01",
                      "2019-12-01")))
  x = yday(in_date)
  x = cut(x, breaks = c(0, br, 366))
  levels(x) = c("Winter", "Spring", "Summer", "Autumn", "Winter")
  x
}

###################################

# get_nc_files(subdir = "202403")
# clean_prelims(subdir = "202403")
# get_nc_files(subdir = "202404")
# clean_prelims(subdir = "202404")


# Get filenames
# All the daily .nc files we downloaded:
all_fnames <- fs::dir_ls(here("raw"), 
                         recurse = TRUE, glob = "*.nc")
fnames2024 <- str_subset(all_fnames, ".2024")
chunked_fnames <- chunk(fnames2024, 25)

crop_bb <- c(110, 170, 10, 60)

df_2024 <- future_map(chunked_fnames, process_raster,
                 crop_area = crop_bb) |>
  list_rbind() |>
  as_tibble() |>
  mutate(date = ymd(date),
         year = lubridate::year(date),
         month = lubridate::month(date),
         day = lubridate::day(date),
         yrday = lubridate::yday(date),
         season = season(date))



FILE <- "aroundJapan1981~2023.csv"
df_org <- read.csv(here("data", FILE))

df <- rbind(df_org, df_2024)




###################################
## World Graph

colors <- ggokabeito::palette_okabe_ito()
colors[2] <- "cornflowerblue"

month_labs <- df |>
  filter(year == 1990,
         day == 15) |>
  dplyr::select(date, year, yrday, month, day) |>
  mutate(month_lab = month(date, label = TRUE, abbr = TRUE))

avg <- df |>
  filter(year > 1981 & year < 2012) |>
  group_by(yrday) |>
  filter(yrday != 366) |>
  summarize(mean_8211 = mean(sst, na.rm = TRUE),
            sd_8211 = sd(sst, na.rm = TRUE)) |>
  mutate(fill = colors[2],
         color = colors[2])

out <- df |>
  mutate(year_flag = case_when(
    year == 2023 ~ "2023",
    year == 2024 ~ "2024",
    .default = "All other years"))

df_2023 <- df |>
  filter(year == 2023)

today <- df_2024 |> 
  filter(date == max(df_2024$date))

LastDate <- gsub("-", "/", gsub("-0", "/", max(df_2024$date)))

g1 <- ggplot() +
  geom_hline(yintercept = 18, color = "gray80") +
  geom_hline(yintercept = 20, color = "gray80") +
  geom_hline(yintercept = 22, color = "gray80") +
  geom_hline(yintercept = 24, color = "gray80") +
  geom_hline(yintercept = 26, color = "gray80") +
  geom_ribbon(data = avg,
              mapping = aes(x = yrday,
                            ymin = mean_8211 - 2*sd_8211,
                            ymax = mean_8211 + 2*sd_8211,
                            fill = fill),
              alpha = 0.3,
              inherit.aes = FALSE) +
  geom_line(data = avg,
            mapping = aes(x = yrday,
                          y = mean_8211,
                          color = color),
            linewidth = 1.5,
            inherit.aes = FALSE) +
  scale_color_identity(name = "Mean Temp. 1982-2011, ±2σ", guide = "legend",
                       breaks = unique(avg$color), labels = "") +
  scale_fill_identity(name = "Mean Temp. 1982-2011, ±2σ", guide = "legend",
                      breaks = unique(avg$fill), labels = "") +
  ggnewscale::new_scale_color() +
  geom_line(data = out,
            mapping = aes(x = yrday, y = sst, group = year, color = year_flag),
            inherit.aes = FALSE, alpha = 0.5) +
  scale_color_manual(values = c("orange", "darkmagenta", "grey50")) +
  scale_x_continuous(breaks = month_labs$yrday, labels = month_labs$month_lab) +
  scale_y_continuous(breaks = seq(16, 26, 2),
                     limits = c(16, 27),
                     expand = expansion(mult = c(-0.05, 0.05))) +
  geom_line(linewidth = rel(1)) +
  geom_line(data = df_2023, aes(x = yrday, y = sst),
            color = "orange", linewidth = 0.8) +
  geom_line(data = df_2024, aes(x = yrday, y = sst),
            color = "darkmagenta", linewidth = 1) +
  geom_point(data = today, aes(x = yrday, y = sst),
            color = "darkmagenta", size = 2) +
  geom_line(data = avg,
            mapping = aes(x = yrday,
                          y = mean_8211),
            linewidth = 1.5, color = "cornflowerblue") +
  theme_classic() +
  guides(
    x = guide_axis(cap = "both"),
    y = guide_axis(minor.ticks = TRUE, cap = "both"),
    color = guide_legend(override.aes = list(linewidth = 1.5))
  ) +
  labs(x = "", y = "Mean Temperature (°Celsius)",
       color = "Year",
       title = "Mean Daily Sea Surface Temperature, Sea areas around Japan",
       subtitle = paste0("Gridded and weighted NOAA OISST v2.1 estimates, 1981/9/1~", 
                         LastDate),
       caption = "https://www.ncei.noaa.gov/data/sea-surface-temperature-optimum-interpolation/v2.1") +
  theme(plot.margin= unit(c(1, 1, 1, 1), "lines"),
        axis.line = element_line(color = "gray30", linewidth = rel(1)),
        plot.title = element_text(size = rel(1.4)),
        plot.subtitle = element_text(size = rel(1.1)),
        plot.caption = element_text(size = rel(1.0)),
        axis.text = element_text(size = rel(0.9)),
        legend.position = "top")


world_map <- ne_countries(scale = "large",
                          returnclass = "sf")
g2 <-  ggplot() +
  theme_void() +
  geom_sf(data = world_map) +
  geom_rect(aes(xmin = 110, xmax = 170, ymin = 10, ymax =60), 
            fill="lightblue", color= NA, alpha = 0.8) +
  geom_sf(data = world_map, linewidth = 0.1, color = "gray60") +
  xlim(c(110, 170)) + ylim(c(10, 60)) 


g <- g1 + annotation_custom(ggplotGrob(g2), xmin = -10, xmax = 80, 
                       ymin = 22, ymax = 27)
print(g)
