#  Carga paquetes -----
library(rgee)
library(stars)
library(ggplot2)

library(tmap)
# Inizialización rgee ----
ee_Initialize(drive = TRUE)

# Definicion funciones ----
cambio_unidades <-  function(img) {
  # bands <- img$selec('total_precipitation')$multiply(1000)
  
  precip = img$expression(
    'pp * 1000', list(
      pp = img$select('total_precipitation')
    )
  )$rename('total_precipitation_mm')
  
  temp = img$expression(
    'temp - 273.15', list(
      temp = img$select('temperature_2m')
    )
  )$rename('temperature_2m_c')

  img$
    addBands(precip)$
    addBands(temp)$
    select(c('temperature_2m_c', 'total_precipitation_mm'))
}


# Archivos ----

## Lectura archivos Provincias ----
my_polygons <- sf::st_read("data/provincia.json")
my_polygons_province <- my_polygons[-c(2, 10, 11, 18, 23) , c("fna", "gna", "nam")]
my_polygons_province <- st_make_valid(my_polygons_province)

## Lectura archivo virus ----

# my_deasease_data <- sf::st_read("data/trigo_maiz_26_09_22.gpkg")
# my_deasease_data <- my_deasease_data[,c(2,4,7,8,11,13)]
# my_deasease_data <- my_deasease_data[!is.na(my_deasease_data$WSMV) & !is.na(my_deasease_data$HPV), ]
# my_deasease_data <- my_deasease_data %>%
#   dplyr::group_by(Año.de.Colecta) %>%
#   dplyr::distinct() %>%
#   dplyr::ungroup() %>%
#   dplyr::mutate(ID = seq_len(dplyr::n()))
# sf::st_write(my_deasease_data, "data/virus.gpkg", delete_layer = TRUE)

my_deasease_data <- sf::st_read("data/virus.gpkg")
plot(my_deasease_data)


## Union Archivos ----
my_deasease_province_data <- 
  st_join(my_deasease_data, my_polygons_province, largest = TRUE)

my_samples_provinces <- table(my_deasease_province_data$fna)


my_provinces_samples <- 
  my_polygons_province[my_polygons_province$fna %in%  names(my_samples_provinces[my_samples_provinces > 5]), ]

my_deasease_province_data <-
  my_deasease_province_data[my_deasease_province_data$fna %in% names(my_samples_provinces[my_samples_provinces > 5]),]
names(my_samples_provinces[my_samples_provinces > 5])


tm_shape(my_provinces_samples) +
  tm_polygons() +
  tm_shape(my_deasease_province_data) +
  tm_dots()



# my_polygons_dpto <- sf::st_read("data/departamento.json")
# my_polygons_dpto <- st_make_valid(my_polygons_dpto)
# # st_intersection(my_polygons_dpto, my_polygons_province)
# my_polygons_dpto <- st_make_valid(my_polygons_dpto[my_polygons_province,])
# #
# # plot(my_polygons_dpto)
# #
# my_polygons_dpto_join <-
#   st_join(my_polygons_dpto, my_polygons_province, join = st_nearest_feature)

# plot(my_polygons_dpto[,1])
# my_polygons_dpto_join <- sf::st_simplify(my_polygons_dpto_join)
# plot(my_polygons_dpto_join[,"fna.y"])

my_polygons_provinces <- sf::st_simplify(my_provinces_samples)


my_polygons <- 
  my_polygons_provinces |>
  sf::st_buffer(0.01) |>
  sf::st_union() |>
  sf::st_boundary() |>
  sf::st_cast("POLYGON") |>
  sf::st_make_valid() |>
  sf::st_simplify() |>
  sf::st_as_sf()
st_geometry(my_polygons) <- "geometry"
plot(my_polygons)


my_polygons_ee <- sf_as_ee(my_polygons)


table <- ee$FeatureCollection('USDOS/LSIB_SIMPLE/2017')



# my_polygons_ee_simple = sf_as_ee(my_polygons_province)$
#   map(function(feature) {
#     feature$bounds()
#   })

era5_ee <- ee$ImageCollection("ECMWF/ERA5_LAND/MONTHLY")$
  select(c('temperature_2m', 'total_precipitation'))$
  filterDate('2021-01-01' , '2021-12-31' )$
  filterBounds(my_polygons_ee)$
  map(cambio_unidades)#$
  # map(function(img){img$clipToCollection(my_polygons_ee)})


# https://github.com/csaybar/rgee/blob/examples//ImageCollection/overview.R
# Resolution 
era5_resolution <- era5_ee$first()$projection()$nominalScale()
cat("ERA5 resolution: ", era5_resolution$getInfo(), '\n')

# Get band names
bandNames = era5_ee$first()$bandNames()$getInfo()
cat("Band Names: ", paste(bandNames, collapse = ", "), '\n')

# Get the number of images.
count <- era5_ee$size()
cat("Count: ", count$getInfo(), '\n')

# Get the date range of images in the collection.
range <- era5_ee$reduceColumns(
  ee$Reducer$minMax(),
  list("system:time_start")
)

col_min <- eedate_to_rdate(range$get("min"))
col_max <- eedate_to_rdate(range$get("max"))
cat("Date range: ", as.character(col_min), as.character(col_max), '\n')


# Image dates
dates = era5_ee$
  map(function(image) {
   ee$Feature(NULL, 
              list('date' = image$date()$format('YYYY-MM-dd')))})$
  distinct('date')$
  aggregate_array('date')$
  getInfo()

cat("Images date (", length(dates), "):", paste(dates, collapse = ", "), '\n')


# Map$addLayers(era5_ee)

# Extract values ----

myErasValues <- lapply(bandNames, function(band) {
  my_stars <- ee_as_stars(
    image = era5_ee$select(band)$toBands(),
    region = my_polygons_ee$geometry(),
    scale = era5_resolution
  )
  st_as_stars(my_stars)
})


ggplot() + geom_stars(data = myErasValues[[1]]) +
  coord_equal() +
  facet_wrap( ~ band) +
  theme_void() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))


ggplot() + geom_stars(data =  myErasValues[[2]]) +
  coord_equal() +
  facet_wrap( ~ band) +
  theme_void() +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0))


#### zonal ----
# Define parameters for the zonalStats function.
# params = list(
#   reducer = ee$Reducer$median(),
#   scale = era5_resolution,
#   crs = era5_ee$toBands()$projection(),
#   bands = c('total_precipitation', 'temperature_2m'),
#   datetimeName = 'date',
#   datetimeFormat = 'YYYY-MM-dd')
# 



#### 1 ----
# myNewPolygon <- sf::st_cast(my_polygons_province, "POLYGON")
# myNewPolygon <- sf::st_simplify(my_polygons_province)myNewPolygon <- sf::st_simplify(my_polygons_province)
# 
# 
# ee_arg_rain <- ee_extract(
#   x = era5_ee,
#   y = myNewPolygon,
#   scale = era5_resolution$getInfo(),
#   fun = ee$Reducer$mean(),
#   via = "drive",
#   lazy = TRUE,
#   sf = TRUE
# )
# 
# 
# 
# plot(myNewPolygon)
# # Extract values
# ee_arg_rain <- ee_extract(
#   x = era5_ee$select('total_precipitation')$toBands(),
#   y = myNewPolygon,
#   scale = era5_resolution,
#   fun = ee$Reducer$mean(),
#   sf = FALSE
# )
# 
# # gganimate
# 
# colnames(ee_arg_rain)[-1] <- sprintf("%02d", 1:12)
# 
# ee_arg_rain %>%
#   pivot_longer(-NAME, names_to = "month", values_to = "pr") %>%
#   ggplot(aes(x = as.integer(month), y = pr, color = pr)) +
#   geom_line(alpha = 0.8, size = 2) +
#   xlab("Month") +
#   ylab("Precipitation (mm)") +
#   theme_minimal() +
#   labs(title = "{closest_state}") +
#   transition_states(NAME) +
#   shadow_mark(size = 0.4, colour = "grey")




## https://github.com/r-spatial/rgee/issues/220

# Statistics in batches ---------------------------------------------------
library(foreach)
nFold <- 5
myNewPolygon <- sf::st_simplify(my_provinces_samples)

myFold <- rep(seq_len(nFold), length = nrow(myNewPolygon))
myNewPolygon$fold <- myFold


# Do something to every element of a collection.
withMoreProperties = sf_as_ee(myNewPolygon)$map(function(f) {
  # Set a property.
  f$set("area_sq_meters", f$area())
})
print(withMoreProperties$first()$get("area_sq_meters")$getInfo())


vectors <- era5_ee$map(function(img) {
  img$reduceToVectors(
  reducer = ee$Reducer$mean(),
  geometry = sf_as_ee(myNewPolygon),
  # Set the scale to the maximum possible given
  # the required precision of the computation.
  scale = 20000
)})

Map$addLayer(vectors$first())
vectors_r <- ee_as_stars(
  image = vectors$toBands(),
  region =  sf_as_ee(myNewPolygon)$geometry()$simplify(5),
  scale = era5_resolution * 10
) 


if (file.exists("cache/statistics_virus.rds")) {
  statistics <- readRDS("cache/statistics_virus.rds")
} else {
  statistics <-
    foreach(index = seq_len(nFold), .combine = rbind) %do% {
      ee_poly <- myNewPolygon[myNewPolygon[['fold']] == index, ] %>%
        st_simplify() %>%
        st_as_sfc() %>%
        sf_as_ee()
      
      ee_poly <- ee_poly$simplify(5)
      
      ee_extract(
        era5_ee,
        ee_poly,
        fun = ee$Reducer$mean(),
        scale = era5_resolution,
        via = "drive",
        sf = TRUE
      )
    }
  statistics
  # plot(myNewPolygon)
  # plot(statistics[,1])
  statistics$nam <- myNewPolygon$nam.y
  statistics$nam_depto <- myNewPolygon$nam.x
  saveRDS(statistics, "cache/statistics_virus.rds")
}






province_stat <- statistics %>%
  tidyr::pivot_longer(cols = starts_with('X'),
                      # names_to = "month_var",
                      # values_to = "value",
                      names_pattern = "X(\\d{6})(.*)$",
                      names_to = c(".value", "variable")) %>% 
  tidyr::pivot_longer(cols = matches('\\d{6}'),
                      names_to = "month",
                      values_to = "value")


ggplot(province_stat) +
  geom_sf(aes(color = value)) +
  facet_wrap(variable ~ month)


province_stat_mean <-
  province_stat %>% 
  st_drop_geometry() %>% 
  dplyr::group_by(nam, variable, month) %>% 
  dplyr::summarise(value = mean(value)) 




library(gganimate)
  ggplot(province_stat_mean, aes(x = as.integer(month), y = value, color = value)) +
  geom_line(alpha = 0.8, size = 2) +
  xlab("Month") +
  ylab("Precipitation (mm)") +
  theme_minimal() +
  facet_wrap(variable ~ ., scales = 'free') +
  labs(title = "{closest_state}") +
  transition_states(nam) +
  shadow_mark(size = 0.4, colour = "grey")


