# Carga paquetes -----
library(rgee)
library(stars)
library(tmap)
library(foreach)
library(dplyr)
library(ggplot2)
library(lubridate)
library(patchwork)

if (!dir.exists('images')) {
  message('Directory images will be created at ', getwd())
  dir.create('images')
}

if (!dir.exists('cache')) {
  message('Directory cache will be created at ', getwd())
  dir.create('cache')
}

# Inizialización rgee ----
ee_Initialize(drive = TRUE)

# Definicion funciones ----
cambio_unidades <-  function(img) {

  precip = img$expression(
    'pp * 1000', list(
      pp = img$select('total_precipitation_sum')
    )
  )$rename('total_precipitation_mm')
  
  temp = img$expression(
    'temp - 273.15', list(
      temp = img$select('temperature_2m')
    )
  )$rename('temperature_2m_c')
  
  temp_pto_rocio = img$expression(
    'dpt - 273.15', list(
      dpt = img$select('dewpoint_temperature_2m')
    )
  )$rename('dewpoint_temperature_2m_c')
  
  img$
    addBands(precip)$
    addBands(temp)$
    addBands(temp_pto_rocio)$
    select(c('total_precipitation_mm', 
             'temperature_2m_c', 
             'dewpoint_temperature_2m_c'))
}

select_polygons <- function(feature){
  filteredGeoms = feature$
    geometry()$
    geometries()$
    map(ee_utils_pyfunc({function(geometry) {
      geometry = ee$Geometry(geometry)
      ee$Algorithms$
        If(geometry$type()$
             compareTo('Polygon'), NULL, geometry)
    }}), dropNulls = TRUE)
  feature$setGeometry(ee$Geometry$MultiPolygon(filteredGeoms))
}
# Archivos ----

## Descarga archivos Provincias ----
my_polygons_country <- ee$FeatureCollection("FAO/GAUL/2015/level1")$
  filter(ee$Filter$eq('ADM0_NAME', 'Argentina'))$
  select(c("ADM1_NAME","Shape_Area", "Shape_Leng"))$
  map(select_polygons)

## Lectura archivo virus ----
my_deasease_data <- sf::st_read("data/virus.gpkg")

# Número de observaciones
cat("Número de observaciones: ", nrow(my_deasease_data), "\n")

# Número de sitios
n_sitios <- 
  my_deasease_data %>%
  sf::st_geometry() %>% 
  sf::st_as_sf() %>% 
  dplyr::distinct() %>% 
  nrow()
cat("Número de sitios observados: ", n_sitios, "\n")

## Union Archivos ----
my_polygons_country_sf <- ee_as_sf(my_polygons_country)

my_polygons_country_sf$n_obs <- 
  st_intersects(my_polygons_country_sf, 
                my_deasease_data) %>%
  lengths()

my_polygons_province_sf <- 
  my_polygons_country_sf %>%
  filter(n_obs > 0)

my_deasease_province_data <- 
  st_join(my_deasease_data, 
          my_polygons_province_sf, 
          largest = TRUE)

# Mapa de sitios de muestreo
muestreo <-
  tm_shape(my_polygons_country_sf) +
  tm_polygons() +
  tm_shape(my_polygons_province_sf) +
  tm_polygons(col = '#bdfcf9') +
  tm_shape(my_deasease_province_data) +
  tm_dots(size = 0.15, 
          # col = 'Año.de.Colecta', 
          # style = 'cat',
          # title = "Año de\nmuestreo"
          ) +
  tm_compass(position = c("left", "top")) +
  tm_scale_bar(text.size = 1) +
  tm_layout(legend.outside = TRUE,
            legend.format = list(big.mark = ".",
                                 decimal.mark = ","),
            frame = FALSE) 

tmap_save(
  muestreo,
  filename = 'images/Figura_1_SitiosMuestreo.png',
  width = 7.5,
  height = 10,
  units = 'cm'
)

# Especificacion de fechas desde la que se descargaran datos
sf_table_dates <- my_deasease_province_data
sf_table <- my_deasease_province_data["geom"]
st_geometry(sf_table) <- 'geometry'

## Manejo de fechas con paquete lubridate
sf_table <-
  my_deasease_province_data %>% 
  mutate(start_date = ifelse(Especie == 'Trigo', 
                             paste0("20/3/", Año.de.Colecta), 
                             paste0("21/9/", Año.de.Colecta)),
         start_date = ymd(as.Date(start_date, "%d/%m/%Y")),
         end_date = start_date %m+% months(3),
         end_date = ymd(as.Date(end_date, "%d/%m/%Y"))) %>% 
  rename(Anio_de_Colecta = Año.de.Colecta)

ee_table <- sf_as_ee(sf_table) # Carga de sf como FeatureCollection

# Otención de imágenes ----
my_polygons_province <- sf_as_ee(my_polygons_province_sf)

## DEM ----
dem_ee <- ee$Image('CGIAR/SRTM90_V4')$
  select('elevation')

if (file.exists("cache/dem.rds")) {
  dem <- readRDS("cache/dem.rds")
} else {
  dem <- rgee::ee_as_stars(
    dem_ee,
    my_polygons_province$geometry(),
    scale = era5_resolution,
    via = "drive"
  ) %>%
    stars::st_as_stars()
  
  # DEM
  saveRDS(dem, "cache/dem.rds")
}

st_get_dimension_values(dem, "band")
plot(dem)

## Datos Meteorológicos ----


min_date <- min(as.Date(sf_table$start_date))
max_date <- max(as.Date(sf_table$end_date))

era5_ee <- ee$ImageCollection("ECMWF/ERA5_LAND/DAILY_RAW")$
  filterBounds(my_polygons_province$geometry())$
  filterDate(as.character(min_date), as.character(max_date))$
  select(c('total_precipitation_sum', 'temperature_2m', 'dewpoint_temperature_2m'))$
  map(cambio_unidades)

# https://github.com/csaybar/rgee/blob/examples//ImageCollection/overview.R
## Resolución del producto 
era5_resolution <- era5_ee$first()$projection()$nominalScale()
cat("ERA5 resolution: ", era5_resolution$getInfo(), '\n')

# Nombre de bandas
bandNames <- era5_ee$first()$bandNames()$getInfo()
cat("Band Names: ", paste(bandNames, collapse = ", "), '\n')

# Número de imágenes
count <- era5_ee$size()
cat("Count: ", count$getInfo(), '\n')

# Obtener el rango de fechas de las imágenes
range <- era5_ee$reduceColumns(
  ee$Reducer$minMax(),
  list("system:time_start")
)

col_min <- eedate_to_rdate(range$get("min"))
col_max <- eedate_to_rdate(range$get("max"))
cat("Date range: ", as.character(col_min), " - ", as.character(col_max), '\n')

# Extracción de variables bioclimáticas ----
temporal_aggregate <- function(monthN, delta = 2, unit = 'week') {
  t1 <- ee$Date(monthN)
  t2 <- t1$advance(delta, unit)
  semanal <- era5_ee$
    filterDate(t1, t2)
  
  pp_semanal <- semanal$select('total_precipitation_mm')$
    sum()$
    set(
      list(
        'system:time_start' = t1$millis()
      ) 
    )
  
  temp_semanal <- semanal$
    select('temperature_2m_c')$
    mean()$
    set(
      list(
        'system:time_start' = t1$millis()
      ) 
    )
  
  ptoRocio_semanal <- semanal$
    select('dewpoint_temperature_2m_c')$
    mean()$
    set(
      list(
        'system:time_start' = t1$millis()
      ) 
    )

  temp_semanal$
    addBands(pp_semanal)$
    addBands(ptoRocio_semanal)
}

## Descarga para sitios y años muestreados ----
if (file.exists("cache/statistics_mean_virus.rds")) {
  statistics <- readRDS("cache/statistics_mean_virus.rds")
} else {
  statistics <-
    foreach(my_date = unique(sf_table$start_date), .combine = rbind) %do% {
      
      my_ss <- sf_table[sf_table[['start_date']] == my_date, ]
      st_geometry(my_ss) <- 'geometry'
      months <- seq(unique(my_ss[['start_date']]), unique(my_ss[['end_date']]), by = '2 week') %>% 
        as.character() %>% 
        lapply(ee$Date) %>% 
        ee$List()
      
      era5_final <- months$map(ee_utils_pyfunc(function(x) {
        temporal_aggregate(x, 2, 'week')
        })) %>% 
        ee$ImageCollection()

      ee_extract(era5_final, 
                 my_ss, 
                 scale = era5_resolution,
                 via = "getInfo",
                 sf = TRUE)
    }
  # statistics
  saveRDS(statistics, "cache/statistics_mean_virus.rds")
}

### Reformateo de base de datos para graficar -----
statistics_long_sf <- 
  statistics %>% 
  tidyr::pivot_longer(
    cols = -c(ID:end_date,geometry),
    names_sep = 3,
    names_to = c('semana', 'variable')
  ) %>% 
  mutate(estacion = ifelse(grepl("*-09-*", start_date), "primavera", "otoño"),
         anio = year(start_date)) %>% 
  tidyr::drop_na(value)

semana_labels <- c(
  X0_ = "Semanas:\n1 y 2",
  X1_ = "Semanas:\n3 y 4",
  X2_ = "Semanas:\n5 y 6",
  X3_ = "Semanas:\n7 y 8",
  X4_ = "Semanas:\n9 y 10",
  X5_ = "Semanas:\n11 y 12",
  X6_ = "Semanas:\n13 y 14"
)
variable_labels <- c(
  total_precipitation_mm = "Precipitación\n(mm)",
  temperature_2m_c = "Temperatura 2m\n(°C)",
  dewpoint_temperature_2m_c = "Punto de rocío 2m\n(°C)"
)

presAus_hpv_plt <-
  statistics_long_sf %>% 
  tidyr::drop_na(HPV) %>% 
  ggplot(aes(as.factor(HPV), value)) +
  stat_summary(fun = 'mean',
               geom = "point") +
  labs(x = "HPWMoV",
       y = "Valor") +
  facet_grid(variable ~ semana,
             scales = "free_y",
             labeller = labeller(
               variable = variable_labels,
               semana = semana_labels
             )) +
  scale_x_discrete(labels = c("0" = "No",
                              "1" = "Si")) +
  theme(
    strip.text = element_text(size = 8),
    legend.position = 'bottom',
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave('images/Figura_4_presenciaAusencia_hpmwov.png',
       presAus_hpv_plt,
       height = 12.0,
       width = 15.0,
       units = 'cm')

presAus_wsmv_plt <-
  statistics_long_sf %>% 
  tidyr::drop_na(WSMV) %>% 
  ggplot(aes(as.factor(WSMV), value)) +
  stat_summary(fun = 'mean',
               geom = "point") +
  labs(x = "WSMV",
       y = "Valor") +
  facet_grid(variable ~ semana,
             scales = "free_y",
             labeller = labeller(
               variable = variable_labels,
               semana = semana_labels
             )) +
  scale_x_discrete(labels = c("0" = "No",
                              "1" = "Si")) +
  theme(
    strip.text = element_text(size = 8),
    legend.position = 'bottom',
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

ggsave('images/Figura_5_presenciaAusencia_wsmv.png',
       presAus_wsmv_plt,
       height = 12.0,
       width = 15.0,
       units = 'cm')


# Ajuste Modelo Logístico ----
## Suma de cuadrados tipo III
options(contrasts = c("contr.sum", "contr.poly"))

## Seleccion de covariables para armar la formula
isCovariable <- agrepl("temperature|dewpoint|total_precipitation", colnames(statistics), fixed = FALSE)
covariable <- colnames(statistics)[isCovariable]
covariable_fml <- paste(covariable, collapse = " + ")


## Ajuste del modelo

mdl_hpv <- 
  glm(
    paste("HPV", "~", covariable_fml),
    family = binomial(link = "logit"),
    data = statistics
  )

deviance(mdl_hpv)
summary(mdl_hpv)

mdl_wsmv <-
  glm(
    paste("WSMV", "~", covariable_fml),
    family = binomial(link = "logit"),
    data = statistics
  )

deviance(mdl_wsmv)
summary(mdl_wsmv)


## Descarga para toda la zona y años muestreados ----

if (file.exists("cache/statistics_mean_provinces.rds")) {
  statistics_province <- readRDS("cache/statistics_mean_provinces.rds")
} else {
  statistics_province <-
    foreach(my_date = unique(sf_table$start_date), 
            .combine = c) %do% {
              my_ss <- sf_table[sf_table[['start_date']] == my_date, ]
              sf::st_geometry(my_ss) <- 'geometry'
              months <- seq(unique(my_ss[['start_date']]), unique(my_ss[['end_date']]), by = '2 week') %>% 
                as.character() %>% 
                lapply(ee$Date) %>%
                ee$List()
              
              
              months <- my_date %>% 
                as.character() %>%
                lapply(ee$Date) %>%
                ee$List()
              
              era5_final <-
                months$map(rgee::ee_utils_pyfunc(function(date) {
                  temporal_aggregate(date, 3, 'month')})) %>%
                ee$ImageCollection()
             
              my_stars_bioclimatic <- rgee::ee_as_stars(
                era5_final$toBands(),
                my_polygons_province$geometry(),
                scale = era5_resolution,
                via = "drive"
              ) %>%
                stars::st_as_stars()
              
            }
  names(statistics_province) <- unique(sf_table$start_date)
  statistics_province <- merge(statistics_province, name = "time")
  names(statistics_province) <- 'value'
  saveRDS(statistics_province, "cache/statistics_mean_provinces.rds")
}


myTimes <- st_get_dimension_values(statistics_province, 'time')
primavera <- grepl("*-09-*", myTimes)

estacion_label <- function(string) {
  my_vector <- unlist(string)
  estacion <- ifelse(grepl("*-09-*", my_vector), "Primavera", "Otoño")
  anio <- gsub("-.*$","\\1", my_vector)
  my_label <- anio#paste0(estacion, ":\n ", anio)
  my_label <- data.frame(my_label)
  colnames(my_label) <- colnames(string) 
  my_label
}

calculate_range <- function(x) {
  my_values <- c(min = purrr::map(x, min, na.rm = TRUE),
                 max = purrr::map(x, max, na.rm = TRUE)) %>% 
    unlist()
  names(my_values) <- c('min', 'max')
  my_values
}


# Media de variables meteoroógica para todos los años muestreados -----
statistics_otono <-
  st_apply(statistics_province[, , , , which(!primavera)], c(1, 2, 3), mean, na.rm = TRUE)

## Descritiva promedio ----
### Temperatura ----
mean(pull(statistics_otono[,,,1]), na.rm = TRUE)
min(pull(statistics_otono[,,,1]), na.rm = TRUE)
max(pull(statistics_otono[,,,1]), na.rm = TRUE)
### Precipitaciones ----
mean(pull(statistics_otono[,,,2]), na.rm = TRUE)
min(pull(statistics_otono[,,,2]), na.rm = TRUE)
max(pull(statistics_otono[,,,2]), na.rm = TRUE)
### Punto Rocio ----
mean(pull(statistics_otono[,,,3]), na.rm = TRUE)
min(pull(statistics_otono[,,,3]), na.rm = TRUE)
max(pull(statistics_otono[,,,3]), na.rm = TRUE)


cat( "El promedio de temperatura para toda la " ,
     "región fue de ",
     round(mean(pull(statistics_otono[,,,1]), na.rm = TRUE), 1),
     "°C, el valor máximo promedio observado en los años en estudio ",
     "fue de ",
     round(max(pull(statistics_otono[,,,1]), na.rm = TRUE), 1),
     "°C, el valor mínimo promedio fue de " ,
     round(min(pull(statistics_otono[,,,1]), na.rm = TRUE), 1),
     "°C observándose en la zona "  ,
     "noroeste de la región. ", 
     "Las precipitaciones promedio para la región en estudio, ",
     "fueron de ",
     round(mean(pull(statistics_otono[,,,2]), na.rm = TRUE), 0),
     "mm para otoño, variando entre ", 
     round(min(pull(statistics_otono[,,,2]), na.rm = TRUE), 0),
     "mm y ",
     round(max(pull(statistics_otono[,,,2]), na.rm = TRUE), 0),
     "mm. La temperatura de punto de rocío promedio fue de ",
     round(mean(pull(statistics_otono[,,,3]), na.rm = TRUE), 1),
     "°C. En el noroeste de la región en estudio se observaron ",
     "las menores temperaturas de punto de rocío (",
     round(min(pull(statistics_otono[,,,3]), na.rm = TRUE), 1),
     "°C), por efecto de la gran altitud de esa zona.",
     "\n",
     sep = "")



## Graficos ----

temp_range <- calculate_range(statistics_otono[,,,1])
rain_range <- calculate_range(statistics_otono[,,,2])
dewpt_range <- calculate_range(statistics_otono[,,,3])


temperatura_promedio_ggplt <- 
  ggplot() + 
  geom_stars(data = statistics_otono[,,,1]) +
  geom_sf(data = my_polygons_province_sf, fill = NA) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  scale_fill_gradientn(limits = temp_range,
                       colours = cptcity::cpt(pal = "arendal_temperature")) +
  labs(fill = "Temperatura\n(°C)",
       y = NULL, x = NULL) +
  theme(axis.text.x = element_text(angle = 90))


pp_promedio_ggplt <-
  ggplot() + 
  geom_stars(data = statistics_otono[,,,2]) +
  geom_sf(data = my_polygons_province_sf, fill = NA) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  scale_fill_gradientn(limits = rain_range,
                       colours = cptcity::cpt(pal = "colo_alpen_sky_rain_brolly")) +
  labs(fill = "Precipitación\n(mm)",
       y = NULL, x = NULL) +
  theme(axis.text.x = element_text(angle = 90))

ptoro_promedio_ggplt <-
  ggplot() + 
  geom_stars(data = statistics_otono[,,,3]) +
  geom_sf(data = my_polygons_province_sf, fill = NA) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  # facet_wrap(. ~ epoca) +
  scale_fill_gradientn(limits = dewpt_range,
                       colours = cptcity::cpt(pal = "idv_relative_humidity")) +
  labs(fill = "Punto de rocío\n(°C)",
       y = NULL, x = NULL) +
  theme(axis.text.x = element_text(angle = 90))


dem_ggplt <- 
  ggplot() + 
  geom_stars(data = dem) +
  geom_sf(data = my_polygons_province_sf, fill = NA) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  scale_fill_gradientn(colours = 
                         cptcity::cpt(pal = "grass_elevation")) +
  labs(fill = "Elevación\n(msnm)",
       y = NULL, x = NULL) +
  theme(axis.text.x = element_text(angle = 90))


temp_pp_ptoro_prom_dem_ggplt <-
  temperatura_promedio_ggplt +
  pp_promedio_ggplt +
  ptoro_promedio_ggplt +
  dem_ggplt +
  plot_annotation(tag_levels = 'a')
  
ggsave('images/Figura_2_TempPrecipitacionesTemperaturaPtoRocioPromedioDem.png',
       temp_pp_ptoro_prom_dem_ggplt,
       height = 11.8,
       width = 15,#10
       units = 'cm')







## Figuras no presentadas ----
# 
# ggsave('images/pp_promedio_roi.png',
#        pp_promedio_ggplt,
#        height = 10,
#        width = 15.0/2,
#        units = 'cm')
# 
# 
# 
# ggsave('images/puntorocio_promedio_roi.png',
#        ptoro_promedio_ggplt,
#        height = 10,
#        width = 15.86/2,
#        units = 'cm')


## Precipitacion promedio ----
precipitacion_ggplt <- 
  ggplot() +
  geom_stars(data = statistics_province[, , , 2, ]) +
  geom_sf(data = my_polygons_province_sf, fill = NA, size = 0.2) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  facet_wrap(. ~ time,
             ncol = 5,
             labeller = estacion_label) +
  scale_fill_gradientn(colours = cptcity::cpt(pal = "colo_alpen_sky_rain_brolly")) +
  labs(fill = "Precipitación\n(mm)",
       y = NULL, x = NULL) +
  # theme_void() +
  theme(
    legend.position = 'bottom',
    legend.key.width = unit(1.3, "cm"),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks.length = unit(0, "pt"),
    axis.ticks.length.x = NULL,
    axis.ticks.length.x.top = NULL,
    axis.ticks.length.x.bottom = NULL,
    axis.ticks.length.y = NULL,
    axis.ticks.length.y.left = NULL,
    axis.ticks.length.y.right = NULL
  )

ggsave('images/Figura_3_VariabilidadEspacioTemporalPrecipitacion.png',
       precipitacion_ggplt,
       height = 18.5,
       width = 15.00,
       units = 'cm')



## Temperatura promedio ----
mean(pull(statistics_province[,,,1,]), na.rm = TRUE)
min(pull(statistics_province[,,,1,]), na.rm = TRUE)
max(pull(statistics_province[,,,1,]), na.rm = TRUE)

## Valor promedio de cada banda en las fechas de muestreo
tempral_mean <- st_apply(statistics_province[,,,,], c(3,4), mean, na.rm = TRUE)
temp_temporal_media <- pull(tempral_mean)
rownames(temp_temporal_media) <- st_get_dimension_values(tempral_mean, 'band')
colnames(temp_temporal_media) <- st_get_dimension_values(tempral_mean, 'time')
temp_temporal_media

# Temperaturas promedio de los años de muestreo
temperaturas_promedio <- temp_temporal_media["0_temperature_2m_c", order(temp_temporal_media[1,])]

# Rango de temperatura
round(range(temperaturas_promedio), 1)
# Anio con temperatura promedio menor
which.min(temperaturas_promedio)
# Anio con temperatura promedio mayor
which.max(temperaturas_promedio)


# Temperaturas promedio de los años de muestreo
precipitaciones_promedio <-
  temp_temporal_media["0_total_precipitation_mm", order(temp_temporal_media[1,])]

# Rango de temperatura
round(range(precipitaciones_promedio), 1)
# Anio con temperatura promedio menor
round(min(precipitaciones_promedio), 0)
which.min(precipitaciones_promedio)
# Anio con temperatura promedio mayor
round(max(precipitaciones_promedio), 0)
which.max(precipitaciones_promedio)



# ggplot por epoca y año  -----

temp_range <- calculate_range(statistics_province[,,,1,])
rain_range <- calculate_range(statistics_province[,,,2,])
dewpt_range <- calculate_range(statistics_province[,,,3,])


pp_ggplt <-
  ggplot() + 
  geom_stars(data = statistics_province[,,,2,]) +
  geom_sf(data = my_polygons_province_sf, fill = NA, size = 0.2) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  facet_wrap(. ~ time,
             ncol = 5, 
             labeller = estacion_label) +
  scale_fill_gradientn(limits = rain_range,
                       colours = cptcity::cpt(pal = "colo_alpen_sky_rain_brolly")) +
  labs(fill = "Precipitación\n(mm)",
       y = NULL, x = NULL) +
  theme_void() +
  theme(legend.position = 'bottom',
        legend.key.width = unit(1.5,"cm"))

ggsave('images/pp_roi.png',
       pp_ggplt,
       height = 14.52,
       width = 15.00,
       units = 'cm')


ptoro_ggplt <-
  ggplot() + 
  geom_stars(data = statistics_province[,,,3,]) +
  geom_sf(data = my_polygons_province_sf, fill = NA, size = 0.2) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  facet_wrap(. ~ time,
             ncol = 4, 
             labeller = estacion_label) +
  scale_fill_gradientn(limits = dewpt_range,
                       colours = cptcity::cpt(pal = "idv_relative_humidity")) +
  labs(fill = "Punto de rocío\n(°C)",
       y = NULL, x = NULL) +
  theme_void() +
  theme(legend.position = 'bottom',
        legend.key.width = unit(1.5,"cm"))

ggsave('images/puntorocio_roi.png',
       ptoro_ggplt,
       height = 14.52,
       width = 15.00,
       units = 'cm')
