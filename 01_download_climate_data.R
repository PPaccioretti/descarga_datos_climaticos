# Carga paquetes -----
library(rgee)
library(stars)
library(tmap)
library(foreach)
library(dplyr)
library(ggplot2)
library(lubridate)

# Inizialización rgee ----
ee_Initialize(drive = TRUE)

# Definicion funciones ----
cambio_unidades <-  function(img) {

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

my_polygons_province <-
  my_polygons_country$
  filter(ee$Filter$inList('ADM1_NAME', 
                          list("La Pampa", "Buenos Aires", "Cordoba",  
                               "Santa Fe", "Entre Rios", "Santiago Del Estero",  
                               "Corrientes", "Salta", "Chaco",
                               "Tucuman", "Jujuy")
  )
  )$map(select_polygons)

my_polygons_province_sf <- ee_as_sf(my_polygons_province)
## Lectura archivo virus ----
my_deasease_data <- sf::st_read("data/virus.gpkg")

# Número de observaciones
cat("Número de observaciones: ", nrow(my_deasease_data))

# Número de sitios
n_sitios <- 
  my_deasease_data %>%
  sf::st_geometry() %>% 
  sf::st_as_sf() %>% 
  dplyr::distinct() %>% 
  nrow()
cat("Número de sitios observados: ", n_sitios)

## Union Archivos ----
my_deasease_province_data <- 
  st_join(my_deasease_data, 
          my_polygons_province_sf, 
          largest = TRUE)


my_polygons_country_sf <- ee_as_sf(my_polygons_country)

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
                             paste0("20/3/",Año.de.Colecta), 
                             paste0("21/9/",Año.de.Colecta)),
         start_date = ymd(as.Date(start_date, "%d/%m/%Y")),
         end_date = start_date %m+% months(3),
         end_date = ymd(as.Date(end_date, "%d/%m/%Y"))) %>% 
  rename(Anio_de_Colecta = Año.de.Colecta)
ee_table <- sf_as_ee(sf_table) # Carga de sf como FeatureCollection

# Otención de imágenes ---- 
min_date <- min(as.Date(sf_table$start_date))
max_date <- max(as.Date(sf_table$end_date))

era5_ee <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY")$
  filterBounds(my_polygons_province$geometry())$
  filterDate(as.character(min_date), as.character(max_date))$
  select(c('total_precipitation', 'temperature_2m', 'dewpoint_temperature_2m'))$
  map(cambio_unidades)

# https://github.com/csaybar/rgee/blob/examples//ImageCollection/overview.R
## Resolución del producto 
era5_resolution <- era5_ee$first()$projection()$nominalScale()
cat("ERA5 resolution: ", era5_resolution$getInfo(), '\n')

# Nombre de bandas
bandNames = era5_ee$first()$bandNames()$getInfo()
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
cat("Date range: ", as.character(col_min), as.character(col_max), '\n')

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

ggsave('images/Figura_5_presenciaAusencia_hpmwov.png',
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

ggsave('images/Figura_6_presenciaAusencia_wsmv.png',
       presAus_wsmv_plt,
       height = 12.0,
       width = 15.0,
       units = 'cm')

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


# Media de temperatura para todos los años muestreados -----
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

ggsave('images/Figura_2_VariabilidadEspacialTemperaturaPromedio.png',
       temperatura_promedio_ggplt,
       height = 10,
       width = 15,#10
       units = 'cm')

library(patchwork)
pp_ptoro_prom_ggplt <-
  pp_promedio_ggplt +
  theme(#legend.position = 'right',
    legend.key.width = unit(0.6, "cm"),
    legend.box = "horizontal") +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5,
                                  
                                direction = "horizontal",
                                # title.position = "top",
                                label.position = "bottom",
                                label.hjust = 0.5,
                                # label.vjust = 1,
                                label.theme = element_text(angle = 90)
  )) +
   plot_spacer() +
  ptoro_promedio_ggplt +
  theme(legend.position = 'right',
    legend.key.width = unit(0.7, "cm"),
    legend.box = "horizontal") +
  scale_fill_continuous(
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      label.position = "bottom",
      label.hjust = 0.5,
      label.vjust = 1,
      label.theme = element_text(angle = 90)
    )
  ) + plot_layout(widths = c(4, 1 ,4))


pp_ptoro_prom_ggplt <-
  pp_promedio_ggplt +
  theme(
    #legend.position = 'right',
    legend.key.width = unit(0.9, "cm"),
    legend.box = "horizontal",
    legend.text = element_text(
      angle = 90,
      hjust = 0.3,
      vjust = 0.4
    )
  ) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) +
  plot_spacer() +
  ptoro_promedio_ggplt +
  theme(#legend.position = 'right',
    legend.key.width = unit(0.9, "cm"),
    legend.box = "horizontal",) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5)) + 
  plot_layout(widths = c(5, 0.5, 5))


ggsave('images/Figura_3_PrecipitacionesTemperaturaPtoRocioPromedio.png',
       pp_ptoro_prom_ggplt,
       height = 12,
       width = 15.0,
       units = 'cm')




## Figuras no presentadas ----
pp_promedio_ggplt <-
  ggplot() + 
  geom_stars(data = statistics_otono[,,,2]) +
  geom_sf(data = my_polygons_province_sf, fill = NA) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  scale_fill_gradientn(limits = rain_range,
                       colours = cptcity::cpt(pal = "colo_alpen_sky_rain_brolly")) +
  labs(fill = "Precipitación\n(mm)",
       y = NULL, x = NULL) +
  theme(legend.position = 'bottom',
        legend.key.width = unit(1.2,"cm"), 
        axis.text.x = element_text(angle = 90))

ggsave('images/pp_promedio_roi.png',
       pp_promedio_ggplt,
       height = 10,
       width = 15.0/2,
       units = 'cm')


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
  theme(legend.position = 'bottom',
        legend.key.width = unit(1.2,"cm"), 
        axis.text.x = element_text(angle = 90))

ggsave('images/puntorocio_promedio_roi.png',
       ptoro_promedio_ggplt,
       height = 10,
       width = 15.86/2,
       units = 'cm')


## Temperatura promedio ----
temperatura_ggplt <- 
  ggplot() +
  geom_stars(data = statistics_province[, , , 1, ]) +
  geom_sf(data = my_polygons_province_sf, fill = NA, size = 0.2) +
  coord_sf(lims_method = "geometry_bbox", expand = FALSE) +
  facet_wrap(. ~ time,
             ncol = 5,
             labeller = estacion_label) +
  scale_fill_gradientn(limits = temp_range,
                       colours = cptcity::cpt(pal = "arendal_temperature")) +
  labs(fill = "Temperatura\n(°C)",
       y = NULL, x = NULL) +
  # theme_void() +
  theme(
    # legend.position = 'bottom',
    # legend.key.width = unit(1.5, "cm"),
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

ggsave('images/Figura_4_VariabilidadEspacioTemporalTemperaturaPromedio.png',
       temperatura_ggplt,
       height = 14.52,
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
temp_temporal_media[1,order(temp_temporal_media[1,])]
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
