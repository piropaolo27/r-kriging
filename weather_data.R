library(riem)
library(purrr)
library(dplyr)

networks = riem_networks()
stations <- riem_stations(network = "FR__ASOS")
december_weather <- map_df(stations$id,
                           riem_measures,
                           date_start = "2019-12-06",
                           date_end = "2019-12-07")
write.csv(december_weather, "csv/december_weather_fr.csv")
