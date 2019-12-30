library(riem)
library(purrr)
library(dplyr)

networks = riem_networks()
north_korea_stations <- riem_stations(network = "KP__ASOS")
south_korea_stations <- riem_stations(network = "KR__ASOS")
japan_stations <- riem_stations(network = "JP__ASOS")
joint_stations <- rbind(north_korea_stations,
                        south_korea_stations,
                        japan_stations)
november_weather <- map_df(joint_stations$id,
                           riem_measures,
                           date_start = "2019-11-01",
                           date_end = "2019-11-02")
write.csv(november_weather, "csv/november_weather_jt.csv")
