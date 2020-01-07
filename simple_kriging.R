library(riem)
library(purrr)
library(dplyr)
library(weathermetrics)
library(data.table)
library(KRIG)
library(plotly)
library(colorRamps)
library(rnaturalearth)
library(uniformly)

# read csv data
november_weather <- read.csv("csv/november_weather_jt.csv", stringsAsFactors = FALSE)
november_weather <- november_weather %>% select(station, lon, lat, tmpf)

# convert into celsius
november_weather <- november_weather %>%
                    mutate(tmpc = convert_temperature(tmpf,
                                                      old_metric = "f",
                                                      new_metric = "c"))

# calc average temperature and count measures
november_weather <- november_weather %>%
                    group_by(station, lon, lat) %>%
                    summarize(avg_tmpc = mean(tmpc, na.rm = TRUE), cnt_tmpc = n())

# convert tibbles to matrices
X <- data.matrix(november_weather[2:3])
Z <- data.matrix(november_weather[4])

japan <- ne_countries(country = "japan")
jpvert <- japan@polygons[[1]]@Polygons[[2]]@coords
points <- runif_in_polygon(n = 1000, vertices = jpvert)

m <- c(10, 10)
lon_lim <- c(min(X[, 1]), max(X[, 1]))
lat_lim <- c(min(X[, 2]), max(X[, 2]))

lon_loc <- c(2, 2)
lat_loc <- c(2, 2)
lon_lim <- c(lon_lim[1] - lon_loc[1], lon_lim[2] + lon_loc[2])
lat_lim <- c(lat_lim[1] - lat_loc[1], lat_lim[2] + lat_loc[2])

# generate regular sequences
Y1 <- seq(lon_lim[1], lon_lim[2], length.out = m[1])
Y2 <- seq(lat_lim[1], lat_lim[2], length.out = m[2])

# create a data frame from all combinations of factor variables
Y <- expand.grid(Y1, Y2)
Y <- as.matrix(Y)
point <- points
Y <- data.matrix(points)
Y1 <- points[,1]
Y2 <- points[,2]
dist <- function(x, y) {
  return(sqrt(sum((x - y) ^ 2)))
}

# calculate sample or residual variogram or variogram cloud
V <- variogram(Z, X, dist)
d <- V$distance[V$sort + 1, 1]

spherical_variogram <- function(d, s, t) {
  return(s - exp_kernel(d, s, t))
}

fit_spherical_kernel <- function(p) {
  FV <- sapply(d, FUN = spherical_variogram, p[1], p[2])
  return(sum((FV - V$variogram[, 1]) ^ 2))
}

NLM <- nlm(fit_spherical_kernel, c(18, 22))
str(NLM)

FV <- sapply(d, FUN = spherical_variogram, NLM$estimate[1], NLM$estimate[2])

# plotting variogram
plot(
  d,
  V$variogram[, 1],
  cex = 0.5,
  pch = 16,
  col = 'purple4',
  xlab = 'd',
  ylab = 'v'
)
points(d,
       FV,
       type = 'l',
       col = 'dodgerblue2' ,
       lwd = 2)

# setting the kernel based on parameters from the variogram
Kern <- function(x, y) {
  h <- sqrt(sum((x - y) ^ 2))
  return(exp_kernel(h, 0.01, 3))
}

# computing covariance matrices
K = Kov(X, X, Kern, TRUE)
k = Kov(Y, X, Kern)

G = matrix(0, 1, 1)
g = matrix(0, 1, 1)

# kriging
KRIG <- Krig(
  Z = Z,
  K = K,
  k = k,
  G = G,
  g = g,
  type = "simple",
  cinv = "syminv"
)

W1 <- KRIG$Z
plot_ly(x = Y1, y = Y2)
