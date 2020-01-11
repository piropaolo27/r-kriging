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
library(sp)
library(dismo)
library(gstat)

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

# sampling
n <- 5
samples <- sample(1:nrow(X), n)
probe_x = X[samples, 1]
probe_y = X[samples, 2]
probe_z = Z[samples, 1]

conv <- convHull(X)
hull_pol <- conv@polygons@polygons[[1]]@Polygons[[1]]@coords

verts_x <- hull_pol[,1]
verts_y <- hull_pol[,2]

# fitting function
fit_into_hull <- function(m, Y1, Y2, verts_x, verts_y, W) {
  for (i in 1:m[1]) {
    for (j in 1:m[2])
      if (!point.in.polygon(Y1[i], Y2[j], verts_x, verts_y))
        W[j, i] = NA
  }
  return(W)
}

m <- c(200, 200)
lon_lim <- c(min(X[, 1]), max(X[, 1]))
lat_lim <- c(min(X[, 2]), max(X[, 2]))

lon_loc <- c(0.1, 0.1)
lat_loc <- c(0.1, 0.1)
lon_lim <- c(lon_lim[1] - lon_loc[1], lon_lim[2] + lon_loc[2])
lat_lim <- c(lat_lim[1] - lat_loc[1], lat_lim[2] + lat_loc[2])

# generate regular sequences
Y1 <- seq(lon_lim[1], lon_lim[2], length.out = m[1])
Y2 <- seq(lat_lim[1], lat_lim[2], length.out = m[2])

# add probe points
m <- c(200 + n, 200 + n)
X <- data.matrix(X[-samples, 1:2])
Z <- data.matrix(Z[-samples])
Y1 <- sort(append(Y1, probe_x))
Y2 <- sort(append(Y2, probe_y))

# verify function
verify <- function(n, W, Y1, Y2, probe_x, probe_y, probe_z) {
  for (i in 1:n) {
    cat("Predicted temperature:", W[which(Y2 == probe_y[i]), which(Y1 == probe_x[i])], "\n")
    cat("Actual temperature:", probe_z[i], "\n")
  }
}

# create a data frame from all combinations of factor variables
grid <- expand.grid(x=Y1, y=Y2)
Y <- as.matrix(grid)

# setting the kernel
Kern <- function(x, y) {
  h <- sqrt(sum((x - y) ^ 2))
  return(exp_kernel(h, 0.01, 3))
}

# computing covariance matrices
K = Kov(X, X, Kern, TRUE)
k = Kov(Y, X, Kern)

G = matrix(0, 1, 1)
g = matrix(0, 1, 1)

# plotting level curves results function
plot_japan <- function(x, y, z) {
  cols <- matlab.like2(40)
  plot_ly(
    x = x,
    y = y,
    z = z,
    type = "contour",
    colors = cols,
    contours = list(
      start = 0,
      size = 1,
      end = 30,
      showlabels = TRUE
    )
  )
}


# simple kriging
KRIG <- Krig(
  Z = Z,
  K = K,
  k = k,
  G = G,
  g = g,
  type = "simple",
  cinv = "syminv"
)

W <- matrix(KRIG$Z, m[1], m[2])
W <- fit_into_hull(m, Y1, Y2, verts_x, verts_y, W)
verify(n, W, Y1, Y2, probe_x, probe_y, probe_z)
# plot_japan(Y1, Y2, W)


# ordinary kriging
KRIG <- Krig(
  Z = Z,
  K = K,
  k = k,
  G = G,
  g = g,
  type = "ordinary",
  cinv = "syminv"
)

W <- matrix(KRIG$Z, m[1], m[2])
W <- fit_into_hull(m, Y1, Y2, verts_x, verts_y, W)
verify(n, W, Y1, Y2, probe_x, probe_y, probe_z)
# plot_japan(Y1, Y2, W)


# kriging universal
G <- rbind(t(X), t(X * X), t(X * X * X))
g <- rbind(t(Y), t(Y * Y), t(Y * Y * Y))

KRIG <- Krig(
  Z = Z,
  K = K,
  k = k,
  G = G,
  g = g,
  type = "universal",
  cinv = "syminv"
)

W <- matrix(KRIG$Z, m[1], m[2])
W <- fit_into_hull(m, Y1, Y2, verts_x, verts_y, W)
# plot_japan(Y1, Y2, W)


# idw
data <- data_frame(x = X[, 1], y = X[, 2], val = Z[, 1])
coordinates(data) =  ~ x + y
coordinates(grid) =  ~ x + y
predictions = idw(formula = val ~ 1,
                  locations = data,
                  newdata = grid)

W <- matrix(predictions$var1.pred, m[1], m[2])
W <- fit_into_hull(m, Y1, Y2, verts_x, verts_y, W)
# plot_japan(Y1, Y2, W)


# linear regression
data <- data_frame(x = X[, 1], y = X[, 2], val = Z[, 1])
model <- lm(val ~ (x + y), data = data)
predictions <- predict(model, new = grid)

W <- matrix(predictions, m[1], m[2])
W <- fit_into_hull(m, Y1, Y2, verts_x, verts_y, W)
# plot_japan(Y1, Y2, W)
