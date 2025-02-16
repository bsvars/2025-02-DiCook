library(bsvarSIGNs)

data <- readxl::read_excel("ecor12749-sup-0001-supinfo/VARData.xlsx")

Y <- as.matrix(data[, 2:5])
Z <- data[, 6:8] |>
  as.matrix() |>
  bsvars::specify_data_matrices$new(p = 4)
Z <- rbind(matrix(0, 4, 12), t(Z$X[-nrow(Z$X), ]))

sign_irf <- matrix(NA, 4, 4)
sign_irf[1, 1] <- sign_irf[4, 1] <- 1
sign_irf[3, 1] <- -1
sign_irf <- array(sign_irf, c(4, 4, 4))

sign_structural <- matrix(NA, 4, 4)
sign_structural[1, ] <- c(NA, -1, -1, 1)

spec <- specify_bsvarSIGN$new(
  Y,
  p = 4,
  exogenous = Z,
  sign_irf = sign_irf,
  sign_structural = sign_structural
)
post <- estimate(spec, S = 5000)
irf <- compute_impulse_responses(post, horizon = 24)
plot(irf, probability = 0.68)
