knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(capivara)

# library(capivara)
# 
# segments <- segment_kinematics(
#   cube_path = "/data/galaxy.fits",
#   redshift = 0.03,
#   emission_line = "halpha",
#   segmentation_mode = "kinematic",
#   knn_k = 50,
#   n_segments = 25,
#   show_plots = TRUE
# )

# result <- run_kinematic_analysis(
#   cube_path = "/data/galaxy.fits",
#   redshift = 0.03,
#   emission_line = "halpha",
#   model = "axisymmetric",
#   segmentation_mode = "kinematic",
#   knn_k = 50,
#   n_segments = 25,
#   show_plots = TRUE
# )
# 
# plot(result, which = "model")

kinematic_models()
