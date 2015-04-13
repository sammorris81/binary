#########################################################################
# A simulation study to determine if a thresholded or skew methods improve
# return-level estimation and threshold exceedance prediction over
# standard kriging methods
#
# data settings:
#   All: rho = 5, s in [0, 10] x [0, 10]
#   1: GEV link
#      a: alpha = 0.3, 100 knots, 1% rareness
#      b: alpha = 0.7, 100 knots, 1% rareness
#      c: alpha = 0.3, 100 knots, 5% rareness
#      d: alpha = 0.7, 100 knots, 5% rareness
#   2: Logit link
#      a: rho = 1, 100 knots, 1% rareness
#      b: rho = 3, 100 knots, 1% rareness
#      c: rho = 1, 100 knots, 5% rareness
#      d: rho = 3, 100 knots, 5% rareness
#   3: Probit link -- Hold off for now
#      a: gamma = 0.9, 100 knots, 1% rareness
#      b: gamma = 0.1, 100 knots, 1% rareness
#      c: gamma = 0.9, 100 knots, 5% rareness
#      d: gamma = 0.1, 100 knots, 5% rareness
#########################################################################