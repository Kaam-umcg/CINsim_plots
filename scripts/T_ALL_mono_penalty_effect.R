#' this script should make a plot where it loads in the
#' best performing simulation from T_ALL_surv_FC and T_ALL_mono_penalty
#' and compares the CnFS performance between the two and also check the amount
#' of non-modal monosomies between those populations.
#' We can also make a plot showing the effect of S_f (fitness) between
#' Having enabled the monosomy penalty or not.
#' We can use the CINsim::get_score() function for this
#' This should serve as proof that the monosomy penalty is applying as we
#' intend it to, and and has the effect of improving the CnFS. We
#' can compare those as well and do statistics on them 

# make this once the other simulations have finsihed running, might just want to generate
# most of the code for this automatically