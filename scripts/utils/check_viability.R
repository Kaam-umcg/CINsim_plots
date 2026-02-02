# Function to check the viability of simulations within a karyoSimParallel object
# viability is defined as the 
check_viability <- function(sim_list, sim_number = 1:100, threshold_value, max_g) {

    # this functions only takes karyoSimParallel 
    if (class(sim_list) != "karyoSimParallel"){
        message("Incorrect type for given object, needs to be karyoSimParallel")
        return(FALSE)
    }

    # Access the last value of the true_cell_count
    last_value <- sim_list[[paste0("sim_", sim_number)]]$gen_measures$true_cell_count[length(sim_list[[paste0("sim_", sim_number)]]$gen_measures$true_cell_count)]

    # Check if it's greater than the threshold
    max_cells_reached <- last_value > threshold_value

    # or if the max amount of generations was reached
    max_g_reached <- max(sim_list[[paste0("sim_", sim_number)]]$gen_measures$g) == max_g

    # if either were true, return TRUE  
    return(max_cells_reached | max_g_reached)
}