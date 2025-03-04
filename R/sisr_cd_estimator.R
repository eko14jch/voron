# Estimates the number of length-n SAWs in Z^d 
sisr_cd_estimator <- function(walk_length, number_of_walks, d){
  
  # Initialize empty matrices to store the weights and sequences. The sequence array will be 3D.
  weights <- numeric(number_of_walks)
  
  sequences <- array(NA, 
                     dim = c(walk_length + 1,  # n
                             d,                # the number of dimensions
                             number_of_walks)) # i
  
  # Set initial weights and positions
  weights[] <- 1
  sequences[1, , ] <- rep(0, d) # All sequences start at the origin 
  
  
  for (k in 2:(walk_length + 1)){ # For each step 
    #cat("Step", k, "\n")
    for (i in seq_len( number_of_walks ) ){  # Iterate through each of the sequences and update for that step
      #cat("Sequence number ", i, "\n") 
      
      current_location <- sequences[k-1, , i] # First check the current location: k-1th step, both x and y, the ith sequence
      
      possible_next_locations <- matrix(rep(current_location, each = 2*d), ncol = d, byrow = FALSE)
      # Adjust each coordinate +1 and -1
      for (j in 1:d) {
        possible_next_locations[(2*j - 1), j] <- possible_next_locations[(2*j - 1), j] + 1  # Move +1 in j-th dimension
        possible_next_locations[(2*j), j] <- possible_next_locations[(2*j), j] - 1  # Move -1 in j-th dimension
      }
      
      
      # Extract visited positions as a matrix
      visited_positions <- matrix(sequences[1:(k-1), , i, drop = TRUE], ncol = d)
      
      # Self-avoiding check by removing possible next locations that we already visited
      all_positions <- rbind(visited_positions, possible_next_locations) # m x n, i.e. 2 dimensions
      allowed_next_locations <- possible_next_locations[ !duplicated(all_positions)[-(1:(k-1))], , drop = FALSE]
      
      number_of_SA_steps <- nrow(allowed_next_locations)
      
      
      # Select next step or stay in place if trapped
      if(number_of_SA_steps == 0){ # If there are no self-avoiding steps, stay at the same location
        next_location <- current_location
      } else { # Else, draw uniformly among the self-avoiding steps
        next_location <- allowed_next_locations[sample(number_of_SA_steps, 1), ]
      }
      
      # Update the next location
      sequences[k, ,i] <- next_location
      
      # Update weights
      #If there are no self-avoiding moves, the weight is set to 0,
      
      updated_weight <- number_of_SA_steps * weights[i]
      weights[i] <- updated_weight
    }
    
    # Compute effective sample size to check if we need to resample
    N_eff <- sum(weights)^2 / sum(weights^2) # 
    if (N_eff < number_of_walks/2){ # If the effective sample size is smaller than half of the sample size, we resample
      
      #cat("Resampling after step", k, "\n")
      
      # Ensure valid probabilities
      if (sum(weights) == 0) {
        probs <- rep(1 / number_of_walks, number_of_walks)
      } else {
        probs <- weights / sum(weights)
        
        # Resample
        resampled_columns <- sample(1:number_of_walks, size = number_of_walks, replace = TRUE, prob = probs)
        sequences <- sequences[ , , resampled_columns, drop = FALSE]
        weights[] <- mean(weights) # Normalize weights instead of resetting to 1
      } 
    }
  }
  estimated_c <- mean(weights)
  return(estimated_c)
}