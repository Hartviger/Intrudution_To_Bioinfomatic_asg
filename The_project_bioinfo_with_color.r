#################################
# Introduction to Bioinformatics
# The Project
# Fall 2024
# Group 3
#################################

# Load necessary libraries
library(Biostrings)
library(DECIPHER)
library(crayon) 

# Set working directory
getwd()
setwd("/Users/jakobhartvig/Desktop/uni/7_semester/DM847 Intrudution To Bioinfomatic/Project/project_sequences")

# Define the query sequence as an AAStringSet
query_seq <- AAStringSet("QVQLQQSGAELARPGASVKMSCKASGYTFTRYTMHWVKQRPGQGLEWIGYINPSRGYTNYNQKFKDKATLTTDKSSSTAYMQLSSLTSEDSAVYYCARYYDDHYCLDYWGQGTTLTVSSAKTTAPSVYPLAPVCGGTTGSSVTLGCLVKGYFPEPVTLTWNSGSLSSGVHTFPAVLQSDLYTLSSSVTVTSSTWPSQSITCNVAHPASSTKVDKKIEPRPKSCDKTHTCPPCPAPELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSRDELTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK")

# Read the .fa files
human_sequences <- readAAStringSet("human.fa")
mouse_sequences <- readAAStringSet("mouse.fa")

# List of datasets to loop through
datasets <- list(human = human_sequences, mouse = mouse_sequences)



# Function to calculate percentage identity
calculate_percentage_identity <- function(aln) {
  # Extract aligned sequences
  aligned_query <- as.character(pattern(aln))
  aligned_target <- as.character(subject(aln))
  
  # Split sequences into individual characters
  query_chars <- strsplit(aligned_query, "")[[1]]
  target_chars <- strsplit(aligned_target, "")[[1]]
  
  # Calculate matches (excluding gaps)
  matches <- sum(query_chars == target_chars & query_chars != "-" & target_chars != "-")
  
  # Calculate total positions compared (excluding gaps)
  total <- sum(query_chars != "-" & target_chars != "-")
  
  # Handle cases where total is zero to avoid division by zero
  if (total == 0) {
    return(0)
  } else {
    return((matches / total) * 100)
  }
}

# Function to display alignments with color
display_alignment_color <- function(aln) {
  aligned_query <- unlist(strsplit(as.character(pattern(aln)), ""))
  aligned_target <- unlist(strsplit(as.character(subject(aln)), ""))
  
  # Generate colored sequences
  colored_query <- mapply(function(q, t) {
    if (q == t) {
      green(q)          # Matches in green
    } else if (q == "-" || t == "-") {
      yellow(q)         # Gaps in yellow
    } else {
      red(q)            # Mismatches in red
    }
  }, aligned_query, aligned_target)
  
  colored_target <- mapply(function(q, t) {
    if (q == t) {
      green(t)
    } else if (q == "-" || t == "-") {
      yellow(t)
    } else {
      red(t)
    }
  }, aligned_query, aligned_target)
  
  # Print colored sequences
  cat("Query  : ", paste(colored_query, collapse = ""), "\n", sep = "")
  cat("Target : ", paste(colored_target, collapse = ""), "\n", sep = "")
}


# Set the alignment type
alignment_type <- "local"  # Change to "global", "overlap", etc., as needed, check out "PairwiseAlignments-class {Biostrings}"



# Loop through each dataset
for (dataset_name in names(datasets)) {
  target_sequences <- datasets[[dataset_name]]
  
  cat("\n\n### Analyzing ", dataset_name, " Dataset ###\n\n", sep = "")
  
  # Combine the query sequence with the target sequences for alignment
  combined_sequences <- c(query_seq, target_sequences)
  alignment <- AlignSeqs(combined_sequences, anchor = NA)
  
  # Perform pairwise alignment to calculate scores for each sequence in the target set
  alignments <- lapply(target_sequences, function(seq) {
    pairwiseAlignment(query_seq, seq, 
                      substitutionMatrix = "BLOSUM62", 
                      gapOpening = -10,  # We need to discuss this value
                      gapExtension = -1, # As well as this
                      type = alignment_type) # Use the specified alignment type
  })
  
  # Calculate scores
  scores <- sapply(alignments, function(aln) score(aln))
  
  # Calculate percentage identities
  percentage_identities <- sapply(alignments, calculate_percentage_identity)
  
  # Extract sequence names if available
  target_names <- names(target_sequences)
  
  # Combine scores, percentage identities, and aligned sequences into a data frame
  alignment_stats <- data.frame(
    Index = 1:length(alignments),
    Name = target_names,
    Score = scores,
    PercentageIdentity = percentage_identities,
    AlignedQuery = sapply(alignments, function(aln) as.character(pattern(aln))),
    AlignedTarget = sapply(alignments, function(aln) as.character(subject(aln))),
    stringsAsFactors = FALSE
  )
  
  # Find top 5 best and top 5 worst alignments based on scores
  best_indices <- order(scores, decreasing = TRUE)[1:5]
  worst_indices <- order(scores, decreasing = FALSE)[1:5]
  
  # Display the top 5 best alignments with color and alignment type
  cat("\nTop 5 Best Alignments for ", dataset_name, " Dataset (Alignment Type: ", alignment_type, "):\n", sep = "")
  lapply(best_indices, function(i) {
    cat("\nAlignment ", " (", target_names[i], ") (Score: ", scores[i], 
        ", Percentage Identity: ", sprintf("%.2f%%", percentage_identities[i]), 
        ", Alignment Type: ", alignment_type, "):\n", sep = "")
    display_alignment_color(alignments[[i]])
  })
  
  # Display the top 5 worst alignments with color and alignment type
  cat("\nTop 5 Worst Alignments for ", dataset_name, " Dataset (Alignment Type: ", alignment_type, "):\n", sep = "")
  lapply(worst_indices, function(i) {
    cat("\nAlignment ", " (", target_names[i], ") (Score: ", scores[i], 
        ", Percentage Identity: ", sprintf("%.2f%%", percentage_identities[i]), 
        ", Alignment Type: ", alignment_type, "):\n", sep = "")
    display_alignment_color(alignments[[i]])
  })
}
